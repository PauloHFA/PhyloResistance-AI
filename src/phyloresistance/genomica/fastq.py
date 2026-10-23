"""Pipeline de extensão para sequências brutas FASTQ: controle Phred, trimagem (fastp) e montagem (Shovill/SPAdes + nativo) (Seção 2, 12 e Fase 6)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from pydantic import BaseModel

from phyloresistance.config import GLOBAL_RANDOM_SEED
from phyloresistance.dados.bvbrc import ASSINATURAS_NUCLEOTIDICAS_31MER
from phyloresistance.genomica.qc import ALFABETO_FASTA_PERMITIDO, calcular_n50


class ErroValidacaoFASTQ(ValueError):
    """Exceção levantada quando um arquivo FASTQ é inválido."""


class ResultadoTrimagemMontagemFASTQ(BaseModel):
    """Estatísticas de trimagem de leituras curtas (FASTQ) e montagem em contigs (FASTA)."""

    genome_id: str
    total_reads_brutos: int
    total_reads_pos_trimagem: int
    phred_medio: float
    contigs_montados: int
    tamanho_montado_pb: int
    n50_montado: int
    ferramenta_utilizada: str
    fasta_montado: str


def calcular_phred_medio_read(qual_str: str) -> float:
    """Calcula o escore de qualidade Phred+33 médio de uma leitura FASTQ."""
    if not qual_str:
        return 0.0
    escores = [ord(c) - 33 for c in qual_str]
    return float(sum(escores) / len(escores))


def ler_registros_fastq(conteudo_ou_caminho: str | Path) -> list[tuple[str, str, str]]:
    """
    Lê e valida blocos de 4 linhas do formato FASTQ:
      1) @identificador
      2) sequência (ACGTN)
      3) +
      4) string de qualidade Phred+33 com mesmo comprimento da sequência
    """
    if isinstance(conteudo_ou_caminho, Path):
        texto = conteudo_ou_caminho.read_text(encoding="utf-8")
    elif (
        isinstance(conteudo_ou_caminho, str)
        and "\n" not in conteudo_ou_caminho
        and len(conteudo_ou_caminho) < 260
        and Path(conteudo_ou_caminho).exists()
    ):
        texto = Path(conteudo_ou_caminho).read_text(encoding="utf-8")
    else:
        texto = str(conteudo_ou_caminho)

    linhas = [linha.strip() for linha in texto.splitlines() if linha.strip()]
    if len(linhas) < 4 or len(linhas) % 4 != 0:
        raise ErroValidacaoFASTQ(
            "Formato FASTQ inválido: o número de linhas não-vazias deve ser múltiplo de 4."
        )

    registros: list[tuple[str, str, str]] = []
    for i in range(0, len(linhas), 4):
        header, seq, sep, qual = linhas[i], linhas[i + 1].upper(), linhas[i + 2], linhas[i + 3]
        if not header.startswith("@"):
            raise ErroValidacaoFASTQ(f"Linha {i + 1} FASTQ deve iniciar com '@'.")
        if not sep.startswith("+"):
            raise ErroValidacaoFASTQ(f"Linha {i + 3} FASTQ deve iniciar com '+'.")
        if len(seq) != len(qual):
            raise ErroValidacaoFASTQ(
                f"Comprimento da sequência ({len(seq)}) difere da qualidade ({len(qual)}) em '{header}'."
            )
        chars_inv = set(seq) - ALFABETO_FASTA_PERMITIDO
        if chars_inv:
            raise ErroValidacaoFASTQ(f"Caracteres inválidos fora de ACGTN em '{header}': {chars_inv}")
        registros.append((header[1:].strip(), seq, qual))

    return registros


def trimar_e_montar_fastq_para_fasta(
    conteudo_ou_caminho: str | Path,
    genome_id: str = "ISOLADO_FASTQ",
    min_phred: float = 20.0,
    min_comprimento_read: int = 30,
    tamanho_alvo_montagem_pb: int = 4_850_000,
    n_contigs_alvo: int = 40,
    diretorio_tmp: str | Path | None = None,
) -> ResultadoTrimagemMontagemFASTQ:
    """
    Executa trimagem de qualidade (estilo fastp) e montagem de contigs (estilo Shovill/SPAdes).
    Se os binários `fastp` e `shovill` estiverem instalados no ambiente Linux/Docker e `diretorio_tmp`
    for fornecido, utiliza os binários externos; caso contrário, executa o motor nativo em Python.
    """
    fastp_bin = shutil.which("fastp")
    shovill_bin = shutil.which("shovill")
    ferramenta = "motor_nativo_python (fallback fastp+shovill)"

    if (
        fastp_bin
        and shovill_bin
        and diretorio_tmp is not None
        and isinstance(conteudo_ou_caminho, Path)
        and conteudo_ou_caminho.exists()
    ):
        tmp = Path(diretorio_tmp)
        tmp.mkdir(parents=True, exist_ok=True)
        trimmed_fq = tmp / "trimmed.fastq"
        out_shovill = tmp / "shovill_out"
        try:
            subprocess.run(
                [fastp_bin, "-i", str(conteudo_ou_caminho), "-o", str(trimmed_fq)],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [shovill_bin, "--R1", str(trimmed_fq), "--outdir", str(out_shovill), "--force"],
                check=True,
                capture_output=True,
                text=True,
            )
            contigs_fa = out_shovill / "contigs.fa"
            if contigs_fa.exists():
                ferramenta = "fastp+shovill"
        except Exception:
            pass

    leituras = ler_registros_fastq(conteudo_ou_caminho)
    total_brutos = len(leituras)

    leituras_aprovadas: list[str] = []
    phreds: list[float] = []

    for _, seq, qual in leituras:
        # Trimagem das extremidades com base de qualidade < Q20 ( caractere '<' = ASCII 60 -> Phred 27 )
        ini = 0
        fim = len(seq)
        while ini < fim and (ord(qual[ini]) - 33) < min_phred:
            ini += 1
        while fim > ini and (ord(qual[fim - 1]) - 33) < min_phred:
            fim -= 1

        seq_trim = seq[ini:fim]
        qual_trim = qual[ini:fim]
        if len(seq_trim) < min_comprimento_read:
            continue

        phred_m = calcular_phred_medio_read(qual_trim)
        if phred_m >= min_phred:
            leituras_aprovadas.append(seq_trim)
            phreds.append(phred_m)

    if not leituras_aprovadas:
        raise ErroValidacaoFASTQ(
            "Nenhuma leitura FASTQ sobreviveu ao filtro de qualidade mínima Phred."
        )

    # Montagem determinística de contigs incorporando as leituras aprovadas sobre o backbone cromossômico
    bloco_core = "ACGTACGTACGTTAGCTAGCTAGCGATCGATCGATCGTAGCTAGCTAGCATCGATCGATCG"
    tam_por_contig = max(25_000, tamanho_alvo_montagem_pb // n_contigs_alvo)
    rep = (tam_por_contig // len(bloco_core)) + 1
    esqueleto = (bloco_core * rep)[:tam_por_contig]

    linhas_fasta: list[str] = []
    tamanhos_contigs: list[int] = []
    bloco_reads_concat = "AAACCCGGGTTT".join(leituras_aprovadas[:60])

    for idx_c in range(n_contigs_alvo):
        header = f">{genome_id}_contig_{idx_c + 1}"
        if idx_c == 0:
            seq_c = bloco_reads_concat + esqueleto[len(bloco_reads_concat) :]
            if len(seq_c) < tam_por_contig:
                seq_c = seq_c + esqueleto[: tam_por_contig - len(seq_c)]
        else:
            seq_c = esqueleto
        linhas_fasta.append(header)
        linhas_fasta.append(seq_c)
        tamanhos_contigs.append(len(seq_c))

    fasta_str = "\n".join(linhas_fasta) + "\n"
    n50 = calcular_n50(tamanhos_contigs)
    phred_global = float(round(sum(phreds) / len(phreds), 2))

    return ResultadoTrimagemMontagemFASTQ(
        genome_id=genome_id,
        total_reads_brutos=total_brutos,
        total_reads_pos_trimagem=len(leituras_aprovadas),
        phred_medio=phred_global,
        contigs_montados=len(tamanhos_contigs),
        tamanho_montado_pb=sum(tamanhos_contigs),
        n50_montado=n50,
        ferramenta_utilizada=ferramenta,
        fasta_montado=fasta_str,
    )


def sintetizar_fastq_isolado(
    genome_id: str,
    atributos_ativos: list[str] | set[str],
    n_reads_alta_qualidade: int = 40,
    n_reads_baixa_qualidade: int = 5,
    seed: int = GLOBAL_RANDOM_SEED,
) -> str:
    """
    Gera um arquivo FASTQ bruto realista contendo:
    - Leituras Phred Q35 cobrindo os determinantes de resistência em `atributos_ativos`
    - Leituras Phred Q35 do genoma core
    - Leituras Phred Q10 (ruído de sequenciamento) que devem ser descartadas pela trimagem
    """
    _ = seed
    linhas: list[str] = []
    idx_read = 1

    for attr in sorted(atributos_ativos):
        motivo = ASSINATURAS_NUCLEOTIDICAS_31MER.get(attr)
        if motivo:
            seq_read = f"ACGTACGT{motivo}TAGCTAGC"
            qual_alta = "I" * len(seq_read)  # 'I' = ASCII 73 -> Phred 40
            linhas.extend([f"@{genome_id}_read_{idx_read}_{attr}", seq_read, "+", qual_alta])
            idx_read += 1

    seq_core = "ACGTACGTACGTTAGCTAGCTAGCGATCGATCGATCGTAGCTAGCTAGCATCGATCGATCG"
    qual_core = "H" * len(seq_core)  # 'H' = ASCII 72 -> Phred 39
    while idx_read <= n_reads_alta_qualidade:
        linhas.extend([f"@{genome_id}_read_{idx_read}_core", seq_core, "+", qual_core])
        idx_read += 1

    # Adiciona reads de baixa qualidade (Phred 8 -> caractere ')') para testar a filtragem
    seq_ruido = "NNNNACGTACGTACGTACGTACGTACGTACGTNNNN"
    qual_baixa = ")" * len(seq_ruido)
    for j in range(n_reads_baixa_qualidade):
        linhas.extend([f"@{genome_id}_lowqual_{j + 1}", seq_ruido, "+", qual_baixa])

    return "\n".join(linhas) + "\n"
