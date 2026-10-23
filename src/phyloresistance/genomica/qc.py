"""Controle de qualidade de genomas montados (FASTA), cálculo de N50 e confirmação de espécie (Seção 3.4, REG-07)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pydantic import BaseModel, Field

from phyloresistance.config import (
    CRITERIOS_QC_POR_ESPECIE,
    ESPECIE_PILOTO,
    CriterioQCEspecie,
)

ALFABETO_FASTA_PERMITIDO: frozenset[str] = frozenset({"A", "C", "G", "T", "N"})


class ErroValidacaoFASTA(ValueError):
    """Exceção levantada quando um arquivo enviado não é um FASTA bacteriano válido."""


class ResultadoQCGenoma(BaseModel):
    """Resultado consolidado do controle de qualidade de um genoma montado."""

    genome_id: str
    tamanho_pb: int
    contigs: int
    n50: int
    fracao_n: float = 0.0
    completude: float = 99.0
    contaminacao: float = 0.5
    especie_declarada: str = ESPECIE_PILOTO
    especie_inferida: str = ESPECIE_PILOTO
    aprovado: bool
    motivos_reprovacao: list[str] = Field(default_factory=list)


def ler_contigs_fasta(conteudo_ou_caminho: str | Path) -> list[tuple[str, str]]:
    """Lê cabeçalhos e sequências de um arquivo FASTA ou string em memória."""
    texto: str
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
    if not linhas or not linhas[0].startswith(">"):
        raise ErroValidacaoFASTA("Formato inválido: arquivo FASTA deve iniciar com cabeçalho '>'.")

    contigs: list[tuple[str, str]] = []
    header_atual: str | None = None
    partes_seq: list[str] = []

    for linha in linhas:
        if linha.startswith(">"):
            if header_atual is not None:
                seq = "".join(partes_seq).upper()
                if not seq:
                    raise ErroValidacaoFASTA(f"Contig '{header_atual}' está vazio.")
                contigs.append((header_atual, seq))
            header_atual = linha[1:].strip() or f"contig_{len(contigs) + 1}"
            partes_seq = []
        else:
            partes_seq.append(linha)

    if header_atual is not None:
        seq = "".join(partes_seq).upper()
        if not seq:
            raise ErroValidacaoFASTA(f"Contig '{header_atual}' está vazio.")
        contigs.append((header_atual, seq))

    if not contigs:
        raise ErroValidacaoFASTA("Nenhum contig encontrado no arquivo FASTA.")

    return contigs


def calcular_n50(tamanhos_contigs: list[int]) -> int:
    """Calcula a estatística N50 de uma montagem a partir dos comprimentos dos contigs."""
    if not tamanhos_contigs:
        return 0
    ordenados = sorted(tamanhos_contigs, reverse=True)
    metade = sum(ordenados) / 2.0
    acumulado = 0
    for tam in ordenados:
        acumulado += tam
        if acumulado >= metade:
            return int(tam)
    return int(ordenados[-1])


def calcular_estatisticas_fasta(conteudo_ou_caminho: str | Path) -> dict[str, int | float]:
    """
    Calcula tamanho total (pb), nº de contigs, N50 e fração de bases N,
    validando o alfabeto nucleotídico (ACGTN).
    """
    contigs = ler_contigs_fasta(conteudo_ou_caminho)
    tamanhos: list[int] = []
    total_n = 0

    for header, seq in contigs:
        caracteres_invalidos = set(seq) - ALFABETO_FASTA_PERMITIDO
        if caracteres_invalidos:
            chars_str = ", ".join(sorted(caracteres_invalidos)[:5])
            raise ErroValidacaoFASTA(
                f"Caractere(s) fora do alfabeto ACGTN no contig '{header}': {chars_str}"
            )
        tamanhos.append(len(seq))
        total_n += seq.count("N")

    tamanho_total = sum(tamanhos)
    n50 = calcular_n50(tamanhos)
    fracao_n = (total_n / tamanho_total) if tamanho_total > 0 else 0.0

    return {
        "tamanho_pb": tamanho_total,
        "contigs": len(contigs),
        "n50": n50,
        "fracao_n": round(fracao_n, 6),
    }


def validar_fasta_entrada(
    conteudo_ou_caminho: str | Path,
    limite_mb: float = 25.0,
    max_contigs: int = 2000,
    min_pb: int = 100,
) -> dict[str, int | float]:
    """
    Valida formato FASTA, alfabeto ACGTN, tamanho máximo e número de contigs na entrada da API (Seção 9.2/9.4).
    """
    if isinstance(conteudo_ou_caminho, Path) or (
        isinstance(conteudo_ou_caminho, str)
        and "\n" not in conteudo_ou_caminho
        and len(conteudo_ou_caminho) < 260
        and Path(conteudo_ou_caminho).exists()
    ):
        caminho = Path(conteudo_ou_caminho)
        tamanho_bytes = caminho.stat().st_size
    else:
        tamanho_bytes = len(str(conteudo_ou_caminho).encode("utf-8"))

    if tamanho_bytes > limite_mb * 1024 * 1024:
        raise ErroValidacaoFASTA(
            f"Arquivo excede o limite máximo de {limite_mb} MB ({tamanho_bytes / (1024 * 1024):.2f} MB)."
        )

    stats = calcular_estatisticas_fasta(conteudo_ou_caminho)
    if int(stats["tamanho_pb"]) < min_pb:
        raise ErroValidacaoFASTA(
            f"Sequência muito curta ({stats['tamanho_pb']} pb < mínimo de {min_pb} pb)."
        )
    if int(stats["contigs"]) > max_contigs:
        raise ErroValidacaoFASTA(
            f"Número de contigs ({stats['contigs']}) excede o máximo permitido ({max_contigs})."
        )
    return stats


def avaliar_qc_genoma(
    genome_id: str,
    tamanho_pb: int,
    contigs: int,
    n50: int,
    fracao_n: float = 0.0,
    completude: float = 99.0,
    contaminacao: float = 0.5,
    especie_declarada: str = ESPECIE_PILOTO,
    especie_inferida: str = ESPECIE_PILOTO,
    criterio: CriterioQCEspecie | None = None,
) -> ResultadoQCGenoma:
    """
    Aplica os gates de controle de qualidade da Seção 3.4 e REG-07:
    - Espécie suportada e espécie declarada == espécie inferida
    - Tamanho total dentro da faixa da espécie
    - Número de contigs <= limite
    - N50 >= limite
    - Completude e contaminação dentro dos limiares CheckM2
    """
    motivos: list[str] = []

    if especie_declarada not in CRITERIOS_QC_POR_ESPECIE and criterio is None:
        motivos.append(f"Espécie '{especie_declarada}' fora do domínio suportado.")
        crit = CriterioQCEspecie(especie=especie_declarada)
    else:
        crit = criterio or CRITERIOS_QC_POR_ESPECIE[especie_declarada]

    if especie_declarada != especie_inferida:
        motivos.append(
            f"Divergência de espécie: declarada '{especie_declarada}' vs inferida '{especie_inferida}'."
        )

    if not (crit.tamanho_min_pb <= tamanho_pb <= crit.tamanho_max_pb):
        motivos.append(
            f"Tamanho total ({tamanho_pb} pb) fora da faixa [{crit.tamanho_min_pb}, {crit.tamanho_max_pb}]."
        )

    if contigs > crit.max_contigs:
        motivos.append(f"Excesso de contigs ({contigs} > {crit.max_contigs}).")

    if n50 < crit.min_n50_pb:
        motivos.append(f"N50 insuficiente ({n50} pb < {crit.min_n50_pb} pb).")

    if fracao_n > crit.max_fracao_n:
        motivos.append(f"Fração de bases ambíguas N elevada ({fracao_n:.4f} > {crit.max_fracao_n}).")

    if completude < crit.min_completude:
        motivos.append(f"Completude baixa ({completude:.1f}% < {crit.min_completude:.1f}%).")

    if contaminacao > crit.max_contaminacao:
        motivos.append(f"Contaminação elevada ({contaminacao:.1f}% > {crit.max_contaminacao:.1f}%).")

    return ResultadoQCGenoma(
        genome_id=genome_id,
        tamanho_pb=tamanho_pb,
        contigs=contigs,
        n50=n50,
        fracao_n=fracao_n,
        completude=completude,
        contaminacao=contaminacao,
        especie_declarada=especie_declarada,
        especie_inferida=especie_inferida,
        aprovado=(len(motivos) == 0),
        motivos_reprovacao=motivos,
    )


def avaliar_qc_coorte(df_genomas_meta: pd.DataFrame) -> tuple[pd.DataFrame, list[ResultadoQCGenoma]]:
    """Executa o QC para todos os isolados da coorte e retorna DataFrame com coluna 'aprovado'."""
    resultados: list[ResultadoQCGenoma] = []
    for row in df_genomas_meta.to_dict(orient="records"):
        res = avaliar_qc_genoma(
            genome_id=str(row["genome_id"]),
            tamanho_pb=int(row["tamanho_pb"]),
            contigs=int(row["contigs"]),
            n50=int(row["n50"]),
            fracao_n=float(row.get("fracao_n", 0.0)),
            completude=float(row.get("completude", 99.0)),
            contaminacao=float(row.get("contaminacao", 0.5)),
            especie_declarada=str(row.get("especie_declarada", ESPECIE_PILOTO)),
            especie_inferida=str(row.get("especie_inferida", ESPECIE_PILOTO)),
        )
        resultados.append(res)

    df_qc = pd.DataFrame([r.model_dump() for r in resultados])
    if "clone_de" in df_genomas_meta.columns:
        df_qc["clone_de"] = df_genomas_meta["clone_de"].values
    if "st" in df_genomas_meta.columns:
        df_qc["st"] = df_genomas_meta["st"].values
    if "pais" in df_genomas_meta.columns:
        df_qc["pais"] = df_genomas_meta["pais"].values
    if "ano" in df_genomas_meta.columns:
        df_qc["ano"] = df_genomas_meta["ano"].values
    return df_qc, resultados

