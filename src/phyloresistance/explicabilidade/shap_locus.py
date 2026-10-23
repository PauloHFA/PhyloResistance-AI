"""TreeSHAP, mapeamento de k-mers/unitigs para loci, agregação por locus e testes de validação XAI (Seção 8, REG-10)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import shap
import xgboost as xgb
from pydantic import BaseModel
from scipy.stats import spearmanr

from phyloresistance.config import (
    AVISO_USO_PESQUISA,
    CATALOGO_DETERMINANTES_ECOLI,
    CATALOGO_DETERMINANTES_POR_ESPECIE,
)
from phyloresistance.dados.bvbrc import MAPEAMENTO_LOCUS_CANDIDATOS


@dataclass(frozen=True)
class InfoLocusGenomico:
    """Metadados de mapeamento de um atributo (gene, mutação ou unitig/k-mer) para um locus genômico."""

    locus_id: str
    tipo: Literal["mutacao_pontual", "gene_amr", "kmer_unitig"]
    gene: str
    alteracao: str | None
    fonte: str
    status: Literal["conhecida", "candidata"]
    farmacos_card: tuple[str, ...]
    mecanismo: str


class EvidenciaGenomica(BaseModel):
    """Evidência genômica individual retornada no relatório por isolado e na API (Seção 9.3, REG-10)."""

    tipo: Literal["mutacao_pontual", "gene_amr", "kmer_unitig"]
    gene: str
    alteracao: str | None = None
    fonte: str
    contribuicao_shap: float
    status: Literal["conhecida", "candidata"]
    nota_associacao: str = (
        "Associação estatística ao modelo (SHAP); candidatos novos são hipóteses para validação experimental."
    )


class AvaliacaoExplicabilidadeQ3(BaseModel):
    """Resultados dos 4 testes obrigatórios de avaliação de explicações (Seção 8.2 e Q3)."""

    farmaco: str
    proporcao_top_k_no_card: float
    determinantes_card_recuperados: list[str]
    candidatos_novos_top_k: list[str]
    correlacao_estabilidade_entre_dobras: float
    concordancia_gwas_populacional: float
    queda_media_probabilidade_remocao_locus: float
    remocao_coerente: bool
    principais_loci_globais: list[dict[str, str | float]]


def construir_mapa_atributo_para_locus(
    nomes_atributos: list[str],
    farmaco: str | None = None,
    especie: str = "Escherichia coli",
) -> dict[str, InfoLocusGenomico]:
    """
    Mapeia cada atributo de F1, F2 e F3 para seu locus genômico e cruza com o catálogo
    CARD / AMRFinderPlus (Seção 8.1 passos 2 e 4).
    """
    todos_determinantes = list(CATALOGO_DETERMINANTES_ECOLI)
    for cat_esp in CATALOGO_DETERMINANTES_POR_ESPECIE.values():
        todos_determinantes.extend(cat_esp)
    catalogo_por_id = {d.id_atributo: d for d in todos_determinantes}
    if especie in CATALOGO_DETERMINANTES_POR_ESPECIE:
        for d in CATALOGO_DETERMINANTES_POR_ESPECIE[especie]:
            catalogo_por_id[d.id_atributo] = d
    mapa: dict[str, InfoLocusGenomico] = {}

    for nome in nomes_atributos:
        if nome in catalogo_por_id:
            det = catalogo_por_id[nome]
            eh_conhecida_farmaco = (farmaco is None) or (farmaco in det.farmacos_associados)
            status_locus: Literal["conhecida", "candidata"] = (
                "conhecida" if eh_conhecida_farmaco else "candidata"
            )
            locus_id = f"{det.gene}:{det.alteracao}" if det.alteracao else det.gene
            mapa[nome] = InfoLocusGenomico(
                locus_id=locus_id,
                tipo=det.tipo,  # type: ignore[arg-type]
                gene=det.gene,
                alteracao=det.alteracao,
                fonte=det.fonte,
                status=status_locus,
                farmacos_card=det.farmacos_associados,
                mecanismo=det.mecanismo,
            )
        elif nome in MAPEAMENTO_LOCUS_CANDIDATOS:
            info_c = MAPEAMENTO_LOCUS_CANDIDATOS[nome]
            locus_id = f"{info_c['gene']}:{info_c['alteracao']}"
            mapa[nome] = InfoLocusGenomico(
                locus_id=locus_id,
                tipo="kmer_unitig",
                gene=info_c["gene"],
                alteracao=info_c["alteracao"],
                fonte=info_c["fonte"],
                status="candidata",
                farmacos_card=(),
                mecanismo=info_c["mecanismo"],
            )
        elif nome.startswith("unitig_core_lineage_"):
            st_nome = nome.replace("unitig_core_lineage_", "").replace("_k31", "")
            mapa[nome] = InfoLocusGenomico(
                locus_id=f"intergenic_lineage_{st_nome}",
                tipo="kmer_unitig",
                gene=f"core_backbone_{st_nome}",
                alteracao="polimorfismo_linhagem",
                fonte="Alinhamento Referência (k=31)",
                status="candidata",
                farmacos_card=(),
                mecanismo="Marcador filogenético de linhagem (não-causal)",
            )
        else:
            mapa[nome] = InfoLocusGenomico(
                locus_id=f"locus_{nome}",
                tipo="kmer_unitig",
                gene=nome,
                alteracao="k31_variant",
                fonte="Alinhamento Referência (k=31)",
                status="candidata",
                farmacos_card=(),
                mecanismo="Variante genômica não catalogada (hipótese)",
            )
    return mapa


def calcular_valores_treeshap(
    modelo_xgb: xgb.XGBClassifier,
    X: np.ndarray,
) -> np.ndarray:
    """
    Calcula valores exatos TreeSHAP para o modelo de árvores XGBoost (Seção 8.1 passo 1).
    Retorna matriz (n_amostras, n_atributos) em escala de log-odds.
    """
    X_arr = np.asarray(X, dtype=np.float64)
    explicador = shap.TreeExplainer(modelo_xgb)
    valores = explicador.shap_values(X_arr)
    if isinstance(valores, list):
        return np.asarray(valores[-1], dtype=np.float64)
    arr = np.asarray(valores, dtype=np.float64)
    if arr.ndim == 3:
        return arr[:, :, -1]
    return arr


def agregar_shap_por_locus(
    valores_shap: np.ndarray,
    nomes_atributos: list[str],
    mapa_atributo_para_locus: dict[str, InfoLocusGenomico],
) -> tuple[np.ndarray, list[str], dict[str, InfoLocusGenomico]]:
    """
    Soma as contribuições SHAP de atributos que mapeiam para o mesmo locus genômico (gene ou mutação),
    pois k-mers vizinhos e atributos correlacionados no mesmo locus dividem o mesmo sinal (Seção 8.1 passo 3, REG-10).
    """
    shap_arr = np.asarray(valores_shap, dtype=np.float64)
    if shap_arr.ndim == 1:
        shap_arr = shap_arr.reshape(1, -1)

    ordem_loci: list[str] = []
    indices_por_locus: dict[str, list[int]] = {}
    info_por_locus: dict[str, InfoLocusGenomico] = {}

    for idx, nome in enumerate(nomes_atributos):
        info = mapa_atributo_para_locus[nome]
        lid = info.locus_id
        if lid not in indices_por_locus:
            ordem_loci.append(lid)
            indices_por_locus[lid] = []
            info_por_locus[lid] = info
        indices_por_locus[lid].append(idx)

    colunas_agregadas = [
        shap_arr[:, indices_por_locus[lid]].sum(axis=1) for lid in ordem_loci
    ]
    matriz_locus = np.column_stack(colunas_agregadas) if colunas_agregadas else np.zeros((shap_arr.shape[0], 0))
    return matriz_locus, ordem_loci, info_por_locus


def extrair_evidencias_isolado(
    vetor_x_isolado: np.ndarray,
    vetor_shap_isolado: np.ndarray,
    nomes_atributos: list[str],
    farmaco: str,
    top_k: int = 5,
    min_contribuicao_abs: float = 0.05,
) -> list[EvidenciaGenomica]:
    """
    Gera a lista ordenada de evidências genômicas que sustentam a predição de um isolado (Seção 8.1 passo 5 e 9.3).
    Inclui apenas atributos presentes no isolado (x_j > 0) ou de impacto relevante, agregados por locus.
    """
    x_arr = np.asarray(vetor_x_isolado, dtype=np.float64).ravel()
    s_arr = np.asarray(vetor_shap_isolado, dtype=np.float64).ravel()
    mapa = construir_mapa_atributo_para_locus(nomes_atributos, farmaco=farmaco)

    # Considera o sinal de atributos presentes no isolado (presença de gene/mutação/unitig)
    s_presentes = np.where(x_arr > 0, s_arr, 0.0)
    matriz_locus, ordem_loci, info_por_locus = agregar_shap_por_locus(
        s_presentes.reshape(1, -1), nomes_atributos, mapa
    )
    vetor_locus = matriz_locus.ravel()

    ordem_idx = np.argsort(np.abs(vetor_locus))[::-1]
    evidencias: list[EvidenciaGenomica] = []

    for idx_l in ordem_idx:
        val = float(vetor_locus[idx_l])
        if abs(val) < min_contribuicao_abs:
            continue
        lid = ordem_loci[idx_l]
        info = info_por_locus[lid]
        evidencias.append(
            EvidenciaGenomica(
                tipo=info.tipo,
                gene=info.gene,
                alteracao=info.alteracao,
                fonte=info.fonte,
                contribuicao_shap=round(val, 4),
                status=info.status,
            )
        )
        if len(evidencias) >= top_k:
            break

    return evidencias


def avaliar_recuperacao_determinantes_card(
    matriz_shap_locus: np.ndarray,
    ordem_loci: list[str],
    info_por_locus: dict[str, InfoLocusGenomico],
    top_k: int = 6,
) -> tuple[float, list[str], list[str]]:
    """
    Teste 1 da Seção 8.2 (Recuperação de determinantes conhecidos):
    Entre os k principais loci por importância média |SHAP|, calcula a proporção catalogada no CARD/AMRFinderPlus
    e lista também os candidatos novos descobertos.
    """
    importancia_media = np.mean(np.abs(matriz_shap_locus), axis=0)
    k_efetivo = min(top_k, len(ordem_loci))
    if k_efetivo == 0:
        return 0.0, [], []

    top_indices = np.argsort(importancia_media)[::-1][:k_efetivo]
    conhecidos: list[str] = []
    candidatos: list[str] = []

    for i in top_indices:
        lid = ordem_loci[int(i)]
        info = info_por_locus[lid]
        if info.status == "conhecida":
            conhecidos.append(lid)
        else:
            candidatos.append(lid)

    proporcao = len(conhecidos) / float(k_efetivo)
    return float(round(proporcao, 4)), conhecidos, candidatos


def avaliar_estabilidade_explicacoes(
    rankings_importancia_por_dobra: list[dict[str, float]],
) -> float:
    """
    Teste 2 da Seção 8.2 (Estabilidade):
    Calcula a correlação média de ranking de Spearman entre as importâncias SHAP por locus
    obtidas em dobras diferentes.
    """
    if len(rankings_importancia_por_dobra) < 2:
        return 1.0

    todos_loci = sorted(
        {lid for r in rankings_importancia_por_dobra for lid in r.keys()}
    )
    if len(todos_loci) < 2:
        return 1.0

    vetores = [
        np.array([r.get(lid, 0.0) for lid in todos_loci], dtype=np.float64)
        for r in rankings_importancia_por_dobra
    ]

    correlacoes: list[float] = []
    for i in range(len(vetores)):
        for j in range(i + 1, len(vetores)):
            rho, _ = spearmanr(vetores[i], vetores[j])
            if not np.isnan(rho):
                correlacoes.append(float(rho))

    return float(round(np.mean(correlacoes), 4)) if correlacoes else 0.0


def avaliar_concordancia_gwas_populacional(
    X: np.ndarray,
    y: np.ndarray,
    grupos_filogeneticos: list[str] | np.ndarray,
    nomes_atributos: list[str],
    mapa_atributo_para_locus: dict[str, InfoLocusGenomico],
    loci_top_shap: list[str],
    top_k: int = 6,
) -> float:
    """
    Teste 3 da Seção 8.2 (Comparação com GWAS bacteriano tipo pyseer):
    Calcula a estatística de associação estratificada por cluster filogenético (controle de estrutura
    populacional estilo Cochran-Mantel-Haenszel / resíduos por linhagem) e mede a concordância com os
    principais loci do TreeSHAP.
    """
    X_arr = np.asarray(X, dtype=np.float64)
    y_arr = np.asarray(y, dtype=np.float64)
    g_arr = np.asarray(grupos_filogeneticos, dtype=str)

    # Centraliza X e y dentro de cada cluster filogenético (fixed-effects lineage control)
    X_res = np.zeros_like(X_arr)
    y_res = np.zeros_like(y_arr)
    for g in np.unique(g_arr):
        m = g_arr == g
        X_res[m] = X_arr[m] - np.mean(X_arr[m], axis=0, keepdims=True)
        y_res[m] = y_arr[m] - np.mean(y_arr[m])

    numerador = np.abs(X_res.T @ y_res)
    denom = np.sqrt(np.sum(X_res**2, axis=0) * max(1e-9, float(np.sum(y_res**2))))
    denom = np.where(denom < 1e-9, 1.0, denom)
    escores_gwas_attr = numerador / denom

    # Agrega escores GWAS por locus
    escore_por_locus: dict[str, float] = {}
    for idx, nome in enumerate(nomes_atributos):
        lid = mapa_atributo_para_locus[nome].locus_id
        escore_por_locus[lid] = max(escore_por_locus.get(lid, 0.0), float(escores_gwas_attr[idx]))

    loci_gwas_ordenados = sorted(escore_por_locus.keys(), key=lambda k: escore_por_locus[k], reverse=True)
    top_gwas = set(loci_gwas_ordenados[:top_k])
    top_shap = set(loci_top_shap[:top_k])
    if not top_shap:
        return 0.0

    concordancia = len(top_shap.intersection(top_gwas)) / float(len(top_shap))
    return float(round(concordancia, 4))


def avaliar_remocao_locus_principal(
    modelo_xgb: xgb.XGBClassifier,
    X_teste: np.ndarray,
    valores_shap: np.ndarray,
) -> tuple[float, bool]:
    """
    Teste 4 da Seção 8.2 (Teste de remoção / ablação de locus principal):
    Remove (zera) o atributo de maior contribuição positiva SHAP em cada isolado predito como resistente
    e verifica se a probabilidade de resistência cai de forma coerente.
    """
    X_arr = np.asarray(X_teste, dtype=np.float64)
    shap_arr = np.asarray(valores_shap, dtype=np.float64)
    probs_orig = modelo_xgb.predict_proba(X_arr)[:, 1]

    mascara_res = probs_orig >= 0.50
    if not np.any(mascara_res):
        mascara_res = np.ones(len(probs_orig), dtype=bool)

    X_pert = X_arr[mascara_res].copy()
    shap_res = shap_arr[mascara_res]

    for i in range(X_pert.shape[0]):
        # Identifica o atributo presente com maior SHAP positivo e o remove
        s_presentes = np.where(X_pert[i] > 0, shap_res[i], -1e9)
        idx_top = int(np.argmax(s_presentes))
        X_pert[i, idx_top] = 0.0

    probs_pos_remocao = modelo_xgb.predict_proba(X_pert)[:, 1]
    quedas = probs_orig[mascara_res] - probs_pos_remocao
    queda_media = float(np.mean(quedas))
    return float(round(queda_media, 4)), bool(queda_media > 0.01)


def avaliar_explicabilidade_completa(
    farmaco: str,
    modelo_xgb: xgb.XGBClassifier,
    X_teste: np.ndarray,
    y_teste: np.ndarray,
    grupos_teste: list[str] | np.ndarray,
    nomes_atributos: list[str],
    rankings_por_dobra: list[dict[str, float]],
    top_k: int = 6,
) -> AvaliacaoExplicabilidadeQ3:
    """Executa todos os passos da Seção 8.1 e os 4 testes da Seção 8.2 para responder Q3."""
    mapa = construir_mapa_atributo_para_locus(nomes_atributos, farmaco=farmaco)
    valores_shap = calcular_valores_treeshap(modelo_xgb, X_teste)
    matriz_locus, ordem_loci, info_por_locus = agregar_shap_por_locus(
        valores_shap, nomes_atributos, mapa
    )

    prop_card, rec_card, cand_novos = avaliar_recuperacao_determinantes_card(
        matriz_locus, ordem_loci, info_por_locus, top_k=top_k
    )
    estabilidade = avaliar_estabilidade_explicacoes(rankings_por_dobra)

    importancia_media = np.mean(np.abs(matriz_locus), axis=0)
    ordem_top = np.argsort(importancia_media)[::-1][:top_k]
    loci_top_shap = [ordem_loci[int(i)] for i in ordem_top]

    conc_gwas = avaliar_concordancia_gwas_populacional(
        X=X_teste,
        y=y_teste,
        grupos_filogeneticos=grupos_teste,
        nomes_atributos=nomes_atributos,
        mapa_atributo_para_locus=mapa,
        loci_top_shap=loci_top_shap,
        top_k=top_k,
    )
    queda_remocao, coerente = avaliar_remocao_locus_principal(
        modelo_xgb=modelo_xgb,
        X_teste=X_teste,
        valores_shap=valores_shap,
    )

    principais_loci: list[dict[str, str | float]] = []
    for i in ordem_top:
        lid = ordem_loci[int(i)]
        info = info_por_locus[lid]
        principais_loci.append(
            {
                "locus": lid,
                "gene": info.gene,
                "alteracao": info.alteracao or "-",
                "tipo": info.tipo,
                "status": info.status,
                "fonte": info.fonte,
                "importancia_shap_media": round(float(importancia_media[int(i)]), 4),
                "mecanismo": info.mecanismo,
            }
        )

    return AvaliacaoExplicabilidadeQ3(
        farmaco=farmaco,
        proporcao_top_k_no_card=prop_card,
        determinantes_card_recuperados=rec_card,
        candidatos_novos_top_k=cand_novos,
        correlacao_estabilidade_entre_dobras=estabilidade,
        concordancia_gwas_populacional=conc_gwas,
        queda_media_probabilidade_remocao_locus=queda_remocao,
        remocao_coerente=coerente,
        principais_loci_globais=principais_loci,
    )


def gerar_relatorio_html_isolado(dados_resposta: dict[str, object]) -> str:
    """Gera relatório visual em HTML para um isolado consultado via GET /v1/predicoes/{job_id}/relatorio."""
    job_id = dados_resposta.get("job_id", "N/A")
    especie = dados_resposta.get("especie", "Escherichia coli")
    qc = dados_resposta.get("qc", {})
    predicoes = dados_resposta.get("predicoes", [])
    versao_modelo = dados_resposta.get("versao_modelo", "ecoli-amr-1.2.0")
    versoes_bancos = dados_resposta.get("versoes_bancos", {})
    aviso = dados_resposta.get("aviso", AVISO_USO_PESQUISA)

    blocos_pred = []
    if isinstance(predicoes, list):
        for p in predicoes:
            if not isinstance(p, dict):
                continue
            linhas_ev = []
            evs = p.get("evidencias", [])
            if isinstance(evs, list):
                for ev in evs:
                    if not isinstance(ev, dict):
                        continue
                    linhas_ev.append(
                        f"<tr>"
                        f"<td>{ev.get('tipo', '')}</td>"
                        f"<td><strong>{ev.get('gene', '')}</strong></td>"
                        f"<td>{ev.get('alteracao') or '-'}</td>"
                        f"<td>{ev.get('fonte', '')}</td>"
                        f"<td>{ev.get('contribuicao_shap', 0.0):+.3f}</td>"
                        f"<td><span class='badge {ev.get('status', '')}'>{ev.get('status', '')}</span></td>"
                        f"</tr>"
                    )
            tabela_ev = (
                "<table><thead><tr><th>Tipo</th><th>Gene/Locus</th><th>Alteração</th>"
                "<th>Fonte</th><th>Contribuição SHAP</th><th>Status</th></tr></thead>"
                f"<tbody>{''.join(linhas_ev) or '<tr><td colspan=6>Nenhuma evidência de resistência ativa detectada</td></tr>'}</tbody></table>"
            )
            blocos_pred.append(
                f"<div class='card'>"
                f"<h3>Fármaco: {p.get('farmaco', '')}</h3>"
                f"<p><strong>Classe Predita:</strong> {p.get('classe', '')} | "
                f"<strong>Probabilidade Calibrada:</strong> {p.get('probabilidade_resistencia', 0.0):.2f} | "
                f"<strong>Confiança:</strong> {p.get('confianca', '')}</p>"
                f"{tabela_ev}"
                f"</div>"
            )

    qc_dict = qc if isinstance(qc, dict) else {}
    return f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
  <meta charset="utf-8" />
  <title>Relatório PhyloResistance AI — Job {job_id}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, sans-serif; margin: 2rem auto; max-width: 960px; color: #1f2937; line-height: 1.5; }}
    .banner-aviso {{ background: #fef3c7; border-left: 5px solid #d97706; padding: 0.8rem 1rem; margin-bottom: 1.5rem; font-weight: 600; }}
    .card {{ border: 1px solid #e5e7eb; border-radius: 8px; padding: 1rem 1.25rem; margin-bottom: 1rem; background: #ffffff; }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 0.75rem; font-size: 0.92rem; }}
    th, td {{ border-bottom: 1px solid #e5e7eb; text-align: left; padding: 0.5rem; }}
    th {{ background: #f9fafb; }}
    .badge {{ padding: 0.15rem 0.5rem; border-radius: 4px; font-size: 0.8rem; font-weight: 600; }}
    .conhecida {{ background: #dcfce7; color: #166534; }}
    .candidata {{ background: #dbeafe; color: #1e40af; }}
    footer {{ margin-top: 2rem; font-size: 0.85rem; color: #6b7280; border-top: 1px solid #e5e7eb; padding-top: 1rem; }}
  </style>
</head>
<body>
  <h1>PhyloResistance AI — Relatório Genômico XAI</h1>
  <div class="banner-aviso">{aviso}</div>
  <div class="card">
    <h2>Identificação e Controle de Qualidade (QC)</h2>
    <p><strong>Job ID:</strong> {job_id} | <strong>Espécie:</strong> <em>{especie}</em></p>
    <p><strong>Tamanho:</strong> {qc_dict.get('tamanho_pb', 0)} pb | <strong>Contigs:</strong> {qc_dict.get('contigs', 0)} | <strong>N50:</strong> {qc_dict.get('n50', 0)} pb | <strong>QC Aprovado:</strong> {qc_dict.get('aprovado', False)}</p>
  </div>
  <h2>Predições e Evidências Genômicas (TreeSHAP Agregado por Locus)</h2>
  {''.join(blocos_pred)}
    <footer>
      <p><strong>Versão do Modelo:</strong> {versao_modelo} | <strong>Versões dos Bancos:</strong> {versoes_bancos}</p>
      <p>Nota científica (REG-10): Contribuições SHAP refletem associação estatística ao modelo preditivo, não causalidade biológica direta.</p>
    </footer>
  </body>
</html>"""


def _escapar_texto_pdf(texto: str) -> str:
    """Escapa parênteses e barras invertidas para strings literais PDF (Latin-1)."""
    limpo = texto.encode("latin-1", errors="replace").decode("latin-1")
    return limpo.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def gerar_relatorio_pdf_isolado(dados_resposta: dict[str, object]) -> bytes:
    """
    Gera relatório visual em formato binário PDF 1.4 (Seção 9.2) para um isolado
    consultado via GET /v1/predicoes/{job_id}/relatorio?formato=pdf.
    """
    job_id = str(dados_resposta.get("job_id", "N/A"))
    especie = str(dados_resposta.get("especie", "Escherichia coli"))
    qc = dados_resposta.get("qc", {})
    qc_dict = qc if isinstance(qc, dict) else {}
    predicoes = dados_resposta.get("predicoes", [])
    versao_modelo = str(dados_resposta.get("versao_modelo", "ecoli-amr-1.2.0"))
    aviso = str(dados_resposta.get("aviso", AVISO_USO_PESQUISA))

    linhas_pdf: list[str] = [
        f"PhyloResistance AI - Relatorio Genomico XAI (Job {job_id})",
        f"AVISO: {aviso}",
        f"Especie: {especie} | Modelo: {versao_modelo}",
        (
            f"QC: Tamanho={qc_dict.get('tamanho_pb', 0)} pb | "
            f"Contigs={qc_dict.get('contigs', 0)} | "
            f"N50={qc_dict.get('n50', 0)} pb | "
            f"Aprovado={qc_dict.get('aprovado', False)}"
        ),
        "-" * 75,
    ]

    if isinstance(predicoes, list):
        for p in predicoes:
            if not isinstance(p, dict):
                continue
            linhas_pdf.append(
                f"Farmaco: {p.get('farmaco', '')} | "
                f"Classe: {p.get('classe', '')} | "
                f"Prob. Calibrada: {float(p.get('probabilidade_resistencia', 0.0)):.2f} | "
                f"Confianca: {p.get('confianca', '')}"
            )
            evs = p.get("evidencias", [])
            if isinstance(evs, list) and evs:
                for ev in evs[:5]:
                    if not isinstance(ev, dict):
                        continue
                    linhas_pdf.append(
                        f"  * [{ev.get('tipo', '')}] {ev.get('gene', '')} "
                        f"({ev.get('alteracao') or '-'}) | "
                        f"SHAP={float(ev.get('contribuicao_shap', 0.0)):+.3f} | "
                        f"Status={ev.get('status', '')} ({ev.get('fonte', '')})"
                    )
            else:
                linhas_pdf.append("  * Nenhuma evidencia de resistencia ativa detectada.")
            linhas_pdf.append("-" * 75)

    linhas_pdf.append(
        "Nota (REG-10): Contribuicoes SHAP refletem associacao estatistica, nao causalidade direta."
    )

    comandos_stream = ["BT", "/F1 10 Tf", "50 790 Td", "14 TL"]
    for idx, linha in enumerate(linhas_pdf[:48]):
        esc = _escapar_texto_pdf(linha)
        if idx == 0:
            comandos_stream.append(f"({esc}) Tj")
        else:
            comandos_stream.append(f"T* ({esc}) Tj")
    comandos_stream.append("ET")
    conteudo_stream = "\n".join(comandos_stream).encode("latin-1")

    objetos: list[bytes] = [
        b"1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n",
        b"2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n",
        (
            b"3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] "
            b"/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>\nendobj\n"
        ),
        b"4 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>\nendobj\n",
        (
            f"5 0 obj\n<< /Length {len(conteudo_stream)} >>\nstream\n".encode("latin-1")
            + conteudo_stream
            + b"\nendstream\nendobj\n"
        ),
    ]

    buffer = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets: list[int] = []
    for obj in objetos:
        offsets.append(len(buffer))
        buffer.extend(obj)

    xref_pos = len(buffer)
    buffer.extend(f"xref\n0 {len(objetos) + 1}\n".encode("latin-1"))
    buffer.extend(b"0000000000 65535 f \n")
    for off in offsets:
        buffer.extend(f"{off:010d} 00000 n \n".encode("latin-1"))

    buffer.extend(
        (
            f"trailer\n<< /Size {len(objetos) + 1} /Root 1 0 R >>\n"
            f"startxref\n{xref_pos}\n%%EOF\n"
        ).encode("latin-1")
    )
    return bytes(buffer)


