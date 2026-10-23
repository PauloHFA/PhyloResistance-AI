"""Gerador de figuras científicas vetoriais (SVG) para o Relatório Quarto e Portfólio (Fase 7)."""

from __future__ import annotations

from pathlib import Path


def _barra_horizontal_svg(
    titulo: str,
    subtitulo: str,
    itens: list[tuple[str, float, str]],
    max_val: float = 1.0,
) -> str:
    """Renderiza gráfico de barras horizontais limpo em SVG."""
    largura = 760
    altura_topo = 75
    altura_barra = 36
    altura = altura_topo + len(itens) * altura_barra + 30
    area_barra_x = 240
    area_barra_w = 430

    linhas_svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {largura} {altura}" width="100%" height="{altura}">',
        '<rect width="100%" height="100%" fill="#ffffff" rx="8" stroke="#e5e7eb" />',
        f'<text x="24" y="34" font-family="system-ui, sans-serif" font-size="16" font-weight="700" fill="#111827">{titulo}</text>',
        f'<text x="24" y="54" font-family="system-ui, sans-serif" font-size="12" fill="#6b7280">{subtitulo}</text>',
    ]

    for idx, (rotulo, valor, cor) in enumerate(itens):
        y = altura_topo + idx * altura_barra
        proporcao = max(0.0, min(1.0, valor / max(1e-6, max_val)))
        w = int(proporcao * area_barra_w)
        linhas_svg.extend(
            [
                f'<text x="24" y="{y + 18}" font-family="monospace" font-size="12" fill="#374151">{rotulo}</text>',
                f'<rect x="{area_barra_x}" y="{y + 4}" width="{area_barra_w}" height="20" fill="#f3f4f6" rx="4" />',
                f'<rect x="{area_barra_x}" y="{y + 4}" width="{w}" height="20" fill="{cor}" rx="4" />',
                f'<text x="{area_barra_x + w + 8}" y="{y + 18}" font-family="system-ui, sans-serif" font-size="12" font-weight="600" fill="#111827">{valor:.4f}</text>',
            ]
        )

    linhas_svg.append("</svg>")
    return "\n".join(linhas_svg)


def gerar_figuras_svg_relatorio(
    relatorio_consolidado: dict[str, object],
    diretorio_figuras: str | Path = "reports/figuras",
) -> list[str]:
    """
    Gera 4 figuras vetoriais SVG em `reports/figuras/` resumindo Q1, Q2, Q3, Q4 e as extensões da Fase 6.
    """
    dir_fig = Path(diretorio_figuras)
    dir_fig.mkdir(parents=True, exist_ok=True)
    caminhos_gerados: list[str] = []

    res_farmacos = relatorio_consolidado.get("resultados_por_farmaco", {})
    if not isinstance(res_farmacos, dict):
        return caminhos_gerados

    # 1) Figura 1: Escada de Modelos (AUPRC em Ciprofloxacina - E. coli)
    if "ciprofloxacina" in res_farmacos:
        escada_cip = res_farmacos["ciprofloxacina"]["escada_modelos_grupo_filogenetico"]
        itens_escada: list[tuple[str, float, str]] = []
        for nome_m, dados_m in escada_cip.items():
            auprc_val = float(dados_m["auprc"])
            cor = "#2563eb" if "B4" in nome_m else ("#16a34a" if "B1" in nome_m else ("#dc2626" if "Controle" in nome_m else "#64748b"))
            itens_escada.append((nome_m[:26], auprc_val, cor))
        svg_1 = _barra_horizontal_svg(
            titulo="Q1: Escada de Modelos — AUPRC em Grupos Filogenéticos Não Vistos (Ciprofloxacina)",
            subtitulo="Verde = Baseline B1 (CARD/AMRFinderPlus) | Azul = B4 (XGBoost F1+F2+F3) | Vermelho = Controle Negativo (k=6)",
            itens=itens_escada,
            max_val=1.0,
        )
        c1 = dir_fig / "escada_auprc_ecoli.svg"
        c1.write_text(svg_1, encoding="utf-8")
        caminhos_gerados.append(str(c1))

    # 2) Figura 2: Redução de VME (Very Major Error) via Predição Conformal / Abstenção (Q4)
    itens_vme: list[tuple[str, float, str]] = []
    for farm, dados_f in res_farmacos.items():
        m_b4 = dados_f["escada_modelos_grupo_filogenetico"]["B4_XGBoost_F1_F2_F3"]
        itens_vme.append((f"{farm[:12]} (VME bruto)", float(m_b4["vme"]), "#f59e0b"))
        itens_vme.append((f"{farm[:12]} (VME pós-abst)", float(m_b4["vme_pos_abstencao"]), "#10b981"))
    svg_2 = _barra_horizontal_svg(
        titulo="Q4: Redução de Erro Muito Grave (VME) por Abstenção Conformal",
        subtitulo="Comparação da taxa de falsos sensíveis antes (laranja) e depois (verde) de abster-se na zona de dúvida",
        itens=itens_vme,
        max_val=0.40,
    )
    c2 = dir_fig / "calibracao_e_abstencao.svg"
    c2.write_text(svg_2, encoding="utf-8")
    caminhos_gerados.append(str(c2))

    # 3) Figura 3: Importância TreeSHAP Agregada por Locus (Q3)
    if "ciprofloxacina" in res_farmacos:
        loci = res_farmacos["ciprofloxacina"]["explicabilidade_q3"]["principais_loci_globais"]
        itens_shap: list[tuple[str, float, str]] = []
        max_shap = max([float(loc["importancia_shap_media"]) for loc in loci] + [1.0])
        for loc in loci:
            cor_l = "#16a34a" if loc["status"] == "conhecida" else "#3b82f6"
            itens_shap.append((f"{loc['locus'][:24]} ({loc['status'][:4]})", float(loc["importancia_shap_media"]), cor_l))
        svg_3 = _barra_horizontal_svg(
            titulo="Q3: Top-6 Loci por Importância TreeSHAP Agregada (E. coli - Ciprofloxacina)",
            subtitulo="Verde = Determinante Conhecido (CARD/AMRFinderPlus) | Azul = Variante Candidata Descoberta em F3 (k=31)",
            itens=itens_shap,
            max_val=max_shap * 1.1,
        )
        c3 = dir_fig / "importancia_shap_loci.svg"
        c3.write_text(svg_3, encoding="utf-8")
        caminhos_gerados.append(str(c3))

    # 4) Figura 4: Extensões da Fase 6 (Multi-espécie, MIC ISO 20776-2 e Validação Externa)
    ext = relatorio_consolidado.get("extensoes_fase6", {})
    if isinstance(ext, dict):
        itens_ext: list[tuple[str, float, str]] = []
        for esp, d_esp in ext.get("multiespecie", {}).items():
            for farm_e, m_e in d_esp.items():
                itens_ext.append((f"{esp[:10]}:{farm_e[:10]} AUPRC", float(m_e["auprc_b4"]), "#6366f1"))
        for farm_mic, d_mic in ext.get("regressao_mic_iso20776", {}).items():
            itens_ext.append((f"MIC {farm_mic[:10]} Conc.Ess.", float(d_mic["concordancia_essencial_1_diluicao"]), "#0d9488"))
        if itens_ext:
            svg_4 = _barra_horizontal_svg(
                titulo="Fase 6: Desempenho Multi-Espécie (K. pneumoniae, S. aureus) e Concordância Essencial MIC",
                subtitulo="Roxo = AUPRC em grupos filogenéticos inéditos | Turquesa = Concordância Essencial ISO 20776-2 (+-1 diluição)",
                itens=itens_ext,
                max_val=1.0,
            )
            c4 = dir_fig / "extensoes_multiespecie_e_mic.svg"
            c4.write_text(svg_4, encoding="utf-8")
            caminhos_gerados.append(str(c4))

    return caminhos_gerados

