"""Regras de curadoria padronizada de rótulos fenotípicos S/I/R e recálculo por MIC (Seção 3.3, REG-06)."""

from __future__ import annotations

from collections import defaultdict
from typing import Literal

import pandas as pd

from phyloresistance.config import (
    ESPECIE_PILOTO,
    FARMACOS_PILOTO,
    PONTOS_CORTE_ECOLI,
    PontoCorteMIC,
)
from phyloresistance.dados.contratos import (
    FunilExclusao,
    RegistroBrutoBVBRC,
    RegistroCurado,
)

PRIORIDADE_METODO: dict[str, int] = {
    "Broth microdilution": 1,
    "Agar dilution": 2,
    "Etest": 3,
    "Disk diffusion": 4,
}


def normalizar_categoria_fenotipo(valor: str | None) -> Literal["R", "S", "I"] | None:
    """Converte strings brutas de fenótipo em categorias padronizadas 'R', 'S', 'I' ou None."""
    if valor is None:
        return None
    limpo = valor.strip().lower()
    if limpo in {"r", "resistant", "resistente", "non-susceptible"}:
        return "R"
    if limpo in {"s", "susceptible", "sensivel", "sensível", "susceptible-dose dependent"}:
        return "S"
    if limpo in {"i", "intermediate", "intermediario", "intermediário"}:
        return "I"
    return None


def determinar_categoria_registro(
    reg: RegistroBrutoBVBRC,
    pontos_corte: dict[str, PontoCorteMIC] | None = None,
) -> tuple[Literal["R", "S", "I"] | None, Literal["recalculado_mic", "fenotipo_declarado"]]:
    """
    Preferir MIC quando existir (Seção 3.3): recalcula a categoria a partir do MIC
    usando um único ponto de corte documentado. Caso contrário, usa o fenótipo bruto.
    """
    tabela_corte = pontos_corte if pontos_corte is not None else PONTOS_CORTE_ECOLI
    if reg.mic_mg_l is not None and reg.mic_mg_l > 0 and reg.farmaco in tabela_corte:
        corte = tabela_corte[reg.farmaco]
        cat_mic = corte.classificar(reg.mic_mg_l)
        return cat_mic, "recalculado_mic"  # type: ignore[return-value]

    cat_bruta = normalizar_categoria_fenotipo(reg.fenotipo_bruto)
    return cat_bruta, "fenotipo_declarado"


def curar_rotulos(
    registros: list[RegistroBrutoBVBRC],
    especie_alvo: str = ESPECIE_PILOTO,
    farmacos_alvo: tuple[str, ...] = FARMACOS_PILOTO,
    modo_intermediario: Literal["excluir", "agrupar_s", "agrupar_r"] = "excluir",
    pontos_corte: dict[str, PontoCorteMIC] | None = None,
) -> tuple[list[RegistroCurado], FunilExclusao]:
    """
    Executa a curadoria padronizada de rótulos (Seção 3.3 e REG-06):
    1. Filtra espécie e fármacos do escopo.
    2. Recalcula categoria via MIC quando disponível.
    3. Trata 'Intermediário (I)' explicitamente (exclui no treino principal ou agrupa em análise de sensibilidade).
    4. Resolve duplicatas/conflitos por par (genome_id, farmaco).
    5. Registra todo o funil de exclusão.
    """
    funil = FunilExclusao(
        total_registros_brutos=len(registros),
        modo_intermediario=modo_intermediario,
        detalhes_por_farmaco={
            f: {"brutos": 0, "intermediarios": 0, "conflitos": 0, "aprovados_R": 0, "aprovados_S": 0}
            for f in farmacos_alvo
        },
    )

    agrupados: dict[tuple[str, str], list[tuple[RegistroBrutoBVBRC, Literal["R", "S"], Literal["recalculado_mic", "fenotipo_declarado"]]]] = defaultdict(list)

    for reg in registros:
        if reg.especie.strip() != especie_alvo:
            funil.removidos_especie_fora_escopo += 1
            continue

        if reg.farmaco not in farmacos_alvo:
            funil.removidos_farmaco_fora_escopo += 1
            continue

        funil.detalhes_por_farmaco[reg.farmaco]["brutos"] += 1

        cat, origem = determinar_categoria_registro(reg, pontos_corte=pontos_corte)
        if cat is None:
            funil.removidos_fenotipo_invalido += 1
            continue

        if cat == "I":
            if modo_intermediario == "excluir":
                funil.removidos_intermediario += 1
                funil.detalhes_por_farmaco[reg.farmaco]["intermediarios"] += 1
                continue
            elif modo_intermediario == "agrupar_s":
                cat_binaria: Literal["R", "S"] = "S"
            else:
                cat_binaria = "R"
        else:
            cat_binaria = cat

        agrupados[(reg.genome_id, reg.farmaco)].append((reg, cat_binaria, origem))

    aprovados: list[RegistroCurado] = []

    for (genome_id, farmaco), leituras in sorted(agrupados.items()):
        categorias_distintas = {item[1] for item in leituras}
        if len(categorias_distintas) > 1:
            # Conflito irreconciliável (ex.: uma leitura R e outra S para o mesmo isolado e fármaco)
            funil.removidos_conflito_rotulo += len(leituras)
            funil.detalhes_por_farmaco[farmaco]["conflitos"] += len(leituras)
            continue

        if len(leituras) > 1:
            funil.removidos_duplicata_concordante += len(leituras) - 1

        # Ordena por prioridade do método laboratorial (Broth microdilution primeiro) e presença de MIC
        leituras_ordenadas = sorted(
            leituras,
            key=lambda item: (
                0 if item[0].mic_mg_l is not None else 1,
                PRIORIDADE_METODO.get(item[0].metodo_laboratorio, 99),
            ),
        )
        melhor_reg, cat_final, origem_final = leituras_ordenadas[0]
        y_bin: Literal[0, 1] = 1 if cat_final == "R" else 0

        corte_ref = (pontos_corte or PONTOS_CORTE_ECOLI).get(farmaco)
        padrao_final = corte_ref.padrao if origem_final == "recalculado_mic" and corte_ref else melhor_reg.padrao
        versao_final = corte_ref.versao if origem_final == "recalculado_mic" and corte_ref else melhor_reg.versao_ponto_corte

        aprovados.append(
            RegistroCurado(
                genome_id=genome_id,
                especie=melhor_reg.especie,
                farmaco=farmaco,
                fenotipo=cat_final,
                y=y_bin,
                origem_rotulo=origem_final,
                mic_mg_l=melhor_reg.mic_mg_l,
                metodo_laboratorio=melhor_reg.metodo_laboratorio,
                padrao=padrao_final,
                versao_ponto_corte=versao_final,
                pais=melhor_reg.pais,
                ano=melhor_reg.ano,
                st=melhor_reg.st,
            )
        )
        if cat_final == "R":
            funil.detalhes_por_farmaco[farmaco]["aprovados_R"] += 1
        else:
            funil.detalhes_por_farmaco[farmaco]["aprovados_S"] += 1

    funil.total_aprovados_curadoria = len(aprovados)
    funil.total_aprovados_pos_qc_dedup = len(aprovados)
    return aprovados, funil


def converter_curados_para_dataframe(curados: list[RegistroCurado]) -> pd.DataFrame:
    """Converte lista de registros curados em DataFrame padronizado."""
    if not curados:
        return pd.DataFrame(
            columns=[
                "genome_id",
                "especie",
                "farmaco",
                "fenotipo",
                "y",
                "origem_rotulo",
                "mic_mg_l",
                "metodo_laboratorio",
                "padrao",
                "versao_ponto_corte",
                "pais",
                "ano",
                "st",
            ]
        )
    return pd.DataFrame([c.model_dump() for c in curados])

