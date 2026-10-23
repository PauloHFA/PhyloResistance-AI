# Model Card — PhyloResistance AI (`ecoli-amr-1.2.0`)

> **Regras de Rigor REG-14 e REG-15**: Documentação transparente de domínio, métricas, calibração, versões de bancos e limitações de uso exclusivo em pesquisa.

---

## 1. Identificação do Modelo

- **Identificador da Suite**: `ecoli-amr-1.2.0`
- **Modelos por Fármaco**:
  - `ecoli-cip-1.2.0` (*Escherichia coli* — Ciprofloxacina)
  - `ecoli-gen-1.2.0` (*Escherichia coli* — Gentamicina)
  - `ecoli-cxm-1.2.0` (*Escherichia coli* — Cefuroxima)
- **Arquitetura**: Ensemble de Árvores (XGBoost sobre atributos F1 + F2 + F3) com calibração probabilística isotônica/Platt, camada de abstenção conformal e explicabilidade TreeSHAP agregada por locus.
- **Versões dos Bancos de Referência**:
  - `AMRFinderPlus`: `3.12.8` (DB `2024-05-02.2`)
  - `CARD`: `3.2.9` (RGI `6.0.3`)
  - `RefSeq E. coli K-12 MG1655`: `NC_000913.3`

---

## 2. Uso Pretendido e Restrições (REG-15)

- **Uso Pretendido**: Pesquisa científica em genômica bacteriana, vigilância epidemiológica retrospectiva, geração de hipóteses sobre determinantes de resistência e demonstração metodológica de Machine Learning Explicável (XAI) livre de vazamento filogenético.
- **Uso Fora de Escopo (Proibido)**:
  - **Não substitui o antibiograma laboratorial (AST)** nem deve orientar prescrição clínica ou decisão terapêutica em pacientes.
  - Não se aplica a espécies diferentes de *Escherichia coli*, nem a metagenomas ou leituras brutas não montadas sem o pipeline prévio de montagem e QC.

---

## 3. Dados de Treinamento e Validação

- **Fonte primária**: Coorte curada no padrão BV-BRC / NCBI de genomas montados de *Escherichia coli* com metadados de MIC, método laboratorial, ano de coleta, país e Sequence Type (ST).
- **Estratégia de validação**:
  - Validação principal por **Grupo Filogenético (`StratifiedGroupKFold`, 5 dobras)** baseada em clusters de distância Mash/k-mers.
  - Validações adicionais: *Leave-One-Lineage-Out* (LOLO), Temporal e Geográfica.
  - Comparação obrigatória contra o **Baseline B1 (regras CARD/AMRFinderPlus)** com intervalos de confiança de 95% via bootstrap em blocos (clusters).

---

## 4. Explicabilidade e Interpretação Biológica (REG-10)

- Cada predição é acompanhada por valores **TreeSHAP** agregados por locus genômico (gene AMR, mutação pontual ou locus de unitig/k-mer mapeado contra a referência).
- Cada evidência recebe a classificação:
  - `conhecida`: determinante catalogado no CARD / AMRFinderPlus para a classe do fármaco.
  - `candidata`: variante genômica ou k-mer mapeado em locus não catalogado, reportado estritamente como **associação estatística / hipótese para validação experimental** (nunca como relação causal direta).

---

## 5. Calibração, Incerteza e Abstenção (REG-11)

- As probabilidades brutas do modelo passam por calibração em conjunto separado.
- Quando a probabilidade calibrada cai na faixa de incerteza conformal próxima ao limiar de decisão, o serviço emite `classe: "indeterminado"` e `confianca: "baixa"`, abstendo-se de uma classificação binária forçada para reduzir erros graves (*Very Major Error* — VME e *Major Error* — ME).
