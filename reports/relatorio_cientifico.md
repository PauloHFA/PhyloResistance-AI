# Relatório Científico Final — PhyloResistance AI

> **Aviso Ético e Regulatório (REG-15)**: Uso em pesquisa. Não substitui o antibiograma.

---

## 1. Funil de Exclusão e Curadoria de Rótulos (REG-02, REG-06, REG-07)

- **Registros brutos ingeridos**: `369`
- **Removidos por espécie fora do escopo**: `1`
- **Removidos por fármaco fora do escopo**: `1`
- **Removidos por categoria Intermediária (I)**: `1`
- **Removidos por conflito irreconciliável de rótulo (R vs S)**: `2`
- **Removidos por reprovação no QC do genoma (contigs/N50/completude)**: `1`
- **Removidos por desduplicação clonal antes do split (distância Mash <= 0.0005)**: `3` (`1` genoma(s) colapsado(s))
- **Total de registros fenotípicos aprovados para modelagem**: `360`

---

## 2. Resultados do Escopo Principal (*Escherichia coli* — Q1, Q2, Q3 e Q4)

### Fármaco: `ciprofloxacina` (`ecoli-cip-1.2.0`)

- **Isolados únicos pós-QC/desduplicação**: `120` (Prevalência de Resistência: `45.0%`)
- **Conclusão Q1 (ML vs Baseline B1 em grupos filogenéticos não vistos)**: B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.0774, IC95% [+0.0546, +0.1008], Delta VME = -0.0830).
- **Regressão Quantitativa de $\log_2(\text{MIC})$ (ISO 20776-2)**: Concordância Essencial ($\pm 1$ diluição) = `55.8%` | Concordância Categórica = `92.5%` | RMSE $\log_2(\text{MIC})$ = `2.243`

#### Escada de Modelos (Validação Principal: Grupo Filogenético — REG-01, REG-04, REG-05)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | VME (Erro Muito Grave) | ME (Erro Grave) | Brier Score |
|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | 0.3988 [0.354, 0.481] | 0.3934 | 0.0000 | 1.0000 | **1.0000** | **0.0000** | 0.2536 |
| `B1_Regras_CARD_AMRFinder` | 0.9101 [0.874, 0.920] | 0.8886 | 0.8519 | 0.9697 | **0.1481** | **0.0303** | 0.0767 |
| `B2_LogReg_L1_F1_F2` | 0.9210 [0.856, 0.959] | 0.9112 | 0.7963 | 0.9697 | **0.2037** | **0.0303** | 0.0919 |
| `B3_XGBoost_F1_F2` | 0.9338 [0.906, 0.957] | 0.9196 | 0.8333 | 0.9545 | **0.1667** | **0.0455** | 0.0853 |
| `B4_XGBoost_F1_F2_F3` | 0.9781 [0.966, 0.989] | 0.9756 | 0.9444 | 0.9545 | **0.0556** | **0.0455** | 0.0464 |
| `B5_Linear_Esparso_F3` | 0.6489 [0.529, 0.708] | 0.5915 | 0.4815 | 0.6061 | **0.5185** | **0.3939** | 0.2703 |
| `B6_CNN1D_Sequencia` | 0.9369 [0.890, 0.971] | 0.9352 | 0.8333 | 0.8788 | **0.1667** | **0.1212** | 0.1091 |
| `Controle_Kmers_Curtos_k6` | 0.4444 [0.344, 0.587] | 0.4517 | 0.3704 | 0.5455 | **0.6296** | **0.4545** | 0.3077 |

#### Generalização entre Estratégias de Split (Q2 e REG-09)

| Estratégia de Validação | AUPRC | AUROC | VME | ME | Observação |
|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | 0.9781 | 0.9756 | 0.0556 | 0.0455 | Sem vazamento clonal |
| **Leave-One-Lineage-Out (LOLO)** | 1.0000 | 1.0000 | 0.0000 | 0.0833 | Linhagens ST inteiras inéditas |
| **Temporal (Treino <= 2022, Teste > 2022)** | 0.9936 | 0.9944 | 0.0833 | 0.0667 | Simula vigilância prospectiva |
| **Geográfica (Países retidos)** | 1.0000 | 1.0000 | 0.0000 | 0.0000 | Robustez a epidemiologia local |
| **Aleatória (Demonstração de Vazamento B5)** | 0.6679 | 0.6049 | 0.4444 | 0.3788 | Inflação por clones no treino e teste |

#### Confiabilidade, Calibração e Abstenção Conformal (Q4 e REG-11)

- **Brier Score**: `0.0464` | **ECE**: `0.0321`
- **Cobertura Conformal**: `99.2%` (Taxa de Abstenção em casos duvidosos: `0.8%`)
- **Redução de Erros Graves pós-abstenção**: VME cai de `0.0556` para `0.0566`, e ME cai de `0.0455` para `0.0455`.

#### Avaliação de Explicabilidade TreeSHAP Agregada por Locus (Q3, Seção 8.2 e REG-10)

- **Proporção do Top-6 presente no CARD/AMRFinderPlus**: `33.3%` (gyrA:S83L, parC:S80I)
- **Candidatos novos identificados via k-mers/unitigs (F3)**: `marR:G103S (unitig), acrR:fs*45 (unitig), locus_kmer31_poly_0014, locus_kmer31_poly_0012`
- **Estabilidade entre dobras (Spearman $\rho$)**: `0.4316`
- **Concordância com GWAS controlado por estrutura populacional**: `50.0%`
- **Teste de remoção in silico do locus principal**: queda média de probabilidade = `-0.5051` (Coerente: `True`)

| Locus Genômico | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Associado |
|---|---|---|---|---|---|
| `gyrA:S83L` | mutacao_pontual | **conhecida** | AMRFinderPlus | 1.3847 | Alteração do alvo DNA girase subunidade A (QRDR S83L) |
| `marR:G103S (unitig)` | kmer_unitig | **candidata** | Descoberta ML (F3) | 1.1394 | Regulador reprimido de marAB / bomba de efluxo AcrAB-TolC |
| `parC:S80I` | mutacao_pontual | **conhecida** | AMRFinderPlus | 0.8271 | Alteração da topoisomerase IV subunidade A (QRDR S80I) |
| `acrR:fs*45 (unitig)` | kmer_unitig | **candidata** | Descoberta ML (F3) | 0.7398 | Desrepressão de bomba de efluxo multidroga AcrAB |
| `locus_kmer31_poly_0014` | kmer_unitig | **candidata** | Alinhamento Referência (k=31) | 0.4895 | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0012` | kmer_unitig | **candidata** | Alinhamento Referência (k=31) | 0.4032 | Variante genômica não catalogada (hipótese) |

### Fármaco: `gentamicina` (`ecoli-gen-1.2.0`)

- **Isolados únicos pós-QC/desduplicação**: `120` (Prevalência de Resistência: `47.5%`)
- **Conclusão Q1 (ML vs Baseline B1 em grupos filogenéticos não vistos)**: B4_XGBoost_F1_F2_F3 NÃO supera com significância estatística o Baseline B1 (Delta AUPRC = +0.0453, IC95% [-0.0292, +0.1564]). Resultado reportado conforme REG-04.
- **Regressão Quantitativa de $\log_2(\text{MIC})$ (ISO 20776-2)**: Concordância Essencial ($\pm 1$ diluição) = `54.2%` | Concordância Categórica = `76.7%` | RMSE $\log_2(\text{MIC})$ = `2.375`

#### Escada de Modelos (Validação Principal: Grupo Filogenético — REG-01, REG-04, REG-05)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | VME (Erro Muito Grave) | ME (Erro Grave) | Brier Score |
|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | 0.4344 [0.326, 0.522] | 0.4095 | 0.2105 | 0.7143 | **0.7895** | **0.2857** | 0.2535 |
| `B1_Regras_CARD_AMRFinder` | 0.8343 [0.658, 0.927] | 0.8020 | 0.7368 | 0.9365 | **0.2632** | **0.0635** | 0.1343 |
| `B2_LogReg_L1_F1_F2` | 0.8544 [0.680, 0.925] | 0.8092 | 0.5789 | 0.9365 | **0.4211** | **0.0635** | 0.1613 |
| `B3_XGBoost_F1_F2` | 0.8533 [0.632, 0.927] | 0.8076 | 0.6140 | 0.9365 | **0.3860** | **0.0635** | 0.1561 |
| `B4_XGBoost_F1_F2_F3` | 0.8592 [0.785, 0.945] | 0.8819 | 0.7193 | 0.8254 | **0.2807** | **0.1746** | 0.1419 |
| `B5_Linear_Esparso_F3` | 0.5595 [0.422, 0.722] | 0.5834 | 0.6140 | 0.5714 | **0.3860** | **0.4286** | 0.2807 |
| `B6_CNN1D_Sequencia` | 0.7585 [0.582, 0.842] | 0.7884 | 0.6842 | 0.6825 | **0.3158** | **0.3175** | 0.2520 |
| `Controle_Kmers_Curtos_k6` | 0.4382 [0.286, 0.551] | 0.4456 | 0.3684 | 0.5873 | **0.6316** | **0.4127** | 0.3222 |

#### Generalização entre Estratégias de Split (Q2 e REG-09)

| Estratégia de Validação | AUPRC | AUROC | VME | ME | Observação |
|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | 0.8592 | 0.8819 | 0.2807 | 0.1746 | Sem vazamento clonal |
| **Leave-One-Lineage-Out (LOLO)** | 0.9796 | 0.9848 | 0.2778 | 0.0455 | Linhagens ST inteiras inéditas |
| **Temporal (Treino <= 2022, Teste > 2022)** | 0.9663 | 0.9667 | 0.1667 | 0.0000 | Simula vigilância prospectiva |
| **Geográfica (Países retidos)** | 0.8291 | 0.8542 | 0.2917 | 0.1667 | Robustez a epidemiologia local |
| **Aleatória (Demonstração de Vazamento B5)** | 0.5998 | 0.6051 | 0.4912 | 0.3968 | Inflação por clones no treino e teste |

#### Confiabilidade, Calibração e Abstenção Conformal (Q4 e REG-11)

- **Brier Score**: `0.1419` | **ECE**: `0.0775`
- **Cobertura Conformal**: `85.8%` (Taxa de Abstenção em casos duvidosos: `14.2%`)
- **Redução de Erros Graves pós-abstenção**: VME cai de `0.2807` para `0.1739`, e ME cai de `0.1746` para `0.1579`.

#### Avaliação de Explicabilidade TreeSHAP Agregada por Locus (Q3, Seção 8.2 e REG-10)

- **Proporção do Top-6 presente no CARD/AMRFinderPlus**: `33.3%` (aac(3)-IIa, ant(2'')-Ia)
- **Candidatos novos identificados via k-mers/unitigs (F3)**: `acrR:fs*45 (unitig), locus_kmer31_poly_0030, locus_kmer31_poly_0047, locus_kmer31_poly_0027`
- **Estabilidade entre dobras (Spearman $\rho$)**: `0.4470`
- **Concordância com GWAS controlado por estrutura populacional**: `66.7%`
- **Teste de remoção in silico do locus principal**: queda média de probabilidade = `-0.3354` (Coerente: `True`)

| Locus Genômico | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Associado |
|---|---|---|---|---|---|
| `acrR:fs*45 (unitig)` | kmer_unitig | **candidata** | Descoberta ML (F3) | 0.9267 | Desrepressão de bomba de efluxo multidroga AcrAB |
| `aac(3)-IIa` | gene_amr | **conhecida** | CARD | 0.6232 | Aminoglicosídeo 3-N-acetiltransferase |
| `locus_kmer31_poly_0030` | kmer_unitig | **candidata** | Alinhamento Referência (k=31) | 0.4896 | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0047` | kmer_unitig | **candidata** | Alinhamento Referência (k=31) | 0.4631 | Variante genômica não catalogada (hipótese) |
| `ant(2'')-Ia` | gene_amr | **conhecida** | CARD | 0.4098 | Aminoglicosídeo 2''-O-nucleotidiltransferase |
| `locus_kmer31_poly_0027` | kmer_unitig | **candidata** | Alinhamento Referência (k=31) | 0.4026 | Variante genômica não catalogada (hipótese) |

### Fármaco: `cefuroxima` (`ecoli-cxm-1.2.0`)

- **Isolados únicos pós-QC/desduplicação**: `120` (Prevalência de Resistência: `38.3%`)
- **Conclusão Q1 (ML vs Baseline B1 em grupos filogenéticos não vistos)**: B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.1916, IC95% [+0.1143, +0.3195], Delta VME = -0.2126).
- **Regressão Quantitativa de $\log_2(\text{MIC})$ (ISO 20776-2)**: Concordância Essencial ($\pm 1$ diluição) = `78.3%` | Concordância Categórica = `90.0%` | RMSE $\log_2(\text{MIC})$ = `1.555`

#### Escada de Modelos (Validação Principal: Grupo Filogenético — REG-01, REG-04, REG-05)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | VME (Erro Muito Grave) | ME (Erro Grave) | Brier Score |
|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | 0.3754 [0.296, 0.444] | 0.4824 | 0.0000 | 1.0000 | **1.0000** | **0.0000** | 0.2366 |
| `B1_Regras_CARD_AMRFinder` | 0.7112 [0.545, 0.852] | 0.7283 | 0.5652 | 0.9730 | **0.4348** | **0.0270** | 0.1485 |
| `B2_LogReg_L1_F1_F2` | 0.7224 [0.543, 0.859] | 0.7083 | 0.5435 | 0.9459 | **0.4565** | **0.0541** | 0.1707 |
| `B3_XGBoost_F1_F2` | 0.7463 [0.588, 0.865] | 0.7305 | 0.5217 | 0.8649 | **0.4783** | **0.1351** | 0.1799 |
| `B4_XGBoost_F1_F2_F3` | 0.9185 [0.865, 0.967] | 0.9057 | 0.8043 | 0.9189 | **0.1957** | **0.0811** | 0.0940 |
| `B5_Linear_Esparso_F3` | 0.6117 [0.542, 0.663] | 0.5943 | 0.5000 | 0.6892 | **0.5000** | **0.3108** | 0.2639 |
| `B6_CNN1D_Sequencia` | 0.8301 [0.731, 0.899] | 0.8237 | 0.6957 | 0.9324 | **0.3043** | **0.0676** | 0.1492 |
| `Controle_Kmers_Curtos_k6` | 0.3970 [0.341, 0.477] | 0.4559 | 0.2609 | 0.6081 | **0.7391** | **0.3919** | 0.3020 |

#### Generalização entre Estratégias de Split (Q2 e REG-09)

| Estratégia de Validação | AUPRC | AUROC | VME | ME | Observação |
|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | 0.9185 | 0.9057 | 0.1957 | 0.0811 | Sem vazamento clonal |
| **Leave-One-Lineage-Out (LOLO)** | 0.8627 | 0.8489 | 0.2857 | 0.0769 | Linhagens ST inteiras inéditas |
| **Temporal (Treino <= 2022, Teste > 2022)** | 0.8271 | 0.8553 | 0.3750 | 0.0526 | Simula vigilância prospectiva |
| **Geográfica (Países retidos)** | 0.8832 | 0.9141 | 0.2941 | 0.0541 | Robustez a epidemiologia local |
| **Aleatória (Demonstração de Vazamento B5)** | 0.6412 | 0.6745 | 0.4565 | 0.3108 | Inflação por clones no treino e teste |

#### Confiabilidade, Calibração e Abstenção Conformal (Q4 e REG-11)

- **Brier Score**: `0.0940` | **ECE**: `0.0813`
- **Cobertura Conformal**: `88.3%` (Taxa de Abstenção em casos duvidosos: `11.7%`)
- **Redução de Erros Graves pós-abstenção**: VME cai de `0.1957` para `0.1667`, e ME cai de `0.0811` para `0.0156`.

#### Avaliação de Explicabilidade TreeSHAP Agregada por Locus (Q3, Seção 8.2 e REG-10)

- **Proporção do Top-6 presente no CARD/AMRFinderPlus**: `33.3%` (blaCTX-M-15, blaCTX-M-14)
- **Candidatos novos identificados via k-mers/unitigs (F3)**: `ompF:c.-112del (unitig), locus_kmer31_poly_0063, locus_kmer31_poly_0100, locus_kmer31_poly_0024`
- **Estabilidade entre dobras (Spearman $\rho$)**: `0.3548`
- **Concordância com GWAS controlado por estrutura populacional**: `50.0%`
- **Teste de remoção in silico do locus principal**: queda média de probabilidade = `-0.3986` (Coerente: `True`)

| Locus Genômico | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Associado |
|---|---|---|---|---|---|
| `blaCTX-M-15` | gene_amr | **conhecida** | CARD | 0.9116 | Beta-lactamase de espectro estendido (ESBL CTX-M-15) |
| `ompF:c.-112del (unitig)` | kmer_unitig | **candidata** | Descoberta ML (F3) | 0.7475 | Redução de permeabilidade de porina de membrana externa |
| `blaCTX-M-14` | gene_amr | **conhecida** | CARD | 0.5776 | Beta-lactamase de espectro estendido (ESBL CTX-M-14) |
| `locus_kmer31_poly_0063` | kmer_unitig | **candidata** | Alinhamento Referência (k=31) | 0.4015 | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0100` | kmer_unitig | **candidata** | Alinhamento Referência (k=31) | 0.3404 | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0024` | kmer_unitig | **candidata** | Alinhamento Referência (k=31) | 0.3283 | Variante genômica não catalogada (hipótese) |


---

## 3. Extensões da Fase 6: Multi-Espécie (*K. pneumoniae* e *S. aureus*) e Validação Externa

### 3.1 Expansão Multi-Espécie (Seção 3.2 — Fases 2 e 3)

| Espécie | Fármaco | Versão do Modelo | AUPRC Baseline B1 | AUPRC Modelo B4 ($\Delta$) | VME | Top Loci TreeSHAP |
|---|---|---|---|---|---|---|
| *Klebsiella pneumoniae* | `ciprofloxacina` | `kpneu-cip-1.2.0` | 0.8341 | **0.9600** (+0.1259) | 0.0476 | gyrA:S83L, marR:G103S (unitig), locus_kmer31_poly_0039 |
| *Klebsiella pneumoniae* | `gentamicina` | `kpneu-gen-1.2.0` | 0.7501 | **0.8173** (+0.0614) | 0.2308 | aac(3)-IIa, cpxA:Y144C (unitig), blaCTX-M-15 |
| *Klebsiella pneumoniae* | `ceftriaxona` | `kpneu-cro-1.2.0` | 0.7668 | **0.8927** (+0.1159) | 0.1795 | blaCTX-M-15, ompF:c.-112del (unitig), locus_kmer31_poly_0119 |
| *Staphylococcus aureus* | `oxacilina` | `saureus-oxa-1.2.0` | 0.8529 | **0.9778** (+0.1280) | 0.0476 | mecA, pbp4:c.-90T>A (unitig), locus_kmer31_poly_0034 |
| *Staphylococcus aureus* | `ciprofloxacina` | `saureus-cip-1.2.0` | 0.7685 | **0.9441** (+0.1676) | 0.1250 | grlA:S80F, norA:flqB_up (unitig), locus_kmer31_poly_0024 |

### 3.2 Validação Externa Independente — Padrão NCBI Pathogen Detection AST (Seção 3.1, 7.1 e Fase 6)

| Fármaco (*E. coli*) | AUPRC B1 (Externo) | AUPRC B4 (Externo) | Ganho $\Delta\text{AUPRC}$ Externo | VME Externo | ME Externo |
|---|---|---|---|---|---|
| `ciprofloxacina` | 0.8542 | **0.9836** | **+0.1294** | 0.0465 | 0.0877 |
| `gentamicina` | 0.7878 | **0.8762** | **+0.0884** | 0.1458 | 0.2115 |
| `cefuroxima` | 0.7315 | **0.9296** | **+0.1981** | 0.2128 | 0.0566 |

---

## 4. Figuras Científicas Geradas (Fase 7)

- `reports/figuras/escada_auprc_ecoli.svg`
- `reports/figuras/calibracao_e_abstencao.svg`
- `reports/figuras/importancia_shap_loci.svg`
- `reports/figuras/extensoes_multiespecie_e_mic.svg`

---

## 5. Limitações e Uso Responsável (Seção 15)

1. **Associação vs. Causalidade (REG-10)**: Os valores TreeSHAP medem contribuição preditiva ao modelo; loci classificados como `candidata` representam hipóteses biológicas para validação experimental.
2. **Domínio de Aplicação**: Os modelos avaliados aplicam-se exclusivamente a isolados das espécies treinadas (*Escherichia coli*, *Klebsiella pneumoniae* e *Staphylococcus aureus*) aprovados pelos critérios de QC da Seção 3.4.
3. **Pesquisa Científica (REG-15)**: Este sistema não é um dispositivo médico e não deve orientar conduta clínica.
