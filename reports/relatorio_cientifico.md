# Relatório Científico Final — PhyloResistance AI

> **Aviso Ético e Regulatório (REG-15)**: Uso em pesquisa. Não substitui o antibiograma.

---

## 1. Funil de Exclusão e Curadoria de Rótulos (REG-02, REG-06, REG-07)

- **Registros brutos ingeridos**: `1089`
- **Removidos por espécie fora do escopo**: `1`
- **Removidos por fármaco fora do escopo**: `1`
- **Removidos por categoria Intermediária (I)**: `1`
- **Removidos por conflito irreconciliável de rótulo (R vs S)**: `2`
- **Removidos por reprovação no QC do genoma (contigs/N50/completude)**: `1`
- **Removidos por desduplicação clonal antes do split (distância Mash <= 0.0005)**: `3` (`1` genoma(s) colapsado(s))
- **Total de registros fenotípicos aprovados para modelagem**: `1080`

---

## 2. Resultados do Escopo Principal (*Escherichia coli* — Q1, Q2, Q3 e Q4)

### Fármaco: `ciprofloxacina` (`ecoli-cip-1.2.0`)

- **Isolados únicos pós-QC/desduplicação**: `360` (Prevalência de Resistência: `48.3%`)
- **Conclusão Q1 (ML vs Baseline B1 em grupos filogenéticos não vistos)**: B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.0875, IC95% [+0.0402, +0.1396], Delta VME = -0.1782).
- **Regressão Quantitativa de $\log_2(\text{MIC})$ (ISO 20776-2)**: Concordância Essencial ($\pm 1$ diluição) = `58.6%` | Concordância Categórica = `88.6%` | RMSE $\log_2(\text{MIC})$ = `2.240`

#### Escada de Modelos (Validação Principal: Grupo Filogenético — REG-01, REG-04, REG-05)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | VME (Erro Muito Grave) | ME (Erro Grave) | Brier Score |
|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | 0.4486 [0.399, 0.508] | 0.4323 | 0.0000 | 1.0000 | **1.0000** | **0.0000** | 0.2525 |
| `B1_Regras_CARD_AMRFinder` | 0.8763 [0.820, 0.922] | 0.8499 | 0.7471 | 0.9677 | **0.2529** | **0.0323** | 0.1143 |
| `B2_LogReg_L1_F1_F2` | 0.9005 [0.852, 0.943] | 0.8691 | 0.7356 | 0.9731 | **0.2644** | **0.0269** | 0.1118 |
| `B3_XGBoost_F1_F2` | 0.8987 [0.844, 0.935] | 0.8665 | 0.7356 | 0.9570 | **0.2644** | **0.0430** | 0.1163 |
| `B4_XGBoost_F1_F2_F3` | 0.9621 [0.946, 0.978] | 0.9644 | 0.9253 | 0.9247 | **0.0747** | **0.0753** | 0.0553 |
| `B5_Linear_Esparso_F3` | 0.7602 [0.706, 0.803] | 0.7249 | 0.5747 | 0.8011 | **0.4253** | **0.1989** | 0.2126 |
| `B6_CNN1D_Sequencia` | 0.9556 [0.942, 0.970] | 0.9529 | 0.8793 | 0.9086 | **0.1207** | **0.0914** | 0.0872 |
| `Controle_Kmers_Curtos_k6` | 0.4862 [0.413, 0.558] | 0.4840 | 0.4483 | 0.4946 | **0.5517** | **0.5054** | 0.2760 |

#### Generalização entre Estratégias de Split (Q2 e REG-09)

| Estratégia de Validação | AUPRC | AUROC | VME | ME | Observação |
|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | 0.9621 | 0.9644 | 0.0747 | 0.0753 | Sem vazamento clonal |
| **Leave-One-Lineage-Out (LOLO)** | 0.9496 | 0.9530 | 0.0909 | 0.0923 | Linhagens ST inteiras inéditas |
| **Temporal (Treino <= 2022, Teste > 2022)** | 0.9921 | 0.9918 | 0.0833 | 0.0217 | Simula vigilância prospectiva |
| **Geográfica (Países retidos)** | 0.9569 | 0.9652 | 0.0513 | 0.0714 | Robustez a epidemiologia local |
| **Aleatória (Demonstração de Vazamento B5)** | 0.7597 | 0.7414 | 0.3736 | 0.3011 | Inflação por clones no treino e teste |

#### Confiabilidade, Calibração e Abstenção Conformal (Q4 e REG-11)

- **Brier Score**: `0.0553` | **ECE**: `0.0205`
- **Cobertura Conformal**: `94.4%` (Taxa de Abstenção em casos duvidosos: `5.6%`)
- **Redução de Erros Graves pós-abstenção**: VME cai de `0.0747` para `0.0370`, e ME cai de `0.0753` para `0.0506`.

#### Avaliação de Explicabilidade TreeSHAP Agregada por Locus (Q3, Seção 8.2 e REG-10)

- **Proporção do Top-6 presente no CARD/AMRFinderPlus**: `50.0%` (gyrA:S83L, gyrA:D87N, parC:S80I)
- **Candidatos novos identificados via k-mers/unitigs (F3)**: `acrR:fs*45 (unitig), marR:G103S (unitig), locus_kmer31_poly_0029`
- **Estabilidade entre dobras (Spearman $\rho$)**: `0.4184`
- **Concordância com GWAS controlado por estrutura populacional**: `83.3%`
- **Teste de remoção in silico do locus principal**: queda média de probabilidade = `-0.5608` (Coerente: `True`)

| Locus Genômico | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Associado |
|---|---|---|---|---|---|
| `gyrA:S83L` | mutacao_pontual | **conhecida** | AMRFinderPlus | 1.4252 | Alteração do alvo DNA girase subunidade A (QRDR S83L) |
| `acrR:fs*45 (unitig)` | kmer_unitig | **candidata** | Descoberta ML (F3) | 1.3396 | Desrepressão de bomba de efluxo multidroga AcrAB |
| `marR:G103S (unitig)` | kmer_unitig | **candidata** | Descoberta ML (F3) | 0.7420 | Regulador reprimido de marAB / bomba de efluxo AcrAB-TolC |
| `gyrA:D87N` | mutacao_pontual | **conhecida** | AMRFinderPlus | 0.4377 | Alteração do alvo DNA girase subunidade A (QRDR D87N) |
| `parC:S80I` | mutacao_pontual | **conhecida** | AMRFinderPlus | 0.4272 | Alteração da topoisomerase IV subunidade A (QRDR S80I) |
| `locus_kmer31_poly_0029` | kmer_unitig | **candidata** | Alinhamento Referência (k=31) | 0.2917 | Variante genômica não catalogada (hipótese) |

### Fármaco: `gentamicina` (`ecoli-gen-1.2.0`)

- **Isolados únicos pós-QC/desduplicação**: `360` (Prevalência de Resistência: `48.6%`)
- **Conclusão Q1 (ML vs Baseline B1 em grupos filogenéticos não vistos)**: B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.0973, IC95% [+0.0393, +0.1689], Delta VME = -0.2181).
- **Regressão Quantitativa de $\log_2(\text{MIC})$ (ISO 20776-2)**: Concordância Essencial ($\pm 1$ diluição) = `61.9%` | Concordância Categórica = `81.7%` | RMSE $\log_2(\text{MIC})$ = `2.236`

#### Escada de Modelos (Validação Principal: Grupo Filogenético — REG-01, REG-04, REG-05)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | VME (Erro Muito Grave) | ME (Erro Grave) | Brier Score |
|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | 0.4397 [0.349, 0.505] | 0.4189 | 0.2686 | 0.6054 | **0.7314** | **0.3946** | 0.2529 |
| `B1_Regras_CARD_AMRFinder` | 0.8174 [0.701, 0.892] | 0.7812 | 0.6457 | 0.9135 | **0.3543** | **0.0865** | 0.1640 |
| `B2_LogReg_L1_F1_F2` | 0.8344 [0.707, 0.889] | 0.7941 | 0.6171 | 0.9081 | **0.3829** | **0.0919** | 0.1674 |
| `B3_XGBoost_F1_F2` | 0.8282 [0.697, 0.888] | 0.7862 | 0.6343 | 0.8811 | **0.3657** | **0.1189** | 0.1728 |
| `B4_XGBoost_F1_F2_F3` | 0.9065 [0.837, 0.944] | 0.9157 | 0.8571 | 0.8649 | **0.1429** | **0.1351** | 0.1086 |
| `B5_Linear_Esparso_F3` | 0.7119 [0.658, 0.771] | 0.6935 | 0.6057 | 0.6811 | **0.3943** | **0.3189** | 0.2304 |
| `B6_CNN1D_Sequencia` | 0.8920 [0.810, 0.930] | 0.8938 | 0.8171 | 0.8054 | **0.1829** | **0.1946** | 0.1567 |
| `Controle_Kmers_Curtos_k6` | 0.5205 [0.385, 0.624] | 0.5112 | 0.4457 | 0.5568 | **0.5543** | **0.4432** | 0.2767 |

#### Generalização entre Estratégias de Split (Q2 e REG-09)

| Estratégia de Validação | AUPRC | AUROC | VME | ME | Observação |
|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | 0.9065 | 0.9157 | 0.1429 | 0.1351 | Sem vazamento clonal |
| **Leave-One-Lineage-Out (LOLO)** | 0.8958 | 0.9186 | 0.1695 | 0.1311 | Linhagens ST inteiras inéditas |
| **Temporal (Treino <= 2022, Teste > 2022)** | 0.9381 | 0.9257 | 0.1250 | 0.0870 | Simula vigilância prospectiva |
| **Geográfica (Países retidos)** | 0.8960 | 0.9046 | 0.1954 | 0.1348 | Robustez a epidemiologia local |
| **Aleatória (Demonstração de Vazamento B5)** | 0.7408 | 0.7518 | 0.3314 | 0.2973 | Inflação por clones no treino e teste |

#### Confiabilidade, Calibração e Abstenção Conformal (Q4 e REG-11)

- **Brier Score**: `0.1086` | **ECE**: `0.0356`
- **Cobertura Conformal**: `93.3%` (Taxa de Abstenção em casos duvidosos: `6.7%`)
- **Redução de Erros Graves pós-abstenção**: VME cai de `0.1429` para `0.1118`, e ME cai de `0.1351` para `0.1257`.

#### Avaliação de Explicabilidade TreeSHAP Agregada por Locus (Q3, Seção 8.2 e REG-10)

- **Proporção do Top-6 presente no CARD/AMRFinderPlus**: `66.7%` (aac(3)-IIa, ant(2'')-Ia, armA, aac(3)-IId)
- **Candidatos novos identificados via k-mers/unitigs (F3)**: `acrR:fs*45 (unitig), cpxA:Y144C (unitig)`
- **Estabilidade entre dobras (Spearman $\rho$)**: `0.4200`
- **Concordância com GWAS controlado por estrutura populacional**: `83.3%`
- **Teste de remoção in silico do locus principal**: queda média de probabilidade = `-0.3262` (Coerente: `True`)

| Locus Genômico | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Associado |
|---|---|---|---|---|---|
| `acrR:fs*45 (unitig)` | kmer_unitig | **candidata** | Descoberta ML (F3) | 0.9693 | Desrepressão de bomba de efluxo multidroga AcrAB |
| `aac(3)-IIa` | gene_amr | **conhecida** | CARD | 0.6578 | Aminoglicosídeo 3-N-acetiltransferase |
| `cpxA:Y144C (unitig)` | kmer_unitig | **candidata** | Descoberta ML (F3) | 0.6226 | Sensor quinase de envelope associado a tolerância a aminoglicosídeos |
| `ant(2'')-Ia` | gene_amr | **conhecida** | CARD | 0.6169 | Aminoglicosídeo 2''-O-nucleotidiltransferase |
| `armA` | gene_amr | **conhecida** | CARD | 0.3557 | Metiltransferase de rRNA 16S |
| `aac(3)-IId` | gene_amr | **conhecida** | AMRFinderPlus | 0.3373 | Aminoglicosídeo 3-N-acetiltransferase variante IId |

### Fármaco: `cefuroxima` (`ecoli-cxm-1.2.0`)

- **Isolados únicos pós-QC/desduplicação**: `360` (Prevalência de Resistência: `39.4%`)
- **Conclusão Q1 (ML vs Baseline B1 em grupos filogenéticos não vistos)**: B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.2023, IC95% [+0.1210, +0.3041], Delta VME = -0.3197).
- **Regressão Quantitativa de $\log_2(\text{MIC})$ (ISO 20776-2)**: Concordância Essencial ($\pm 1$ diluição) = `76.4%` | Concordância Categórica = `90.0%` | RMSE $\log_2(\text{MIC})$ = `1.592`

#### Escada de Modelos (Validação Principal: Grupo Filogenético — REG-01, REG-04, REG-05)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | VME (Erro Muito Grave) | ME (Erro Grave) | Brier Score |
|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | 0.3870 [0.340, 0.457] | 0.4845 | 0.0000 | 1.0000 | **1.0000** | **0.0000** | 0.2390 |
| `B1_Regras_CARD_AMRFinder` | 0.7364 [0.610, 0.825] | 0.7382 | 0.5352 | 0.9358 | **0.4648** | **0.0642** | 0.1699 |
| `B2_LogReg_L1_F1_F2` | 0.7328 [0.560, 0.851] | 0.7052 | 0.4930 | 0.9358 | **0.5070** | **0.0642** | 0.1763 |
| `B3_XGBoost_F1_F2` | 0.7379 [0.583, 0.853] | 0.7115 | 0.5423 | 0.8853 | **0.4577** | **0.1147** | 0.1834 |
| `B4_XGBoost_F1_F2_F3` | 0.9335 [0.904, 0.958] | 0.9471 | 0.8592 | 0.9450 | **0.1408** | **0.0550** | 0.0726 |
| `B5_Linear_Esparso_F3` | 0.6887 [0.643, 0.740] | 0.6884 | 0.5986 | 0.6468 | **0.4014** | **0.3532** | 0.2275 |
| `B6_CNN1D_Sequencia` | 0.9176 [0.881, 0.945] | 0.9345 | 0.8380 | 0.8853 | **0.1620** | **0.1147** | 0.1142 |
| `Controle_Kmers_Curtos_k6` | 0.3799 [0.332, 0.439] | 0.4753 | 0.3239 | 0.6284 | **0.6761** | **0.3716** | 0.2778 |

#### Generalização entre Estratégias de Split (Q2 e REG-09)

| Estratégia de Validação | AUPRC | AUROC | VME | ME | Observação |
|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | 0.9335 | 0.9471 | 0.1408 | 0.0550 | Sem vazamento clonal |
| **Leave-One-Lineage-Out (LOLO)** | 0.8852 | 0.9159 | 0.1500 | 0.0250 | Linhagens ST inteiras inéditas |
| **Temporal (Treino <= 2022, Teste > 2022)** | 0.9495 | 0.9488 | 0.1111 | 0.0517 | Simula vigilância prospectiva |
| **Geográfica (Países retidos)** | 0.9261 | 0.9332 | 0.1515 | 0.0545 | Robustez a epidemiologia local |
| **Aleatória (Demonstração de Vazamento B5)** | 0.6842 | 0.7072 | 0.4437 | 0.2615 | Inflação por clones no treino e teste |

#### Confiabilidade, Calibração e Abstenção Conformal (Q4 e REG-11)

- **Brier Score**: `0.0726` | **ECE**: `0.0209`
- **Cobertura Conformal**: `95.3%` (Taxa de Abstenção em casos duvidosos: `4.7%`)
- **Redução de Erros Graves pós-abstenção**: VME cai de `0.1408` para `0.1212`, e ME cai de `0.0550` para `0.0474`.

#### Avaliação de Explicabilidade TreeSHAP Agregada por Locus (Q3, Seção 8.2 e REG-10)

- **Proporção do Top-6 presente no CARD/AMRFinderPlus**: `50.0%` (blaCTX-M-14, blaCTX-M-15, ampC:c.-42C>T)
- **Candidatos novos identificados via k-mers/unitigs (F3)**: `ompF:c.-112del (unitig), ftsI (PBP3):ins333YRIN (unitig), locus_kmer31_poly_0032`
- **Estabilidade entre dobras (Spearman $\rho$)**: `0.4167`
- **Concordância com GWAS controlado por estrutura populacional**: `83.3%`
- **Teste de remoção in silico do locus principal**: queda média de probabilidade = `-0.3545` (Coerente: `True`)

| Locus Genômico | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Associado |
|---|---|---|---|---|---|
| `ompF:c.-112del (unitig)` | kmer_unitig | **candidata** | Descoberta ML (F3) | 0.9043 | Redução de permeabilidade de porina de membrana externa |
| `blaCTX-M-14` | gene_amr | **conhecida** | CARD | 0.7801 | Beta-lactamase de espectro estendido (ESBL CTX-M-14) |
| `blaCTX-M-15` | gene_amr | **conhecida** | CARD | 0.6904 | Beta-lactamase de espectro estendido (ESBL CTX-M-15) |
| `ftsI (PBP3):ins333YRIN (unitig)` | kmer_unitig | **candidata** | Descoberta ML (F3) | 0.5812 | Redução de afinidade de PBP3 a cefalosporinas |
| `ampC:c.-42C>T` | mutacao_pontual | **conhecida** | AMRFinderPlus | 0.5351 | Superexpressão de beta-lactamase cromossômica AmpC por mutação promotora |
| `locus_kmer31_poly_0032` | kmer_unitig | **candidata** | Alinhamento Referência (k=31) | 0.2641 | Variante genômica não catalogada (hipótese) |


---

## 3. Extensões da Fase 6: Multi-Espécie (*K. pneumoniae* e *S. aureus*) e Validação Externa

### 3.1 Expansão Multi-Espécie (Seção 3.2 — Fases 2 e 3)

| Espécie | Fármaco | Versão do Modelo | AUPRC Baseline B1 | AUPRC Modelo B4 ($\Delta$) | VME | Top Loci TreeSHAP |
|---|---|---|---|---|---|---|
| *Klebsiella pneumoniae* | `ciprofloxacina` | `kpneu-cip-1.2.0` | 0.9356 | **0.9939** (+0.0588) | 0.0508 | gyrA:S83L, marR:G103S (unitig), qnrS1 |
| *Klebsiella pneumoniae* | `gentamicina` | `kpneu-gen-1.2.0` | 0.7784 | **0.8864** (+0.1030) | 0.1429 | aac(3)-IIa, cpxA:Y144C (unitig), armA |
| *Klebsiella pneumoniae* | `ceftriaxona` | `kpneu-cro-1.2.0` | 0.8789 | **0.9206** (+0.0464) | 0.0492 | blaCTX-M-15, ompF:c.-112del (unitig), locus_kmer31_poly_0049 |
| *Staphylococcus aureus* | `oxacilina` | `saureus-oxa-1.2.0` | 0.7910 | **0.8977** (+0.1066) | 0.0667 | mecA, pbp4:c.-90T>A (unitig), locus_kmer31_poly_0041 |
| *Staphylococcus aureus* | `ciprofloxacina` | `saureus-cip-1.2.0` | 0.8067 | **0.9097** (+0.1048) | 0.0862 | grlA:S80F, norA:flqB_up (unitig), locus_kmer31_poly_0118 |

### 3.2 Validação Externa Independente — Padrão NCBI Pathogen Detection AST (Seção 3.1, 7.1 e Fase 6)

| Fármaco (*E. coli*) | AUPRC B1 (Externo) | AUPRC B4 (Externo) | Ganho $\Delta\text{AUPRC}$ Externo | VME Externo | ME Externo |
|---|---|---|---|---|---|
| `ciprofloxacina` | 0.8542 | **0.9778** | **+0.1236** | 0.0465 | 0.1228 |
| `gentamicina` | 0.7878 | **0.9395** | **+0.1517** | 0.1875 | 0.1154 |
| `cefuroxima` | 0.7315 | **0.9215** | **+0.1900** | 0.2553 | 0.0000 |

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
