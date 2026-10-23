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

## 2. Resultados por Fármaco (Perguntas Científicas Q1, Q2, Q3 e Q4)

### Fármaco: `ciprofloxacina` (`ecoli-cip-1.2.0`)

- **Isolados únicos pós-QC/desduplicação**: `360` (Prevalência de Resistência: `48.3%`)
- **Conclusão Q1 (ML vs Baseline B1 em grupos filogenéticos não vistos)**: B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.0875, IC95% [+0.0402, +0.1396], Delta VME = -0.1782).

#### Escada de Modelos (Validação Principal: Grupo Filogenético — REG-01, REG-04, REG-05)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | VME (Erro Muito Grave) | ME (Erro Grave) | Brier Score |
|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | 0.4486 [0.404, 0.510] | 0.4323 | 0.0000 | 1.0000 | **1.0000** | **0.0000** | 0.2525 |
| `B1_Regras_CARD_AMRFinder` | 0.8763 [0.820, 0.922] | 0.8499 | 0.7471 | 0.9677 | **0.2529** | **0.0323** | 0.1143 |
| `B2_LogReg_L1_F1_F2` | 0.9005 [0.843, 0.944] | 0.8691 | 0.7356 | 0.9731 | **0.2644** | **0.0269** | 0.1118 |
| `B3_XGBoost_F1_F2` | 0.8987 [0.837, 0.941] | 0.8665 | 0.7356 | 0.9570 | **0.2644** | **0.0430** | 0.1163 |
| `B4_XGBoost_F1_F2_F3` | 0.9621 [0.946, 0.978] | 0.9644 | 0.9253 | 0.9247 | **0.0747** | **0.0753** | 0.0553 |
| `B5_Linear_Esparso_F3` | 0.7602 [0.704, 0.816] | 0.7249 | 0.5747 | 0.8011 | **0.4253** | **0.1989** | 0.2126 |
| `Controle_Kmers_Curtos_k6` | 0.4862 [0.417, 0.562] | 0.4840 | 0.4483 | 0.4946 | **0.5517** | **0.5054** | 0.2760 |

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

#### Escada de Modelos (Validação Principal: Grupo Filogenético — REG-01, REG-04, REG-05)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | VME (Erro Muito Grave) | ME (Erro Grave) | Brier Score |
|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | 0.4397 [0.365, 0.531] | 0.4189 | 0.2686 | 0.6054 | **0.7314** | **0.3946** | 0.2529 |
| `B1_Regras_CARD_AMRFinder` | 0.8174 [0.701, 0.892] | 0.7812 | 0.6457 | 0.9135 | **0.3543** | **0.0865** | 0.1640 |
| `B2_LogReg_L1_F1_F2` | 0.8344 [0.712, 0.909] | 0.7941 | 0.6171 | 0.9081 | **0.3829** | **0.0919** | 0.1674 |
| `B3_XGBoost_F1_F2` | 0.8282 [0.705, 0.896] | 0.7862 | 0.6343 | 0.8811 | **0.3657** | **0.1189** | 0.1728 |
| `B4_XGBoost_F1_F2_F3` | 0.9065 [0.837, 0.944] | 0.9157 | 0.8571 | 0.8649 | **0.1429** | **0.1351** | 0.1086 |
| `B5_Linear_Esparso_F3` | 0.7119 [0.646, 0.780] | 0.6935 | 0.6057 | 0.6811 | **0.3943** | **0.3189** | 0.2304 |
| `Controle_Kmers_Curtos_k6` | 0.5205 [0.407, 0.642] | 0.5112 | 0.4457 | 0.5568 | **0.5543** | **0.4432** | 0.2767 |

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

#### Escada de Modelos (Validação Principal: Grupo Filogenético — REG-01, REG-04, REG-05)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | VME (Erro Muito Grave) | ME (Erro Grave) | Brier Score |
|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | 0.3870 [0.327, 0.458] | 0.4845 | 0.0000 | 1.0000 | **1.0000** | **0.0000** | 0.2390 |
| `B1_Regras_CARD_AMRFinder` | 0.7364 [0.610, 0.825] | 0.7382 | 0.5352 | 0.9358 | **0.4648** | **0.0642** | 0.1699 |
| `B2_LogReg_L1_F1_F2` | 0.7328 [0.595, 0.846] | 0.7052 | 0.4930 | 0.9358 | **0.5070** | **0.0642** | 0.1763 |
| `B3_XGBoost_F1_F2` | 0.7379 [0.601, 0.846] | 0.7115 | 0.5423 | 0.8853 | **0.4577** | **0.1147** | 0.1834 |
| `B4_XGBoost_F1_F2_F3` | 0.9335 [0.904, 0.958] | 0.9471 | 0.8592 | 0.9450 | **0.1408** | **0.0550** | 0.0726 |
| `B5_Linear_Esparso_F3` | 0.6887 [0.639, 0.744] | 0.6884 | 0.5986 | 0.6468 | **0.4014** | **0.3532** | 0.2275 |
| `Controle_Kmers_Curtos_k6` | 0.3799 [0.322, 0.453] | 0.4753 | 0.3239 | 0.6284 | **0.6761** | **0.3716** | 0.2778 |

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

## 3. Limitações e Uso Responsável (Seção 15)

1. **Associação vs. Causalidade (REG-10)**: Os valores TreeSHAP medem contribuição preditiva ao modelo; loci classificados como `candidata` representam hipóteses biológicas para validação experimental (ex.: mutagênese sítio-dirigida ou ensaio funcional).
2. **Domínio de Aplicação**: Os modelos avaliados aplicam-se exclusivamente a isolados montados de *Escherichia coli* aprovados pelos critérios de QC da Seção 3.4 para os fármacos `ciprofloxacina`, `gentamicina` e `cefuroxima`.
3. **Pesquisa Científica (REG-15)**: Este sistema não é um dispositivo médico e não deve orientar conduta clínica.
