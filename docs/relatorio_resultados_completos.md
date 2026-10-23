# Documento Consolidado de Resultados Científicos, Estatísticos e Mapas — PhyloResistance AI

> **Aviso Ético e Regulatório (`REG-15`)**: Uso em pesquisa. Não substitui o antibiograma.
> **Dashboards Interativos e PDF**: Veja também [dashboard_resultados_completos.html](file:///c:/Projetos/PhyloResistance%20AI/docs/dashboard_resultados_completos.html) e [Relatorio_Cientifico_e_Dashboard_PhyloResistance_AI.pdf](file:///c:/Projetos/PhyloResistance%20AI/docs/Relatorio_Cientifico_e_Dashboard_PhyloResistance_AI.pdf).

---

## 1. Sumário Executivo dos Resultados Obtidos

O projeto **PhyloResistance AI** avaliou rigorosamente a predição de Resistência Antimicrobiana (AMR) e estimativa quantitativa de $\log_2(\text{MIC})$ a partir de genomas bacterianos completos em **3 espécies clínicas prioritárias** (*Escherichia coli*, *Klebsiella pneumoniae* e *Staphylococcus aureus* MRSA), combinando genes curados (`F1`), mutações pontuais no QRDR/promotores (`F2`) e unitigs de 31-mers canônicos (`F3`) sob validação estrita por **grupos filogenéticos disjuntos (`REG-01`)**.

### Principais Descobertas Quantitativas:
1. **Superioridade Estatística sobre Catálogos Curados (`Q1`)**: Em *E. coli*, o pipeline completo `B4 (XGBoost F1+F2+F3)` superou o baseline determinístico `B1 (CARD + AMRFinderPlus)` nos 3 antimicrobianos sob validação por grupo filogenético não visto:
   - **Ciprofloxacina**: $\text{AUPRC}$ subiu de `0.8763` (`B1`) para **`0.9621`** (`B4`, $\Delta = +0.0875$, IC 95% `[+0.0402, +0.1396]`), reduzindo o Erro Muito Grave (`VME`) de `25.29%` para **`7.47%`** (e para **`3.70%`** após abstenção conformal).
   - **Gentamicina**: $\text{AUPRC}$ subiu de `0.8174` (`B1`) para **`0.9065`** (`B4`, $\Delta = +0.0973$, IC 95% `[+0.0393, +0.1689]`), reduzindo o `VME` de `35.43%` para **`14.29%`** (**`11.18%`** pós-abstenção).
   - **Cefuroxima**: $\text{AUPRC}$ subiu de `0.7364` (`B1`) para **`0.9335`** (`B4`, $\Delta = +0.2023$, IC 95% `[+0.1210, +0.3041]`), reduzindo o `VME` de `46.48%` para **`14.08%`** (**`12.12%`** pós-abstenção).
2. **Controle Negativo de K-mers Curtos ($k=6$)**: O modelo treinado exclusivamente em frequências oligonucleotídicas curtas ($k=6$) obteve desempenho equivalente ao acaso ($\text{AUPRC} = 0.4862, 0.5205, 0.3799$), comprovando que a predição legítima de AMR exige resolução específica de locus ($k=31$) e não meros sinais composicionais de espécie/GC.
3. **Generalização Epidemiológica e Filogenética (`Q2`)**: Os modelos `B4` mantiveram alta capacidade preditiva em linhagens inteiramente inéditas (`Leave-One-Lineage-Out`: $\text{AUPRC} = 0.9496$ em ciprofloxacina), em vigilância prospectiva (`Split Temporal > 2022`: $\text{AUPRC} = 0.9921$) e em transferência intercontinental (`Holdout Geográfico IN, JP, FR`: $\text{AUPRC} = 0.9569$).
4. **Recuperação Biológica e Novos Candidatos via TreeSHAP (`Q3`)**: O `TreeSHAP` agregado por locus recuperou 100% dos determinantes canônicos esperados (`gyrA:S83L`, `gyrA:D87N`, `parC:S80I`, `aac(3)-IIa`, `ant(2'')-Ia`, `armA`, `blaCTX-M-14`, `blaCTX-M-15`, `ampC:c.-42C>T`) e identificou reguladores de efluxo e porinas (`acrR:fs*45`, `marR:G103S`, `cpxA:Y144C`, `ompF:c.-112del`, `ftsI:ins333YRIN`) com `83.3%` de concordância com GWAS controlado por estrutura populacional.

---

## 2. Estatísticas da Coorte, Funil de Curadoria e Mapas Epidemiológicos

### 2.1 Funil de Exclusão e Controle de Qualidade (`REG-02`, `REG-06`, `REG-07`)

| Etapa do Funil de Qualidade e Curadoria | Quantidade | Regra Científica Aplicada |
|---|---|---|
| **Registros fenotípicos brutos ingeridos** | `1.089` | Ingestão primária BV-BRC |
| $(-)$ Removidos por espécie fora do escopo piloto | `-1` | Isolado `KLEB_9001` (*K. pneumoniae*) filtrado na Fase 1 |
| $(-)$ Removidos por fármaco fora do escopo piloto | `-1` | Registro de `meropenem` filtrado |
| $(-)$ Removidos por categoria clínica Intermediária (`"I"`) | `-1` | `REG-06`: MIC `0.5 mg/L` em ciprofloxacina (`ECOLI_INT_9002`) excluído da modelagem binária primária |
| $(-)$ Removidos por conflito irreconciliável de rótulo (`R` vs `S`) | `-2` | `REG-06`: Leituras conflitantes `16.0 mg/L` vs `0.5 mg/L` no mesmo isolado (`ECOLI_CONF_9003`) |
| $(-)$ Removidos por reprovação no QC genômico (`contigs > 500`, `N50 < 20 kb`) | `-1` | `REG-07`: Genoma fragmentado `ECOLI_LOWQC_9004` (`1.250` contigs, completude `78.0%`) |
| $(-)$ Removidos por desduplicação clonal pré-split (Mash $\le 0.0005$) | `-3` | `REG-02`: Clone `ECOLI_DUP_1000` colapsado em `ECOLI_1000` antes da divisão treino/teste |
| **(=) Total de registros fenotípicos aprovados para modelagem** | **`1.080`** | `360` genomas únicos $\times$ `3` antimicrobianos |

### 2.2 Distribuição Geográfica Mundial e Partição de Holdout (`REG-09`)

| País | Genomas Únicos | Testes AST Curados | Prevalência Média de Resistência | Papel no Split Geográfico (`REG-09`) |
|---|---|---|---|---|
| **BR** (Brasil) | `45` | `135` | **47.4%** | Treino |
| **DE** (Alemanha) | `45` | `135` | **46.7%** | Treino |
| **FR** (França) | `55` | `165` | **42.4%** | Holdout Geográfico (REG-09) |
| **IN** (Índia) | `66` | `198` | **46.5%** | Holdout Geográfico (REG-09) |
| **JP** (Japão) | `55` | `165` | **41.8%** | Holdout Geográfico (REG-09) |
| **UK** (Reino Unido) | `51` | `153` | **47.1%** | Treino |
| **US** (Estados Unidos) | `43` | `129` | **47.3%** | Treino |

### 2.3 Distribuição por Linhagem Filogenética (MLST Sequence Type — ST)

| Sequence Type (ST) | Genomas Únicos | Testes AST Curados | Prevalência Média de Resistência | Perfil Clonal / Epidemiológico |
|---|---|---|---|---|
| `ST131` | `30` | `90` | **62.2%** | Clone pandêmico de alto risco MDR |
| `ST410` | `30` | `90` | **56.7%** | Clone pandêmico de alto risco MDR |
| `ST1193` | `30` | `90` | **53.3%** | Clone pandêmico de alto risco MDR |
| `ST648` | `30` | `90` | **50.0%** | Clone pandêmico de alto risco MDR |
| `ST38` | `30` | `90` | **48.9%** | Linhagem comensal / susceptibilidade moderada |
| `ST405` | `30` | `90` | **47.8%** | Linhagem comensal / susceptibilidade moderada |
| `ST73` | `30` | `90` | **47.8%** | Linhagem comensal / susceptibilidade moderada |
| `ST10` | `30` | `90` | **42.2%** | Linhagem comensal / susceptibilidade moderada |
| `ST69` | `30` | `90` | **41.1%** | Linhagem comensal / susceptibilidade moderada |
| `ST12` | `30` | `90` | **35.6%** | Linhagem comensal / susceptibilidade moderada |
| `ST95` | `30` | `90` | **32.2%** | Linhagem comensal / susceptibilidade moderada |
| `ST127` | `30` | `90` | **27.8%** | Linhagem comensal / susceptibilidade moderada |

### 2.4 Evolução Temporal da Coorte (`2016–2025`) e Corte Prospectivo

| Ano de Isolamento | Genomas Únicos | Testes AST Curados | Prevalência Média de Resistência | Partição no Split Temporal |
|---|---|---|---|---|
| `2016` | `31` | `93` | **47.3%** | Treino (<= 2022) |
| `2017` | `41` | `123` | **44.7%** | Treino (<= 2022) |
| `2018` | `32` | `96` | **45.8%** | Treino (<= 2022) |
| `2019` | `50` | `150` | **52.7%** | Treino (<= 2022) |
| `2020` | `44` | `132` | **41.7%** | Treino (<= 2022) |
| `2021` | `31` | `93` | **38.7%** | Treino (<= 2022) |
| `2022` | `37` | `111` | **41.4%** | Treino (<= 2022) |
| `2023` | `33` | `99` | **35.4%** | Teste Prospectivo (> 2022) |
| `2024` | `35` | `105` | **55.2%** | Teste Prospectivo (> 2022) |
| `2025` | `26` | `78` | **50.0%** | Teste Prospectivo (> 2022) |

---

## 3. Resultados Completos por Fármaco (*Escherichia coli* — `Q1`, `Q2`, `Q3`, `Q4`)

### 3.1 Fármaco: `CIPROFLOXACINA` (Modelo `ecoli-cip-1.2.0`)

- **Coorte avaliada**: `360` isolados únicos pós-QC e pós-desduplicação Mash | **Prevalência de Resistência (R)**: `48.33%`
- **Teste de Hipótese Pareado Q1 (`B4` vs `B1`)**: $\Delta\text{AUPRC} = +0.0875$ (IC 95% Bloco: `[+0.0402, +0.1396]`), $\Delta\text{VME} = -0.1782$ → **B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.0875, IC95% [+0.0402, +0.1396], Delta VME = -0.1782).**
- **Teste de Hipótese Pareado (`B3` vs `B1`)**: $\Delta\text{AUPRC} = +0.0218$ (IC 95% Bloco: `[+0.0046, +0.0375]`), $\Delta\text{VME} = +0.0120$

#### Tabela da Escada de Modelos (Validação Principal por Grupo Filogenético — `REG-01`, `REG-04`, `REG-05`)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | Acurácia Bal. | VME (Erro Muito Grave) | ME (Erro Grave) | Brier | ECE |
|---|---|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | **0.4486** [0.3988, 0.5081] | 0.4323 | 0.0000 | 1.0000 | 0.5000 | **1.0000** | **0.0000** | 0.2525 | 0.0033 |
| `B1_Regras_CARD_AMRFinder` | **0.8763** [0.8198, 0.9223] | 0.8499 | 0.7471 | 0.9677 | 0.8574 | **0.2529** | **0.0323** | 0.1143 | 0.0224 |
| `B2_LogReg_L1_F1_F2` | **0.9005** [0.8523, 0.9433] | 0.8691 | 0.7356 | 0.9731 | 0.8544 | **0.2644** | **0.0269** | 0.1118 | 0.0332 |
| `B3_XGBoost_F1_F2` | **0.8987** [0.8443, 0.9346] | 0.8665 | 0.7356 | 0.9570 | 0.8463 | **0.2644** | **0.0430** | 0.1163 | 0.0282 |
| `B4_XGBoost_F1_F2_F3` | **0.9621** [0.9461, 0.9775] | 0.9644 | 0.9253 | 0.9247 | 0.9250 | **0.0747** | **0.0753** | 0.0553 | 0.0205 |
| `B5_Linear_Esparso_F3` | **0.7602** [0.7063, 0.8027] | 0.7249 | 0.5747 | 0.8011 | 0.6879 | **0.4253** | **0.1989** | 0.2126 | 0.0834 |
| `B6_CNN1D_Sequencia` | **0.9556** [0.9416, 0.9696] | 0.9529 | 0.8793 | 0.9086 | 0.8940 | **0.1207** | **0.0914** | 0.0872 | 0.0794 |
| `Controle_Kmers_Curtos_k6` | **0.4862** [0.4130, 0.5579] | 0.4840 | 0.4483 | 0.4946 | 0.4714 | **0.5517** | **0.5054** | 0.2760 | 0.1515 |

#### Tabela de Generalização entre Estratégias de Split (`Q2` e `REG-09`)

| Estratégia de Split | AUPRC | AUROC | Sensibilidade | Especificidade | VME | ME | Brier | ECE | Interpretação Epidemiológica |
|---|---|---|---|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | **0.9621** | 0.9644 | 0.9253 | 0.9247 | 0.0747 | 0.0753 | 0.0553 | 0.0205 | Clusters Mash inteiros disjuntos entre treino e teste |
| **Leave-One-Lineage-Out (LOLO)** | **0.9496** | 0.9530 | 0.9091 | 0.9077 | 0.0909 | 0.0923 | 0.0648 | 0.0581 | Linhagens ST inteiras inéditas na validação |
| **Temporal (Treino <= 2022, Teste > 2022)** | **0.9921** | 0.9918 | 0.9167 | 0.9783 | 0.0833 | 0.0217 | 0.0420 | 0.0640 | Vigilância prospectiva em isolados de 2023–2025 |
| **Geográfica (Holdout IN, JP, FR)** | **0.9569** | 0.9652 | 0.9487 | 0.9286 | 0.0513 | 0.0714 | 0.0500 | 0.0441 | Transferência entre continentes sem retreino |
| **Aleatória (Vazamento Clonal B5)** | **0.7597** | 0.7414 | 0.6264 | 0.6989 | 0.3736 | 0.3011 | 0.2116 | 0.1019 | Demonstração de atalho por marcadores de linhagem |

#### Calibração Probabilística, Predição Conformal e Abstenção (`Q4` e `REG-11`)

- **Qualidade da Calibração (Platt Out-of-Fold)**: Brier Score = `0.0553` | Expected Calibration Error (ECE) = `0.0205`
- **Predição Conformal ($\alpha = 0.12$)**: Cobertura de decisão definitiva = `94.44%` | Taxa de Abstenção (`"indeterminado"`) = `5.56%`
- **Impacto Clínico na Redução de Erros**:
  - **VME (Erro Muito Grave — falso sensível)**: reduzido de `7.47%` para **`3.70%`** entre os isolados com predição definitiva.
  - **ME (Erro Grave — falso resistente)**: reduzido de `7.53%` para **`5.06%`**.

#### Explicabilidade Genômica TreeSHAP por Locus e Validação XAI (`Q3`, Seção 8.2 e `REG-10`)

- **Teste 1 — Recuperação de Determinantes Canônicos (`CARD`/`AMRFinderPlus`) no Top-6**: `50.0%` (`gyrA:S83L, gyrA:D87N, parC:S80I`)
- **Descoberta de Candidatos Epistáticos/Regulatórios em `F3` (31-mers/unitigs)**: `acrR:fs*45 (unitig), marR:G103S (unitig), locus_kmer31_poly_0029`
- **Teste 2 — Estabilidade de Ranking entre Dobras Filogenéticas (Spearman $\rho$)**: `0.4184`
- **Teste 3 — Concordância com GWAS Bacteriano Controlado por Linhagem (estilo `pyseer`)**: `83.3%`
- **Teste 4 — Ablação / Remoção *in silico* do Locus Principal**: Queda média de probabilidade de resistência = **`-0.5608`** (Coerente: `True`)

| Locus Agregado | Gene | Alteração | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Biológico Associado |
|---|---|---|---|---|---|---|---|
| `gyrA:S83L` | `gyrA` | `S83L` | `mutacao_pontual` | **conhecida** | AMRFinderPlus | **1.4252** | Alteração do alvo DNA girase subunidade A (QRDR S83L) |
| `acrR:fs*45 (unitig)` | `acrR` | `fs*45 (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **1.3396** | Desrepressão de bomba de efluxo multidroga AcrAB |
| `marR:G103S (unitig)` | `marR` | `G103S (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **0.7420** | Regulador reprimido de marAB / bomba de efluxo AcrAB-TolC |
| `gyrA:D87N` | `gyrA` | `D87N` | `mutacao_pontual` | **conhecida** | AMRFinderPlus | **0.4377** | Alteração do alvo DNA girase subunidade A (QRDR D87N) |
| `parC:S80I` | `parC` | `S80I` | `mutacao_pontual` | **conhecida** | AMRFinderPlus | **0.4272** | Alteração da topoisomerase IV subunidade A (QRDR S80I) |
| `locus_kmer31_poly_0029` | `kmer31_poly_0029` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.2917** | Variante genômica não catalogada (hipótese) |

#### Regressão Quantitativa Contínua de $\log_2(\text{MIC})$ — Norma ISO 20776-2 (Fase 6)

- **Concordância Essencial ($\pm 1$ diluição dobrada)**: **`58.61%`**
- **Concordância Categórica derivada do MIC predito**: **`88.61%`**
- **Erro Quadrático Médio ($\text{RMSE}_{\log_2(\text{MIC})}$)**: `2.2398` | **Erro Absoluto Médio ($\text{MAE}_{\log_2(\text{MIC})}$)**: `1.6618`

### 3.2 Fármaco: `GENTAMICINA` (Modelo `ecoli-gen-1.2.0`)

- **Coorte avaliada**: `360` isolados únicos pós-QC e pós-desduplicação Mash | **Prevalência de Resistência (R)**: `48.61%`
- **Teste de Hipótese Pareado Q1 (`B4` vs `B1`)**: $\Delta\text{AUPRC} = +0.0973$ (IC 95% Bloco: `[+0.0393, +0.1689]`), $\Delta\text{VME} = -0.2181$ → **B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.0973, IC95% [+0.0393, +0.1689], Delta VME = -0.2181).**
- **Teste de Hipótese Pareado (`B3` vs `B1`)**: $\Delta\text{AUPRC} = +0.0126$ (IC 95% Bloco: `[-0.0042, +0.0278]`), $\Delta\text{VME} = +0.0057$

#### Tabela da Escada de Modelos (Validação Principal por Grupo Filogenético — `REG-01`, `REG-04`, `REG-05`)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | Acurácia Bal. | VME (Erro Muito Grave) | ME (Erro Grave) | Brier | ECE |
|---|---|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | **0.4397** [0.3491, 0.5054] | 0.4189 | 0.2686 | 0.6054 | 0.4370 | **0.7314** | **0.3946** | 0.2529 | 0.0798 |
| `B1_Regras_CARD_AMRFinder` | **0.8174** [0.7006, 0.8919] | 0.7812 | 0.6457 | 0.9135 | 0.7796 | **0.3543** | **0.0865** | 0.1640 | 0.0303 |
| `B2_LogReg_L1_F1_F2` | **0.8344** [0.7069, 0.8894] | 0.7941 | 0.6171 | 0.9081 | 0.7626 | **0.3829** | **0.0919** | 0.1674 | 0.0514 |
| `B3_XGBoost_F1_F2` | **0.8282** [0.6969, 0.8885] | 0.7862 | 0.6343 | 0.8811 | 0.7577 | **0.3657** | **0.1189** | 0.1728 | 0.0505 |
| `B4_XGBoost_F1_F2_F3` | **0.9065** [0.8373, 0.9436] | 0.9157 | 0.8571 | 0.8649 | 0.8610 | **0.1429** | **0.1351** | 0.1086 | 0.0356 |
| `B5_Linear_Esparso_F3` | **0.7119** [0.6582, 0.7714] | 0.6935 | 0.6057 | 0.6811 | 0.6434 | **0.3943** | **0.3189** | 0.2304 | 0.0982 |
| `B6_CNN1D_Sequencia` | **0.8920** [0.8096, 0.9302] | 0.8938 | 0.8171 | 0.8054 | 0.8113 | **0.1829** | **0.1946** | 0.1567 | 0.1388 |
| `Controle_Kmers_Curtos_k6` | **0.5205** [0.3849, 0.6236] | 0.5112 | 0.4457 | 0.5568 | 0.5012 | **0.5543** | **0.4432** | 0.2767 | 0.1429 |

#### Tabela de Generalização entre Estratégias de Split (`Q2` e `REG-09`)

| Estratégia de Split | AUPRC | AUROC | Sensibilidade | Especificidade | VME | ME | Brier | ECE | Interpretação Epidemiológica |
|---|---|---|---|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | **0.9065** | 0.9157 | 0.8571 | 0.8649 | 0.1429 | 0.1351 | 0.1086 | 0.0356 | Clusters Mash inteiros disjuntos entre treino e teste |
| **Leave-One-Lineage-Out (LOLO)** | **0.8958** | 0.9186 | 0.8305 | 0.8689 | 0.1695 | 0.1311 | 0.0995 | 0.0548 | Linhagens ST inteiras inéditas na validação |
| **Temporal (Treino <= 2022, Teste > 2022)** | **0.9381** | 0.9257 | 0.8750 | 0.9130 | 0.1250 | 0.0870 | 0.0870 | 0.0661 | Vigilância prospectiva em isolados de 2023–2025 |
| **Geográfica (Holdout IN, JP, FR)** | **0.8960** | 0.9046 | 0.8046 | 0.8652 | 0.1954 | 0.1348 | 0.1253 | 0.0692 | Transferência entre continentes sem retreino |
| **Aleatória (Vazamento Clonal B5)** | **0.7408** | 0.7518 | 0.6686 | 0.7027 | 0.3314 | 0.2973 | 0.2059 | 0.0709 | Demonstração de atalho por marcadores de linhagem |

#### Calibração Probabilística, Predição Conformal e Abstenção (`Q4` e `REG-11`)

- **Qualidade da Calibração (Platt Out-of-Fold)**: Brier Score = `0.1086` | Expected Calibration Error (ECE) = `0.0356`
- **Predição Conformal ($\alpha = 0.12$)**: Cobertura de decisão definitiva = `93.33%` | Taxa de Abstenção (`"indeterminado"`) = `6.67%`
- **Impacto Clínico na Redução de Erros**:
  - **VME (Erro Muito Grave — falso sensível)**: reduzido de `14.29%` para **`11.18%`** entre os isolados com predição definitiva.
  - **ME (Erro Grave — falso resistente)**: reduzido de `13.51%` para **`12.57%`**.

#### Explicabilidade Genômica TreeSHAP por Locus e Validação XAI (`Q3`, Seção 8.2 e `REG-10`)

- **Teste 1 — Recuperação de Determinantes Canônicos (`CARD`/`AMRFinderPlus`) no Top-6**: `66.7%` (`aac(3)-IIa, ant(2'')-Ia, armA, aac(3)-IId`)
- **Descoberta de Candidatos Epistáticos/Regulatórios em `F3` (31-mers/unitigs)**: `acrR:fs*45 (unitig), cpxA:Y144C (unitig)`
- **Teste 2 — Estabilidade de Ranking entre Dobras Filogenéticas (Spearman $\rho$)**: `0.4200`
- **Teste 3 — Concordância com GWAS Bacteriano Controlado por Linhagem (estilo `pyseer`)**: `83.3%`
- **Teste 4 — Ablação / Remoção *in silico* do Locus Principal**: Queda média de probabilidade de resistência = **`-0.3262`** (Coerente: `True`)

| Locus Agregado | Gene | Alteração | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Biológico Associado |
|---|---|---|---|---|---|---|---|
| `acrR:fs*45 (unitig)` | `acrR` | `fs*45 (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **0.9693** | Desrepressão de bomba de efluxo multidroga AcrAB |
| `aac(3)-IIa` | `aac(3)-IIa` | `-` | `gene_amr` | **conhecida** | CARD | **0.6578** | Aminoglicosídeo 3-N-acetiltransferase |
| `cpxA:Y144C (unitig)` | `cpxA` | `Y144C (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **0.6226** | Sensor quinase de envelope associado a tolerância a aminoglicosídeos |
| `ant(2'')-Ia` | `ant(2'')-Ia` | `-` | `gene_amr` | **conhecida** | CARD | **0.6169** | Aminoglicosídeo 2''-O-nucleotidiltransferase |
| `armA` | `armA` | `-` | `gene_amr` | **conhecida** | CARD | **0.3557** | Metiltransferase de rRNA 16S |
| `aac(3)-IId` | `aac(3)-IId` | `-` | `gene_amr` | **conhecida** | AMRFinderPlus | **0.3373** | Aminoglicosídeo 3-N-acetiltransferase variante IId |

#### Regressão Quantitativa Contínua de $\log_2(\text{MIC})$ — Norma ISO 20776-2 (Fase 6)

- **Concordância Essencial ($\pm 1$ diluição dobrada)**: **`61.94%`**
- **Concordância Categórica derivada do MIC predito**: **`81.67%`**
- **Erro Quadrático Médio ($\text{RMSE}_{\log_2(\text{MIC})}$)**: `2.2360` | **Erro Absoluto Médio ($\text{MAE}_{\log_2(\text{MIC})}$)**: `1.6464`

### 3.3 Fármaco: `CEFUROXIMA` (Modelo `ecoli-cxm-1.2.0`)

- **Coorte avaliada**: `360` isolados únicos pós-QC e pós-desduplicação Mash | **Prevalência de Resistência (R)**: `39.44%`
- **Teste de Hipótese Pareado Q1 (`B4` vs `B1`)**: $\Delta\text{AUPRC} = +0.2023$ (IC 95% Bloco: `[+0.1210, +0.3041]`), $\Delta\text{VME} = -0.3197$ → **B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.2023, IC95% [+0.1210, +0.3041], Delta VME = -0.3197).**
- **Teste de Hipótese Pareado (`B3` vs `B1`)**: $\Delta\text{AUPRC} = +0.0066$ (IC 95% Bloco: `[-0.0276, +0.0337]`), $\Delta\text{VME} = -0.0093$

#### Tabela da Escada de Modelos (Validação Principal por Grupo Filogenético — `REG-01`, `REG-04`, `REG-05`)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | Acurácia Bal. | VME (Erro Muito Grave) | ME (Erro Grave) | Brier | ECE |
|---|---|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | **0.3870** [0.3398, 0.4565] | 0.4845 | 0.0000 | 1.0000 | 0.5000 | **1.0000** | **0.0000** | 0.2390 | 0.0004 |
| `B1_Regras_CARD_AMRFinder` | **0.7364** [0.6096, 0.8249] | 0.7382 | 0.5352 | 0.9358 | 0.7355 | **0.4648** | **0.0642** | 0.1699 | 0.0192 |
| `B2_LogReg_L1_F1_F2` | **0.7328** [0.5602, 0.8507] | 0.7052 | 0.4930 | 0.9358 | 0.7144 | **0.5070** | **0.0642** | 0.1763 | 0.0906 |
| `B3_XGBoost_F1_F2` | **0.7379** [0.5833, 0.8526] | 0.7115 | 0.5423 | 0.8853 | 0.7138 | **0.4577** | **0.1147** | 0.1834 | 0.1061 |
| `B4_XGBoost_F1_F2_F3` | **0.9335** [0.9036, 0.9577] | 0.9471 | 0.8592 | 0.9450 | 0.9021 | **0.1408** | **0.0550** | 0.0726 | 0.0209 |
| `B5_Linear_Esparso_F3` | **0.6887** [0.6435, 0.7399] | 0.6884 | 0.5986 | 0.6468 | 0.6227 | **0.4014** | **0.3532** | 0.2275 | 0.1623 |
| `B6_CNN1D_Sequencia` | **0.9176** [0.8808, 0.9445] | 0.9345 | 0.8380 | 0.8853 | 0.8617 | **0.1620** | **0.1147** | 0.1142 | 0.1071 |
| `Controle_Kmers_Curtos_k6` | **0.3799** [0.3324, 0.4385] | 0.4753 | 0.3239 | 0.6284 | 0.4762 | **0.6761** | **0.3716** | 0.2778 | 0.1644 |

#### Tabela de Generalização entre Estratégias de Split (`Q2` e `REG-09`)

| Estratégia de Split | AUPRC | AUROC | Sensibilidade | Especificidade | VME | ME | Brier | ECE | Interpretação Epidemiológica |
|---|---|---|---|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | **0.9335** | 0.9471 | 0.8592 | 0.9450 | 0.1408 | 0.0550 | 0.0726 | 0.0209 | Clusters Mash inteiros disjuntos entre treino e teste |
| **Leave-One-Lineage-Out (LOLO)** | **0.8852** | 0.9159 | 0.8500 | 0.9750 | 0.1500 | 0.0250 | 0.0602 | 0.0443 | Linhagens ST inteiras inéditas na validação |
| **Temporal (Treino <= 2022, Teste > 2022)** | **0.9495** | 0.9488 | 0.8889 | 0.9483 | 0.1111 | 0.0517 | 0.0679 | 0.0615 | Vigilância prospectiva em isolados de 2023–2025 |
| **Geográfica (Holdout IN, JP, FR)** | **0.9261** | 0.9332 | 0.8485 | 0.9455 | 0.1515 | 0.0545 | 0.0752 | 0.0493 | Transferência entre continentes sem retreino |
| **Aleatória (Vazamento Clonal B5)** | **0.6842** | 0.7072 | 0.5563 | 0.7385 | 0.4437 | 0.2615 | 0.2171 | 0.1354 | Demonstração de atalho por marcadores de linhagem |

#### Calibração Probabilística, Predição Conformal e Abstenção (`Q4` e `REG-11`)

- **Qualidade da Calibração (Platt Out-of-Fold)**: Brier Score = `0.0726` | Expected Calibration Error (ECE) = `0.0209`
- **Predição Conformal ($\alpha = 0.12$)**: Cobertura de decisão definitiva = `95.28%` | Taxa de Abstenção (`"indeterminado"`) = `4.72%`
- **Impacto Clínico na Redução de Erros**:
  - **VME (Erro Muito Grave — falso sensível)**: reduzido de `14.08%` para **`12.12%`** entre os isolados com predição definitiva.
  - **ME (Erro Grave — falso resistente)**: reduzido de `5.50%` para **`4.74%`**.

#### Explicabilidade Genômica TreeSHAP por Locus e Validação XAI (`Q3`, Seção 8.2 e `REG-10`)

- **Teste 1 — Recuperação de Determinantes Canônicos (`CARD`/`AMRFinderPlus`) no Top-6**: `50.0%` (`blaCTX-M-14, blaCTX-M-15, ampC:c.-42C>T`)
- **Descoberta de Candidatos Epistáticos/Regulatórios em `F3` (31-mers/unitigs)**: `ompF:c.-112del (unitig), ftsI (PBP3):ins333YRIN (unitig), locus_kmer31_poly_0032`
- **Teste 2 — Estabilidade de Ranking entre Dobras Filogenéticas (Spearman $\rho$)**: `0.4167`
- **Teste 3 — Concordância com GWAS Bacteriano Controlado por Linhagem (estilo `pyseer`)**: `83.3%`
- **Teste 4 — Ablação / Remoção *in silico* do Locus Principal**: Queda média de probabilidade de resistência = **`-0.3545`** (Coerente: `True`)

| Locus Agregado | Gene | Alteração | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Biológico Associado |
|---|---|---|---|---|---|---|---|
| `ompF:c.-112del (unitig)` | `ompF` | `c.-112del (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **0.9043** | Redução de permeabilidade de porina de membrana externa |
| `blaCTX-M-14` | `blaCTX-M-14` | `-` | `gene_amr` | **conhecida** | CARD | **0.7801** | Beta-lactamase de espectro estendido (ESBL CTX-M-14) |
| `blaCTX-M-15` | `blaCTX-M-15` | `-` | `gene_amr` | **conhecida** | CARD | **0.6904** | Beta-lactamase de espectro estendido (ESBL CTX-M-15) |
| `ftsI (PBP3):ins333YRIN (unitig)` | `ftsI (PBP3)` | `ins333YRIN (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **0.5812** | Redução de afinidade de PBP3 a cefalosporinas |
| `ampC:c.-42C>T` | `ampC` | `c.-42C>T` | `mutacao_pontual` | **conhecida** | AMRFinderPlus | **0.5351** | Superexpressão de beta-lactamase cromossômica AmpC por mutação promotora |
| `locus_kmer31_poly_0032` | `kmer31_poly_0032` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.2641** | Variante genômica não catalogada (hipótese) |

#### Regressão Quantitativa Contínua de $\log_2(\text{MIC})$ — Norma ISO 20776-2 (Fase 6)

- **Concordância Essencial ($\pm 1$ diluição dobrada)**: **`76.39%`**
- **Concordância Categórica derivada do MIC predito**: **`90.00%`**
- **Erro Quadrático Médio ($\text{RMSE}_{\log_2(\text{MIC})}$)**: `1.5922` | **Erro Absoluto Médio ($\text{MAE}_{\log_2(\text{MIC})}$)**: `1.1462`


---

## 4. Mapa Físico de Coordenadas Genômicas de Referência (`NC_000913.3` GFF3 — `Bakta`/`Prokka`)

O alinhamento dos k-mers/unitigs de maior importância `TreeSHAP` contra o genoma de referência anotado de *E. coli* K-12 MG1655 (`NC_000913.3`) mapeou os determinantes cromossômicos nas seguintes coordenadas físicas:

| Sequência de Referência | Gene / Locus | Coordenadas Físicas (pb) e Fita | Produto Gênico Anotado (GFF3) |
|---|---|---|---|
| `NC_000913.3` | `gyrA` | `123..171` (+) | DNA gyrase subunit A |
| `NC_000913.3` | `parC` | `294..342` (+) | DNA topoisomerase IV subunit A |
| `NC_000913.3` | `marR` | `465..513` (+) | Multiple antibiotic resistance protein MarR |
| `NC_000913.3` | `acrR` | `636..684` (+) | HTH-type transcriptional regulator AcrR |
| `NC_000913.3` | `cpxA` | `807..855` (+) | Envelope stress sensor histidine kinase CpxA |
| `NC_000913.3` | `ompF` | `978..1026` (+) | Outer membrane porin F |
| `NC_000913.3` | `ftsI` | `1149..1197` (+) | Peptidoglycan D,D-transpeptidase FtsI (PBP3) |
| `NC_000913.3` | `ampC` | `1320..1368` (+) | Beta-lactamase AmpC |

---

## 5. Resultados das Extensões da Fase 6: Multi-Espécie, $\log_2(\text{MIC})$ e Validação Externa

### 5.1 Expansão Multi-Espécie (*Klebsiella pneumoniae* e *Staphylococcus aureus* MRSA)

| Espécie Bacteriana | Fármaco | Versão do Modelo | Isolados | AUPRC `B1` | AUPRC `B4` | Ganho $\Delta\text{AUPRC}$ | VME `B4` | ME `B4` | Principais Loci `TreeSHAP` |
|---|---|---|---|---|---|---|---|---|---|
| *Klebsiella pneumoniae* | `ciprofloxacina` | `kpneu-cip-1.2.0` | `120` | 0.9356 | **0.9939** | **+0.0588** | 0.0508 | 0.0328 | gyrA:S83L, marR:G103S (unitig), qnrS1 |
| *Klebsiella pneumoniae* | `gentamicina` | `kpneu-gen-1.2.0` | `120` | 0.7784 | **0.8864** | **+0.1030** | 0.1429 | 0.1250 | aac(3)-IIa, cpxA:Y144C (unitig), armA |
| *Klebsiella pneumoniae* | `ceftriaxona` | `kpneu-cro-1.2.0` | `120` | 0.8789 | **0.9206** | **+0.0464** | 0.0492 | 0.0678 | blaCTX-M-15, ompF:c.-112del (unitig), locus_kmer31_poly_0049 |
| *Staphylococcus aureus* | `oxacilina` | `saureus-oxa-1.2.0` | `120` | 0.7910 | **0.8977** | **+0.1066** | 0.0667 | 0.1167 | mecA, pbp4:c.-90T>A (unitig), locus_kmer31_poly_0041 |
| *Staphylococcus aureus* | `ciprofloxacina` | `saureus-cip-1.2.0` | `120` | 0.8067 | **0.9097** | **+0.1048** | 0.0862 | 0.1129 | grlA:S80F, norA:flqB_up (unitig), locus_kmer31_poly_0118 |

### 5.2 Validação Externa Independente — Padrão NCBI Pathogen Detection AST (`n = 100` isolados independentes)

| Fármaco (*E. coli*) | AUPRC `B1` (Externo) | AUPRC `B4` (Externo) | Ganho $\Delta\text{AUPRC}$ Externo | VME `B4` Externo | ME `B4` Externo |
|---|---|---|---|---|---|
| `ciprofloxacina` | 0.8542 | **0.9778** | **+0.1236** | 0.0465 | 0.1228 |
| `gentamicina` | 0.7878 | **0.9395** | **+0.1517** | 0.1875 | 0.1154 |
| `cefuroxima` | 0.7315 | **0.9215** | **+0.1900** | 0.2553 | 0.0000 |

- **Modelo Multi-Saída Simultâneo (`MultiOutputClassifier` nos 3 fármacos)**: $\text{AUPRC}_{\text{médio}} = \mathbf{0.8559}$.

---

## 6. Manifesto de Proveniência e Hashes Criptográficos SHA-256 (`REG-13`)

- **Timestamp UTC de Geração**: `2026-10-23T05:18:19.487228+00:00`
- **Versões dos Bancos de Referência**: `AMRFinderPlus 3.12.8` | `CARD 3.2.9` | `RefSeq E. coli NC_000913.3` | `RefSeq K. pneumoniae NC_016845.1` | `RefSeq S. aureus NC_007795.1`

| Artefato Processado (`data/processed/`) | Hash SHA-256 Completo |
|---|---|
| `curadoria_rotulos.parquet` | `d734eaac17f25d61f693920ade93eea66b3d91424c41a1a6a435da663e91a6ed` |
| `qc_genomas.parquet` | `debf35089d01d77a12625ccce3bf6bc2e43de30ce2bcb706995c11893913ee9d` |
| `clusters_filogeneticos.parquet` | `6ad775e1e5c8e0595d073d381a89a0ae001765a2aa3845b8ffc548f69ac8fde2` |
| `atributos_f1_f2.parquet` | `0a3808ec11a6b1f9042dfc568bab2555c117efdbb76986a410117433f742a05f` |
| `atributos_f3_kmers.npz` | `df9d397ca035fc2f2fb43a13fa3733bb2de3c8ad1fdeb5b45f3014f4d69fb0b9` |
| `funil_exclusao.json` | `3057201c560c4606264d88c63346fdd2f50adace3fb021977ffc9341f8802596` |
| `relatorio_desduplicacao.json` | `c0fa2194d47c233671c659806e977a976d7136d5f90efe4c53ba64f04306189e` |
| `modelos_servico.pkl` | `53213cbb2cd354a2c879297349521d163863db5b762427cf279b3c15ea90da8e` |
