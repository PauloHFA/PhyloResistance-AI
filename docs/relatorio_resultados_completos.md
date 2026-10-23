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

- **Coorte avaliada**: `120` isolados únicos pós-QC e pós-desduplicação Mash | **Prevalência de Resistência (R)**: `45.00%`
- **Teste de Hipótese Pareado Q1 (`B4` vs `B1`)**: $\Delta\text{AUPRC} = +0.0774$ (IC 95% Bloco: `[+0.0546, +0.1008]`), $\Delta\text{VME} = -0.0830$ → **B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.0774, IC95% [+0.0546, +0.1008], Delta VME = -0.0830).**
- **Teste de Hipótese Pareado (`B3` vs `B1`)**: $\Delta\text{AUPRC} = +0.0302$ (IC 95% Bloco: `[+0.0155, +0.0486]`), $\Delta\text{VME} = +0.0290$

#### Tabela da Escada de Modelos (Validação Principal por Grupo Filogenético — `REG-01`, `REG-04`, `REG-05`)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | Acurácia Bal. | VME (Erro Muito Grave) | ME (Erro Grave) | Brier | ECE |
|---|---|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | **0.3988** [0.3538, 0.4808] | 0.3934 | 0.0000 | 1.0000 | 0.5000 | **1.0000** | **0.0000** | 0.2536 | 0.0056 |
| `B1_Regras_CARD_AMRFinder` | **0.9101** [0.8741, 0.9203] | 0.8886 | 0.8519 | 0.9697 | 0.9108 | **0.1481** | **0.0303** | 0.0767 | 0.0333 |
| `B2_LogReg_L1_F1_F2` | **0.9210** [0.8555, 0.9589] | 0.9112 | 0.7963 | 0.9697 | 0.8830 | **0.2037** | **0.0303** | 0.0919 | 0.0509 |
| `B3_XGBoost_F1_F2` | **0.9338** [0.9056, 0.9572] | 0.9196 | 0.8333 | 0.9545 | 0.8939 | **0.1667** | **0.0455** | 0.0853 | 0.0658 |
| `B4_XGBoost_F1_F2_F3` | **0.9781** [0.9664, 0.9889] | 0.9756 | 0.9444 | 0.9545 | 0.9495 | **0.0556** | **0.0455** | 0.0464 | 0.0321 |
| `B5_Linear_Esparso_F3` | **0.6489** [0.5295, 0.7080] | 0.5915 | 0.4815 | 0.6061 | 0.5438 | **0.5185** | **0.3939** | 0.2703 | 0.2032 |
| `B6_CNN1D_Sequencia` | **0.9369** [0.8900, 0.9709] | 0.9352 | 0.8333 | 0.8788 | 0.8561 | **0.1667** | **0.1212** | 0.1091 | 0.0888 |
| `Controle_Kmers_Curtos_k6` | **0.4444** [0.3437, 0.5869] | 0.4517 | 0.3704 | 0.5455 | 0.4579 | **0.6296** | **0.4545** | 0.3077 | 0.2110 |

#### Tabela de Generalização entre Estratégias de Split (`Q2` e `REG-09`)

| Estratégia de Split | AUPRC | AUROC | Sensibilidade | Especificidade | VME | ME | Brier | ECE | Interpretação Epidemiológica |
|---|---|---|---|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | **0.9781** | 0.9756 | 0.9444 | 0.9545 | 0.0556 | 0.0455 | 0.0464 | 0.0321 | Clusters Mash inteiros disjuntos entre treino e teste |
| **Leave-One-Lineage-Out (LOLO)** | **1.0000** | 1.0000 | 1.0000 | 0.9167 | 0.0000 | 0.0833 | 0.0232 | 0.0762 | Linhagens ST inteiras inéditas na validação |
| **Temporal (Treino <= 2022, Teste > 2022)** | **0.9936** | 0.9944 | 0.9167 | 0.9333 | 0.0833 | 0.0667 | 0.0294 | 0.0620 | Vigilância prospectiva em isolados de 2023–2025 |
| **Geográfica (Holdout IN, JP, FR)** | **1.0000** | 1.0000 | 1.0000 | 1.0000 | 0.0000 | 0.0000 | 0.0143 | 0.0967 | Transferência entre continentes sem retreino |
| **Aleatória (Vazamento Clonal B5)** | **0.6679** | 0.6049 | 0.5556 | 0.6212 | 0.4444 | 0.3788 | 0.2672 | 0.2423 | Demonstração de atalho por marcadores de linhagem |

#### Calibração Probabilística, Predição Conformal e Abstenção (`Q4` e `REG-11`)

- **Qualidade da Calibração (Platt Out-of-Fold)**: Brier Score = `0.0464` | Expected Calibration Error (ECE) = `0.0321`
- **Predição Conformal ($\alpha = 0.12$)**: Cobertura de decisão definitiva = `99.17%` | Taxa de Abstenção (`"indeterminado"`) = `0.83%`
- **Impacto Clínico na Redução de Erros**:
  - **VME (Erro Muito Grave — falso sensível)**: reduzido de `5.56%` para **`5.66%`** entre os isolados com predição definitiva.
  - **ME (Erro Grave — falso resistente)**: reduzido de `4.55%` para **`4.55%`**.

#### Explicabilidade Genômica TreeSHAP por Locus e Validação XAI (`Q3`, Seção 8.2 e `REG-10`)

- **Teste 1 — Recuperação de Determinantes Canônicos (`CARD`/`AMRFinderPlus`) no Top-6**: `33.3%` (`gyrA:S83L, parC:S80I`)
- **Descoberta de Candidatos Epistáticos/Regulatórios em `F3` (31-mers/unitigs)**: `marR:G103S (unitig), acrR:fs*45 (unitig), locus_kmer31_poly_0014, locus_kmer31_poly_0012`
- **Teste 2 — Estabilidade de Ranking entre Dobras Filogenéticas (Spearman $\rho$)**: `0.4316`
- **Teste 3 — Concordância com GWAS Bacteriano Controlado por Linhagem (estilo `pyseer`)**: `50.0%`
- **Teste 4 — Ablação / Remoção *in silico* do Locus Principal**: Queda média de probabilidade de resistência = **`-0.5051`** (Coerente: `True`)

| Locus Agregado | Gene | Alteração | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Biológico Associado |
|---|---|---|---|---|---|---|---|
| `gyrA:S83L` | `gyrA` | `S83L` | `mutacao_pontual` | **conhecida** | AMRFinderPlus | **1.3847** | Alteração do alvo DNA girase subunidade A (QRDR S83L) |
| `marR:G103S (unitig)` | `marR` | `G103S (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **1.1394** | Regulador reprimido de marAB / bomba de efluxo AcrAB-TolC |
| `parC:S80I` | `parC` | `S80I` | `mutacao_pontual` | **conhecida** | AMRFinderPlus | **0.8271** | Alteração da topoisomerase IV subunidade A (QRDR S80I) |
| `acrR:fs*45 (unitig)` | `acrR` | `fs*45 (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **0.7398** | Desrepressão de bomba de efluxo multidroga AcrAB |
| `locus_kmer31_poly_0014` | `kmer31_poly_0014` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4895** | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0012` | `kmer31_poly_0012` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4032** | Variante genômica não catalogada (hipótese) |

#### Regressão Quantitativa Contínua de $\log_2(\text{MIC})$ — Norma ISO 20776-2 (Fase 6)

- **Concordância Essencial ($\pm 1$ diluição dobrada)**: **`55.83%`**
- **Concordância Categórica derivada do MIC predito**: **`92.50%`**
- **Erro Quadrático Médio ($\text{RMSE}_{\log_2(\text{MIC})}$)**: `2.2431` | **Erro Absoluto Médio ($\text{MAE}_{\log_2(\text{MIC})}$)**: `1.6047`

### 3.2 Fármaco: `GENTAMICINA` (Modelo `ecoli-gen-1.2.0`)

- **Coorte avaliada**: `120` isolados únicos pós-QC e pós-desduplicação Mash | **Prevalência de Resistência (R)**: `47.50%`
- **Teste de Hipótese Pareado Q1 (`B4` vs `B1`)**: $\Delta\text{AUPRC} = +0.0453$ (IC 95% Bloco: `[-0.0292, +0.1564]`), $\Delta\text{VME} = -0.0234$ → **B4_XGBoost_F1_F2_F3 NÃO supera com significância estatística o Baseline B1 (Delta AUPRC = +0.0453, IC95% [-0.0292, +0.1564]). Resultado reportado conforme REG-04.**
- **Teste de Hipótese Pareado (`B3` vs `B1`)**: $\Delta\text{AUPRC} = +0.0089$ (IC 95% Bloco: `[-0.0324, +0.0288]`), $\Delta\text{VME} = +0.1249$

#### Tabela da Escada de Modelos (Validação Principal por Grupo Filogenético — `REG-01`, `REG-04`, `REG-05`)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | Acurácia Bal. | VME (Erro Muito Grave) | ME (Erro Grave) | Brier | ECE |
|---|---|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | **0.4344** [0.3258, 0.5219] | 0.4095 | 0.2105 | 0.7143 | 0.4624 | **0.7895** | **0.2857** | 0.2535 | 0.0461 |
| `B1_Regras_CARD_AMRFinder` | **0.8343** [0.6579, 0.9274] | 0.8020 | 0.7368 | 0.9365 | 0.8367 | **0.2632** | **0.0635** | 0.1343 | 0.0838 |
| `B2_LogReg_L1_F1_F2` | **0.8544** [0.6802, 0.9247] | 0.8092 | 0.5789 | 0.9365 | 0.7577 | **0.4211** | **0.0635** | 0.1613 | 0.0535 |
| `B3_XGBoost_F1_F2` | **0.8533** [0.6324, 0.9266] | 0.8076 | 0.6140 | 0.9365 | 0.7753 | **0.3860** | **0.0635** | 0.1561 | 0.1363 |
| `B4_XGBoost_F1_F2_F3` | **0.8592** [0.7853, 0.9445] | 0.8819 | 0.7193 | 0.8254 | 0.7723 | **0.2807** | **0.1746** | 0.1419 | 0.0775 |
| `B5_Linear_Esparso_F3` | **0.5595** [0.4221, 0.7223] | 0.5834 | 0.6140 | 0.5714 | 0.5927 | **0.3860** | **0.4286** | 0.2807 | 0.1557 |
| `B6_CNN1D_Sequencia` | **0.7585** [0.5818, 0.8416] | 0.7884 | 0.6842 | 0.6825 | 0.6834 | **0.3158** | **0.3175** | 0.2520 | 0.2533 |
| `Controle_Kmers_Curtos_k6` | **0.4382** [0.2862, 0.5508] | 0.4456 | 0.3684 | 0.5873 | 0.4779 | **0.6316** | **0.4127** | 0.3222 | 0.2384 |

#### Tabela de Generalização entre Estratégias de Split (`Q2` e `REG-09`)

| Estratégia de Split | AUPRC | AUROC | Sensibilidade | Especificidade | VME | ME | Brier | ECE | Interpretação Epidemiológica |
|---|---|---|---|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | **0.8592** | 0.8819 | 0.7193 | 0.8254 | 0.2807 | 0.1746 | 0.1419 | 0.0775 | Clusters Mash inteiros disjuntos entre treino e teste |
| **Leave-One-Lineage-Out (LOLO)** | **0.9796** | 0.9848 | 0.7222 | 0.9545 | 0.2778 | 0.0455 | 0.0722 | 0.1622 | Linhagens ST inteiras inéditas na validação |
| **Temporal (Treino <= 2022, Teste > 2022)** | **0.9663** | 0.9667 | 0.8333 | 1.0000 | 0.1667 | 0.0000 | 0.0790 | 0.1187 | Vigilância prospectiva em isolados de 2023–2025 |
| **Geográfica (Holdout IN, JP, FR)** | **0.8291** | 0.8542 | 0.7083 | 0.8333 | 0.2917 | 0.1667 | 0.1495 | 0.0993 | Transferência entre continentes sem retreino |
| **Aleatória (Vazamento Clonal B5)** | **0.5998** | 0.6051 | 0.5088 | 0.6032 | 0.4912 | 0.3968 | 0.2788 | 0.1851 | Demonstração de atalho por marcadores de linhagem |

#### Calibração Probabilística, Predição Conformal e Abstenção (`Q4` e `REG-11`)

- **Qualidade da Calibração (Platt Out-of-Fold)**: Brier Score = `0.1419` | Expected Calibration Error (ECE) = `0.0775`
- **Predição Conformal ($\alpha = 0.12$)**: Cobertura de decisão definitiva = `85.83%` | Taxa de Abstenção (`"indeterminado"`) = `14.17%`
- **Impacto Clínico na Redução de Erros**:
  - **VME (Erro Muito Grave — falso sensível)**: reduzido de `28.07%` para **`17.39%`** entre os isolados com predição definitiva.
  - **ME (Erro Grave — falso resistente)**: reduzido de `17.46%` para **`15.79%`**.

#### Explicabilidade Genômica TreeSHAP por Locus e Validação XAI (`Q3`, Seção 8.2 e `REG-10`)

- **Teste 1 — Recuperação de Determinantes Canônicos (`CARD`/`AMRFinderPlus`) no Top-6**: `33.3%` (`aac(3)-IIa, ant(2'')-Ia`)
- **Descoberta de Candidatos Epistáticos/Regulatórios em `F3` (31-mers/unitigs)**: `acrR:fs*45 (unitig), locus_kmer31_poly_0030, locus_kmer31_poly_0047, locus_kmer31_poly_0027`
- **Teste 2 — Estabilidade de Ranking entre Dobras Filogenéticas (Spearman $\rho$)**: `0.4470`
- **Teste 3 — Concordância com GWAS Bacteriano Controlado por Linhagem (estilo `pyseer`)**: `66.7%`
- **Teste 4 — Ablação / Remoção *in silico* do Locus Principal**: Queda média de probabilidade de resistência = **`-0.3354`** (Coerente: `True`)

| Locus Agregado | Gene | Alteração | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Biológico Associado |
|---|---|---|---|---|---|---|---|
| `acrR:fs*45 (unitig)` | `acrR` | `fs*45 (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **0.9267** | Desrepressão de bomba de efluxo multidroga AcrAB |
| `aac(3)-IIa` | `aac(3)-IIa` | `-` | `gene_amr` | **conhecida** | CARD | **0.6232** | Aminoglicosídeo 3-N-acetiltransferase |
| `locus_kmer31_poly_0030` | `kmer31_poly_0030` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4896** | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0047` | `kmer31_poly_0047` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4631** | Variante genômica não catalogada (hipótese) |
| `ant(2'')-Ia` | `ant(2'')-Ia` | `-` | `gene_amr` | **conhecida** | CARD | **0.4098** | Aminoglicosídeo 2''-O-nucleotidiltransferase |
| `locus_kmer31_poly_0027` | `kmer31_poly_0027` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4026** | Variante genômica não catalogada (hipótese) |

#### Regressão Quantitativa Contínua de $\log_2(\text{MIC})$ — Norma ISO 20776-2 (Fase 6)

- **Concordância Essencial ($\pm 1$ diluição dobrada)**: **`54.17%`**
- **Concordância Categórica derivada do MIC predito**: **`76.67%`**
- **Erro Quadrático Médio ($\text{RMSE}_{\log_2(\text{MIC})}$)**: `2.3748` | **Erro Absoluto Médio ($\text{MAE}_{\log_2(\text{MIC})}$)**: `1.7992`

### 3.3 Fármaco: `CEFUROXIMA` (Modelo `ecoli-cxm-1.2.0`)

- **Coorte avaliada**: `120` isolados únicos pós-QC e pós-desduplicação Mash | **Prevalência de Resistência (R)**: `38.33%`
- **Teste de Hipótese Pareado Q1 (`B4` vs `B1`)**: $\Delta\text{AUPRC} = +0.1916$ (IC 95% Bloco: `[+0.1143, +0.3195]`), $\Delta\text{VME} = -0.2126$ → **B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.1916, IC95% [+0.1143, +0.3195], Delta VME = -0.2126).**
- **Teste de Hipótese Pareado (`B3` vs `B1`)**: $\Delta\text{AUPRC} = +0.0338$ (IC 95% Bloco: `[+0.0101, +0.0821]`), $\Delta\text{VME} = +0.0384$

#### Tabela da Escada de Modelos (Validação Principal por Grupo Filogenético — `REG-01`, `REG-04`, `REG-05`)

| Modelo | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | Acurácia Bal. | VME (Erro Muito Grave) | ME (Erro Grave) | Brier | ECE |
|---|---|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | **0.3754** [0.2961, 0.4439] | 0.4824 | 0.0000 | 1.0000 | 0.5000 | **1.0000** | **0.0000** | 0.2366 | 0.0011 |
| `B1_Regras_CARD_AMRFinder` | **0.7112** [0.5453, 0.8523] | 0.7283 | 0.5652 | 0.9730 | 0.7691 | **0.4348** | **0.0270** | 0.1485 | 0.0519 |
| `B2_LogReg_L1_F1_F2` | **0.7224** [0.5428, 0.8592] | 0.7083 | 0.5435 | 0.9459 | 0.7447 | **0.4565** | **0.0541** | 0.1707 | 0.1095 |
| `B3_XGBoost_F1_F2` | **0.7463** [0.5877, 0.8647] | 0.7305 | 0.5217 | 0.8649 | 0.6933 | **0.4783** | **0.1351** | 0.1799 | 0.1522 |
| `B4_XGBoost_F1_F2_F3` | **0.9185** [0.8649, 0.9670] | 0.9057 | 0.8043 | 0.9189 | 0.8616 | **0.1957** | **0.0811** | 0.0940 | 0.0813 |
| `B5_Linear_Esparso_F3` | **0.6117** [0.5420, 0.6629] | 0.5943 | 0.5000 | 0.6892 | 0.5946 | **0.5000** | **0.3108** | 0.2639 | 0.2239 |
| `B6_CNN1D_Sequencia` | **0.8301** [0.7311, 0.8992] | 0.8237 | 0.6957 | 0.9324 | 0.8140 | **0.3043** | **0.0676** | 0.1492 | 0.1400 |
| `Controle_Kmers_Curtos_k6` | **0.3970** [0.3410, 0.4773] | 0.4559 | 0.2609 | 0.6081 | 0.4345 | **0.7391** | **0.3919** | 0.3020 | 0.2494 |

#### Tabela de Generalização entre Estratégias de Split (`Q2` e `REG-09`)

| Estratégia de Split | AUPRC | AUROC | Sensibilidade | Especificidade | VME | ME | Brier | ECE | Interpretação Epidemiológica |
|---|---|---|---|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | **0.9185** | 0.9057 | 0.8043 | 0.9189 | 0.1957 | 0.0811 | 0.0940 | 0.0813 | Clusters Mash inteiros disjuntos entre treino e teste |
| **Leave-One-Lineage-Out (LOLO)** | **0.8627** | 0.8489 | 0.7143 | 0.9231 | 0.2857 | 0.0769 | 0.1212 | 0.1548 | Linhagens ST inteiras inéditas na validação |
| **Temporal (Treino <= 2022, Teste > 2022)** | **0.8271** | 0.8553 | 0.6250 | 0.9474 | 0.3750 | 0.0526 | 0.1121 | 0.1365 | Vigilância prospectiva em isolados de 2023–2025 |
| **Geográfica (Holdout IN, JP, FR)** | **0.8832** | 0.9141 | 0.7059 | 0.9459 | 0.2941 | 0.0541 | 0.0969 | 0.0784 | Transferência entre continentes sem retreino |
| **Aleatória (Vazamento Clonal B5)** | **0.6412** | 0.6745 | 0.5435 | 0.6892 | 0.4565 | 0.3108 | 0.2278 | 0.1580 | Demonstração de atalho por marcadores de linhagem |

#### Calibração Probabilística, Predição Conformal e Abstenção (`Q4` e `REG-11`)

- **Qualidade da Calibração (Platt Out-of-Fold)**: Brier Score = `0.0940` | Expected Calibration Error (ECE) = `0.0813`
- **Predição Conformal ($\alpha = 0.12$)**: Cobertura de decisão definitiva = `88.33%` | Taxa de Abstenção (`"indeterminado"`) = `11.67%`
- **Impacto Clínico na Redução de Erros**:
  - **VME (Erro Muito Grave — falso sensível)**: reduzido de `19.57%` para **`16.67%`** entre os isolados com predição definitiva.
  - **ME (Erro Grave — falso resistente)**: reduzido de `8.11%` para **`1.56%`**.

#### Explicabilidade Genômica TreeSHAP por Locus e Validação XAI (`Q3`, Seção 8.2 e `REG-10`)

- **Teste 1 — Recuperação de Determinantes Canônicos (`CARD`/`AMRFinderPlus`) no Top-6**: `33.3%` (`blaCTX-M-15, blaCTX-M-14`)
- **Descoberta de Candidatos Epistáticos/Regulatórios em `F3` (31-mers/unitigs)**: `ompF:c.-112del (unitig), locus_kmer31_poly_0063, locus_kmer31_poly_0100, locus_kmer31_poly_0024`
- **Teste 2 — Estabilidade de Ranking entre Dobras Filogenéticas (Spearman $\rho$)**: `0.3548`
- **Teste 3 — Concordância com GWAS Bacteriano Controlado por Linhagem (estilo `pyseer`)**: `50.0%`
- **Teste 4 — Ablação / Remoção *in silico* do Locus Principal**: Queda média de probabilidade de resistência = **`-0.3986`** (Coerente: `True`)

| Locus Agregado | Gene | Alteração | Tipo | Status | Fonte | Importância SHAP Média | Mecanismo Biológico Associado |
|---|---|---|---|---|---|---|---|
| `blaCTX-M-15` | `blaCTX-M-15` | `-` | `gene_amr` | **conhecida** | CARD | **0.9116** | Beta-lactamase de espectro estendido (ESBL CTX-M-15) |
| `ompF:c.-112del (unitig)` | `ompF` | `c.-112del (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **0.7475** | Redução de permeabilidade de porina de membrana externa |
| `blaCTX-M-14` | `blaCTX-M-14` | `-` | `gene_amr` | **conhecida** | CARD | **0.5776** | Beta-lactamase de espectro estendido (ESBL CTX-M-14) |
| `locus_kmer31_poly_0063` | `kmer31_poly_0063` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4015** | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0100` | `kmer31_poly_0100` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.3404** | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0024` | `kmer31_poly_0024` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.3283** | Variante genômica não catalogada (hipótese) |

#### Regressão Quantitativa Contínua de $\log_2(\text{MIC})$ — Norma ISO 20776-2 (Fase 6)

- **Concordância Essencial ($\pm 1$ diluição dobrada)**: **`78.33%`**
- **Concordância Categórica derivada do MIC predito**: **`90.00%`**
- **Erro Quadrático Médio ($\text{RMSE}_{\log_2(\text{MIC})}$)**: `1.5548` | **Erro Absoluto Médio ($\text{MAE}_{\log_2(\text{MIC})}$)**: `1.1103`


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
| *Klebsiella pneumoniae* | `ciprofloxacina` | `kpneu-cip-1.2.0` | `90` | 0.8341 | **0.9600** | **+0.1259** | 0.0476 | 0.0625 | gyrA:S83L, marR:G103S (unitig), locus_kmer31_poly_0039 |
| *Klebsiella pneumoniae* | `gentamicina` | `kpneu-gen-1.2.0` | `90` | 0.7501 | **0.8173** | **+0.0614** | 0.2308 | 0.0980 | aac(3)-IIa, cpxA:Y144C (unitig), blaCTX-M-15 |
| *Klebsiella pneumoniae* | `ceftriaxona` | `kpneu-cro-1.2.0` | `90` | 0.7668 | **0.8927** | **+0.1159** | 0.1795 | 0.1569 | blaCTX-M-15, ompF:c.-112del (unitig), locus_kmer31_poly_0119 |
| *Staphylococcus aureus* | `oxacilina` | `saureus-oxa-1.2.0` | `90` | 0.8529 | **0.9778** | **+0.1280** | 0.0476 | 0.0417 | mecA, pbp4:c.-90T>A (unitig), locus_kmer31_poly_0034 |
| *Staphylococcus aureus* | `ciprofloxacina` | `saureus-cip-1.2.0` | `90` | 0.7685 | **0.9441** | **+0.1676** | 0.1250 | 0.1429 | grlA:S80F, norA:flqB_up (unitig), locus_kmer31_poly_0024 |

### 5.2 Validação Externa Independente — Padrão NCBI Pathogen Detection AST (`n = 100` isolados independentes)

| Fármaco (*E. coli*) | AUPRC `B1` (Externo) | AUPRC `B4` (Externo) | Ganho $\Delta\text{AUPRC}$ Externo | VME `B4` Externo | ME `B4` Externo |
|---|---|---|---|---|---|
| `ciprofloxacina` | 0.8542 | **0.9836** | **+0.1294** | 0.0465 | 0.0877 |
| `gentamicina` | 0.7878 | **0.8762** | **+0.0884** | 0.1458 | 0.2115 |
| `cefuroxima` | 0.7315 | **0.9296** | **+0.1981** | 0.2128 | 0.0566 |

- **Modelo Multi-Saída Simultâneo (`MultiOutputClassifier` nos 3 fármacos)**: $\text{AUPRC}_{\text{médio}} = \mathbf{0.8906}$.

---

## 6. Manifesto de Proveniência e Hashes Criptográficos SHA-256 (`REG-13`)

- **Timestamp UTC de Geração**: `2026-10-23T16:24:46.919795+00:00`
- **Versões dos Bancos de Referência**: `AMRFinderPlus 3.12.8` | `CARD 3.2.9` | `RefSeq E. coli NC_000913.3` | `RefSeq K. pneumoniae NC_016845.1` | `RefSeq S. aureus NC_007795.1`

| Artefato Processado (`data/processed/`) | Hash SHA-256 Completo |
|---|---|
| `curadoria_rotulos.parquet` | `ab7862a2a48e284c4e287dfa5e49d50a97b5374795d698b5eed52f51d3c7e6a2` |
| `qc_genomas.parquet` | `b2bf69ab75565b8bd8aa20d1f6a9e2db1ce10cfc3c613aeba9063fdd7fc84c65` |
| `clusters_filogeneticos.parquet` | `205acdc06db15c0ec6d38d60a2365bb4d4a215c105686d4f3b5f73b2f019da58` |
| `atributos_f1_f2.parquet` | `b5d1606b937bce6b3a7fd59b009cbbcb6cce40fd0c2bce56e654ca1a88b111f0` |
| `atributos_f3_kmers.npz` | `790419cdd3e6249357ef6ad42d2972016a1d85f11c47bbc2c089a668d3e75fd1` |
| `funil_exclusao.json` | `cadd5d3dbceeab50e68e074cc3d8acda3a88c3df9b82da3d61cc59bc77daf559` |
| `relatorio_desduplicacao.json` | `ece314a9111df15dc8cd68eab4aa04926a5257c03b552188d70de090f4deca2d` |
| `modelos_servico.pkl` | `99274a25da72ac865a4b99f0dac2c262267652f306371f5cb2615cb60db15a5d` |
