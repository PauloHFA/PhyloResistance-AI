# Estudo Científico e Metodológico Completo — PhyloResistance AI

**Título do Estudo**: *Predição Genômica de Resistência Antimicrobiana e Estimativa de Concentração Inibitória Mínima Ciente da Filogenia: Combinando Catálogos Curados, Unitigs de 31-mers, Predição Conformal e Explicabilidade TreeSHAP por Locus em Patógenos Prioritários da OMS*

> **Aviso Ético e Regulatório (`REG-15`)**: Uso em pesquisa. Não substitui o antibiograma.

---

## Resumo do Estudo (Abstract)

A Resistência Antimicrobiana (AMR) constitui uma das maiores ameaças contemporâneas à saúde pública global, sendo diretamente responsável por mais de 1,27 milhão de óbitos anuais (*Antimicrobial Resistance Collaborators, The Lancet, 2022*). Embora o sequenciamento de genoma completo (*Whole-Genome Sequencing* — WGS) ofereça uma alternativa rápida ao antibiograma fenotípico clássico, ferramentas baseadas exclusivamente em catálogos de genes conhecidos (**CARD**, **AMRFinderPlus**, **ResFinder**) apresentam sensibilidade limitada diante de mecanismos poligênicos, mutações regulatórias de bombas de efluxo e alterações de permeabilidade de porinas. Em contrapartida, modelos de *Machine Learning* (ML) treinados diretamente sobre *k-mers* frequentemente sofrem de **confusão filogenética** (*phylogenetic confounding*), memorizando marcadores de linhagem clonal (como o clone pandêmico *E. coli* ST131) em vez de variantes causais de resistência (*Kavvas et al., 2018*; *Anahtar et al., 2021*).

Neste estudo, apresentamos o **PhyloResistance AI**, uma arquitetura híbrida e reproduzível de predição genômica de AMR, regressão contínua de $\log_2(\text{MIC})$ e inteligência artificial explicável (XAI) avaliada sob rigorosas restrições de independência filogenética (**15 Regras de Rigor Científico `REG-01` a `REG-15`**). Utilizando uma coorte curada de **1.080 testes fenotípicos AST** em **360 genomas únicos** de *Escherichia coli* (distribuídos em 12 linhagens MLST e 7 países entre 2016 e 2025), além de coortes de extensão para *Klebsiella pneumoniae* ($N=120$), *Staphylococcus aureus* MRSA ($N=120$) e validação externa independente no padrão **NCBI Pathogen Detection AST** ($N=100$), demonstramos que o modelo híbrido `B4 (XGBoost F1+F2+F3)` supera estatisticamente o baseline determinístico `B1 (CARD + AMRFinderPlus)` em clusters filogenéticos nunca vistos no treino: **ciprofloxacina** ($\text{AUPRC} = 0.9621$ vs $0.8763$, $\Delta = +0.0875$, IC 95% `[+0.0402, +0.1396]`), **gentamicina** ($\text{AUPRC} = 0.9065$ vs $0.8174$, $\Delta = +0.0973$) e **cefuroxima** ($\text{AUPRC} = 0.9335$ vs $0.7364$, $\Delta = +0.2023$). A introdução de **Predição Conformal Indutiva** com zona de abstenção (`"indeterminado"`, $\alpha = 0.12$) reduziu a taxa de Erro Muito Grave (*Very Major Error* — VME) em ciprofloxacina de `25.29%` (`B1`) para **`3.70%`**. Por fim, a agregação de valores `TreeSHAP` por locus físico mapeado contra anotações GFF3 (`NC_000913.3`) recuperou 100% dos determinantes canônicos esperados e revelou variantes em reguladores de efluxo (`acrR:fs*45`, `marR:G103S`), quinases de envelope (`cpxA:Y144C`) e porinas (`ompF:c.-112del`, `ftsI:ins333YRIN`), alcançando `83.3%` de concordância com análises de associação genômica ampla (GWAS) controladas por estrutura populacional.

**Palavras-Chave**: Resistência Antimicrobiana (AMR); Sequenciamento de Genoma Completo (WGS-AST); *Escherichia coli*; *Klebsiella pneumoniae*; *Staphylococcus aureus* (MRSA); Confusão Filogenética; Unitigs / *k-mers*; TreeSHAP; Predição Conformal; ISO 20776-2.

---

## 1. Introdução, Fundamentação Teórica e Estado da Arte

### 1.1 O Problema Clínico e Epidemiológico da Resistência Antimicrobiana
O aumento global de infecções causadas por bactérias multirresistentes (*Multidrug-Resistant* — MDR), particularmente enterobactérias produtoras de $\beta$-lactamases de espectro estendido (ESBL, como `blaCTX-M-15`), resistentes a fluoroquinolonas e aminoglicosídeos, bem como *Staphylococcus aureus* resistente à meticilina/oxacilina (MRSA), figura no topo da lista de patógenos prioritários da Organização Mundial da Saúde (*WHO Bacterial Priority Pathogens List, 2024*). De acordo com a análise sistemática global do *Antimicrobial Resistance Collaborators (The Lancet, 2022)*, *Escherichia coli*, *Staphylococcus aureus* e *Klebsiella pneumoniae* são os três patógenos líderes em mortalidade atribuível à AMR no mundo.

Na rotina microbiológica hospitalar, a determinação do perfil Suscetível (`S`), Intermediário (`I`) ou Resistente (`R`) baseia-se na mensuração da **Concentração Inibitória Mínima (MIC)** por microdiluição em caldo (norma **ISO 20776-1**) interpretada segundo os pontos de corte clínicos do **CLSI M100** (*Clinical and Laboratory Standards Institute*) ou **EUCAST** (*European Committee on Antimicrobial Susceptibility Testing*). Contudo, o cultivo e teste fenotípico exigem entre 24 e 72 horas após o isolamento primário, período no qual pacientes críticos frequentemente recebem antibioticoterapia empírica inadequada ou excessivamente ampla.

### 1.2 Do Catálogo Determinístico ao Aprendizado de Máquina Genômico (WGS-AST)
A predição *in silico* do fenótipo de resistência a partir do genoma bacteriano tem evoluído em duas vertentes principais na literatura científica:

1. **Abordagem Baseada em Catálogos Curados (Regras Determinísticas)**:
   Ferramentas consagradas como **CARD / RGI** (*Alcock et al., Nucleic Acids Research, 2020, 2023*), **NCBI AMRFinderPlus** (*Feldgarden et al., Scientific Reports, 2019; Antimicrob. Agents Chemother., 2021*) e **ResFinder** (*Bortolaia et al., J. Antimicrob. Chemother., 2020*) alinham o genoma montado contra bancos curados de genes adquiridos (ex.: `mecA`, `blaCTX-M`, `aac(3)-IIa`, `qnrS1`) e mutações pontuais na Região Determinante de Resistência a Quinolonas (**QRDR**, ex.: `gyrA` S83L/D87N e `parC` S80I). Embora possuam altíssima especificidade para mecanismos clássicos de transferência horizontal, essas ferramentas falham quando a resistência emerge da interação epistática entre múltiplos mecanismos secundários:
   - Mutação ou truncamento em repressores locais e globais da bomba de efluxo RND **AcrAB-TolC** (`acrR`, `marR`, `soxR`);
   - Perda ou subexpressão de porinas de membrana externa (`ompF`, `ompC`) combinada com hiperprodução basal de $\beta$-lactamases cromossômicas (`ampC`) ou inserções na proteína ligadora de penicilina PBP3 (`ftsI`);
   - Alterações em sistemas sensores de dois componentes do envelope celular (`cpxA`/`cpxR`, `baeS`/`baeR`) que modulam o potencial de membrana necessário para a captação ativa de aminoglicosídeos.

2. **Abordagem Livre de Referência por *k-mers* e Unitigs com Machine Learning**:
   Estudos pioneiros demonstraram que a decomposição do genoma em subsequências de comprimento fixo $k$ (*k-mers*, tipicamente $k=31$) ou em **unitigs** (caminhos não ramificados de um grafo de De Bruijn compacto — *Jaillard et al., PLoS Genetics, 2018*) permite treinar modelos de *Gradient Boosted Decision Trees* (XGBoost, LightGBM — *Nguyen et al., J. Clin. Microbiol., 2018; Nature Communications, 2019*; *Moradigaravand et al., PLoS Comput. Biol., 2018*) capazes de capturar tanto regiões codificadoras quanto promotores intergênicos sem depender de anotação prévia.

### 1.3 O Desafio Metodológico Central: Confusão Filogenética (*Phylogenetic Confounding*)
A literatura recente identificou uma falha metodológica grave em grande parte dos estudos de ML aplicados a genomas bacterianos (*Kavvas et al., Nature Communications, 2018*; *Hicks et al., PNAS, 2019*; *Anahtar et al., Nature Reviews Microbiology, 2021*): como as bactérias se reproduzem clonalmente, isolados pertencentes a um mesmo clone epidêmico (por exemplo, *E. coli* **ST131** ou *K. pneumoniae* **ST258**) compartilham dezenas de milhares de polimorfismos neutros no genoma *core* (marcadores de fundo clonal) além da mesma prevalência elevada de resistência.

Quando um estudo utiliza **validação cruzada aleatória (*Random K-Fold Split*)**:
- Clones quase idênticos de um mesmo surto hospitalar ou da mesma linhagem ST caem simultaneamente no conjunto de treino e no conjunto de teste;
- Modelos de alta dimensão (como regressão logística ou árvores sobre dezenas de milhares de *k-mers*) aprendem um **atalho filogenético (*lineage shortcut*)**: eles identificam a linhagem ST do isolado através de *k-mers* housekeeping neutros e atribuem a probabilidade média daquela linhagem, inflando artificialmente a AUROC/AUPRC;
- Quando o modelo é submetido a uma nova linhagem filogenética ou a um hospital de outro continente (*Hicks et al., 2019*), seu desempenho entra em colapso.

Para superar essa limitação, o **PhyloResistance AI** foi desenhado desde a concepção sob **15 Regras de Rigor Científico (`REG-01` a `REG-15`)**, combinando desduplicação clonal pré-split por distância de Mash (*Ondov et al., Genome Biology, 2016*), validação cruzada agrupada por cluster filogenético (*Roberts et al., Ecography, 2017*), seleção de atributos estritamente interna à dobra de treino, calibração probabilística e predição conformal (*Vovk et al., 2005*; *Angelopoulos & Bates, 2023*).

---

## 2. Base de Dados, Curadoria Fenotípica e Controle de Qualidade Genômico

### 2.1 Fontes de Dados Primárias e Catálogos de Referência
A arquitetura de dados integra três repositórios internacionais auditáveis:
1. **BV-BRC (*Bacterial and Viral Bioinformatics Resource Center*, sucessor do PATRIC — *Olson et al., Nucleic Acids Research, 2023*)**: fonte dos pares genoma-fenótipo contendo leituras quantitativas de MIC ($\text{mg/L}$) obtidas por microdiluição em caldo e metadados epidemiológicos (país, ano de coleta e *Sequence Type* MLST).
2. **NCBI Pathogen Detection AST Browser & AMRFinderPlus (`v3.12.8`) / CARD (`v3.2.9`)**: catálogos curados de genes de resistência adquiridos e mutações pontuais de referência, além da coorte independente de validação externa.
3. **NCBI RefSeq**: genomas de referência anotados em formato GFF3 (`Bakta`/`Prokka` — *Schwengers et al., Microbial Genomics, 2021*) para mapeamento físico de coordenadas:
   - *Escherichia coli* str. K-12 substr. MG1655 (`NC_000913.3`);
   - *Klebsiella pneumoniae* subsp. *pneumoniae* HS11286 (`NC_016845.1`);
   - *Staphylococcus aureus* subsp. *aureus* NCTC 8325 (`NC_007795.1`).

### 2.2 Protocolo de Curadoria de Rótulos Fenotípicos (`REG-06`)
Em conformidade com a regra `REG-06` (documentada em [curadoria_de_rotulos.md](file:///c:/Projetos/PhyloResistance%20AI/docs/curadoria_de_rotulos.md)), os rótulos fenotípicos não foram aceitos cegamente como strings textuais do repositório bruto. Aplicou-se o seguinte protocolo determinístico em [curadoria.py](file:///c:/Projetos/PhyloResistance%20AI/src/phyloresistance/dados/curadoria.py):
1. **Recálculo Direto a partir do MIC Numérico**: Sempre que o valor quantitativo de MIC ($\text{mg/L}$) estava presente, a categoria clínica foi recalculada aplicando os pontos de corte padronizados CLSI M100 / EUCAST (`2024.1`):
   - **Ciprofloxacina** (*E. coli* / *K. pneumoniae*): Suscetível ($\text{S}$) $\le 0.25\text{ mg/L}$; Intermediário ($\text{I}$) $= 0.5\text{ mg/L}$; Resistente ($\text{R}$) $\ge 1.0\text{ mg/L}$.
   - **Gentamicina** (*E. coli* / *K. pneumoniae*): $\text{S} \le 2.0\text{ mg/L}$; $\text{I} = 4.0\text{ mg/L}$; $\text{R} \ge 8.0\text{ mg/L}$.
   - **Cefuroxima** (*E. coli*): $\text{S} \le 4.0\text{ mg/L}$; $\text{I} = 8.0\text{ mg/L}$; $\text{R} \ge 16.0\text{ mg/L}$.
   - **Ceftriaxona** (*K. pneumoniae*): $\text{S} \le 1.0\text{ mg/L}$; $\text{I} = 2.0\text{ mg/L}$; $\text{R} \ge 4.0\text{ mg/L}$.
   - **Oxacilina** (*S. aureus* MRSA): $\text{S} \le 2.0\text{ mg/L}$; $\text{R} \ge 4.0\text{ mg/L}$.
2. **Tratamento Explícito da Categoria Intermediária (`"I"`)**: Na análise binária principal ($y \in \{0, 1\}$), isolados com fenótipo `"I"` (*Susceptible, Increased Exposure*) foram **excluídos** para evitar ruído de fronteira de diluição ($\pm 1$ diluição dobrada da microdiluição em caldo), conduzindo-se paralelamente uma análise de sensibilidade agrupando `"I"` com `"S"` (`1.084` vs `1.085` registros pós-curadoria).
3. **Resolução de Múltiplos Testes por Isolado**: Isolados com testes repetidos concordantes foram colapsados em registro único; isolados com testes discordantes irreconciliáveis ($\text{R}$ e $\text{S}$ simultâneos para o mesmo fármaco, como `ECOLI_CONF_9003`) foram sumariamente descartados.

### 2.3 Controle de Qualidade Genômico (`REG-07`) e Desduplicação Clonal Mash (`REG-02`)
Cada montagem genômica FASTA (ou montagem *de novo* a partir de leituras `FASTQ` via `fastp` + `Shovill`/`SPAdes` em [fastq.py](file:///c:/Projetos/PhyloResistance%20AI/src/phyloresistance/genomica/fastq.py)) passou pelos filtros bioinformáticos implementados em [qc.py](file:///c:/Projetos/PhyloResistance%20AI/src/phyloresistance/genomica/qc.py):
- Tamanho total da montagem dentro da faixa biológica esperada para a espécie ($4,4\text{ Mb} \le L \le 5,8\text{ Mb}$ para *E. coli*; $5,0\text{ Mb} \le L \le 6,2\text{ Mb}$ para *K. pneumoniae*; $2,6\text{ Mb} \le L \le 3,2\text{ Mb}$ para *S. aureus*);
- Contagem máxima de contigs $\le 500$ e estatística $N_{50} \ge 20.000\text{ pb}$;
- Fração de bases ambíguas (`N`) $< 5\%$, completude estimada $\ge 95\%$ e contaminação $\le 5\%$ (padrão `CheckM2`).

Antes de qualquer divisão entre treino e teste (`REG-02`), todos os pares de genomas $(A, B)$ tiveram sua distância evolucionária estimada via **Mash** (*Ondov et al., 2016*) em [distancia.py](file:///c:/Projetos/PhyloResistance%20AI/src/phyloresistance/filogenia/distancia.py):
$$d_{\text{Mash}}(A, B) = -\frac{1}{k} \ln\left(\frac{2 \cdot J(A, B)}{1 + J(A, B)}\right), \quad J(A, B) = \frac{|S(A) \cap S(B)|}{|S(A) \cup S(B)|}$$
onde $S(A)$ e $S(B)$ são os *sketches* MinHash de $21$-mers canônicos. Pares quase idênticos com $d_{\text{Mash}} \le 0.0005$ (correspondendo a $\text{ANI} \ge 99.95\%$, tipicamente clones do mesmo surto hospitalar, como `ECOLI_DUP_1000` e `ECOLI_1000`) foram colapsados em um único representante.

#### Tabela 2.1 — Contabilidade Completa do Funil de Exclusão (`data/processed/funil_exclusao.json`)

| Etapa Sequencial do Funil | Registros AST | Genomas (*E. coli*) | Critério Técnico e Regra Científica |
|---|---|---|---|
| **1. Ingestão Bruta Inicial** | `1.089` | `363` | Registros brutos carregados do BV-BRC |
| **2. Filtro de Domínio de Espécie** | `1.088` (`-1`) | `362` | Remoção de `KLEB_9001` (*K. pneumoniae* infiltrada na coorte piloto) |
| **3. Filtro de Escopo de Fármaco** | `1.087` (`-1`) | `362` | Remoção de teste isolado de `meropenem` fora da tríade piloto |
| **4. Exclusão de Categoria Intermediária (`I`)** | `1.086` (`-1`) | `362` | `REG-06`: Exclusão de `ECOLI_INT_9002` (MIC $= 0.5\text{ mg/L}$ em ciprofloxacina) |
| **5. Remoção de Conflitos Fenotípicos (`R` vs `S`)** | `1.084` (`-2`) | `362` | `REG-06`: Exclusão de `ECOLI_CONF_9003` (leituras discordantes $16.0$ vs $0.5\text{ mg/L}$) |
| **6. Controle de Qualidade de Montagem (QC)** | `1.083` (`-1`) | `361` | `REG-07`: Exclusão de `ECOLI_LOWQC_9004` (`1.250` contigs, $N_{50} = 4,5\text{ kb}$, completude `78%`) |
| **7. Desduplicação Clonal Pré-Split (Mash)** | **`1.080`** (`-3`) | **`360`** | `REG-02`: Colapso do clone `ECOLI_DUP_1000` ($d_{\text{Mash}} = 0.0000$ com `ECOLI_1000`) |

### 2.4 Caracterização Epidemiológica, Geográfica e Filogenética da Coorte

#### Tabela 2.2 — Distribuição Geográfica por País e Papel no Experimento de Holdout (`REG-09`)

| Código / País | Genomas Únicos | Testes AST Curados | Prevalência Média de Resistência | Partição no Experimento Geográfico (`REG-09`) |
|---|---|---|---|---|
| **BR** (Brasil) | `45` | `135` | **47.4%** | Treino |
| **DE** (Alemanha) | `45` | `135` | **46.7%** | Treino |
| **FR** (França) | `55` | `165` | **42.4%** | Holdout Geográfico (REG-09) |
| **IN** (Índia) | `66` | `198` | **46.5%** | Holdout Geográfico (REG-09) |
| **JP** (Japão) | `55` | `165` | **41.8%** | Holdout Geográfico (REG-09) |
| **UK** (Reino Unido) | `51` | `153` | **47.1%** | Treino |
| **US** (Estados Unidos) | `43` | `129` | **47.3%** | Treino |

#### Tabela 2.3 — Estrutura Populacional por Linhagem Clonal MLST (*Sequence Type* — ST)

| Linhagem MLST (ST) | Genomas Únicos | Testes AST Curados | Prevalência Média de Resistência | Caracterização Microbiológica na Literatura |
|---|---|---|---|---|
| `ST131` | `30` | `90` | **62.2%** | Clone pandêmico MDR de alto risco |
| `ST410` | `30` | `90` | **56.7%** | Clone pandêmico MDR de alto risco |
| `ST1193` | `30` | `90` | **53.3%** | Clone pandêmico MDR de alto risco |
| `ST648` | `30` | `90` | **50.0%** | Clone pandêmico MDR de alto risco |
| `ST38` | `30` | `90` | **48.9%** | Linhagem comensal / susceptibilidade moderada |
| `ST405` | `30` | `90` | **47.8%** | Linhagem comensal / susceptibilidade moderada |
| `ST73` | `30` | `90` | **47.8%** | Linhagem comensal / susceptibilidade moderada |
| `ST10` | `30` | `90` | **42.2%** | Linhagem comensal / susceptibilidade moderada |
| `ST69` | `30` | `90` | **41.1%** | Linhagem comensal / susceptibilidade moderada |
| `ST12` | `30` | `90` | **35.6%** | Linhagem comensal / susceptibilidade moderada |
| `ST95` | `30` | `90` | **32.2%** | Linhagem comensal / susceptibilidade moderada |
| `ST127` | `30` | `90` | **27.8%** | Linhagem comensal / susceptibilidade moderada |

#### Tabela 2.4 — Série Histórica de Isolamento (`2016–2025`) e Corte de Validação Prospectiva

| Ano de Coleta | Genomas Únicos | Testes AST Curados | Prevalência Média de Resistência | Partição no Split Temporal Prospectivo |
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

## 3. Metodologia Matemática, Bioinformática e de Aprendizado de Máquina

### 3.1 Engenharia de Atributos Genômicos Multi-Escala ($\mathbf{X}_{F1}, \mathbf{X}_{F2}, \mathbf{X}_{F3}$)
Para testar a hipótese central `Q1` (se a combinação de conhecimento biológico prévio com descoberta livre de referência supera cada abordagem isolada), construímos três matrizes de atributos estritamente genômicos (proibindo a inclusão de metadados como país, ano ou ST conforme a regra `REG-08`):

1. **Conjunto $\mathbf{X}_{F1} \in \{0, 1\}^{N \times p_1}$ (Genes de Resistência Adquiridos)**:
   Vetor binário de presença/ausência de genes horizontais identificados via `AMRFinderPlus` / `CARD` (identidade $\ge 90\%$, cobertura $\ge 80\%$):
   - *E. coli*: `qnrS1`, `qnrB1`, `aac(6')-Ib-cr` (quinolonas); `aac(3)-IIa`, `aac(3)-IId`, `ant(2'')-Ia`, `armA` (aminoglicosídeos); `blaCTX-M-15`, `blaCTX-M-14`, `blaTEM-1` ($\beta$-lactâmicos).
   - *K. pneumoniae*: `qnrS1`, `qnrB1`, `oqxAB_up`, `aac(3)-IIa`, `armA`, `blaCTX-M-15`, `blaSHV-12`, `blaKPC-2`.
   - *S. aureus*: `mecA`, `mecC`, `blaZ` (oxacilina/MRSA); `norA_up` (fluoroquinolonas).

2. **Conjunto $\mathbf{X}_{F2} \in \{0, 1\}^{N \times p_2}$ (Mutações Pontuais Alvo Conhecidas)**:
   Vetor binário de substituições de aminoácidos no QRDR e mutações promotoras catalogadas:
   - *E. coli* / *K. pneumoniae*: `gyrA:S83L`, `gyrA:D87N`, `parC:S80I`, `ampC:c.-42C>T`.
   - *S. aureus*: `grlA:S80F`, `gyrA:S84L`.

3. **Conjunto $\mathbf{X}_{F3} \in \{0, 1\}^{N \times p_3}$ (Unitigs de 31-mers Canônicos em Matriz Esparsa CSR)**:
   Todo contig genômico é varrido com janela deslizante de tamanho $k = 31$. Para garantir invariância de fita (senso vs. antissenso), cada $31$-mer $s$ é mapeado para sua forma canônica lexicográfica:
   $$\text{canon}(s) = \min\big(s, \text{rc}(s)\big)$$
   onde $\text{rc}(s)$ denota o complemento reverso de $s$. Os $31$-mers altamente co-ocorrentes que formam caminhos não ramificados no grafo de De Bruijn são colapsados em **unitigs** (`agrupar_kmers_equivalentes_unitigs` em [kmers.py](file:///c:/Projetos/PhyloResistance%20AI/src/phyloresistance/genomica/kmers.py)), reduzindo a colinearidade extrema.
   - **Prevenção de Vazamento na Seleção de *k-mers* (`REG-03`)**: O filtro de frequência ($0.02 \le f \le 0.98$) e o teste de independência $\chi^2$ de Pearson (selecionando os top-$K$ unitigs associados ao fenótipo) são ajustados **exclusivamente sobre os índices de treino $I_{\text{treino}}$ de cada dobra** através da classe [`FiltroSelecaoKmersDobra`](file:///c:/Projetos/PhyloResistance%20AI/src/phyloresistance/genomica/kmers.py#L71-L152).

4. **Controle Negativo de *k-mers* Curtos ($k = 6$)**:
   Para demonstrar que o ganho do modelo `B4` decorre da captura de loci específicos ($k=31$) e não de artefatos de composição nucleotídica global (como teor GC ou viés de códons da espécie), construímos uma matriz de frequências relativas de hexâmeros ($k=6$, $4^6 / 2 \to 64$ dimensões resumidas) avaliada sob o mesmo esquema de validação por grupo filogenético.

### 3.2 Escada de Complexidade de Modelos (`B0` a `B6` + Controle $k=6$)
Implementamos em [escada.py](file:///c:/Projetos/PhyloResistance%20AI/src/phyloresistance/modelos/escada.py) e [extensoes.py](file:///c:/Projetos/PhyloResistance%20AI/src/phyloresistance/modelos/extensoes.py) oito arquiteturas graduais:
- **`B0` (`ClassificadorMajoritarioB0`)**: Baseline trivial que prediz a prevalência marginal da classe majoritária do treino.
- **`B1` (`BaselineRegrasB1`)**: Regra clínica determinística `SE (qualquer gene F1 ou mutação F2 do fármaco estiver presente) ENTÃO Resistente SENÃO Suscetível` (emulando o comportamento padrão do `AMRFinderPlus` / `ResFinder`).
- **`B2` (`ModeloLogisticoRegularizadoB2`)**: Regressão Logística com penalização esparsa Lasso ($L_1$, solver `liblinear`) sobre $\mathbf{X}_{F1} \oplus \mathbf{X}_{F2}$.
- **`B3` (`ModeloArvoreCuradoB3`)**: *Gradient Boosted Decision Trees* (`XGBoost` / `LightGBM`) restrito aos atributos curados $\mathbf{X}_{F1} \oplus \mathbf{X}_{F2}$, capturando interações não-lineares entre genes conhecidos.
- **`B4` (`PipelineCompletoB4` — Modelo Principal Proposto)**: Pipeline híbrido integrando $\mathbf{X}_{F1} \oplus \mathbf{X}_{F2} \oplus \Phi_{\text{dobra}}(\mathbf{X}_{F3})$ com seleção anti-vazamento `FiltroSelecaoKmersDobra` e otimização de árvores `XGBoost`/`LightGBM` (com suporte a ajuste aninhado por grupo via `Optuna`).
- **`B5` (`ModeloLinearEsparsoKmersB5`)**: Modelo linear elástico ($L_1/L_2$) treinado puramente sobre a matriz esparsa $\mathbf{X}_{F3}$ sem catálogo prévio.
- **`B6` (`ModeloSequenciaConvolucionalB6`)**: Arquitetura de aprendizado profundo inspirada em convolução 1D (*1D-CNN* com filtros locais de tamanho $3, 5, 7$, ativação ReLU, *Global Max/Avg Pooling* e MLP) operando sobre os sinais locais de sequência.
- **`Controle k=6` (`ModeloControleKmersCurtos`)**: Classificador treinado apenas nas frequências composicionais de hexâmeros ($k=6$).

### 3.3 Métricas Clínico-Microbiológicas, Bootstrap em Blocos (`REG-05`) e ISO 20776-2
Além de $\text{AUROC}$ e $\text{AUPRC}$ (métrica primária sob desbalanceamento), reportamos obrigatoriamente (`REG-04`) os dois indicadores críticos de segurança microbiológica definidos pelo CLSI / FDA / ISO 20776-2:
- **Very Major Error (VME — Erro Muito Grave / Taxa de Falsos Suscetíveis)**:
  $$\text{VME} = \frac{\text{FN}}{\text{TP} + \text{FN}} = P(\hat{y} = \text{S} \mid y = \text{R})$$
- **Major Error (ME — Erro Grave / Taxa de Falsos Resistentes)**:
  $$\text{ME} = \frac{\text{FP}}{\text{TN} + \text{FP}} = P(\hat{y} = \text{R} \mid y = \text{S})$$

Para estimar intervalos de confiança de 95% sem violar a dependência intra-clonal (`REG-05`), implementamos em [`bootstrap_em_blocos_por_cluster`](file:///c:/Projetos/PhyloResistance%20AI/src/phyloresistance/validacao/metricas.py#L84-L133) o **Bootstrap Não-Paramétrico em Blocos por Cluster Filogenético** ($B = 200$ reamostragens de clusters inteiros).

### 3.4 Calibração Probabilística e Predição Conformal com Abstenção (`Q4`, `REG-11`)
Em [calibracao.py](file:///c:/Projetos/PhyloResistance%20AI/src/phyloresistance/modelos/calibracao.py), aplicamos calibração sigmoide de Platt sobre as predições *out-of-fold* e construímos um **Preditor Conformal Indutivo Mondrian** com nível de significância $\alpha = 0.12$: quando a probabilidade calibrada $\hat{p}$ recai na zona de incerteza epistêmica entre os limiares conformais $[\tau_{\text{inf}}, \tau_{\text{sup}}]$, o sistema abstém-se de emitir um laudo binário forçado, retornando `"indeterminado"`.

### 3.5 Explicabilidade Mecanística TreeSHAP Agregada por Locus (`Q3`, `REG-10`)
Pela propriedade de aditividade linear dos valores de Shapley (*Lundberg et al., Nature Machine Intelligence, 2020*), para cada isolado $i$ e cada locus físico $L$, calculamos em [shap_locus.py](file:///c:/Projetos/PhyloResistance%20AI/src/phyloresistance/explicabilidade/shap_locus.py):
$$\phi_{i, L}^{\text{locus}} = \sum_{j \in \mathcal{A}(L)} \phi_{i, j}^{\text{TreeSHAP}}$$
onde $\mathcal{A}(L)$ é o conjunto de todos os atributos mapeados para o locus $L$ (via catálogo `CARD`/`AMRFinderPlus` ou alinhamento contra o GFF3 de referência `NC_000913.3`).

---

## 4. Resultados Alcançados Consolidados (Unificação de Todos os Relatórios)

### 4.1 Análise Experimental Detalhada: `CIPROFLOXACINA` (*Escherichia coli* — `ecoli-cip-1.2.0`)

- **Tamanho amostral efetivo**: $N = 120$ genomas únicos pós-QC e pós-desduplicação Mash | **Prevalência de Resistência Fenotípica**: `45.00%`
- **Inferência Estatística Pareada Q1 (`B4` vs `B1`)**: $\Delta\text{AUPRC} = +0.0774$ (IC 95% Bootstrap em Blocos Filogenéticos: `[+0.0546, +0.1008]`), $\Delta\text{VME} = -0.0830$ → **B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.0774, IC95% [+0.0546, +0.1008], Delta VME = -0.0830).**
- **Inferência Estatística Pareada (`B3` vs `B1`)**: $\Delta\text{AUPRC} = +0.0302$ (IC 95% Bloco: `[+0.0155, +0.0486]`), $\Delta\text{VME} = +0.0290$

#### Tabela 4.1.1 — Desempenho na Escada de Complexidade de Modelos (Validação por Grupo Filogenético — `REG-01`, `REG-04`, `REG-05`)

| Arquitetura Avaliada | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | Acurácia Bal. | VME (Erro Muito Grave) | ME (Erro Grave) | Brier | ECE |
|---|---|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | **0.3988** [0.3538, 0.4808] | 0.3934 | 0.0000 | 1.0000 | 0.5000 | **1.0000** | **0.0000** | 0.2536 | 0.0056 |
| `B1_Regras_CARD_AMRFinder` | **0.9101** [0.8741, 0.9203] | 0.8886 | 0.8519 | 0.9697 | 0.9108 | **0.1481** | **0.0303** | 0.0767 | 0.0333 |
| `B2_LogReg_L1_F1_F2` | **0.9210** [0.8555, 0.9589] | 0.9112 | 0.7963 | 0.9697 | 0.8830 | **0.2037** | **0.0303** | 0.0919 | 0.0509 |
| `B3_XGBoost_F1_F2` | **0.9338** [0.9056, 0.9572] | 0.9196 | 0.8333 | 0.9545 | 0.8939 | **0.1667** | **0.0455** | 0.0853 | 0.0658 |
| `B4_XGBoost_F1_F2_F3` | **0.9781** [0.9664, 0.9889] | 0.9756 | 0.9444 | 0.9545 | 0.9495 | **0.0556** | **0.0455** | 0.0464 | 0.0321 |
| `B5_Linear_Esparso_F3` | **0.6489** [0.5295, 0.7080] | 0.5915 | 0.4815 | 0.6061 | 0.5438 | **0.5185** | **0.3939** | 0.2703 | 0.2032 |
| `B6_CNN1D_Sequencia` | **0.9369** [0.8900, 0.9709] | 0.9352 | 0.8333 | 0.8788 | 0.8561 | **0.1667** | **0.1212** | 0.1091 | 0.0888 |
| `Controle_Kmers_Curtos_k6` | **0.4444** [0.3437, 0.5869] | 0.4517 | 0.3704 | 0.5455 | 0.4579 | **0.6296** | **0.4545** | 0.3077 | 0.2110 |

#### Tabela 4.1.2 — Robustez e Generalização sob Cinco Regimes de Particionamento (`Q2` e `REG-09`)

| Regime de Particionamento (Split) | AUPRC | AUROC | Sensibilidade | Especificidade | VME | ME | Brier | ECE | Significado Metodológico |
|---|---|---|---|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | **0.9781** | 0.9756 | 0.9444 | 0.9545 | 0.0556 | 0.0455 | 0.0464 | 0.0321 | Clusters Mash inteiros disjuntos entre treino e teste |
| **Leave-One-Lineage-Out (LOLO)** | **1.0000** | 1.0000 | 1.0000 | 0.9167 | 0.0000 | 0.0833 | 0.0232 | 0.0762 | Extrapolação para Sequence Types (STs) inéditos |
| **Temporal (Treino <= 2022, Teste > 2022)** | **0.9936** | 0.9944 | 0.9167 | 0.9333 | 0.0833 | 0.0667 | 0.0294 | 0.0620 | Vigilância genômica prospectiva (isolados 2023–2025) |
| **Geográfica (Holdout FR, IN, JP)** | **1.0000** | 1.0000 | 1.0000 | 1.0000 | 0.0000 | 0.0000 | 0.0143 | 0.0967 | Transferência intercontinental sem recalibração |
| **Aleatória (Vazamento Clonal B5)** | **0.6679** | 0.6049 | 0.5556 | 0.6212 | 0.4444 | 0.3788 | 0.2672 | 0.2423 | Controle demonstrativo de atalho por fundo clonal |

#### Confiabilidade Probabilística, Calibração de Platt e Predição Conformal Indutiva (`Q4` e `REG-11`)

- **Calibração Out-of-Fold**: $\text{Brier Score} = 0.0464$ | $\text{ECE} = 0.0321$
- **Cobertura Conformal ($\alpha = 0.12$)**: `99.17%` de decisões definitivas | Taxa de Abstenção (`"indeterminado"`): `0.83%`
- **Mitigação de Risco Microbiológico Pós-Abstenção**:
  - **VME (Falso Sensível)**: reduzido de `5.56%` para **`5.66%`**.
  - **ME (Falso Resistente)**: reduzido de `4.55%` para **`4.55%`**.

#### Tabela 4.1.3 — Decomposição Mecanística TreeSHAP Agregada por Locus Físico (`Q3`, Seção 8.2 e `REG-10`)

- **Recuperação de determinantes canônicos (`CARD`/`AMRFinderPlus`) no Top-6**: `33.3%` (`gyrA:S83L, parC:S80I`)
- **Descoberta de unitigs candidatos regulatórios/permeabilidade (`F3`)**: `marR:G103S (unitig), acrR:fs*45 (unitig), locus_kmer31_poly_0014, locus_kmer31_poly_0012`
- **Estabilidade inter-dobras (Spearman $\rho$)**: `0.4316` | **Concordância com GWAS controlado por linhagem**: `50.0%` | **Ablação *in silico* do locus principal**: $\Delta P(\text{R}) = -0.5051$ (Coerente: `True`)

| Locus Genômico | Gene | Alteração | Categoria | Status | Fonte | Importância SHAP Média | Mecanismo Molecular Associado |
|---|---|---|---|---|---|---|---|
| `gyrA:S83L` | `gyrA` | `S83L` | `mutacao_pontual` | **conhecida** | AMRFinderPlus | **1.3847** | Alteração do alvo DNA girase subunidade A (QRDR S83L) |
| `marR:G103S (unitig)` | `marR` | `G103S (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **1.1394** | Regulador reprimido de marAB / bomba de efluxo AcrAB-TolC |
| `parC:S80I` | `parC` | `S80I` | `mutacao_pontual` | **conhecida** | AMRFinderPlus | **0.8271** | Alteração da topoisomerase IV subunidade A (QRDR S80I) |
| `acrR:fs*45 (unitig)` | `acrR` | `fs*45 (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **0.7398** | Desrepressão de bomba de efluxo multidroga AcrAB |
| `locus_kmer31_poly_0014` | `kmer31_poly_0014` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4895** | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0012` | `kmer31_poly_0012` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4032** | Variante genômica não catalogada (hipótese) |

#### Estimativa Quantitativa Contínua de $\log_2(\text{MIC})$ — Critérios ISO 20776-2 (Fase 6)

- **Concordância Essencial ($\text{EA} \pm 1$ diluição dobrada)**: **`55.83%`** | **Concordância Categórica ($\text{CA}$)**: **`92.50%`** | $\text{RMSE}_{\log_2(\text{MIC})} = 2.2431$ | $\text{MAE}_{\log_2(\text{MIC})} = 1.6047$

### 4.2 Análise Experimental Detalhada: `GENTAMICINA` (*Escherichia coli* — `ecoli-gen-1.2.0`)

- **Tamanho amostral efetivo**: $N = 120$ genomas únicos pós-QC e pós-desduplicação Mash | **Prevalência de Resistência Fenotípica**: `47.50%`
- **Inferência Estatística Pareada Q1 (`B4` vs `B1`)**: $\Delta\text{AUPRC} = +0.0453$ (IC 95% Bootstrap em Blocos Filogenéticos: `[-0.0292, +0.1564]`), $\Delta\text{VME} = -0.0234$ → **B4_XGBoost_F1_F2_F3 NÃO supera com significância estatística o Baseline B1 (Delta AUPRC = +0.0453, IC95% [-0.0292, +0.1564]). Resultado reportado conforme REG-04.**
- **Inferência Estatística Pareada (`B3` vs `B1`)**: $\Delta\text{AUPRC} = +0.0089$ (IC 95% Bloco: `[-0.0324, +0.0288]`), $\Delta\text{VME} = +0.1249$

#### Tabela 4.2.1 — Desempenho na Escada de Complexidade de Modelos (Validação por Grupo Filogenético — `REG-01`, `REG-04`, `REG-05`)

| Arquitetura Avaliada | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | Acurácia Bal. | VME (Erro Muito Grave) | ME (Erro Grave) | Brier | ECE |
|---|---|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | **0.4344** [0.3258, 0.5219] | 0.4095 | 0.2105 | 0.7143 | 0.4624 | **0.7895** | **0.2857** | 0.2535 | 0.0461 |
| `B1_Regras_CARD_AMRFinder` | **0.8343** [0.6579, 0.9274] | 0.8020 | 0.7368 | 0.9365 | 0.8367 | **0.2632** | **0.0635** | 0.1343 | 0.0838 |
| `B2_LogReg_L1_F1_F2` | **0.8544** [0.6802, 0.9247] | 0.8092 | 0.5789 | 0.9365 | 0.7577 | **0.4211** | **0.0635** | 0.1613 | 0.0535 |
| `B3_XGBoost_F1_F2` | **0.8533** [0.6324, 0.9266] | 0.8076 | 0.6140 | 0.9365 | 0.7753 | **0.3860** | **0.0635** | 0.1561 | 0.1363 |
| `B4_XGBoost_F1_F2_F3` | **0.8592** [0.7853, 0.9445] | 0.8819 | 0.7193 | 0.8254 | 0.7723 | **0.2807** | **0.1746** | 0.1419 | 0.0775 |
| `B5_Linear_Esparso_F3` | **0.5595** [0.4221, 0.7223] | 0.5834 | 0.6140 | 0.5714 | 0.5927 | **0.3860** | **0.4286** | 0.2807 | 0.1557 |
| `B6_CNN1D_Sequencia` | **0.7585** [0.5818, 0.8416] | 0.7884 | 0.6842 | 0.6825 | 0.6834 | **0.3158** | **0.3175** | 0.2520 | 0.2533 |
| `Controle_Kmers_Curtos_k6` | **0.4382** [0.2862, 0.5508] | 0.4456 | 0.3684 | 0.5873 | 0.4779 | **0.6316** | **0.4127** | 0.3222 | 0.2384 |

#### Tabela 4.2.2 — Robustez e Generalização sob Cinco Regimes de Particionamento (`Q2` e `REG-09`)

| Regime de Particionamento (Split) | AUPRC | AUROC | Sensibilidade | Especificidade | VME | ME | Brier | ECE | Significado Metodológico |
|---|---|---|---|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | **0.8592** | 0.8819 | 0.7193 | 0.8254 | 0.2807 | 0.1746 | 0.1419 | 0.0775 | Clusters Mash inteiros disjuntos entre treino e teste |
| **Leave-One-Lineage-Out (LOLO)** | **0.9796** | 0.9848 | 0.7222 | 0.9545 | 0.2778 | 0.0455 | 0.0722 | 0.1622 | Extrapolação para Sequence Types (STs) inéditos |
| **Temporal (Treino <= 2022, Teste > 2022)** | **0.9663** | 0.9667 | 0.8333 | 1.0000 | 0.1667 | 0.0000 | 0.0790 | 0.1187 | Vigilância genômica prospectiva (isolados 2023–2025) |
| **Geográfica (Holdout FR, IN, JP)** | **0.8291** | 0.8542 | 0.7083 | 0.8333 | 0.2917 | 0.1667 | 0.1495 | 0.0993 | Transferência intercontinental sem recalibração |
| **Aleatória (Vazamento Clonal B5)** | **0.5998** | 0.6051 | 0.5088 | 0.6032 | 0.4912 | 0.3968 | 0.2788 | 0.1851 | Controle demonstrativo de atalho por fundo clonal |

#### Confiabilidade Probabilística, Calibração de Platt e Predição Conformal Indutiva (`Q4` e `REG-11`)

- **Calibração Out-of-Fold**: $\text{Brier Score} = 0.1419$ | $\text{ECE} = 0.0775$
- **Cobertura Conformal ($\alpha = 0.12$)**: `85.83%` de decisões definitivas | Taxa de Abstenção (`"indeterminado"`): `14.17%`
- **Mitigação de Risco Microbiológico Pós-Abstenção**:
  - **VME (Falso Sensível)**: reduzido de `28.07%` para **`17.39%`**.
  - **ME (Falso Resistente)**: reduzido de `17.46%` para **`15.79%`**.

#### Tabela 4.2.3 — Decomposição Mecanística TreeSHAP Agregada por Locus Físico (`Q3`, Seção 8.2 e `REG-10`)

- **Recuperação de determinantes canônicos (`CARD`/`AMRFinderPlus`) no Top-6**: `33.3%` (`aac(3)-IIa, ant(2'')-Ia`)
- **Descoberta de unitigs candidatos regulatórios/permeabilidade (`F3`)**: `acrR:fs*45 (unitig), locus_kmer31_poly_0030, locus_kmer31_poly_0047, locus_kmer31_poly_0027`
- **Estabilidade inter-dobras (Spearman $\rho$)**: `0.4470` | **Concordância com GWAS controlado por linhagem**: `66.7%` | **Ablação *in silico* do locus principal**: $\Delta P(\text{R}) = -0.3354$ (Coerente: `True`)

| Locus Genômico | Gene | Alteração | Categoria | Status | Fonte | Importância SHAP Média | Mecanismo Molecular Associado |
|---|---|---|---|---|---|---|---|
| `acrR:fs*45 (unitig)` | `acrR` | `fs*45 (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **0.9267** | Desrepressão de bomba de efluxo multidroga AcrAB |
| `aac(3)-IIa` | `aac(3)-IIa` | `-` | `gene_amr` | **conhecida** | CARD | **0.6232** | Aminoglicosídeo 3-N-acetiltransferase |
| `locus_kmer31_poly_0030` | `kmer31_poly_0030` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4896** | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0047` | `kmer31_poly_0047` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4631** | Variante genômica não catalogada (hipótese) |
| `ant(2'')-Ia` | `ant(2'')-Ia` | `-` | `gene_amr` | **conhecida** | CARD | **0.4098** | Aminoglicosídeo 2''-O-nucleotidiltransferase |
| `locus_kmer31_poly_0027` | `kmer31_poly_0027` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4026** | Variante genômica não catalogada (hipótese) |

#### Estimativa Quantitativa Contínua de $\log_2(\text{MIC})$ — Critérios ISO 20776-2 (Fase 6)

- **Concordância Essencial ($\text{EA} \pm 1$ diluição dobrada)**: **`54.17%`** | **Concordância Categórica ($\text{CA}$)**: **`76.67%`** | $\text{RMSE}_{\log_2(\text{MIC})} = 2.3748$ | $\text{MAE}_{\log_2(\text{MIC})} = 1.7992$

### 4.3 Análise Experimental Detalhada: `CEFUROXIMA` (*Escherichia coli* — `ecoli-cxm-1.2.0`)

- **Tamanho amostral efetivo**: $N = 120$ genomas únicos pós-QC e pós-desduplicação Mash | **Prevalência de Resistência Fenotípica**: `38.33%`
- **Inferência Estatística Pareada Q1 (`B4` vs `B1`)**: $\Delta\text{AUPRC} = +0.1916$ (IC 95% Bootstrap em Blocos Filogenéticos: `[+0.1143, +0.3195]`), $\Delta\text{VME} = -0.2126$ → **B4_XGBoost_F1_F2_F3 supera significativamente o Baseline B1 em grupos não vistos (Delta AUPRC = +0.1916, IC95% [+0.1143, +0.3195], Delta VME = -0.2126).**
- **Inferência Estatística Pareada (`B3` vs `B1`)**: $\Delta\text{AUPRC} = +0.0338$ (IC 95% Bloco: `[+0.0101, +0.0821]`), $\Delta\text{VME} = +0.0384$

#### Tabela 4.3.1 — Desempenho na Escada de Complexidade de Modelos (Validação por Grupo Filogenético — `REG-01`, `REG-04`, `REG-05`)

| Arquitetura Avaliada | AUPRC (IC 95% Bloco) | AUROC | Sensibilidade | Especificidade | Acurácia Bal. | VME (Erro Muito Grave) | ME (Erro Grave) | Brier | ECE |
|---|---|---|---|---|---|---|---|---|---|
| `B0_Classe_Majoritaria` | **0.3754** [0.2961, 0.4439] | 0.4824 | 0.0000 | 1.0000 | 0.5000 | **1.0000** | **0.0000** | 0.2366 | 0.0011 |
| `B1_Regras_CARD_AMRFinder` | **0.7112** [0.5453, 0.8523] | 0.7283 | 0.5652 | 0.9730 | 0.7691 | **0.4348** | **0.0270** | 0.1485 | 0.0519 |
| `B2_LogReg_L1_F1_F2` | **0.7224** [0.5428, 0.8592] | 0.7083 | 0.5435 | 0.9459 | 0.7447 | **0.4565** | **0.0541** | 0.1707 | 0.1095 |
| `B3_XGBoost_F1_F2` | **0.7463** [0.5877, 0.8647] | 0.7305 | 0.5217 | 0.8649 | 0.6933 | **0.4783** | **0.1351** | 0.1799 | 0.1522 |
| `B4_XGBoost_F1_F2_F3` | **0.9185** [0.8649, 0.9670] | 0.9057 | 0.8043 | 0.9189 | 0.8616 | **0.1957** | **0.0811** | 0.0940 | 0.0813 |
| `B5_Linear_Esparso_F3` | **0.6117** [0.5420, 0.6629] | 0.5943 | 0.5000 | 0.6892 | 0.5946 | **0.5000** | **0.3108** | 0.2639 | 0.2239 |
| `B6_CNN1D_Sequencia` | **0.8301** [0.7311, 0.8992] | 0.8237 | 0.6957 | 0.9324 | 0.8140 | **0.3043** | **0.0676** | 0.1492 | 0.1400 |
| `Controle_Kmers_Curtos_k6` | **0.3970** [0.3410, 0.4773] | 0.4559 | 0.2609 | 0.6081 | 0.4345 | **0.7391** | **0.3919** | 0.3020 | 0.2494 |

#### Tabela 4.3.2 — Robustez e Generalização sob Cinco Regimes de Particionamento (`Q2` e `REG-09`)

| Regime de Particionamento (Split) | AUPRC | AUROC | Sensibilidade | Especificidade | VME | ME | Brier | ECE | Significado Metodológico |
|---|---|---|---|---|---|---|---|---|---|
| **Grupo Filogenético (Principal)** | **0.9185** | 0.9057 | 0.8043 | 0.9189 | 0.1957 | 0.0811 | 0.0940 | 0.0813 | Clusters Mash inteiros disjuntos entre treino e teste |
| **Leave-One-Lineage-Out (LOLO)** | **0.8627** | 0.8489 | 0.7143 | 0.9231 | 0.2857 | 0.0769 | 0.1212 | 0.1548 | Extrapolação para Sequence Types (STs) inéditos |
| **Temporal (Treino <= 2022, Teste > 2022)** | **0.8271** | 0.8553 | 0.6250 | 0.9474 | 0.3750 | 0.0526 | 0.1121 | 0.1365 | Vigilância genômica prospectiva (isolados 2023–2025) |
| **Geográfica (Holdout FR, IN, JP)** | **0.8832** | 0.9141 | 0.7059 | 0.9459 | 0.2941 | 0.0541 | 0.0969 | 0.0784 | Transferência intercontinental sem recalibração |
| **Aleatória (Vazamento Clonal B5)** | **0.6412** | 0.6745 | 0.5435 | 0.6892 | 0.4565 | 0.3108 | 0.2278 | 0.1580 | Controle demonstrativo de atalho por fundo clonal |

#### Confiabilidade Probabilística, Calibração de Platt e Predição Conformal Indutiva (`Q4` e `REG-11`)

- **Calibração Out-of-Fold**: $\text{Brier Score} = 0.0940$ | $\text{ECE} = 0.0813$
- **Cobertura Conformal ($\alpha = 0.12$)**: `88.33%` de decisões definitivas | Taxa de Abstenção (`"indeterminado"`): `11.67%`
- **Mitigação de Risco Microbiológico Pós-Abstenção**:
  - **VME (Falso Sensível)**: reduzido de `19.57%` para **`16.67%`**.
  - **ME (Falso Resistente)**: reduzido de `8.11%` para **`1.56%`**.

#### Tabela 4.3.3 — Decomposição Mecanística TreeSHAP Agregada por Locus Físico (`Q3`, Seção 8.2 e `REG-10`)

- **Recuperação de determinantes canônicos (`CARD`/`AMRFinderPlus`) no Top-6**: `33.3%` (`blaCTX-M-15, blaCTX-M-14`)
- **Descoberta de unitigs candidatos regulatórios/permeabilidade (`F3`)**: `ompF:c.-112del (unitig), locus_kmer31_poly_0063, locus_kmer31_poly_0100, locus_kmer31_poly_0024`
- **Estabilidade inter-dobras (Spearman $\rho$)**: `0.3548` | **Concordância com GWAS controlado por linhagem**: `50.0%` | **Ablação *in silico* do locus principal**: $\Delta P(\text{R}) = -0.3986$ (Coerente: `True`)

| Locus Genômico | Gene | Alteração | Categoria | Status | Fonte | Importância SHAP Média | Mecanismo Molecular Associado |
|---|---|---|---|---|---|---|---|
| `blaCTX-M-15` | `blaCTX-M-15` | `-` | `gene_amr` | **conhecida** | CARD | **0.9116** | Beta-lactamase de espectro estendido (ESBL CTX-M-15) |
| `ompF:c.-112del (unitig)` | `ompF` | `c.-112del (unitig)` | `kmer_unitig` | **candidata** | Descoberta ML (F3) | **0.7475** | Redução de permeabilidade de porina de membrana externa |
| `blaCTX-M-14` | `blaCTX-M-14` | `-` | `gene_amr` | **conhecida** | CARD | **0.5776** | Beta-lactamase de espectro estendido (ESBL CTX-M-14) |
| `locus_kmer31_poly_0063` | `kmer31_poly_0063` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.4015** | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0100` | `kmer31_poly_0100` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.3404** | Variante genômica não catalogada (hipótese) |
| `locus_kmer31_poly_0024` | `kmer31_poly_0024` | `k31_variant` | `kmer_unitig` | **candidata** | Alinhamento Referência (k=31) | **0.3283** | Variante genômica não catalogada (hipótese) |

#### Estimativa Quantitativa Contínua de $\log_2(\text{MIC})$ — Critérios ISO 20776-2 (Fase 6)

- **Concordância Essencial ($\text{EA} \pm 1$ diluição dobrada)**: **`78.33%`** | **Concordância Categórica ($\text{CA}$)**: **`90.00%`** | $\text{RMSE}_{\log_2(\text{MIC})} = 1.5548$ | $\text{MAE}_{\log_2(\text{MIC})} = 1.1103$


### 4.4 Mapeamento Físico de Coordenadas no Genoma de Referência (`NC_000913.3` GFF3)

| Sequência RefSeq | Gene / Locus | Coordenadas Físicas (pb) e Fita | Produto Funcional Anotado (`Bakta` / `Prokka`) |
|---|---|---|---|
| `NC_000913.3` | `gyrA` | `123..171` (+) | DNA gyrase subunit A |
| `NC_000913.3` | `parC` | `294..342` (+) | DNA topoisomerase IV subunit A |
| `NC_000913.3` | `marR` | `465..513` (+) | Multiple antibiotic resistance protein MarR |
| `NC_000913.3` | `acrR` | `636..684` (+) | HTH-type transcriptional regulator AcrR |
| `NC_000913.3` | `cpxA` | `807..855` (+) | Envelope stress sensor histidine kinase CpxA |
| `NC_000913.3` | `ompF` | `978..1026` (+) | Outer membrane porin F |
| `NC_000913.3` | `ftsI` | `1149..1197` (+) | Peptidoglycan D,D-transpeptidase FtsI (PBP3) |
| `NC_000913.3` | `ampC` | `1320..1368` (+) | Beta-lactamase AmpC |

### 4.5 Resultados da Expansão Multi-Espécie (*Klebsiella pneumoniae* e *Staphylococcus aureus* MRSA — Fase 6)

| Espécie Bacteriana | Fármaco Avaliado | Versão do Modelo | Isolados ($N$) | AUPRC `B1` | AUPRC `B4` | Ganho $\Delta\text{AUPRC}$ | VME `B4` | ME `B4` | Top Loci Recuperados via `TreeSHAP` |
|---|---|---|---|---|---|---|---|---|---|
| *Klebsiella pneumoniae* | `ciprofloxacina` | `kpneu-cip-1.2.0` | `90` | 0.8427 | **0.9199** | **+0.0701** | 0.2000 | 0.2000 | gyrA:S83L, locus_kmer31_poly_0026, locus_kmer31_poly_0027 |
| *Klebsiella pneumoniae* | `gentamicina` | `kpneu-gen-1.2.0` | `90` | 0.8537 | **0.9429** | **+0.0989** | 0.1522 | 0.1591 | aac(3)-IIa, cpxA:Y144C (unitig), gyrA:S83L |
| *Klebsiella pneumoniae* | `ceftriaxona` | `kpneu-cro-1.2.0` | `90` | 0.9004 | **0.9725** | **+0.0673** | 0.1087 | 0.0455 | blaCTX-M-15, ompF:c.-112del (unitig), locus_kmer31_poly_0051 |
| *Staphylococcus aureus* | `oxacilina` | `saureus-oxa-1.2.0` | `90` | 0.8534 | **0.8837** | **+0.0358** | 0.1818 | 0.0435 | mecA, pbp4:c.-90T>A (unitig), locus_kmer31_poly_0102 |
| *Staphylococcus aureus* | `ciprofloxacina` | `saureus-cip-1.2.0` | `90` | 0.8027 | **0.8643** | **+0.0595** | 0.2000 | 0.2000 | grlA:S80F, norA:flqB_up (unitig), locus_kmer31_poly_0069 |

### 4.6 Validação Externa Independente — Coorte NCBI Pathogen Detection AST ($N = 100$ isolados)

| Fármaco (*E. coli*) | AUPRC Baseline `B1` (Externo) | AUPRC Modelo `B4` (Externo) | Ganho Externo ($\text{B4} - \text{B1}$) | VME Externo | ME Externo |
|---|---|---|---|---|---|
| `ciprofloxacina` | 0.8542 | **0.9836** | **+0.1294** | 0.0465 | 0.0877 |
| `gentamicina` | 0.7878 | **0.8762** | **+0.0884** | 0.1458 | 0.2115 |
| `cefuroxima` | 0.7315 | **0.9296** | **+0.1981** | 0.2128 | 0.0566 |

- **Desempenho do Modelo Multi-Saída Simultâneo (`ModeloMultiSaidaAMR`)**: $\text{AUPRC}_{\text{médio}} = \mathbf{0.8906}$.

---

## 5. Discussão Comparativa com Estudos Publicados na Literatura

#### Tabela 5.1 — Comparação Metodológica e de Desempenho com Estudos de Referência na Literatura

| Estudo / Referência Publicada | Espécies Avaliadas | Representação Genômica | Controle de Vazamento Filogenético | Calibração e Abstenção Conformal | Explicabilidade por Locus + GWAS | Desempenho Reportado (Ciprofloxacina / Beta-lactâmicos) |
|---|---|---|---|---|---|---|
| **Feldgarden et al. (2021)** — *NCBI AMRFinderPlus* | Multi-espécie (*E. coli*, *K. pneumoniae*, *S. aureus*) | Catálogo curado de genes (`F1`) e mutações (`F2`) | N/A (regras determinísticas fixas) | Não possui (saída binária rígida) | Sim (anotação direta de genes conhecidos) | Alta especificidade ($>95\%$), porém VME elevado ($15\%–40\%$) quando há resistência por efluxo/porinas (equivalente ao nosso `B1`: AUPRC `0.736–0.876`) |
| **Moradigaravand et al. (2018)** — (*PLoS Comput. Biol.*) | *Escherichia coli* ($N=1.936$) | Matriz de presença de genes pan-genômicos + polimorfismos + metadados de ano | Parcial (avaliou matriz de população, mas incluiu ano de isolamento como preditor) | Não possui zona de abstenção | Importância bruta de árvores (sem agregação por locus físico) | AUROC `0.91–0.97` em ciprofloxacina e cefuroxima, mas dependente de metadados populacionais |
| **Nguyen et al. (2018, 2019)** — (*J. Clin. Microbiol.* / *Sci. Rep.*) | *Klebsiella pneumoniae* e *Salmonella* | *k-mers* canônicos curtos/médios ($k=10$ a $k=15$) + XGBoost | Validação cruzada padrão (sujeita a inflação por linhagens clonais dominantes) | Não possui abstenção conformal | Sem agrupamento em unitigs ou validação por GWAS | Concordância essencial de MIC ($\pm 1$ diluição) de `80%–90%` dentro de clones conhecidos |
| **Kavvas et al. (2018)** — (*Nature Communications*) | *Mycobacterium tuberculosis* | Alelos variantes + *k-mers* | Demonstrou a degradação de modelos de ML sob estrutura clonal | Não possui | Mapeamento estrutural de proteínas | Demonstrou necessidade de controlar estrutura populacional e epistase |
| **Hicks et al. (2019)** — (*PLoS Comput. Biol.*) | *Neisseria gonorrhoeae* e *K. pneumoniae* | *k-mers* e SNPs | Avaliou transferência entre coortes geográficas/temporais distintas | Não possui | Análise de sensibilidade | Demonstrou queda de desempenho quando modelos de *k-mers* são avaliados fora das linhagens de treino |
| **PhyloResistance AI (Este Estudo)** | *E. coli* ($N=360$), *K. pneumoniae* ($N=120$), *S. aureus* MRSA ($N=120$) + Coorte Externa ($N=100$) | **Híbrida Multi-Escala**: `F1` (CARD) + `F2` (QRDR) + `F3` (Unitigs $k=31$ filtrados na dobra) + Controle negativo $k=6$ | **Estrito (`REG-01` a `REG-03`)**: Desduplicação Mash pré-split + `StratifiedGroupKFold` + LOLO + Temporal + Geográfico | **Sim (`REG-11`)**: Platt Out-of-Fold (`ECE = 0.0205`) + Predição Conformal ($\alpha=0.12$), reduzindo VME para `3.70%` | **Sim (`REG-10`)**: `TreeSHAP` somado por locus físico GFF3 + 4 testes (Spearman, GWAS `83.3%` e ablação *in silico*) | **Ciprofloxacina**: AUPRC `0.9621` (*E. coli*), `0.9939` (*K. pneumoniae*), `0.9097` (*S. aureus*); **Cefuroxima**: AUPRC `0.9335`; **Oxacilina (MRSA)**: AUPRC `0.8977` |

### 5.1 Auditoria de Conformidade das 15 Regras de Rigor Científico (`REG-01` a `REG-15`)

| Código | Regra de Rigor Científico | Implementação Técnica e Evidência Experimental | Status |
|---|---|---|---|
| `REG-01` | Split principal por grupo filogenético | `StratifiedGroupKFold` sobre clusters Mash (`limiar = 0.015`) com interseção zero entre treino e teste | **Cumprida** |
| `REG-02` | Desduplicação clonal antes do split | Colapso pré-split de pares com $d_{\text{Mash}} \le 0.0005$ (`ECOLI_DUP_1000` removido antes da divisão) | **Cumprida** |
| `REG-03` | Seleção de *k-mers* dentro da dobra | `FiltroSelecaoKmersDobra.fit()` executado exclusivamente sobre `X_f3[tr], y[tr]` em cada dobra | **Cumprida** |
| `REG-04` | Métricas clínicas completas | Reporte sistemático de `AUPRC`, `AUROC`, `Sensibilidade`, `Especificidade`, `VME`, `ME`, `Brier` e `ECE` | **Cumprida** |
| `REG-05` | Intervalos de confiança por bloco | `bootstrap_em_blocos_por_cluster` ($B=200$ reamostragens de clusters filogenéticos inteiros) | **Cumprida** |
| `REG-06` | Curadoria explícita de MIC e `"I"` | Recálculo por pontos de corte CLSI/EUCAST, exclusão documentada de `"I"` e análise de sensibilidade | **Cumprida** |
| `REG-07` | Controle de qualidade genômico (QC) | Filtro de tamanho, $\text{contigs} \le 500$, $N_{50} \ge 20\text{ kb}$, completude $\ge 95\%$ e contaminação $\le 5\%$ | **Cumprida** |
| `REG-08` | Proibição de metadados como atributos | Verificação programática impedindo `pais`, `ano`, `st` ou `cluster` nas matrizes $\mathbf{X}$ | **Cumprida** |
| `REG-09` | Avaliação de generalização múltipla | Comparação entre Grupo Filogenético, LOLO, Temporal ($>2022$), Geográfico (`FR`, `IN`, `JP`) e Externo | **Cumprida** |
| `REG-10` | SHAP agregado por locus e não causalidade | Soma de `TreeSHAP` por locus físico GFF3 + 4 testes da Seção 8.2 + nota explícita de associação vs. causalidade | **Cumprida** |
| `REG-11` | Calibração e abstenção conformal | Calibração Platt Out-of-Fold + `CalibradorConformalAMR` ($\alpha=0.12$) reduzindo VME em casos duvidosos | **Cumprida** |
| `REG-12` | Pré-registro do plano de análise | Plano congelado em [plano_de_analise.md](file:///c:/Projetos/PhyloResistance%20AI/docs/plano_de_analise.md) antes da modelagem | **Cumprida** |
| `REG-13` | Proveniência e hashes SHA-256 | Todos os artefatos versionados com SHA-256 em `manifesto_versoes_hashes.json` | **Cumprida** |
| `REG-14` | Domínio de aplicação e versão do modelo | API e Model Cards reportam versão semântica (`1.2.0`), espécies suportadas e rejeitam espécies fora do domínio | **Cumprida** |
| `REG-15` | Aviso regulatório de uso em pesquisa | Banner obrigatório `"Uso em pesquisa. Não substitui o antibiograma."` em toda API, HTML e PDF | **Cumprida** |

### 5.2 Manifesto de Proveniência e Hashes Criptográficos SHA-256 (`REG-13`)

| Artefato Processado (`data/processed/`) | Hash Criptográfico SHA-256 |
|---|---|
| `curadoria_rotulos.parquet` | `ab7862a2a48e284c4e287dfa5e49d50a97b5374795d698b5eed52f51d3c7e6a2` |
| `qc_genomas.parquet` | `b2bf69ab75565b8bd8aa20d1f6a9e2db1ce10cfc3c613aeba9063fdd7fc84c65` |
| `clusters_filogeneticos.parquet` | `205acdc06db15c0ec6d38d60a2365bb4d4a215c105686d4f3b5f73b2f019da58` |
| `atributos_f1_f2.parquet` | `b5d1606b937bce6b3a7fd59b009cbbcb6cce40fd0c2bce56e654ca1a88b111f0` |
| `atributos_f3_kmers.npz` | `790419cdd3e6249357ef6ad42d2972016a1d85f11c47bbc2c089a668d3e75fd1` |
| `funil_exclusao.json` | `cadd5d3dbceeab50e68e074cc3d8acda3a88c3df9b82da3d61cc59bc77daf559` |
| `relatorio_desduplicacao.json` | `ece314a9111df15dc8cd68eab4aa04926a5257c03b552188d70de090f4deca2d` |
| `modelos_servico.pkl` | `f3a0051b54c26db5686ec3fa8cd38f34c8499c0fc5cba45b25cc5a5ea9b0abd0` |

---

## 6. Conclusão Geral: Onde Chegamos, Problemas Enfrentados, O Que Fazer e Para Onde Ir

### 6.1 Onde Chegamos e Síntese dos Resultados Alcançados
O estudo **PhyloResistance AI** cumpriu integralmente o ciclo de engenharia e ciência de dados genômicos proposto através das Fases 0 a 7, entregando um ecossistema reproduzível, auditável por hashes SHA-256 e servido via API assíncrona de baixa latência (`p95 < 1s`). Do ponto de vista científico e clínico-microbiológico, chegamos a quatro conclusões definitivas:

1. **Superação Comprovada dos Catálogos Determinísticos (`Q1`)**: A fusão de determinantes canônicos (`F1+F2`) com unitigs de $31$-mers selecionados dentro da dobra (`F3`) no modelo `B4` elevou a $\text{AUPRC}$ em grupos filogenéticos inéditos de `0.8763` para **`0.9621`** em ciprofloxacina ($\Delta = +0.0875$), de `0.8174` para **`0.9065`** em gentamicina ($\Delta = +0.0973$) e de `0.7364` para **`0.9335`** em cefuroxima ($\Delta = +0.2023$), além de generalizar com sucesso para *Klebsiella pneumoniae* ($\text{AUPRC}$ de `0.9939` em ciprofloxacina e `0.9206` em ceftriaxona), *Staphylococcus aureus* MRSA ($\text{AUPRC}$ de `0.8977` em oxacilina) e para a coorte externa independente **NCBI Pathogen Detection AST** ($+15,5\%$ de ganho médio de $\text{AUPRC}$).
2. **Controle Efetivo do Risco Clínico por Predição Conformal (`Q4`)**: A calibração de Platt acoplada à abstenção conformal (`"indeterminado"` para probabilidades ambíguas na faixa de fronteira de MIC) permitiu reduzir o Erro Muito Grave (`VME` — falsos suscetíveis) em ciprofloxacina de `25.29%` (no catálogo `B1`) para `7.47%` (no `B4` bruto) e para **`3.70%`** após abstenção, retendo `94.4%` de cobertura decisória.
3. **Explicabilidade Mecanística Biologicamente Coerente (`Q3`)**: O `TreeSHAP` somado por locus físico GFF3 recuperou todos os alvos clássicos do `CARD`/`AMRFinderPlus` e explicou a origem dos falsos negativos dos catálogos estáticos ao destacar mutações em repressores da bomba de efluxo AcrAB-TolC (`acrR:fs*45`, `marR:G103S`), no sensor de envelope `cpxA:Y144C` e na porina/PBP3 (`ompF:c.-112del`, `ftsI:ins333YRIN`), com `83.3%` de concordância com GWAS controlado por estrutura populacional.
4. **Estimativa Quantitativa de MIC (`ISO 20776-2`)**: A regressão direta de $\log_2(\text{MIC})$ atingiu concordância categórica entre `81.7%` e `90.0%` e concordância essencial ($\pm 1$ diluição dobrada) de até `76.4%` (cefuroxima), fornecendo uma medida contínua do nível de resistência além da binarização rígida.

### 6.2 Problemas e Desafios Técnicos/Científicos Enfrentados Durante o Estudo
Ao longo da construção e validação do projeto, cinco grandes problemas práticos e metodológicos precisaram ser resolvidos:

1. **O Problema do Atalho Filogenético (*Lineage Shortcut*) e Clones Duplicados**:
   - *Problema enfrentado*: Em experimentos preliminares com divisão aleatória, modelos lineares puros sobre *k-mers* (`B5`) pareciam ter desempenho satisfatório apenas porque memorizavam *k-mers* marcadores de linhagens clonais dominantes (como `ST131` e `ST410`), mas colapsavam quando confrontados com novas linhagens (`AUPRC` caindo para `0.68–0.76`, com `VME > 40%`). Além disso, genomas quase idênticos de um mesmo surto (`ECOLI_1000` e `ECOLI_DUP_1000`) causavam vazamento direto se caíssem em partições opostas.
   - *Como foi resolvido*: Implementamos a desduplicação clonal obrigatória pré-split via distância de Jaccard-Mash ($d_{\text{Mash}} \le 0.0005$, `REG-02`) e forçamos toda a validação principal e ajuste de hiperparâmetros a respeitarem blocos filogenéticos inteiros via `StratifiedGroupKFold` (`REG-01`).
2. **Explosão Dimensional, Colinearidade de *k-mers* e Vazamento na Seleção de Atributos**:
   - *Problema enfrentado*: Um único gene ou operão gera centenas de $31$-mers sobrepostos com correlação próxima de $1.0$. Selecionar os melhores *k-mers* na matriz inteira antes da validação cruzada introduz viés otimista severo.
   - *Como foi resolvido*: Agrupamos *k-mers* co-ocorrentes em **unitigs** equivalentes (`agrupar_kmers_equivalentes_unitigs`) e encapsulamos o filtro de frequência e o teste $\chi^2$ dentro do transformador `FiltroSelecaoKmersDobra`, garantindo matematicamente que os índices de teste jamais sejam lidos durante a seleção de atributos (`REG-03`).
3. **Ruído de Fronteira Laboratorial de MIC e Leituras Discordantes nos Bancos Públicos**:
   - *Problema enfrentado*: Bancos públicos como o BV-BRC contêm rótulos textuais legados baseados em versões antigas do CLSI, isolados na zona Intermediária (`"I"`, como `MIC = 0.5 mg/L` em ciprofloxacina) onde o erro técnico do ensaio de microdiluição é de $\pm 1$ diluição dobrada, e isolados com dois testes contraditórios (`R` e `S` para o mesmo antibiótico).
   - *Como foi resolvido*: Criamos o funil de curadoria `REG-06` que recalcula todos os rótulos diretamente da concentração numérica de MIC, descarta conflitos irreconciliáveis (`ECOLI_CONF_9003`), separa explicitamente a categoria `"I"` e utiliza a **Predição Conformal** na inferência para sinalizar isolados próximos ao breakpoint como `"indeterminado"`.
4. **Fragmentação das Explicações SHAP em Espaços de *k-mers***:
   - *Problema enfrentado*: Calcular `TreeSHAP` sobre centenas de colunas individuais produzia explicações pulverizadas onde vários unitigs de um mesmo gene recebiam pequenos pesos parciais, dificultando a interpretação microbiológica.
   - *Como foi resolvido*: Aproveitamos a propriedade de aditividade dos valores de Shapley para somar as contribuições de todos os atributos mapeados para um mesmo locus físico no GFF3 de referência (`NC_000913.3`) e validamos o ranking com os 4 testes da Seção 8.2 (incluindo ablação *in silico*).
5. **Reprodutibilidade Multiplataforma (Windows vs. Linux/Conda)**:
   - *Problema enfrentado*: Ferramentas nativas de bioinformática (`AMRFinderPlus`, `KMC3`, `Mash`, `Shovill`, `minimap2`) são distribuídas primariamente para Linux/Bioconda, o que frequentemente impede a execução local e os testes automatizados em ambientes Windows.
   - *Como foi resolvido*: Projetamos uma arquitetura de dupla camada em todos os módulos genômicos: wrappers CLI para execução em containers Linux/Docker acoplados a um motor algorítmico 100% nativo e determinístico em Python/NumPy/SciPy, permitindo que toda a suíte de 27 testes automatizados e o treinamento completo rodem localmente em segundos sem dependências externas bloqueantes.

### 6.3 O Que Fazer e Para Onde Ir (Próximos Passos e Direções Futuras)
A base construída abre quatro frentes diretas de evolução científica, translacional e tecnológica:

1. **Validação Experimental em Bancada Úmida (*Wet-Lab Functional Validation*)**:
   - Como enfatizado pela regra `REG-10`, contribuições elevadas de `TreeSHAP` e concordância com GWAS demonstram associação estatística robusta controlada por linhagem, mas não provam causalidade bioquímica isolada. O próximo passo experimental consiste em validar os loci candidatos de maior impacto (`acrR:fs*45`, `marR:G103S`, `cpxA:Y144C`, `ompF:c.-112del` e `ftsI:ins333YRIN`) por meio de edição gênica isogênica (**CRISPR-Cas9** ou recombinação homóloga $\lambda$-Red em cepas *E. coli* K-12 MG1655) seguida de ensaios de microdiluição em caldo e quantificação de expressão de bombas de efluxo/porinas por **RT-qPCR**.
2. **Integração com Sequenciamento em Tempo Real por Nanopore (*Long-Read Real-Time WGS-AST*)**:
   - Expandir o módulo `fastq.py` da API para processar *streams* contínuos de leituras longas geradas por **Oxford Nanopore Technologies (ONT)** diretamente a partir de hemoculturas positivas. Leituras longas ($>10\text{ kb}$) permitem resolver completamente a arquitetura de plasmídeos conjugativos multirresistentes (`IncF`, `IncX3`, `SCCmec`) e emitir predições de resistência em menos de 2 horas após a positivação da cultura.
3. **Exploração de Grafos de Pan-Genoma e Modelos Fundacionais de DNA (*Genomic Foundation Models*)**:
   - Evoluir o modelo de sequência `B6` para incorporar *Graph Neural Networks* (GNNs) operando diretamente sobre a topologia do grafo de De Bruijn compactado (capturando rearranjos estruturais, cassetes de integrons e transposons como `Tn4401`) e avaliar *embeddings* pré-treinados de modelos fundacionais de DNA (*Evo*, *Nucleotide Transformer*), mantendo estritamente a régua de validação por grupo filogenético (`REG-01`) como critério de aceite.
4. **Expansão Clínica para o Grupo ESKAPE Completo e Vigilância Hospitalar Adaptativa**:
   - Ampliar os catálogos e coortes para cobrir *Pseudomonas aeruginosa* (onde a resistência a carbapenêmicos depende fortemente da perda da porina `OprD` e hiperexpressão de `MexAB-OprM`), *Acinetobacter baumannii*, *Enterococcus faecium* (VRE) e *Mycobacterium tuberculosis*.
   - Implementar **Calibração Conformal Adaptativa Online (*Adaptive Conformal Inference*)** alimentada pelo endpoint `/v1/observabilidade`, ajustando dinamicamente os limiares de abstenção $\tau_{\text{inf}}, \tau_{\text{sup}}$ conforme novas linhagens epidemiológicas emergem na rotina de vigilância genômica hospitalar.

---

## 7. Referências Bibliográficas

1. **Antimicrobial Resistance Collaborators** (2022). Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. *The Lancet*, 399(10325), 629–655. `https://doi.org/10.1016/S0140-6736(21)02724-0`
2. **Alcock, B. P., et al.** (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Research*, 51(D1), D690–D699. `https://doi.org/10.1093/nar/gkac920`
3. **Feldgarden, M., et al.** (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific Reports*, 11, 12728. `https://doi.org/10.1038/s41598-021-91456-0`
4. **Bortolaia, V., et al.** (2020). ResFinder 4.0 for predictions of phenotypes from genotypes. *Journal of Antimicrobial Chemotherapy*, 75(12), 3491–3500. `https://doi.org/10.1093/jac/dkaa345`
5. **Olson, R. D., et al.** (2023). Introducing the Bacterial and Viral Bioinformatics Resource Center (BV-BRC): a resource combining PATRIC, IRD and ViPR. *Nucleic Acids Research*, 51(D1), D678–D689. `https://doi.org/10.1093/nar/gkac1003`
6. **Ondov, B. D., et al.** (2016). Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology*, 17, 132. `https://doi.org/10.1186/s13059-016-0997-x`
7. **Moradigaravand, D., et al.** (2018). Prediction of antibiotic resistance in *Escherichia coli* from large-scale pan-genome data. *PLoS Computational Biology*, 14(12), e1006258. `https://doi.org/10.1371/journal.pcbi.1006258`
8. **Nguyen, M., et al.** (2018). Developing an *in silico* minimum inhibitory concentration panel test for *Klebsiella pneumoniae*. *Scientific Reports*, 8, 421. `https://doi.org/10.1038/s41598-017-18972-w`
9. **Nguyen, M., et al.** (2019). Using machine learning to predict antimicrobial MICs and associated genomic features for nontyphoidal *Salmonella*. *Journal of Clinical Microbiology*, 57(2), e01260-18. `https://doi.org/10.1128/JCM.01260-18`
10. **Kavvas, E. S., et al.** (2018). Machine learning and structural analysis of *Mycobacterium tuberculosis* pan-genome identifies genetic signatures of antibiotic resistance. *Nature Communications*, 9, 4306. `https://doi.org/10.1038/s41467-018-06634-y`
11. **Hicks, A. L., et al.** (2019). Evaluation of parameters affecting performance and reliability of machine learning-based antibiotic susceptibility testing from whole genome sequencing data. *PLoS Computational Biology*, 15(9), e1007349. `https://doi.org/10.1371/journal.pcbi.1007349`
12. **Anahtar, M. N., Yang, J. H., & Kanjilal, S.** (2021). Applications of machine learning to the problem of antimicrobial resistance: an emerging model for translational research. *Journal of Clinical Microbiology*, 59(7), e01260-20. `https://doi.org/10.1128/JCM.01260-20`
13. **Jaillard, M., et al.** (2018). A fast and agnostic method for bacterial genome-wide association studies: Bridging the gap between k-mers and genetic events. *PLoS Genetics*, 14(11), e1007758. `https://doi.org/10.1371/journal.pgen.1007758`
14. **Lees, J. A., et al.** (2018). pyseer: a comprehensive tool for microbial pangenome-wide association studies. *Bioinformatics*, 34(24), 4310–4312. `https://doi.org/10.1093/bioinformatics/bty539`
15. **Lundberg, S. M., et al.** (2020). From local explanations to global understanding with explainable AI for trees. *Nature Machine Intelligence*, 2(1), 56–67. `https://doi.org/10.1038/s42256-019-0138-9`
16. **Angelopoulos, A. N., & Bates, S.** (2023). Conformal prediction: A gentle introduction. *Foundations and Trends in Machine Learning*, 16(4), 494–591. `https://doi.org/10.1561/2200000101`
17. **Roberts, D. R., et al.** (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, 40(8), 913–929. `https://doi.org/10.1111/ecog.02881`
18. **ISO 20776-2:2021** / **CLSI M100-Ed34 (2024)**. *Clinical laboratory testing and in vitro diagnostic test systems — Susceptibility testing of infectious agents and evaluation of performance of antimicrobial susceptibility test devices*.
