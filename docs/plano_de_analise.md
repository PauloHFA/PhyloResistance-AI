# Plano de Análise Pré-Registrado — PhyloResistance AI

> **Regra de Rigor REG-12**: Este plano de análise define as perguntas científicas, hipóteses, conjuntos de atributos, estratégias de divisão (splits), modelos, métricas e critérios de sucesso **antes** da abertura do conjunto de teste final.

- **Data de Registro**: 2026-10-23
- **Versão do Plano**: `1.0.0`
- **Escopo Piloto**: *Escherichia coli* — Fármacos: `ciprofloxacina`, `gentamicina`, `cefuroxima`
- **Semente Global (Random Seed)**: `42`

---

## 1. Perguntas Científicas e Hipóteses

| ID | Pergunta Científica | Tipo | Hipótese Operacional | Critério Primário de Decisão |
|---|---|---|---|---|
| **Q1** | Um modelo de ML sobre atributos genômicos (F1+F2 e F1+F2+F3) prevê resistência melhor do que a detecção de genes e mutações já conhecidos (CARD/AMRFinderPlus — Baseline B1)? | Preditiva (Principal) | Modelos B3 (XGBoost F1+F2) e/ou B4 (XGBoost F1+F2+F3) apresentam AUPRC superior e menor VME que o Baseline B1 em grupos filogenéticos não vistos. | $\Delta \text{AUPRC}(\text{ML} - \text{B1}) > 0$ com intervalo de confiança (IC 95%) por bootstrap em blocos (clusters). Se não superar, o resultado negativo é reportado como válido (REG-04). |
| **Q2** | O desempenho se mantém em linhagens, anos e países não vistos no treino? | Generalização | A queda de AUPRC e sensibilidade entre validação por grupo filogenético, *leave-one-lineage-out* (LOLO), split temporal e split geográfico permanece controlada e acima do piso B0. | Comparação entre splits: Aleatório (apenas demonstração de vazamento), Grupo Filogenético (principal), LOLO, Temporal (anos antigos $\rightarrow$ recentes) e Geográfico (REG-01, REG-09). |
| **Q3** | As explicações do modelo (TreeSHAP agregado por locus) recuperam determinantes biologicamente conhecidos e sugerem candidatos novos? | Interpretabilidade | Entre os $k$ principais loci por importância SHAP agregada, a maioria corresponde a determinantes catalogados no CARD/AMRFinderPlus para o respectivo fármaco, com estabilidade entre dobras. | Taxa de recuperação no top-$k$ contra o CARD, correlação de Spearman de ranking entre dobras ($\rho \ge 0.6$) e queda coerente de probabilidade no teste de remoção in silico (Seção 8.2, REG-10). |
| **Q4** | Quanta incerteza o modelo carrega, e ele consegue abster-se nos casos duvidosos? | Confiabilidade | A calibração probabilística (isotônica/Platt) reduz o Brier score, e a predição conformal (abstenção em zona de incerteza) reduz erros graves (VME e ME) nos isolados cobertos. | Curva de cobertura vs. erro (VME/ME) com faixa de abstenção configurada por risco clínico ($\alpha = 0.10$, REG-11). |

---

## 2. Escopo Piloto e Domínio de Aplicação

- **Espécie alvo**: *Escherichia coli* (NCBI Taxonomy ID: `562`)
- **Fármacos avaliados (1 modelo por par espécie–fármaco)**:
  1. **Ciprofloxacina** (Fluoroquinolona): mecanismos principais mediados por mutações pontuais na região QRDR (`gyrA` S83L, D87N; `parC` S80I) e genes plasmidiais/efluxo (`qnrS1`, `qnrB`, `aac(6')-Ib-cr`, `oqxAB`, `acrR`).
  2. **Gentamicina** (Aminoglicosídeo): mecanismos principais mediados por enzimas modificadoras de aminoglicosídeos (`aac(3)-IIa`, `aac(3)-IId`, `aac(3)-IV`, `ant(2'')-Ia`, `aph(3')-Ia`) e metilases 16S rRNA (`armA`, `rmtB`).
  3. **Cefuroxima** (Cefalosporina de 2ª geração / $\beta$-lactâmico): mecanismos mediados por $\beta$-lactamases e ESBLs/AmpC (`blaCTX-M-15`, `blaCTX-M-14`, `blaTEM-1B`, `blaSHV-12`, `blaCMY-2`, `blaOXA-1`, mutações promotoras de `ampC`).

---

## 3. Controle de Qualidade (QC) e Desduplicação

Antes de qualquer divisão ou extração de atributos (REG-02, REG-07):
1. **Critérios de montagem (`Escherichia coli`)**:
   - Tamanho total do genoma: `4.4 Mb` a `5.8 Mb` (`4_400_000` a `5_800_000` pb)
   - Número máximo de contigs: $\le 500$
   - $N_{50} \ge 20\,000$ pb
   - Fração de bases ambíguas (`N`): $\le 2\%$
   - Completude estimada $\ge 95\%$ e contaminação $\le 5\%$
   - Espécie declarada idêntica à inferida por distância genômica / assinatura de k-mers.
2. **Desduplicação genômica (REG-02)**:
   - Pares de isolados com similaridade de Jaccard/Mash de k-mers $\ge 0.9995$ (distância $\le 0.0005$) e mesmo fenótipo são colapsados em um único representante por cluster clonal quase-idêntico antes do particionamento.

---

## 4. Conjuntos de Atributos e Prevenção de Vazamento

- **F1 (Genes AMR)**: Matriz binária ($0/1$) de presença/ausência de genes de resistência (identidade $\ge 90\%$, cobertura $\ge 80\%$).
- **F2 (Mutações Pontuais)**: Matriz binária ($0/1$) de mutações conhecidas em genes-alvo (QRDR `gyrA`, `parC`, `parE`, `ampC`, `acrR`) e SNPs contra referência.
- **F3 (K-mers Longos Canônicos)**: Presença/ausência ($0/1$) de k-mers canônicos longos ($k = 21..31$, padrão operacional $k=31$ ou $k=21$), agrupados em *unitigs*/padrões equivalentes e armazenados em formato esparso (`scipy.sparse.csr_matrix`).
- **Controle Negativo**: Frequência relativa de k-mers curtos ($k = 6$) para demonstrar que composição nucleotídica global não localiza determinantes específicos.
- **Exclusão de Metadados (REG-08)**: `pais`, `ano`, `st` (Sequence Type) e `cluster_filogenetico` **jamais** entram na matriz de atributos $X$.
- **Filtragem Estrita na Dobra (REG-03)**:
  - Filtro de frequência mínima ($\ge 2\%$ dos isolados de treino) e máxima ($\le 98\%$ dos isolados de treino) ajustado **exclusivamente** nos índices de treino de cada dobra.
  - Seleção supervisionada ($\chi^2$ / Informação Mútua) e ajuste de hiperparâmetros (Optuna) executados **exclusivamente** nos dados de treino da dobra.

---

## 5. Escada de Modelos

| Nível | Modelo | Atributos | Papel |
|---|---|---|---|
| **B0** | Classificador de Classe Majoritária | Nenhum | Piso trivial |
| **B1** | Regra Determinística (CARD/AMRFinderPlus) | F1 + F2 | **Baseline científico obrigatório (REG-04)**: prediz Resistente (1) se houver pelo menos 1 determinante catalogado para o fármaco |
| **B2** | Regressão Logística Regularizada (L1 / ElasticNet) | F1 + F2 | Modelo linear interpretável |
| **B3** | XGBoost / LightGBM | F1 + F2 | Modelo não-linear sobre determinantes curados |
| **B4** | XGBoost / LightGBM | F1 + F2 + F3 | Modelo completo para captura de variantes não catalogadas e interações epistáticas |
| **B5** | Regressão Logística Regularizada Esparsa | F3 | Modelo baseado apenas em k-mers (comparação com literatura) |
| **Controle** | XGBoost em Frequência de K-mers Curtos ($k=6$) | Composição ($k=6$) | Controle negativo (ablação de resolução genômica) |

---

## 6. Estratégias de Validação e Métricas

### 6.1 Divisões (Splits)
1. **Por Grupo Filogenético (`StratifiedGroupKFold`, $k=5$)**: Validação principal (REG-01). Todos os isolados de um mesmo cluster filogenético/ST permanecem juntos (treino ou teste).
2. **Leave-One-Lineage-Out (LOLO)**: Reserva uma linhagem/ST inteira para teste.
3. **Temporal**: Treino em isolados até o ano de corte ($t \le t_{\text{corte}}$) e teste em isolados posteriores ($t > t_{\text{corte}}$).
4. **Geográfica**: Treino em um conjunto de países/regiões e teste em países retidos.
5. **Aleatória (`StratifiedKFold`)**: Executada **apenas** para quantificar e demonstrar a inflação por vazamento filogenético (reportada explicitamente como limite superior enganoso).

### 6.2 Métricas Obrigatórias (REG-04, REG-05, REG-11)
- **Discriminação**: AUPRC (principal para classes desbalanceadas) e AUROC.
- **Desempenho Clínico-Microbiológico (ISO 20776-2)**:
  - **VME (*Very Major Error*)**: $\frac{\text{FN}}{\text{TP} + \text{FN}}$ (predizer Sensível quando o isolado é verdadeiramente Resistente — erro crítico).
  - **ME (*Major Error*)**: $\frac{\text{FP}}{\text{TN} + \text{FP}}$ (predizer Resistente quando o isolado é verdadeiramente Sensível).
  - Sensibilidade, Especificidade, Acurácia Balanceada e Concordância Categórica.
- **Calibração e Incerteza**:
  - Brier Score e ECE (*Expected Calibration Error*).
  - Taxa de Cobertura e Taxa de Abstenção (predição conformal / limiar de confiança).
- **Inferência Estatística**:
  - Intervalos de confiança de 95% via **bootstrap em blocos por cluster filogenético** ($B = 200$ reamostragens de clusters) e diferença pareada $\Delta(\text{Modelo} - \text{B1})$.
