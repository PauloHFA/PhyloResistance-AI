# Datasheet do Conjunto de Dados — PhyloResistance AI

> **Entregável da Seção 16 (Datasheet for Datasets — padrão Gebru et al.)** e conformidade ética da Seção 15 (`REG-06`, `REG-13`, `REG-15`).

---

## 1. Motivação

- **Com qual finalidade o conjunto de dados foi criado?**
  Permitir o treinamento e a avaliação rigorosa (livre de vazamento clonal/filogenético) de modelos de Machine Learning Explicável (XAI) para predição de Resistência Antimicrobiana (AMR) e regressão de $\log_2(\text{MIC})$ a partir de genomas bacterianos montados.
- **Quem financiou/desenvolveu a curadoria?**
  Equipe de pesquisa do projeto **PhyloResistance AI**, integrando esquemas e ontologias públicas do **BV-BRC (PATRIC)**, **CARD (Comprehensive Antibiotic Resistance Database)**, **NCBI RefSeq / AMRFinderPlus** e **NCBI Pathogen Detection**.

---

## 2. Composição

- **O que representam as instâncias?**
  Cada instância representa um isolado bacteriano sequenciado (genoma montado em contigs FASTA ou leituras FASTQ pós-montagem) pareado com testes laboratoriais de sensibilidade antimicrobiana (AST: categoria $S/I/R$ e valor de MIC em $\text{mg/L}$).
- **Quais espécies e fármacos compõem o dataset?**
  1. **Fase 1 (*Escherichia coli*)**: `ciprofloxacina`, `gentamicina`, `cefuroxima`.
  2. **Fase 2 (*Klebsiella pneumoniae*)**: `ciprofloxacina`, `gentamicina`, `ceftriaxona`.
  3. **Fase 3 (*Staphylococcus aureus*)**: `oxacilina` (MRSA) e `ciprofloxacina`.
- **Quais atributos compõem as matrizes processadas (`data/processed/`)?**
  - `F1`: Presença/ausência binária de genes de resistência adquiridos (`CARD` / `AMRFinderPlus`).
  - `F2`: Presença/ausência binária de mutações pontuais cromossômicas (`gyrA`, `parC`, `grlA`, promotor `ampC`).
  - `F3`: Matriz esparsa (`scipy.sparse.csr_matrix` / `.npz`) de presença/ausência de k-mers canônicos longos ($k=31$) e *unitigs*.
  - **Metadados de estratificação (excluídos de $X$ por `REG-08`)**: `pais`, `ano`, `st` (Sequence Type), `cluster_filogenetico`, `metodo_laboratorio`, `padrao` (`CLSI`/`EUCAST`).

---

## 3. Processo de Coleta e Controle de Qualidade (QC)

- **Critérios de inclusão/exclusão (`REG-06` e `REG-07`)**:
  - Genomas reprovados nos limiares de montagem da espécie (tamanho total fora da faixa esperada, excesso de contigs, $N_{50} < 20\,\text{kb}$, completude $< 95\%$, contaminação $> 5\%$ ou divergência entre espécie declarada e inferida) são removidos antes de qualquer análise.
  - Clones quase idênticos (distância Mash $\le 0.0005$) são colapsados na etapa de desduplicação pré-split (`REG-02`).
  - Rótulos fenotípicos contraditórios ($R$ e $S$ para o mesmo isolado e fármaco) são excluídos e auditados em `data/processed/funil_exclusao.json`.

---

## 4. Privacidade, Ética e Regulação (Seção 15)

- **Privacidade**: As sequências correspondem exclusivamente a genomas bacterianos sem DNA humano. Metadados geográficos restringem-se ao nível de país e ano de coleta. No serviço FastAPI, arquivos enviados para predição são apagados imediatamente após o término do job.
- **Restrição de Uso Clínico (`REG-15`)**: Os dados e modelos destinam-se exclusivamente à pesquisa científica e demonstração metodológica. Não substituem o antibiograma laboratorial nem devem orientar prescrição médica.

---

## 5. Versionamento e Reprodutibilidade (`REG-13`)

- Todos os arquivos gerados em `data/processed/` possuem seus hashes criptográficos **SHA-256** registrados automaticamente em `data/processed/manifesto_versoes_hashes.json`.

