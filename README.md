# PhyloResistance AI

Predição de **Resistência Antimicrobiana (AMR)** e mapeamento de determinantes genômicos (genes de resistência e mutações pontuais) por **Machine Learning Explicável (XAI)** a partir de genomas bacterianos montados (FASTA).

> **Aviso de Uso Responsável (REG-15)**: Projeto de pesquisa científica e demonstração técnica. **Não substitui o antibiograma laboratorial (AST)** nem deve orientar prescrição clínica.

---

## 1. Perguntas Científicas (Pré-Registradas em `docs/plano_de_analise.md`)

- **Q1 (Preditiva)**: Um modelo de ML sobre atributos genômicos (`F1+F2` e `F1+F2+F3`) prevê resistência melhor do que a detecção por regras de genes e mutações já conhecidos (`CARD`/`AMRFinderPlus` — Baseline `B1`) em grupos filogenéticos não vistos?
- **Q2 (Generalização)**: O desempenho se mantém em linhagens (*leave-one-lineage-out*), anos (split temporal) e países (split geográfico) não vistos no treino?
- **Q3 (Interpretabilidade)**: As explicações `TreeSHAP` agregadas por locus recuperam os determinantes biologicamente conhecidos no `CARD` e sugerem candidatos novos com estabilidade entre dobras?
- **Q4 (Confiabilidade)**: A calibração probabilística e a predição conformal permitem abster-se (`indeterminado`) nos casos duvidosos, reduzindo erros graves (*Very Major Error* — `VME` e *Major Error* — `ME`)?

---

## 2. Arquitetura do Pipeline

```text
BV-BRC / NCBI / CARD (download versionado por SHA-256)
   │
   ▼
[1] Ingestão e curadoria padronizada de rótulos (MIC -> S/I/R, CLSI/EUCAST, funil de exclusão)
   │
   ▼
[2] Controle de Qualidade (QC: tamanho, contigs, N50, completude, espécie) + Desduplicação (REG-02)
   │
   ▼
[3] Estrutura Filogenética (distâncias Mash/MinHash -> clusters clonais e ST)
   │
   ▼
[4] Engenharia de Atributos:
    • F1: Genes AMR (AMRFinderPlus / CARD)
    • F2: Mutações pontuais (ex.: gyrA S83L/D87N, parC S80I, ampC)
    • F3: K-mers longos canônicos (k=21..31) em matriz esparsa CSR + unitigs
   │
   ▼
[5] Divisão Anti-Vazamento (Grupo Filogenético StratifiedGroupKFold, LOLO, Temporal, Geográfico)
   │
   ▼
[6] Escada de Modelos (B0 -> B1 Regra CARD -> B2 LogReg L1 -> B3 XGBoost F1+F2 -> B4 XGBoost F1+F2+F3 -> B5 Linear F3 + Controle k=6)
   │
   ▼
[7] Calibração Probabilística + Predição Conformal (Abstenção para incerteza)
   │
   ▼
[8] Explicabilidade (TreeSHAP -> mapeamento k-mer/unitig para locus -> agregação por gene/mutação -> validação CARD)
   │
   ▼
[9] Serviço FastAPI Assíncrono (/v1/predicoes, /v1/modelos, /health, /ready) + Relatórios HTML/JSON
```

---

## 3. Instalação e Execução Rápida (Comando Único)

### Instalar dependências com `uv`
```bash
uv sync --all-extras
```

### Executar o pipeline completo ponta a ponta (Fases 1 a 4 + Relatórios)
```bash
uv run python -m phyloresistance.cli run-all --output-dir data/processed --reports-dir reports
```

### Executar todos os testes automatizados (incluindo Regras de Rigor REG-01 a REG-15)
```bash
uv run pytest -v
```

### Iniciar o serviço de predição FastAPI
```bash
uv run uvicorn phyloresistance.api.app:app --host 0.0.0.0 --port 8000
```

---

## 4. Estrutura do Repositório

```text
phyloresistance/
├── README.md
├── pyproject.toml                  # Dependências Python travadas
├── Makefile                        # Atalhos para reprodução com 1 comando
├── docs/
│   ├── plano_de_analise.md         # Pré-registro de hipóteses, splits e métricas (REG-12)
│   ├── curadoria_de_rotulos.md     # Regras de padronização S/I/R, MIC e funil (REG-06)
│   ├── model_card_ecoli_cip.md     # Model Card e restrições de domínio (REG-14, REG-15)
│   └── decisoes/                   # Architecture Decision Records (ADRs)
├── workflow/                       # Snakemake + ambientes Conda
├── data/
│   ├── raw/                        # Dados brutos versionados por hash SHA-256
│   ├── interim/                    # Artefatos intermediários
│   └── processed/                  # Matrizes esparsas, tabelas curadas e modelos treinados
├── src/phyloresistance/
│   ├── dados/                      # Ingestão BV-BRC, curadoria de rótulos e contratos
│   ├── genomica/                   # QC de montagem, k-mers canônicos, anotação AMR e mutações
│   ├── filogenia/                  # Distâncias Mash/MinHash, clustering e splits anti-vazamento
│   ├── modelos/                    # Escada B0-B5, controle k-curto, calibração e abstenção conformal
│   ├── validacao/                  # Métricas (AUROC, AUPRC, VME, ME, Brier) e bootstrap por bloco
│   ├── explicabilidade/            # TreeSHAP, agregação por locus, validação CARD e relatórios
│   └── api/                        # FastAPI assíncrono, workers e esquemas Pydantic
├── notebooks/                      # Exploração e avaliação de explicabilidade (XAI)
├── tests/                          # Testes unitários, anti-vazamento (REG-01 a REG-15) e API
├── docker/                         # Dockerfiles da API/worker e docker-compose.yml
└── reports/                        # Relatório científico final e métricas exportadas
```
