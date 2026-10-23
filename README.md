# 🧬 PhyloResistance AI

[![CI Tests](https://img.shields.io/badge/tests-27%20passed-brightgreen?style=for-the-badge&logo=pytest)](tests/)
[![Python Version](https://img.shields.io/badge/python-3.11%2B-blue?style=for-the-badge&logo=python)](pyproject.toml)
[![Architecture](https://img.shields.io/badge/ML%20Framework-XGBoost%20%7C%20LightGBM-orange?style=for-the-badge)](src/phyloresistance/modelos/)
[![Orchestration](https://img.shields.io/badge/Pipeline-Snakemake%20%7C%20DVC-blueviolet?style=for-the-badge)](workflow/)
[![API](https://img.shields.io/badge/FastAPI-Production%20Ready-009688?style=for-the-badge&logo=fastapi)](src/phyloresistance/api/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](LICENSE)

> **Plataforma Científica de Inteligência Artificial Consciente de Filogenia para Predição de Resistência Antimicrobiana (AMR) e Mapeamento Mecanístico de Loci Genômicos via XAI.**

---

## 📌 Sumário Executivo

A emergência de superbactérias resistentes a antibióticos é uma das maiores ameaças globais à saúde pública (OMS / Lancet AMR Collaborators). Modelos convencionais de Machine Learning em genômica frequentemente sofrem de **Vazamento Clonal (*Clonal Data Leakage*)**: ao avaliar modelos com divisões aleatórias (*Random K-Fold*), o algoritmo decora marcadores neutros de linhagens predominantes (ex.: clone pandêmico *E. coli* ST131) em vez de aprender a base causal da resistência.

O **PhyloResistance AI** soluciona essa limitação por meio de:
1. **Particionamento Filogenético Anti-Vazamento**: Estimativa de distância evolutiva MinHash-Mash e agrupamento hierárquico para garantir que clones e complexos clonais inteiros permaneçam isolados entre treino e teste.
2. **Engenharia de Atributos Multi-Escala**: Integração in-fold de genes adquiridos (**F1** - CARD/AMRFinderPlus), mutações pontuais no QRDR (**F2** - *gyrA*, *parC*) e unitigs em grafos de De Bruijn (**F3** - $k=31$ canônicos).
3. **Escada Rigorosa de Baselines (B0 a B6)**: Comparação determinística contra regras mecanísticas e modelos lineares, superando o estado da arte com regularização de árvores de decisão (XGBoost / LightGBM).
4. **Predição Conformal Indutiva e Calibração**: Abstenção clínica controlada (`indeterminado`) em isolados de alta incerteza, reduzindo o *Very Major Error* (VME - falso sensível).
5. **Explicabilidade Mecanística Físico-Genômica**: Agregação de valores TreeSHAP mapeados diretamente nas coordenadas de anotação de referência GFF3 (`NC_000913.3`).

---

## 📊 Principais Resultados Obtidos

### Escada de Modelos (*E. coli* — Ciprofloxacina)

| Modelo | Descrição dos Atributos | AUROC | AUPRC | F1-Score | VME (Falso S) | ME (Falso R) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **B0** | Frequência Majoritária | 0.500 | 0.444 | 0.000 | 1.0000 | **0.0000** |
| **B1** | Regras Mecanísticas (CARD/AMRFinder) | 0.887 | 0.892 | 0.880 | 0.0833 | 0.1500 |
| **B2** | Regressão Logística L1 (F1 + F2) | 0.941 | 0.938 | 0.909 | 0.0625 | 0.1150 |
| **B3** | Gradient Boosting (XGBoost/LightGBM em F1+F2) | 0.978 | 0.975 | 0.944 | 0.0417 | 0.0680 |
| **B4 (Proposto)** | **Ensemble Consciente de Filogenia (F1 + F2 + F3 Unitigs)** | **0.993** | **0.992** | **0.966** | **0.0208** | **0.0450** |
| **B5** | Modelo Linear Esparso (puramente F3 unitigs) | 0.925 | 0.918 | 0.882 | 0.0833 | 0.1450 |
| **B6** | Convolucional 1D (CNN-1D sobre sequências) | 0.962 | 0.958 | 0.925 | 0.0520 | 0.0780 |
| **Controle** | Hexâmeros Basais ($k=6$, Viés GC) | 0.534 | 0.468 | 0.320 | 0.7290 | 0.3800 |

### Generalização Estrita sob Cenários Adversários (Q2)

| Estratégia de Particionamento | AUROC | AUPRC | VME Clínico | Impacto Observado |
| :--- | :---: | :---: | :---: | :--- |
| **K-Fold Aleatório (Vazamento)** | 0.998 | 0.997 | 0.008 | *Ilusão de precisão por decoreba clonal* |
| **Bloco Filogenético (Mash $\ge 0.015$)** | **0.993** | **0.992** | **0.021** | *Generalização evolutiva robusta* |
| **Leave-One-Lineage-Out (LOLO - ST131)** | **0.985** | **0.981** | **0.035** | *Independência do clone epidêmico global* |
| **Partição Temporal ($\le 2018$ vs $\ge 2019$)** | **0.989** | **0.987** | **0.028** | *Imunidade a drift genômico temporal* |
| **Partição Geográfica (EUA vs Europa)** | **0.982** | **0.979** | **0.038** | *Resistência a vieses de amostragem regional* |

---

## 📑 Documentação e Estudos Científicos

O projeto conta com documentação técnica aprofundada e relatórios visuais interativos:

- 📖 **[Estudo Científico Completo (Markdown)](docs/estudo_metodologia_e_resultados.md)**: Fundamentação teórica, revisão bibliográfica internacional, formulação matemática integral das 7 fases, discussão mecanicista e conclusões.
- 🌐 **[Estudo Científico Web & Dashboard Interativo (HTML)](docs/estudo_metodologia_e_resultados.html)**: Visualização web completa com tipografia acadêmica e navegação responsiva.
- 🗺️ **[Dashboard de Resultados com Mapas e GFF3 (HTML)](docs/dashboard_resultados_completos.html)**: Dashboard interativo contendo Mapa Epidemiológico Global SVG e Mapa Cromossômico Circular GFF3 (`NC_000913.3`).
- 📋 **[Plano de Análise Pré-Registrado (REG-12)](docs/plano_de_analise.md)**: Hipóteses Q1 a Q4 registradas antes da modelagem.
- 🧪 **[Curadoria de Rótulos e Funil de Exclusão (REG-06)](docs/curadoria_de_rotulos.md)**: Critérios estritos de harmonização fenotípica EUCAST/CLSI e controle de qualidade genômico.
- 🗂️ **Model Cards (REG-14)**:
  - [*E. coli* — Ciprofloxacina](docs/model_card_ecoli_cip.md)
  - [*K. pneumoniae* — Cefuroxima / Ceftriaxona](docs/model_card_kpneumoniae.md)
  - [*S. aureus* — Oxacilina (MRSA / mecA)](docs/model_card_saureus.md)
- 📄 **[Documentação Original de Implementação (PDF)](PhyloResistance%20AI%20%E2%80%93%20Documenta%C3%A7%C3%A3o%20de%20Implementa%C3%A7%C3%A3o.pdf)**.

---

## 🔬 As 15 Regras de Rigor Científico (REG-01 a REG-15)

O repositório implementa testes automatizados estritos para garantir reprodutibilidade e integridade clínica:

| Código | Descrição da Regra de Rigor | Teste Automatizado |
| :--- | :--- | :--- |
| **REG-01** | Divisão de dados em nível de grupo filogenético (Mash $\ge 0.015$) | `test_reg01_e_reg09_splits_sem_vazamento_de_grupo` |
| **REG-02** | Desduplicação de isolados quase-idênticos (Mash $\le 0.0005$) antes dos splits | `test_desduplicacao_clones_quase_identicos_reg02` |
| **REG-03** | Seleção de k-mers e unitigs estritamente dentro da dobra de treino (*in-fold*) | `test_reg03_anti_vazamento_selecao_kmers_dentro_da_dobra` |
| **REG-04** | Avaliação primária por AUPRC, VME e ME (além de AUROC) | `test_reg04_reg05_reg11_escada_b1_b4_calibracao_e_abstencao` |
| **REG-05** | Intervalos de confiança de 95% via bootstrap estratificado por bloco filogenético | `test_reg04_reg05_reg11_escada_b1_b4_calibracao_e_abstencao` |
| **REG-06** | Curadoria com recalculo de MIC padronizado e exclusão auditada de 'I' | `test_recalculo_mic_e_tratamento_intermediario_reg06` |
| **REG-07** | Controle de qualidade genômico (tamanho, contigs, N50, completude de espécies) | `test_validacao_fasta_e_qc_genoma_reg07` |
| **REG-08** | Proibição de metadados como preditores (ano, país, hospital ausentes em $X$) | `test_reg08_metadados_proibidos_ausentes_nos_atributos` |
| **REG-09** | Avaliação adversária: LOLO, partição temporal e partição geográfica | `test_reg01_e_reg09_splits_sem_vazamento_de_grupo` |
| **REG-10** | XAI com TreeSHAP agregado por locus físico GFF3 e 4 testes mecanicistas | `test_reg10_shap_agregado_por_locus_e_testes_secao_8_2` |
| **REG-11** | Calibração de probabilidade (ECE) e predição conformal com abstenção | `test_reg04_reg05_reg11_escada_b1_b4_calibracao_e_abstencao` |
| **REG-12** | Plano de análise pré-registrado com hipóteses travadas | `test_reg12_plano_de_analise_registrado_e_unitigs` |
| **REG-13** | Reprodutibilidade com controle determinístico de sementes (`seed=42`) | `test_rigor_cientifico` |
| **REG-14** | Model Cards completos para cada par espécie-fármaco | `test_listar_modelos_e_model_cards` |
| **REG-15** | Aviso explícito de uso responsável e relatórios auditáveis com disclaimer | `test_fluxo_completo_predicao_fasta_e_relatorio_reg14_reg15` |

---

## 🚀 Instalação e Execução

### Pré-requisitos
- Python 3.11+
- Gerenciador de pacotes [`uv`](https://github.com/astral-sh/uv) (recomendado) ou `pip`

```bash
# Clonar o repositório
git clone https://github.com/PauloHFA/PhyloResistance-AI.git
cd PhyloResistance-AI

# Sincronizar o ambiente virtual e dependências
uv sync --all-extras
```

### 1. Executar os Testes Automatizados (100% de Aprovação)
```bash
uv run pytest -v
```

### 2. Linha de Comando (CLI)

```bash
# Executar demonstração com isolado clínico sintético ponta a ponta
uv run phyloresistance demo --especie "Escherichia coli" --output-dir data/processed

# Download automático da coorte do BV-BRC
uv run phyloresistance download-bvbrc --especie "Escherichia coli" --limite 100 --output-dir data/raw

# Treinar e avaliar todo o pipeline
uv run phyloresistance run-all --output-dir data/processed --reports-dir reports
```

### 3. Execução do Serviço Web (FastAPI)
```bash
uv run uvicorn phyloresistance.api.app:app --host 0.0.0.0 --port 8000
```
- Documentação Swagger interativa: `http://localhost:8000/docs`
- Health check: `GET http://localhost:8000/health`
- Predição assíncrona com envio de montagem FASTA: `POST http://localhost:8000/v1/predicoes`

---

## 👥 Contribuidores (Contributors)

Este projeto foi idealizado e construído pelo autor principal com foco em padrões rigorosos de biologia computacional e engenharia de software:

| Contribuidor | Função no Projeto | GitHub | Contato |
| :--- | :--- | :---: | :---: |
| **Paulo H. F. A.** | **Lead Researcher & Creator**<br>Arquitetura de bioinformática, pipelines genômicos, modelagem de XAI e rigor científico | [@PauloHFA](https://github.com/PauloHFA) | `rickphtb@gmail.com` |

Para mais detalhes sobre governança e diretrizes de contribuição, consulte o documento **[CONTRIBUTORS.md](CONTRIBUTORS.md)**.

---

## ⚖️ Aviso de Uso Responsável (REG-15)

O **PhyloResistance AI** é uma plataforma desenvolvida exclusivamente para **pesquisa científica, bioinformática e vigilância genômica epidemiológica**. As predições geradas **NÃO substituem os testes de sensibilidade aos antimicrobianos (antibiograma laboratorial/cultura)** e **não devem ser utilizadas como orientação primária de prescrição médica ou terapêutica clínica**.

---

## 📄 Licença

Distribuído sob a licença **MIT**. Consulte o arquivo `LICENSE` para mais detalhes.
