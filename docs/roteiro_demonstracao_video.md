# Roteiro de Demonstração em Vídeo Curto (Fase 7 — Comunicação e Portfólio)

> **Entregável da Seção 13 (Fase 7) e Seção 16**: Roteiro enxuto (3 a 4 minutos) para apresentação técnica e científica do **PhyloResistance AI**.

---

## Cena 1 — O Problema Científico e o Risco do Vazamento Filogenético (0:00 – 0:45)

1. Apresentar o objetivo: prever resistência antimicrobiana (AMR) e explicar quais genes/mutações sustentam a decisão a partir do genoma montado (FASTA ou FASTQ).
2. Destacar os dois diferenciais científicos do projeto (Seção 2 e 11):
   - **Prevenção de vazamento clonal (`REG-01`, `REG-02`, `REG-03`)**: bactérias são clonais; dividir treino e teste aleatoriamente infla artificialmente as métricas. O projeto desduplica clones por distância Mash e valida com `StratifiedGroupKFold` por linhagem/ST, *Leave-One-Lineage-Out*, split temporal, split geográfico e coorte externa (*NCBI Pathogen Detection*).
   - **Superação do Baseline por Regras (`B1`, `REG-04`)**: o modelo completo `B4` (`F1 + F2 + F3` com k-mers canônicos $k=31$ e *unitigs*) é comparado contra o catálogo `CARD`/`AMRFinderPlus` com bootstrap em blocos por cluster.

---

## Cena 2 — Reprodução Completa com Um Único Comando (0:45 – 1:45)

No terminal, executar:
```bash
uv run python -m phyloresistance.cli run-all --output-dir data/processed --reports-dir reports
```
Mostrar os artefatos gerados:
- `data/processed/funil_exclusao.json` e `manifesto_versoes_hashes.json` (hashes SHA-256).
- `reports/relatorio_cientifico.md` e as 4 figuras vetoriais em `reports/figuras/`.
- Destacar a queda do **Erro Muito Grave (VME)** graças à calibração e abstenção conformal (`Q4`) e a recuperação dos determinantes do CARD (`gyrA S83L`, `parC S80I`, `blaCTX-M-15`, `mecA`) + candidatos novos via `TreeSHAP` agregado por locus (`Q3`).

---

## Cena 3 — Demonstração ao Vivo da API Assíncrona e Observabilidade (1:45 – 2:45)

1. Subir o serviço FastAPI:
   ```bash
   uv run uvicorn phyloresistance.api.app:app --port 8000
   ```
2. Abrir `http://localhost:8000/docs` (OpenAPI automática):
   - Consultar `GET /v1/modelos?todas_especies=true` (*E. coli*, *K. pneumoniae*, *S. aureus*).
   - Enviar um arquivo FASTA ou FASTQ em `POST /v1/predicoes` $\rightarrow$ receber `202 Accepted` com `job_id`.
   - Consultar `GET /v1/predicoes/{job_id}` (mostrando probabilidade calibrada, estimativa de $\text{MIC mg/L}$, evidências `TreeSHAP` e versões travadas dos bancos) e abrir o relatório visual em `GET /v1/predicoes/{job_id}/relatorio`.
   - Consultar `GET /v1/observabilidade` para exibir latência `p50`/`p95` e monitoramento de *drift*.

---

## Cena 4 — Suíte de Testes de Rigor Científico (2:45 – 3:15)

Executar:
```bash
uv run pytest -v
```
Mostrar todos os testes automatizados verificando `REG-01` a `REG-15`, testes de carga `p95` e extensões da Fase 6 passando em segundos.

