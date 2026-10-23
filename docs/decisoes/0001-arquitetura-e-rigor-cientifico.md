# ADR 0001: Arquitetura Híbrida e Garantias de Rigor Científico Anti-Vazamento

- **Status**: Aceito
- **Data**: 2026-10-23

## Contexto

O projeto **PhyloResistance AI** exige:
1. Rigor científico estrito para evitar os problemas clássicos de Machine Learning em genômica bacteriana: vazamento clonal/filogenético (REG-01, REG-02), vazamento de seleção de k-mers fora da dobra (REG-03), ausência de baseline de regras (REG-04) e uso de k-mers curtos que medem apenas composição nucleotídica (Seção 2).
2. Portabilidade e reprodutibilidade: o pipeline precisa orquestrar ferramentas bioinformáticas padrão (Snakemake, AMRFinderPlus, KMC3, Mash, QUAST, CheckM2) em ambientes Linux/containers, mas também executar 100% nativamente em estações Windows para desenvolvimento, testes automatizados de CI e demonstração ponta a ponta sem bloqueios de binários externos.

## Decisão

1. **Camada Genômica Híbrida (`src/phyloresistance/genomica` e `filogenia`)**:
   - Implementamos *Adapters* para as ferramentas CLI externas (`kmc`, `amrfinder`, `mash`, `quast`, `checkm2`) que são acionados automaticamente quando os binários estão presentes no `PATH` ou via container, e um **motor nativo determinístico em Python/NumPy/SciPy** que executa exatamente a mesma lógica matemática (extração de k-mers canônicos $k=21..31$ em matriz esparsa CSR, esboço MinHash/Jaccard para distância filogenética, cálculo de N50/contigs/GC/ambiguidade FASTA, alinhamento/varredura de genes AMR e mutações QRDR/promotoras) quando executado localmente sem os binários Linux.
2. **Encapsulamento Anti-Vazamento (`FoldFeaturePipeline`)**:
   - O filtro de frequência de k-mers e o seletor supervisionado ($\chi^2$ / Informação Mútua) seguem a interface `fit(X_train, y_train)` / `transform(X_test)`, com asserções em tempo de execução e testes automatizados (`tests/test_rigor_cientifico.py`) que garantem que nenhum atributo é filtrado ou selecionado usando o conjunto de teste (REG-03).
3. **Fila Assíncrona Resiliente na API (`src/phyloresistance/api`)**:
   - O serviço FastAPI suporta Redis + RQ quando `REDIS_URL` está configurado e provê um executor assíncrono em *thread pool* com armazenamento persistente em memória/disco quando executado localmente sem Redis, mantendo exatamente o contrato assíncrono `202 Accepted` + `job_id` da Seção 9.

## Consequências

- Qualquer pesquisador ou avaliador consegue reproduzir todo o pipeline (`phyloresistance run-all` ou `pytest`) com um único comando em Windows, macOS, Linux ou Docker.
- As 15 regras de rigor científico (REG-01 a REG-15) são verificadas automaticamente pela suíte de testes.
