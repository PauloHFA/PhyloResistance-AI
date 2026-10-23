# Model Card — PhyloResistance AI (*Staphylococcus aureus* `saureus-amr-1.2.0`)

> **Extensão Fase 6 (Seção 3.2 — Fase 3: Caso de Controle MRSA `mecA`)** | Conformidade `REG-14` e `REG-15`.

---

## 1. Identificação dos Modelos

- **Espécie**: *Staphylococcus aureus* (NCBI Taxonomy ID: `1280`)
- **Modelos por Fármaco**:
  - `saureus-oxa-1.2.0` (Oxacilina — fenótipo MRSA)
  - `saureus-cip-1.2.0` (Ciprofloxacina)
- **Papel Científico no Projeto**: Caso de controle positivo com mecanismo genético canônico altamente definido (`mecA` / `mecC` codificando PBP2a para resistência à oxacilina, mais candidato regulatório `pbp4` para fenótipo BORSA).

---

## 2. Critérios de QC e Pontos de Corte (CLSI M100)

- **Critérios de QC da montagem**: Tamanho entre `2.6 Mb` e `3.2 Mb`, $\le 300$ contigs, $N_{50} \ge 25\,000\text{ pb}$.
- **Pontos de corte de MIC**:
  - *Oxacilina*: Sensível $\le 2.0\text{ mg/L}$, Resistente (MRSA) $\ge 4.0\text{ mg/L}$.
  - *Ciprofloxacina*: Sensível $\le 1.0\text{ mg/L}$, Resistente $\ge 4.0\text{ mg/L}$.

---

## 3. Aviso Regulatório (`REG-15`)

Uso exclusivo em pesquisa científica. Não substitui o antibiograma laboratorial (AST).

