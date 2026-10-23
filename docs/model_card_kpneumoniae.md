# Model Card — PhyloResistance AI (*Klebsiella pneumoniae* `kpneu-amr-1.2.0`)

> **Extensão Fase 6 (Seção 3.2 — Fase 2)** | Conformidade `REG-14` e `REG-15`.

---

## 1. Identificação dos Modelos

- **Espécie**: *Klebsiella pneumoniae* (NCBI Taxonomy ID: `573`)
- **Modelos por Fármaco**:
  - `kpneu-cip-1.2.0` (Ciprofloxacina)
  - `kpneu-gen-1.2.0` (Gentamicina)
  - `kpneu-cro-1.2.0` (Ceftriaxona)
- **Arquitetura**: XGBoost sobre `F1 + F2 + F3` (`PipelineCompletoB4`) com seleção supervisionada de k-mers dentro da dobra (`REG-03`), calibração Platt e abstenção conformal (`CalibradorConformalAMR`), acompanhado de `RegressorMICLog2` e `TreeSHAP` agregado por locus.

---

## 2. Determinantes Biológicos Principais e Critérios de QC

- **Critérios de QC da montagem**: Tamanho entre `5.0 Mb` e `6.2 Mb`, $\le 500$ contigs, $N_{50} \ge 20\,000\text{ pb}$, completude $\ge 95\%$, contaminação $\le 5\%$.
- **Determinantes catalogados (`F1 + F2`)**:
  - *Ciprofloxacina*: `gyrA` S83I, `qnrB1`.
  - *Gentamicina*: `aac(3)-IIa`, `armA`.
  - *Ceftriaxona*: `blaCTX-M-15`, `blaKPC-2`.

---

## 3. Aviso Regulatório (`REG-15`)

Uso exclusivo em pesquisa científica. Não substitui o antibiograma laboratorial (AST).

