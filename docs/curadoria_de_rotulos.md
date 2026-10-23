# Protocolo de Curadoria Padronizada de Rótulos Fenotípicos

> **Regra de Rigor REG-06**: Todos os rótulos de teste de sensibilidade antimicrobiana (AST) passam por padronização determinística e rastreável, com funil de exclusão documentado.

---

## 1. Fontes e Metadados Exigidos

Cada registro de fenótipo oriundo do **BV-BRC** ou **NCBI Pathogen Detection** deve conter:
- `genome_id`: Identificador único do isolado.
- `especie`: Nome científico padronizado (ex.: *Escherichia coli*).
- `farmaco`: Nome padronizado em português/inglês normalizado (`ciprofloxacina`, `gentamicina`, `cefuroxima`).
- `fenotipo_bruto`: Rótulo original (`Resistant`, `Susceptible`, `Intermediate`, `R`, `S`, `I`, `Susceptible-dose dependent`).
- `mic_mg_l`: Concentração Inibitória Mínima (MIC em $\text{mg/L}$), quando disponível.
- `metodo_laboratorio`: Método de teste (`Broth microdilution`, `Agar dilution`, `Disk diffusion`, `Etest`).
- `padrao`: Comitê de referência (`CLSI` ou `EUCAST`) e versão do documento de ponto de corte.

---

## 2. Pontos de Corte Clínicos Documentados (CLSI M100 / EUCAST) para *Escherichia coli*

Quando o valor numérico de **MIC ($\text{mg/L}$)** estiver disponível, a categoria é **recalculada** a partir de um único ponto de corte padronizado para evitar mistura de interpretações históricas:

| Fármaco | Sensível ($S$, $\text{MIC} \le$) | Intermediário ($I$) | Resistente ($R$, $\text{MIC} \ge$) | Padrão de Referência |
|---|---|---|---|---|
| **Ciprofloxacina** | $\le 0.25\text{ mg/L}$ | $0.5\text{ mg/L}$ | $\ge 1.0\text{ mg/L}$ | CLSI M100 / EUCAST Harmonizado |
| **Gentamicina** | $\le 2.0\text{ mg/L}$ | $4.0\text{ mg/L}$ | $\ge 8.0\text{ mg/L}$ | CLSI M100 / EUCAST |
| **Cefuroxima** | $\le 8.0\text{ mg/L}$ | $16.0\text{ mg/L}$ | $\ge 32.0\text{ mg/L}$ | CLSI M100 (Parenteral) |

---

## 3. Regras de Padronização e Resolução

1. **Preferência ao MIC**: Se `mic_mg_l` for válido ($> 0$), aplica-se a tabela da Seção 2 para obter `S`, `I` ou `R`. Caso contrário, normaliza-se `fenotipo_bruto`.
2. **Tratamento da Categoria Intermediária (`I`)**:
   - **Análise Principal (`modo_intermediario="excluir"`)**: Isolados com fenótipo `I` são excluídos do treino/avaliação binária principal e contabilizados no funil de exclusão.
   - **Análise de Sensibilidade (`modo_intermediario="agrupar_s"` ou `"agrupar_r"`)**: Permite avaliar o impacto de agrupar `I` com `S` (ou `R` para definição não-suscetível).
3. **Resolução de Conflitos**:
   - Se um mesmo `genome_id` possuir múltiplos registros para o mesmo `farmaco`:
     - Se todos concordarem após recálculo por MIC, mantém-se um único registro consolidado (preferindo `Broth microdilution`).
     - Se houver discordância entre `R` e `S` para o mesmo par `(genome_id, farmaco)`, o par é marcado como **conflito irreconciliável** e excluído no funil de exclusão.

---

## 4. Funil de Exclusão Auditável

Toda execução da etapa de curadoria gera um objeto estruturado `FunilExclusao` contendo:
- `total_registros_brutos`: Total inicial ingerido.
- `removidos_especie_fora_escopo`: Registros de outras espécies.
- `removidos_farmaco_fora_escopo`: Registros de outros antibióticos.
- `removidos_fenotipo_invalido`: Sem fenótipo interpretável nem MIC numérico.
- `removidos_intermediario`: Registros classificados como `I` (no modo principal).
- `removidos_conflito_rotulo`: Isolados com leituras contraditórias (`R` vs `S`).
- `removidos_qc_genoma`: Isolados reprovados no controle de qualidade do FASTA.
- `removidos_desduplicacao`: Clones quase idênticos removidos (REG-02).
- `total_aprovados_finais`: Isolados únicos prontos para modelagem.
