# Proteinbase Protein Design Competition — Challenge 1 submission

Anti-EGFR binder designs for the three Challenge 1 objectives: binding
affinity against human EGFR, mouse cross-reactivity, and pH-selective binding
(active at pH 6.5, silent at pH 7.4).

## Submission file

**`designs/submission.csv`** — upload this directly. 20 designs, columns
`name,sequence,molecule_class`, all `nanobody`, all 118–122 aa.

| Prefix | Objective | Count |
|---|---|---|
| `EGFR_D3_AFF_*` | 1 — binding affinity, human EGFR domain III | 8 |
| `EGFR_D3_XRC_*` | 2 — human + mouse cross-reactivity | 6 |
| `EGFR_D3_PHS_*` | 3 — pH-selective (binds 6.5, not 7.4) | 6 |

Competition rules checked: 20 proteins (limit 20), every chain 118–122 aa
(limit 10–250), single-chain nanobodies so no `{VH}:{VL}` formatting needed.
The novelty score is checked by the platform after upload.

## Files

| File | What it is |
|---|---|
| `designs/submission.csv` | The deliverable |
| `designs/DESIGN_RATIONALE.md` | Epitope analysis, scaffold choice, per-objective strategy, limitations |
| `designs/build_designs.py` | Regenerates the CSV and re-runs developability QC |
| `designs/qc_report.tsv` | Per-design CDRs, charge at both assay pH values, pI, GRAVY, liability flags |

## Reproduce

```bash
python3 designs/build_designs.py
```

No dependencies beyond the standard library. Prints a summary and rewrites
`submission.csv` and `qc_report.tsv`.

## Honest status

These are knowledge-based designs built on published EGFR structural biology,
not structure-validated ones. The container has no GPU and its network policy
blocks UniProt, NCBI, RCSB and EBI, so no RFdiffusion / ProteinMPNN /
AlphaFold2 filtering was possible. The nanobody scaffold was chosen
specifically because it confines that uncertainty to the three CDR loops while
the framework folds regardless. See the limitations section of
`designs/DESIGN_RATIONALE.md` before relying on these.
