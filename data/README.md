# Input data

| File | Status | Source |
|---|---|---|
| `mouse_secretome.tsv` | included | UniProt export, mouse (N = 4975 proteins) |
| `Supplemental_Table_2.xlsx` | included | published with the paper |
| `peak_boundaries.csv` | **you add** | Skyline report, see `viewer/README.md` |
| `mzml/` | not in git | Mendeley Data, DOI 10.17632/9wm5x5ncfb.1 |

`mouse_secretome.tsv` is a UniProt export requiring the keywords "secreted" or
"extracellular" and the absence of "transmembrane", restricted to *Mus
musculus*. Columns used, in order: Entry (0), Gene Names (2), Modified residue
(4), Signal peptide (7), Peptide (8), Sequence (9), Propeptide (10), Chain (11).
