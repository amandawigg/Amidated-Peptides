# Global profiling of peptide amidation

Code accompanying *Global profiling of peptide amidation to discover bioactive
fragments of the secretome* (Wiggenhorn & Long).

Two steps, one script each:

| | script | input | output |
|---|---|---|---|
| 1. prediction | `predict_amidated_peptides.py` | proteome FASTA | the candidate amidated peptides |
| 2. detection | `amidated_peptide_detection.py` | Supplemental Table 2 | dot products, tiers and detection calls |

Step 2 reproduces Supplemental Table 2 from the table's own integrations, so
anyone can re-derive every call in the paper without the raw files.

## Install

```
pip install pandas numpy openpyxl
```

Python 3.8 or newer.

## Detection and scoring

```
python amidated_peptide_detection.py Supplemental_Table_2.xlsx
```

Reads two sheets:

- **Ion Integrations** — the five fingerprint ions per peptide and their
  Skyline peak areas in every replicate, plus the monoisotopic precursor and
  its mass error
- **Spectral Validation** — the MS1 apex difference from the synthetic
  standard

and writes three files:

- `fingerprints.csv` — the ion panel used for each peptide
- `detection_by_replicate.csv` — one row per peptide-replicate: dot product,
  ions above background, precursor mass accuracy, apex difference, tier
- `detection_by_peptide.csv` — one row per peptide, with the tissues it was
  detected in

`--out-prefix results/` puts them somewhere other than the working directory.

Expected output on the published table:

```
56 peptides read from the table
56 peptides reported (>=1 tier 1 replicate)
275 peptide-tissue detections
```

## Criteria

| | |
|---|---|
| tier 1 (anchoring) | dot product >= 0.70 across >= 4 of the 5 fingerprint ions |
| tier 2 (distribution) | dot product >= 0.70 across >= 3 fingerprint ions |
| precursor mass accuracy | <= 10 ppm, monoisotopic precursor only |
| MS1 apex | within 0.5 min of the 1 uM synthetic standard |

Both tiers share the dot product threshold; they differ only in how many
fingerprint ions were above background. A peptide is reported only if at
least one replicate reaches tier 1. Tier 2 is then applied to the remaining
replicates of that peptide to map its tissue distribution, and never
establishes a detection on its own.

A window with no integrated monoisotopic precursor has no mass error and no
apex, and cannot meet either tier.

The dot product is the normalized cosine between the fingerprint areas in the
replicate and in the 1 uM standard:

```
dp = sum(A_sample * A_standard) / (||A_sample|| * ||A_standard||)
```

All thresholds are in the CONFIG block at the top of the script.

## Two tolerances, doing different jobs

Easy to conflate, so stated plainly:

- **0.5 Da** is the fragment mass tolerance Skyline used to extract the
  chromatograms. It is spent before this code runs and appears nowhere in it.
- **`MZ_TOL = 0.02` Da** is the tolerance for matching a fingerprint ion to
  its exported transition. Both m/z come from the same file, so this only has
  to absorb rounding. Widening it lets a neighbouring transition stand in for
  a fingerprint ion that is absent.

## Re-running from the raw export

```
python amidated_peptide_detection.py --transitions Transition_Results_final.csv
```

This redoes fingerprint **selection** as well as scoring: the five most
abundant fragment ions in the 1 uM standard, excluding any ion whose blank
area exceeds 10% of its standard area, and skipping any ion within 0.5 Da of
one already chosen, since transitions closer than that are the same
chromatographic peak. Selection used the standard and blank runs only, before
any tissue data were examined.

Only needed to reproduce the selection step. The table path gives identical
calls because the sheet lists the ions that selection chose.

Two cosmetic differences from the published table, neither affecting any
call: the table blanks the dot product where fewer than three fingerprint
ions were above background, and zeroes the areas for windows with no MS2
scan at that precursor.

## Notes

- `FINGERPRINT_OVERRIDE` holds one manual panel, GGFSFRF (PEP-QRFP), whose
  standard is dominated by a co-eluting contaminant; its ions were chosen
  from manually inspected spectra.
- PEP-ADAMTS4 yields only four separable fragment ions and is scored on a
  four-ion fingerprint requiring 4 of 4. Peptides yielding fewer than four
  are not assessable and are reported as undetected.
- Replicate names: Skyline's `Gut_*` and `plasma_*` are the manuscript's
  `Ileum_*` and `Plasma_*`. Both spellings are accepted.
- `blank_01`, `blank_02` and `blank_03` are standard carryover and are
  excluded; the blank filter uses `Blanks` and `Blanks1`.
- z ions are reported as y-NH3. Measured against the matching y ion in the
  1 uM standard, all 743 z/y pairs differ by 17.026549 Da to within 8e-7, so
  these are classical z = y - NH3. The relabelling changes no m/z.

## Citing

Please cite the paper. Supplemental Table 2 holds the data these scripts
read; the raw files are deposited separately.
