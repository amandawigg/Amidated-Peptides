# Global profiling of peptide amidation

Code accompanying *Global profiling of peptide amidation to discover bioactive
fragments of the secretome* (Wiggenhorn & Long).

Two steps, one script each:

| | script | input | output |
|---|---|---|---|
| 1. prediction | `predict_amidated_peptides.py` | proteome FASTA | the candidate amidated peptides |
| 2. detection | `amidated_peptide_detection.py` | Supplemental Table 2 | dot products, tiers and detection calls |

Step 2 reads nothing but the published table, so anyone can re-derive every
call in the paper without the raw files.

## Install

```
pip install pandas numpy openpyxl
```

Python 3.8 or newer.

## Detection and scoring

```
python amidated_peptide_detection.py Supplemental_Table_2.xlsx
```

Reads three sheets:

- **Detected Peptides** — the fingerprint ion panel per peptide, from the
  `Ions for Dot Product` column, as published and in the published order
- **Ion Integrations** — each fingerprint ion's m/z and Skyline peak area in
  every replicate, plus the monoisotopic precursor and its mass error
- **Spectral Validation** — the MS1 apex difference from the synthetic
  standard

If `Ions for Dot Product` is missing, the panel falls back to the ions listed
for the 1 uM standard on Ion Integrations; both give the same calls.

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
