# Amidated Peptides

Code accompanying **Wiggenhorn, Lone, Svensson & Long, "Global profiling of peptide amidation to discover bioactive fragments of the secretome."**

Candidate amidated peptides are predicted from primary sequence across the mouse secretome, then tested for endogenous presence against authentic synthetic standards by LC-MS/MS.

---

## Contents

| File | Language | Produces |
|---|---|---|
| `predict_amidated_peptides.py` | Python | predicted amidated peptides → Table S2, "All Predicted Peptides" |
| `tissue_panel_analysis.Rmd` | R | dot products, detection calls and Fig. 1e → Table S2, "spectral valid" |
| `viewer/` | Python | `viewer_MSMS.html`, the interactive spectral viewer (Supplemental Data) |

## Requirements

**Python** ≥3.9 — `openpyxl`, `pyteomics`
**R** ≥4.2 — `dplyr`, `tidyr`, `readxl`, `gplots`

## Running

```bash
# 1. Predict candidate amidated peptides
python predict_amidated_peptides.py

# 2. Score detections and draw Fig. 1e
Rscript -e 'rmarkdown::render("tissue_panel_analysis.Rmd")'

# 3. Build the spectral viewer
cd viewer && make
```

Steps 1 and 2 are independent. Step 3 needs the raw acquisitions from Mendeley
(see `data/README.md`).

## Data

Inputs live in `data/`; see `data/README.md`. Raw LC-MS files are deposited on
Mendeley Data (DOI: 10.17632/9wm5x5ncfb.1) and are not versioned here.

---

## Method summary

### Prediction

Classically secreted proteins are scanned from the signal-peptide cleavage site
to the C-terminus for a site of PAM-mediated amidation. Two cases are searched:

- a glycine-dibasic tripeptide within the sequence — GKR, GKK, GRK or GRR
- a terminal GR or GK, where the protein sequence ends

For each site, an N-terminal cleavage position is taken as the most 3' dibasic
pair (KR, KK, RK, RR) lying downstream of the signal peptide, or the annotated
signal-peptide cleavage site. The predicted peptide is the sequence between
that position and the residue preceding the glycine.

`MAX_LENGTH` sets the upper bound on peptide length. The supplied value of 15
gives the set examined experimentally; raising it gives the full predicted set.

### Detection scoring

Each tissue replicate is compared to an authentic synthetic standard of the
identical sequence run in the same batch. For each peptide, 3–4 fragment ions
are selected from the standard: the most intense qualifying ions, each more
than 0.5 Da — the fragment matching tolerance — from every ion already chosen.
The normalised dot product between sample and standard intensities across those
fixed channels is:

```
dp = Σ(A_sample · A_standard) / ( √Σ(A_sample²) · √Σ(A_standard²) )
```

Ion assignments cover a, b and y ions and their ammonia-loss forms, consistent
with HCD fragmentation. Ion-type labels annotate matched m/z channels and are
not inputs to the score.

Detection uses two tiers. A peptide is confirmed if any replicate reaches dot
product ≥ 0.85 with precursor mass accuracy < 15 ppm and ≥3 fingerprint ions
above background. Once confirmed, further replicates of that peptide qualify at
≥ 0.70 with ≥2 ions. A peptide never confirmed has no detections. Dot products
are reported, and criteria applied, at two decimal places.

For TRH, a tripeptide, the a1, b1 and y1 ions fall below the m/z 120
acquisition floor and the y2 ion gives no signal in the synthetic standard;
b2 and a2 are the observable ions, and detection required dot product ≥ 0.85
across both. This is set in `MIN_IONS_OVERRIDE`.

The analysis gives 62 peptides across 786 tissue replicate detections.

---

## Citation

If you use this code, please cite the paper. Correspondence: jzlong@stanford.edu
