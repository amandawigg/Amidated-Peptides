# MS/MS spectral viewer

Builds `viewer_MSMS.html`, the interactive viewer deposited as Supplemental
Data: for every confirmed peptide, the endogenous MS/MS spectrum mirrored
against its 1 µM synthetic standard, with matched fragment ions annotated.

## Build

```bash
make          # extract spectra from mzML, then assemble the viewer
```

## What you need

| Input | Where it comes from | Goes in |
|---|---|---|
| `../data/mzml/*.mzML` | Mendeley Data, DOI 10.17632/9wm5x5ncfb.1 | not in git |
| `../data/peak_boundaries.csv` | one Skyline report, see below | git |
| `../data/Supplemental_Table_2.xlsx` | published with the paper | git |

Fragment ion m/z values are computed from the peptide sequence at build time.
Peptide names, fingerprint ions, dot products, precursor accuracy and detection
calls are read from Supplemental Table 2.

### The Skyline report

Export > Report > Edit List > Add, with these fields:

| Field | Category |
|---|---|
| Protein Name | Proteins |
| Peptide Modified Sequence | Peptides |
| Precursor Mz, Precursor Charge | Precursors |
| Replicate Name | Replicates |
| File Name, Min Start Time, Max End Time | Precursor Results |

`Min Start Time` and `Max End Time` are the integrated peak boundaries; the
extractor averages the MS2 scans inside that window. `File Name` maps each
replicate to its mzML.

## Steps

```
mzML + peak boundaries  --extract_spectra.py-->  spectra.json
spectra.json + Table S2 --build_viewer.py---->  viewer_MSMS.html
```

## Layout

```
viewer/
├── Makefile
├── extract_spectra.py     raw acquisitions -> spectra.json
├── build_viewer.py        spectra.json + Table S2 -> viewer
├── viewer_template.html   interface only, 33 KB
└── data/
    ├── mzml/              from Mendeley, not versioned
    ├── peak_boundaries.csv
    └── Supplemental_Table_2.xlsx
```

The template is versioned; the spectra are regenerated.

## Requirements

Python ≥3.9 with `pyteomics` (mzML reading) and `openpyxl` (Excel).
