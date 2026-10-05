"""
amidated_peptide_detection.py
-----------------------------
Reproduce the detection calls in Supplemental Table 2 from the table itself.

  Detected Peptides     the fingerprint ion panel per peptide, from the
                        'Ions for Dot Product' column
  Ion Integrations      each fingerprint ion's m/z and peak area in every
                        replicate, plus the monoisotopic precursor and its
                        mass error
  Spectral Validation   the MS1 apex difference from the 1 uM standard

For each peptide and replicate it computes the normalized dot product
against the 1 uM synthetic standard, counts how many fingerprint ions were
above background, applies the detection criteria, and writes three CSVs.

Run:

    pip install pandas numpy openpyxl
    python amidated_peptide_detection.py Supplemental_Table_2.xlsx
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# ══════════════════════════════════════════════════════════════════════════
# CRITERIA
#
# Both tiers share the dot product threshold; they differ only in how many
# fingerprint ions had to be above background. A peptide is reported only if
# at least one replicate reaches tier 1. Tier 2 is then applied to the
# remaining replicates of that peptide to map its tissue distribution, and
# never establishes a detection on its own.
TIER1 = {'dp': 0.70, 'ions': 4}
TIER2 = {'dp': 0.70, 'ions': 3}
PPM_TOL = 10.0        # precursor mass accuracy, absolute ppm
APEX_TOL = 0.5        # MS1 apex difference from the standard, minutes
MIN_PANEL = 4         # fewer fingerprint ions than this -> not assessable

STANDARD = '1uMStd6'  # the reference replicate every window is scored against

TISSUES = ['BAT', 'Brain', 'Ileum', 'Heart', 'Kidney', 'Liver', 'Lung',
           'Pancreas', 'Quad', 'Spleen', 'iWAT', 'Plasma']

# Skyline named the ileum and plasma runs Gut_* and plasma_*; the manuscript
# and the table use Ileum_* and Plasma_*. Accept either spelling so a sheet
# written under the older convention still joins correctly.
ALIAS = {}
for _i in (1, 2, 3):
    ALIAS[f'Gut_{_i}'] = f'Ileum_{_i}'
    ALIAS[f'plasma_{_i}'] = f'Plasma_{_i}'

SHEETS = {'panel': 'Detected Peptides',
          'areas': 'Ion Integrations',
          'meta': ('Spectral Validation', 'spectral valid')}
# ══════════════════════════════════════════════════════════════════════════


def replicate(name):
    name = str(name).strip()
    return ALIAS.get(name, name)


def tissue_of(rep):
    """Tissue for a replicate, or None for a standard or blank."""
    stem = replicate(rep).rsplit('_', 1)[0]
    return stem if stem in TISSUES else None


def cosine(a, b):
    """Normalized dot product. None when either vector is empty."""
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return None if na == 0 or nb == 0 else float(a @ b / (na * nb))


def ok(x):
    return x is not None and not (isinstance(x, float) and np.isnan(x))


def sheet(book, which):
    want = SHEETS[which]
    want = (want,) if isinstance(want, str) else want
    for name in want:
        if name in book.sheetnames:
            return book[name]
    raise SystemExit(f'the workbook has no {want[0]} sheet; it has '
                     f'{book.sheetnames}')


def header_at(ws, row):
    head = next(ws.iter_rows(min_row=row, max_row=row, values_only=True))
    return {h: i for i, h in enumerate(head) if h}


# ── reading the workbook ──────────────────────────────────────────────────

def read_areas(book):
    """{(peptide, replicate): {ion: area}} and {(peptide, ion): m/z}."""
    ws = sheet(book, 'areas')
    col = header_at(ws, 1)
    areas, mz_of = {}, {}
    for r in ws.iter_rows(min_row=2, values_only=True):
        if r[col['PEP Name']] is None:
            continue
        ion = str(r[col['Fragment Ion']]).strip()
        if ion == 'precursor':
            continue
        key = (r[col['PEP Name']], replicate(r[col['Replicate']]))
        areas.setdefault(key, {})[ion] = float(r[col['Area']] or 0)
        mz_of.setdefault((r[col['PEP Name']], ion),
                         round(float(r[col['Product m/z']]), 4))
    if not areas:
        raise SystemExit('no fragment rows on the Ion Integrations sheet')
    return areas, mz_of


def read_panels(book, mz_of):
    """{peptide: [ion, ...]} from 'Ions for Dot Product', in published order."""
    ws = sheet(book, 'panel')
    col = header_at(ws, 2)
    if 'Ions for Dot Product' not in col:
        raise SystemExit("Detected Peptides has no 'Ions for Dot Product' "
                         "column")
    panels, seq, missing = {}, {}, []
    for r in ws.iter_rows(min_row=3, values_only=True):
        pep = r[col['PEP Name']]
        if pep is None or r[col['Peptide Sequence']] is None:
            continue
        seq[pep] = r[col['Peptide Sequence']]
        ions = [i.strip() for i in
                str(r[col['Ions for Dot Product']]).split(',') if i.strip()]
        for ion in ions:
            if (pep, ion) not in mz_of:
                missing.append(f'{pep} {ion}')
        panels[pep] = ions
    if missing:
        raise SystemExit('listed on Detected Peptides but absent from Ion '
                         'Integrations: ' + ', '.join(missing[:10]))
    return panels, seq


def read_meta(book):
    """{(peptide, replicate): (ppm, apex difference)}."""
    ws = sheet(book, 'meta')
    col = header_at(ws, 2)
    meta = {}
    for r in ws.iter_rows(min_row=3, values_only=True):
        if r[col['PEP Name']] is None:
            continue
        meta[(r[col['PEP Name']], replicate(r[col['Replicate']]))] = (
            r[col['Precursor Accuracy (ppm)']],
            r[col['Apex Difference (min)']])
    return meta


# ── scoring ───────────────────────────────────────────────────────────────

def score(panels, seq, areas, meta):
    """One row per peptide-replicate: dot product, ion count, tier."""
    rows = []
    for peptide, ions in panels.items():
        reference = np.array([areas.get((peptide, STANDARD), {}).get(i, 0.0)
                              for i in ions])
        if np.linalg.norm(reference) == 0:
            continue
        for (pep, rep), found in areas.items():
            if pep != peptide:
                continue
            observed = np.array([found.get(i, 0.0) for i in ions])
            dp = cosine(observed, reference)
            n_above = int((observed > 0).sum())
            error, drift = meta.get((peptide, rep), (None, None))

            # The precursor requirement is carried by the mass accuracy:
            # with no integrated monoisotopic precursor there is no ppm and
            # no apex, and such a window cannot meet either tier.
            passes = (ok(error) and abs(float(error)) <= PPM_TOL
                      and ok(drift) and abs(float(drift)) <= APEX_TOL)

            if len(ions) < MIN_PANEL:
                tier = 'not assessable'
            elif not ok(dp):
                tier = 'none'
            elif passes and dp >= TIER1['dp'] and n_above >= TIER1['ions']:
                tier = 'tier 1'
            elif passes and dp >= TIER2['dp'] and n_above >= TIER2['ions']:
                tier = 'tier 2'
            else:
                tier = 'none'

            rows.append({'peptide': peptide, 'sequence': seq.get(peptide),
                         'replicate': rep, 'tissue': tissue_of(rep),
                         'fingerprint_ions': ', '.join(ions),
                         'panel_size': len(ions),
                         'dot_product': None if not ok(dp) else round(dp, 4),
                         'ions_above_background': n_above,
                         'precursor_ppm': error,
                         'apex_difference_min': drift,
                         'tier': tier})
    return pd.DataFrame(rows)


def summarise(calls):
    """One row per peptide: best tier, replicate counts, per-tissue calls.

    A peptide is reported when at least one replicate reaches tier 1. A
    tissue is Yes when any of its replicates reaches tier 1 or tier 2,
    since tier 2 maps distribution for peptides already anchored elsewhere.
    """
    tissue_calls = calls[calls['tissue'].notna()]
    rows = []
    for peptide, g in tissue_calls.groupby('peptide'):
        t1 = int((g['tier'] == 'tier 1').sum())
        t2 = int((g['tier'] == 'tier 2').sum())
        na = int((g['tier'] == 'not assessable').sum())
        row = {'peptide': peptide,
               'detection_tier': ('tier 1' if t1 else 'tier 2' if t2 else
                                  'not assessable' if na else 'not detected'),
               'tier1_replicates': t1, 'tier2_replicates': t2,
               'reported': 'yes' if t1 else 'no'}
        for tissue in TISSUES:
            sub = g[g['tissue'] == tissue]
            row[tissue] = ('Yes' if sub['tier'].isin(['tier 1', 'tier 2']).any()
                           else 'No')
        rows.append(row)
    return pd.DataFrame(rows).sort_values('peptide')


def panel_table(panels, mz_of, areas):
    return pd.DataFrame(
        [{'peptide': p, 'rank': i, 'ion': ion,
          'product_mz': mz_of[(p, ion)],
          'standard_area': areas.get((p, STANDARD), {}).get(ion, 0.0)}
         for p, ions in panels.items()
         for i, ion in enumerate(ions, start=1)])


# ── entry point ───────────────────────────────────────────────────────────

def main(argv=None):
    global STANDARD
    ap = argparse.ArgumentParser(
        description='Reproduce the detection calls in Supplemental Table 2.')
    ap.add_argument('table', type=Path, help='Supplemental Table 2 (.xlsx)')
    ap.add_argument('--out-prefix', default='',
                    help='prefix for the output files, e.g. results/')
    ap.add_argument('--standard', default=STANDARD,
                    help=f'reference replicate (default {STANDARD})')
    args = ap.parse_args(argv)
    STANDARD = args.standard

    try:
        import openpyxl
    except ImportError:
        raise SystemExit('needs openpyxl: pip install openpyxl')
    book = openpyxl.load_workbook(args.table, read_only=True, data_only=True)

    areas, mz_of = read_areas(book)
    panels, seq = read_panels(book, mz_of)
    meta = read_meta(book)
    print(f'{len(panels)} peptides read from the table')

    calls = score(panels, seq, areas, meta)
    summary = summarise(calls)

    prefix = args.out_prefix
    if prefix and not prefix.endswith(('/', '\\')):
        Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    elif prefix:
        Path(prefix).mkdir(parents=True, exist_ok=True)

    out = [f'{prefix}fingerprints.csv',
           f'{prefix}detection_by_replicate.csv',
           f'{prefix}detection_by_peptide.csv']
    panel_table(panels, mz_of, areas).to_csv(out[0], index=False)
    calls.to_csv(out[1], index=False)
    summary.to_csv(out[2], index=False)

    reported = int((summary['reported'] == 'yes').sum())
    cells = int((summary[TISSUES] == 'Yes').values.sum())
    counts = calls['tier'].value_counts().to_dict()
    print(f'{reported} peptides reported (>=1 tier 1 replicate)')
    print(f'{cells} peptide-tissue detections')
    print('tiers: ' + ', '.join(f'{k} {v}' for k, v in sorted(counts.items())))
    print('wrote ' + ', '.join(out))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
