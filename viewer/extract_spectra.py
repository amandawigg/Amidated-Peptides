"""
Extract MS/MS spectra from raw acquisitions for the spectral viewer.

Wiggenhorn et al., "Global profiling of peptide amidation to discover
bioactive fragments of the secretome."

For every peptide and replicate, collects the MS2 scans falling inside the
integrated chromatographic peak, averages them, and writes the result as
spectra.json for build_viewer.py.

Inputs
------
1. mzML files              raw acquisitions (Mendeley DOI: 10.17632/9wm5x5ncfb.1)
2. peak_boundaries.csv     one Skyline report giving, per peptide and replicate:

       Peptide Modified Sequence, Protein Name,
       Precursor Mz, Precursor Charge,
       Replicate Name, File Name,
       Min Start Time, Max End Time

   In Skyline: Export > Report > Edit List > Add, then take Precursor Mz and
   Precursor Charge from Precursors, and Min Start Time / Max End Time /
   File Name from Precursor Results.

Usage
-----
    python extract_spectra.py \\
        --peaks    data/peak_boundaries.csv \\
        --mzml-dir data/mzml \\
        --output   data/spectra.json
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

from pyteomics import mzml

PRECURSOR_TOL = 0.7    # m/z, isolation window half-width for matching MS2 scans
MZ_BIN        = 0.1    # m/z bin width when averaging scans
BASE_PEAK     = 10000  # intensities are normalised to this


def parse_replicate(name: str) -> tuple[str, str]:
    """'BAT_1' -> ('BAT', '1'). Standards return ('standard', '')."""
    if "Std" in name:
        return "standard", ""
    m = re.match(r"^(.*)_(\d+)$", name)
    return (m.group(1), m.group(2)) if m else (name, "1")


def read_peaks(path: Path) -> list[dict]:
    """Skyline peak-boundary report, one row per peptide per replicate."""
    rows = []
    with path.open() as handle:
        for r in csv.DictReader(handle):
            try:
                rows.append({
                    "peptide":   r.get("Protein Name") or r["Peptide Modified Sequence"],
                    "replicate": r["Replicate Name"],
                    "file":      r["File Name"],
                    "mz":        float(r["Precursor Mz"]),
                    "charge":    int(float(r["Precursor Charge"])),
                    "start":     float(r["Min Start Time"]),
                    "end":       float(r["Max End Time"]),
                })
            except (KeyError, ValueError):
                continue          # rows without an integrated peak
    return rows


def average_scans(scans: list[tuple[list, list]]) -> dict:
    """
    Bin and sum scans across the peak, then normalise to the base peak.

Scans are summed across the chromatographic peak rather than taken singly.
    """
    binned: dict[int, float] = defaultdict(float)
    for mzs, ints in scans:
        for mz, it in zip(mzs, ints):
            binned[round(mz / MZ_BIN)] += it
    if not binned:
        return {"mz": [], "intensity": []}

    top = max(binned.values())
    order = sorted(binned)
    return {
        "mz":        [round(k * MZ_BIN, 1) for k in order],
        "intensity": [round(binned[k] / top * BASE_PEAK) for k in order],
    }


def extract(peaks: list[dict], mzml_dir: Path) -> dict:
    by_file: dict[str, list[dict]] = defaultdict(list)
    for row in peaks:
        by_file[row["file"]].append(row)

    collected: dict[tuple[str, str], dict] = {}

    for filename, targets in sorted(by_file.items()):
        path = mzml_dir / filename
        if not path.exists():                       # tolerate .raw -> .mzML naming
            alt = mzml_dir / (Path(filename).stem + ".mzML")
            if not alt.exists():
                print(f"  missing: {filename}")
                continue
            path = alt

        scans: dict[tuple[str, str], list] = defaultdict(list)
        rts:   dict[tuple[str, str], list] = defaultdict(list)

        with mzml.read(str(path)) as reader:
            for spec in reader:
                if spec.get("ms level") != 2:
                    continue
                rt = float(spec["scanList"]["scan"][0]["scan start time"])
                pmz = float(spec["precursorList"]["precursor"][0]
                            ["selectedIonList"]["selectedIon"][0]["selected ion m/z"])
                for t in targets:
                    if (t["start"] <= rt <= t["end"]
                            and abs(pmz - t["mz"]) <= PRECURSOR_TOL):
                        key = (t["peptide"], t["replicate"])
                        scans[key].append((spec["m/z array"], spec["intensity array"]))
                        rts[key].append(rt)

        for key, sc in scans.items():
            t = next(x for x in targets
                     if (x["peptide"], x["replicate"]) == key)
            collected[key] = {
                **average_scans(sc),
                "n_scans":      len(sc),
                "mean_rt":      round(sum(rts[key]) / len(rts[key]), 2),
                "precursor_mz": t["mz"],
                "charge":       t["charge"],
            }
        print(f"  {path.name}: {len(scans)} peptide peaks")

    # Reshape to the viewer's schema
    data: dict[str, dict] = {}
    for (peptide, replicate), spec in collected.items():
        entry = data.setdefault(peptide, {"standard": None, "tissues": {}})
        tissue, rep = parse_replicate(replicate)
        if tissue == "standard":
            entry["standard"] = spec
        else:
            entry["tissues"].setdefault(tissue, {})[rep] = spec

    incomplete = [p for p, v in data.items() if v["standard"] is None]
    if incomplete:
        print(f"  no standard spectrum for: {', '.join(incomplete)}")
    return {p: v for p, v in data.items() if v["standard"] is not None}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--peaks",    type=Path, default=Path("../data/peak_boundaries.csv"))
    p.add_argument("--mzml-dir", type=Path, default=Path("../data/mzml"))
    p.add_argument("--output",   type=Path, default=Path("../data/spectra.json"))
    args = p.parse_args()

    for path in (args.peaks, args.mzml_dir):
        if not path.exists():
            p.error(f"missing input: {path}")

    peaks = read_peaks(args.peaks)
    print(f"peak boundaries: {len(peaks)} peptide-replicate pairs")

    data = extract(peaks, args.mzml_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, separators=(",", ":")))

    n_rep = sum(len(r) for v in data.values() for r in v["tissues"].values())
    print(f"peptides   : {len(data)}")
    print(f"replicates : {n_rep}")
    print(f"output     : {args.output} "
          f"({args.output.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
