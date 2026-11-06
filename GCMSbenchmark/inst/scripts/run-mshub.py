"""Run the standalone MSHub pipeline and export its aligned areas and spectra."""

import argparse
import csv
import os
from pathlib import Path
import runpy
import sys
from time import perf_counter

import h5py
import numpy as np


def step(path, *args):
    sys.argv = [path, *map(str, args)]
    try:
        runpy.run_path(path, run_name="__main__")
    except SystemExit as error:
        if error.code not in (None, 0):
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    options = parser.parse_args()
    source = options.source.resolve()
    inputs = options.input.resolve()
    output = options.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    os.chdir(source)
    sys.path.insert(0, str(source))
    # MSHub predates removal of these builtin aliases in NumPy.
    for name, base in (("float", float), ("int", int), ("bool", bool),
                       ("complex", complex), ("str", str), ("object", object),
                       ("string_", np.bytes_), ("unicode_", np.str_)):
        if name not in np.__dict__:
            setattr(np, name, base)
    db = output / "chromatograms.h5"
    if db.exists():
        raise RuntimeError("Use a fresh output folder for MSHub.")
    started = perf_counter()
    step("proc/io/importmsdata.py", inputs, db, "-f", "netcdf")
    step("proc/preproc/intrapalign.py", db, "--h5writepath", "sp2D")
    step("proc/preproc/noisefilter.py", db, "--h5readpath", "sp2D",
         "--h5writepath", "spproc2D", "--window", 6, "--frame", 50)
    step("proc/preproc/interpalign.py", db, "--h5readpath", "spproc2D",
         "--h5writepath", "spal2D")
    step("proc/preproc/peakdetect.py", db, "--frag_pattern", "deconvolution",
         "--h5readpath", "spal2D", "--individual", "no")
    elapsed = perf_counter() - started
    from proc.io.manageh5db import h5read_strings
    with h5py.File(db) as store:
        result = store["spal2D_peakdetect"]
        area = result["quantity_integrals"][:]
        rt = result["grouped_rts"][:] / 60
        mass = result["grouped_cmz"][:]
        spectra = result["integral_MS_spectra"][:]
        ids = h5read_strings(result["dataset_names"], result["utf_8"])
        if len(ids) != area.shape[0] or len(rt) != area.shape[1]:
            raise ValueError("MSHub export axes do not match the sample and component counts.")
        if not result.attrs.get("deconvolved", False):
            raise ValueError("MSHub did not produce deconvolved spectra.")
    with (output / "peaks.tsv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(("sample", "id", "rt", "area", "spectrum"))
        for j in range(len(rt)):
            spectrum = " ".join(f"{m:g},{v:g}" for m, v in zip(mass, spectra[j]) if v > 0)
            for i, sample in enumerate(ids):
                if area[i, j] > 0:
                    writer.writerow((sample.lstrip("/"), j + 1, rt[j], area[i, j], spectrum))
    (output / "seconds.txt").write_text(f"{elapsed}\n")


if __name__ == "__main__":
    main()
