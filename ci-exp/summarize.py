"""Print a one-line fingerprint of a 5MW_Land_BD_Linear result directory.

usage: summarize.py <case_dir> <label>
"""
import hashlib
import os
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "reg_tests", "lib"))
from fast_io import load_output  # noqa: E402
from fast_linearization_file import FASTLinearizationFile  # noqa: E402

case_dir, label = sys.argv[1], sys.argv[2]
case = "5MW_Land_BD_Linear"


def h(arrs):
    m = hashlib.sha256()
    for a in arrs:
        m.update(np.ascontiguousarray(np.asarray(a, dtype=float)).tobytes())
    return m.hexdigest()[:10]


def linhash(suffix):
    f = FASTLinearizationFile(os.path.join(case_dir, f"{case}.{suffix}.lin"))
    keys = [k for k in ["x", "xdot", "u", "y", "A", "B", "C", "D"] if k in f.keys() and f[k] is not None]
    return f, h([f[k] for k in keys]), h([f[k] for k in ["x", "u", "y"] if k in f.keys() and f[k] is not None])


try:
    data, info, _ = load_output(os.path.join(case_dir, f"{case}.outb"))
    ch = info["attribute_names"]
    ci = data[:, ch.index("ConvIter")]
    nj = data[:, ch.index("NumUJac")]
    f1, h1, hop = linhash("1")
    _, hed, _ = linhash("1.ED")
    _, hbd, _ = linhash("1.BD1")
    C = np.asarray(f1["C"], dtype=float)
    print(
        f"FPRINT {label} outb={h([data])} op={hop} lin={h1} ED={hed} BD1={hbd} "
        f"C18_4={C[17, 3]:+.10e} C18_6={C[17, 5]:+.10e} "
        f"ConvIterSum={int(ci.sum())} ConvIterMax={int(ci.max())} NumUJac={int(nj[-1])}"
    )
except Exception as e:  # keep the loop going
    print(f"FPRINT {label} ERROR {type(e).__name__}: {e}")
