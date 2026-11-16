"""Build small tampered certificate files and run the REAL penta/verify.py on them (--dir)."""
import gzip, json, os, subprocess, sys, time
PY = r"C:\Projects\Claude4\.venv\Scripts\python.exe"
VERIFY = r"C:\Projects\Claude4\penta\verify.py"
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tamper")

def write(case, Z, P, raw=False):
    d = os.path.join(HERE, case); os.makedirs(d, exist_ok=True)
    for name, obj in (("Z", Z), ("P", P)):
        with gzip.open(os.path.join(d, f"certs_{name}.json.gz"), "wt") as fh:
            fh.write(obj if raw else json.dumps(obj))
    return d

R = 0.9375   # dyadic rho0 > 7/8 so that the cell [1/2,1]^3 of a 2x2x2 grid of chart Z is a legal skip cell
emptyP = {"chart": "P", "rho0": R, "grid": [1, 1, 1], "cells": []}
Zskip = {"chart": "Z", "rho0": R, "grid": [2, 2, 2], "cells": [{"index": [1, 1, 1], "node": {"skip": True}}]}
cases = {
    # 1. only a skip cell in Z, nothing in P: must FAIL (grid incomplete)
    "missing_cells": (Zskip, emptyP, []),
    "missing_cells_with_-O": (Zskip, emptyP, ["-O"]),
    # 2. illegal skip leaf
    "illegal_skip": ({"chart": "Z", "rho0": 0.5, "grid": [1, 1, 1], "cells": [{"index": [0, 0, 0], "node": {"skip": True}}]},
                     {"chart": "P", "rho0": 0.5, "grid": [1, 1, 1], "cells": [{"index": [0, 0, 0], "node": {"skip": True}}]}, []),
    # 3. NaN multipliers (json.load accepts the NaN token)
    "nan_multiplier": ('{"chart":"Z","rho0":0.5,"grid":[1,1,1],"cells":[{"index":[0,0,0],"node":{"certs":[[0,[NaN,NaN,NaN,NaN],[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]]]}}]}',
                       json.dumps({"chart": "P", "rho0": 0.5, "grid": [1, 1, 1], "cells": []}), [], True),
    # 4. a single child
    "one_child": ({"chart": "Z", "rho0": 0.5, "grid": [1, 1, 1], "cells": [{"index": [0, 0, 0], "node": {"axis": 0, "children": [{"certs": []}]}}]},
                  emptyP | {"rho0": 0.5}, []),
    # 5. node that is not a dict
    "nondict_node": ({"chart": "Z", "rho0": 0.5, "grid": [1, 1, 1], "cells": [{"index": [0, 0, 0], "node": []}]}, emptyP | {"rho0": 0.5}, []),
    # 6. rho0 mismatch
    "rho0_mismatch": (Zskip, emptyP | {"rho0": 0.5}, []),
}
for name, spec in cases.items():
    Z, P, flags = spec[0], spec[1], spec[2]
    raw = len(spec) > 3
    d = write(name.replace("_with_-O", ""), Z, P, raw)
    t0 = time.time()
    r = subprocess.run([PY] + flags + [VERIFY, "--dir", d, "--workers", "1"], capture_output=True, text=True,
                       env=dict(os.environ, OMP_NUM_THREADS="1"))
    tail = (r.stdout.strip().splitlines() or [""])[-1][:160]
    err = (r.stderr.strip().splitlines() or [""])[-1][:160]
    print(f"{name:24s} exit={r.returncode}  last stdout: {tail!r}  last stderr: {err!r}  [{time.time()-t0:.0f}s]", flush=True)
