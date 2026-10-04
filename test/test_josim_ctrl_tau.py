#!/usr/bin/env python3
"""JoSIM (aether fork) -- the run-level coupling time constant ``.ctrltau`` (aether_sims T35, W0.2 b).

The D8 ``CTRL=`` law makes a junction's Ic follow its control junction's branch current; the card's ``CTRLLAG`` (T3) puts a
first-order low-pass on that input.  ``.ctrltau <seconds>`` is the same filter as a run-level field: every coupled junction whose
card does not carry ``CTRLLAG`` takes it; a card value wins.  These checks run through ``josim-cli`` and compare output files
byte for byte:

  absent_vs_zero   ``.ctrltau 0`` is the memoryless law of record: identical to no ``.ctrltau`` (noise off and noise on)
  field_vs_card    ``.ctrltau X`` on an untagged card is identical to ``CTRLLAG=X`` on the card (noise off and noise on)
  card_wins        a card's ``CTRLLAG=Y`` is kept under ``.ctrltau X``: identical to the card alone
  field_acts       ``.ctrltau X`` (X > 0) changes the run (the field is not a no-op)
  bad_value        ``.ctrltau -1p`` is refused (non-zero exit, the control error)
  run_log          the engine id and the 'CTRL coupling:' line (counts of field / card / memoryless junctions)

The analytic synthetic truth (the step response at tau_c, the rectified mean Ic cut of a THz noise current) needs the per-step Ic
read-back and lives in the superproject (aether_sims/tests/test_t35_ctrl_tau.py, pyjosim).
Usage: test_josim_ctrl_tau.py --josim PATH [--out result.json].  Exit 0 when every check passes, 1 otherwise.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile

ENGINE = "aether-josim-v2-ctrl-tau"


def netlist(card_lag: str | None = None, field: str | None = None, temp: bool = False) -> str:
    """A D8 coupled pair (the 4TJJ M1 card form, 77 K weak-link law) with a driven control arm: the control junction BC carries a
    PWL current step through its own branch; the sample BS is DC-biased near its Ic, so its phase follows Ic(Ictrl)."""
    lag = f", ctrllag={card_lag}" if card_lag is not None else ""
    lines = [
        "* aether T35 ctrl-tau semantics",
        ".model ms jj(rtype=0, rn=0.7, cap=5e-14, icrit=0.006, tc=89, t=77, ictemp=weaklink, n=1.5, "
        f"icctrl=table, ctrltab={{0,1,1,0}}, ctrln=0.0006{lag})",
        ".model mc jj(rtype=0, rn=0.7, cap=5e-14, icrit=0.006, tc=89, t=77, ictemp=weaklink, n=1.5)",
        "BS 1 0 ms ctrl=BC",
        "BC 2 0 mc",
        "IBS 0 1 pwl(0 0 20p 0.0009)",
        "IBC 0 2 pwl(0 0 100p 0 100.1p 0.00025 300p 0.00025 300.1p 0)",
    ]
    if temp:
        lines += [".temp 77", ".neb 1e13"]
    if field is not None:
        lines += [f".ctrltau {field}"]
    lines += [".tran 0.1p 500p 0", ".print phase BS", ".print devv BS", ".print devi BC", ".end", ""]
    return "\n".join(lines)


def run(exe: str, text: str, seed: int = 7) -> tuple[int, bytes, str]:
    with tempfile.TemporaryDirectory() as d:
        cir, csv = os.path.join(d, "t.cir"), os.path.join(d, "t.csv")
        with open(cir, "w") as fh:
            fh.write(text)
        p = subprocess.run([exe, "-o", csv, cir], capture_output=True, text=True,
                           env=dict(os.environ, JOSIM_SEED=str(seed)))
        data = open(csv, "rb").read() if (p.returncode == 0 and os.path.exists(csv)) else b""
        return p.returncode, data, p.stdout + p.stderr


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--josim", required=True)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    res: dict = {"engine": ENGINE, "checks": {}}

    def check(name: str, ok: bool, **info) -> None:
        res["checks"][name] = dict(ok=bool(ok), **info)
        print(f"{'PASS' if ok else 'FAIL'}  {name}  {json.dumps(info)}")

    for temp in (False, True):
        tag = "noise_on" if temp else "noise_off"
        rc0, base, log0 = run(a.josim, netlist(temp=temp))
        rcz, zero, _ = run(a.josim, netlist(field="0", temp=temp))
        check(f"absent_vs_zero[{tag}]", rc0 == 0 and rcz == 0 and base == zero and len(base) > 0, bytes=len(base))
        rcf, fld, logf = run(a.josim, netlist(field="10p", temp=temp))
        rcc, crd, _ = run(a.josim, netlist(card_lag="10p", temp=temp))
        check(f"field_vs_card[{tag}]", rcf == 0 and rcc == 0 and fld == crd and len(fld) > 0, bytes=len(fld))
        rcw, win, logw = run(a.josim, netlist(card_lag="5p", field="10p", temp=temp))
        rc5, c5, _ = run(a.josim, netlist(card_lag="5p", temp=temp))
        check(f"card_wins[{tag}]", rcw == 0 and rc5 == 0 and win == c5 and len(win) > 0, bytes=len(win))
        check(f"field_acts[{tag}]", fld != base, differs=fld != base)
        if not temp:
            mf = re.search(r"CTRL coupling: tau_c = (\S+) s \(\.ctrltau\) on 1 junction\(s\); 0 card CTRLLAG kept; "
                           r"0 memoryless", logf)
            mw = re.search(r"CTRL coupling: tau_c = (\S+) s \(\.ctrltau\) on 0 junction\(s\); 1 card CTRLLAG kept; "
                           r"0 memoryless", logw)
            ok_log = bool(f"Engine: {ENGINE}" in logf and mf and mw and abs(float(mf.group(1)) - 1e-11) < 1e-22
                          and "CTRL coupling" not in log0)
            check("run_log", ok_log, field_line=[ln for ln in logf.splitlines() if "CTRL coupling" in ln],
                  card_line=[ln for ln in logw.splitlines() if "CTRL coupling" in ln])
    rcb, _, logb = run(a.josim, netlist(field="-1p"))
    check("bad_value", rcb != 0 and "CTRLTAU" in logb.upper(), rc=rcb)
    ok = all(c["ok"] for c in res["checks"].values())
    res["ok"] = ok
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(res, fh, indent=1)
    print("ALL PASS" if ok else "FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
