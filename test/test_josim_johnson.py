#!/usr/bin/env python3
"""JoSIM (aether fork) -- the Johnson source against synthetic truth (aether_sims T33-R1, R0.3).

The engine's thermal noise is a current source across every resistor and junction shunt at temperature T (``.temp`` or the instance
``TEMP=``): a HELD process on the absolute time grid k / neb, per-sample variance (4kT / R)(neb / 2), whose exact charge the engine
receives in every step. These checks run circuits whose answer is known in closed form:

  equipartition  R || L carries <i_L^2> = kT / L in the inductor, whatever R, grounded or two-node
  psd            R shorted by 1 fH: the short carries the source's whole current, one-sided PSD 4kT / R
  step           the equipartition reading does not depend on the step the engine takes (0.1 / 0.05 / 0.025 / 0.0125 ps requested)
  ah             an overdamped junction's time-averaged voltage below Ic: the Ambegaokar-Halperin / Stratonovich closed form
  escape         the record line junction (678.8 uA at 65 K) at 0.95 Ic: its slip rate against thermal activation at 65 K

Two backends: ``--josim <josim-cli>`` (subprocess; JOSIM_SEED seeds the run) or ``--backend pyjosim`` (in process).
Usage: test_josim_johnson.py [--josim PATH | --backend pyjosim] [--neb 1e13] [--plateau] [--quick] [--out result.json]
Exit 0 when every bar at the tested neb passes, 1 otherwise, 77 when numpy is not available (ctest SKIP_RETURN_CODE).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
import time

try:
    import numpy as np
except ImportError:  # pragma: no cover
    print("SKIP: numpy is not available")
    sys.exit(77)

KB = 1.380649e-23
PHI0 = 2.067833848e-15
NEBS = (2e12, 5e12, 1e13, 4e13)


# ---------------------------------------------------------------------------------------------------------------- backends
class CliBackend:
    def __init__(self, exe: str):
        self.exe = exe

    def run(self, text: str, seed: int) -> tuple[dict, dict]:
        with tempfile.TemporaryDirectory() as d:
            cir, csv = os.path.join(d, "t.cir"), os.path.join(d, "t.csv")
            with open(cir, "w") as fh:
                fh.write(text)
            env = dict(os.environ, JOSIM_SEED=str(int(seed)))
            p = subprocess.run([self.exe, "-o", csv, cir], capture_output=True, text=True, env=env)
            if p.returncode != 0:
                raise RuntimeError(f"josim-cli failed ({p.returncode}): {p.stderr[-2000:]}")
            info = {}
            for ln in p.stdout.splitlines():
                i = ln.find("Transient step:")
                if i >= 0:
                    toks = ln[i:].replace("(", " ").replace(",", " ").split()
                    info = {"engine_step": float(toks[2]), "requested_step": float(toks[5]), "halvings": int(toks[7])}
            with open(csv) as fh:
                head = fh.readline().strip().split(",")
            data = np.loadtxt(csv, delimiter=",", skiprows=1, ndmin=2)
            tr = {h.strip('"'): data[:, j] for j, h in enumerate(head)}
            return tr, info


class PyjosimBackend:
    def __init__(self):
        import pyjosim  # noqa: F401

        self.pj = pyjosim

    def run(self, text: str, seed: int) -> tuple[dict, dict]:
        import contextlib
        import io

        pj = self.pj
        fd, path = tempfile.mkstemp(suffix=".cir")
        try:
            with os.fdopen(fd, "w") as fh:
                fh.write(text)
            with contextlib.redirect_stdout(io.StringIO()):
                cli = pj.CliOptions([path])
                inp = pj.Input(cli.analysis_type, cli.verbose, cli.minimal)
                inp.parse_input(cli.cir_file_name)
                if len(inp.parameters) > 0:
                    inp.parameters.parse()
                inp.parse_models()
                inp.netlist.expand_subcircuits()
                inp.netlist.expand_maindesign()
                inp.identify_simulation()
                pj.seed_noise(int(seed))
                m = pj.Matrix(inp)
                sim = pj.Simulation(inp, m)
                tr = {t.name: np.asarray(t.data) for t in pj.Output(inp, m, sim).traces}
                info = {"engine_step": sim.step_size(), "requested_step": sim.requested_step(), "halvings": sim.halvings()}
        finally:
            os.unlink(path)
        return tr, info


# ---------------------------------------------------------------------------------------------------------------- closed forms
def ah_velocity(i: float, D: float, n: int = 2048) -> float:
    """Overdamped RSJ with white thermal noise: <V> / (Ic R) = 2 pi D (1 - exp(-2 pi i / D)) / I,
    I = int_0^2pi dx int_0^2pi dy exp{[-i y - cos x + cos(x - y)] / D}, D = kT / E_J = 2 pi k T / (Phi0 Ic)
    (Stratonovich 1967; Ambegaokar & Halperin, PRL 22, 1364 (1969))."""
    x = (np.arange(n) + 0.5) * (2 * np.pi / n)
    X, Y = np.meshgrid(x, x, indexing="ij")
    E = (-i * Y - np.cos(X) + np.cos(X - Y)) / D
    m = E.max()
    s = np.exp(E - m).sum() * (2 * np.pi / n) ** 2
    lognum = math.log(2 * math.pi * D) + math.log1p(-math.exp(-2 * math.pi * i / D))
    return math.exp(lognum - (m + math.log(s)))


def kramers_overdamped(ic: float, rn: float, i: float, temp: float) -> float:
    """Kramers-Smoluchowski escape rate of an overdamped junction: (omega_c / 2 pi) sqrt(1 - i^2) exp(-dU / kT),
    omega_c = 2 pi Ic Rn / Phi0, dU = E_J [2 sqrt(1 - i^2) - 2 i arccos i]."""
    ej = PHI0 * ic / (2 * math.pi)
    du = ej * (2 * math.sqrt(1 - i * i) - 2 * i * math.acos(i))
    wc = 2 * math.pi * ic * rn / PHI0
    return wc / (2 * math.pi) * math.sqrt(1 - i * i) * math.exp(-du / (KB * temp))


def equipartition_theory(r: float, l: float, h: float, tau: float, nf: int = 200001) -> float:
    """<i_L^2> L / kT the discrete system delivers: the source's PSD (held, box-averaged over the step, folded to the step grid)
    through the BDF2 transfer function of the inductor current. The deficit from 1 is the finite noise bandwidth and the step."""
    alpha = 2 * h * r / l
    f = np.linspace(0, 1 / (2 * h), nf)
    z = np.exp(1j * 2 * np.pi * f * h)
    H = alpha / ((3 + alpha) - 4 / z + 1 / z ** 2)
    S = np.zeros_like(f)
    for mm in range(-100, 101):
        fm = f + mm / h
        S += np.sinc(fm * tau) ** 2 * np.sinc(fm * h) ** 2
    return float(np.trapezoid(S * np.abs(H) ** 2, f) * 4 / r * l)


# ---------------------------------------------------------------------------------------------------------------- circuits
def _controls(temp, neb, tstep, tstop, prstep):
    return [f".temp {temp:g}", f".neb {neb:g}", f".tran {tstep:g} {tstop:g} 0 {prstep:g}"]


def rl_netlist(pairs, *, temp, neb, tstep, tstop, prstep):
    ln = ["* R || L equipartition"]
    for nm, r, l, mode in pairs:
        a, b = f"a{nm}", ("0" if mode == "g" else f"b{nm}")
        ln += [f"R{nm} {a} {b} {r:g}", f"L{nm} {a} {b} {l:g}"]
        if mode != "g":
            ln.append(f"LG{nm} {b} 0 1e-15")
    ln += _controls(temp, neb, tstep, tstop, prstep)
    ln += [f".print DEVI L{nm}" for nm, *_ in pairs] + [".end"]
    return "\n".join(ln) + "\n"


def trace(tr, key):
    for k in (key, key.upper()):
        if k in tr:
            return np.asarray(tr[k], float)
    raise KeyError(f"{key} not in {list(tr)[:6]}")


PAIRS = [(1.0, 10e-12), (0.7, 13e-12), (2.0, 20e-12), (10.0, 20e-12)]


def check_equipartition(be, neb, *, temp=65.0, tstep=1e-13, tstop=400e-9, seeds=(1, 2)):
    pairs = [(f"{j}{m}", r, l, m) for j, (r, l) in enumerate(PAIRS) for m in ("g", "n")]
    acc = {p[0]: [] for p in pairs}
    info = {}
    for sd in seeds:
        tr, info = be.run(rl_netlist(pairs, temp=temp, neb=neb, tstep=tstep, tstop=tstop, prstep=0.5e-12), sd)
        for nm, r, l, m in pairs:
            i = trace(tr, f"I(L{nm})")[200:]
            acc[nm].append(float(np.var(i) * l / (KB * temp)))
    rows = []
    for nm, r, l, m in pairs:
        v = float(np.mean(acc[nm]))
        rows.append({"R_ohm": r, "L_pH": l * 1e12, "reference": "grounded" if m == "g" else "two-node", "ratio": v,
                     "theory": equipartition_theory(r, l, info.get("engine_step", tstep), 1 / neb), "pass": bool(abs(v - 1) <= 0.05)})
    return {"neb": neb, "engine": info, "rows": rows, "pass": all(r_["pass"] for r_ in rows)}


def check_psd(be, neb, *, temp=65.0, tstep=1e-13, tstop=200e-9, seed=3):
    bands = ((1e9, 10e9), (10e9, 50e9), (50e9, 200e9), (200e9, 500e9))
    out = {"neb": neb, "rows": []}
    for r in (1.0, 10.0):
        spec = {}
        for mode in ("g", "n"):
            # the same seed and the same single noise element: the two references carry the SAME sample stream
            tr, info = be.run(rl_netlist([("0", r, 1e-15, mode)], temp=temp, neb=neb, tstep=tstep, tstop=tstop, prstep=tstep), seed)
            i = trace(tr, "I(L0)")[100:]
            i = i - i.mean()
            n = len(i)
            f = np.fft.rfftfreq(n, tstep)
            p = (np.abs(np.fft.rfft(i)) ** 2) * 2 * tstep / n
            s0 = 4 * KB * temp / r
            spec[mode] = {f"{lo / 1e9:g}-{hi / 1e9:g} GHz": float(p[(f >= lo) & (f < hi)].mean() / s0) for lo, hi in bands}
            spec[mode]["1-500 GHz"] = float(p[(f >= 1e9) & (f < 500e9)].mean() / s0)
        ok_b = all(abs(spec[m_][k] - 1) <= 0.10 for m_ in spec for k in spec[m_] if k != "1-500 GHz")
        ratio = spec["n"]["1-500 GHz"] / spec["g"]["1-500 GHz"]
        out["rows"].append({"R_ohm": r, "grounded": spec["g"], "two_node": spec["n"], "two_node_over_grounded": ratio,
                            "pass": bool(ok_b and abs(ratio - 1) <= 0.03)})
    out["pass"] = all(r_["pass"] for r_ in out["rows"])
    return out


def check_step(be, neb, *, temp=65.0, tstop=100e-9, seed=4):
    pairs = [("0g", 1.0, 10e-12, "g"), ("0n", 1.0, 10e-12, "n"), ("3g", 10.0, 20e-12, "g")]
    rows = []
    for h in (1e-13, 5e-14, 2.5e-14, 1.25e-14):
        tr, info = be.run(rl_netlist(pairs, temp=temp, neb=neb, tstep=h, tstop=tstop, prstep=0.5e-12), seed)
        rows.append({"requested_ps": h * 1e12, "engine": info,
                     **{nm: float(np.var(trace(tr, f"I(L{nm})")[200:]) * l / (KB * temp)) for nm, r, l, m in pairs}})
    ref = rows[0]
    for r_ in rows:
        r_["pass"] = bool(all(abs(r_[nm] / ref[nm] - 1) <= 0.05 for nm, *_ in pairs))
    return {"neb": neb, "rows": rows, "pass": all(r_["pass"] for r_ in rows)}


def ah_netlist(*, temp, neb, ic, rn, cap, fracs, tstep, tstop, prstep):
    ln = ["* overdamped junctions at a DC bias", f".model jah jj(rtype=0, rn={rn:g}, cap={cap:g}, icrit={ic:g})"]
    for j, fr in enumerate(fracs):
        ln += [f"B{j} n{j} 0 jah", f"I{j} 0 n{j} pwl(0 0 100p {fr * ic:.12g})"]
    ln += _controls(temp, neb, tstep, tstop, prstep)
    ln += [f".print PHASE B{j}" for j in range(len(fracs))] + [".end"]
    return "\n".join(ln) + "\n"


def check_ah(be, neb, *, ic=2e-6, rn=20.0, cap=1e-18, tstep=1e-13, tstop=2e-6, seed=5, temps=(4.0, 40.0, 65.0, 77.0), fracs=(0.8, 0.9, 0.95)):
    rows = []
    for T in temps:
        tr, info = be.run(ah_netlist(temp=T, neb=neb, ic=ic, rn=rn, cap=cap, fracs=fracs, tstep=tstep, tstop=tstop, prstep=1e-12), seed)
        t = trace(tr, "time")
        m = t >= 0.2e-9
        D = 2 * math.pi * KB * T / (PHI0 * ic)
        for j, fr in enumerate(fracs):
            ph = trace(tr, f"P(B{j})")
            v = (ph[m][-1] - ph[m][0]) / (2 * math.pi) * PHI0 / (t[m][-1] - t[m][0])
            v_ah = ah_velocity(fr, D) * ic * rn
            slips = (ph[m][-1] - ph[m][0]) / (2 * math.pi)
            rows.append({"T_K": T, "i": fr, "D": D, "v_sim_over_IcR": v / (ic * rn), "v_ah_over_IcR": v_ah / (ic * rn), "ratio": v / v_ah, "slips": slips,
                         "pass": bool(abs(v / v_ah - 1) <= 0.10)})
    return {"neb": neb, "junction": {"Ic_uA": ic * 1e6, "Rn_ohm": rn, "C_F": cap, "beta_c": 2 * math.pi * ic * rn * rn * cap / PHI0}, "rows": rows,
            "pass": all(r_["pass"] for r_ in rows)}


# the record line-class junction: icrit is the T = 0 base of the weak-link law pinned to 678.8 uA at 65 K (tc 89 K, n 1.5)
REC_ICRIT, REC_TC, REC_N, REC_RN, REC_CAP, REC_T = 4.847573276049231e-3, 89.0, 1.5, 0.7, 5e-14, 65.0


def rec_ic(T=REC_T):
    return REC_ICRIT * (1 - T / REC_TC) ** REC_N


def check_escape(be, neb, *, frac=0.95, n=8, tstep=1e-13, tstop=60e-9, seeds=(11, 12)):
    ic = rec_ic()
    card = f".model jrec jj(rtype=0, rn={REC_RN:g}, cap={REC_CAP:g}, icrit={REC_ICRIT:.15g}, tc={REC_TC:g}, t={REC_T:g}, ictemp=weaklink, n={REC_N:g})"
    tot, t_tot, info = 0.0, 0.0, {}
    for sd in seeds:
        ln = ["* escape", card]
        for j in range(n):
            ln += [f"B{j} n{j} 0 jrec", f"I{j} 0 n{j} pwl(0 0 0.2n 0 1n {frac * ic:.12g})"]
        ln += _controls(REC_T, neb, tstep, tstop, 1e-12) + [f".print PHASE B{j}" for j in range(n)] + [".end"]
        tr, info = be.run("\n".join(ln) + "\n", sd)
        t = trace(tr, "time")
        m = t > 2e-9
        for j in range(n):
            ph = trace(tr, f"P(B{j})")
            tot += round((ph[m][-1] - ph[m][0]) / (2 * math.pi))
        t_tot += n * (t[m][-1] - t[m][0])
    rate = tot / t_tot
    D = 2 * math.pi * KB * REC_T / (PHI0 * ic)
    r_ah = ah_velocity(frac, D) * ic * REC_RN / PHI0
    r_kr = kramers_overdamped(ic, REC_RN, frac, REC_T)
    return {"neb": neb, "engine": info, "Ic_uA": ic * 1e6, "frac": frac, "beta_c": 2 * math.pi * ic * REC_RN ** 2 * REC_CAP / PHI0, "slips": tot,
            "junction_ns": t_tot * 1e9, "rate_per_ns": rate * 1e-9, "stratonovich_per_ns": r_ah * 1e-9, "kramers_per_ns": r_kr * 1e-9,
            "ratio_to_stratonovich": rate / r_ah, "ratio_to_kramers": rate / r_kr, "pass": bool(0.5 <= rate / r_kr <= 2.0 and 0.5 <= rate / r_ah <= 2.0)}


# ---------------------------------------------------------------------------------------------------------------- main
def run_all(be, neb, quick=False):
    t0 = time.time()
    res = {"neb": neb}
    res["equipartition"] = check_equipartition(be, neb, tstop=150e-9 if quick else 400e-9)
    res["psd"] = check_psd(be, neb, tstop=100e-9 if quick else 200e-9)
    res["step"] = check_step(be, neb, tstop=60e-9 if quick else 100e-9)
    res["ah"] = check_ah(be, neb, tstop=0.8e-6 if quick else 2e-6)
    res["escape"] = check_escape(be, neb, tstop=40e-9 if quick else 60e-9)
    res["pass"] = all(res[k]["pass"] for k in ("equipartition", "psd", "step", "ah", "escape"))
    res["elapsed_s"] = round(time.time() - t0, 1)
    return res


def summarize(res) -> str:
    L = [f"neb {res['neb']:.0e}: {'PASS' if res['pass'] else 'FAIL'} [{res['elapsed_s']} s]"]
    e = res["equipartition"]
    L.append("  equipartition <i^2>L/kT: " + ", ".join(f"{r['R_ohm']:g}/{r['L_pH']:g} {r['reference'][0]} {r['ratio']:.3f} (theory {r['theory']:.3f})" for r in e["rows"])
             + f" -> {'PASS' if e['pass'] else 'FAIL'}")
    for r in res["psd"]["rows"]:
        L.append(f"  psd R {r['R_ohm']:g}: grounded " + ", ".join(f"{k} {v:.3f}" for k, v in r["grounded"].items()) + f"; two-node / grounded {r['two_node_over_grounded']:.4f}"
                 + f" -> {'PASS' if r['pass'] else 'FAIL'}")
    L.append("  step (<i^2>L/kT at each requested step; bar: within 5 % of the 0.1 ps reading): " + "; ".join(f"{r['requested_ps']:g} ps (engine {r['engine'].get('engine_step', 0) * 1e12:g} ps) "
                                                     + " ".join(f"{k} {v:.3f}" for k, v in r.items() if k in ('0g', '0n', '3g')) for r in res["step"]["rows"])
             + f" -> {'PASS' if res['step']['pass'] else 'FAIL'}")
    L.append("  AH <V>/(IcR) sim / closed form: " + "; ".join(f"{r['T_K']:g} K i {r['i']:g}: {r['v_sim_over_IcR']:.4f} / {r['v_ah_over_IcR']:.4f}" for r in res["ah"]["rows"])
             + f" -> {'PASS' if res['ah']['pass'] else 'FAIL'}")
    x = res["escape"]
    L.append(f"  escape {x['Ic_uA']:.1f} uA at {x['frac']} Ic, 65 K: {x['slips']} slips in {x['junction_ns']:.0f} junction-ns = {x['rate_per_ns']:.3f}/ns; "
             f"Stratonovich {x['stratonovich_per_ns']:.3f}/ns, Kramers {x['kramers_per_ns']:.3f}/ns -> {'PASS' if x['pass'] else 'FAIL'}")
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--josim", help="path to josim-cli")
    ap.add_argument("--backend", default=None, choices=("cli", "pyjosim"))
    ap.add_argument("--neb", type=float, default=1e13)
    ap.add_argument("--plateau", action="store_true", help=f"run every check at each neb of {NEBS}")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    if a.backend == "pyjosim" or (a.backend is None and not a.josim):
        be = PyjosimBackend()
    else:
        be = CliBackend(a.josim)
    nebs = NEBS if a.plateau else (a.neb,)
    results = [run_all(be, nb, quick=a.quick) for nb in nebs]
    for r in results:
        print(summarize(r), flush=True)
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(results, fh, indent=1, default=lambda o: bool(o) if isinstance(o, np.bool_) else float(o))
    ok = results[-1]["pass"] if a.plateau else results[0]["pass"]
    if a.plateau:
        # the convention: the lowest neb >= 2e12 at which every bar passes (and every higher neb passes too)
        conv = next((r["neb"] for j, r in enumerate(results) if all(x["pass"] for x in results[j:])), None)
        print(f"convention neb (lowest plateau neb at which every bar passes): {conv}", flush=True)
        ok = conv is not None
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
