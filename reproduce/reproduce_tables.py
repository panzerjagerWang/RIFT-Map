#!/usr/bin/env python3
"""Recompute Tables 3 to 6 of the RIFT-Map manuscript from the released measurements.

    python reproduce_tables.py            # print the four tables and compare them with the manuscript
    python reproduce_tables.py --quiet    # print only the comparison

Inputs (data/):
    table3_registration.csv           sizes of the connected models of each registration run
    table4_withheld_intervals.csv     signed relative length discrepancy (%) of the 31 withheld tape intervals,
                                      one column per row of Table 4
    table5_rendering_agreement.csv    image coverage and structural similarity of the six withheld photographs
    table6_simulation.csv             per-seed results of the controlled simulations

Rows and arms are numbered in the order in which the manuscript prints them (row01 is the first row of Table 4,
tape1 and arm1 the first tape and the first arm of Table 6). The manuscript defines them.

The script exits with status 0 when every recomputed value equals the value printed in the manuscript
(expected_tables.json) and with status 1 otherwise.

Tested with Python 3.12, NumPy 1.26.4 and SciPy 1.11.4 and 1.13.1.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

BOOTSTRAP_RESAMPLES = 20000
BOOTSTRAP_SEED = 1337
PROPOSED_ROW = 6                       # row of Table 4 against which the paired tests are made
# The original analysis drew its bootstrap samples from one random stream: first for the leave-one-out intervals,
# then for the check intervals, then for all 31 intervals, and inside each subset in this order of rows.
# Replaying the same order reproduces the printed intervals to the last digit.
STREAM_ROWS = [1, 2, 3, 4, 5, 7, 8, 6, 9, 11, 10]
ROWS = list(range(1, 12))


def read_csv(name):
    with open(DATA / name, newline="") as f:
        return list(csv.DictReader(f))


def rmse(x):
    return float(np.sqrt(np.mean(np.square(x))))


def mae(x):
    return float(np.mean(np.abs(x)))


def fmt_p(p):
    return f"{p:.4f}" if p < 0.001 else (f"{p:.3f}" if p < 0.1 else f"{p:.2f}")


# ------------------------------------------------------------------------------------------------ Table 3
def table3():
    out = []
    for r in read_csv("table3_registration.csv"):
        sizes = [int(v) for v in r["component_sizes"].split()]
        out.append({"row": int(r["row"]), "cameras": r["cameras"], "registered": sum(sizes), "components": len(sizes),
                    "largest": max(sizes), "frames": int(r["frames"])})
    return out


# ------------------------------------------------------------------------------------------------ Table 4
def table4():
    rows = read_csv("table4_withheld_intervals.csv")
    for r in rows:
        r["reading"] = float(r["midpoint_reading_m"])
        for k in ROWS:
            r[k] = float(r[f"row{k:02d}"])
    subsets = {"leave_one_out": [r for r in rows if r["subset"] == "leave_one_out"],
               "check": [r for r in rows if r["subset"] == "check"]}
    subsets["all"] = subsets["leave_one_out"] + subsets["check"]
    rng = np.random.default_rng(BOOTSTRAP_SEED)

    def interval(x, f):
        x = np.asarray(x)
        idx = rng.integers(0, len(x), size=(BOOTSTRAP_RESAMPLES, len(x)))
        v = np.array([f(x[i]) for i in idx])
        return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))

    res = {}
    for name, srows in subsets.items():
        res[name] = {}
        for k in STREAM_ROWS:
            e = np.array([r[k] for r in srows])
            res[name][k] = {"n": len(e), "mae": mae(e), "rmse": rmse(e), "max_abs": float(np.max(np.abs(e))),
                            "rmse_interval": interval(e, rmse), "mae_interval": interval(e, mae)}
        a = np.abs(np.array([r[PROPOSED_ROW] for r in srows]))
        for k in STREAM_ROWS:
            if k == PROPOSED_ROW:
                continue
            b = np.abs(np.array([r[k] for r in srows]))
            interval(b - a, lambda x: float(np.mean(x)))        # drawn by the original analysis, not printed
            res[name][k]["p"] = float(wilcoxon(b, a, alternative="greater").pvalue)
            res[name][k]["smaller"] = int(np.sum(a < b))
    allr = subsets["all"]
    text = {
        # numbers quoted in the text of Section 5.3
        "row01_above_41m": sorted(-r[1] for r in allr if r["reading"] > 41.0),
        "row01_39_to_21m_max_abs": max(abs(r[1]) for r in allr if 21.0 <= r["reading"] <= 39.0),
        "row01_below_steep": sorted(-r[1] for r in allr if r["reading"] < 12.0),
        "interval_41m": {k: next(r[k] for r in allr if r["subset"] == "check" and abs(r["reading"] - 41.1) < 0.06)
                         for k in ROWS},
        "rmse_without_41m": {k: rmse([r[k] for r in allr
                                      if not (r["subset"] == "check" and abs(r["reading"] - 41.1) < 0.06)])
                             for k in (5, 6)},
    }
    return res, text


# ------------------------------------------------------------------------------------------------ Table 5
def table5():
    rows = read_csv("table5_rendering_agreement.csv")
    frames = sorted({int(r["frame"]) for r in rows})
    out = {}
    for fr in frames:
        sub = {r["surface"]: r for r in rows if int(r["frame"]) == fr}
        out[fr] = {"reading": float(sub["sparse"]["reading_m"]),
                   "coverage": {s: 100 * float(sub[s]["coverage"]) for s in ("sparse", "dense", "textured")},
                   "ssim": {s: float(sub[s]["ssim"]) for s in ("sparse", "dense", "textured")},
                   "tape": 100 * float(sub["sparse"]["tape_points_in_mask"])}
    mean = {"coverage": {s: float(np.mean([out[f]["coverage"][s] for f in frames])) for s in ("sparse", "dense", "textured")},
            "ssim": {s: float(np.mean([out[f]["ssim"][s] for f in frames])) for s in ("sparse", "dense", "textured")},
            "tape": float(np.mean([out[f]["tape"] for f in frames]))}
    return out, mean


# ------------------------------------------------------------------------------------------------ Table 6
def table6():
    rows = read_csv("table6_simulation.csv")
    cols = ("scale_error_pct", "camera_error_mm", "check_mae_pct", "check_rmse_pct")
    main = {}
    for tape in ("tape1", "tape2", "tape3"):
        for arm in ("arm1", "arm2", "arm3"):
            sub = [r for r in rows if r["study"] == "main" and r["tape"] == tape and r["arm"] == arm]
            main[(tape, arm)] = {"n": len(sub), **{c: (float(np.mean([float(r[c]) for r in sub])),
                                                      float(np.std([float(r[c]) for r in sub], ddof=1))) for c in cols}}

    def per_seed(tape, arm, col="scale_error_pct"):
        return {int(r["seed"]): float(r[col]) for r in rows
                if r["study"] == "main" and r["tape"] == tape and r["arm"] == arm}

    c1, c3 = per_seed("tape2", "arm1"), per_seed("tape2", "arm3")
    m2, m3 = per_seed("tape3", "arm2"), per_seed("tape3", "arm3")
    sweep = {}
    for amp in sorted({float(r["sweep_value_m"]) for r in rows if r["study"] == "sweep"}):
        sweep[amp] = {arm: float(np.mean([float(r["scale_error_pct"]) for r in rows if r["study"] == "sweep"
                                         and float(r["sweep_value_m"]) == amp and r["arm"] == arm]))
                      for arm in ("arm1", "arm2", "arm3")}
    text = {"tape2_arm3_below_arm1": sum(c3[s] < c1[s] for s in c1), "tape3_arm3_below_arm2": sum(m3[s] < m2[s] for s in m2),
            "sweep": sweep}
    return main, text


# ------------------------------------------------------------------------------------------------ report
def build():
    got = {}
    t3 = table3()
    got["table3"] = {f"row{r['row']}": f"{r['registered']} & {r['components']} & {r['largest']}" for r in t3}
    t4, x4 = table4()
    a = t4["all"]
    got["table4"] = {}
    for k in ROWS:
        v = a[k]
        s = (f"{v['mae']:.2f} & {v['rmse']:.2f} & {v['rmse_interval'][0]:.2f} to {v['rmse_interval'][1]:.2f} & "
             f"{v['max_abs']:.2f}")
        if k != PROPOSED_ROW:
            s += f" & {fmt_p(v['p'])} & {v['smaller']} of {v['n']}"
        got["table4"][f"row{k:02d}"] = s
    got["text_5_3"] = {
        "row01 above the 41 m reading": f"{len(x4['row01_above_41m'])} intervals, {x4['row01_above_41m'][0]:.1f}% to "
                                        f"{x4['row01_above_41m'][-1]:.1f}%",
        "row01 between 39 and 21 m": f"within {x4['row01_39_to_21m_max_abs']:.1f}%",
        "row01 beyond the steep section": f"{len(x4['row01_below_steep'])} intervals, {x4['row01_below_steep'][0]:.1f}% to "
                                          f"{x4['row01_below_steep'][-1]:.1f}%",
        "interval near 41.1 m, rows 3 to 6": f"{max(x4['interval_41m'][k] for k in range(3, 7)):.1f}% to "
                                             f"{min(x4['interval_41m'][k] for k in range(3, 7)):.1f}%",
        "RMSE without that interval": f"row06 {x4['rmse_without_41m'][6]:.2f}%, row05 {x4['rmse_without_41m'][5]:.2f}%",
    }
    t5, m5 = table5()
    got["table5"] = {}
    for fr, v in t5.items():
        got["table5"][str(fr)] = (f"{v['reading']:.1f} & " + " & ".join(f"{v['coverage'][s]:.1f}" for s in ("sparse", "dense", "textured"))
                                  + " & " + " & ".join(f"{v['ssim'][s]:.2f}" for s in ("sparse", "dense", "textured"))
                                  + f" & {v['tape']:.0f}")
    got["table5"]["mean"] = (" & ".join(f"{m5['coverage'][s]:.1f}" for s in ("sparse", "dense", "textured")) + " & "
                             + " & ".join(f"{m5['ssim'][s]:.2f}" for s in ("sparse", "dense", "textured")) + f" & {m5['tape']:.1f}")
    t6, x6 = table6()
    got["table6"] = {}
    for (tape, arm), v in t6.items():
        got["table6"][f"{tape} {arm}"] = (f"{v['scale_error_pct'][0]:.2f} ± {v['scale_error_pct'][1]:.2f} & "
                                          f"{v['camera_error_mm'][0]:.1f} ± {v['camera_error_mm'][1]:.1f} & "
                                          f"{v['check_mae_pct'][0]:.2f} ± {v['check_mae_pct'][1]:.2f} & "
                                          f"{v['check_rmse_pct'][0]:.2f} ± {v['check_rmse_pct'][1]:.2f}")
    got["text_5_5"] = {
        "tape2, seeds in which arm3 is below arm1": f"{x6['tape2_arm3_below_arm1']} of 100",
        "tape3, seeds in which arm3 is below arm2": f"{x6['tape3_arm3_below_arm2']} of 100",
        "arm1 at the sweep value of 2 cm": f"{x6['sweep'][0.02]['arm1']:.2f}%",
        "arm2 and arm3 over the sweep, largest mean": f"{max(x6['sweep'][a][k] for a in x6['sweep'] for k in ('arm2', 'arm3')):.2f}%",
    }
    return got


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quiet", action="store_true", help="print only the comparison with the manuscript")
    ap.add_argument("--write-expected", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()
    got = build()
    exp_path = HERE / "expected_tables.json"
    if args.write_expected:
        exp_path.write_text(json.dumps(got, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        print("wrote", exp_path)
        return 0
    expected = json.loads(exp_path.read_text(encoding="utf-8"))
    bad = 0
    for block, items in expected.items():
        ok = sum(got.get(block, {}).get(k) == v for k, v in items.items())
        bad += len(items) - ok
        print(f"{block:9s} {ok} of {len(items)} entries equal the manuscript")
        for k, v in items.items():
            g = got.get(block, {}).get(k)
            if not args.quiet or g != v:
                mark = "ok  " if g == v else "DIFF"
                print(f"   {mark} {k:48s} {g}" + ("" if g == v else f"    (manuscript: {v})"))
    print("ALL VALUES REPRODUCED" if bad == 0 else f"{bad} VALUES DIFFER")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
