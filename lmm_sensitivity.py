# -*- coding: utf-8 -*-
"""
측점 수 대비 자유도 민감도 — 자문 요구 항목(「15점 최소 구성」)
=================================================================

무작위 부분표본을 반복 추출해 측점 수별 LOO 를 잰다. 적합·검증은 생산 설정
(lmm_build 의 현재 상수) 그대로 쓴다 — 지각 벡터·External 포함.

    python lmm_sensitivity.py [--reps 12]

⚠ 표본을 줄이면 inlier 판정도 함께 흔들린다. 여기의 LOO 는 각 부분표본이
   스스로 고른 inlier 위에서 계산되므로, 절대값보다 **무너지는 지점**을 보라.
"""
import argparse
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd


def main(reps=12, seed=20260819):
    sys.stdout.reconfigure(encoding="utf-8")
    import lmm_build as LB

    pts = LB.load_all_points(include_2019=LB.INCLUDE_2019)
    res = LB.igrf_residuals(pts)
    sites0 = LB.aggregate_sites(pts, res)
    lons, lats, grid = LB.load_kigam_grid()
    crustal = LB.CrustalGrid(lons, lats, grid)
    cvec = None
    if LB.CRUSTAL_VECTOR:
        cvec = LB.CrustalVector(lons, lats, *LB.crustal_vector(lons, lats, grid))
    sites0 = LB.attach_crustal_di(sites0, cvec)

    n_all = len(sites0)
    rng = np.random.default_rng(seed)
    rows = []
    for n in [n_all, 14, 12, 10, 8]:
        if n > n_all:
            continue
        for deg in (1, 2):
            n_terms = (deg + 1) * (deg + 2) // 2
            if n < n_terms + 3:
                rows.append(dict(n=n, degree=deg, D=np.nan, I=np.nan, F=np.nan,
                                 note="자유도 부족"))
                continue
            got = {"D": [], "I": [], "F": []}
            for _ in range(1 if n == n_all else reps):
                idx = (np.arange(n_all) if n == n_all
                       else rng.choice(n_all, n, replace=False))
                sub = sites0.iloc[sorted(idx)].reset_index(drop=True)
                try:
                    cv = LB.loo_cv(sub, crustal, deg)
                except Exception:      # noqa: BLE001
                    continue
                for k in got:
                    got[k].append(cv[k])
            rows.append(dict(n=n, degree=deg, note="",
                             **{k: float(np.median(v)) if v else np.nan
                                for k, v in got.items()}))
            print(f"  측점 {n:2d} · degree {deg}: "
                  f"LOO D {rows[-1]['D']:.4f}°  I {rows[-1]['I']:.4f}°  "
                  f"F {rows[-1]['F']:.1f} nT")

    df = pd.DataFrame(rows)
    piv = df.pivot(index="n", columns="degree", values=["D", "F"])
    print("\n" + "=" * 62)
    print(f"{'측점':>4} {'1차 LOO D':>11} {'2차 LOO D':>11} {'1차 F':>9} {'2차 F':>9}")
    for n in sorted(df.n.unique(), reverse=True):
        r = piv.loc[n]
        print(f"{n:4d} {r[('D',1)]:11.4f} {r[('D',2)]:11.4f} "
              f"{r[('F',1)]:9.1f} {r[('F',2)]:9.1f}")
    print("=" * 62)

    from pathlib import Path
    import datetime as dt
    out = Path(__file__).parent / "docs" / "output"
    out.mkdir(parents=True, exist_ok=True)
    p = out / f"{dt.datetime.now():%Y%m%d_%H%M%S}_측점수_민감도.csv"
    df.to_csv(p, index=False, encoding="utf-8-sig")
    # 보고서·발표자료가 하드코딩 없이 읽도록 JSON 도 남긴다
    j = Path(__file__).parent / "docs" / "data" / "lmm_sensitivity.json"
    df.to_json(j, orient="records", force_ascii=False)
    print(f"[저장] {p}\n[저장] {j}")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=12)
    main(ap.parse_args().reps)
