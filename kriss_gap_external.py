# -*- coding: utf-8 -*-
"""
재현 안 되는 쌍이 청양 보정으로 설명되는가
==========================================

2020 보고서는 「IGRF 자료와 청양 상시관측소의 데이터를 참조하였다」고 적었다.
`kriss_reproduce2020.py` 가 어떤 세션 조합으로도 재현하지 못한 쌍이 **청양
보정을 걸면 설명되는지** 본다.

방법은 `lmm_external.py` 와 같다 — 관측일 ±10일 · Kp≤2 · 지방시 23~03시의
청양 중앙값을 기준선으로 잡고, 세션 시간창의 청양 편차를 보정량으로 본다.
그 보정량이 **메워야 할 간격과 크기·부호가 맞는지**가 판정이다.

⚠️ **보정량이 간격보다 작으면 청양으로 설명되지 않는다.** 그때는 원인을
   다른 데서 찾아야 한다(방위 기준·전사·기준시점 환산).

⚠️ 2012·2013 구간에는 청양 자료가 없다(2013-10-30 부터이고 그 해 11·12월은
   HAPI 가 거부한다). 그 쌍은 **검정할 수 없다**고 적고 넘어간다.

    python kriss_gap_external.py
"""
from __future__ import annotations

import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
REP = ROOT / "docs" / "data" / "kriss_reproduce2020.csv"
SESS = ROOT / "docs" / "data" / "kriss_legacy_sessions.csv"

from lmm_cyg import fetch_day, fetch_kp, kp_for          # noqa: E402
from lmm_external import derive_dif, to_utc              # noqa: E402

BASE_DAYS = 10          # 기준선 창 (관측일 ±)
QUIET_KP = 2.0
NIGHT = (23, 3)         # 지방시


def quiet_base(day, kp):
    """관측일 ±10일의 정온야간 청양 편각 중앙값."""
    vals = []
    for k in range(-BASE_DAYS, BASE_DAYS + 1):
        d = day + dt.timedelta(days=k)
        try:
            g = fetch_day(d)
        except Exception:
            continue
        if g is None or not len(g):
            continue
        g = derive_dif(g).dropna(subset=["X", "Y"])
        if not len(g):
            continue
        t = pd.to_datetime(g["time"], utc=True)
        loc = t.dt.tz_convert("Asia/Seoul").dt.hour
        night = (loc >= NIGHT[0]) | (loc < NIGHT[1])
        quiet = kp_for(t, kp).to_numpy() <= QUIET_KP
        sel = g[night.to_numpy() & quiet]
        if len(sel):
            vals.append(sel["D_deg"])
    if not vals:
        return None
    return float(pd.concat(vals).median())


def main():
    if not REP.exists():
        print("⚠ kriss_reproduce2020.csv 가 없다 — 먼저 그 스크립트를 돌릴 것")
        return 1
    T = pd.read_csv(REP, encoding="utf-8-sig")
    df = pd.read_csv(SESS, encoding="utf-8-sig")
    df = df[df["잔재"].isna() | (df["잔재"].astype(str).str.strip() == "")]

    gaps = T[T["편차초"].isna()].copy()
    print(f"■ 재현 안 되는 쌍 {len(gaps)}")

    kp = fetch_kp()
    print(f"{'측점':<6}{'연도':<6}{'원자료':>11}{'보고서':>11}"
          f"{'메울 간격':>11}{'청양 보정량':>13}{'보정량/간격':>12}  판정")
    out = []
    for _, r in gaps.iterrows():
        g = df[(df["측점"] == r["측점"]) & (df["연도"] == r["연도"])]
        d2 = -pd.to_numeric(g["편각"], errors="coerce").abs()
        raw = float(d2.mean())
        need = (r["보고서_편각"] - raw) * 60          # 분
        day = dt.date.fromisoformat(str(g["날짜"].iloc[0])[:10])

        try:
            got = fetch_day(day)
        except Exception:
            got = None
        if got is None or not len(got) or not got["X"].notna().any():
            print(f"{r['측점']:<6}{int(r['연도']):<6}{raw:>11.5f}"
                  f"{r['보고서_편각']:>11.5f}{need:>+10.2f}′"
                  f"{'청양 자료 없음':>13}{'-':>12}  검정 불가")
            out.append(dict(측점=r["측점"], 연도=int(r["연도"]),
                            간격분=round(need, 2), 보정량분=None,
                            판정="검정 불가 — 청양 자료 없음"))
            continue

        base = quiet_base(day, kp)
        gg = derive_dif(got)
        gg["time"] = pd.to_datetime(gg["time"], utc=True)
        gg = gg.set_index("time").sort_index()

        devs = []
        for _, s in g.iterrows():
            if pd.isna(s.get("편각시작")) or pd.isna(s.get("편각종료")):
                continue
            try:
                t0 = to_utc(day, dt.time.fromisoformat(str(s["편각시작"])[:8]))
                t1 = to_utc(day, dt.time.fromisoformat(str(s["편각종료"])[:8]))
            except ValueError:
                continue
            seg = gg.loc[min(t0, t1):max(t0, t1), "D_deg"].dropna()
            if len(seg) and base is not None:
                devs.append(float(seg.median()) - base)
        if not devs:
            print(f"{r['측점']:<6}{int(r['연도']):<6}   시간창이 안 맞는다")
            continue
        corr = -float(np.mean(devs)) * 60      # 보정은 편차를 «빼는» 것이다
        ratio = corr / need if abs(need) > 1e-9 else np.nan
        verdict = ("청양이 설명한다" if 0.5 <= ratio <= 1.5 else
                   "부호가 반대다" if ratio < 0 else
                   "보정량이 너무 작다" if abs(ratio) < 0.5 else
                   "보정량이 너무 크다")
        print(f"{r['측점']:<6}{int(r['연도']):<6}{raw:>11.5f}"
              f"{r['보고서_편각']:>11.5f}{need:>+10.2f}′{corr:>+12.2f}′"
              f"{ratio:>12.2f}  {verdict}")
        out.append(dict(측점=r["측점"], 연도=int(r["연도"]),
                        간격분=round(need, 2), 보정량분=round(corr, 2),
                        비=round(float(ratio), 3), 판정=verdict))

    print("\n  보정량은 «기준선 대비 세션 편차의 평균을 뺀 값»이다."
          " 간격과 비가 1 에 가까우면 청양으로 설명된다.")
    print("  ⚠️ 청양 한 곳의 편차를 전국에 같이 적용하는 균일 근사다"
          " — 거리가 멀면 성립하지 않는다(12절).")
    P = ROOT / "docs" / "data" / "kriss_gap_external.csv"
    pd.DataFrame(out).to_csv(P, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {P}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
