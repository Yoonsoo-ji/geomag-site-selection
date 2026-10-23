# -*- coding: utf-8 -*-
"""
2012~2016 절대관측 ↔ 청양(CYG) 상시관측 시계열 대조 — 자료 신뢰성 판정
========================================================================

`parse_kriss_legacy.py` 가 낸 세션표에는 **분 단위 관측시각**이 있다. 같은
시각의 청양 1분 자료를 붙여, **세션 사이에 값이 흔들린 것이 실제 자기장 변화
때문인지 관측 오차 때문인지** 가른다.

핵심 검정은 12-E 가 세운 **β** 다.

    (세션 관측값 − 그 방문의 평균)  =  β · (세션 CYG 값 − 그 방문의 평균)

    β 가 1 에 가까우면 흔들림은 **외부장**이다. 보정하면 걷힌다.
    β 가 0 에 가까우면 흔들림은 **관측 오차**다. 보정해도 잡음만 는다.

⚠️ **청양 자료는 2013-10-30 부터 있다.** 2012 년 세션은 청양으로 대조할 수
   없다. 그 파일에는 **현장 변화계(Flare) 동시관측**이 있으므로 그쪽으로
   따로 본다(`--flare`).

⚠️ **야장 시각은 KST 이고 청양 자료는 UT 다.** `to_utc()` 가 −9h 한다.

    python kriss_cyg_timeseries.py
"""
from __future__ import annotations

import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
SESS = ROOT / "docs" / "data" / "kriss_legacy_sessions.csv"
OUT = ROOT / "docs" / "data" / "kriss_cyg_match.csv"

from lmm_cyg import fetch_day                      # noqa: E402
from lmm_external import derive_dif, to_utc        # noqa: E402

# 성분마다 (관측 열, 청양 열, 단위, 분 환산 배율)
COMPS = (("D", "편각", "D_deg", "′", 60.0),
         ("I", "복각", "I_deg", "′", 60.0),
         ("F", "총자력", "F", " nT", 1.0))


def window(row, kind):
    """세션의 성분별 시간창을 UTC 로 돌려준다. 못 읽으면 None."""
    a, b = ("편각시작", "편각종료") if kind == "D" else ("복각시작", "복각종료")
    if pd.isna(row.get(a)) or pd.isna(row.get(b)) or pd.isna(row.get("날짜")):
        return None
    day = dt.date.fromisoformat(str(row["날짜"])[:10])
    try:
        t0 = dt.time.fromisoformat(str(row[a])[:8])
        t1 = dt.time.fromisoformat(str(row[b])[:8])
    except ValueError:
        return None
    u0, u1 = to_utc(day, t0), to_utc(day, t1)
    if u1 < u0:                                     # 자정을 넘긴 세션
        u1 += pd.Timedelta(days=1)
    if u1 == u0:                                    # 시각이 한 점이면 ±5분
        u0, u1 = u0 - pd.Timedelta(minutes=5), u1 + pd.Timedelta(minutes=5)
    return u0, u1


def load_cyg(days):
    """관측일마다 청양 1분 자료를 읽어 시각 색인으로 돌려준다."""
    got, miss = {}, []
    for d in days:
        try:
            g = fetch_day(dt.date.fromisoformat(d))
        except Exception:
            g = None
        if g is None or not len(g) or not g["X"].notna().any():
            miss.append(d)
            continue
        g = derive_dif(g)
        g["time"] = pd.to_datetime(g["time"], utc=True)
        got[d] = g.set_index("time").sort_index()
    return got, miss


def match(df, cyg):
    """세션마다 같은 시간창의 청양 중앙값·변동폭을 붙인다."""
    keep = ("측점", "연도", "날짜", "회차", "편각", "복각", "총자력",
            "D_판정", "I_판정", "F_판정", "변화계")
    rows = []
    for _, r in df.iterrows():
        d = str(r.get("날짜"))[:10]
        if d not in cyg:
            continue
        g = cyg[d]
        rec = {k: r.get(k) for k in keep}
        for comp, _obs, col, _u, _s in COMPS:
            # F 는 야장이 세션 단위로 적지 않으므로 편각 구간을 쓴다
            w = window(r, "D" if comp == "F" else comp)
            if not w:
                continue
            seg = g.loc[w[0]:w[1], col].dropna()
            if not len(seg):
                continue
            rec[f"CYG_{comp}"] = seg.median()
            rec[f"CYG_{comp}_폭"] = seg.max() - seg.min()
            rec[f"CYG_{comp}_n"] = len(seg)
        rows.append(rec)
    return pd.DataFrame(rows)


def beta_test(M, label="청양"):
    """방문 안에서 중심화한 뒤 관측 편차를 외부장 편차에 회귀한다."""
    print(f"\n■ β 검정 — 세션 간 흔들림이 {label} 변동으로 설명되는가")
    print(f"{'성분':<5}{'n':>4}{'방문':>5}{'관측 산포':>12}{'외부 산포':>12}"
          f"{'상관 r':>9}{'β':>8}   판정")
    out = {}
    for comp, obs, _col, unit, scale in COMPS:
        ext = f"CYG_{comp}"
        if ext not in M.columns:
            continue
        sub = M.dropna(subset=[obs, ext]).copy()
        sub = sub[sub[f"{comp}_판정"] == "정상"]
        # 방문 안에 세션이 둘 이상이어야 중심화가 뜻을 갖는다
        n_in = sub.groupby(["측점", "연도"])[obs].transform("size")
        sub = sub[n_in >= 2]
        if len(sub) < 6:
            print(f"{comp:<5}{len(sub):>4}   표본 부족")
            continue
        g = sub.groupby(["측점", "연도"])
        y = (sub[obs] - g[obs].transform("mean")).to_numpy() * scale
        x = (sub[ext] - g[ext].transform("mean")).to_numpy() * scale
        if np.std(x) < 1e-9:
            print(f"{comp:<5}{len(sub):>4}   외부 변동이 없다")
            continue
        beta = float(np.polyfit(x, y, 1)[0])
        r = float(np.corrcoef(x, y)[0, 1])
        lo, hi = beta - 1.96 * se(x, y, beta), beta + 1.96 * se(x, y, beta)
        verdict = ("외부장이 설명한다" if lo > 0.5 else
                   "관측 오차가 지배한다" if hi < 0.5 else "판정 보류")
        print(f"{comp:<5}{len(sub):>4}{sub.groupby(['측점','연도']).ngroups:>5}"
              f"{np.std(y):>11.2f}{unit}{np.std(x):>11.2f}{unit}"
              f"{r:>9.3f}{beta:>8.3f}   {verdict}  [{lo:+.2f}, {hi:+.2f}]")
        out[comp] = (len(sub), r, beta, lo, hi)
    return out


def se(x, y, beta):
    """회귀 기울기의 표준오차."""
    n = len(x)
    if n < 3:
        return float("inf")
    a = float(np.mean(y) - beta * np.mean(x))
    resid = y - (a + beta * x)
    sxx = float(np.sum((x - np.mean(x)) ** 2))
    if sxx <= 0:
        return float("inf")
    return float(np.sqrt(np.sum(resid ** 2) / (n - 2) / sxx))


def per_visit(M):
    print("\n■ 방문별 — 관측 산포와 같은 구간 청양 변동을 나란히")
    print(f"{'측점':<6}{'연도':<6}{'세션':>4}{'관측 D 산포':>13}"
          f"{'청양 D 산포':>13}{'청양 D 폭':>11}{'설명되는 몫':>12}")
    for (s, y), g in M.groupby(["측점", "연도"]):
        gg = g[(g["D_판정"] == "정상") & g["CYG_D"].notna()]
        if len(gg) < 2:
            continue
        obs = (gg["편각"].max() - gg["편각"].min()) * 60
        ext = (gg["CYG_D"].max() - gg["CYG_D"].min()) * 60
        share = f"{ext / obs * 100:>10.0f} %" if obs > 1e-9 else "          -"
        print(f"{s:<6}{int(y):<6}{len(gg):>4}{obs:>12.2f}′{ext:>12.2f}′"
              f"{gg['CYG_D_폭'].median() * 60:>10.2f}′{share}")


def main():
    df = pd.read_csv(SESS, encoding="utf-8-sig")
    df = df[df["잔재"].isna() | (df["잔재"].astype(str).str.strip() == "")]
    days = sorted({str(d)[:10] for d in df["날짜"].dropna().unique()})

    cyg, miss = load_cyg(days)
    print(f"■ 관측일 {len(days)} · 청양 자료가 붙는 날 {len(cyg)}")
    if miss:
        print(f"  못 붙는 날 {len(miss)} — {', '.join(miss)}")

    M = match(df, cyg)
    if not len(M):
        print("대조된 세션이 없다")
        return 1
    print(f"■ 대조된 세션 {len(M)} · 측점 {M['측점'].nunique()}"
          f" · 연도 {sorted(set(int(y) for y in M['연도']))}")

    beta_test(M)
    per_visit(M)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    M.to_csv(OUT, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {OUT}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
