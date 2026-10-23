# -*- coding: utf-8 -*-
"""
2012~2016 절대관측 신뢰성 판정 — 방문마다 쓸 수 있는지 가른다
================================================================

`parse_kriss_legacy.py` 가 읽은 세션표와, `kriss_cyg_timeseries.py` 가 붙인
같은 시각의 청양 값을 함께 본다. 묻는 것은 하나다 —
**방문(측점,연도)마다 그 성과를 LMM 에 넣을 수 있는가.**

세 단계다.

  1단계 **블런더 솎기** — 방문 안에서 중앙값에서 크게 벗어난 «한 세션»을 찾는다.
     세션 하나가 10° 틀어지는 일이 실제로 있다(봉평 2012 4회차 D=+2.06°).
     그대로 평균을 내면 그 하나가 방문 전체를 끌고 간다.

  2단계 **β 검정** — 블런더를 뺀 나머지 세션의 흔들림을 청양 변동에 회귀한다.
     β 가 1 에 가까우면 그 흔들림은 외부장이고, 세션 평균을 쓰면 걷힌다.
     β 가 0 에 가까우면 관측 오차이고, 평균을 내도 남는다.

  3단계 **방문 등급** — 남은 세션의 산포에서 청양으로 설명되는 몫을 뺀 뒤,
     남는 몫으로 등급을 준다.

⚠️ **블런더 표시는 모든 연도에 건다.** 청양 자료는 2014~2016 분만 있지만
   블런더 판정에는 청양이 필요하지 않다 — 방문 안의 중앙값만 쓴다.
   (처음에 청양이 붙은 행에만 걸었더니 2012·2013 의 블런더를 통째로 놓쳐
   2020 보고서 재현이 19/29 에서 멈췄다.)

⚠️ **블런더를 지우지 않는다.** 표시만 하고 원 세션은 남긴다. 원자료를 버리면
   나중에 되짚을 수 없다(현장 탈락 금지 규약과 같은 어법이다).

    python kriss_reliability.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
SESS = ROOT / "docs" / "data" / "kriss_legacy_sessions.csv"
MATCH = ROOT / "docs" / "data" / "kriss_cyg_match.csv"
OUT_S = ROOT / "docs" / "data" / "kriss_session_flags.csv"
OUT_V = ROOT / "docs" / "data" / "kriss_reliability.csv"

# 성분 (이름, 관측 열, 청양 열, 단위, 분 환산 배율)
COMPS = (("D", "편각", "CYG_D", "′", 60.0),
         ("I", "복각", "CYG_I", "′", 60.0),
         ("F", "총자력", "CYG_F", " nT", 1.0))

# 블런더 문턱 — 방문 안 중앙값에서 이만큼 벗어나면 그 세션을 표시한다.
# 근거: 같은 날 세션 산포가 2019 야장 최대 2.00′ · 2022~25 중앙 1.41′ 다
# (13절). 그 열 배를 넘으면 산포가 아니라 블런더로 본다.
BLUNDER = {"D": 20.0, "I": 20.0, "F": 300.0}

# 방문 등급 경계 — 청양으로 설명되지 않고 남는 편각 산포
GRADE = {"A": 6.0, "B": 20.0}


def flag_blunders(df):
    """방문 안에서 중앙값에서 크게 벗어난 세션을 성분마다 표시한다."""
    df = df.copy()
    for comp, obs, _cyg, _u, scale in COMPS:
        col = f"{comp}_블런더"
        df[col] = False
        v = pd.to_numeric(df[obs], errors="coerce")
        for _, g in df.groupby(["측점", "연도"]):
            ok = g[g[f"{comp}_판정"] == "정상"]
            if len(ok) < 3:
                continue
            vv = v.loc[ok.index]
            dev = (vv - vv.median()).abs() * scale
            df.loc[dev[dev > BLUNDER[comp]].index, col] = True
        # 범위 밖도 못 쓰는 세션이라 같이 표시한다
        df.loc[df[f"{comp}_판정"] != "정상", col] = True
    return df


def report_blunders(df):
    print("■ 1단계 — 블런더로 표시된 세션")
    print(f"{'측점':<6}{'연도':<6}{'회차':>4}{'성분':>5}{'값':>14}"
          f"{'방문 중앙과의 차':>18}  사유")
    n = 0
    for comp, obs, _c, unit, scale in COMPS:
        v = pd.to_numeric(df[obs], errors="coerce")
        for _, g in df.groupby(["측점", "연도"]):
            ok = g[g[f"{comp}_판정"] == "정상"]
            med = v.loc[ok.index].median() if len(ok) else np.nan
            for idx, r in g[g[f"{comp}_블런더"]].iterrows():
                why = ("범위 밖" if r[f"{comp}_판정"] != "정상"
                       else "방문 중앙에서 이탈")
                d = (v.loc[idx] - med) * scale if np.isfinite(med) else np.nan
                ds = f"{d:>+16.2f}{unit}" if np.isfinite(d) else f"{'-':>18}"
                print(f"{r['측점']:<6}{int(r['연도']):<6}{int(r['회차']):>4}"
                      f"{comp:>5}{v.loc[idx]:>14.4f}{ds}  {why}")
                n += 1
    tot = sum(int(df[f"{c}_블런더"].any()) for c, *_ in COMPS)  # noqa: F841
    print(f"  표시 {n} 건 / 세션 {len(df)}")
    return n


def slope(x, y):
    """기울기·상관·95% 구간."""
    n = len(x)
    if n < 4 or np.std(x) < 1e-12:
        return None
    b, a = np.polyfit(x, y, 1)
    resid = y - (a + b * x)
    sxx = float(np.sum((x - np.mean(x)) ** 2))
    se = float(np.sqrt(np.sum(resid ** 2) / (n - 2) / sxx)) if sxx > 0 else np.inf
    r = float(np.corrcoef(x, y)[0, 1])
    return float(b), r, float(b - 1.96 * se), float(b + 1.96 * se)


def beta_test(M, title, drop_blunder=True):
    print(f"\n■ 2단계 — β 검정 ({title})")
    print(f"{'성분':<5}{'n':>4}{'방문':>5}{'관측 산포':>12}{'청양 산포':>12}"
          f"{'r':>8}{'β':>8}{'95% 구간':>20}   판정")
    for comp, obs, cyg, unit, scale in COMPS:
        if cyg not in M.columns:
            continue
        s = M[M[f"{comp}_판정"] == "정상"].dropna(subset=[obs, cyg]).copy()
        if drop_blunder:
            s = s[~s[f"{comp}_블런더"]]
        s = s[s.groupby(["측점", "연도"])[obs].transform("size") >= 2]
        if len(s) < 6:
            print(f"{comp:<5}{len(s):>4}   표본 부족")
            continue
        g = s.groupby(["측점", "연도"])
        y = (s[obs] - g[obs].transform("mean")).to_numpy() * scale
        x = (s[cyg] - g[cyg].transform("mean")).to_numpy() * scale
        res = slope(x, y)
        if res is None:
            continue
        b, r, lo, hi = res
        verdict = ("외부장이 설명한다" if lo > 0.5 else
                   "관측 오차가 지배한다" if hi < 0.5 else "판정 보류")
        print(f"{comp:<5}{len(s):>4}{g.ngroups:>5}{np.std(y):>11.2f}{unit}"
              f"{np.std(x):>11.2f}{unit}{r:>8.3f}{b:>8.3f}"
              f"{f'[{lo:+.2f}, {hi:+.2f}]':>20}   {verdict}")


def grade_visits(M):
    """방문마다 블런더를 뺀 산포에서 청양 몫을 걷어내고 등급을 준다."""
    print("\n■ 3단계 — 방문 등급 (편각 기준)")
    print(f"{'측점':<6}{'연도':<6}{'쓸 세션':>7}{'블런더':>6}{'남은 산포':>11}"
          f"{'청양 산포':>11}{'설명 안 되는 몫':>16}  등급")
    rows = []
    has_cyg = "CYG_D" in M.columns
    for (s, y), g in M.groupby(["측점", "연도"]):
        ok = g[(g["D_판정"] == "정상") & (~g["D_블런더"])]
        nb = int(g["D_블런더"].sum())
        if not len(ok):
            rows.append(dict(측점=s, 연도=int(y), 쓸세션=0, 블런더=nb,
                             등급="X", 사유="쓸 세션이 없다"))
            print(f"{s:<6}{int(y):<6}{0:>7}{nb:>6}   쓸 세션이 없다")
            continue
        obs = ((ok["편각"].max() - ok["편각"].min()) * 60 if len(ok) > 1
               else 0.0)
        ext = np.nan
        if has_cyg and ok["CYG_D"].notna().any() and len(ok) > 1:
            e = ok["CYG_D"].dropna()
            if len(e) > 1:
                ext = float((e.max() - e.min()) * 60)
        # 설명되지 않고 남는 몫 — 직교분해로 본다
        un = float(np.sqrt(max(obs ** 2 - (ext if np.isfinite(ext) else 0.0) ** 2, 0)))
        gr = ("A" if un <= GRADE["A"] else "B" if un <= GRADE["B"] else "C")
        es = f"{ext:>10.2f}′" if np.isfinite(ext) else f"{'-':>10} "
        rows.append(dict(측점=s, 연도=int(y), 쓸세션=len(ok), 블런더=nb,
                         산포=round(obs, 2),
                         청양산포=round(ext, 2) if np.isfinite(ext) else None,
                         미설명=round(un, 2), 등급=gr,
                         편각대표=round(float(ok["편각"].median()), 6),
                         사유=("청양 대조 없음" if not np.isfinite(ext)
                               else "")))
        print(f"{s:<6}{int(y):<6}{len(ok):>7}{nb:>6}{obs:>10.2f}′{es}"
              f"{un:>15.2f}′  {gr}")
    print(f"\n  A 미설명 {GRADE['A']:.0f}′ 이하 · B {GRADE['B']:.0f}′ 이하"
          f" · C 그 밖 · X 쓸 세션 없음")
    print("  ⚠️ 6′ 는 설명자료의 공학적 목표치이고 발주기관이 정한 값이 아니다"
          " — 등급 경계로만 쓴다")
    print("  ⚠️ 청양 대조가 없는 방문(2012·2013)은 산포 전부를 미설명으로 본다"
          " — 그만큼 보수적이다")
    return pd.DataFrame(rows)


def main():
    df = pd.read_csv(SESS, encoding="utf-8-sig")
    df = df[df["잔재"].isna() | (df["잔재"].astype(str).str.strip() == "")]
    df = flag_blunders(df)
    report_blunders(df)

    # 청양 값을 붙인다 — 없는 방문은 등급에서 보수적으로 다룬다
    if MATCH.exists():
        C = pd.read_csv(MATCH, encoding="utf-8-sig")
        cols = ["측점", "연도", "회차"] + [c for c in C.columns
                                          if c.startswith("CYG_")]
        df = df.merge(C[cols].drop_duplicates(subset=["측점", "연도", "회차"]),
                      on=["측점", "연도", "회차"], how="left")
        beta_test(df, "블런더 포함 — 비교용", drop_blunder=False)
        beta_test(df, "블런더 제외 — 정본")
    else:
        print(f"\n⚠ {MATCH.name} 이 없다 — β 검정을 건너뛴다")

    V = grade_visits(df)

    OUT_S.parent.mkdir(parents=True, exist_ok=True)
    keep = (["측점", "연도", "날짜", "회차", "편각", "복각", "총자력",
             "D_판정", "I_판정", "F_판정"]
            + [f"{c}_블런더" for c, *_ in COMPS]
            + [c for c in df.columns if c.startswith("CYG_")])
    df[[c for c in keep if c in df.columns]].to_csv(
        OUT_S, index=False, encoding="utf-8-sig")
    V.to_csv(OUT_V, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {OUT_S} — 세션별 블런더 표시")
    print(f"[저장] {OUT_V} — 방문 등급")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
