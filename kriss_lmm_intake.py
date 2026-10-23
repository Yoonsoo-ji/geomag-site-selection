# -*- coding: utf-8 -*-
"""
2012~2016 절대관측을 LMM 에 어떻게 넣을 것인가 — 투입 판정
============================================================

앞선 네 단계의 결과를 모아 **방문마다 투입 여부를 정한다.**

    ① parse_kriss_legacy.py    세션별 원시 판독값 151
    ② kriss_cyg_timeseries.py  같은 시각 청양 값 (2014~2016)
    ③ kriss_reliability.py     블런더 표시 · β 검정 · 방문 등급
    ④ kriss_reproduce2020.py   2020 보고서 재현 24/29
    ⑤ 이 스크립트              투입 판정

판정은 네 관문을 모두 통과해야 한다.

  관문 1 **좌표** — 파일 좌표가 관측망 좌표와 맞거나, 어긋나더라도 점조서
     이설 이력으로 설명돼야 한다. 설명되지 않으면 **그 자리를 모른다**는
     뜻이라 Regional 적합에 넣을 수 없다.
  관문 2 **블런더를 뺀 세션이 3개 이상** — 적으면 평균이 한 판독에 끌려간다.
  관문 3 **방문 등급 A 또는 B** — C 는 청양으로 설명되지 않는 산포가 20′ 를
     넘는다. 2022~25 성과의 방위 기준 오류와 같은 크기다.
  관문 4 **IGRF 잔차가 상식 범위** — 편각 60′ 초과면 모형이 아니라 자료 문제다.

성분마다 따로 판정한다. **편각은 되고 총자력은 안 되는 방문이 많다.**

⚠️ **투입은 이 스크립트가 하지 않는다.** 판정표만 내고, `lmm_build.py` 에
   거는 것은 동결(Midterm Baseline v1) 해제 뒤의 별도 결정이다.

    python kriss_lmm_intake.py
"""
from __future__ import annotations

import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
FLAGS = ROOT / "docs" / "data" / "kriss_session_flags.csv"
REL = ROOT / "docs" / "data" / "kriss_reliability.csv"
REP = ROOT / "docs" / "data" / "kriss_reproduce2020.csv"
LMM = ROOT / "docs" / "data" / "lmm_model.json"
OUT = ROOT / "docs" / "data" / "kriss_lmm_intake.csv"

import existing_network as EN                        # noqa: E402

MIN_SESSIONS = 3
MAX_RD_ARCMIN = 60.0
MAX_RF_NT = 400.0

# 파일 좌표가 관측망과 멀어도 **점조서 이설 이력으로 설명되는** 방문.
# 폐지점(차수가 다른 옛 표석)에서 관측한 것이라 파일 좌표가 맞다.
RELOCATED = {
    ("삼척", 2012): "폐지된 임계 1-8-3 (2004 폐지) 에서 관측",
    ("영천", 2015): "폐지된 경주 1-22-3 (2005 폐지) 에서 관측",
    ("언양", 2012): "폐지된 양산 1-23-3 (2005 폐지) 에서 관측",
    ("여주", 2012): "2023 양평 이설 전 위치",
    ("여주", 2016): "2016 재설 · 2023 양평 이설 전 위치",
    ("남양", 2016): "2016 장지저수지 이설 — 15-B 에서 확정",
    ("남지", 2015): "2015 망실·재설 (이력 등급 C)",
}

# 파일 머리글이 **다른 측점의 것**인 방문. 관측값은 제 측점의 것이나
# 좌표·전자력은 못 쓴다. 양식을 복사해 쓴 흔적이다(Code 칸과 같은 병).
TEMPLATE_BAD = {
    ("제천", 2015): ("머리글 좌표 34.62731/126.97783 이 **장흥과 0.00 km** 다."
                    " 같은 파일의 전자력도 66,38x nT 로 한반도 값의 1.31배다."),
}

COORD_TOL_KM = 5.0


def haversine(a, b, c, d):
    la, lo, lc, ld = map(np.radians, (a, b, c, d))
    return 6371.0 * 2 * np.arcsin(np.sqrt(
        np.sin((lc - la) / 2) ** 2
        + np.cos(la) * np.cos(lc) * np.sin((ld - lo) / 2) ** 2))


def igrf_resid(lat, lon, day, D, I, F):
    """IGRF-14 잔차. 좌표나 날짜가 없으면 NaN."""
    import ppigrf
    if not all(np.isfinite(v) for v in (lat, lon)):
        return (np.nan,) * 3
    e, n, u = ppigrf.igrf(lon, lat, 0.0, day)
    e, n, u = float(np.ravel(e)[0]), float(np.ravel(n)[0]), float(np.ravel(u)[0])
    d0 = np.degrees(np.arctan2(e, n))
    h = np.hypot(e, n)
    i0 = np.degrees(np.arctan2(-u, h))      # ppigrf 의 u 는 상향이다
    f0 = float(np.sqrt(e * e + n * n + u * u))
    return ((D - d0) * 60 if np.isfinite(D) else np.nan,
            (I - i0) * 60 if np.isfinite(I) else np.nan,
            (F - f0) if np.isfinite(F) else np.nan)


def lmm_sites():
    import json
    if not LMM.exists():
        return set()
    j = json.loads(LMM.read_text(encoding="utf-8"))
    for key in ("sites", "stations", "points"):
        v = j.get(key)
        if isinstance(v, list) and v:
            out = set()
            for s in v:
                nm = (s.get("name") or s.get("site") or s.get("측점")
                      if isinstance(s, dict) else None)
                if nm:
                    out.add(str(nm).split("#")[0].strip())
            if out:
                return out
    return set()


def main():
    S = pd.read_csv(FLAGS, encoding="utf-8-sig")
    V = pd.read_csv(REL, encoding="utf-8-sig")
    rep = (pd.read_csv(REP, encoding="utf-8-sig") if REP.exists()
           else pd.DataFrame())
    sess = pd.read_csv(ROOT / "docs" / "data" / "kriss_legacy_sessions.csv",
                       encoding="utf-8-sig")
    sess = sess[sess["잔재"].isna()
                | (sess["잔재"].astype(str).str.strip() == "")]

    reg = EN.load_register()
    reg["지점명"] = reg["지점명"].astype(str).str.strip()
    R = {r["지점명"]: (float(r["위도"]), float(r["경도"]))
         for _, r in reg.iterrows()}
    inlmm = lmm_sites()

    rows = []
    for _, v in V.iterrows():
        s, y = v["측점"], int(v["연도"])
        g = sess[(sess["측점"] == s) & (sess["연도"] == y)]
        fl = S[(S["측점"] == s) & (S["연도"] == y)]
        lat = float(pd.to_numeric(g["위도"], errors="coerce").median())
        lon = float(pd.to_numeric(g["경도"], errors="coerce").median())
        day = dt.datetime.fromisoformat(str(g["날짜"].iloc[0])[:10])

        # ── 관문 1 좌표 ─────────────────────────────────────────
        why, coord_ok, dist = "", False, np.nan
        if (s, y) in TEMPLATE_BAD:
            why = "머리글이 다른 측점의 것 — " + TEMPLATE_BAD[(s, y)]
        elif not np.isfinite(lat):
            why = "파일에 좌표가 없다"
        elif s not in R:
            why = "관측망 명단에 이름이 없다"
        else:
            dist = float(haversine(lat, lon, *R[s]))
            if dist <= COORD_TOL_KM:
                coord_ok, why = True, f"관측망과 {dist:.2f} km"
            elif (s, y) in RELOCATED:
                coord_ok, why = True, RELOCATED[(s, y)]
            else:
                why = f"관측망과 {dist:.1f} km 어긋나고 이설 근거가 없다"

        # ── 성분별 대표값과 잔차 ────────────────────────────────
        reps, resid = {}, {}
        for comp, col in (("D", "편각"), ("I", "복각"), ("F", "총자력")):
            keep = fl[(fl[f"{comp}_판정"] == "정상")
                      & (~fl[f"{comp}_블런더"].astype(bool))]
            val = pd.to_numeric(keep[col], errors="coerce")
            # ⚠️ 편각 부호 교정 — 한반도 편각은 전부 서편(음)이다. 파일 일부가
            # 부호 없이 적혀 있다(설악 2012 는 +7.9989 로 기재). 2020 보고서도
            # 같은 교정을 했음이 재현에서 확인됐다(kriss_reproduce2020).
            if comp == "D":
                val = -val.abs()
            reps[comp] = float(val.median()) if len(val) else np.nan
            reps[comp + "_n"] = len(keep)
        rD, rI, rF = igrf_resid(lat if coord_ok else np.nan,
                                lon if coord_ok else np.nan, day,
                                reps["D"], reps["I"], reps["F"])
        resid = {"D": rD, "I": rI, "F": rF}

        # ── 관문 2~4 성분별 판정 ────────────────────────────────
        verdict = {}
        for comp, lim in (("D", MAX_RD_ARCMIN), ("I", MAX_RD_ARCMIN),
                          ("F", MAX_RF_NT)):
            if not coord_ok:
                verdict[comp] = "불가 — 좌표"
            elif reps[comp + "_n"] < MIN_SESSIONS:
                verdict[comp] = f"불가 — 세션 {reps[comp + '_n']}"
            elif v["등급"] in ("C", "X"):
                verdict[comp] = f"보류 — 등급 {v['등급']}"
            elif not np.isfinite(resid[comp]):
                verdict[comp] = "불가 — 잔차 계산 불가"
            elif abs(resid[comp]) > lim:
                verdict[comp] = "불가 — 잔차 과대"
            else:
                verdict[comp] = "투입"

        rows.append(dict(
            측점=s, 연도=y, 등급=v["등급"], 쓸세션=int(v["쓸세션"]),
            블런더=int(v["블런더"]), 좌표=("OK" if coord_ok else "NG"),
            관측망거리km=(round(dist, 2) if np.isfinite(dist) else None),
            좌표사유=why,
            편각=round(reps["D"], 6) if np.isfinite(reps["D"]) else None,
            복각=round(reps["I"], 6) if np.isfinite(reps["I"]) else None,
            총자력=round(reps["F"], 2) if np.isfinite(reps["F"]) else None,
            rD분=round(rD, 2) if np.isfinite(rD) else None,
            rI분=round(rI, 2) if np.isfinite(rI) else None,
            rF_nT=round(rF, 1) if np.isfinite(rF) else None,
            D판정=verdict["D"], I판정=verdict["I"], F판정=verdict["F"],
            LMM보유=("있음" if s in inlmm else "없음")))

    T = pd.DataFrame(rows).sort_values(["측점", "연도"])

    print("■ 방문별 투입 판정")
    print(f"{'측점':<6}{'연도':<6}{'등급':>3}{'세션':>4}{'좌표':>5}"
          f"{'rD':>9}{'rI':>9}{'rF':>10}  {'D':<14}{'I':<14}{'F':<14}LMM")
    for _, r in T.iterrows():
        f = lambda v, w, u="": (f"{v:>{w}.1f}{u}" if v is not None  # noqa: E731
                                else f"{'-':>{w}}{u}")
        print(f"{r['측점']:<6}{int(r['연도']):<6}{r['등급']:>3}"
              f"{r['쓸세션']:>4}{r['좌표']:>5}"
              f"{f(r['rD분'], 7, '′')}{f(r['rI분'], 7, '′')}"
              f"{f(r['rF_nT'], 8, 'nT')}  {r['D판정']:<14}"
              f"{r['I판정']:<14}{r['F판정']:<14}{r['LMM보유']}")

    print("\n■ 요약")
    for comp in ("D", "I", "F"):
        c = T[f"{comp}판정"].value_counts()
        print(f"  {comp}  " + " · ".join(f"{k} {v}" for k, v in c.items()))

    ok = T[T["D판정"] == "투입"]
    new = sorted(set(ok.loc[ok["LMM보유"] == "없음", "측점"]))
    have = sorted(set(ok.loc[ok["LMM보유"] == "있음", "측점"]))
    print(f"\n  편각 투입 가능 {len(ok)} 방문 · {ok['측점'].nunique()} 측점")
    print(f"  그중 LMM 에 없는 측점 {len(new)} — {', '.join(new) or '없음'}")
    print(f"  이미 있는 측점 {len(have)} — {', '.join(have) or '없음'}"
          "  (같은 측점의 다른 시기라 Track A 시기간 대조에 쓴다)")

    if len(ok):
        print(f"\n  투입 가능분의 IGRF 편각 잔차 RMS "
              f"{np.sqrt(np.mean(ok['rD분'].dropna() ** 2)):.1f}′")
    okf = T[T["F판정"] == "투입"]
    if len(okf):
        print(f"  총자력 투입 가능 {len(okf)} 방문 · 잔차 RMS "
              f"{np.sqrt(np.mean(okf['rF_nT'].dropna() ** 2)):.0f} nT")

    print("\n  ⚠️ 이 표는 판정이고 투입이 아니다. 동결(Midterm Baseline v1)"
          " 해제 뒤 lmm_build 에 거는 것은 별도 결정이다.")
    T.to_csv(OUT, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {OUT}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
