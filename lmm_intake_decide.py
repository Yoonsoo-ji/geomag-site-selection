# -*- coding: utf-8 -*-
"""
LMM 투입 판정 — 통합 시계열에서 «쓸 수 있는 자료»를 가른다
============================================================

`build_obs_archive.py` 가 모은 다섯 출처의 통합 시계열을 읽어 두 가지를 낸다.

  ① **측점별 시계열 변동 정리** — 같은 표석을 여러 번 쟀을 때 IGRF 잔차가
     얼마나 흔들렸는가. 영년변화가 빠지므로 **0 에 가까워야 정상**이고,
     흔들리면 그 측점·그 방문에 문제가 있다는 뜻이다.

  ② **LMM 초기 구성요소 투입 판정** — 성분마다 네 관문을 모두 통과해야 한다.

관문
----
  1 **중복 제거** — 같은 (측점, 연도) 가 여러 출처에 겹쳐 실린다. 출처마다
    정보량이 달라 **우선순위**로 하나만 고른다. 겹친 출처끼리 값이 어긋나면
    그 자체가 경고다(`출처불일치`).
  2 **좌표** — 그 관측의 좌표를 알아야 한다. 현재점 좌표를 이설 이전 관측에
    붙이면 그 방문의 잔차가 통째로 틀어진다.
  3 **세션** — 블런더를 뺀 유효 세션이 기준 이상. 방문 단위 출처는 세션
    산포를 알 수 없으므로 **보조 등급**으로 내린다.
  4 **잔차** — IGRF-14 잔차가 상식 범위. 넘으면 모형이 아니라 자료 문제다.

⚠️ **투입은 이 스크립트가 하지 않는다.** 판정표만 내고, `lmm_build` 에 거는
   것은 Midterm Baseline v1 동결 해제 뒤의 별도 결정이다.

    python build_obs_archive.py     # ← 먼저
    python lmm_intake_decide.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DOCS = ROOT / "docs" / "data"
SESS = DOCS / "obs_archive_sessions.csv"
VIS = DOCS / "obs_archive_visits.csv"
OUT_T = DOCS / "lmm_intake_decision.csv"
OUT_V = DOCS / "obs_timeseries_summary.csv"

import existing_network as EN                        # noqa: E402

# 출처 우선순위 — 숫자가 작을수록 먼저 쓴다.
#
# 근거는 **정보량**이다. 세션별 원시 판독값과 실제 관측시각·좌표가 있는 쪽이
# 방문 대표값 하나만 있는 쪽보다 낫다 — 블런더를 솎을 수 있고 External 보정량도
# 낼 수 있다. 성과표는 세션이 없지만 좌표가 실측이고 가장 최근 자료다.
SRC_RANK = {"kriss": 1, "survey2225": 2, "fieldbook": 3,
            "register": 4, "legacy1019": 5}
SRC_NOTE = {
    "kriss": "2012~2016 절대관측 원자료 — 세션·시각·참방위각",
    "survey2225": "22~25 성과정리 — 실측 좌표·표고",
    "fieldbook": "2019~2025 야장 — 세션·시각(원시 역산)",
    "register": "지자기 점조서 관측이력 — 1984~ 방문 대표값",
    "legacy1019": "'10~'19 관측현황 — 방문 대표값",
}

MIN_SESSIONS = 3            # 세션 출처의 최소 유효 세션
MAX_RD = 60.0               # 편각·복각 잔차 한계 (분)
MAX_RF = 400.0              # 총자력 잔차 한계 (nT)
BLUNDER = {"D": 20.0, "I": 20.0, "F": 300.0}   # 방문 안 중앙값 이탈 문턱

# 같은 (측점, 연도) 를 두 출처가 담았을 때 「같은 관측」으로 볼 허용차
AGREE = {"D": 3.0, "I": 3.0, "F": 20.0}        # 분 · 분 · nT


def haversine(a, b, c, d):
    la, lo, lc, ld = map(np.radians, (a, b, c, d))
    return 6371.0 * 2 * np.arcsin(np.sqrt(
        np.sin((lc - la) / 2) ** 2
        + np.cos(la) * np.cos(lc) * np.sin((ld - lo) / 2) ** 2))


# 관측 좌표가 그 측점의 관측망 좌표에서 이만큼 넘게 떨어지면 따져 본다.
COORD_TOL_KM = 5.0

# 관측망 좌표와 멀어도 **이설 이력으로 설명되는** 방문. 폐지된 옛 표석에서
# 관측한 것이라 파일 좌표가 맞고, 그 방문을 「다른 측점」으로 다뤄야 한다.
#
# ⚠️ 이 목록에 넣는 것은 «좌표가 설명된다»는 뜻일 뿐 «투입한다»는 뜻이 아니다.
#    이설 전후를 한 측점으로 묶으면 안 되므로(12-A 0단계) 판정은 「이설 이전」로
#    남기고, Regional 적합에 넣으려면 `site_key` 를 분리해야 한다.
RELOCATED = {
    ("삼척", 2012): "폐지된 임계 1-8-3 (2004 폐지)",
    ("영천", 2015): "폐지된 경주 1-22-3 (2005 폐지)",
    ("언양", 2012): "폐지된 양산 1-23-3 (2005 폐지)",
    ("여주", 2012): "2023 양평 이설 전",
    ("여주", 2016): "2016 재설 · 2023 양평 이설 전",
    ("남양", 2012): "2016 장지저수지 이설 전",
    ("남지", 2015): "2015 망실·재설",
}

# 파일 머리글이 **다른 측점의 것**인 방문. 관측값은 제 측점의 것이나
# 좌표·전자력은 못 쓴다. 양식을 복사해 쓴 흔적이다.
TEMPLATE_BAD = {
    ("제천", 2015): ("머리글 좌표가 장흥과 0.00 km · 같은 파일 전자력도 "
                    "66,38x nT (한반도 값의 1.31배)"),
}


def coord_gate(T):
    """
    관측 좌표가 그 측점의 것인지 가른다.

    ⚠️ **이 관문이 없으면 양식 복사 오염이 통과한다.** 제천 2015 는 머리글
       좌표가 장흥의 것인데 편각 잔차만 보면 −11.7′ 로 멀쩡해 보였다(복각이
       182′ 로 걸렸을 뿐이다). 좌표를 보지 않으면 **다른 자리의 자료를 그
       측점에 넣는다.**
    """
    reg = EN.load_register()
    R = {str(r["지점명"]).strip(): (float(r["위도"]), float(r["경도"]))
         for _, r in reg.iterrows()}
    dist, judge, why = [], [], []
    for _, r in T.iterrows():
        s = r["측점"]
        y = int(r["연도"]) if pd.notna(r["연도"]) else None
        la, lo = r["위도"], r["경도"]
        if (s, y) in TEMPLATE_BAD:
            dist.append(np.nan); judge.append("오염")
            why.append("머리글이 다른 측점의 것 — " + TEMPLATE_BAD[(s, y)])
            continue
        if not (np.isfinite(la) and np.isfinite(lo)) or s not in R:
            dist.append(np.nan); judge.append("미상")
            why.append("좌표 없음" if s in R else "관측망 명단에 없음")
            continue
        km = float(haversine(la, lo, *R[s]))
        dist.append(round(km, 2))
        if km <= COORD_TOL_KM:
            judge.append("현재점"); why.append(f"관측망과 {km:.2f} km")
        elif (s, y) in RELOCATED:
            judge.append("이설 이전"); why.append(RELOCATED[(s, y)])
        else:
            judge.append("불일치")
            why.append(f"관측망과 {km:.1f} km · 이설 근거 없음")
    T["관측망거리km"] = dist
    T["좌표판정"] = judge
    T["좌표사유"] = why
    return T


def flag_blunders(S):
    """방문 안에서 중앙값에서 크게 벗어난 세션을 성분마다 표시한다."""
    S = S.copy()
    for comp, col in (("D", "D"), ("I", "I"), ("F", "F")):
        S[f"{comp}_블런더"] = S[f"{comp}_판정"] != "정상"
        v = pd.to_numeric(S[col], errors="coerce")
        for _, g in S.groupby(["측점", "일자", "출처"], dropna=False):
            ok = g[g[f"{comp}_판정"] == "정상"]
            if len(ok) < 3:
                continue
            vv = v.loc[ok.index]
            sc = 60.0 if comp != "F" else 1.0
            dev = (vv - vv.median()).abs() * sc
            S.loc[dev[dev > BLUNDER[comp]].index, f"{comp}_블런더"] = True
    return S


def visit_table(S):
    """세션을 (측점, 일자, 출처) 로 묶어 블런더를 뺀 대표값을 낸다."""
    rows = []
    for (s, day, src), g in S.groupby(["측점", "일자", "출처"], dropna=False):
        rec = dict(측점=s, 일자=day, 출처=src, 단위=g["단위"].iloc[0],
                   연도=(int(g["연도"].iloc[0])
                         if pd.notna(g["연도"].iloc[0]) else None),
                   세션=len(g), 위도=g["위도"].median(), 경도=g["경도"].median(),
                   좌표출처=g["좌표출처"].iloc[0],
                   좌표신뢰=g["좌표신뢰"].iloc[0],
                   일자정밀도=g["일자정밀도"].iloc[0])
        for comp in ("D", "I", "F"):
            ok = g[~g[f"{comp}_블런더"]]
            sc = 60.0 if comp != "F" else 1.0
            rec[comp] = (round(float(pd.to_numeric(ok[comp]).median()), 6)
                         if len(ok) else None)
            rec[f"{comp}_유효"] = len(ok)
            rec[f"{comp}_블런더"] = int(g[f"{comp}_블런더"].sum())
            rec[f"{comp}_산포"] = (
                round(float((pd.to_numeric(ok[comp]).max()
                             - pd.to_numeric(ok[comp]).min()) * sc), 2)
                if len(ok) > 1 else None)
            r = g.loc[~g[f"{comp}_블런더"],
                      "rD분" if comp == "D" else
                      "rI분" if comp == "I" else "rF_nT"].dropna()
            rec[f"r{comp}"] = round(float(r.median()), 2) if len(r) else None
        rows.append(rec)
    return pd.DataFrame(rows)


def dedupe(V):
    """
    같은 (측점, 연도) 를 여러 출처가 담았을 때 하나만 고른다.

    ⚠️ **겹친 쪽을 조용히 버리지 않는다.** 고른 출처와 버린 출처의 값이
       `AGREE` 를 넘게 어긋나면 `출처불일치` 에 적는다 — 같은 관측을 두
       문서가 다르게 적었다는 뜻이고, 그것 자체가 자료 문제의 신호다.
    """
    V = V.copy()
    V["순위"] = V["출처"].map(SRC_RANK).fillna(9)
    out = []
    for (s, y), g in V.groupby(["측점", "연도"], dropna=False):
        g = g.sort_values(["순위", "세션"], ascending=[True, False])
        best = g.iloc[0].to_dict()
        others, bad = [], []
        for _, o in g.iloc[1:].iterrows():
            others.append(f"{o['출처']}")
            for comp, sc in (("D", 60.0), ("I", 60.0), ("F", 1.0)):
                a, b = best.get(comp), o.get(comp)
                if a is None or b is None or not np.isfinite(a) or not np.isfinite(b):
                    continue
                dd = abs(a - b) * sc
                if dd > AGREE[comp]:
                    bad.append(f"{comp} {o['출처']} {dd:.1f}"
                               f"{'′' if comp != 'F' else 'nT'}")
        best["겹친출처"] = ",".join(others)
        best["출처불일치"] = " · ".join(bad)
        out.append(best)
    return pd.DataFrame(out).sort_values(["측점", "연도"])


def timeseries_summary(T):
    """
    ① 측점별 시계열 변동 — 같은 표석을 다시 재면 IGRF 잔차가 같아야 한다.

    ⚠️ 측점 안이라도 **이설 전후는 다른 표석**이다. `좌표신뢰` 가 「이설 이전」인
       방문은 따로 세어 두고 산포 계산에서 뺀다 — 안 그러면 이설이 자료 품질
       문제로 보인다.
    """
    rows = []
    for s, g in T.groupby("측점"):
        g = g.sort_values("연도")
        cur = g[g["좌표신뢰"] != "이설 이전 — 좌표 불확실"]
        rec = dict(측점=s, 방문=len(g), 이설이전=len(g) - len(cur),
                   최초=int(g["연도"].min()) if g["연도"].notna().any() else None,
                   최종=int(g["연도"].max()) if g["연도"].notna().any() else None,
                   출처="+".join(sorted(set(g["출처"]))))
        for comp, unit in (("D", "분"), ("I", "분"), ("F", "nT")):
            v = cur[f"r{comp}"].dropna()
            rec[f"r{comp}_중앙"] = round(float(v.median()), 2) if len(v) else None
            rec[f"r{comp}_폭"] = (round(float(v.max() - v.min()), 2)
                                 if len(v) > 1 else None)
            rec[f"r{comp}_n"] = len(v)
        # 세션 출처에서 같은 날 산포(정밀도)의 중앙값
        sp = cur["D_산포"].dropna()
        rec["같은날_D산포_중앙"] = round(float(sp.median()), 2) if len(sp) else None
        rows.append(rec)
    return pd.DataFrame(rows).sort_values("측점")


def _gap_arcmin(tok):
    """「D register 10.2′」 꼴에서 크기만. nT 는 분으로 환산하지 않는다."""
    import re
    m = re.search(r"([0-9.]+)[ ]*(′|nT)", tok)
    if not m:
        return 0.0
    v = float(m.group(1))
    return v if m.group(2) == "′" else 0.0     # nT 는 이 문턱에 넣지 않는다


def decide(T):
    """② 성분별 네 관문."""
    T = T.copy()
    for comp, lim in (("D", MAX_RD), ("I", MAX_RD), ("F", MAX_RF)):
        verd, why = [], []
        for _, r in T.iterrows():
            v = r.get(comp)
            res = r.get(f"r{comp}")
            n = r.get(f"{comp}_유효") or 0
            if v is None or not np.isfinite(v):
                verd.append("불가"); why.append("값 없음"); continue
            # ── 관문 2 좌표 ───────────────────────────────────────
            cj = r["좌표판정"]
            if cj == "오염":
                verd.append("불가"); why.append("좌표 오염 — " + r["좌표사유"])
                continue
            if cj in ("불일치", "미상"):
                verd.append("불가"); why.append("좌표 — " + r["좌표사유"])
                continue
            if cj == "이설 이전" or r["좌표신뢰"] == "이설 이전 — 좌표 불확실":
                # 좌표는 설명되지만 **다른 표석**이다. 같은 측점으로 묶으면
                # 안 되므로(12-A 0단계) 투입하려면 site_key 를 분리해야 한다.
                verd.append("별도측점"); why.append("이설 이전 — " + r["좌표사유"])
                continue
            # ── 관문 4 잔차 ───────────────────────────────────────
            if res is None or not np.isfinite(res):
                verd.append("불가"); why.append("잔차 계산 불가"); continue
            if abs(res) > lim:
                verd.append("불가")
                why.append(f"잔차 {res:+.0f} 과대"); continue
            # ── 관문 3 세션 ───────────────────────────────────────
            if r["단위"] == "세션" and n < MIN_SESSIONS:
                verd.append("보류"); why.append(f"유효 세션 {int(n)}"); continue
            if r["단위"] == "방문":
                verd.append("보조")
                why.append("방문 대표값 — 세션 산포 미확인"); continue
            # ⚠️ 출처불일치는 **투입을 막지 않는다.** 우리가 고른 쪽이 정보량이
            #    많은 출처이고, 어긋남은 «버린 쪽»의 전사 문제일 수 있다.
            #    다만 30′ 넘게 벌어지면 같은 관측이 아닐 소지가 있어 보류한다.
            big = [x for x in str(r["출처불일치"]).split(" · ")
                   if x and _gap_arcmin(x) > 30.0]
            if big:
                verd.append("보류")
                why.append("출처 간 큰 불일치 — " + ", ".join(big)); continue
            verd.append("투입")
            why.append(("출처 간 작은 불일치 있음" if r["출처불일치"] else ""))
        T[f"{comp}판정"] = verd
        T[f"{comp}사유"] = why
    return T


def main():
    if not SESS.exists():
        print("⚠ obs_archive_sessions.csv 가 없다 —"
              " build_obs_archive.py 를 먼저 돌릴 것")
        return 1
    S = pd.read_csv(SESS, encoding="utf-8-sig")
    S = flag_blunders(S)
    V = visit_table(S)
    T = dedupe(V)
    T = coord_gate(T)
    print(f"■ 통합 {len(S)} 행 → 방문 {len(V)} → 중복 제거 {len(T)}"
          f" · 측점 {T['측점'].nunique()}")

    nb = {c: int(S[f"{c}_블런더"].sum()) for c in ("D", "I", "F")}
    print(f"  블런더 표시 — D {nb['D']} · I {nb['I']} · F {nb['F']} / {len(S)} 행")

    print("\n■ 출처별 채택 — 중복 제거 뒤")
    print(T.groupby("출처").agg(
        방문=("측점", "size"), 측점=("측점", "nunique")).to_string())
    mism = T[T["출처불일치"].astype(str).str.len() > 0]
    print(f"\n■ 출처 간 값이 어긋난 방문 {len(mism)}")
    for _, r in mism.head(14).iterrows():
        print(f"  {r['측점']:<5}{int(r['연도']) if r['연도'] else '':<6}"
              f"{r['출처']:<12}겹침 {r['겹친출처']:<24}{r['출처불일치']}")

    Z = timeseries_summary(T)
    print(f"\n■ ① 측점별 시계열 변동 (IGRF 잔차 · 이설 이전 제외)")
    print(f"{'측점':<6}{'방문':>4}{'이설전':>5}{'기간':>12}"
          f"{'rD 중앙':>9}{'rD 폭':>8}{'rF 중앙':>10}{'rF 폭':>9}  출처")
    for _, r in Z.iterrows():
        f = lambda v, w, u="": (f"{v:>{w}.1f}{u}" if v is not None  # noqa: E731
                                and np.isfinite(v) else f"{'-':>{w}}{u}")
        per = (f"{r['최초']}~{r['최종']}" if r["최초"] else "")
        print(f"{r['측점']:<6}{r['방문']:>4}{r['이설이전']:>5}{per:>12}"
              f"{f(r['rD_중앙'], 8, '′')}{f(r['rD_폭'], 7, '′')}"
              f"{f(r['rF_중앙'], 9)}{f(r['rF_폭'], 9)}  "
              f"{r['출처'].replace('legacy1019', 'L').replace('survey2225', 'S').replace('fieldbook', 'B').replace('register', 'R').replace('kriss', 'K')}")

    T = decide(T)
    print("\n■ ② 투입 판정 요약")
    for comp in ("D", "I", "F"):
        c = T[f"{comp}판정"].value_counts()
        print(f"  {comp}  " + " · ".join(f"{k} {v}" for k, v in c.items()))

    ok = T[T["D판정"] == "투입"]
    aux = T[T["D판정"] == "보조"]
    # ── LMM 현재 측점과 대조 ────────────────────────────────────────
    # ⚠️ `lmm_model.json` 은 **옛 이름**을 담고 있다(「양산」 = 언양). 그대로
    #    비교하면 언양이 「새 측점」으로 세어진다 — SITE_CANON 을 거칠 것.
    import json

    from build_obs_archive import SITE_CANON
    cur = set()
    jf = DOCS / "lmm_model.json"
    if jf.exists():
        j = json.loads(jf.read_text(encoding="utf-8"))
        for k in ("sites", "stations", "points"):
            v = j.get(k)
            if isinstance(v, list) and v:
                for x in v:
                    if isinstance(x, dict):
                        n = x.get("name") or x.get("site") or x.get("측점")
                        if n:
                            n = str(n).split("#")[0].strip()
                            cur.add(SITE_CANON.get(n, n))
                break
    T["LMM보유"] = T["측점"].map(lambda s: "있음" if s in cur else "없음")
    new_s = sorted(set(ok["측점"]) - cur)
    have_s = sorted(set(ok["측점"]) & cur)
    aux_s = sorted(set(aux["측점"]) - cur - set(ok["측점"]))

    print(f"\n  편각 «투입» {len(ok)} 방문 · {ok['측점'].nunique()} 측점")
    print(f"    LMM 현재 {len(cur)} 측점에 없는 측점 {len(new_s)} — "
          + (", ".join(new_s) or "없음"))
    print(f"    이미 있는 측점 {len(have_s)} — " + (", ".join(have_s) or "없음")
          + "  (같은 측점의 다른 시기 — Track A 시기간 대조)")
    print(f"    «보조»로만 더해지는 측점 {len(aux_s)} — "
          + (", ".join(aux_s) or "없음"))
    print(f"  편각 «보조» {len(aux)} 방문 · {aux['측점'].nunique()} 측점"
          f"  (방문 대표값 — 세션 산포를 못 재는 자료)")
    if len(ok):
        print(f"  투입분 IGRF 편각 잔차 RMS "
              f"{np.sqrt(np.mean(ok['rD'].dropna() ** 2)):.1f}′"
              f" · 총자력 "
              f"{np.sqrt(np.mean(T[T['F판정'] == '투입']['rF'].dropna() ** 2)):.0f} nT")

    print("\n■ 투입 가능 방문 목록 (편각)")
    print(f"{'측점':<6}{'연도':<6}{'출처':<12}{'세션':>4}{'rD':>8}{'rI':>8}"
          f"{'rF':>9}  {'I':<6}{'F':<6}")
    for _, r in ok.sort_values(["측점", "연도"]).iterrows():
        f = lambda v, w: (f"{v:>{w}.1f}" if v is not None  # noqa: E731
                          and np.isfinite(v) else f"{'-':>{w}}")
        print(f"{r['측점']:<6}{int(r['연도']):<6}{r['출처']:<12}"
              f"{int(r['D_유효']):>4}{f(r['rD'], 7)}′{f(r['rI'], 7)}′"
              f"{f(r['rF'], 8)}  {r['I판정']:<6}{r['F판정']:<6}")

    T.to_csv(OUT_T, index=False, encoding="utf-8-sig")
    Z.to_csv(OUT_V, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {OUT_T}")
    print(f"[저장] {OUT_V}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
