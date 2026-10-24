# -*- coding: utf-8 -*-
"""
점검 대시보드의 자료층 — 측점·방문·세션·플래그를 한 묶음으로
================================================================

`build_obs_archive.py` → `lmm_intake_decide.py` 가 낸 CSV 넷을 읽어
**대시보드가 그대로 쓸 수 있는 형태**로 가공한다. HTML 판과 엑셀 판이
같은 함수를 호출하므로 **두 산출물의 수치가 갈라지지 않는다.**

내는 것
-------
  `stations()`  측점 30 — 좌표·기간·방문수·잔차 변동·플래그 수
  `visits()`    방문 — 연도별 D·I·F 와 잔차, 출처, 원본 파일
  `sessions()`  세션 원자료 — 회차별 판독값
  `flags()`     **확인이 필요한 것** — 종류·심각도·대상·근거·원본 위치

플래그 설계
-----------
점검용이므로 **「무엇이 이상한가」가 아니라 「무엇을 확인해야 하는가」**로
적는다. 각 플래그는 네 가지를 갖춘다 — ① 무엇이 ② 얼마나 ③ 왜 문제인가
④ 어느 파일을 열어 보면 되는가. ④가 없으면 점검을 시작할 수 없다.

⚠️ **자동으로 버리지 않는다.** 이 저장소의 「현장 탈락 금지」 규약과 같다 —
   플래그는 사람이 확인할 목록이지 배제 목록이 아니다.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DOCS = ROOT / "docs" / "data"

S_CSV = DOCS / "obs_archive_sessions.csv"
V_CSV = DOCS / "obs_archive_visits.csv"
T_CSV = DOCS / "lmm_intake_decision.csv"
Z_CSV = DOCS / "obs_timeseries_summary.csv"

import existing_network as EN                        # noqa: E402
from build_obs_archive import resolve_path            # noqa: E402

# 출처를 표에 적을 짧은 이름
SRC_LABEL = {
    "register": "점조서",
    "legacy1019": "관측현황 10~19",
    "survey2225": "성과표 22~25",
    "kriss": "절대관측 원자료 12~16",
    "fieldbook": "야장 19~25",
}

# 플래그 심각도 — 숫자가 클수록 먼저 본다
SEV = {"확인필요": 3, "주의": 2, "참고": 1}

# 재방문 잔차 폭 문턱(분). 13절 재방문 감사가 RMS 33.7′ 를 「불량」으로
# 판정했고 법정 §20 정수차 한계가 30′ 다. 그 둘을 넘는 60′ 를 1차 선으로,
# 20′ 를 2차 선으로 둔다.
SPREAD_HI, SPREAD_MID = 60.0, 20.0


def _load():
    for f in (S_CSV, V_CSV, T_CSV, Z_CSV):
        if not f.exists():
            raise SystemExit(
                f"⚠ {f.name} 이 없다 — build_obs_archive.py →"
                f" lmm_intake_decide.py 순으로 먼저 돌릴 것")
    S = pd.read_csv(S_CSV, encoding="utf-8-sig")
    V = pd.read_csv(V_CSV, encoding="utf-8-sig")
    T = pd.read_csv(T_CSV, encoding="utf-8-sig")
    Z = pd.read_csv(Z_CSV, encoding="utf-8-sig")
    for d in (S, V, T, Z):
        d.columns = [str(c).strip() for c in d.columns]
    return S, V, T, Z


def _f(v, nd=2):
    try:
        x = float(v)
        return None if not np.isfinite(x) else round(x, nd)
    except Exception:
        return None


def build():
    S, V, T, Z = _load()
    reg = EN.load_register()
    R = {str(r["지점명"]).strip(): r for _, r in reg.iterrows()}

    # ── 세션 원자료 ────────────────────────────────────────────────
    S = S.copy()
    S["출처명"] = S["출처"].map(lambda k: SRC_LABEL.get(k, k))
    # ⚠️ CSV 에는 파일 **이름만** 실려 있다(절대 경로는 커밋하지 않는다).
    #    점검에는 전체 경로가 필요하므로 여기서 되만든다 — 로컬 전용이다.
    S["원본파일"] = [resolve_path(a, b) for a, b in
                    zip(S["출처"], S["원본파일"].fillna(""))]

    # ── 방문 — 중복 제거 전 전부. 출처가 겹치면 그대로 둔다 ─────────
    #    (점검 대상이므로 «어느 출처가 무엇을 적었는가»를 다 봐야 한다)
    V = V.copy()
    V["출처명"] = V["출처"].map(lambda k: SRC_LABEL.get(k, k))
    prov = (S.groupby(["측점", "일자", "출처"], dropna=False)
             .agg(원본파일=("원본파일", "first"),
                  원본시트=("원본시트", "first"),
                  일자정밀도=("일자정밀도", "first"),
                  좌표신뢰=("좌표신뢰", "first")).reset_index())
    V = V.merge(prov, on=["측점", "일자", "출처"], how="left")

    # 판정 결과(중복 제거 뒤)를 방문에 붙인다 — 없으면 「미판정」
    jt = T[["측점", "일자", "출처", "D판정", "I판정", "F판정", "D사유",
            "좌표판정", "관측망거리km", "출처불일치", "LMM보유"]]
    V = V.merge(jt, on=["측점", "일자", "출처"], how="left")
    for c in ("D판정", "I판정", "F판정"):
        V[c] = V[c].fillna("중복(미채택)")

    F = flags(S, V, T, Z)

    # ── 측점 ───────────────────────────────────────────────────────
    rows = []
    for s, g in V.groupby("측점"):
        z = Z[Z["측점"] == s]
        z = z.iloc[0] if len(z) else {}
        rr = R.get(s)
        nf = F[F["측점"] == s]
        rows.append(dict(
            측점=s,
            위도=_f(rr["위도"], 6) if rr is not None else None,
            경도=_f(rr["경도"], 6) if rr is not None else None,
            소재지=(str(rr.get("소재지", "")) if rr is not None else ""),
            점의부호=(str(rr.get("지점코드", "")) if rr is not None else ""),
            방문=len(g),
            최초=int(g["연도"].min()) if g["연도"].notna().any() else None,
            최종=int(g["연도"].max()) if g["연도"].notna().any() else None,
            출처수=g["출처"].nunique(),
            rD_중앙=_f(z.get("rD_중앙") if len(z) else None),
            rD_폭=_f(z.get("rD_폭") if len(z) else None),
            rF_중앙=_f(z.get("rF_중앙") if len(z) else None, 0),
            rF_폭=_f(z.get("rF_폭") if len(z) else None, 0),
            같은날산포=_f(z.get("같은날_D산포_중앙") if len(z) else None),
            투입D=int((g["D판정"] == "투입").sum()),
            플래그=len(nf),
            확인필요=int((nf["심각도"] == "확인필요").sum()) if len(nf) else 0,
            LMM보유=("있음" if (g["LMM보유"] == "있음").any() else "없음")))
    ST = pd.DataFrame(rows).sort_values("측점").reset_index(drop=True)
    return ST, V, S, F


def _site_files(S, site):
    f = S[S["측점"] == site]["원본파일"].dropna().astype(str)
    f = sorted({x for x in f if x})
    return " · ".join(f) if len(f) <= 3 else f"{len(f)}개 파일 — ⑥ 시트 참조"


def flags(S, V, T, Z):
    """확인이 필요한 항목을 한 표로. 종류마다 근거와 원본 위치를 붙인다."""
    out = []

    # ⚠️ **플래그마다 원본 파일을 반드시 붙인다.** 「어느 파일을 열면 되는가」가
    #    없으면 점검을 시작할 수 없다. 판정표(T)에는 그 열이 없으므로
    #    세션표에서 (측점, 일자, 출처) 로 끌어온다.
    PROV = {}
    for (st_, dy_, sr_), g in S.groupby(["측점", "일자", "출처"], dropna=False):
        f_ = g["원본파일"].dropna().astype(str)
        f_ = f_[f_.str.len() > 0]
        if len(f_):
            PROV[(st_, str(dy_), sr_)] = f_.iloc[0]

    def prov(r):
        """그 방문의 원본 파일. 없으면 같은 측점·출처의 아무 파일이라도."""
        k = (r.get("측점"), str(r.get("일자")), r.get("출처"))
        if k in PROV:
            return PROV[k]
        for (st_, _d, sr_), v in PROV.items():
            if st_ == r.get("측점") and sr_ == r.get("출처"):
                return v
        return ""

    def add(kind, sev, site, year, what, why, where, src="", detail=""):
        out.append(dict(종류=kind, 심각도=sev, 측점=site,
                        연도=(int(year) if pd.notna(year) else None),
                        내용=what, 왜=why, 원본파일=where, 출처=src,
                        상세=detail))

    # ① 재방문 잔차가 크게 흔들린다 — 같은 표석이면 0 이어야 한다
    for _, r in Z.iterrows():
        w = r.get("rD_폭")
        if w is None or not np.isfinite(w):
            continue
        if w > SPREAD_HI:
            add("재방문 편각 불일치", "확인필요", r["측점"], None,
                f"편각 잔차 폭 {w:.0f}′",
                "같은 표석을 다시 재면 IGRF 잔차가 같아야 한다. "
                "방문마다 달라지는 방위표지 기준이 주된 후보다.",
                _site_files(S, r["측점"]), "",
                f"기간 {r.get('최초')}~{r.get('최종')} · 방문 {r.get('방문')}")
        elif w > SPREAD_MID:
            add("재방문 편각 불일치", "주의", r["측점"], None,
                f"편각 잔차 폭 {w:.0f}′",
                "법정 §20 정수차 한계 30′ 수준이다. 방위표지 이력을 본다.",
                _site_files(S, r["측점"]), "",
                f"기간 {r.get('최초')}~{r.get('최종')}")
        wf = r.get("rF_폭")
        if wf is not None and np.isfinite(wf) and wf > 1000:
            add("재방문 총자력 불일치", "확인필요", r["측점"], None,
                f"총자력 잔차 폭 {wf:,.0f} nT",
                "단위·기기·기록 위치 가운데 무엇인지 아직 가리지 못했다.",
                _site_files(S, r["측점"]), "", "")

    # ② 값이 한반도 범위 밖 — 전사·단위 오류를 먼저 의심한다
    for comp, nm, unit in (("D", "편각", "°"), ("I", "복각", "°"),
                           ("F", "총자력", " nT")):
        bad = S[S[f"{comp}_판정"] != "정상"]
        for _, r in bad.iterrows():
            v = r.get(comp)
            add("범위 밖 값", "확인필요", r["측점"], r.get("연도"),
                f"{nm} {v:,.4f}{unit}" if pd.notna(v) else f"{nm} 결측",
                "한반도에서 나올 수 없는 값이다. 원본 칸을 직접 볼 것.",
                r.get("원본파일", ""), r.get("출처명", ""),
                (f"{r.get('일자')} 회차 {int(r['회차'])}"
                 if pd.notna(r.get("회차")) else f"{r.get('일자')}"))

    # ③ 세션 블런더 — 방문 안 중앙값에서 크게 벗어난 «한 세션»
    for comp, nm in (("D", "편각"), ("I", "복각"), ("F", "총자력")):
        col = f"{comp}_블런더"
        if col not in T.columns:
            continue
        for _, r in T[T[col].fillna(0) > 0].iterrows():
            add("세션 블런더", "주의", r["측점"], r.get("연도"),
                f"{nm} 블런더 {int(r[col])}세션 / {int(r['세션'])}",
                "한 세션만 크게 벗어났다. 그 회차의 원시 판독을 확인한다.",
                prov(r), SRC_LABEL.get(r["출처"], r["출처"]),
                f"유효 {int(r.get(comp + '_유효') or 0)}세션")

    # ④ 좌표 — 오염·불일치는 다른 자리의 자료를 그 측점에 넣게 된다
    for _, r in T.iterrows():
        cj = str(r.get("좌표판정", ""))
        if cj in ("오염", "불일치", "미상"):
            add("좌표", "확인필요", r["측점"], r.get("연도"),
                f"좌표 {cj}", str(r.get("좌표사유", "")),
                prov(r), SRC_LABEL.get(r["출처"], r["출처"]),
                (f"관측망과 {r['관측망거리km']:.1f} km"
                 if pd.notna(r.get("관측망거리km")) else ""))
        elif cj == "이설 이전":
            add("이설 이전 관측", "참고", r["측점"], r.get("연도"),
                "폐지된 옛 표석에서 관측",
                "같은 측점으로 묶으면 안 된다. 측점 식별자를 분리해야 쓴다.",
                prov(r), SRC_LABEL.get(r["출처"], r["출처"]),
                str(r.get("좌표사유", "")))

    # ⑤ 출처 두 곳이 같은 관측을 다르게 적었다
    #
    # ⚠️ **한 덩어리로 묶지 않는다.** 한 방문에서 D·I·F 가 각각 다른 출처와
    #    다른 크기로 어긋나는데, 묶어 놓으면 「2′ 차이」와 「2,000 nT 차이」가
    #    같은 줄에 같은 등급으로 앉는다(처음에 192건이 전부 «주의» 였다).
    #    성분마다 한 줄로 쪼개고 **크기로 등급**을 매긴다.
    import re
    for _, r in T[T["출처불일치"].astype(str).str.len() > 0].iterrows():
        for tok in str(r["출처불일치"]).split(" · "):
            m = re.match(r"([DIF])\s+(\S+)\s+([0-9.]+)(′|nT)", tok.strip())
            if not m:
                continue
            comp, other, val, unit = m.group(1), m.group(2), float(m.group(3)), m.group(4)
            nm = {"D": "편각", "I": "복각", "F": "총자력"}[comp]
            if unit == "′":
                sev = ("확인필요" if val > 30 else "주의" if val > 10 else "참고")
                shown = f"{val:,.1f}′"
            else:
                sev = ("확인필요" if val > 500 else "주의" if val > 100 else "참고")
                shown = f"{val:,.0f} nT"
            add("출처 간 값 불일치", sev, r["측점"], r.get("연도"),
                f"{nm} {shown} 차이",
                "같은 관측을 두 문서가 다르게 적었다. 어느 쪽이 전사 오류인지 "
                "원본에서 가린다.",
                prov(r), f"{SRC_LABEL.get(r['출처'], r['출처'])} ↔ "
                         f"{SRC_LABEL.get(other, other)}", "")

    # ⑥ 잔차 과대 — 모형이 아니라 자료 문제다
    for _, r in T.iterrows():
        for comp, nm, lim, unit in (("D", "편각", 60.0, "′"),
                                    ("I", "복각", 60.0, "′"),
                                    ("F", "총자력", 400.0, " nT")):
            v = r.get(f"r{comp}")
            if v is None or not np.isfinite(v) or abs(v) <= lim:
                continue
            add("IGRF 잔차 과대", "확인필요", r["측점"], r.get("연도"),
                f"{nm} 잔차 {v:+,.1f}{unit} (한계 {lim:.0f}{unit})",
                "모형 문제가 아니라 자료 문제다. 관측값·좌표·날짜를 본다.",
                prov(r), SRC_LABEL.get(r["출처"], r["출처"]), "")

    # ⑦ 날짜를 연도까지만 안다 — 영년변화를 반년 어긋나게 건다
    for (s, src), g in S[S["일자정밀도"] == "연"].groupby(["측점", "출처"]):
        add("관측일자 미상", "참고", s, g["연도"].min(),
            f"{len(g)}행이 연도까지만 기재",
            "7월 1일로 대체하므로 영년변화가 최대 반년 어긋난다"
            "(편각 약 ±2.5′).",
            g["원본파일"].iloc[0], SRC_LABEL.get(src, src),
            f"연도 {sorted(set(int(y) for y in g['연도'].dropna()))}")

    # ⑧ 편각이 양수로 적혀 있었다 — 한반도는 전부 서편(음)이다
    for (s, src), g in S[S["D부호교정"].astype(str).str.lower()
                         .isin(["true", "1"])].groupby(["측점", "출처"]):
        add("편각 부호 교정", "참고", s, g["연도"].min(),
            f"{len(g)}행이 양수로 기재 — 음수로 교정",
            "그대로 쓰면 IGRF 잔차가 1,000′ 가까이 나온다. 원본 기재가 "
            "부호를 뺀 것인지 확인한다.",
            g["원본파일"].iloc[0], SRC_LABEL.get(src, src), "")

    # ⑨ 좌표를 교정한 측점 — 점조서 소재지와 실제 표석 자리가 어긋난다
    #
    # ⚠️ 대시보드가 보이는 소재지는 **점조서 기재 그대로**이고 좌표는 교정값
    #    이다. 둘이 다른 곳을 가리키는 것이 정상 상태이므로 그대로 두면
    #    보는 사람이 혼동한다 — 플래그로 명시한다.
    for nm, ov in getattr(EN, "COORD_OVERRIDE", {}).items():
        add("좌표 교정됨", "참고", nm, None,
            f"표석 좌표를 {ov['lat']:.6f} / {ov['lon']:.6f} 로 교정",
            str(ov.get("why", "")),
            "", "existing_network.COORD_OVERRIDE",
            "화면의 소재지는 점조서 기재 그대로다 — 좌표와 다른 곳을 "
            "가리키는 것이 현재 상태다")

    F = pd.DataFrame(out)
    if len(F):
        F["정렬"] = F["심각도"].map(SEV).fillna(0)
        F = F.sort_values(["정렬", "측점", "연도"],
                          ascending=[False, True, True]).drop(columns="정렬")
    return F.reset_index(drop=True)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    ST, V, S, F = build()
    print(f"■ 측점 {len(ST)} · 방문 {len(V)} · 세션 {len(S)} · 플래그 {len(F)}")
    print("\n■ 플래그 종류별")
    print(F.groupby(["심각도", "종류"]).size().to_string())
    print(f"\n■ 플래그가 가장 많은 측점")
    print(ST.sort_values("플래그", ascending=False)
          [["측점", "방문", "플래그", "확인필요", "rD_폭", "rF_폭"]]
          .head(10).to_string(index=False))
