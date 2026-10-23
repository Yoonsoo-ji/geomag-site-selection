# -*- coding: utf-8 -*-
"""
2012~2016 지자기 절대관측 원자료 파서 — 주성중 박사 제공분
================================================================

`20260908_주성중 박사_2020년 사용 데이터` 는 **2020 연구가 실제로 쓴 원 관측
자료**다. 성과값만 있는 `'10~'19년 관측현황.xls` 와 달리 **세션별 원시 판독값
(ED/WD/EU/WU · SU/ND/SD/NU) · 관측 시각 · 참방위각 · 변화계 동시관측**이
그대로 들어 있다.

⚠️ **이 자료가 푸는 것** — 13절이 「자료만으로 구분 불가」라 적어 둔 세 가지:

  ① **방위 기준** — 세션마다 `FM True`(참방위각)가 적혀 있어, 재방문 사이에
     방위각이 얼마나 옮겨 갔는지를 **2012~16 구간에서 직접** 잴 수 있다.
  ② **관측 시각** — 2022~25 야장과 같은 분 단위 시각이 있어 External 보정량을
     낼 수 있다.
  ③ **변화계 동시관측(Flare)** — 당시 일변화 보정을 어떻게 했는지 재현된다.
     ⚠️ **연도·측점마다 있고 없다** — 2015 남지는 Flare 가 전부 0 이다.

양식이 두 갈래다. **고정 좌표가 아니라 라벨 탐색**으로 읽는다(행이 밀린다).

| | 시트 | 세션 | 특징 |
|---|---|---|---|
| **A형** 2012·2013·2014 | `종합` + `측정_1..N` | 4 또는 6 | `GPS방위각`·`Mark1방위각` 행이 따로 있다 |
| **B형** 2015·2016 | `<측점>1..6` + `<측점>결과` | 6 | 좌표가 1·2행 머리글에 있다 |

    python parse_kriss_legacy.py            # 요약
    python parse_kriss_legacy.py --csv      # docs/data/kriss_legacy_sessions.csv
"""
from __future__ import annotations

import argparse
import datetime as dt
import glob
import os
import re
import sys
from pathlib import Path

import xlrd

ROOT = Path(__file__).parent
SRC = Path(r"D:\LX_yoons\2026_research\2026_지자기 연구"
           r"\20260908_주성중 박사_2020년 사용 데이터")
OUT = ROOT / "docs" / "data" / "kriss_legacy_sessions.csv"

# 파일명 → 측점명. 괄호 안 옛 이름·점의부호는 떼고 **현재 이름**으로 맞춘다
# (`existing_network.SUPERSEDED` 와 같은 대응).
FILE_SITE = re.compile(r"(?:\[결과\]\s*)?\d*\s*년?\s*([가-힣]+)")
CODE_SITE = {           # 점의부호 → 측점명 (13년 폴더가 부호로만 돼 있다)
    "1-2-2": "강화", "1-12-02": "상주", "1-13-2": "미원", "1-14-3": "서산",
    "1-15-3": "청양", "1-17-2": "부안", "1-19-4": "순창", "1-20-3": "함양",
    "1-25-2": "거제", "1-26-2": "남해", "1-27-3": "순천",
}
RENAME = {"임계": "삼척", "경주": "영천", "양산": "언양", "곡성": "순창",
          "인제": "설악"}


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _dms(d, m, s):
    """도·분·초 → 십진도. 부호는 «도» 칸에만 붙는다(분·초는 절댓값)."""
    d, m, s = _f(d), _f(m), _f(s)
    if d is None:
        return None
    sign = -1 if d < 0 else 1
    return sign * (abs(d) + (m or 0) / 60 + (s or 0) / 3600)


def _clock(v):
    """엑셀 시각 → 'HH:MM:SS'. 문자열 `1;39:43` 같은 오타도 받는다."""
    if isinstance(v, str):
        t = re.findall(r"\d+", v)
        if len(t) >= 2:
            h, m = int(t[0]), int(t[1])
            s = int(t[2]) if len(t) > 2 else 0
            return f"{h:02d}:{m:02d}:{s:02d}"
        return None
    f = _f(v)
    if f is None or not (0 <= f < 1):
        return None
    sec = round(f * 86400)
    return f"{sec//3600:02d}:{(sec//60)%60:02d}:{sec%60:02d}"


def _date(v, datemode):
    if isinstance(v, str):
        t = re.findall(r"\d+", v)
        if len(t) >= 3:
            return f"{int(t[0]):04d}-{int(t[1]):02d}-{int(t[2]):02d}"
        return None
    f = _f(v)
    if f is None or f < 1000:
        return None
    y, m, d = xlrd.xldate_as_tuple(f, datemode)[:3]
    return f"{y:04d}-{m:02d}-{d:02d}"


class Grid:
    """시트를 (행,열) 격자로 들고 라벨로 찾는다."""

    def __init__(self, sh):
        self.sh = sh
        self.n, self.m = sh.nrows, sh.ncols

    def v(self, r, c):
        if 0 <= r < self.n and 0 <= c < self.m:
            x = self.sh.cell_value(r, c)
            return x if x != "" else None
        return None

    def find(self, text, col=None, exact=False):
        """라벨이 있는 (행,열). `col` 을 주면 그 열만 본다."""
        t = str(text).strip()
        for r in range(self.n):
            for c in ([col] if col is not None else range(self.m)):
                x = self.v(r, c)
                if x is None:
                    continue
                s = str(x).strip()
                if (s == t) if exact else (t in s):
                    return r, c
        return None

    def row_after(self, text, col=None):
        hit = self.find(text, col)
        return hit[0] if hit else None


def parse_session(sh, datemode, style):
    """세션 시트 한 장 → dict. 양식 A·B 공통 처리."""
    g = Grid(sh)
    out = {}

    # ── 날짜·관측자·코드 ────────────────────────────────────────────
    hit = g.find("Date")
    if hit:
        r, c = hit
        out["날짜"] = _date(g.v(r, c + 1), datemode)
    hit = g.find("Observer")
    if hit:
        r, c = hit
        out["관측자"] = str(g.v(r, c + 1) or "").strip()
    hit = g.find("Code")
    if hit:
        r, c = hit
        out["코드"] = str(g.v(r, c + 1) or "").strip()

    # ── 좌표 ────────────────────────────────────────────────────────
    for key, label in (("위도", "위도"), ("경도", "경도")):
        hit = g.find(label)
        if hit:
            r, c = hit
            out[key] = _dms(g.v(r, c + 1), g.v(r, c + 2), g.v(r, c + 3))

    # ── 방위 기준 ───────────────────────────────────────────────────
    # `FM True` 는 두 양식 모두에 있다 — **재방문 사이 방위각 비교의 핵심**
    for key, label in (("FM_True", "FM True"), ("TN_Circle", "TN Circle"),
                       ("GPS방위각", "GPS방위각"), ("Mark1방위각", "Mark1방위각")):
        hit = g.find(label)
        if hit:
            r, c = hit
            # A형 GPS·Mark1 은 라벨 오른쪽 3칸이 도분초
            out[key] = _dms(g.v(r, c + 1), g.v(r, c + 2), g.v(r, c + 3))

    # ── 편각 세션 (ED/WD/EU/WU) ────────────────────────────────────
    times_d, times_i = [], []
    for lab in ("ED", "WD", "EU", "WU"):
        hit = g.find(lab, exact=True)
        if hit:
            r, _ = hit
            t = _clock(g.v(r, 1))
            if t:
                times_d.append(t)
    for lab in ("SU", "ND", "SD", "NU"):
        hit = g.find(lab, exact=True)
        if hit:
            r, _ = hit
            t = _clock(g.v(r, 1))
            if t:
                times_i.append(t)
    if times_d:
        out["편각시작"], out["편각종료"] = min(times_d), max(times_d)
    if times_i:
        out["복각시작"], out["복각종료"] = min(times_i), max(times_i)

    # ── 성과 (D·I·F) ───────────────────────────────────────────────
    hit = g.find("Declination")
    if hit:
        r, c = hit
        # A형은 도분초(C:E), B형은 F열 십진도
        d = _dms(g.v(r, 2), g.v(r, 3), g.v(r, 4))
        out["편각"] = d if d is not None else _f(g.v(r, 5))
    hit = g.find("Inclination")
    if hit is None:
        hit = g.find("Mean (Inc.)")
    if hit:
        r, c = hit
        i = _f(g.v(r, 5))
        if i is None:
            i = _dms(g.v(r, 2), g.v(r, 3), g.v(r, 4))
        out["복각"] = i
    hit = g.find("Corrected")
    if hit:
        r, c = hit
        out["총자력"] = _f(g.v(r, c + 1))

    # ── 변화계(Flare) 동시관측 — 있는지 없는지가 중요하다 ───────────
    fl = []
    for lab in ("Flare H", "Flare Z", "Flare F", "Flare D"):
        hit = g.find(lab)
        if hit:
            r, c = hit
            for rr in range(r + 1, min(r + 6, g.n)):
                x = _f(g.v(rr, c))
                if x is not None:
                    fl.append(abs(x))
    out["변화계"] = "있음" if any(x > 0 for x in fl) else "없음"

    # ── 시준오차 ───────────────────────────────────────────────────
    for key, label in (("Dec_delta", "Dec. δ"), ("Dec_eps", "Dec. ε"),
                       ("Inc_eps", "Inc. ε")):
        hit = g.find(label)
        if hit:
            r, c = hit
            out[key] = _f(g.v(r, c + 1))
    return out


def site_from(path, code):
    """측점명은 **파일명에서** 읽는다.

    ⚠️ **시트의 `Code` 칸을 믿으면 안 된다** — 양식을 복사해 쓴 흔적이 남아
    있다. `12년 봉평.xls`·`12년 제천.xls` 의 Code 가 둘 다 서산의 부호
    `1-14-3` 이라, 코드를 우선하면 **두 파일이 통째로 서산이 된다.**
    파일명 쪽은 좌표와도 맞는다(제천 37.168/128.016 = 백마저수지).
    """
    base = os.path.basename(path).replace(".xls", "")
    base = re.sub(r"\(.*?\)", "", base)                  # 「설악(인제)」 → 「설악」
    runs = re.findall(r"[가-힣]+", base)
    runs = [x for x in runs if x not in ("년", "결과", "직영출장결과")]
    if runs:
        nm = max(runs, key=len)
    elif code in CODE_SITE:                              # 파일명에 이름이 없을 때만
        nm = CODE_SITE[code]
    else:
        nm = base
    return RENAME.get(nm, nm)


def parse_file(path):
    wb = xlrd.open_workbook(path)
    names = wb.sheet_names()
    style = "A" if any(n.startswith("측정") for n in names) else "B"
    sheets = ([n for n in names if n.startswith("측정")] if style == "A"
              else [n for n in names
                    if re.search(r"[\(]?\d+[\)]?$", n) and "결과" not in n
                    and "그래프" not in n and not n.startswith("Sheet")])
    rows, dropped = [], []
    for k, sn in enumerate(sheets, 1):
        rec = parse_session(wb.sheet_by_name(sn), wb.datemode, style)
        if not rec.get("편각") and not rec.get("복각"):
            continue
        rec["파일"] = os.path.basename(path)
        rec["양식"] = style
        rec["회차"] = k
        rec["측점"] = site_from(path, rec.get("코드"))
        rows.append(rec)

    # ── 버리지 않고 **성분마다 표시**한다 ─────────────────────────
    # ⚠️ 세션을 통째로 떨어뜨리면 멀쩡한 성분까지 잃는다. 실제로 그랬다 —
    #    15년 제천은 D·I 가 정상인데 F 만 66,410 nT 라 전부 날아갔었다.
    #    「값이 이상해도 버리지 않고 그대로 보여 준다」가 이 저장소 규약이다.
    #
    # 한반도 범위는 `lmm_build.parse_results()` 와 같은 규약을 쓴다.
    for r in rows:
        d, i, f = r.get("편각"), r.get("복각"), r.get("총자력")
        r["D_판정"] = ("결측" if d is None else
                       "정상" if abs(d) <= 20 else "범위밖")
        r["I_판정"] = ("결측" if i is None else
                       "정상" if 30 <= i <= 70 else "범위밖")
        r["F_판정"] = ("결측" if f is None else
                       "정상" if 30000 <= f <= 60000 else "범위밖")

    # ⚠️ **양식 복사 잔재** — 13년 파일 6개에 값이 **완전히 같은**
    #    2010-10-22 시트가 붙어 있다(D −8.6269 · I **90.0°** · F 50730.5858).
    #    복각 90° 는 자극이라 한반도에서 나올 수 없다. 14년 거제도 6세션이
    #    전부 I≈90° 라 **그 파일에는 쓸 자료가 없다.**
    days = [r["날짜"] for r in rows if r.get("날짜")]
    base_day = max(set(days), key=days.count) if days else None
    for r in rows:
        why = []
        if r["I_판정"] == "범위밖" and r.get("복각") and abs(r["복각"] - 90) < 0.1:
            why.append("복각 90° (빈 양식)")
        if base_day and r.get("날짜"):
            gap = abs((dt.date.fromisoformat(r["날짜"])
                       - dt.date.fromisoformat(base_day)).days)
            if gap > 3:
                why.append(f"본 날짜에서 {gap}일")
        r["잔재"] = " · ".join(why)
        if why:
            dropped.append((r["회차"], r.get("날짜"), r["잔재"]))
    for k, d, w in dropped:
        print(f"  · 잔재 표시 {os.path.basename(path)} 회차{k} {d} — {w}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", action="store_true")
    a = ap.parse_args()

    files = sorted(glob.glob(str(SRC / "*.xls"))
                   + glob.glob(str(SRC / "13년 직영출장결과" / "*.xls")))
    if not files:
        print(f"[오류] 원자료를 못 찾음: {SRC}")
        return 1

    rows = []
    for p in files:
        try:
            rows += parse_file(p)
        except Exception as exc:                      # 한 파일이 깨져도 계속
            print(f"  ⚠ {os.path.basename(p)} — {exc}")
    if not rows:
        print("[오류] 세션을 하나도 못 읽음")
        return 1

    import pandas as pd
    df = pd.DataFrame(rows)
    df["연도"] = pd.to_datetime(df["날짜"], errors="coerce").dt.year
    cols = ["측점", "연도", "날짜", "회차", "관측자", "코드", "위도", "경도",
            "편각", "복각", "총자력", "D_판정", "I_판정", "F_판정", "잔재",
            "편각시작", "편각종료", "복각시작", "복각종료", "변화계",
            "FM_True", "TN_Circle", "GPS방위각", "Mark1방위각",
            "Dec_delta", "Dec_eps", "Inc_eps", "양식", "파일"]
    df = df[[c for c in cols if c in df.columns]]
    ok = df["잔재"].eq("")

    print(f"\n■ 파일 {len(files)} · 세션 {len(df)}"
          f" (잔재 {(~ok).sum()} 빼면 {ok.sum()}) · 측점 {df.loc[ok,'측점'].nunique()}")
    print("\n연도별 — 잔재를 뺀 세션 기준")
    print(df[ok].groupby("연도").agg(
        측점=("측점", "nunique"), 세션=("측점", "size"),
        D정상=("D_판정", lambda s: (s == "정상").sum()),
        I정상=("I_판정", lambda s: (s == "정상").sum()),
        F정상=("F_판정", lambda s: (s == "정상").sum()),
        시각=("편각시작", lambda s: s.notna().sum()),
        변화계=("변화계", lambda s: (s == "있음").sum())).to_string())

    bad = df[ok & ((df.D_판정 != "정상") | (df.F_판정 != "정상")
                   | (df.I_판정 != "정상"))]
    if len(bad):
        print("\n⚠ 잔재는 아니나 값이 범위 밖 — 버리지 않고 표시만 한다")
        for _, r in bad.iterrows():
            print(f"   {r['측점']:<4} {r['날짜']} 회차{r['회차']}  "
                  f"D={r['편각']} ({r['D_판정']}) · "
                  f"I={r['복각']} ({r['I_판정']}) · "
                  f"F={r['총자력']} ({r['F_판정']})")
    print("\n측점별 (방문연도)")
    g = (df.groupby("측점")["연도"].apply(lambda s: "·".join(map(str, sorted(set(s)))))
         .sort_index())
    for k, v in g.items():
        print(f"  {k:<5} {v}")

    if a.csv:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(OUT, index=False, encoding="utf-8-sig")
        print(f"\n[저장] {OUT}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
