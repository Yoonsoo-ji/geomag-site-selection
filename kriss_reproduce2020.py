# -*- coding: utf-8 -*-
"""
2020 보고서 표4-12 재현 — 원 관측자료에서 공표 성과가 나오는가
================================================================

2020 연구보고서 <표4-12>「1등자기점 절대측정 결과」는 2012~2016 절대관측을
측점·연도마다 한 값으로 공표했다. `parse_kriss_legacy.py` 가 읽은 **세션별
원시 판독값**에서 그 값이 나오는지 확인한다. 나오면 **이 파일들이 그
보고서의 원자료**임이 확정되고, 보고서가 어느 세션을 뺐는지도 드러난다.

두 단계다.

  ① **고정 규칙으로 평균** — 전 세션 / 부호 교정 / 범위 밖 제외 / 블런더 제외

  ② **세션 조합 되짚기** — 고정 규칙으로 안 맞는 쌍은 **어느 세션을 평균했을
     때 보고서 값이 나오는지** 부분집합을 전수로 훑어 찾는다. 세션이 6개
     이하라 조합이 63개뿐이다.

②가 요점이다. **보고서는 하나의 규칙을 쓰지 않았다** — 쌍마다 뺀 세션이
다르다. 그러니 「규칙 하나로 몇 개 맞는가」로는 재현 여부를 판정할 수 없고,
**각 쌍이 어느 세션 조합에서 나오는가**를 보여야 한다.

⚠️ **보고서 값을 손으로 적지 않는다.** 2020 보고서 PDF 의 표를 긁어 쓴다.
   ⚠️ 본문 글줄로는 못 읽는다 — 표라서 `find_tables()` 가 필요하다.

⚠️ **「조합을 찾았다」가 「그 세션이 틀렸다」는 뜻은 아니다.** 보고서가 왜
   그 세션을 뺐는지는 보고서에 적혀 있지 않다. 되짚은 조합은 **보고서의
   처리를 복원한 것**이고, 그 판단이 옳았는지는 별개의 물음이다.

    python kriss_reliability.py      # ← 먼저. 블런더 표시를 만든다
    python kriss_reproduce2020.py
"""
from __future__ import annotations

import itertools
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
SESS = ROOT / "docs" / "data" / "kriss_legacy_sessions.csv"
FLAGS = ROOT / "docs" / "data" / "kriss_session_flags.csv"
OUT = ROOT / "docs" / "data" / "kriss_reproduce2020.csv"

PDF = Path(r"D:\LX_yoons\2026_research\2026_지자기 연구\20251210_선행연구"
           r"\200125 2020 지구물리측량 연구사업 (최종).pdf")

# 보고서는 옛 이름을 쓴다 — 세션표는 현재 이름이다
RENAME = {"임계": "삼척", "경주": "영천", "양산": "언양",
          "곡성": "순창", "인제": "설악"}

TOL = 1.0        # 재현 판정 문턱 (초)
DMS = re.compile(r"(-?\d+)\s*[°]\s*(\d+)\s*['′]\s*([\d.]+)\s*[\"″]?")


def dms(txt):
    """도분초 문자열을 십진도로. 부호는 「도」 칸에만 붙는다."""
    m = DMS.search(str(txt))
    if not m:
        return None
    d, mi, s = float(m.group(1)), float(m.group(2)), float(m.group(3))
    sign = -1.0 if m.group(1).strip().startswith("-") else 1.0
    return sign * (abs(d) + mi / 60 + s / 3600)


def report_table():
    """표4-12 에서 측점·연도별 공표 편각·복각을 긁는다."""
    try:
        import fitz
    except ImportError:
        print("⚠ pymupdf 가 없다 — 보고서 표를 읽을 수 없다")
        return None
    if not PDF.exists():
        print(f"⚠ 보고서 PDF 가 없다: {PDF}")
        return None

    doc = fitz.open(PDF)
    rows = []
    for pno in range(len(doc)):
        txt = doc[pno].get_text()
        if "점의 번호" not in txt and "점의번호" not in txt:
            continue
        for t in doc[pno].find_tables().tables:
            tab = t.extract()
            if not tab or "편각" not in " ".join(str(c) for c in tab[0]):
                continue
            year = None   # 연도 칸은 그 해 첫 줄에만 있다 — 이어서 쓴다
            for r in tab[1:]:
                cell = [("" if c is None else str(c).strip()) for c in r]
                if len(cell) < 5:
                    continue
                m = re.search(r"(20[0-2]\d)", cell[0])
                if m:
                    year = int(m.group(1))
                if year is None or not (2012 <= year <= 2016):
                    continue
                nm, d, i = cell[2], dms(cell[3]), dms(cell[4])
                if not nm or d is None or i is None:
                    continue
                rows.append(dict(측점=RENAME.get(nm, nm), 연도=year,
                                 보고서_편각=d, 보고서_복각=i,
                                 보고서_비고=cell[-1], 쪽=pno + 1))
    doc.close()
    if not rows:
        print("⚠ 보고서에서 표4-12 를 찾지 못했다 — 쪽 추출 방식을 확인할 것")
        return None
    R = pd.DataFrame(rows).drop_duplicates(subset=["측점", "연도"])
    print(f"■ 보고서에서 읽은 성과 {len(R)} 쌍 (쪽 {sorted(set(R['쪽']))})")
    return R


def blunder_map(flags):
    bl = {}
    if flags is not None and "D_블런더" in flags.columns:
        for _, r in flags[flags["D_블런더"].astype(bool)].iterrows():
            bl.setdefault((r["측점"], int(r["연도"])), set()).add(int(r["회차"]))
    return bl


def fixed_rules(df, bl):
    """고정 규칙 넷으로 방문 대표값을 낸다."""
    rows = []
    for (s, y), g in df.groupby(["측점", "연도"]):
        d = pd.to_numeric(g["편각"], errors="coerce")
        if not d.notna().any():
            continue
        d2 = -d.abs()                                    # 한반도 편각은 음수
        inrange = g["D_판정"] == "정상"
        keep = ~g["회차"].astype(int).isin(bl.get((s, int(y)), set()))
        rows.append(dict(
            측점=s, 연도=int(y), 세션=len(g),
            그대로=float(d.mean()),
            부호교정=float(d2.mean()),
            범위밖제외=(float(d2[inrange].mean()) if inrange.any() else np.nan),
            블런더제외=(float(d2[keep].mean()) if keep.any() else np.nan)))
    return pd.DataFrame(rows)


def find_subset(d2, rounds, target):
    """
    보고서 값이 나오는 세션 조합을 찾는다. 작게 빼는 쪽을 먼저 보므로
    가장 적게 뺀 조합이 돌아온다.
    """
    n = len(d2)
    best = None
    for k in range(n, 0, -1):                 # 많이 쓰는 조합부터
        for idx in itertools.combinations(range(n), k):
            diff = abs(float(np.mean([d2[j] for j in idx])) - target) * 3600
            if diff <= TOL:
                used = sorted(rounds[j] for j in idx)
                drop = sorted(set(rounds) - set(used))
                if best is None or len(drop) < len(best[1]):
                    best = (used, drop, diff)
        if best is not None:
            return best
    return None


def main():
    df = pd.read_csv(SESS, encoding="utf-8-sig")
    df = df[df["잔재"].isna() | (df["잔재"].astype(str).str.strip() == "")]
    flags = pd.read_csv(FLAGS, encoding="utf-8-sig") if FLAGS.exists() else None
    if flags is None:
        print("⚠ kriss_session_flags.csv 가 없다 — kriss_reliability.py 를"
              " 먼저 돌릴 것 (블런더 규칙 열이 빠진다)")
    bl = blunder_map(flags)

    R = report_table()
    V = fixed_rules(df, bl)
    print(f"■ 원자료에서 낸 방문 대표값 {len(V)} 쌍")
    if R is None:
        V.to_csv(OUT, index=False, encoding="utf-8-sig")
        print(f"\n보고서 값을 못 읽었으므로 대조를 멈춘다.\n[저장] {OUT}")
        return 1

    T = V.merge(R, on=["측점", "연도"], how="inner")
    print(f"■ 대조 가능한 쌍 {len(T)}\n")

    print(f"■ ① 고정 규칙 — {TOL:.0f}″ 이내로 재현되는 쌍")
    for col in ("그대로", "부호교정", "범위밖제외", "블런더제외"):
        e = (T[col] - T["보고서_편각"]).abs() * 3600
        print(f"  {col:<12}{(e <= TOL).sum():>3} / {e.notna().sum()}"
              f"   중앙 편차 {e.median():>8.1f}″")
    print("  → 규칙 하나로는 모두 재현되지 않는다. 쌍마다 뺀 세션이 다르다.")

    print(f"\n■ ② 세션 조합 되짚기 — 보고서 값이 나오는 조합")
    print(f"{'측점':<6}{'연도':<6}{'세션':>4}{'쓴 회차':>14}{'뺀 회차':>10}"
          f"{'편차':>9}  블런더 표시와")
    rec, miss = [], []
    for _, r in T.iterrows():
        g = df[(df["측점"] == r["측점"]) & (df["연도"] == r["연도"])]
        d2 = (-pd.to_numeric(g["편각"], errors="coerce").abs()).tolist()
        rounds = [int(x) for x in g["회차"]]
        res = find_subset(d2, rounds, r["보고서_편각"])
        flagged = bl.get((r["측점"], int(r["연도"])), set())
        if res is None:
            miss.append((r, flagged))
            continue
        used, drop, diff = res
        agree = ("일치" if set(drop) == flagged else
                 "블런더 없음" if not drop and not flagged else
                 f"표시 {sorted(flagged) or '없음'}")
        print(f"{r['측점']:<6}{int(r['연도']):<6}{len(rounds):>4}"
              f"{','.join(map(str, used)):>14}"
              f"{(','.join(map(str, drop)) or '-'):>10}{diff:>8.2f}″  {agree}")
        rec.append(dict(측점=r["측점"], 연도=int(r["연도"]),
                        쓴회차=",".join(map(str, used)),
                        뺀회차=",".join(map(str, drop)),
                        편차초=round(diff, 3),
                        블런더표시=",".join(map(str, sorted(flagged)))))
    print(f"\n  재현 {len(rec)} / {len(T)} 쌍")

    if miss:
        print("\n■ 어떤 조합으로도 재현되지 않는 쌍")
        print(f"{'측점':<6}{'연도':<6}{'세션':>4}{'원자료 범위':>24}"
              f"{'보고서':>12}  비고")
        for r, flagged in miss:
            g = df[(df["측점"] == r["측점"]) & (df["연도"] == r["연도"])]
            d2 = -pd.to_numeric(g["편각"], errors="coerce").abs()
            print(f"{r['측점']:<6}{int(r['연도']):<6}{len(g):>4}"
                  f"{f'{d2.min():.5f} ~ {d2.max():.5f}':>24}"
                  f"{r['보고서_편각']:>12.5f}  {r['보고서_비고'] or ''}")
            rec.append(dict(측점=r["측점"], 연도=int(r["연도"]),
                            쓴회차="", 뺀회차="", 편차초=None,
                            블런더표시=",".join(map(str, sorted(flagged)))))

    pd.DataFrame(rec).merge(T, on=["측점", "연도"], how="right").to_csv(
        OUT, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {OUT}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
