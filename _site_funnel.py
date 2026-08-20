# -*- coding: utf-8 -*-
"""
Regional 입력 16측점의 근거 추적 — 「자료가 얼마나 있었고 무엇이 왜 빠졌나」
=============================================================================

세 갈래를 한자리에 놓는다.

  ① 1등 지자기점 30점                  — 관측망 전체
  ② 실제 관측 자료가 있는 (측점,연도)   — 성과표 + 지리원 야장
  ③ 모델에 실제로 들어간 것             — lmm_build.load_all_points()

②와 ③의 차이가 「빠진 것」이고, 그 사유를 단계별로 센다.

    python _site_funnel.py
"""
import sys
import warnings
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))


def canon(n):
    """개명·별칭을 하나로 — 언양은 양산의 다른 이름이라 별도 측점이 아니다."""
    return {"임계": "삼척", "경주": "영천", "언양": "양산"}.get(n, n)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    warnings.filterwarnings("ignore")
    import lmm_build as L
    import legacy2020_sets as LS

    # ── ① 관측망 ─────────────────────────────────────────────────────────
    net = {canon(n) for n in LS.pts30()}

    # ── ② 자료가 존재하는 (측점, 연도) ───────────────────────────────────
    #    성과표
    sv = pd.read_excel(ROOT / "data" / "지자기측량 성과정리(22_25).xlsx",
                       sheet_name="Sheet1")
    sv = sv.rename(columns={"도엽명": "name", "관측연도": "year",
                            "실.2": "D", "실.3": "I", "총자력": "F"})
    sv["name"] = sv["name"].ffill()
    sv_ok = sv.dropna(subset=["D", "I", "F"])
    sv_pairs = {(canon(r["name"]), int(r["year"])) for _, r in sv_ok.iterrows()}

    #    야장 (지리원 원본 트리 전체)
    fb = pd.read_csv(ROOT / "docs" / "data" / "fieldbook_sessions.csv",
                     encoding="utf-8-sig")
    fb["year"] = fb["날짜"].astype(str).str[:4].astype(int)
    fb_pairs = {(canon(r["측점"]), int(r["year"])) for _, r in fb.iterrows()}

    have = sv_pairs | fb_pairs

    # ── ③ 모델 투입분 ────────────────────────────────────────────────────
    used = L.load_all_points(include_2019=L.INCLUDE_2019)
    used_pairs = {(canon(n), int(y))
                  for n, y in zip(used["name"], used["year"])}

    P = lambda s: " · ".join(sorted(s))            # noqa: E731

    print("═" * 76)
    print("① 관측망 — 1등 지자기점")
    print("═" * 76)
    print(f"  {len(net)}점")

    print("\n" + "═" * 76)
    print("② 자료가 존재하는 (측점, 연도)")
    print("═" * 76)
    print(f"  성과표(’22~25)      {len(sv_pairs):>3}쌍 / {len({p[0] for p in sv_pairs})}측점")
    print(f"  지리원 야장(’19~25) {len(fb_pairs):>3}쌍 / {len({p[0] for p in fb_pairs})}측점"
          f" · 세션 {len(fb)}건")
    print(f"  합집합              {len(have):>3}쌍 / {len({p[0] for p in have})}측점")
    print(f"\n  연도별 야장 세션")
    for y, g in fb.groupby("year"):
        print(f"    {y}  {len(g):>3}세션 / {g['측점'].nunique()}측점  "
              f"{P({canon(x) for x in g['측점']})}")

    print("\n" + "═" * 76)
    print("③ 모델에 실제로 들어간 것")
    print("═" * 76)
    print(f"  {len(used)}행 / {len({p[0] for p in used_pairs})}측점")

    print("\n" + "═" * 76)
    print("④ 빠진 것 — ② 에는 있으나 ③ 에 없는 (측점, 연도)")
    print("═" * 76)
    miss = sorted(have - used_pairs)
    print(f"  총 {len(miss)}쌍\n")

    # ── 사유 분류 — 코드가 실제로 어디서 잘라냈는지에 맞춘다 ─────────────
    #   ㄱ 성과표에 있으나 신뢰도 감사에서 배제 (REGIONAL_FILTER="strict")
    #   ㄴ 2019 야장에는 있으나 투입 규칙에서 배제
    #   ㄷ 야장만 있고 성과가 아직 없음 — 투입 경로 자체가 열려 있지 않음
    strict = {("부안", 2023), ("포천", 2023)}
    fb2019 = {("미원", 2019), ("성산", 2019)}
    rest = set(miss) - strict - fb2019

    print(f"  ㄱ. 성과표에 있으나 신뢰도 감사로 배제 : {len(strict & set(miss))}쌍")
    print("       부안 2023 · 포천 2023 — 재방문 잔여 부호 반전(+62.2′→−70.6′,")
    print("       −32.2′→+22.3′)으로 그 방문의 방위 기준 오류가 특정됨")

    print(f"\n  ㄴ. 2019 야장에 있으나 투입 규칙에서 배제 : "
          f"{len(fb2019 & set(miss))}쌍")
    print("       미원 2019 — 측점 이력 등급 B (2012~19 편각 신뢰도 문제)")
    print("       성산 2019 — 야장에 관측시각 미기입 (212세션 중 유일)")

    print(f"\n  ㄷ. 야장만 있고 성과가 없음 — 투입 경로 미개통 : {len(rest)}쌍")
    for y in sorted({q[1] for q in rest}):
        print(f"       {y}  {P({q[0] for q in rest if q[1] == y})}")

    print("\n" + "═" * 76)
    print("⑤ 측점 단위 — 자료는 있으나 모델에 한 행도 못 들어간 측점")
    print("═" * 76)
    have_s = {p[0] for p in have}
    used_s = {p[0] for p in used_pairs}
    print(f"  자료 보유 {len(have_s)}측점 → 모델 {len(used_s)}측점")
    print(f"  탈락 : {P(have_s - used_s)}")
    print(f"\n  관측망 30점 중 ’19~25 자료가 아예 없는 측점 ({len(net - have_s)}점)")
    print(f"    {P(net - have_s)}")


if __name__ == "__main__":
    main()
