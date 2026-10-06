# -*- coding: utf-8 -*-
"""야장 세션 ↔ 성과표 대조 — 성과표 각 행에 실제 관측일시를 붙일 수 있는지 본다.

성과표(지자기측량 성과정리)에는 연도만 있어 lmm_build 가 7월 1일로 대체하고 있다
([lmm_build.py] df["date"] = dt.datetime(y, 7, 1)). 야장에서 복원한 일시로
이 대체값을 실제 관측일시로 바꿀 수 있는지 (측점, 연도) 로 맞춰 확인한다.
"""
import sys, warnings
warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd
from pathlib import Path
from lmm_build import load_survey_points

ROOT = Path(__file__).parent
fb = pd.read_csv(ROOT / "docs/data/fieldbook_sessions.csv", encoding="utf-8-sig")
fb["날짜"] = pd.to_datetime(fb["날짜"])
fb["연도"] = fb["날짜"].dt.year

sp = load_survey_points()
sp["year"] = sp["year"].astype(int)

print(f"성과표 {len(sp)}행 · {sp['name'].nunique()}측점 "
      f"({sp['year'].min()}~{sp['year'].max()})")
print(f"야장   {len(fb)}세션 · {fb['측점'].nunique()}측점 "
      f"({fb['연도'].min()}~{fb['연도'].max()})\n")

# 야장 측점·연도별 대표 일시(첫 세션)와 세션 수
g = (fb.sort_values(["측점", "날짜", "시작"])
       .groupby(["측점", "연도"])
       .agg(일자수=("날짜", "nunique"), 세션=("날짜", "size"),
            첫날=("날짜", "first"), 첫시각=("시작", "first"),
            시각있음=("시작", lambda s: s.notna().sum()))
       .reset_index())

m = sp.merge(g, left_on=["name", "year"], right_on=["측점", "연도"], how="left")
hit = m["첫날"].notna()
print(f"■ 성과표 {len(m)}행 중 야장 매칭 {hit.sum()}행 · 미매칭 {(~hit).sum()}행\n")

print(f"{'측점':8} {'연도':5} {'성과표대체일':12} {'야장 첫 관측':12} {'시각':9} {'세션':>4}  판정")
print("-" * 78)
for _, r in m.sort_values(["name", "year"]).iterrows():
    stub = r["date"].strftime("%Y-07-01")
    if pd.notna(r["첫날"]):
        real = pd.Timestamp(r["첫날"]).strftime("%Y-%m-%d")
        t = str(r["첫시각"])[:8] if pd.notna(r["첫시각"]) else "—"
        gap = abs((pd.Timestamp(r["첫날"]) - r["date"]).days)
        v = f"복원 가능 (7/1 대체와 {gap}일 차)"
        print(f"{r['name'][:7]:8} {r['year']:5} {stub:12} {real:12} {t:9} "
              f"{int(r['세션']):4}  {v}")
    else:
        print(f"{r['name'][:7]:8} {r['year']:5} {stub:12} {'—':12} {'—':9} "
              f"{'—':>4}  야장 없음")

print("\n■ 야장에는 있으나 성과표에 없는 (측점, 연도)")
back = g.merge(sp[["name", "year"]], left_on=["측점", "연도"],
               right_on=["name", "year"], how="left")
for _, r in back[back["name"].isna()].sort_values(["측점", "연도"]).iterrows():
    print(f"   {r['측점']:8} {int(r['연도'])}  세션 {int(r['세션'])}건 "
          f"({pd.Timestamp(r['첫날']).strftime('%Y-%m-%d')}~)")
