# -*- coding: utf-8 -*-
"""
약도(略圖) 생성 — 기준점 + 방위표지 1·2 를 위성지도 위에 표시
=============================================================

조사 카드 ③ 방위표지 조사의 기준점·방위표지 좌표로 위성영상 배경 약도를 만든다.
지자기점(중심)에서 방위표지 1·2 로 향하는 진북기준 방위각·거리를 라벨로 표시하고
북 화살표를 넣는다. 결과 PNG 는 merge 단계에서 카드 ⑥ 약도 칸에 삽입된다.

    python make_sketches.py            # 전체(좌표 완비 카드)
    python make_sketches.py --only 12  # DS-012 만 (테스트)

출력: docs/output/_sketches/{DS-0xx}.png
"""
import argparse
import math
import re
import socket
import sys
from pathlib import Path

socket.setdefaulttimeout(20)   # 타일 다운로드 무한 대기 방지

import contextily as ctx
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

from aggregate_survey_xlsx import parse_workbook, survey_files

ROOT = Path(__file__).parent
OUT = ROOT / "docs" / "output" / "_sketches"

# 한글 폰트
for _cand in ("Malgun Gothic", "맑은 고딕", "NanumGothic"):
    if any(_cand.lower() in f.name.lower() for f in fm.fontManager.ttflist):
        plt.rcParams["font.family"] = _cand
        break
plt.rcParams["axes.unicode_minus"] = False

R = 6378137.0


def to3857(lon, lat):
    x = R * math.radians(lon)
    y = R * math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
    return x, y


def fmt_az(s):
    """방위각 문자열 정규화 → 'ddd°mm′ss″'. 실패 시 원문."""
    nums = re.findall(r"\d+\.?\d*", str(s))
    if len(nums) >= 3:
        return f"{int(float(nums[0]))}°{int(float(nums[1]))}′{float(nums[2]):.0f}″"
    return str(s)


def _valid(ll):
    return ll and 124 <= ll[0] <= 132 and 33 <= ll[1] <= 39   # (lon, lat) 한반도


def make_sketch(d, out_path):
    base = d.get("기준점ll")
    m1 = d.get("표지1ll")
    m2 = d.get("표지2ll")
    if not (base and m1 and m2):
        return False
    if not all(_valid(p) for p in (base, m1, m2)):
        print(f"    ! 좌표 범위 이상 건너뜀 {d['관리번호']} {d['후보지명']}: "
              f"기준점={base} 표지1={m1} 표지2={m2}", file=sys.stderr)
        return False
    det = d.get("방위표지상세", {})
    pts = {
        "기준점": (base, "#E8531F", "*", 320),
        "방위표지 1": (m1, "#1D4ED8", "^", 150),
        "방위표지 2": (m2, "#0E8A6B", "^", 150),
    }
    xs, ys = [], []
    xy = {}
    for k, (ll, *_rest) in pts.items():
        x, y = to3857(*ll)
        xy[k] = (x, y)
        xs.append(x)
        ys.append(y)

    cx, cy = xy["기준점"]
    span = max(max(xs) - min(xs), max(ys) - min(ys), 40)
    pad = span * 0.55 + 15
    half = span / 2 + pad
    xlim = (cx - half, cx + half)
    ylim = (cy - half, cy + half)

    fig, ax = plt.subplots(figsize=(6.2, 5.0), dpi=140)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")

    # 라인: 기준점 → 각 방위표지
    for k, col in (("방위표지 1", "#1D4ED8"), ("방위표지 2", "#0E8A6B")):
        x, y = xy[k]
        ax.plot([cx, x], [cy, y], color=col, lw=2.0, zorder=5,
                solid_capstyle="round")
        tag = "표지1" if k.endswith("1") else "표지2"
        az = fmt_az(det.get(tag, {}).get("방위각", ""))
        dist = det.get(tag, {}).get("거리", "")
        mx, my = (cx + x) / 2, (cy + y) / 2
        ax.annotate(f"{az}\n{dist} m", (mx, my), color="#111", fontsize=9,
                    ha="center", va="center", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=col, alpha=0.9))

    # 점 + 라벨
    for k, (ll, col, mk, sz) in pts.items():
        x, y = xy[k]
        ax.scatter([x], [y], c=col, marker=mk, s=sz, ec="white", lw=1.4, zorder=8)
        ax.annotate(k, (x, y), color=col, fontsize=9.5, fontweight="bold",
                    xytext=(0, 11), textcoords="offset points", ha="center", zorder=9,
                    bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=col, alpha=0.85))

    # 위성 배경
    for z in (18, 17, 16):
        try:
            ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery,
                            crs="EPSG:3857", attribution=False, zoom=z)
            break
        except Exception as e:   # noqa: BLE001
            print(f"    ! basemap z{z} 실패 {d['관리번호']}: {e}", file=sys.stderr)

    # 북 화살표
    ax.annotate("N", xy=(0.94, 0.93), xytext=(0.94, 0.80),
                xycoords="axes fraction", textcoords="axes fraction",
                ha="center", va="center", fontsize=12, fontweight="bold", color="white",
                arrowprops=dict(arrowstyle="-|>", color="white", lw=2.2))

    # 축척 막대 (하단 좌측, 20 m 또는 span 에 맞춰)
    bar_m = 20 if span < 120 else 50
    x0 = xlim[0] + (xlim[1] - xlim[0]) * 0.06
    y0 = ylim[0] + (ylim[1] - ylim[0]) * 0.06
    ax.plot([x0, x0 + bar_m], [y0, y0], color="white", lw=3, zorder=10)
    ax.annotate(f"{bar_m} m", (x0 + bar_m / 2, y0), color="white", fontsize=8.5,
                ha="center", va="bottom", zorder=10, fontweight="bold")

    ax.set_title(f"약도  |  {d['관리번호']}  {d['후보지명']}", fontsize=11, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout(pad=0.6)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="DS 번호(예: 12) 하나만")
    ap.add_argument("--force", action="store_true", help="기존 약도도 재생성")
    a = ap.parse_args()
    recs = []
    for f in survey_files():
        recs += parse_workbook(f)
    if a.only:
        want = f"DS-{int(a.only):03d}"
        recs = [d for d in recs if d["관리번호"] == want]

    ok = skip = exist = 0
    for d in recs:
        out = OUT / f"{d['관리번호']}.png"
        if out.exists() and not a.force and not a.only:
            exist += 1
            continue
        if make_sketch(d, out):
            ok += 1
            print(f"  약도 저장 {d['관리번호']} {d['후보지명']}", flush=True)
        else:
            skip += 1
    print(f"\n신규 {ok}건 / 기존유지 {exist}건 / 좌표미비 {skip}건 → {OUT}")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
