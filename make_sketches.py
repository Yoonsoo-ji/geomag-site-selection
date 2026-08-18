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

from aggregate_survey_xlsx import (parse_workbook, survey_files,
                                   _geo_dist, _card_dist)

ROOT = Path(__file__).parent
OUT = ROOT / "docs" / "output" / "_sketches"

# 한글 폰트
for _cand in ("Malgun Gothic", "맑은 고딕", "NanumGothic"):
    if any(_cand.lower() in f.name.lower() for f in fm.fontManager.ttflist):
        plt.rcParams["font.family"] = _cand
        break
plt.rcParams["axes.unicode_minus"] = False

R = 6378137.0

# 카드 기재값이 좌표와 어긋난 건 — 실행 끝에 요약해 재확인 대상으로 낸다.
FLAGS = []


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


def az_deg(s):
    """카드 방위각 문자열 → 도(십진). 실패 시 None."""
    nums = re.findall(r"\d+\.?\d*", str(s))
    if not nums:
        return None
    v = float(nums[0])
    if len(nums) >= 2:
        v += float(nums[1]) / 60
    if len(nums) >= 3:
        v += float(nums[2]) / 3600
    return v % 360


def az_from_ll(a, b):
    """a→b 진북기준 방위각(도). a, b = (lon, lat)."""
    if not a or not b:
        return None
    lon1, lat1 = a
    lon2, lat2 = b
    la = math.radians((lat1 + lat2) / 2)
    m_lat = 111132.92 - 559.82 * math.cos(2 * la) + 1.175 * math.cos(4 * la)
    m_lon = 111412.84 * math.cos(la) - 93.5 * math.cos(3 * la)
    return math.degrees(math.atan2((lon2 - lon1) * m_lon,
                                   (lat2 - lat1) * m_lat)) % 360


def az_label(base, mark_ll, det, tag):
    """방위각 라벨 — 카드 기재값이 원칙이나 **좌표와 어긋나면 함께 적는다**.

    약도의 삼각형은 좌표로 찍히므로 라벨만 카드값이면 그림과 모순된다.
    196쌍 중 29쌍이 5° 넘게 어긋나고 그중 16쌍은 정확히 ±180°(역방위 —
    표지→기준점 방향을 적었거나 두 좌표를 맞바꿔 적은 것으로 보인다).
    측정 방위각이 좌표(휴대GPS)보다 정확할 수 있어 카드값을 지우지 않고,
    5° 이내면 카드값만, 넘으면 좌표값을 병기한다.

    반환 (라벨문자열, 사유) — 사유는 "" / "역방위" / "불일치".
    """
    card_s = det.get(tag, {}).get("방위각", "")
    card = az_deg(card_s)
    geo = az_from_ll(base, mark_ll)
    if card is None:
        return (f"{geo:.1f}° (좌표)" if geo is not None else "-"), ""
    if geo is None:
        return fmt_az(card_s), ""
    diff = abs((geo - card + 180) % 360 - 180)
    if diff <= 5:
        return fmt_az(card_s), ""
    why = "역방위" if abs(diff - 180) <= 5 else "불일치"
    return f"{fmt_az(card_s)} (카드)\n{geo:.1f}° (좌표·{why})", why


def dist_label(base, mark_ll, det, tag):
    """거리 라벨 — **좌표 계산 측지거리**를 쓴다.

    카드 기재 거리는 오기입 사례가 있어(경남 DS-004·045·084 은 좌표상 20~70m 를
    100~117m 로 적음) 삼각형 위치와 라벨이 어긋난다. `mark_max_dist()` 의 등급
    판정과 같은 기준(`_geo_dist` 우선, 좌표 없으면 카드값 폴백)을 쓰되, 카드값이
    유의하게 다르면 괄호로 함께 적어 기재 이력을 남긴다.

    반환 (라벨문자열, 카드값과 불일치 여부).
    """
    geo = _geo_dist(base, mark_ll)
    card = _card_dist(det, tag)
    if geo is None:                       # 좌표 없음 → 카드값 폴백
        return (f"{card:.0f} m (카드)" if card is not None else "- m"), False
    lab = f"{geo:.0f} m"
    # 불일치 기준: 5 m 넘게 & 10% 넘게 차이 (반올림·소수 표기 차이는 무시)
    if card is not None and abs(card - geo) > 5 and abs(card - geo) > geo * 0.10:
        return f"{lab}  (카드 {card:.0f} m)", True
    return lab, False


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
        mll = m1 if tag == "표지1" else m2
        az, az_why = az_label(base, mll, det, tag)
        dl, d_bad = dist_label(base, mll, det, tag)
        bad = bool(az_why) or d_bad
        if bad:
            FLAGS.append((d["관리번호"], d["후보지명"], tag,
                          az_why, dl if d_bad else ""))
        mx, my = (cx + x) / 2, (cy + y) / 2
        ax.annotate(f"{az}\n{dl}", (mx, my), color="#111",
                    fontsize=8 if bad else 9,
                    ha="center", va="center", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white",
                              ec="#CC3333" if bad else col,
                              lw=1.6 if bad else 1.0, alpha=0.9))

    # 점 + 라벨
    for k, (ll, col, mk, sz) in pts.items():
        x, y = xy[k]
        ax.scatter([x], [y], c=col, marker=mk, s=sz, ec="white", lw=1.4, zorder=8)
        ax.annotate(k, (x, y), color=col, fontsize=9.5, fontweight="bold",
                    xytext=(0, 11), textcoords="offset points", ha="center", zorder=9,
                    bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=col, alpha=0.85))

    # 위성 배경
    basemap_ok = False
    for z in (18, 17, 16):
        try:
            ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery,
                            crs="EPSG:3857", attribution=False, zoom=z)
            basemap_ok = True
            break
        except Exception as e:   # noqa: BLE001
            print(f"    ! basemap z{z} 실패 {d['관리번호']}: {e}", file=sys.stderr)
    if not basemap_ok:
        # 위성배경 없이 저장하면 '빈 약도'가 남는다 → 저장하지 않고 다음 실행에서 재시도
        plt.close(fig)
        print(f"    ! 위성배경 실패로 저장 생략(재시도 대상): {d['관리번호']}", file=sys.stderr)
        return False

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
    if FLAGS:
        print(f"\n■ 카드 기재값이 좌표와 불일치 — 재확인 대상 {len(FLAGS)}건")
        for mid, nm, tag, why, dl in FLAGS:
            bits = []
            if why:
                bits.append(f"방위각 {why}")
            if dl:
                bits.append(f"거리 {dl}")
            print(f"   {mid:8} {nm[:14]:15} {tag}  " + " · ".join(bits))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
