# -*- coding: utf-8 -*-
"""
지자기 모델(LMM) 시범 구축(안) — 학술 조판 발표자료 생성.

    python create_lmm_slides.py   -> docs/output/YYYYMMDD_HHMMSS_LMM_시범구축_발표자료.pptx

조판 원칙 (논문 어법)
  · 세리프 전면 사용 — 라틴 Cambria, 한글 Noto Serif KR (run 단위로 latin/ea 분리 지정)
  · 절 번호(1~6) · 번호 붙인 식 · 「표 n」「그림 n」 캡션 체계
  · booktabs 규칙 — 가로 괘선 3줄(top/mid/bottom)만, 세로 괘선·음영 없음
  · 무채색 본문 + 잉크블루 단색 강조. 라운드 카드·배지·그라디언트 일절 배제
  · 표 캡션은 위, 그림 캡션은 아래

도해 위주 구성 — 표제를 뺀 모든 면에 그림 또는 표를 둔다
  그림 1  층 분해 모식도      (matplotlib, 합성 파형)
  그림 2  층별 공간 파장 대역  (IGRF degree 13 -> λ≈3,000 km 계산)
  그림 3  구축 절차 흐름도     (pptx 네이티브 도형 — 편집 가능)
  그림 4  측점 분포           (lmm_model.json + korea_boundary.geojson)
  그림 5  산출 체계           (pptx 네이티브 도형)
  그림 6  요약 모식도         (matplotlib)

그림 1·6 은 개념을 보이는 모식도이므로 캡션에 그 사실을 명시한다 (실측이 아님).
그림 2·4 는 실제 자료·계산에서 생성한다.

내용은 4-층 결합 구조와 시범 구축 6단계 절차를 재구성한 것으로,
현 시점의 적합 결과(LOO 등)는 싣지 않는다 — 이 문서는 결과 보고가 아니라 절차(안)이다.
"""
import json
import math
import tempfile
from datetime import datetime
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

HERE = Path(__file__).parent
DOCS_OUT = HERE / "docs" / "output"
MODEL_JSON = HERE / "docs" / "data" / "lmm_model.json"
BOUNDARY = HERE / "data" / "korea_boundary.geojson"
ASSETS = Path(tempfile.gettempdir()) / "lmm_slide_figs"

# ── 팔레트 (무채색 + 단색 강조) ───────────────────────────────────────────
INK      = RGBColor(0x1A, 0x1A, 0x1A)
BODY     = RGBColor(0x3A, 0x3A, 0x3A)
MUTED    = RGBColor(0x70, 0x70, 0x70)
ACCENT   = RGBColor(0x1F, 0x3A, 0x5F)   # 잉크블루
MARK     = RGBColor(0x8C, 0x2F, 0x27)   # 강조용 적갈
RULE     = RGBColor(0x1A, 0x1A, 0x1A)
RULE_LT  = RGBColor(0xB8, 0xB8, 0xB8)
BOXLINE  = RGBColor(0x5A, 0x6B, 0x80)
BOXTINT  = RGBColor(0xF2, 0xF4, 0xF7)
WHITE    = RGBColor(0xFF, 0xFF, 0xFF)

HEX_ACCENT = "#1F3A5F"
HEX_MARK = "#8C2F27"

LATIN = "Cambria"
KR = "Noto Serif KR"

SW, SH = 13.333, 7.5
M = 1.00
CW = SW - 2 * M
MEASURE = 7.9


# ══════════════════════════════════════════════════════════════════════
# 도판 (matplotlib)
# ══════════════════════════════════════════════════════════════════════
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    # Noto Serif KR 은 가변 글꼴(VF)이라 matplotlib 의 기본 스캔에 잡히지 않는다.
    # 파일 경로로 직접 등록해야 본문 조판과 도판의 글꼴이 일치한다.
    for path in (r"C:\Windows\Fonts\NotoSerifKR-VF.ttf",
                 r"C:\Windows\Fonts\NanumMyeongjo.ttf"):
        if Path(path).exists():
            try:
                font_manager.fontManager.addfont(path)
            except Exception:
                pass

    have = {f.name for f in font_manager.fontManager.ttflist}
    for cand in (KR, "NanumMyeongjo", "Batang", "serif"):
        if cand in have or cand == "serif":
            plt.rcParams["font.family"] = cand
            break
    plt.rcParams.update({
        "axes.unicode_minus": False,
        "axes.edgecolor": "#4A4A4A",
        "axes.linewidth": 0.7,
        "xtick.color": "#4A4A4A",
        "ytick.color": "#4A4A4A",
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })
    return plt


def fig_decomposition(path):
    """그림 1 — 층 분해 모식도. 파장이 다른 성분이 합쳐지는 구조를 보인다."""
    plt = _mpl()
    import numpy as np

    # 진폭은 파장이 짧아질수록 작아진다 — 이래야 합에서 장파장 추세가 살아난다.
    # (실제로도 주자기장 수만 nT, 지역 수백 nT, 지각 이상 수십~수백 nT 규모)
    x = np.linspace(0, 1, 800)
    core = 5.0 + 3.6 * x + 0.45 * np.sin(2 * np.pi * 0.6 * x)
    regional = (1.05 * np.sin(2 * np.pi * 1.3 * x + 0.5)
                + 0.45 * np.sin(2 * np.pi * 2.2 * x))
    crustal = (0.40 * np.sin(2 * np.pi * 11 * x)
               + 0.25 * np.sin(2 * np.pi * 23 * x + 1.0)
               + 0.15 * np.sin(2 * np.pi * 41 * x + 2.2))
    total = core + regional + crustal

    rows = [("① Core", core, 0.62), ("② Regional", regional, 0.80),
            ("③ Crustal", crustal, 1.00), ("합  ①+②+③", total, 1.00)]

    fig, axes = plt.subplots(len(rows), 1, figsize=(10.6, 3.05), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1, 1, 1.45]})
    for ax, (label, y, alpha) in zip(axes, rows):
        is_sum = label.startswith("합")
        ax.plot(x, y, color=HEX_ACCENT, lw=1.7 if is_sum else 1.3, alpha=alpha)
        ax.set_ylabel(label, rotation=0, ha="right", va="center",
                      labelpad=12, fontsize=9.5,
                      color="#1A1A1A" if is_sum else "#3A3A3A")
        ax.set_yticks([])
        ax.margins(y=0.28)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color("#DDDDDD")
        ax.tick_params(axis="x", length=0)

    # 합 패널 위 가로선 = 총합 기호의 가로줄 역할
    axes[-1].spines["top"].set_visible(True)
    axes[-1].spines["top"].set_color("#4A4A4A")
    axes[-1].spines["bottom"].set_color("#4A4A4A")
    axes[-1].set_xlabel("거리 (모식적 표현, 눈금 없음)", fontsize=9.5, labelpad=6)
    axes[-1].set_xticks([])
    fig.subplots_adjust(left=0.105, right=0.99, top=0.99, bottom=0.14, hspace=0.42)
    fig.savefig(path, dpi=300)
    plt.close(fig)


def fig_wavelength(path):
    """그림 2 — 층별 공간 파장 대역. IGRF 차단차수에서 파장을 직접 계산한다."""
    plt = _mpl()

    R = 6371.0
    n_core = 13
    lam_core = 2 * math.pi * R / math.sqrt(n_core * (n_core + 1))
    lam_round = round(lam_core, -2)          # 2,967 -> 3,000 (유효숫자 정리)

    bands = [
        # alpha 하한 0.70 — 그 아래로는 막대 위 흰 라벨의 명도 대비가 부족하다
        ("③ Crustal", 0.05, 50.0, "0.05 – 50 km", 1.00),
        ("② Regional", 50.0, lam_core, "50 – 3,000 km", 0.85),
        ("① Core", lam_core, 20000.0, "3,000 km 이상", 0.70),
    ]
    fig, ax = plt.subplots(figsize=(9.4, 2.15))
    for i, (label, lo, hi, tag, alpha) in enumerate(bands):
        ax.barh(i, hi - lo, left=lo, height=0.46,
                color=HEX_ACCENT, alpha=alpha, edgecolor="none")
        ax.text(math.sqrt(lo * hi), i, tag,
                ha="center", va="center", color="white", fontsize=9)

    ax.set_xscale("log")
    ax.set_xlim(0.02, 30000)
    ax.set_ylim(-0.7, 3.25)
    ax.set_yticks(range(len(bands)))
    ax.set_yticklabels([b[0] for b in bands], fontsize=10)
    ax.set_xlabel("공간 파장 (km, 로그 눈금)", fontsize=10, labelpad=7)
    ax.axvline(lam_core, color=HEX_MARK, lw=0.9, ls=(0, (4, 3)), ymax=0.86)
    # 「차단차수」는 발표 자리에서 걸리는 용어라 평이한 말을 앞세우고
    # 전문 용어는 괄호로만 남긴다
    ax.annotate(f"IGRF 표현 한계 ≈ {lam_round:,.0f} km  (차수 {n_core})",
                xy=(lam_core, 2.62), fontsize=8.5, color=HEX_MARK,
                ha="right", va="center",
                xytext=(-6, 0), textcoords="offset points")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#DDDDDD", lw=0.5, which="major")
    ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.135, right=0.985, top=0.98, bottom=0.26)
    fig.savefig(path, dpi=300)
    plt.close(fig)


def fig_sites(path):
    """그림 4 — 지표 절대측정 측점 분포. 원 크기는 반복 관측 횟수."""
    plt = _mpl()
    model = json.loads(MODEL_JSON.read_text(encoding="utf-8"))
    sites = model["sites"]
    gj = json.loads(BOUNDARY.read_text(encoding="utf-8"))

    rings = []
    for feat in gj["features"]:
        geom = feat["geometry"]
        polys = (geom["coordinates"] if geom["type"] == "MultiPolygon"
                 else [geom["coordinates"]])
        for poly in polys:
            rings.append(poly[0])

    fig, ax = plt.subplots(figsize=(4.3, 5.0))
    for ring in rings:
        ax.fill([p[0] for p in ring], [p[1] for p in ring],
                facecolor="#F0F0F0", edgecolor="#9A9A9A", lw=0.5)

    ax.scatter([s["lon"] for s in sites], [s["lat"] for s in sites],
               s=[16 + 13 * s["n_visit"] for s in sites], facecolor=HEX_ACCENT,
               edgecolor="white", linewidth=0.7, zorder=5)

    ax.set_xlim(125.0, 130.2)
    ax.set_ylim(33.0, 38.8)
    ax.set_aspect(1.0 / math.cos(math.radians(36.0)))
    ax.set_xlabel("경도 (°E)", fontsize=9.5, labelpad=5)
    ax.set_ylabel("위도 (°N)", fontsize=9.5, labelpad=3)
    ax.grid(color="#E4E4E4", lw=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    for v in (2, 4):
        ax.scatter([], [], s=16 + 13 * v, facecolor=HEX_ACCENT,
                   edgecolor="white", linewidth=0.7, label=f"{v}회")
    leg = ax.legend(title="반복 관측", fontsize=8.5, title_fontsize=8.5,
                    loc="upper left", frameon=False, labelspacing=0.9,
                    borderpad=0.2)   # 좌하단은 제주가 있어 범례를 겹치게 한다
    leg._legend_box.align = "left"
    fig.subplots_adjust(left=0.19, right=0.99, top=0.99, bottom=0.11)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return len(sites)


def fig_summary(path):
    """그림 6 — 분리·결합·검증 모식도."""
    plt = _mpl()
    import numpy as np

    x = np.linspace(0, 1, 500)
    a = 1.05 * np.sin(2 * np.pi * 1.1 * x + 0.4)
    b = 0.62 * np.sin(2 * np.pi * 4.5 * x)
    c = 0.34 * np.sin(2 * np.pi * 13 * x + 1.1)

    fig, axes = plt.subplots(1, 3, figsize=(10.6, 1.95))

    # 분리 — 성분을 따로 그린다
    for k, (y, off) in enumerate(((a, 2.1), (b, 0.0), (c, -1.9))):
        axes[0].plot(x, y + off, color=HEX_ACCENT, lw=1.3,
                     alpha=[1.0, 0.75, 0.5][k])
    axes[0].set_title("분 리", fontsize=10.5, pad=9, color=INK_HEX)

    # 결합 — 합쳐 하나의 면으로
    axes[1].plot(x, a + b + c, color=HEX_ACCENT, lw=1.8)
    axes[1].fill_between(x, a + b + c, -2.6, color=HEX_ACCENT, alpha=0.07)
    axes[1].set_ylim(-2.6, 2.6)
    axes[1].set_title("결 합", fontsize=10.5, pad=9, color=INK_HEX)

    # 검증 — 한 점을 빼고 나머지로 예측
    rng = np.random.default_rng(7)
    px = rng.uniform(0.08, 0.92, 12)
    py = rng.uniform(-1.5, 1.5, 12)
    axes[2].scatter(px, py, s=22, facecolor=HEX_ACCENT, edgecolor="none")
    axes[2].scatter([0.5], [0.25], s=90, facecolor="white",
                    edgecolor=HEX_MARK, linewidth=1.5, zorder=5)
    axes[2].annotate("제외한 1점", xy=(0.5, 0.25), xytext=(0.5, -2.05),
                     fontsize=8.5, color=HEX_MARK, ha="center",
                     arrowprops=dict(arrowstyle="->", color=HEX_MARK, lw=0.9))
    axes[2].set_ylim(-2.6, 2.6)
    axes[2].set_xlim(0, 1)
    axes[2].set_title("검 증", fontsize=10.5, pad=9, color=INK_HEX)

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(False)
    # 아래 3단 본문과 패널을 정확히 맞춘다 (본문 단 간격 0.5", 폭 CW)
    gap_frac = 0.5 / CW
    fig.subplots_adjust(left=0.0, right=1.0, top=0.84, bottom=0.04,
                        wspace=gap_frac / ((1 - 2 * gap_frac) / 3))
    fig.savefig(path, dpi=300)
    plt.close(fig)


INK_HEX = "#1A1A1A"


def build_figures():
    ASSETS.mkdir(parents=True, exist_ok=True)
    paths = {
        "decomp": ASSETS / "fig1_decomposition.png",
        "wave": ASSETS / "fig2_wavelength.png",
        "sites": ASSETS / "fig4_sites.png",
        "summary": ASSETS / "fig6_summary.png",
    }
    fig_decomposition(paths["decomp"])
    fig_wavelength(paths["wave"])
    n = fig_sites(paths["sites"])
    fig_summary(paths["summary"])
    return paths, n


# ══════════════════════════════════════════════════════════════════════
# 조판 헬퍼
# ══════════════════════════════════════════════════════════════════════
def font_pair(run, latin=LATIN, ea=KR):
    """라틴·한글에 서로 다른 세리프를 물린다 (Cambria 에 한글 글리프가 없다)."""
    run.font.name = latin                      # <a:latin> 을 스키마 위치에 삽입
    rPr = run.font._rPr
    lat = rPr.find(qn("a:latin"))
    ea_el = rPr.find(qn("a:ea"))
    if ea_el is None:
        ea_el = rPr.makeelement(qn("a:ea"), {"typeface": ea})
        lat.addnext(ea_el)
    else:
        ea_el.set("typeface", ea)


def textbox(slide, x, y, w, h, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = anchor
    tf.paragraphs[0].alignment = align
    return tf


def para(tf, text, size, color, *, bold=False, italic=False, line=1.55,
         space_before=0, space_after=0, align=None, first=False):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    if align is not None:
        p.alignment = align
    p.line_spacing = line
    p.space_before = Pt(space_before)
    p.space_after = Pt(space_after)
    r = p.add_run()
    r.text = text
    r.font.size = Pt(size)
    r.font.bold = bold
    r.font.italic = italic
    r.font.color.rgb = color
    font_pair(r)
    return p


def run_in(p, text, size, color, *, bold=False, italic=False, sub=False):
    r = p.add_run()
    r.text = text
    r.font.size = Pt(size)
    r.font.bold = bold
    r.font.italic = italic
    r.font.color.rgb = color
    font_pair(r)
    if sub:
        r.font._rPr.set("baseline", "-25000")
    return r


def rule(slide, x, y, w, weight=0.75, color=RULE):
    """가로 괘선 — booktabs 의 toprule/midrule/bottomrule."""
    ln = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x),
                                    Inches(y), Inches(x + w), Inches(y))
    ln.line.color.rgb = color
    ln.line.width = Pt(weight)
    return ln


def arrow(slide, x1, y1, x2, y2, *, color=BOXLINE, weight=1.0):
    """화살표 — 흐름도 연결선. 커넥터에 tailEnd 를 직접 물린다."""
    cn = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1),
                                    Inches(y1), Inches(x2), Inches(y2))
    cn.line.color.rgb = color
    cn.line.width = Pt(weight)
    ln = cn.line._get_or_add_ln()
    tail = ln.makeelement(qn("a:tailEnd"),
                          {"type": "triangle", "w": "med", "len": "med"})
    ln.append(tail)
    return cn


def flowbox(slide, x, y, w, h, *, fill=WHITE, line=BOXLINE, weight=0.9):
    """흐름도 상자 — 직각 모서리, 얇은 테두리, 그림자 없음."""
    sh = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y),
                                Inches(w), Inches(h))
    sh.fill.solid()
    sh.fill.fore_color.rgb = fill
    sh.line.color.rgb = line
    sh.line.width = Pt(weight)
    sh.shadow.inherit = False
    style = sh._element.find(qn("p:style"))
    if style is not None:
        sh._element.remove(style)          # 테마 effectRef 까지 떼어야 그림자가 사라진다
    sh.text_frame.word_wrap = True
    return sh


def section_head(slide, number, title, *, y=0.68):
    tf = textbox(slide, M, y, 0.62, 0.5)
    para(tf, number, 25, ACCENT, bold=True, line=1.0, first=True)
    tf2 = textbox(slide, M + 0.62, y, CW - 0.62, 0.5)
    para(tf2, title, 25, INK, bold=True, line=1.0, first=True)


def caption(slide, x, y, w, text, *, align=PP_ALIGN.LEFT, h=0.3):
    tf = textbox(slide, x, y, w, h, align=align)
    p = tf.paragraphs[0]
    p.line_spacing = 1.35
    head, _, rest = text.partition(". ")
    run_in(p, head + ".", 10, MUTED, bold=True)
    run_in(p, " " + rest, 10, MUTED)
    return tf


def footnote(slide, text, y=6.74):
    tf = textbox(slide, M, y, CW, 0.34)
    para(tf, text, 9.5, MUTED, line=1.4, first=True)


def booktabs(slide, x, y, w, fracs, header, rows, *, size=11.5,
             head_h=0.40, row_h=0.40, aligns=None):
    """가로 괘선 3줄만 쓰는 학술 표. 세로 괘선·행 음영 없음."""
    aligns = aligns or ["l"] * len(fracs)
    amap = {"l": PP_ALIGN.LEFT, "r": PP_ALIGN.RIGHT, "c": PP_ALIGN.CENTER}
    xs, acc = [], 0.0
    for f in fracs:
        xs.append(x + acc * w)
        acc += f

    rule(slide, x, y, w, weight=1.25)
    for i, cell in enumerate(header):
        tf = textbox(slide, xs[i], y + 0.09, fracs[i] * w - 0.12,
                     head_h, align=amap[aligns[i]])
        para(tf, cell, size, INK, bold=True, line=1.2, first=True)

    ymid = y + head_h + 0.06
    rule(slide, x, ymid, w, weight=0.6)

    yy = ymid
    for row in rows:
        for i, cell in enumerate(row):
            tf = textbox(slide, xs[i], yy + 0.11, fracs[i] * w - 0.12,
                         row_h, align=amap[aligns[i]])
            para(tf, cell, size, BODY, line=1.2, first=True)
        yy += row_h

    rule(slide, x, yy + 0.08, w, weight=1.25)
    return yy + 0.08


def blank(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


# ══════════════════════════════════════════════════════════════════════
# 슬라이드
# ══════════════════════════════════════════════════════════════════════
def s1_title(prs, n_sites):
    s = blank(prs)

    tf = textbox(s, M, 1.62, CW * 0.86, 0.34)
    para(tf, "지구자기장 기준장 모형", 12.5, MUTED, line=1.3, first=True)

    tf = textbox(s, M, 2.10, CW * 0.88, 1.5)
    para(tf, "지자기 모델(LMM) 시범 구축 절차(안)", 36, INK, bold=True,
         line=1.22, first=True)

    tf = textbox(s, M, 3.30, MEASURE, 0.4)
    para(tf, "IGRF · Regional · Crustal · External 4층 결합 모형", 15, ACCENT,
         line=1.4, first=True)

    rule(s, M, 3.96, MEASURE, weight=0.75, color=RULE_LT)

    tf = textbox(s, M, 4.24, MEASURE, 1.5)
    para(tf, "개요", 11.5, INK, bold=True, line=1.3, first=True)
    para(tf, "지구자기장을 파장과 시간 규모가 다른 네 개의 물리 성분으로 "
             "분리하여 각각 산출한 뒤 다시 결합하고, 관측점 잔차면으로 "
             f"층 사이의 공간 편차를 보정한다. 지표 절대측정 {n_sites}개 측점을 "
             "기준으로 Leave-One-Out 교차검증을 수행하여 채택 여부를 판단하며, "
             "최종 산출물은 도엽별 자침편차와 웹·엑셀 계산기이다.",
         12.5, BODY, line=1.62, space_before=6)

    tf = textbox(s, M, 6.42, MEASURE, 0.3)
    para(tf, "2026년 7월", 11.5, MUTED, line=1.3, first=True)

    s.notes_slide.notes_text_frame.text = (
        "이 문서는 결과 보고가 아니라 절차(안)이다. 현 시점 적합 결과는 싣지 않는다.")
    return s


def s2_structure(prs, figs):
    s = blank(prs)
    section_head(s, "1", "모형 구조")

    tf = textbox(s, M, 1.42, CW, 0.4)
    para(tf, "발생 원인이 다른 네 개 장(場)의 합으로 기술한다. 층마다 담당하는 "
             "공간 파장이 달라 서로 섞이지 않는다.",
         12.5, BODY, line=1.6, first=True)

    # 번호 붙인 표시식
    tfe = textbox(s, M, 2.06, CW - 0.7, 0.5, align=PP_ALIGN.CENTER)
    p = tfe.paragraphs[0]
    p.line_spacing = 1.2
    run_in(p, "B", 19, INK, bold=True, italic=True)
    run_in(p, "LMM", 12, INK, sub=True)
    run_in(p, "(r, t)  =  ", 19, INK)
    for k, name in enumerate(("IGRF", "Regional", "Crustal", "External")):
        run_in(p, "B", 19, ACCENT, bold=True, italic=True)
        run_in(p, name, 12, ACCENT, sub=True)
        if k < 3:
            run_in(p, "  +  ", 19, INK)
    # 19pt 식과 12.5pt 번호의 베이스라인을 맞춘다
    tfn = textbox(s, SW - M - 0.5, 2.16, 0.5, 0.4, align=PP_ALIGN.RIGHT)
    para(tfn, "(1)", 12.5, MUTED, line=1.2, first=True)

    s.shapes.add_picture(str(figs["decomp"]), Inches(M), Inches(2.80),
                         width=Inches(CW))
    caption(s, M, 6.14, CW,
            "그림 1. 층 분해 모식도. 파장이 짧을수록 진폭이 작아 합에서는 "
            "장파장 성분이 지배한다. 개념도이며 실측 자료가 아니다.")
    footnote(s, "④ External 은 공간이 아니라 시간에 따라 변하는 항이므로 "
                "그림 1 에 함께 싣지 않았다.", y=6.62)
    return s


def s3_layers(prs, figs):
    s = blank(prs)
    section_head(s, "2", "층별 구성과 자료원")

    tf = textbox(s, M, 1.40, CW, 0.34)
    para(tf, "표 1. 층별 자료원과 담당 대역", 11.5, INK, bold=True,
         line=1.3, first=True)

    booktabs(
        s, M, 1.76, CW, [0.155, 0.245, 0.20, 0.175, 0.225],
        ["층", "자료원", "공간 파장", "시간 규모", "비고"],
        [["① Core", "IGRF-14", "3,000 km 이상", "5년 주기",
          "차수 13 이하 전지구 기준장"],
         ["② Regional", "지표 절대측정", "50 – 3,000 km", "수년",
          "반복관측점 기반 지역 편차"],
         ["③ Crustal", "KIGAM 항공자력", "0.05 – 50 km", "정적",
          "암상·지질구조 기원 국소 이상"],
         ["④ External", "상시관측(청양, CYG)", "전역 근사", "1분 – 수시간",
          "전리층·자기권 기원 시간변화"]],
        size=11.5, row_h=0.44)

    s.shapes.add_picture(str(figs["wave"]), Inches(M + 0.35), Inches(4.18),
                         width=Inches(CW - 0.7))
    caption(s, M, 6.66, CW,
            "그림 2. 층별 공간 파장 대역. 차단차수 13 에 대응하는 파장 "
            "약 3,000 km 가 ①과 ②의 경계이다.")
    return s


def s4_process(prs):
    """그림 3 — 절차 흐름도. 네이티브 도형이라 발표자가 직접 편집할 수 있다."""
    s = blank(prs)
    section_head(s, "3", "시범 구축 절차")

    tf = textbox(s, M, 1.40, CW, 0.34)
    para(tf, "여섯 단계로 구성하며, 제5단계에서 기준에 미달하면 제1단계로 되돌아간다.",
         12.5, BODY, line=1.6, first=True)

    steps = [
        ("제1단계", "입력자료\n수집", "IGRF-14 · 절대측정 17점\n항공자력 · CYG 1분"),
        ("제2단계", "층별 성분\n산출", "기준장 · 지역 · 지각\n시간변화로 분리"),
        ("제3단계", "결합 및\n잔차면 구성", "식 (1) 결합 후\n관측점 잔차 내삽"),
        ("제4단계", "최적계산", "편각 D · 복각 I\n총자력 F 산출"),
        ("제5단계", "교차검증", "Leave-One-Out\n표 2 기준 점검"),
        ("제6단계", "산출물 작성", "도엽별 자침편차\n웹 · 엑셀 계산기"),
    ]
    bw, gap = 1.62, 0.32
    by, bh = 2.22, 2.32
    centers = []
    for i, (num, title, desc) in enumerate(steps):
        x = M + i * (bw + gap)
        centers.append(x + bw / 2)
        last = (i == len(steps) - 1)
        flowbox(s, x, by, bw, bh, fill=BOXTINT if last else WHITE)

        tfn = textbox(s, x + 0.14, by + 0.20, bw - 0.28, 0.24,
                      align=PP_ALIGN.CENTER)
        para(tfn, num, 9, ACCENT, bold=True, line=1.15, first=True)

        tft = textbox(s, x + 0.12, by + 0.56, bw - 0.24, 0.66,
                      align=PP_ALIGN.CENTER)
        para(tft, title, 12, INK, bold=True, line=1.22, first=True)

        tfd = textbox(s, x + 0.12, by + 1.46, bw - 0.24, 0.72,
                      align=PP_ALIGN.CENTER)
        para(tfd, desc, 9, BODY, line=1.4, first=True)

        if i < len(steps) - 1:
            arrow(s, x + bw + 0.06, by + bh / 2, x + bw + gap - 0.06, by + bh / 2)

    # 되먹임 경로 — 제5단계에서 제1단계로
    fy = by + bh + 0.62
    rule(s, centers[0], fy, centers[4] - centers[0], weight=1.0, color=MARK)
    arrow(s, centers[0], fy, centers[0], by + bh + 0.08, color=MARK)
    ln = s.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(centers[4]),
                                Inches(by + bh + 0.08), Inches(centers[4]),
                                Inches(fy))
    ln.line.color.rgb = MARK
    ln.line.width = Pt(1.0)

    tf = textbox(s, centers[0] + 0.25, fy + 0.10, 5.6, 0.3)
    para(tf, "기준 미달 시 입력자료를 보강하여 재수행", 10, MARK, line=1.3,
         first=True)

    caption(s, M, fy + 0.62, CW,
            "그림 3. 시범 구축 절차 흐름도. 붉은 경로는 기준 미달 시의 되먹임을 "
            "나타낸다.")
    return s


def s5_validation(prs, figs, n_sites):
    s = blank(prs)
    section_head(s, "4", "검증 설계")

    tw = 6.55
    tf = textbox(s, M, 1.42, tw, 0.8)
    para(tf, f"측점 {n_sites}개를 대상으로 Leave-One-Out 교차검증을 수행한다. "
             "측점을 하나씩 제외하고 나머지로 재적합한 뒤, 제외한 측점의 "
             "예측값을 실측과 대조한다.",
         12.5, BODY, line=1.6, first=True)

    tf = textbox(s, M, 2.72, tw, 0.34)
    para(tf, "표 2. 성분별 점검 기준", 11.5, INK, bold=True, line=1.3, first=True)

    booktabs(
        s, M, 3.08, tw, [0.26, 0.24, 0.50],
        ["성분", "기준", "근거"],
        [["편각 D", "< 0.1°", "도폭 단위 자침편차의 단일 값 표기"],
         ["총자력 F", "< 50 nT", "산출값이 실측을 대신할 수 있는 범위"],
         ["복각 I", "기준 없음", "참고 지표로 병기"]],
        size=11.5, row_h=0.44)

    s.shapes.add_picture(str(figs["sites"]), Inches(M + tw + 0.95),
                         Inches(1.38), height=Inches(4.72))
    caption(s, M + tw + 0.75, 6.22, CW - tw - 0.75,
            f"그림 4. 지표 절대측정 {n_sites}개 측점 분포. "
            "원의 크기는 반복 관측 횟수이다.")

    footnote(s, "기준값은 공학적 목표치이며 「지구물리측량 작업규정」 제20조의 "
                "측정오차 한계(정수차 30′)와는 성격이 다르다.", y=6.74)
    return s


def s6_outputs(prs):
    """그림 5 — 산출 체계. 하나의 모형 파일에서 두 산출물이 파생됨을 보인다."""
    s = blank(prs)
    section_head(s, "5", "산출물")

    tf = textbox(s, M, 1.40, CW, 0.34)
    para(tf, "두 산출물은 동일한 모형 파일을 참조하므로 재적합 시 함께 갱신된다.",
         12.5, BODY, line=1.6, first=True)

    # 원천 — 모형 파일
    src_w, src_h = 3.05, 1.02
    src_x, src_y = M + 0.30, 3.10
    flowbox(s, src_x, src_y, src_w, src_h, fill=BOXTINT)
    tf = textbox(s, src_x + 0.14, src_y + 0.20, src_w - 0.28, 0.66,
                 align=PP_ALIGN.CENTER)
    para(tf, "LMM 모형 파일", 13, INK, bold=True, line=1.2, first=True)
    para(tf, "층별 계수 · 잔차면", 9.5, MUTED, line=1.3, space_before=3)

    # 파생 — 두 산출물
    out_x = src_x + src_w + 1.45
    out_w, out_h = 5.55, 1.34
    outs = [
        ("가.", "도엽별 자침편차", 2.00,
         "지형도 도폭 단위 단일 값 산출\n도폭 내 편각 변화가 표기 정밀도보다 작아야 성립"),
        ("나.", "웹 · 엑셀 계산기", 3.88,
         "좌표·표고 입력 → 편각 D · 복각 I · 총자력 F\n웹은 단일 파일, 엑셀은 수식 공개형"),
    ]
    for mark, title, oy, desc in outs:
        flowbox(s, out_x, oy, out_w, out_h)
        tfm = textbox(s, out_x + 0.22, oy + 0.20, 0.42, 0.28)
        para(tfm, mark, 11.5, ACCENT, bold=True, line=1.2, first=True)
        tft = textbox(s, out_x + 0.66, oy + 0.18, out_w - 0.9, 0.32)
        para(tft, title, 14, INK, bold=True, line=1.2, first=True)
        tfd = textbox(s, out_x + 0.66, oy + 0.62, out_w - 0.9, 0.82)
        para(tfd, desc, 10.5, BODY, line=1.48, first=True)

        arrow(s, src_x + src_w + 0.10, src_y + src_h / 2,
              out_x - 0.10, oy + out_h / 2)

    caption(s, M, 6.00, CW,
            "그림 5. 산출 체계. 두 산출물이 하나의 모형 파일에서 파생되므로 "
            "수치를 개별 관리하지 않는다.")
    footnote(s, "입력 높이는 표고로 통일한다. 타원체고와의 차이는 지오이드고에 "
                "해당하나 완만하여 Regional 상수항에 흡수된다.", y=6.52)
    return s


def s7_summary(prs, figs):
    s = blank(prs)
    section_head(s, "6", "요약")

    s.shapes.add_picture(str(figs["summary"]), Inches(M), Inches(1.42),
                         width=Inches(CW))
    caption(s, M, 3.62, CW,
            "그림 6. 분리 · 결합 · 검증의 세 단계 모식도.")

    points = [
        ("층 분리", "파장과 시간 규모가 다른 성분을 섞지 않는다. 정밀도의 제약 "
                    "요인이 어느 층에 있는지 식별하기 위한 구조이다."),
        ("결합과 보정", "식 (1)로 결합한 뒤 남는 관측점 잔차를 내삽하여, "
                        "단순 내삽으로는 닿지 않는 공간 정합성을 확보한다."),
        ("검증", "기준 미달 시 모형 구조가 아니라 입력자료의 밀도와 관측시각 "
                 "기록을 먼저 점검한다."),
    ]
    colw = (CW - 2 * 0.5) / 3
    for i, (title, desc) in enumerate(points):
        x = M + i * (colw + 0.5)
        tfn = textbox(s, x, 4.22, colw, 0.3)
        para(tfn, f"{i + 1})  {title}", 13, INK, bold=True, line=1.25,
             first=True)
        tfd = textbox(s, x, 4.66, colw, 1.3)
        para(tfd, desc, 11, BODY, line=1.55, first=True)

    rule(s, M, 6.64, CW, weight=0.5, color=RULE_LT)
    tf = textbox(s, M, 6.80, CW, 0.3)
    para(tf, "지자기 모델(LMM) 시범 구축 절차(안) · 2026년 7월", 9.5, MUTED,
         line=1.3, first=True)
    return s


def main():
    figs, n_sites = build_figures()

    prs = Presentation()
    prs.slide_width = Inches(SW)
    prs.slide_height = Inches(SH)

    s1_title(prs, n_sites)
    s2_structure(prs, figs)
    s3_layers(prs, figs)
    s4_process(prs)
    s5_validation(prs, figs, n_sites)
    s6_outputs(prs)
    s7_summary(prs, figs)

    DOCS_OUT.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = DOCS_OUT / f"{stamp}_LMM_시범구축_발표자료.pptx"
    prs.save(out)
    print("saved:", out)
    return out


if __name__ == "__main__":
    main()
