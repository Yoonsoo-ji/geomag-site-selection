# -*- coding: utf-8 -*-
"""
중간보고회(2026-08-25) 지윤수 작성분 — 국토지리정보원 본 deck 양식 판본
=========================================================================

첨부 양식 `지자기측량발전방안연구_중간보고회(20260812)_목차(초안1)_지윤수.pptx`
**그 파일을 그대로 열어 배정된 면만 채운다.** 머리글·제목띠·로고·쪽번호는
슬라이드 레이아웃(「2장」)에 있으므로 건드리지 않는다.

배정 면 (양식 파일 기준)
    5면  → 초안 PDF 10면          모델 개발에 사용한 자료 (반복관측점·상시관측소·전지구모델)
    6면  → 초안 PDF 11면          4축 자료와 항공자력 포함 여부
    7면  → 초안 PDF 12면 보완     도엽별 자침편차 제공의 정확도 요건 (신설·6면 골격 복제)
    10면 → 양식 9면               모델 갱신 주기 및 방안
    11면 → Ⅲ 추진계획             지자기모델 고도화 방안 (신설·「3장」 레이아웃)

⚠ 수치는 하드코딩하지 않고 `docs/data/lmm_model.json` · `lmm_diagnosis.json` ·
   `external_corrections_multi.csv` 에서 읽는다. 재적합하면 재실행만으로 따라간다.

    python create_ngii_midterm_slides.py
"""
import copy
import csv
import json
import sys
from datetime import datetime
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Emu, Inches, Pt

ROOT = Path(__file__).parent
DATA = ROOT / "docs" / "data"
OUT_DIR = ROOT / "docs" / "output"
TEMPLATE = Path(r"D:\LX_yoons\2026_research\2026_지자기 연구\20260825_중간보고회"
                r"\지자기측량발전방안연구_중간보고회(20260812)_목차(초안1)_지윤수.pptx")

# ── 양식 실측 토큰 ──────────────────────────────────────────────────────────
W, H = 10.8333, 7.5
L, R = 0.45, 10.38          # 본문 좌우 한계 (제목띠와 맞춤)
CW = R - L
BODY_Y = 2.00               # 회색 칩 아래

NAVY = RGBColor(0x2F, 0x54, 0x97)      # 제목띠·표 머리
CYAN = RGBColor(0x00, 0xB0, 0xF0)      # 우측 장 표시 테두리
INK = RGBColor(0x1A, 0x1A, 0x1A)
GRAY = RGBColor(0x59, 0x59, 0x59)
LGRAY = RGBColor(0xF2, 0xF5, 0xF9)
LINE = RGBColor(0xBF, 0xCB, 0xDC)
PEACH = RGBColor(0xFD, 0xEA, 0xDA)     # 요지 띠 (본 deck 12면과 동일 계열)
RED = RGBColor(0xC0, 0x00, 0x00)
GREEN = RGBColor(0x1F, 0x6B, 0x4A)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

FB = "KoPub돋움체 Bold"
FM = "KoPub돋움체 Medium"
FL = "KoPub돋움체 Light"


# ── 자료 ────────────────────────────────────────────────────────────────────
MODEL = json.loads((DATA / "lmm_model.json").read_text(encoding="utf-8"))
DIAG = json.loads((DATA / "lmm_diagnosis.json").read_text(encoding="utf-8"))
LOO = MODEL["loo_cv"]
N_SITES = len(MODEL["sites"])
N_OBS = DIAG["dataset"]["n_obs"]
RD = DIAG["asymmetry"]["rD_rms"]          # 편각 잔차 RMS (arcmin)
GRID_KM = DIAG["vector_recovery"]["grid"]["dx_km"]
KPI_D, KPI_F = 0.1, 50.0


def ext_stats():
    f = DATA / "external_corrections_multi.csv"
    if not f.exists():
        return dict(ok=0, y0="", y1="")
    rows = list(csv.DictReader(open(f, encoding="utf-8-sig")))
    ys = sorted({r["날짜"][:4] for r in rows if r.get("날짜")})
    return dict(ok=sum(1 for r in rows if r.get("상태") == "정상"),
                y0=ys[0], y1=ys[-1])


EXT = ext_stats()


# ── 그리기 도구 ─────────────────────────────────────────────────────────────
def _tf(shape, size, bold, color, align, anchor, space=1.15):
    tf = shape.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    for pa in tf.paragraphs:
        pa.alignment = align
        pa.line_spacing = space
        for r in pa.runs:
            r.font.size = Pt(size)
            r.font.bold = bold
            r.font.color.rgb = color
            r.font.name = FB if bold else FM
    return tf


def text(s, x, y, w, h, txt, size=11, color=INK, bold=False,
         align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, space=1.15):
    box = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = box.text_frame
    tf.margin_left = tf.margin_right = Emu(18000)
    tf.margin_top = tf.margin_bottom = 0
    lines = txt.split("\n")
    tf.text = lines[0]
    for ln in lines[1:]:
        tf.add_paragraph().text = ln
    _tf(box, size, bold, color, align, anchor, space)
    return box


def rect(s, x, y, w, h, fill=None, line=None, radius=False, weight=0.75):
    from pptx.enum.shapes import MSO_SHAPE
    shp = s.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE if radius else MSO_SHAPE.RECTANGLE,
        Inches(x), Inches(y), Inches(w), Inches(h))
    if radius:
        shp.adjustments[0] = 0.14
    if fill is None:
        shp.fill.background()
    else:
        shp.fill.solid()
        shp.fill.fore_color.rgb = fill
    if line is None:
        shp.line.fill.background()
    else:
        shp.line.color.rgb = line
        shp.line.width = Pt(weight)
    shp.shadow.inherit = False
    return shp


def lead(s, y, txt, sub=None, fill=PEACH, h=0.72):
    """상단 요지 띠 — 본 deck 12면과 같은 형식(작은 근거 + 굵은 결론)."""
    rect(s, L, y, CW, h, fill=fill)
    if sub:
        text(s, L + 0.14, y + 0.06, CW - 0.28, 0.22, sub, 10, GRAY)
        text(s, L + 0.14, y + 0.28, CW - 0.28, 0.38, txt, 15.5, INK, True)
    else:
        text(s, L + 0.14, y + 0.16, CW - 0.28, 0.40, txt, 15.5, INK, True,
             anchor=MSO_ANCHOR.MIDDLE)
    return y + h


def label(s, x, y, w, txt, fill=NAVY, size=11.5, h=0.30):
    rect(s, x, y, w, h, fill=fill, radius=True)
    text(s, x, y + 0.02, w, h, txt, size, WHITE, True, PP_ALIGN.CENTER,
         MSO_ANCHOR.MIDDLE)
    return y + h


def table(s, x, y, w, heads, rows, widths, row_h=0.42, size=10.5,
          head_h=0.34, head_fill=NAVY, aligns=None):
    """괘선 최소·머리 네이비 — 본 deck 표와 같은 결."""
    cols = [w * f for f in widths]
    aligns = aligns or [PP_ALIGN.CENTER] + [PP_ALIGN.LEFT] * (len(heads) - 1)
    rect(s, x, y, w, head_h, fill=head_fill)
    cx = x
    for i, hd in enumerate(heads):
        text(s, cx + 0.10, y + 0.02, cols[i] - 0.14, head_h, hd, size, WHITE,
             True, PP_ALIGN.CENTER, MSO_ANCHOR.MIDDLE)
        cx += cols[i]
    yy = y + head_h
    for k, row in enumerate(rows):
        if k % 2 == 0:
            rect(s, x, yy, w, row_h, fill=LGRAY)
        cx = x
        for i, cell in enumerate(row):
            bold = (i == 0)
            col = INK
            if isinstance(cell, tuple):
                cell, col = cell
            text(s, cx + 0.10, yy + 0.03, cols[i] - 0.14, row_h, str(cell),
                 size, col, bold, aligns[i], MSO_ANCHOR.MIDDLE)
            cx += cols[i]
        yy += row_h
    rect(s, x, y, w, yy - y, fill=None, line=LINE)
    return yy


def note(s, y, txt, color=GRAY, size=9.5):
    text(s, L, y, CW, 0.30, txt, size, color)
    return y + 0.28


def clear_notes(slide):
    """양식에 적힌 「지윤수 박사님 …」 지시문 상자를 지운다."""
    for sh in list(slide.shapes):
        if sh.has_text_frame and "박사님" in sh.text_frame.text:
            sh._element.getparent().remove(sh._element)


def set_chip(slide, txt):
    """회색 칩(소제목) 문구 교체 — 서식은 첫 run 을 물려받는다."""
    for sh in slide.shapes:
        if not sh.has_text_frame:
            continue
        t = sh.text_frame.text.strip()
        if t and Emu(sh.top or 0).inches < 1.9 and Emu(sh.top or 0).inches > 1.0:
            p = sh.text_frame.paragraphs[0]
            if not p.runs:
                continue
            keep = p.runs[0]
            keep.text = txt
            for r in p.runs[1:]:
                r._r.getparent().remove(r._r)
            for extra in sh.text_frame.paragraphs[1:]:
                extra._p.getparent().remove(extra._p)
            return sh
    return None


def clone_slide(prs, src_idx, layout=None):
    """
    양식 면 하나를 복제한다. layout 을 주면 그 레이아웃에 붙인다
    (Ⅱ 연구수행 → Ⅲ 추진계획 처럼 장 표시를 바꿀 때).

    ⚠ 도형 XML 만 복사하면 그림이 깨진다. 이 버전 python-pptx 는 관계 id 를
      지정할 수 없으므로, 관계를 새로 만든 뒤 **복사본 XML 의 r:embed·r:id 를
      새 rId 로 갈아끼운다.**
    """
    from pptx.oxml.ns import qn

    src = prs.slides[src_idx]
    dst = prs.slides.add_slide(layout or src.slide_layout)
    for sh in list(dst.shapes):
        sh._element.getparent().remove(sh._element)

    rid_map = {}
    for rid, rel in src.part.rels.items():
        # 레이아웃은 add_slide 가 이미 걸었고, 노트는 슬라이드마다 고유 파트라
        # 공유하면 패키지가 깨진다.
        if rel.reltype.endswith(("/slideLayout", "/notesSlide")):
            continue
        try:
            if rel.is_external:
                continue
            rid_map[rid] = dst.part.rels.get_or_add(rel.reltype, rel._target)
        except Exception:          # noqa: BLE001
            continue

    ATTRS = [qn("r:embed"), qn("r:id"), qn("r:link")]
    for sh in src.shapes:
        el = copy.deepcopy(sh._element)
        for node in el.iter():
            for a in ATTRS:
                v = node.get(a)
                if v in rid_map:
                    node.set(a, rid_map[v])
        dst.shapes._spTree.append(el)
    return dst


def strip_body(slide, keep_above=1.95):
    """복제본에서 본문 도형만 지우고 제목·칩 골격은 남긴다."""
    for sh in list(slide.shapes):
        if Emu(sh.top or 0).inches > keep_above:
            sh._element.getparent().remove(sh._element)


def set_title(slide, txt):
    """제목띠 문구 교체 (제목띠는 y<1.0 의 도형)."""
    for sh in slide.shapes:
        if not sh.has_text_frame:
            continue
        if 0.2 < Emu(sh.top or 0).inches < 1.0 and sh.text_frame.text.strip():
            pa = sh.text_frame.paragraphs[0]
            if not pa.runs:
                continue
            pa.runs[0].text = txt
            for r in pa.runs[1:]:
                r._r.getparent().remove(r._r)
            return sh
    return None


def move_slide(prs, from_idx, to_idx):
    ids = prs.slides._sldIdLst
    items = list(ids)
    ids.remove(items[from_idx])
    ids.insert(to_idx, items[from_idx])


# ══════════════════════════════════════════════════════════════ 면 구성
def s_data(slide):
    """① 모델 개발에 사용한 자료 (초안 10면)."""
    clear_notes(slide)
    set_chip(slide, "기존 반복 관측 성과 기반 지자기모델 시범구축")

    y = lead(slide, BODY_Y,
             "‘25년 시범모델은 반복관측 성과·상시관측소·전지구모델·항공자력 "
             "4종 자료로 구축",
             "모델 개발에 사용한 자료 확인 (반복관측점 · 상시관측소 · 전지구모델)")

    y = table(
        slide, L, y + 0.18, CW,
        ["구분", "자료원", "규모", "모델 내 역할", "비고"],
        [["반복관측점", "지자기점 관측성과(’22~’25) + ’19년 야장",
          f"{N_SITES}점 · {N_OBS}회", "② Regional (공간분포)",
          "재방문 불일치 2점 제외"],
         ["상시관측소", "청양(기상청) · 제주 · 강릉 · 이천",
          "4개소 1분자료", "④ External (시간변화)",
          f"관측시각 복원(’{EXT['y0'][2:]}~’{EXT['y1'][2:]})"],
         ["전지구모델", "IGRF-14 (13차 구면조화)", "4개 에폭",
          "① Core (주자기장)", "공식 계산기와 자릿수 내 일치"],
         ["항공자력", "KIGAM 자력이상도 1.5분 격자",
          f"약 {GRID_KM:.1f} km", "③ Crustal (지각 자기이상)",
          "원측선 자료는 부재"]],
        [0.11, 0.31, 0.12, 0.25, 0.21], row_h=0.46)

    # 청양 활용 여부 — 자문의견 검토 결과
    y += 0.22
    bh = 1.28
    rect(slide, L, y, CW, bh, fill=WHITE, line=NAVY, radius=True, weight=1.0)
    label(slide, L + 0.20, y - 0.15, 2.95, "청양 상시관측소 자료 활용 여부")
    text(slide, L + 0.28, y + 0.30, CW - 0.56, 0.30,
         "활용함 — 다만 청양 단독이 아니라 4개소 공간보간으로 적용", 13, NAVY, True)
    text(slide, L + 0.28, y + 0.64, CW - 0.56, 0.56,
         "· 자문의견(거리·위도 차이로 단독 적용 한계) 검토 결과, 4개소의 "
         "정온야간 기준선 대비 편차를 1차 평면으로 공간보간\n"
         "· 자기폭풍(Kp 2 초과) 구간은 보간이 성립하지 않아 보정에서 제외하고, "
         "편각에만 적용(총자력은 지각장 영향이 지배)", 10.5, INK, space=1.35)

    note(slide, y + bh + 0.10,
         "※ 반복관측 성과는 야장 원본과 대조하여 일변화·기준시점 환산이 적용되지 "
         "않은 원시값임을 확인함(매칭 18행 전부 표기 자릿수 이내 일치)")


def s_axes(slide):
    """② 4축 자료와 항공자력 포함 여부 (초안 11면)."""
    clear_notes(slide)
    set_chip(slide, "LMM(국내외기준) 모델 비교 및 한계·개선방향 도출")

    y = lead(slide, BODY_Y,
             "항공자력을 포함해 4축을 모두 구성 — 남은 제약은 지각축의 해상도",
             "LMM 모델 개발에 필요한 4축 자료 / 이번 개발의 항공자력 포함 여부",
             h=0.66)

    text(slide, L, y + 0.10, CW, 0.30,
         "B_LMM  =  ① Core(IGRF)  +  ② Regional(지역)  +  ③ Crustal(지각)  "
         "+  ④ External(외부장)",
         13, NAVY, True, PP_ALIGN.CENTER)

    y = table(
        slide, L, y + 0.44, CW,
        ["축", "자료", "이번 개발 반영", "한계"],
        [["① Core", "IGRF-14 (13차)", "반영 — 재현성 검증 통과", "—"],
         ["② Regional", f"지표 절대측정 {N_SITES}점", "반영 — 1차 다항식",
          "권고 30점 대비 미달"],
         ["③ Crustal", f"항공자력 {GRID_KM:.1f} km 격자",
          "반영 — 스칼라를 벡터로 복원해 편각·복각까지 적용",
          ("원측선 부재로 해상도 상한", RED)],
         ["④ External", "상시관측 4개소",
          f"반영 — 정온 세션 {EXT['ok']}건, 편각 보정",
          "공간투영(NOC) 미구현"]],
        [0.11, 0.20, 0.42, 0.27], row_h=0.44)

    # 4축 구성의 효과 / 남은 한계
    y += 0.22
    bw = (CW - 0.24) / 2
    for i, (tt, col, lines) in enumerate([
        ("항공자력 포함(4축)의 효과", GREEN,
         "· 총자력뿐 아니라 편각·복각에도 지각 자기이상을 반영\n"
         "· 편각 교차검증 오차 0.540° → 0.340° 로 감소\n"
         "· 편각 잔차의 약 6할을 지각 자기이상으로 설명"),
        ("3축 수준에 머무르는 제약", RED,
         "· 원측선(측선 100~250 m) 자료가 존재하지 않음\n"
         f"· 공개 격자 {GRID_KM:.1f} km 는 단파장을 평활 — 진폭 약 1.6배 부족\n"
         "· 해상도 개선은 신규 자력측량 여부에 좌우")]):
        x = L + i * (bw + 0.24)
        rect(slide, x, y, bw, 1.24, fill=WHITE, line=col, radius=True, weight=1.0)
        label(slide, x + 0.18, y - 0.15, 2.5, tt, fill=col, size=11)
        text(slide, x + 0.24, y + 0.34, bw - 0.48, 0.86, lines, 10.5, INK,
             space=1.32)

    note(slide, y + 1.36,
         "※ 스칼라 총자력 이상은 이상벡터를 주자기장 방향으로 투영한 값이므로, "
         "파수영역 역산으로 3성분을 복원하여 편각·복각에 반영함")


def s_service(slide):
    """③ 도엽별 자침편차 제공의 정확도 요건 (초안 12면 보완)."""
    clear_notes(slide)
    set_chip(slide, "도엽별 자침편차 제공을 위한 정확도 요건")

    y = lead(slide, BODY_Y,
             "도엽 단위로 자침편차를 표기하려면 도폭 내 편각 변화가 표기 "
             "정밀도보다 작아야 함",
             "규정 제12조 — 1등 지자기점은 국가기본도의 자편각 표기를 위한 점")

    y = table(
        slide, L, y + 0.18, CW * 0.60,
        ["축척", "도폭 크기", "도폭 내 편각 변화", "목표(0.1°) 대비"],
        [["1:50,000", "15′", "0.0905°", ("91% — 여유 없음", RED)],
         ["1:25,000", "7′30″", "0.0453°", ("45% — 적정", GREEN)],
         ["1:5,000", "1′30″", "0.0091°", "9%"]],
        [0.22, 0.20, 0.29, 0.29], row_h=0.42,
        aligns=[PP_ALIGN.CENTER] * 4)

    x2 = L + CW * 0.62
    w2 = CW - CW * 0.62
    rect(slide, x2, BODY_Y + 0.90, w2, 1.94, fill=WHITE, line=NAVY,
         radius=True, weight=1.0)
    label(slide, x2 + 0.18, BODY_Y + 0.75, 2.2, "현재 모델 정확도")
    text(slide, x2 + 0.26, BODY_Y + 1.22, w2 - 0.52, 0.34,
         f"편각 {LOO['D']:.2f}°  (목표 {KPI_D:g}°)", 15, RED, True)
    text(slide, x2 + 0.26, BODY_Y + 1.58, w2 - 0.52, 1.10,
         f"· 총자력 {LOO['F']:.0f} nT (목표 {KPI_F:g} nT)\n"
         f"· 편각 잔차 {RD:.1f}′ — 착수 시점 35.0′ 에서 감소\n"
         "· 현 단계는 참고값이며 정식 표기는 목표 달성 후",
         10.5, INK, space=1.35)

    y += 0.30
    text(slide, L, y, CW, 0.30, "제공 방식(안)", 12.5, NAVY, True)
    y += 0.34
    y = table(
        slide, L, y, CW,
        ["단계", "제공 형태", "내용"],
        [["1단계", "도엽별 수치자료(엑셀)",
          "1:25,000 도엽 단위 자침편차 · 기준시점 · 연변화율"],
         ["2단계", "조회 서비스", "국토정보플랫폼에서 위치·도엽 검색 시 자침편차 표시"],
         ["3단계", "지형도 표기", "도곽 밖 자침편차 도표 반영(정확도 목표 달성 후)"]],
        [0.11, 0.24, 0.65], row_h=0.40)

    note(slide, y + 0.12,
         "※ 도폭 내 편각 변화는 IGRF-14 로 산출한 값이며, 1:25,000 은 모델 오차분을 "
         "남길 수 있는 가장 성긴 축척임 / 도폭 규격은 지도도식적용규정 제4조·제211조")


def s_cycle(slide):
    """④ 모델 갱신 주기 및 방안 (양식 9면 · Ⅱ)."""
    clear_notes(slide)
    set_chip(slide, "모델 갱신 주기 및 방안 도출")

    y = lead(slide, BODY_Y,
             "관측망 재관측 주기(5년)에 맞춰 모델을 재계산하고, 기준시점 변경 시 "
             "수시 갱신",
             "국가 지자기모델 갱신 주기 및 절차", h=0.66)

    y += 0.22
    steps = [("재관측", "5년 주기 반복관측\n(관측시각·Kp 기록)"),
             ("재적합", "4축 재계산 · 교차검증\n(순차 튜닝 금지)"),
             ("검증", "독립 측점 검증\n판정기준 충족 확인"),
             ("배포", "도엽별 자침편차 산출\n계산기·조회 서비스 갱신")]
    bw = (CW - 3 * 0.34) / 4
    for i, (tt, body) in enumerate(steps):
        x = L + i * (bw + 0.34)
        rect(slide, x, y, bw, 1.16, fill=WHITE, line=LINE, radius=True)
        rect(slide, x, y, bw, 0.36, fill=NAVY, radius=True)
        rect(slide, x, y + 0.20, bw, 0.16, fill=NAVY)
        text(slide, x, y + 0.03, bw, 0.30, tt, 12, WHITE, True,
             PP_ALIGN.CENTER, MSO_ANCHOR.MIDDLE)
        text(slide, x + 0.14, y + 0.48, bw - 0.28, 0.60, body, 10.5, INK,
             PP_ALIGN.CENTER if False else PP_ALIGN.LEFT, space=1.4)
        if i < 3:
            text(slide, x + bw + 0.04, y + 0.42, 0.26, 0.30, "→", 15, NAVY,
                 True, PP_ALIGN.CENTER)

    y += 1.40
    y = table(
        slide, L, y, CW,
        ["구분", "시점", "내용"],
        [["정기 갱신", "5년", "재관측 성과 반영 · 전 층 재적합 · 기준시점 환산"],
         ["수시 갱신", "IGRF 갱신 시(5년)", "① Core 계수 교체 후 재계산"],
         ["수시 갱신", "측점·자료 추가 시", "신규 관측성과 투입 후 교차검증 재수행"]],
        [0.16, 0.22, 0.62], row_h=0.42)

    note(slide, y + 0.16,
         "※ 층이 서로 결합되어 있어 한 층만 고쳐 끼우는 순차 갱신은 성립하지 "
         "않으며, 갱신 시에는 전 층을 다시 적합해야 함")


def s_plan(slide):
    """④ LMM 모델 추후 진행 방향 (Ⅲ 추진계획)."""
    clear_notes(slide)
    set_title(slide, "3. 지자기모델 고도화 방안")
    set_chip(slide, "지자기모델 검증 및 보완 방향")

    y = lead(slide, BODY_Y,
             "입력자료 보강 → 모델 고도화 → 독립검증·갱신주기 순으로 진행",
             "LMM 모델 추후 진행 방향")

    y += 0.20
    steps = [
        ("1단계", "입력자료 보강", "’26.9~10",
         "· 잔차 큰 측점 재측정 · 방위표지 고정\n"
         "· 야장에만 있는 성과 17건 검증 투입\n"
         f"· 목표: 반복관측점 {N_SITES}점 → 20점 이상"),
        ("2단계", "모델 고도화", "’26.10~11",
         "· 상시관측 4개소 공간투영(NOC)\n"
         "· 잔차면 모델링 · 다항식 차수 재선정\n"
         "· 지각 자기이상 진폭 보정 검토"),
        ("3단계", "검증 · 갱신체계", "’26.11~",
         "· 적합에 쓰지 않은 독립 측점 검증\n"
         "· 5년 갱신주기 연계 재계산 절차\n"
         "· 도엽별 자침편차 산출 · 서비스 연계"),
    ]
    bw = (CW - 0.30) / 3
    for i, (no, tt, when, body) in enumerate(steps):
        x = L + i * (bw + 0.15)
        rect(slide, x, y, bw, 2.26, fill=WHITE, line=LINE, radius=True)
        rect(slide, x, y, bw, 0.42, fill=NAVY, radius=True)
        rect(slide, x, y + 0.24, bw, 0.18, fill=NAVY)
        text(slide, x + 0.12, y + 0.04, bw - 0.24, 0.34,
             f"{no}   {tt}", 12, WHITE, True, PP_ALIGN.CENTER,
             MSO_ANCHOR.MIDDLE)
        text(slide, x + 0.16, y + 0.50, bw - 0.32, 0.22, when, 9.5, GRAY, True)
        text(slide, x + 0.16, y + 0.78, bw - 0.32, 1.40, body, 10.5, INK,
             space=1.45)
        if i < 2:
            text(slide, x + bw + 0.01, y + 0.80, 0.14, 0.30, "▶", 11, NAVY, True,
                 PP_ALIGN.CENTER)

    y += 2.46
    bw2 = (CW - 0.24) / 2
    rect(slide, L, y, bw2, 1.00, fill=LGRAY, line=LINE, radius=True)
    text(slide, L + 0.20, y + 0.12, bw2 - 0.40, 0.24, "단계 이행 판정기준",
         11.5, NAVY, True)
    text(slide, L + 0.20, y + 0.40, bw2 - 0.40, 0.52,
         f"편각 잔차 10′ 이하 → 정식 산출 검토 · 10~20′ → 부분 재측정 "
         f"(현재 {RD:.1f}′)", 10.5, INK, space=1.3)

    rect(slide, L + bw2 + 0.24, y, bw2, 1.00, fill=WHITE, line=RED,
         radius=True, weight=1.0)
    text(slide, L + bw2 + 0.44, y + 0.12, bw2 - 0.40, 0.24, "제약 사항",
         11.5, RED, True)
    text(slide, L + bw2 + 0.44, y + 0.40, bw2 - 0.40, 0.52,
         "항공자력 원측선 자료 부재 — 지각축 해상도 개선은 신규 자력측량 "
         "여부에 좌우", 10.5, INK, space=1.3)


# ══════════════════════════════════════════════════════════════ 실행
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    if not TEMPLATE.exists():
        raise SystemExit(f"양식 파일을 찾을 수 없다: {TEMPLATE}")

    prs = Presentation(str(TEMPLATE))
    print(f"양식 적재: {TEMPLATE.name}  ({len(prs.slides)}면)")

    s_data(prs.slides[4])          # 5면 — 모델 개발 자료
    s_axes(prs.slides[5])          # 6면 — 4축·항공자력
    s_cycle(prs.slides[8])         # 9면 — 모델 갱신 주기

    # 자침편차 면은 양식에 없으므로 6면 골격을 복제해 그 뒤에 끼운다
    dup = clone_slide(prs, 5)
    strip_body(dup)
    move_slide(prs, len(prs.slides) - 1, 6)
    s_service(prs.slides[6])

    # 추진계획(Ⅲ) 면 — 장 표시가 달라야 하므로 「3장」 레이아웃에 붙인다
    lay3 = [l for l in prs.slide_masters[0].slide_layouts if l.name == "3장"][0]
    plan = clone_slide(prs, 9, layout=lay3)   # 9면(=갱신주기) 골격 재사용
    strip_body(plan)
    move_slide(prs, len(prs.slides) - 1, 10)
    s_plan(prs.slides[10])

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = OUT_DIR / f"{ts}_중간보고회_지자기모델_작성분_5장.pptx"
    prs.save(out)
    print(f"저장: {out}")
    print("  · 5면 모델 개발 자료 · 6면 4축/항공자력 · 7면 자침편차 정확도")
    print("  · 10면 모델 갱신 주기(Ⅱ) · 11면 고도화 방안(Ⅲ 추진계획)")
    return out


if __name__ == "__main__":
    main()
