# -*- coding: utf-8 -*-
"""
현장 야장 기입 안내서 — `docs/output/*_시험탐사_현장야장_기입안내.docx`
========================================================================

    python create_field_guide.py

`make_trial_field_card.py` 가 내는 현장 카드를 «무엇을 어떻게 채우는가» 설명하는
인쇄용 안내서다. 카드와 함께 들고 다닐 수 있도록 A4 몇 장으로 끝낸다.

⚠️ **수치와 규격은 여기서 다시 적지 않는다.** 측선 간격·판정 참고값·발주자 결정은
전부 `trial_survey_spec.py` 에서 읽어 온다. 안내서와 카드가 갈라지면 현장에서
어느 쪽을 믿어야 할지 알 수 없게 된다.

말투는 카드와 같이 **설명하는 투**로 쓴다 — 「~한다」로 끊지 않는다.
"""
from __future__ import annotations

import datetime as dt
import sys
from pathlib import Path

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

import trial_survey_spec as SP

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "docs" / "output"

NAVY = RGBColor(0x1F, 0x38, 0x64)
RED = RGBColor(0xB2, 0x18, 0x2B)
GREY = RGBColor(0x55, 0x55, 0x55)


def shade(cell, hex_color):
    el = OxmlElement("w:shd")
    el.set(qn("w:fill"), hex_color)
    cell._tc.get_or_add_tcPr().append(el)


def table(doc, headers, rows, widths=None, head_fill="1F3864"):
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(headers):
        c = t.rows[0].cells[i]
        c.text = ""
        run = c.paragraphs[0].add_run(h)
        run.font.bold = True
        run.font.size = Pt(9.5)
        run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        c.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
        shade(c, head_fill)
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            run = cells[i].paragraphs[0].add_run(str(v))
            run.font.size = Pt(9.5)
            if i == 0:
                run.font.bold = True
    if widths:
        for r in t.rows:
            for i, w in enumerate(widths):
                r.cells[i].width = Cm(w)
    return t


def para(doc, text, size=10.5, color=None, bold=False, space=6, indent=0):
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(space)
    if indent:
        p.paragraph_format.left_indent = Cm(indent)
    run = p.add_run(text)
    run.font.size = Pt(size)
    run.font.bold = bold
    if color is not None:
        run.font.color.rgb = color
    return p


def callout(doc, text):
    """짚고 넘어갈 것 — 표 한 칸으로 만들어 눈에 띄게."""
    t = doc.add_table(rows=1, cols=1)
    t.style = "Table Grid"
    c = t.rows[0].cells[0]
    c.text = ""
    run = c.paragraphs[0].add_run(text)
    run.font.size = Pt(9.5)
    run.font.color.rgb = RED
    shade(c, "FDE9E7")
    doc.add_paragraph().paragraph_format.space_after = Pt(4)
    return t


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    doc = Document()

    st = doc.styles["Normal"]
    st.font.name = "맑은 고딕"
    st._element.rPr.rFonts.set(qn("w:eastAsia"), "맑은 고딕")
    st.font.size = Pt(10.5)
    st.paragraph_format.space_after = Pt(6)
    st.paragraph_format.line_spacing = 1.35
    for lv, size in ((0, 16), (1, 13), (2, 11.5)):
        hs = doc.styles[f"Heading {lv + 1}"]
        hs.font.name = "맑은 고딕"
        hs._element.rPr.rFonts.set(qn("w:eastAsia"), "맑은 고딕")
        hs.font.size = Pt(size)
        hs.font.color.rgb = NAVY

    for s in doc.sections:
        s.left_margin = s.right_margin = Cm(2.0)
        s.top_margin = s.bottom_margin = Cm(1.8)

    # ── 표지 ─────────────────────────────────────────────
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run("지자기 시험 탐사 현장 야장")
    run.font.size = Pt(22)
    run.font.bold = True
    run.font.color.rgb = NAVY
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run("기입 안내서")
    run.font.size = Pt(15)
    run.font.color.rgb = GREY
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run(f"{dt.date.today():%Y년 %m월}   ·   선점 검토 50지점")
    run.font.size = Pt(10)
    run.font.color.rgb = GREY

    para(doc, "")
    para(doc,
         "이 안내서는 현장 야장 카드를 어떻게 채우는지 설명합니다. 카드는 점마다 "
         "한 장씩 있고, 한 장을 다 채우면 그 지점의 기록이 끝납니다. 카드와 함께 "
         "들고 다니시면 됩니다.")

    callout(doc,
            "이 탐사는 법정 선점이 아니고, 얼마 이하여야 합격이라는 기준도 아직 "
            "없습니다. 오히려 그 기준을 만들 자료를 모으는 것이 이번 일의 목적입니다. "
            "그래서 값이 크게 나오더라도 현장에서 그 자리를 접지 마시고, 나온 대로 "
            "적어 주시는 것이 가장 중요합니다.")

    # ── 1. 하루의 흐름 ───────────────────────────────────
    doc.add_heading("1. 하루의 흐름", level=1)
    para(doc, "한 지점에서 하는 일을 순서대로 적으면 이렇습니다.")
    table(doc,
          ["순서", "무엇을", "카드 구역", "대략 걸리는 시간"],
          [["1", "몸에 지닌 자기성 물품 정리, 차량·전자기기 확인", "0", "5분"],
           ["2", "도착 기록, Kp·예보 확인, 장비 준비", "1", "10분"],
           ["3", "GNSS 로 기준점 좌표 실측", "2", "10분"],
           ["4", "수평 자기구배 — 동·서·남·북 네 방향", "3", "1시간 남짓"],
           ["5", "수직 자기구배 — 중심점에서 높이별", "4", "20분"],
           ["6", "방위표지 시준", "5", "20분"],
           ["7", "사진 촬영", "6", "15분"],
           ["8", "중심점 결정, 소견 작성, 완전성 점검", "7", "10분"]],
          widths=[1.4, 8.2, 2.2, 4.0])
    para(doc, "")
    para(doc,
         f"수평은 한 방향에 {len(SP.H_OFFSETS_M)}개 지점"
         f"({', '.join(str(x) for x in SP.H_OFFSETS_M)} m)을 재고, 방향이 바뀔 "
         f"때마다 중심점으로 돌아와 한 번 더 잽니다. 네 방향이니 "
         f"{len(SP.H_DIRECTIONS) * (len(SP.H_OFFSETS_M) + 2)}번, 수직까지 하면 "
         f"한 지점에서 {len(SP.H_DIRECTIONS) * (len(SP.H_OFFSETS_M) + 2) + len(SP.V_HEIGHTS_CM) + 2}"
         f"번쯤 재게 됩니다. 한 자리에서 한 번씩만 재시면 됩니다 — 기기가 안에서 "
         f"여러 번 재어 평균을 내 줍니다.")

    # ── 2. 절대 빠뜨리면 안 되는 것 ──────────────────────
    doc.add_heading("2. 이것만은 꼭 — 빠지면 되살릴 수 없는 것", level=1)
    para(doc,
         "아래 항목들은 나중에 사무실에서 어떤 방법으로도 복원할 수 없습니다. "
         "실제로 예전 야장 68개를 열어 보니 총자력 측정시각이 한 건도 적혀 있지 "
         "않았고, 그 때문에 시간에 따른 보정을 아예 할 수 없게 된 일이 있었습니다.")
    table(doc,
          ["항목", "어디에", "왜 없으면 안 되는가"],
          [["측정 시각 (시·분·초)", "3·4 구역",
            "이 시각으로 재는 동안 자기장이 흘러간 만큼을 빼냅니다. "
            "분까지만 적으면 계산이 되지 않습니다."],
           ["P0 앞·뒤 두 번", "3 구역",
            "한 방향을 다 재고 중심점으로 돌아와 한 번 더 재야 그 사이를 "
            "이을 수 있습니다. 한쪽이 비면 그 방향 전체가 계산되지 않습니다."],
           ["차수", "3·4 구역",
            "다시 잰 것인지 처음 잰 것인지를 가르는 유일한 표시입니다. "
            "안 적으면 1차와 2차가 섞여 둘 다 못 쓰게 됩니다."],
           ["장비 원시파일명·레코드 범위", "1 구역",
            "종이 기록과 장비에 저장된 원본을 잇는 유일한 끈입니다."],
           ["실제 센서 높이", "4 구역",
            "목표 높이가 아니라 센서 가운데가 실제로 몇 cm 였는지를 적어야 "
            "나중에 같은 조건으로 다시 재어 볼 수 있습니다."],
           ["표지 좌표와 취득방법", "5 구역",
            "폐합차가 작아도 참방위각이 틀릴 수 있습니다. 같은 점을 다시 "
            "찾아가 재었을 때 편각이 33.7분이나 어긋난 일이 있었는데, "
            "그 원인이 방문할 때마다 달라진 참방위각이었습니다."]],
          widths=[3.6, 2.0, 10.2])

    doc.add_page_break()

    # ── 3. 구역별 안내 ───────────────────────────────────
    doc.add_heading("3. 구역별로 무엇을 적는가", level=1)

    doc.add_heading("0 구역 — 측정 전 확인", level=2)
    para(doc,
         "측정을 시작하기 전에 몸에 지닌 자기성 물품부터 걷어내 주세요. 시계나 펜, "
         "금속 혁대 하나만 있어도 네 방향 값이 모두 같은 쪽으로 밀립니다. 그러면 "
         "나중에 보았을 때 마치 그 자리의 지형이 만든 값인 것처럼 보이게 되어, "
         "멀쩡한 자리를 나쁘게 판단하게 됩니다. 차량은 측정 범위 밖에 두시고, "
         "무전기나 노트북 같은 전자기기는 꺼 주세요.")
    para(doc,
         "사정상 빼지 못한 것이 있다면 무엇을 왜 못 뺐는지 적어 주시면 됩니다. "
         "적어 두시면 나중에 그 방향 값을 해석할 때 참고할 수 있습니다.",
         color=GREY, size=10)

    doc.add_heading("1 구역 — 도착 · 관측 조건 · 원시자료 연결", level=2)
    para(doc,
         "관측일자와 도착 시각, 관측자, 기상을 적습니다. Kp 지수와 우주기상 예보도 "
         "확인해서 적어 주시는데, 값이 크게 나왔다고 해서 그 자료를 빼지는 마세요. "
         "몇 이상이면 걸러야 하는지가 아직 정해지지 않았고, 사실 그 기준을 만들려고 "
         "이 탐사를 하는 것이기 때문입니다.")
    para(doc,
         "야간 관측은 하지 않습니다. 다만 눈에 띄게 교란이 있었다면 그 구간은 다시 "
         "재시고, 어떤 상황이었는지 소견에 남겨 주세요.")
    callout(doc,
            "장비 원시파일명과 레코드 범위를 꼭 적어 주세요. 장비에 저장된 원본은 "
            "나중에 따로 불러올 예정인데, 그때 어느 파일의 어느 구간이 이 카드에 "
            "해당하는지 알아낼 방법은 여기 적힌 정보밖에 없습니다.")

    doc.add_heading("2 구역 — 기준점 좌표", level=2)
    para(doc,
         "GNSS 로 실측한 위도·경도를 십진도로 적습니다. 미리 알려 드린 사전좌표와 "
         "얼마나 차이가 나는지는 카드가 자동으로 계산해 줍니다. 많이 차이가 난다면 "
         "어떤 사정으로 자리를 옮기셨는지 소견에 적어 주시면 됩니다. 도상에서 고른 "
         "자리라 실제로 가 보면 진입이 어렵거나 나무가 우거진 경우가 있습니다.")

    doc.add_heading("3 구역 — 수평 자기구배", level=2)
    para(doc,
         "가장 시간이 많이 걸리는 부분입니다. 순서는 이렇습니다.")
    para(doc,
         "① 중심점(P0)에서 한 번 재고 시각과 값을 적습니다.  "
         f"② 동쪽으로 {SP.H_OFFSETS_M[0]} m, {SP.H_OFFSETS_M[1]} m, 그다음부터는 "
         f"1 m 씩 {SP.H_OFFSETS_M[-1]} m 까지 옮겨 가며 잽니다.  "
         "③ 다시 중심점으로 돌아와 한 번 더 잽니다.  "
         "④ 서·남·북 방향도 같은 방식으로 반복합니다.",
         indent=0.5)
    para(doc,
         "중심점을 앞뒤로 두 번 재는 이유는, 그 사이에 자기장이 저절로 변한 만큼을 "
         "빼내기 위해서입니다. 자력계가 한 대뿐이라 위치를 옮기는 동안 시간도 함께 "
         "흘러가는데, 앞뒤 두 값을 직선으로 이으면 각 측점을 잰 순간의 중심점 값을 "
         "추정할 수 있습니다. 이것을 빼지 않으면 시간이 흘러간 몫까지 그 자리의 "
         "값인 것처럼 계산됩니다.")
    para(doc,
         "카드가 방향마다 두 가지를 자동으로 보여 줍니다. 「P0 전후차」가 수십 nT 로 "
         "벌어졌다면 그 방향을 재는 동안 자기장이 꽤 흔들렸다는 뜻이니, 여유가 되면 "
         "차수를 올려 다시 재 주세요. 「관측값 max−min」은 그 방향에서 값이 얼마나 "
         "벌어졌는지를 보여 줍니다. 두 값 모두 시간변화를 빼기 전의 참고치라서, "
         "이것만 보고 좋고 나쁨을 가르지는 않습니다.")

    doc.add_heading("4 구역 — 수직 자기구배", level=2)
    para(doc,
         f"중심점에서 센서 높이만 바꿔 가며 잽니다. 지상 {SP.V_HEIGHTS_CM[0]} cm 부터 "
         f"{SP.V_HEIGHTS_CM[-1]} cm 까지 "
         f"{SP.V_HEIGHTS_CM[1] - SP.V_HEIGHTS_CM[0]} cm 간격이고, 시작할 때와 끝날 "
         f"때 기준높이에서 한 번씩 더 재 주시면 됩니다. 수평과 같은 이유입니다.")
    para(doc,
         "높이는 목표값이 아니라 센서 가운데가 실제로 몇 cm 였는지를 적어 주세요. "
         "그리고 어디를 0 으로 보았는지(지면인지 표석 윗면인지)도 함께 적어 주셔야 "
         "나중에 같은 조건을 재현할 수 있습니다.")

    doc.add_heading("5 구역 — 방위표지 시준", level=2)
    para(doc,
         "기준점에서 방위표지를 바라본 방향을 잽니다. 진북을 기준으로 북쪽이 0도, "
         "시계 방향입니다. 방향은 언제나 기준점에서 표지 쪽이며, 반대로 적으면 그 "
         "방문의 편각이 통째로 틀어집니다.")
    para(doc,
         "각도는 곤(gon)으로 적습니다. 한 바퀴가 400 이고 반대쪽이 200 입니다. "
         "정·반 시준의 폐합차는 카드가 계산해 주는데, 0.02 gon(약 65초)을 넘으면 "
         "다시 시준해 주시는 편이 좋습니다.")
    callout(doc,
            "폐합차가 작다고 해서 참방위각이 맞다는 뜻은 아닙니다. 그래서 표지의 "
            "좌표와 그 좌표를 어떻게 얻었는지, 예전에 쓰던 값이 있으면 그것도 함께 "
            "적어 두시는 것입니다. 표지를 바꾸셨다면 반드시 그 사실과 이유를 "
            "남겨 주세요.")

    doc.add_heading("6 구역 — 현장 사진", level=2)
    para(doc,
         "아홉 칸이 있습니다. 중심점 전경, 측정 모습, 방위표지 1·2, 동·서·남·북 "
         "네 방향 전경, 그리고 자기교란 요소입니다. 방위별 전경은 중심점에 서서 그 "
         "방향을 보고, 측선과 10 m 끝점이 함께 담기도록 찍어 주세요. 나중에 어느 "
         "방향 값이 크게 나왔을 때 그 사진에서 원인을 찾게 됩니다.")
    para(doc,
         "사진마다 원본 파일명과 촬영시각을 적어 주세요. 촬영시각까지 적는 이유는, "
         "카메라가 여럿이면 파일명이 겹치거나 옮기는 과정에서 이름이 바뀌는 일이 "
         "있어 파일명만으로는 어느 사진인지 되짚기 어렵기 때문입니다.")

    doc.add_heading("7 구역 — 중심점 결정과 소견", level=2)
    para(doc,
         "어느 한 방향의 값이 유난히 크게 나온다면, 후보지 안에서 더 조용한 자리를 "
         "찾아 중심점을 옮기셔도 됩니다. 다만 옮기셨다면 옮긴 자리에서 3·4 구역을 "
         "다시 재 주셔야 합니다. 옮기기 전 자리에서 잰 값이 새 자리가 조용하다는 "
         "근거가 되지는 못하기 때문입니다.")
    para(doc,
         "옮기실 때는 이 카드를 고쳐 쓰지 마시고, 새 Location ID 로 카드를 한 장 더 "
         "쓰시면 됩니다. 앞 카드도 그대로 두시고요.")
    para(doc,
         "마지막 소견 칸에는 현장에서 느낀 것을 자유롭게 적어 주세요. 값만으로는 "
         "알 수 없는 것들 — 주변에 무엇이 있었는지, 접근이 어땠는지, 다시 온다면 "
         "무엇을 준비해야 하는지 — 이 나중에 큰 도움이 됩니다.")

    doc.add_page_break()

    # ── 4. 자주 생기는 상황 ──────────────────────────────
    doc.add_heading("4. 이럴 때는 어떻게 하나", level=1)
    table(doc,
          ["상황", "이렇게 하시면 됩니다"],
          [["재는 도중 값이 이상하다",
            "그 방향을 처음부터 다시 재시고 차수를 2 로 올려 주세요. 앞서 적은 "
            "값은 지우지 마시고 그대로 두시면 됩니다. 사무실에서 차수로 갈라 "
            "보기 때문에 1차와 2차가 섞이지 않습니다."],
           ["P0 전후차가 크게 벌어졌다",
            "재는 동안 자기장이 흔들렸다는 뜻입니다. 여유가 되면 차수를 올려 "
            "다시 재시고, 어려우면 소견에 적어 주세요."],
           ["중간에 중단했다",
            "중단한 시각과 사유를 적어 두시고, 다시 시작하실 때 차수를 올려 "
            "주세요."],
           ["값이 참고값보다 크게 나온다",
            f"그대로 적어 주세요. 반경 10 m 안에서 {SP.IAGA_RANGE_NT:.0f} nT, "
            f"구배 {SP.EURO_GRAD_NT_PER_M:.0f} nT/m 라는 국제 권고가 있긴 하지만 "
            "국내 기준은 아직 없습니다. 현장에서 그 자리를 접지 마세요."],
           ["방위표지를 바꿔야 했다",
            "바꾸신 것 자체는 문제가 되지 않습니다. 다만 새 표지의 좌표와 "
            "예전 값, 바꾼 이유를 반드시 적어 주세요."],
           ["기준점이 사전좌표와 많이 다르다",
            "실측한 좌표를 적으시고 왜 옮기셨는지 소견에 남겨 주세요. 도상에서 "
            "고른 자리라 진입이 어려운 경우가 있습니다."],
           ["Variometer 를 못 돌렸다",
            "「미운용」으로 표시해 주세요. 값을 못 쓰는 것은 아니지만, 시간에 "
            "따른 변화를 검증할 수 없다는 표시가 남습니다."]],
          widths=[4.2, 11.6])

    # ── 5. 철수 전 점검 ──────────────────────────────────
    doc.add_heading("5. 철수하기 전에", level=1)
    para(doc,
         "카드 3 구역 아래에 「완전성 점검」 줄이 있습니다. 다섯 가지를 스스로 "
         "확인하시고 「예 / 아니오」를 적어 주세요. 하나라도 「아니오」가 있으면 "
         "현장을 떠나기 전에 채우시는 편이 좋습니다. 돌아온 뒤에는 채울 방법이 "
         "없습니다.")
    table(doc,
          ["점검 항목", "확인할 것"],
          [["필수행 다 채움", "빈 칸으로 남은 측점이 없는지"],
           ["시각 시:분:초", "분까지만 적은 곳이 없는지"],
           ["P0 전·후 있음", "네 방향 모두 앞뒤로 중심점을 쟀는지"],
           ["측정시각 P0 사이", "측점 시각이 P0 두 시각 사이에 들어오는지"],
           ["Variometer 연결", "파일명을 적었는지 (미운용이면 그렇게 표시)"]],
          widths=[4.6, 11.2])
    para(doc, "")
    para(doc,
         "그리고 사진 아홉 칸의 파일명과 촬영시각, 장비 원시파일명과 레코드 범위를 "
         "적으셨는지 한 번 더 봐 주세요. 이 둘이 빠지면 나중에 자료를 이을 수 "
         "없습니다.")

    # ── 6. 발주자 결정 ───────────────────────────────────
    doc.add_heading("6. 참고 — 이번 탐사의 성격", level=1)
    para(doc,
         "야장을 채우다 보면 「이건 왜 이렇게 하지」 싶은 곳이 있을 텐데, 대개 "
         "아래 네 가지 결정에서 나옵니다.")
    # ⚠️ 「정해진 것」은 SP.DECISIONS 에서 그대로 가져와 갈라지지 않게 하고,
    #    「현장에서는」 칸만 이 안내서의 말투로 따로 쓴다. 원문을 마침표로 잘라
    #    쓰면 「2.6 h」 같은 숫자에서 끊겨 버린다(실제로 그랬다).
    FIELD_MEANING = {
        "법정 성격": "값이 크게 나와도 현장에서 그 자리를 접지 마시고 그대로 "
                   "적어 주세요. 합격선이 아직 없습니다.",
        "관측 조건": "밤에 나가실 필요는 없습니다. Kp 는 적어 두기만 하시고 "
                   "그 값으로 자료를 빼지는 마세요.",
        "현장 원본 매체": "이 카드가 원본입니다. 장비 원시파일명과 레코드 범위를 "
                       "적어 두셔야 나중에 기기 자료와 이을 수 있습니다.",
        "작업량": "한 지점에 몇 시간이 걸리는지는 이번 탐사로 재어 보는 중입니다. "
                "서두르시다 기록을 빠뜨리는 편이 더 손해입니다.",
    }
    table(doc,
          ["항목", "정해진 것", "현장에서는"],
          [[k, v.replace("**", ""), FIELD_MEANING.get(k, "")]
           for k, v, _ in SP.DECISIONS],
          widths=[3.0, 5.4, 7.4])

    para(doc, "")
    para(doc,
         "측정 설계와 참고값은 IAGA 「Guide for Magnetic Repeat Station Surveys」"
         "(1996) 제4.2절과 유럽 반복관측점 권고를 따랐습니다. 자세한 근거는 "
         "중간보고 자료를 보시면 됩니다.",
         color=GREY, size=9.5)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"{dt.datetime.now():%Y%m%d_%H%M%S}_시험탐사_현장야장_기입안내.docx"
    doc.save(path)
    print(f"[저장] {path}  ({path.stat().st_size/1000:.0f} KB)")
    return path


if __name__ == "__main__":
    main()
