# -*- coding: utf-8 -*-
r"""
시험 탐사 현장 야장 카드 — 2판 (부안 시범 측정 반영)
=====================================================

    python make_trial_field_card_v2.py --into <NOAA Kp 야장.xlsx>            # S01 한 장만
    python make_trial_field_card_v2.py --into <NOAA Kp 야장.xlsx> --all      # 50장 전부

`20260917_시험탐사_현장야장(초안)_수정.xlsx`(부안_점검결과)에 적힌 항목대로 카드를
다시 짠다. 그 시트는 부안에서 실제로 잰 값을 새 양식에 채워 둔 것이라, 계산 규칙을
**그 숫자에서 역산해 맞췄다**(아래 표).

## 무엇이 바뀌었나

| 구역 | 1판 | 2판 |
|---|---|---|
| 기본 정보 | 도착시각·Variometer·원시파일명·GSM 설정 8칸 | **관측일자·관측자·기상·Kp(자동)·우주기상 예보(붙여넣기)·기기번호** |
| 수평 구배 | 동·서 / 남·북 두 표 · 시각·F·품질 · 0.5 m 포함 | **네 방향 한 표 · 시각·평균 F** · 1~10 m · 안 잰 거리는 「-」 |
| 수평 결과 | P0 전후차 · max−min | **P0 전후차 · 관측구간 F 변화폭 · 최대 구간 구배와 발생구간 · P0→10 m ΔF · 평균구배** (자동) |
| 수직 구배 | 20~200 cm 20 cm 간격 12줄 | **막대 마디 4단(185·140·95·50 cm)** · 인접 높이 차·F 차이·구배 자동 |
| 사진 | 9칸(파일명·촬영시각) | **3×3 넓은 칸(한 칸 = 세 열 폭)** — 전경(방향 표시)·중심점·표지 설치 / 동·서·남 / 북·교란 요소·추가 · 칸마다 촬영 방향·파일명 · 한 쪽을 통째로 쓴다 |
| 소견 | 조정·최종좌표·표지설치·인상 | **현황조사 위치(자동)·조정 여부·사유·특이사항·사진 2** + 소견 |
| 끝 | — | **유의사항 4 · 조치사항 5** |
| 뺀 것 | 0구역(측정 전 확인) · 완전성 점검 · GSM 설정·원시파일·Variometer 칸 | 수정본에 없어서 뺐다 — 되살릴지는 발주자 확인 |

2판에 없는 기준점 좌표(GNSS)·방위표지 시준 구역은 **수정본이 다루지 않았으므로 1판
그대로** 둔다.

## 계산 규칙 — 부안 수정본의 숫자로 역산

| 칸 | 정의 | 부안 재현 |
|---|---|---|
| P0 전후차 | P0(후) − P0(전) | 동 0.28 (수정본 0.27 — 원값 소수 셋째 자리 차) |
| 관측구간 F 변화폭 | P0(전)~10 m 가운데 잰 값의 max − min | 서 140.19 ✔ |
| 최대 구간 수평구배 | **바로 앞에 잰 측점**과의 \|ΔF\| ÷ 거리차 중 최대. P0(전)=0 m | 동 12.09 · 6~8 m ✔ / 서 28.04 · P0~5 m ✔ |
| P0→10 m ΔF · 평균구배 | F(10 m) − F(P0 전) · 그 절댓값 ÷ 10 | 동 −38.51 · 3.85 ✔ |
| 수직 F 차이 | (그 높이) − (바로 위 높이) · + 면 땅에 가까울수록 크다 | 140 cm +0.788 ✔ |
| 수직구배 | F 차이 ÷ 인접 높이 차 | 1.75 · 6.57 · 9.85 ✔ · 전체 +8.18 nT / 6.06 ✔ |

⚠️ **시간변화 보정 전 값이다.** 수정본이 보정 없이 계산했으므로 그대로 따랐다.
P0 전후차가 곧 보정 크기의 눈금이다(부안은 0.3 nT 이하).

⚠️ 구간 구배는 **숨긴 보조열(K~AE)** 이 계산한다. 안 잰 거리(「-」·빈칸)를 건너뛰고
직전에 잰 측점과 잇는다 — 부안은 동·남을 2 m 간격, 서·북을 5 m 간격으로 쟀다.

## ⚠️ 사용자 파일은 openpyxl 로 저장하지 않는다

대상 파일(`…_NOAA_Kp예보_자동그래프.xlsx`)에는 **NOAA 웹 쿼리 두 개**(connections ·
queryTables, 열 때 자동 새로 고침)와 **그래프 두 개**가 있다. openpyxl 로 열어 저장하면
웹 쿼리가 사라진다. 그래서 카드만 임시 통합문서에 만들고 **엑셀(COM)이 시트를 복사해
끼워 넣는다.** Kp 수식은 다른 시트를 가리키므로 복사 «뒤에» 엑셀로 넣는다(먼저 넣으면
임시 파일을 가리키는 외부 링크가 된다).
"""
from __future__ import annotations

import argparse
import datetime as dt
import shutil
import sys
import tempfile
from pathlib import Path

from openpyxl import Workbook
from openpyxl.comments import Comment
from openpyxl.worksheet.datavalidation import DataValidation
from openpyxl.formatting.rule import FormulaRule
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.pagebreak import Break

import make_trial_field_card as MC
import trial_survey_points as TP
import trial_survey_spec as SP

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "docs" / "output"
KP_SHEET = "⑧a Kp 예보"

# ── Kp 두 칸 — 카드마다 받아 두는 웹 쿼리 ─────────────────────────────
# ⚠️ 「⑧a Kp 예보」는 열 때마다 새로 받으므로 거기를 가리키면 카드 값이 여는 날마다
#    바뀐다(3일 → 27일 일 최대 → 범위 밖). 그래서 카드가 «자기 사본»을 받아 둔다.
NOAA_3DAY = "https://services.swpc.noaa.gov/text/3-day-forecast.txt"
NOAA_27DAY = "https://services.swpc.noaa.gov/text/27-day-outlook.txt"
GFZ_JSON = "https://kp.gfz.de/app/json/"   # ⚠️ kp.gfz-potsdam.de 는 리디렉션에서 엑셀이 멈췄다
MONTHS = '{"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"}'
SLOTS = ["00-03UT", "03-06UT", "06-09UT", "09-12UT", "12-15UT", "15-18UT", "18-21UT", "21-00UT"]
# 숨긴 열: AG 매개변수 · AH 3일 예보 원문 · AI 27일 예보 원문 · AJ GFZ 원문 · AK~AN 풀이
HID = dict(par="AG", r3="AH", r27="AI", gfz="AJ", p3="AK", g="AL", gk="AM", gs="AN")
DATE_MEMO = (
    "관측일자 입력 양식\n"
    "· 2026-09-14 처럼 「연-월-일」로 적습니다(2026/9/14 도 됩니다).\n"
    "· 2026.09.14 · 9월14일 · 0914 는 날짜로 인식되지 않습니다.\n"
    "· 날짜만 적고 시각은 넣지 않습니다.\n"
    "· 날짜를 넣는 순간 옆 칸 Kp 예보를 NOAA 에서 받아 고정합니다(인터넷 필요).\n"
    "· 예보를 다시 받으려면 날짜를 지웠다가 다시 넣으세요.\n"
    "· Kp(확정)은 관측일이 지난 뒤 「데이터 > 모두 새로 고침」으로 들어옵니다.")


COLS = "ABCDEFGHI"
WIDTHS = [11, 10, 13, 10, 13, 10, 13, 10, 13]
DIRS = ["동", "서", "남", "북"]
DISTS = list(range(1, 11))          # ⚠️ 수정본은 0.5 m 를 뺐다 (규격 H_OFFSETS_M 에는 있다)
VERT = [(4, 185), (3, 140), (2, 95), (1, 50)]   # (막대 마디, 센서 높이 cm) — 수정본 기본값
F_FMT = "#,##0.00"
T_FMT = "h:mm"
RED = Font(name=MC.FONT, size=10, bold=True, color="C00000")
F_HELP = Font(name=MC.FONT, size=8, color="7F7F7F")
ROW_IN = 20          # 손으로 적는 줄의 높이(pt)
PHOTO_ROWS, PHOTO_ROW_H = 12, 16     # 사진 칸 높이 = 12 × 16 pt (85% 에서 세 줄이 한 쪽에 들어가게)

CAUTIONS = [
    "센서는 작업자 및 콘솔 등 자성체와 충분히 이격하여 측정",
    "센서 방향은 측정 전 과정에서 동일하게 유지",
    "센서는 항상 동일한 자세로 사용",
    "철제 구조물, 차량, 휴대품 등 자기장 영향 물체 접근 최소화",
]
ACTIONS = [
    "주변 철제물·매설물·케이블 등 자기장 영향원 제거",
    "동일 지점 재측정으로 재현성 확인",
    "이상구간 전후를 더 촘촘한 간격으로 추가 측정",
    "주변 철제물·매설물·케이블 등 자기장 영향원 확인",
    "이상이 지속되면 위치·변화량 기록 후 중심점 조정 또는 후보지 재검토",
]


def _note(ws, r, text, font=None, fill=None, per=92):
    ws.merge_cells(f"A{r}:I{r}")
    c = ws[f"A{r}"]
    c.value, c.font, c.alignment = text, font or MC.F_SM, MC.AL_TL
    if fill is not None:
        c.fill = fill
    ws.row_dimensions[r].height = MC._text_h(text, per=per)


def _sub(ws, r, text):
    """구역 안의 작은 머리(① ②)."""
    ws.merge_cells(f"A{r}:I{r}")
    c = ws[f"A{r}"]
    c.value, c.font, c.alignment = text, Font(name=MC.FONT, size=10, bold=True,
                                              color="1F3864"), MC.AL_L
    ws.row_dimensions[r].height = 17


# ══════════════════════════════════════════════════════════════
def build_card(wb, p, idx):
    """카드 한 장. Kp 수식이 들어갈 칸 주소를 돌려준다(엑셀이 나중에 채운다)."""
    ws = wb.create_sheet(f"{idx:02d}_{MC._safe(p['지점명'])}")
    for col, w in zip(COLS, WIDTHS):
        ws.column_dimensions[col].width = w
    # ⚠️ 가로로 두면 한 방향 표와 결과표가 쪽 사이에서 끊긴다 — 수정본처럼 세로로.
    ws.page_setup.orientation = "portrait"
    ws.page_setup.paperSize = ws.PAPERSIZE_A4
    ws.page_margins.left = ws.page_margins.right = 0.4
    ws.page_margins.top = ws.page_margins.bottom = 0.55
    ws.page_margins.header = ws.page_margins.footer = 0.25
    # ⚠️ 「폭에 맞춤」은 엑셀이 파일마다 배율을 조금씩 다르게 잡아(0.85~0.87) 꽉 찬 쪽이
    #    경계에서 넘쳤다. 배율을 고정한다 — 쪽 구성은 ①기본·수평 ②수직·좌표·방위표지
    #    ③사진 ④소견.
    ws.page_setup.scale = 85
    ws.print_title_rows = "1:3"
    ws.oddHeader.right.text = f"{p['_sid']} {p['지점명']}"
    ws.oddHeader.right.size = 9
    ws.oddFooter.center.text = "&P / &N"
    ws.oddFooter.center.size = 9
    sid = p["_sid"]

    # ── 머리 (1판 그대로) ────────────────────────────────────
    ws.merge_cells("A1:I1")
    c = ws["A1"]
    c.value = f"[{sid}] {p['지점명']} — 지자기 시험 탐사 현장야장"
    c.font, c.fill, c.alignment = MC.F_TITLE, MC.FILL_TITLE, MC.AL_L
    ws.row_dimensions[1].height = 28
    ws.merge_cells("A2:I2")
    c = ws["A2"]
    caution = SP.SITE_CAUTION.get(p["지점명"], "")
    c.value = (f"구분 {p['구분']} · 원번호 {p['번호']} · 도엽 {p['도엽번호']} "
               f"{p['도엽명']} · 사전좌표 {p['위도']:.6f} / {p['경도']:.6f} · "
               f"표고 {p['표고'] if p['표고'] else '-'} m"
               + (f"   ⚠ {caution}" if caution else ""))
    c.font, c.fill, c.alignment = MC.F_SUB, MC.FILL_TITLE, MC.AL_L
    ws.row_dimensions[2].height = 16
    ws.merge_cells("A3:E3")
    ws["A3"].value, ws["A3"].font, ws["A3"].alignment = (
        f"소재지: {p['소재지']}", MC.F_SM, MC.AL_L)
    MC._lbl(ws, "F3", "Visit ID")
    MC._fld(ws, "G3")
    MC._lbl(ws, "H3", "Location ID")
    MC._fld(ws, "I3")
    ws["I3"].value, ws["I3"].font, ws["I3"].fill = f"{sid}-P0a", MC.F_CALC, MC.FILL_CALC
    r = 5

    # ── 1. 기본 정보 ─────────────────────────────────────────
    r = MC._s(ws, r, "1. 기본 정보")
    for ref, t in (("A", "관측일자"), ("B", "관측자"), ("C", "기상"),
                   ("D", "Kp 예보\n(참고·자동,\n기록 아님)"),
                   ("E{r}:G{r}", "Kp (확정)\n관측일에 실제로 기록된 값 · 자동"),
                   ("H{r}:I{r}", "GSM-19T 기기번호")):
        MC._lbl(ws, ref.format(r=r) if "{" in ref else f"{ref}{r}", t)
    ws.row_dimensions[r].height = 46
    r += 1
    info = r
    MC._fld(ws, f"A{r}", "yyyy-mm-dd")
    memo = Comment(DATE_MEMO, "현장야장")
    memo.width, memo.height = 330, 170
    ws[f"A{r}"].comment = memo
    dv = DataValidation(type="date", operator="between", formula1="DATE(2026,1,1)",
                        formula2="DATE(2027,12,31)", allow_blank=True,
                        showInputMessage=True, showErrorMessage=True,
                        promptTitle="관측일자", prompt="예: 2026-09-14 (연-월-일, 시각 없이)",
                        errorTitle="날짜 양식", error="2026-09-14 처럼 연-월-일로 적어 주세요. "
                        "2026.09.14 는 날짜로 인식되지 않습니다.")
    ws.add_data_validation(dv)
    dv.add(f"A{r}")
    MC._fld(ws, f"B{r}")
    MC._fld(ws, f"C{r}")
    MC._dv(ws, f"C{r}", ["맑음", "구름 조금", "흐림", "비", "눈", "바람 강함"])
    MC._fld(ws, f"D{r}", "0.00", calc="(자동)")
    ws[f"D{r}"].font = Font(name=MC.FONT, size=10, bold=True, color="7F7F7F")
    MC._fld(ws, f"E{r}:G{r}", "0.00", calc="(자동)")
    MC._fld(ws, f"H{r}:I{r}")
    ws.row_dimensions[r].height = 36
    r += 1
    MC._lbl(ws, f"A{r}", "Kp 예보 기준")
    MC._fld(ws, f"B{r}:I{r}", calc="(자동)")
    pre_src = r
    r += 1
    MC._lbl(ws, f"A{r}", "Kp 확정 기준")
    MC._fld(ws, f"B{r}:I{r}", calc="(자동)")
    fix_src = r
    for rr in (pre_src, fix_src):
        ws[f"B{rr}"].font = Font(name=MC.FONT, size=9, color="1F3864")
        ws[f"B{rr}"].alignment = MC.AL_L
        ws.row_dimensions[rr].height = 28
    r += 1
    _note(ws, r, "「Kp 예보」는 현장에 나가기 전에 보는 참고값입니다. 관측일자를 넣는 순간 NOAA "
                 "예보를 받아 그날(UTC) 하루 최대값을 보여 주고, 관측일자를 바꾸기 전까지는 파일을 "
                 "다시 열거나 새로 고쳐도 바뀌지 않습니다. 「Kp (확정)」은 관측일에 실제로 기록된 "
                 "GFZ 관측값(그날 UTC 하루 최대)이며, 관측일이 지난 뒤 「데이터 > 모두 새로 고침」을 "
                 "하면 들어옵니다. 확정값은 한 달쯤 뒤에 나오므로 그 전에는 「잠정」으로 표시됩니다.")
    r += 2

    # ── 2. 수평 구배 ─────────────────────────────────────────
    r = MC._s(ws, r, "2. 수평 자기구배",
              "중심점(P0)에서 한 번 재고, 한 방향으로 거리를 옮겨 가며 잰 뒤 P0 로 돌아와 "
              "한 번 더 잽니다. 네 방향 모두 같습니다. 한 자리에서 여러 번 읽었다면 그 "
              "평균을 적어 주세요. 재지 않은 거리는 「-」로 두시면 계산에서 알아서 빠집니다. "
              "자력계는 Base 모드로 두고 한 자리에서 자동 반복 판독을 받습니다. 측점을 옮길 "
              "때마다 기록을 끊어 새 블록이 생기게 하고, 거리는 블록 순서로만 남으니 "
              "순서를 건너뛰지 마세요.")
    _sub(ws, r, "① 방향별 관측자료")
    r += 1
    MC._lbl(ws, f"A{r}", "거리")
    for k, d in enumerate(DIRS):
        MC._lbl(ws, f"{COLS[1 + 2 * k]}{r}", f"{d} 시각")
        MC._lbl(ws, f"{COLS[2 + 2 * k]}{r}", f"{d} 평균 F(nT)")
    ws.row_dimensions[r].height = 30
    r += 1
    r0 = r                               # P0(전)
    labels = ["P0(전)"] + [f"{x} m" for x in DISTS] + ["P0(후)"]
    for i, lab in enumerate(labels):
        rr = r0 + i
        isP0 = lab.startswith("P0")
        c = ws[f"A{rr}"]
        c.value, c.font, c.border, c.alignment = lab, MC.F_LBL, MC.BOX, MC.AL_C
        c.fill = MC.FILL_P0 if isP0 else MC.FILL_LBL
        ws.row_dimensions[rr].height = ROW_IN
        for k in range(4):
            MC._fld(ws, f"{COLS[1 + 2 * k]}{rr}", T_FMT)
            MC._fld(ws, f"{COLS[2 + 2 * k]}{rr}", F_FMT)
            if isP0:
                ws[f"{COLS[1 + 2 * k]}{rr}"].fill = MC.FILL_P0
                ws[f"{COLS[2 + 2 * k]}{rr}"].fill = MC.FILL_P0
            ws[f"{COLS[2 + 2 * k]}{rr}"].font = Font(name=MC.FONT, size=10, bold=True)
    r1 = r0 + len(DISTS)                 # 10 m
    r2 = r1 + 1                          # P0(후)
    r = r2 + 1

    # 숨긴 보조열 — K: 거리, 방향마다 [직전 거리, 직전 F, 구간 구배, |구배|, 구간 이름]
    ws.column_dimensions["K"].width = 5
    for i in range(len(DISTS) + 1):
        ws[f"K{r0 + i}"] = i
    helper = {}
    for k in range(4):
        F = COLS[2 + 2 * k]
        base = 12 + 5 * k
        cD, cF, cS, cA, cL = (get_column_letter(base + j) for j in range(5))
        helper[k] = (cA, cL)
        ws[f"{cD}{r0}"] = f'=IF(ISNUMBER({F}{r0}),0,"")'
        ws[f"{cF}{r0}"] = f'=IF(ISNUMBER({F}{r0}),{F}{r0},"")'
        for rr in range(r0 + 1, r1 + 1):
            ws[f"{cD}{rr}"] = f'=IF(ISNUMBER({F}{rr}),K{rr},{cD}{rr-1})'
            ws[f"{cF}{rr}"] = f'=IF(ISNUMBER({F}{rr}),{F}{rr},{cF}{rr-1})'
            ws[f"{cS}{rr}"] = (f'=IF(AND(ISNUMBER({F}{rr}),ISNUMBER({cD}{rr-1})),'
                               f'({F}{rr}-{cF}{rr-1})/(K{rr}-{cD}{rr-1}),"")')
            ws[f"{cA}{rr}"] = f'=IF(ISNUMBER({cS}{rr}),ABS({cS}{rr}),"")'
            ws[f"{cL}{rr}"] = (f'=IF(ISNUMBER({cS}{rr}),IF({cD}{rr-1}=0,"P0",{cD}{rr-1})'
                               f'&"~"&K{rr}&" m","")')
    for j in list(range(11, 12 + 5 * 4)) + list(range(33, 41)):   # K~AE 구배 · AG~AN Kp
        cd = ws.column_dimensions[get_column_letter(j)]
        cd.hidden = True
        cd.outlineLevel = 1
    for rr in range(r0, r1 + 1):
        for j in range(11, 12 + 5 * 4):
            ws.cell(rr, j).font = F_HELP

    _sub(ws, r, "② 방향별 수평구배 결과 (자동 계산)")
    r += 1
    for col, t in zip("ABCDEFG", ["방향", "P0 전후차\n(nT)", "관측구간\nF 변화폭 (nT)",
                                   "최대 구간\n수평구배\n(nT/m)", "발생구간",
                                   "P0→10 m\nΔF (nT)", "P0→10 m\n평균구배\n(nT/m)"]):
        MC._lbl(ws, f"{col}{r}", t)
    MC._lbl(ws, f"H{r}:I{r}", "비고")
    ws.row_dimensions[r].height = 44
    r += 1
    res0 = r
    for k, d in enumerate(DIRS):
        F = COLS[2 + 2 * k]
        cA, cL = helper[k]
        c = ws[f"A{r}"]
        c.value, c.font, c.border, c.alignment, c.fill = d, MC.F_LBL, MC.BOX, MC.AL_C, MC.FILL_LBL
        MC._fld(ws, f"B{r}", "0.00",
                calc=f'=IF(COUNT({F}{r0},{F}{r2})<2,"",{F}{r2}-{F}{r0})')
        MC._fld(ws, f"C{r}", "0.00",
                calc=f'=IF(COUNT({F}{r0}:{F}{r1})=0,"",MAX({F}{r0}:{F}{r1})-MIN({F}{r0}:{F}{r1}))')
        MC._fld(ws, f"D{r}", "0.00",
                calc=f'=IF(COUNT({cA}{r0+1}:{cA}{r1})=0,"",MAX({cA}{r0+1}:{cA}{r1}))')
        MC._fld(ws, f"E{r}",
                calc=f'=IF(D{r}="","",INDEX({cL}{r0+1}:{cL}{r1},MATCH(D{r},{cA}{r0+1}:{cA}{r1},0)))')
        MC._fld(ws, f"F{r}", "0.00",
                calc=f'=IF(COUNT({F}{r0},{F}{r1})<2,"",{F}{r1}-{F}{r0})')
        MC._fld(ws, f"G{r}", "0.00", calc=f'=IF(F{r}="","",ABS(F{r})/10)')
        MC._fld(ws, f"H{r}:I{r}")
        ws.row_dimensions[r].height = ROW_IN
        r += 1
    # 참고값(IAGA 50 nT)을 넘은 변화폭만 붉게 — 판정이 아니라 눈에 띄게 하는 표시
    ws.conditional_formatting.add(
        f"C{res0}:C{r-1}",
        FormulaRule(formula=[f'AND(ISNUMBER(C{res0}),C{res0}>{SP.IAGA_RANGE_NT:.0f})'],
                    font=RED))
    _note(ws, r, f"참고 1) IAGA: 반복관측점 후보지의 총자기장 변화 확인 시 반경 "
                 f"{SP.IAGA_RADIUS_M:.0f} m 내 약 {SP.IAGA_RANGE_NT:.0f} nT 참고 "
                 f"(넘으면 변화폭이 붉게 표시됩니다)")
    r += 1
    _note(ws, r, f"참고 2) 유럽 반복관측 권고: 측점의 자기장 구배 "
                 f"{SP.EURO_GRAD_NT_PER_M:.0f} nT/m 미만 권고, PPM cross-profile 및 "
                 f"수직구배 확인")
    r += 1
    _note(ws, r, "구간 구배는 바로 앞에 잰 측점과의 F 차이를 거리 차로 나눈 값이고, P0(전)을 "
                 "0 m 로 봅니다. 모두 시간변화를 빼기 «전» 값이라, P0 전후차가 크게 벌어졌다면 "
                 "그만큼 섞여 있다는 뜻입니다. 국내 기준은 아직 없으니 값이 커도 현장에서 "
                 "자리를 접지 마시고 나온 대로 적어 주세요.", font=MC.F_WARN, fill=MC.FILL_WARN,
          per=80)
    r += 2

    # ── 3. 수직 구배 ─────────────────────────────────────────
    ws.row_breaks.append(Break(id=r - 1))
    r = MC._s(ws, r, "3. 수직 자기구배 — 중심점(P0)에서 높이별",
              "막대 마디 수만 바꿔 가며 중심점 위에서 잽니다. 높이는 센서 가운데 높이이고 "
              "기본값을 채워 두었으니, 실제로 다르면 고쳐 적어 주세요. 한 높이에서 여러 번 "
              "읽었다면 평균을 적고, 읽은 시각 범위는 비고에 남겨 주세요. F 차이는 «그 높이 "
              "− 바로 위 높이»라서 + 면 땅에 가까울수록 값이 커진다는 뜻입니다.")
    for col, t in zip("ABCDEF", ["봉 번호", "높이\n(cm)", "평균 F\n(nT)", "인접\n높이 차",
                                  "F 차이\n(nT)", "수직구배\n(nT/m)"]):
        MC._lbl(ws, f"{col}{r}", t)
    MC._lbl(ws, f"G{r}:I{r}", "비고 (측정 시각 · 삭제한 값 등)")
    ws.row_dimensions[r].height = 30
    r += 1
    v0 = r
    for i, (bong, h) in enumerate(VERT):
        rr = v0 + i
        ws.row_dimensions[rr].height = ROW_IN
        c = ws[f"A{rr}"]
        c.value, c.font, c.border, c.alignment, c.fill = bong, MC.F_LBL, MC.BOX, MC.AL_C, MC.FILL_LBL
        MC._fld(ws, f"B{rr}", '0" cm"')
        ws[f"B{rr}"].value = h
        MC._fld(ws, f"C{rr}", F_FMT)
        ws[f"C{rr}"].font = Font(name=MC.FONT, size=10, bold=True)
        if i == 0:
            for cc in "DEF":
                MC._lbl(ws, f"{cc}{rr}", "")
        else:
            MC._fld(ws, f"D{rr}", '0.00" m"',
                    calc=f'=IF(COUNT(B{rr},B{rr-1})<2,"",(B{rr-1}-B{rr})/100)')
            MC._fld(ws, f"E{rr}", '+0.000;-0.000;0.000',
                    calc=f'=IF(COUNT(C{rr},C{rr-1})<2,"",C{rr}-C{rr-1})')
            MC._fld(ws, f"F{rr}", '+0.00;-0.00;0.00',
                    calc=f'=IF(OR(D{rr}="",E{rr}="",D{rr}=0),"",E{rr}/D{rr})')
        MC._fld(ws, f"G{rr}:I{rr}")
    vt, vb = v0, v0 + len(VERT) - 1
    r = vb + 1
    MC._lbl(ws, f"A{r}", "전체 높이구간")
    MC._fld(ws, f"B{r}", '0.00" m"',
            calc=f'=IF(COUNT(B{vt},B{vb})<2,"",(B{vt}-B{vb})/100)')
    MC._lbl(ws, f"C{r}", "")
    MC._lbl(ws, f"D{r}", "")
    MC._fld(ws, f"E{r}", '+0.000;-0.000;0.000',
            calc=f'=IF(COUNT(C{vt},C{vb})<2,"",C{vb}-C{vt})')
    MC._fld(ws, f"F{r}", '+0.00;-0.00;0.00',
            calc=f'=IF(OR(B{r}="",E{r}="",B{r}=0),"",E{r}/B{r})')
    MC._fld(ws, f"G{r}:I{r}", calc='="맨 위 ~ 맨 아래 전체 구간"')
    r += 2

    # ── 4. 기준점 좌표 (1판 그대로) ──────────────────────────
    r = MC._s(ws, r, "4. 기준점 좌표 (GNSS 실측)",
              "십진도로 적어 주세요. 미리 알려 드린 사전좌표와 많이 차이가 난다면, 어떤 "
              "사정으로 자리를 옮기셨는지 소견에 적어 주시면 됩니다.")
    MC._lbl(ws, f"A{r}", "위도");   MC._fld(ws, f"B{r}:C{r}", "0.000000")
    MC._lbl(ws, f"D{r}", "경도");   MC._fld(ws, f"E{r}:F{r}", "0.000000")
    MC._lbl(ws, f"G{r}", "취득방법"); MC._fld(ws, f"H{r}:I{r}")
    MC._dv(ws, f"H{r}", ["RTK-GNSS", "네트워크RTK", "정적GNSS", "휴대GNSS"])
    r += 1
    MC._lbl(ws, f"A{r}", "사전값 대비(m)")
    MC._fld(ws, f"B{r}:C{r}", "0.0",
            calc=f'=IF(COUNT(B{r-1},E{r-1})<2,"",ROUND(SQRT((({p["경도"]}-E{r-1})*'
                 f'111320*COS(RADIANS({p["위도"]})))^2+(({p["위도"]}-B{r-1})*110574)^2),1))')
    r += 2

    # ── 5. 방위표지 (1판 그대로) ─────────────────────────────
    r = MC._s(ws, r, "5. 방위표지 시준 — 대상 지점만 (진북 기준 · 북 0° · 시계방향)",
              "대상 지점이 아니면 「해당없음」으로 두시면 됩니다. 각도는 곤(gon)으로 적어 "
              "주세요(한 바퀴 400, 반대쪽 200). 정·반 시준의 폐합차가 0.02 gon 을 넘으면 "
              "다시 시준해 주시고, 표지 좌표와 그 좌표를 얻은 방법도 함께 적어 주세요.")
    MC._lbl(ws, f"A{r}", "표지");          MC._lbl(ws, f"B{r}", "표지 위도")
    MC._lbl(ws, f"C{r}", "표지 경도");     MC._lbl(ws, f"D{r}", "좌표 취득")
    MC._lbl(ws, f"E{r}", "정 시준(gon)");  MC._lbl(ws, f"F{r}", "반 시준(gon)")
    MC._lbl(ws, f"G{r}", "폐합차(gon)");   MC._lbl(ws, f"H{r}:I{r}", "표지 변경 여부 · 사유")
    r += 1
    for k in (1, 2):
        c = ws[f"A{r}"]
        c.value, c.font, c.fill, c.border, c.alignment = (
            f"방위표지{k}", MC.F_LBL, MC.FILL_LBL, MC.BOX, MC.AL_C)
        MC._fld(ws, f"B{r}", "0.000000")
        MC._fld(ws, f"C{r}", "0.000000")
        MC._fld(ws, f"D{r}")
        MC._dv(ws, f"D{r}", ["RTK-GNSS", "네트워크RTK", "정적GNSS", "휴대GNSS", "기존 성과 인용"])
        MC._fld(ws, f"E{r}", "0.0000")
        MC._fld(ws, f"F{r}", "0.0000")
        MC._fld(ws, f"G{r}", "0.0000",
                calc=f'=IF(COUNT(E{r}:F{r})<2,"",ABS(ABS(F{r}-E{r})-200))')
        MC._fld(ws, f"H{r}:I{r}")
        r += 1
    r += 1

    # ── 6. 현장 사진 기록 ───────────────────────────────────
    # ⚠️ 한 줄에 일곱 칸을 두었더니 사진 칸이 한 열 폭밖에 안 됐다(발주자 지적).
    #    세 칸씩 세 줄 — 한 칸이 세 열 폭이고 한 쪽을 통째로 쓴다.
    ws.row_breaks.append(Break(id=r - 1))
    r = MC._s(ws, r, "6. 현장 사진 기록",
              "전경은 사진 방향이 드러나게 찍어 주세요(나침반 화면을 함께 담으면 좋습니다). "
              "방향 사진은 중심점에 서서 그 방향을 보고 찍습니다. 사진 아래에 촬영 방향과 "
              "파일명을 적어 주세요.")
    SPANS = [("A", "C"), ("D", "F"), ("G", "I")]
    SHOTS = [
        [("전경 (방향 표시)", "[ 촬영 방향이 보이게 · 나침반 화면 함께 ]"),
         ("중심점 사진", "[ 표석·마커가 보이게 ]"),
         ("표지 설치 사진", "[ 설치한 표지 ]")],
        [("동 사진", "[ 중심점에서 동쪽을 보고 ]"),
         ("서 사진", "[ 중심점에서 서쪽을 보고 ]"),
         ("남 사진", "[ 중심점에서 남쪽을 보고 ]")],
        [("북 사진", "[ 중심점에서 북쪽을 보고 ]"),
         ("자기교란 요소 · 기타", "[ 없으면 비워 두세요 ]"),
         ("추가 사진", "[ 필요할 때 ]")],
    ]
    for band in SHOTS:
        for (a, b), (t, _) in zip(SPANS, band):
            MC._lbl(ws, f"{a}{r}:{b}{r}", t)
        ws.row_dimensions[r].height = 20
        r += 1
        for (a, b), (_, hint) in zip(SPANS, band):
            MC._photo(ws, r, (a, b), hint, rows=PHOTO_ROWS)
        for rr in range(r, r + PHOTO_ROWS):
            ws.row_dimensions[rr].height = PHOTO_ROW_H
        r += PHOTO_ROWS
        for lab in ("촬영 방향", "파일명"):
            for (a, b), _ in zip(SPANS, band):
                MC._lbl(ws, f"{a}{r}", lab)
                MC._fld(ws, f"{chr(ord(a) + 1)}{r}:{b}{r}")
            ws.row_dimensions[r].height = 20
            r += 1
    r += 1

    # ── 7. 현장 소견 ────────────────────────────────────────
    ws.row_breaks.append(Break(id=r - 1))
    r = MC._s(ws, r, "7. 현장 소견")
    MC._lbl(ws, f"A{r}", "현황조사 위치")
    MC._lbl(ws, f"B{r}", "위치 조정\n여부")
    MC._lbl(ws, f"C{r}:D{r}", "위치 조정 사유")
    MC._lbl(ws, f"E{r}:G{r}", "관측 특이사항")
    MC._lbl(ws, f"H{r}", "사진")
    MC._lbl(ws, f"I{r}", "사진")
    ws.row_dimensions[r].height = 30
    r += 1
    MC._fld(ws, f"A{r}", calc="=$I$3")
    MC._fld(ws, f"B{r}")
    MC._dv(ws, f"B{r}", ["동일", "조정"])
    MC._fld(ws, f"C{r}:D{r}")
    MC._fld(ws, f"E{r}:G{r}")
    for cc in "CE":
        ws[f"{cc}{r}"].alignment = MC.AL_TL
    for cc in "HI":
        MC._fld(ws, f"{cc}{r}")
        ws[f"{cc}{r}"].value = "[ 사진 ]"
        ws[f"{cc}{r}"].font = Font(name=MC.FONT, size=9, color="999999")
    ws.row_dimensions[r].height = 120     # 7구역은 한 쪽을 혼자 쓰므로 손글씨 칸을 넉넉히
    r += 1
    MC._lbl(ws, f"A{r}", "소견")
    MC._fld(ws, f"B{r}:I{r}")
    ws[f"B{r}"].alignment = MC.AL_TL
    ws.row_dimensions[r].height = 150
    r += 2

    # ── 유의사항 · 조치사항 ─────────────────────────────────
    MC._lbl(ws, f"A{r}:D{r}", "유의사항")
    MC._lbl(ws, f"E{r}:I{r}", "조치사항")
    r += 1
    for i in range(max(len(CAUTIONS), len(ACTIONS))):
        for (n, t0, t1), items in (((("A", "B", "D")), CAUTIONS), ((("E", "F", "I")), ACTIONS)):
            if i < len(items):
                c = ws[f"{n}{r}"]
                c.value, c.font, c.alignment, c.border = i + 1, MC.F_LBL, MC.AL_C, MC.BOX
                ws.merge_cells(f"{t0}{r}:{t1}{r}")
                c = ws[f"{t0}{r}"]
                c.value, c.font, c.alignment = items[i], MC.F_VAL, MC.AL_L
                MC._span(ws, f"{t0}{r}:{t1}{r}", None)
        ws.row_dimensions[r].height = 32      # 두 줄로 넘어가는 문장이 있다
        r += 1
    ws.print_area = f"A1:I{r}"
    return ws.title, {"info": info, "pre_src": pre_src, "fix_src": fix_src, "p0": r0}


def kp_formulas(info, pre_src, fix_src, p0=None):
    """Kp 예보(카드에 받아 둔 NOAA 사본) · Kp 확정(GFZ) 수식."""
    A = f"$A${info}"
    H = HID
    f = {}
    # 매개변수 — 관측일자가 바뀌면 이 값이 바뀌고, 그때만 웹 쿼리가 다시 받는다
    f[f"{H['par']}1"] = f'=IF({A}="","none",TEXT({A},"yyyymmdd"))'
    f[f"{H['par']}2"] = f'=IF({A}="","2026-09-01",TEXT({A},"yyyy-mm-dd"))'
    # 3일 예보 사본 풀이 — :Issued: 줄(2행)에서 발행일, 날짜 열은 발행일 + 0·1·2
    r3 = f"${H['r3']}$1:${H['r3']}$70"
    p = H["p3"]
    f[f"{p}1"] = (f'=IFERROR(DATE(VALUE(MID(${H["r3"]}$2,10,4)),'
                  f'MATCH(MID(${H["r3"]}$2,15,3),{MONTHS},0),VALUE(MID(${H["r3"]}$2,19,2))),"")')
    f[f"{p}2"] = f'=IF(OR({A}="",{p}1=""),"",INT({A})-{p}1)'
    for k, slot in enumerate(SLOTS):
        f[f"{p}{3 + k}"] = (f'=IF(OR({p}2="",N({p}2)<0,N({p}2)>2),"",'
                            f'IFERROR(VALUE(MID(INDEX({r3},MATCH("{slot}*",{r3},0)),15+13*{p}2,4)),""))')
    f[f"{p}11"] = f'=IF(COUNT({p}3:{p}10)=0,"",MAX({p}3:{p}10))'
    # 27일 예보 사본 — 「2026 Sep 14 …」 줄의 42번째 두 글자가 Kp
    r27 = f"${H['r27']}$1:${H['r27']}$60"
    f[f"{p}12"] = f'=IF({A}="","",YEAR({A})&" "&INDEX({MONTHS},MONTH({A}))&" "&TEXT(DAY({A}),"00"))'
    f[f"{p}13"] = f'=IF({p}12="","",IFERROR(VALUE(MID(INDEX({r27},MATCH({p}12&"*",{r27},0)),42,2)),""))'
    # GFZ JSON 풀이 — {"Kp":[…],"datetime":[…],…,"status":[…]}
    g, gk, gs, J = H["g"], H["gk"], H["gs"], f"${H['gfz']}$1"
    f[f"{g}1"] = (f'=IFERROR(MID({J},SEARCH("""Kp"":[",{J})+6,'
                  f'SEARCH("]",{J},SEARCH("""Kp"":[",{J}))-SEARCH("""Kp"":[",{J})-6),"")')
    f[f"{g}2"] = (f'=IFERROR(MID({J},SEARCH("""status"":[",{J})+10,'
                  f'SEARCH("]",{J},SEARCH("""status"":[",{J}))-SEARCH("""status"":[",{J})-10),"")')
    f[f"{g}3"] = f'=AND({A}<>"",ISNUMBER(SEARCH(TEXT({A},"yyyy-mm-dd")&"T00",{J})))'
    for n in range(1, 9):
        f[f"{gk}{n}"] = (f'=IF(OR(NOT(${g}$3),${g}$1=""),"",IFERROR(VALUE(TRIM(MID('
                         f'SUBSTITUTE(${g}$1,",",REPT(" ",99)),{(n - 1) * 99 + 1},99))),""))')
        f[f"{gs}{n}"] = (f'=IF({gk}{n}="","",SUBSTITUTE(TRIM(MID('
                         f'SUBSTITUTE(${g}$2,",",REPT(" ",99)),{(n - 1) * 99 + 1},99)),"""",""))')
    f[f"{g}4"] = f'=COUNT({gk}1:{gk}8)'
    f[f"{g}5"] = f'=COUNTIF({gs}1:{gs}8,"def")'
    f[f"{g}6"] = f'=IF({g}4=0,"",MAX({gk}1:{gk}8))'
    # 보이는 칸
    f[f"D{info}"] = (f'=IF({A}="","",IF(ISNUMBER({p}11),{p}11,'
                     f'IF(ISNUMBER({p}13),{p}13,"예보범위 밖")))')
    f[f"B{pre_src}"] = (
        f'=IF({A}="","관측일자를 넣으면 그 순간의 NOAA 예보로 채워지고, 관측일자를 바꾸기 전까지 그대로 남습니다",'
        f'IF(ISNUMBER({p}11),"NOAA 3일 예보 ("&MID(${H["r3"]}$2,10,20)&" 발행) · "&TEXT({A},"m/d")&" UTC 하루 최대",'
        f'IF(ISNUMBER({p}13),"NOAA 27일 예보 ("&MID(${H["r27"]}$2,10,20)&" 발행) · "&TEXT({A},"m/d")&" UTC 하루 최대",'
        f'"예보 범위 밖이거나 받지 못했습니다 — 인터넷 연결 후 관측일자를 지웠다가 다시 넣어 보세요"))'
        f'&" · 관측일자를 넣을 때 받은 값(이후 고정)")')
    f[f"E{info}"] = f'=IF({A}="","",IF({g}4=0,"미수신",{g}6))'
    f[f"B{fix_src}"] = (
        f'=IF({A}="","관측일이 지난 뒤 「데이터 > 모두 새로 고침」을 하면 GFZ 관측 Kp 가 들어옵니다",'
        f'IF({g}4=0,"아직 없음 — 관측일이 지난 뒤 「데이터 > 모두 새로 고침」",'
        f'"GFZ 관측 Kp · "&TEXT({A},"m/d")&" UTC 하루 최대 · "&'
        f'IF({g}5={g}4,"확정","잠정 ("&{g}5&"/"&{g}4&" 확정 — 한 달쯤 뒤 새로 고침하면 확정값)")))')
    return f


def add_kp_queries(ws, tag):
    """카드 한 장에 웹 쿼리 셋 — NOAA 3일·27일 사본(관측일자 바뀔 때만) · GFZ(모두 새로 고침)."""
    for i in range(ws.QueryTables.Count, 0, -1):
        ws.QueryTables(i).Delete()
    H = HID
    specs = [
        (f"KpPre3_{tag}", f'URL;{NOAA_3DAY}?d=["d","관측일"]', f"{H['r3']}1", [f"{H['par']}1"], False),
        (f"KpPre27_{tag}", f'URL;{NOAA_27DAY}?d=["d","관측일"]', f"{H['r27']}1", [f"{H['par']}1"], False),
        (f"KpFix_{tag}", f'URL;{GFZ_JSON}?start=["s","시작일"]T00:00:00Z'
                         f'&end=["e","종료일"]T23:59:59Z&index=Kp',
         f"{H['gfz']}1", [f"{H['par']}2", f"{H['par']}2"], True),
    ]
    for name, url, dest, cells, with_all in specs:
        qt = ws.QueryTables.Add(Connection=url, Destination=ws.Range(dest))
        qt.Name = name
        qt.WebSelectionType = 1           # 전체 페이지
        qt.WebFormatting = 1              # 서식 없음
        qt.RefreshStyle = 0               # 덮어쓰기
        qt.BackgroundQuery = False
        qt.SaveData = True                # 받은 원문을 파일에 남긴다
        qt.RefreshOnFileOpen = False      # ⚠️ 열 때 받으면 고정이 깨진다
        qt.AdjustColumnWidth = False
        qt.PreserveFormatting = True
        for i, cell in zip(range(1, qt.Parameters.Count + 1), cells):
            prm = qt.Parameters(i)
            prm.SetParam(2, ws.Range(cell))   # xlRange
            prm.RefreshOnChange = True        # 관측일자가 바뀌면 다시 받는다
        qt.Refresh(False)
        # ⚠️ 예보 사본은 「모두 새로 고침」에서 빼야 고정된다. GFZ 는 거기서 받는다.
        qt.WorkbookConnection.RefreshWithRefreshAll = with_all


# ══════════════════════════════════════════════════════════════
def inject(target: Path, out: Path, cards: Path, jobs):
    """엑셀로 카드 시트를 바꿔 끼운다 — 웹 쿼리·그래프를 지키기 위해서다."""
    import pythoncom
    import win32com.client as win32

    shutil.copy(target, out)
    pythoncom.CoInitialize()
    xl = win32.DispatchEx("Excel.Application")
    xl.Visible = False
    xl.DisplayAlerts = False
    xl.AskToUpdateLinks = False
    try:
        dst = xl.Workbooks.Open(str(out), UpdateLinks=0)
        src = xl.Workbooks.Open(str(cards), UpdateLinks=0, ReadOnly=True)
        for name, cells in jobs:
            old = dst.Worksheets(name)
            src.Worksheets(name).Copy(Before=old)
            new = dst.Worksheets(old.Index - 1)
            old.Delete()
            new.Name = name
            for addr, f in cells.items():
                new.Range(addr).Formula = f
            add_kp_queries(new, name.split("_")[0])
            new.Activate()
            new.Range("A1").Select()
        src.Close(False)
        dst.Worksheets(jobs[0][0]).Activate()
        dst.Save()
        dst.Close(True)
    finally:
        xl.Quit()
        pythoncom.CoUninitialize()


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--into", required=True, help="NOAA Kp 예보가 붙은 현장야장 xlsx (원본은 건드리지 않음)")
    ap.add_argument("--sites", default="S01", help="쉼표로 구분한 Site ID (기본 S01)")
    ap.add_argument("--all", action="store_true", help="50장 전부")
    ap.add_argument("--keep-cards", help="임시 카드 통합문서를 이 경로에 남긴다(검증용)")
    a = ap.parse_args()

    pts = TP.load_points()
    for i, p in enumerate(pts, 1):
        p["_sid"] = f"S{i:02d}"
    want = {p["_sid"] for p in pts} if a.all else set(a.sites.split(","))
    wb = Workbook()
    wb.remove(wb.active)
    jobs = []
    for i, p in enumerate(pts, 1):
        if p["_sid"] in want:
            name, rows = build_card(wb, p, i)
            jobs.append((name, kp_formulas(**rows)))
    tmpdir = Path(tempfile.mkdtemp())
    cards = tmpdir / "cards_v2.xlsx"
    wb.save(cards)
    if a.keep_cards:
        shutil.copy(cards, a.keep_cards)

    target = Path(a.into).resolve()
    tag = "전체50" if a.all else "_".join(sorted(want))
    out = OUT_DIR / f"{dt.datetime.now():%Y%m%d_%H%M%S}_시험탐사_현장야장_2판_{tag}.xlsx"
    inject(target, out, cards, jobs)
    print(f"카드 {len(jobs)}장 교체: {', '.join(n for n, _ in jobs)}")
    print(f"[저장] {out}  ({out.stat().st_size / 1e6:.2f} MB)")
    return out


if __name__ == "__main__":
    main()
