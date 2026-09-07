# -*- coding: utf-8 -*-
"""
지자기 시험 탐사 표준 야장 (자기교란 확인표) — `docs/output/*_시험탐사_표준야장.xlsx`
======================================================================================

선점 검토 50점(등급 A 5 · B 도엽대표 29 · 기존 관측망 선점대상 16)의 시험 탐사용
**사무실 정리 야장**. 중간보고 자료 23면의 기록양식 3종 중 셋째
**「자기교란 확인표」**다 — 앞의 둘(도상선점표 103점 · 현황조사표 34점)은 끝났다.
현장 기록은 `make_trial_field_card.py` 의 점별 카드로 받는다.

⚠️ **절대측정(D·I)은 이 단계에 없다**(2026-09-07 사용자 확인). 원 계획 31면의
선점실시 절차는 「GNSS 수신기 → 오버하우저 자기구배 확인 → 중심점 선정 및 표지
설치」이며 DI-flux 절대관측은 선점이 «끝난 뒤»의 별도 단계다. 따라서 §19①·§20 도
이 야장의 게이트가 아니고, 재는 것은 **오버하우저 총자력뿐**이다.

    python make_trial_survey_book.py

## 규격은 내가 정하지 않는다 — `trial_survey_spec.py`

⚠️ **원 계획을 읽지 않고 야장을 설계했다가 통째로 다시 만든 일이 있다**(2026-09-07).
IAGA 측선 규격·수직 구배·정량 판정기준·이동식 Variometer·장비 검교정 주기가 이미
중간보고 자료에 확정돼 있었는데 그것을 모르고 임의 설계했다. 측선 간격도 근거 없이
「±30 m·2 m 간격」으로 잡았었다. 새 항목을 넣기 전에 **반드시 원 자료를 먼저 볼 것.**

## Codex 독립 검토 2회 (`2026_nscenter/docs/codex/*_Geomag_TrialSurvey*`)

1차 Critical 6·Major 6 → 2차 Critical 6·Major 8. 2차의 요지는 **「측정 횟수가
부족한 게 아니라, 보정에 쓸 관측행을 유일하게 연결하는 키와 판정지표의 연산 정의가
완결되지 않았다」**는 것이다. 반영한 것:

| 지적 | 반영 |
|---|---|
| **C1** Set 식별자가 키에 없고 `SUMIFS` 가 누락·중복을 0/합계로 삼킨다 | ⑧⑨ 에 **Set 검증 열** 신설 — 유효 P0 가 전·후 «각각 정확히 1개»이고 `t0 < t < t1` 일 때만 계산한다. 아니면 수치 대신 **오류 사유**를 낸다 |
| **C2** P0 선형성은 두 끝점으로 검증 못 하고 Variometer 조건이 충돌 | Set 지속시간·P0 전후차를 자동 산출. Variometer 파일 ID·시간구간을 세트마다 기록하고, 없으면 **「선형성 미검증」**으로 구분. 임의 시간상한은 두지 않는다 |
| **C3** `ΔF/거리` 와 「그 점의 자기구배」를 같은 것으로 다룸 | **부호를 보존**한 중심점기준 변화율·절댓값·**인접구간 변화율**을 각각 다른 열로. 지점 대표값은 **확정하지 않는다**(국내 적용방안 결정 전) |
| **C4** 4방향 측선으로 「반경 10 m 전체」를 판정 못 함 | 표현을 **「반경 10 m 내 4개 방위 측선의 관측점에서 확인된 변화」**로 한정. 50 nT 연산도 `max|F−P0|` 와 `max−min` 둘 다 두고 «정의 미확정» 표기 |
| **C5** 옮긴 중심점의 적합성 재확인 절차 없음 | 최초 후보 중심과 최종 중심에 **다른 Location ID** · 재측정 여부를 명시 필드로 |
| **C6** D·I 원시관측과 F 판독의 연결이 불완전 | ⚠️ **해당 없음으로 바뀌었다** — 절대측정을 이 단계에서 하지 않으므로 D·I 관측행 자체가 없다(2026-09-07) |
| **M1** 수직 기준높이·산출물 미정의 | 기준높이·실제 센서 중심높이·측정기준면 기록. 자동 판정하지 않고 기술지표만 |
| **M2** 작업량이 «행 수»로 과소 표현 | 행마다 3회 판독이면 점당 **64 회**다(1 판독 기준). ⑭ 에 산정표를 싣는다 |
| **M3** PPM 1대인데 기기대조 시트가 2대를 전제 | 두지점 F 시트 자체를 **뺐다** — 절대측정 절차라 이 단계에 없다 |
| **M4** §13② 시설별 증빙 필드 없음 | 시설 «유형마다 한 행» — ID·출처·기준일·좌표·도상·현장·방법·사진·판정·확인자 |
| **M5** 미원 HOLD 범위 불명확 | 「현 물리점 식별 / 공식 좌표 확인 / 과거점 동일성 / 성분별 대조」 넷으로 분리 |
| **M6** 중단·결측 시 세트 처리 규칙 없음 | Set 상태·중단시각·사유·재시작 Set ID |
| **M7** 「사람이 적는 판정」에 근거·기준 버전 없음 | 관측사실 / 국제 참고기준 대비 / 국내 판정상태 / 사유 / 기준명·버전 / 판단자·일자 분리 |
| **M8** 해소 주장을 입증할 추적표 없음 | ⑭ **검토 반영 추적** 시트 신설 |
| **m1** 예시에 거리가 없어 1.75 nT/m 재현 불가 | 거리 **2 m** 명시 |
| **m2** 상태코드에 사유·이력 없음 | 사유·입력자·입력시각 열 |

⚠️ Codex 는 작업공간에 실제 파일이 없어 셀 단위 구현을 검증하지 못했다. 그래서
⑭ 추적 시트에 **검증 사례와 실제 결과**를 함께 적는다.
"""
from __future__ import annotations

import datetime as dt
import sys
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.datavalidation import DataValidation

import trial_survey_points as TP
import trial_survey_spec as SP

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "docs" / "output"

# ⚠️ HOLD 는 «범위»를 나눠야 한다(Codex M5). 물리점 자체가 특정되지 않으면 전부
#    보류지만, 현 표석이 확실하고 과거 좌표와의 연결만 불명확하면 신규 측정은
#    별도 Location ID 로 진행하고 «과거 성과 대조만» 보류하는 것이 옳다.
HOLD_SITES = {
    "미원": ("과거대조", "'10~'19 성과표와 2019 야장 좌표가 353 m 어긋난다. "
                       "현 표석에서의 신규 측정은 별도 Location ID 로 진행하되, "
                       "과거 성과와의 대조는 표석 동일성 확인 전까지 보류."),
}
STATUS = ["NOT AVAILABLE", "NOT OBSERVED", "INVALID", "NOT APPLICABLE", "PENDING"]
FACILITIES = [("직류철도", 5.0), ("교류·일반철도", 2.0),
              ("고압철탑", 1.0), ("송전탑", 0.5)]

FONT = "맑은 고딕"
F_TITLE = Font(name=FONT, size=13, bold=True, color="FFFFFF")
F_SEC = Font(name=FONT, size=10.5, bold=True, color="FFFFFF")
F_LBL = Font(name=FONT, size=10, bold=True)
F_VAL = Font(name=FONT, size=10)
F_SM = Font(name=FONT, size=9, color="555555")
F_WARN = Font(name=FONT, size=9, bold=True, color="B2182B")
F_HDR = Font(name=FONT, size=9, bold=True, color="FFFFFF")
F_KEY = Font(name=FONT, size=9, bold=True, color="FFFFFF")

FILL_TITLE = PatternFill("solid", fgColor="1F3864")
FILL_KEY = PatternFill("solid", fgColor="2E5C8A")
FILL_SEC = PatternFill("solid", fgColor="4472C4")
FILL_LBL = PatternFill("solid", fgColor="F2F2F2")
FILL_FIELD = PatternFill("solid", fgColor="FFF8E1")
FILL_AUTO = PatternFill("solid", fgColor="E8F0E8")
FILL_WARN = PatternFill("solid", fgColor="FDE9E7")
FILL_LOCK = PatternFill("solid", fgColor="EDEDED")
FILL_P0 = PatternFill("solid", fgColor="DDEBF7")

AL_L = Alignment(horizontal="left", vertical="center", wrap_text=True)
AL_C = Alignment(horizontal="center", vertical="center", wrap_text=True)
AL_TL = Alignment(horizontal="left", vertical="top", wrap_text=True)
_thin = Side(style="thin", color="BBBBBB")
BORDER = Border(left=_thin, right=_thin, top=_thin, bottom=_thin)
BAND = {"A": "E2F0E4", "B": "E4EDF7", "기존점": "FDF2E0"}
DTFMT = "yyyy-mm-dd hh:mm:ss"


def title(ws, span, text, row=1):
    ws.merge_cells(f"A{row}:{span}{row}")
    c = ws[f"A{row}"]
    c.value, c.font, c.fill, c.alignment = text, F_TITLE, FILL_TITLE, AL_L
    ws.row_dimensions[row].height = 26


def note(ws, row, span, text, warn=False):
    ws.merge_cells(f"A{row}:{span}{row}")
    c = ws[f"A{row}"]
    c.value, c.alignment = text, AL_TL
    c.font = F_WARN if warn else F_SM
    if warn:
        c.fill = FILL_WARN
    segs = str(text).split("\n")
    ws.row_dimensions[row].height = 14 * sum(max(1, -(-len(s) // 84)) for s in segs) + 5


def table(ws, row, headers, widths, keys=()):
    for j, h in enumerate(headers, 1):
        c = ws.cell(row=row, column=j, value=h)
        c.font = F_KEY if h in keys else F_HDR
        c.fill = FILL_KEY if h in keys else FILL_TITLE
        c.alignment, c.border = AL_C, BORDER
    ws.row_dimensions[row].height = 36
    for j, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = ws.cell(row=row + 1, column=1)


def cell(ws, r, c, v=None, fill=FILL_FIELD, fmt=None, al=AL_C, font=F_VAL):
    x = ws.cell(row=r, column=c, value=v)
    x.border, x.fill, x.font, x.alignment = BORDER, fill, font, al
    if fmt:
        x.number_format = fmt
    return x


def blank_rows(ws, start, n, ncol, fmts=None):
    for r in range(start, start + n):
        for j in range(1, ncol + 1):
            cell(ws, r, j, fmt=(fmts or {}).get(j))


def dv(ws, ref, items):
    d = DataValidation(type="list", formula1='"' + ",".join(items) + '"',
                       allow_blank=True, showDropDown=False)
    ws.add_data_validation(d)
    d.add(ref)


def dv_warn(ws, ref, lo, hi, msg):
    """⚠️ 저장을 차단하지 않는다 — 실제 이상값을 잃지 않기 위해서다."""
    d = DataValidation(type="decimal", operator="between", formula1=lo, formula2=hi,
                       allow_blank=True, showErrorMessage=True, errorStyle="warning",
                       errorTitle="예상 범위 밖", error=msg)
    ws.add_data_validation(d)
    d.add(ref)


# ═══════════════════════════════════════════════ ① 기입규약
def sheet_intro(wb):
    ws = wb.create_sheet("① 기입규약")
    for col, w in zip("ABCDEFGH", [24, 20, 20, 20, 20, 20, 20, 20]):
        ws.column_dimensions[col].width = w
    title(ws, "H", "지자기 시험 탐사 표준 야장 (자기교란 확인표) — 기입규약")
    r = 2
    note(ws, r, "H", f"측정 설계·판정기준은 「{SP.DECK}」에 확정된 것을 그대로 쓴다. "
                     f"현행 규정과의 관계: {SP.REG_GAP}"); r += 1
    note(ws, r, "H", "야장에 칸이 없으면 그 값은 나중에 어떤 방법으로도 되살릴 수 없습니다. "
                     "실제로 과거 야장 68개를 전부 열어 보았더니 F 측정일시가 한 건도 적혀 "
                     "있지 않았고, 그 때문에 세션별 외부장 보정을 아예 할 수 없게 되었습니다. "
                     "이 양식이 이렇게 꼼꼼한 이유가 여기에 있습니다.", warn=True); r += 2

    def block(head, rows):
        nonlocal r
        ws.merge_cells(f"A{r}:H{r}")
        c = ws[f"A{r}"]
        c.value, c.font, c.fill, c.alignment = head, F_SEC, FILL_SEC, AL_L
        r += 1
        for k, v in rows:
            c1 = ws.cell(row=r, column=1, value=k)
            c1.font, c1.fill, c1.alignment, c1.border = F_LBL, FILL_LBL, AL_C, BORDER
            ws.merge_cells(f"B{r}:H{r}")
            c2 = ws.cell(row=r, column=2, value=v)
            c2.font, c2.alignment, c2.border = F_VAL, AL_L, BORDER
            # ⚠️ 병합하면 좌상단 셀에만 테두리가 남아 오른쪽 끝이 열려 보인다.
            #    범위 전체에 둘러야 인쇄물에서 칸이 닫힌다.
            for cc in range(2, 9):
                ws.cell(row=r, column=cc).border = BORDER
            ws.row_dimensions[r].height = 14 * max(1, -(-len(v) // 98)) + 6
            r += 1
        r += 1

    block("기록양식 3종 중 이 야장의 자리 (원 자료 23면)",
          [(f"{i+1}. {s}", f"{c} → 「{f}」 · {st}")
           for i, (s, c, f, st) in enumerate(SP.RECORD_FORMS)])
    block(f"측정 설계 — {SP.IAGA_SOURCE}", [
        ("장비 기준", "PPM(GSM-19T) **1대**. 자력계가 하나뿐이라 위치와 시각이 함께 움직인다 "
                    "— 그래서 중심점 전·후 반복관측이 필수다."),
        ("수평", f"P0 → 방향별 측정 → P0 재측정. 한 방향 "
                f"{', '.join(str(x) for x in SP.H_OFFSETS_M)} m → "
                f"{'·'.join(SP.H_DIRECTIONS)} 네 방향 반복."),
        ("수직", f"중심점(P0)에서 지상 {SP.V_HEIGHTS_CM[0]}~{SP.V_HEIGHTS_CM[-1]} cm 를 "
                f"{SP.V_HEIGHTS_CM[1]-SP.V_HEIGHTS_CM[0]} cm 간격으로 순차 측정."),
        ("시간변화 보정", "P0추정(t) = F0 + (F1−F0)·(t−t0)/(t1−t0) · ΔF = F측정 − P0추정 · "
                       "변화율 = ΔF ÷ 거리.  " + SP.TIME_CORR_NOTE),
        ("원 계획 예시", "P0 100 nT(10:00) → 동쪽 **거리 2 m** 지점 106 nT(10:02) → "
                      "P0 105 nT(10:04) ⇒ P0추정 102.5 · ΔF 3.5 · 변화율 "
                      "**1.75 nT/m**. 보정하지 않으면 3.0 nT/m 로 71 % 과대평가된다."),
        ("⚠ Set 검증", "유효한 P0(전)·P0(후)가 «각각 정확히 1개»이고 측정시각이 그 사이일 "
                     "때만 계산한다. 아니면 수치 대신 오류 사유가 뜬다 — 다른 방향의 P0 를 "
                     "물거나 중복 P0 를 합산하는 사고를 구조적으로 막는다."),
    ])
    block("판정기준 — 국제 권고이며 **국내 기준은 확정되지 않았다**", [
        ("총자기장 변화", SP.CRITERIA_SOURCE["range"]),
        ("자기구배", SP.CRITERIA_SOURCE["grad"]),
        ("⚠ 공간 범위", "4방향 측선은 원 내부의 «일부 선»만 잰다. 측선 사이와 대각선은 "
                      "관측하지 않으므로 「반경 10 m 전체가 50 nT 이내」라고 말할 수 없다. "
                      "「반경 10 m 내 4개 방위 측선의 관측점에서 확인된 변화」로 적는다."),
        ("⚠ 연산 미확정", "50 nT 를 max|F−P0| 로 볼지 관측값 max−min 으로 볼지, 지점 대표 "
                       "구배를 최대·중앙·특정 거리 중 무엇으로 할지 «정해지지 않았다». "
                       "야장은 셋 다 산출해 두고 확정은 국내 적용방안에 맡긴다."),
        ("⚠ 현실성", "시범관측 두 곳이 이미 네 방향 모두 3 nT/m 를 넘었다(이원 여의저수지 "
                   "최대 250.7 · 다른 기존점 최대 38.5, 1 m 기준). 두 지점만으로는 "
                   "「임계가 비현실적」·「1 m 가 전제와 다름」·「그 점이 실제 부적합」 셋을 "
                   "구분할 수 없다 — 자동 탈락시키지 말 것."),
        ("산출물", "합격/불합격만이 아니다. 함양 임도처럼 **중심점을 가장 균질한 지점으로 "
                 "조정**하는 것이 목적의 하나다(⑩). 옮겼으면 새 중심에서 재측정할지를 "
                 "명시해야 한다 — 옛 중심 자료가 새 중심을 입증하지는 않는다."),
    ])
    block("단위 — 논쟁 없이 고정한다", [
        ("각도(현장 원시)", "곤 gon (1회전 400 · 대척 200 gon). 도(°) 환산은 읽기전용 계산 칸."),
        ("총자력 F", "nT, 소수 1자리 (§17 0.1 nT 판독)"),
        ("좌표", "십진도. 위도·경도를 각각 다른 열에 숫자로 + 좌표계·취득방법·정확도"),
        ("거리·높이", "수평 거리 m, 수직 높이 cm. **명목값과 실측값을 따로** 적는다"),
        ("부호", "ΔF·변화율은 **부호를 보존**한다. 절댓값은 별도 열이다"),
    ])
    block("시각", [
        ("기입", "날짜 + 시:분:초 (예 2026-10-14 13:24:05). 분까지만 적으면 선형보간이 깨진다."),
        ("시간대", "야장은 KST. CYG·Kp 는 UT, KASA 캐시는 KST. UTC 환산은 계산 칸이 만든다."),
    ])
    block("방위각 — 정·역 혼동이 29쌍 있었다", [
        ("규약", "진북 기준 · 북 = 0° · 시계방향 · **FROM = 기준점, TO = 방위표지**"),
        ("참방위각", "결정방법·근거좌표·정확도와 **이전 확정값 대비 차이·재결정 사유** — "
                   "재방문 잔여 RMS 33.7′ 의 원인이 이것이다."),
    ])
    block("관측 조건 (원 자료 31면)", [(f"{i+1}", t) for i, t in enumerate(SP.OBS_CONDITION)])
    block("상태코드 — 빈칸을 0 으로 세지 않기 위해",
          list(zip(STATUS, ["자료가 존재하지 않음(예: 항공자력 미측선)", "관측 미실시",
                            "관측했으나 무효 — 사유 필수", "그 점에는 해당 없음",
                            "확인 대기"])) +
          [("공통", "상태코드에는 «사유·입력자·입력시각»을 반드시 함께 적는다. 같은 코드도 "
                  "장비오류·접근불가·기상·조사자 판단으로 원인이 다르다.")])
    block("색 규약", [("연한 노랑", "현장 기입란"), ("연한 초록", "사전 채움·계산값"),
                    ("연한 파랑", "중심점(P0) 반복관측 행"), ("회색", "편집 금지"),
                    ("연한 적색", "주의·경고")])
    return ws


# ═══════════════════════════════════════════════ ② 대상 50점
def sheet_points(wb, pts):
    ws = wb.create_sheet("② 대상 50점")
    n = {"A": 0, "B": 0, "기존점": 0}
    for p in pts:
        n[p["구분"]] += 1
    title(ws, "T", f"시험 탐사 대상 {len(pts)}점 — 등급 A {n['A']} · "
                   f"B 도엽대표 {n['B']} · 기존점 {n['기존점']}")
    note(ws, 2, "T", "선점 검토 결과 그대로다(export_selection_50 단일 출처). "
                     "Site ID 는 이 야장 전체의 연결키이며 다른 시트가 그대로 쓴다.")
    note(ws, 3, "T", "「예측 |∇ΔT|」는 KIGAM 항공자력 1.5분(약 2.8 km) 격자에서 계산한 값이라, "
                     "현장에서 수 m 간격으로 재는 값과는 물리적으로 다른 양입니다. 그래서 "
                     "합격 기준으로 쓰지 않고 사전 참고용으로만 둡니다. 오호·문산·강화·포천·"
                     "화천 다섯 곳은 항공자력 측선이 지나가지 않아 예측값이 아예 없는데, "
                     "이 경우 NOT AVAILABLE 로 두시면 됩니다. 값이 없다고 감점하거나 "
                     "이웃 값으로 채워 넣지는 마세요.", warn=True)
    hdr = ["Site ID", "연번", "구분", "원 번호", "지점명", "관할본부", "도엽번호", "도엽명",
           "위도(십진)", "경도(십진)", "표고(m)", "소재지",
           "예측 |∇ΔT|\n(nT/km)\n※판정 미사용", "예측 상태", "최근접\n관측소",
           "관측소\n거리(km)", "표지1 방위각\n(기재값)", "표지1\n거리(m)",
           "원 계획 유의사항", "HOLD 범위"]
    table(ws, 4, hdr, [11, 5, 7, 10, 16, 12, 13, 10, 11, 11, 9, 30, 12, 13, 9, 9,
                       14, 9, 18, 34], keys={"Site ID"})
    for i, p in enumerate(pts, 1):
        r, gap = 4 + i, p["예측구배"] is None
        hold = HOLD_SITES.get(p["지점명"])
        p["_sid"] = sid = f"S{i:02d}"
        vals = [sid, i, p["구분"], p["번호"], p["지점명"], p["관할본부"], p["도엽번호"],
                p["도엽명"], p["위도"], p["경도"], p["표고"], p["소재지"],
                None if gap else round(p["예측구배"], 1),
                "NOT AVAILABLE" if gap else "OK",
                p["최근접관측소"], round(p["관측소거리"]),
                p["표지1_방위각"] or "NOT APPLICABLE", p["표지1_거리"],
                SP.SITE_CAUTION.get(p["지점명"], ""),
                f"{hold[0]} 보류 — {hold[1]}" if hold else "없음 — 전 항목 진행"]
        fill = PatternFill("solid", fgColor=BAND[p["구분"]])
        for j, v in enumerate(vals, 1):
            cell(ws, r, j, v, fill=fill, al=AL_L if j in (12, 19, 20) else AL_C)
        for j, f in ((9, "0.000000"), (10, "0.000000"), (11, "0.0"),
                     (13, "0.0"), (18, "0.0")):
            ws.cell(row=r, column=j).number_format = f
        if gap:
            ws.cell(row=r, column=14).font = F_WARN
        if SP.SITE_CAUTION.get(p["지점명"]):
            ws.cell(row=r, column=19).font = F_WARN
        if hold:
            ws.cell(row=r, column=20).font = F_WARN
            ws.cell(row=r, column=20).fill = FILL_WARN
    return ws


# ═══════════════════════════════════════════════ ③ 점정보·좌표
def sheet_site(wb, pts):
    ws = wb.create_sheet("③ 점정보·좌표")
    title(ws, "M", "점 정보 · 좌표 검증")
    note(ws, 2, "M", "좌표는 «십진도 숫자»로 위도·경도를 각각 적고 좌표계·취득방법·정확도·"
                     "출처를 함께 남긴다. 자릿수만 맞추면 허위 정밀도가 생기고 위·경도 "
                     "뒤바뀜도 못 막는다. 선점 절차: " + " → ".join(SP.FIELD_FLOW))
    hdr = ["Site ID", "지점명", "구분", "사전 위도", "사전 경도",
           "현장 실측 위도", "현장 실측 경도", "좌표계", "측지기준", "취득방법",
           "수평정확도(m)", "좌표 출처", "좌표 검증상태"]
    table(ws, 4, hdr, [11, 15, 7, 12, 12, 13, 13, 12, 12, 14, 12, 18, 16],
          keys={"Site ID"})
    for i, p in enumerate(pts, 1):
        r = 4 + i
        for j, v in enumerate([p["_sid"], p["지점명"], p["구분"], p["위도"], p["경도"]], 1):
            cell(ws, r, j, v, fill=FILL_AUTO, fmt="0.000000" if j in (4, 5) else None)
        for j in range(6, len(hdr) + 1):
            cell(ws, r, j, fmt=("0.000000" if j in (6, 7) else
                                "0.00" if j == 11 else None))
    last = 4 + len(pts)
    dv(ws, f"H5:H{last}", ["WGS84", "KGD2002", "Bessel1841", "PENDING"])
    dv(ws, f"J5:J{last}", ["RTK-GNSS", "네트워크RTK", "정적GNSS", "휴대GNSS",
                           "기존 성과 인용", "PENDING"])
    dv(ws, f"M5:M{last}", ["검증완료 — 사전값과 일치", "검증완료 — 사전값 정정",
                           "불일치 — HOLD", "PENDING"])
    dv_warn(ws, f"F5:F{last}", 33, 39.5, "한반도 위도 범위 밖이다. 위·경도를 바꿔 적지 않았는지.")
    dv_warn(ws, f"G5:G{last}", 124, 132, "한반도 경도 범위 밖이다. 위·경도를 바꿔 적지 않았는지.")
    return ws


# ═══════════════════════════════════════════════ ④ §13② 이격
def sheet_clearance(wb, pts):
    """시설 «유형마다 한 행» — 단일 거리칸으로는 법정 판정을 재현할 수 없다(M4)."""
    ws = wb.create_sheet("④ §13② 이격거리")
    title(ws, "N", "§13② 이격거리 — 시설 유형마다 한 행")
    note(ws, 2, "N", "법정 기준: " + " · ".join(f"{n} {v} km" for n, v in FACILITIES))
    note(ws, 3, "N", "⚠ 도상 확인만으로 확정하지 않는다. 시설 분류·자료 기준일·거리 산정점·"
                     "최근 설치 여부에 따라 결과가 달라지므로 시설 ID·출처·좌표·현장 확인값·"
                     "확인방법·사진·확인자를 각각 남긴다. 분류가 불명확하면 임의 판정하지 말고 "
                     "PENDING.", warn=True)
    hdr = ["Site ID", "지점명", "시설 유형", "기준(km)", "최근접 시설 ID", "시설 명칭",
           "자료 출처", "자료 기준일", "시설 위도", "시설 경도", "도상거리(km)",
           "현장 확인거리(km)", "확인 방법", "사진 파일명", "판정", "확인자", "비고"]
    table(ws, 4, hdr, [11, 14, 15, 9, 15, 20, 16, 12, 12, 12, 12, 13, 14, 18, 10,
                       11, 24], keys={"Site ID", "시설 유형"})
    r = 5
    for p in pts:
        for nm, lim in FACILITIES:
            for j, v in enumerate([p["_sid"], p["지점명"], nm, lim], 1):
                cell(ws, r, j, v, fill=FILL_AUTO, fmt="0.0" if j == 4 else None)
            for j in range(5, len(hdr) + 1):
                cell(ws, r, j, fmt={8: "yyyy-mm-dd", 9: "0.000000", 10: "0.000000",
                                    11: "0.000", 12: "0.000"}.get(j),
                     al=AL_L if j in (6, 17) else AL_C)
            r += 1
    dv(ws, f"M5:M{r-1}", ["도상만", "현장 실측", "지자체 확인", "PENDING"])
    dv(ws, f"O5:O{r-1}", ["적합", "부적합", "PENDING", "NOT APPLICABLE"])
    return ws


# ═══════════════════════════════════════════════ ⑤ 기존점 대조가능성
def sheet_existing(wb, pts):
    ws = wb.create_sheet("⑤ 기존점 대조가능성")
    ex = [p for p in pts if p["구분"] == "기존점"]
    title(ws, "Q", f"기존 관측점 {len(ex)}점 — 과거 성과와의 «성분별» 대조 가능성")
    note(ws, 2, "Q", "기존 성과가 남아 있다는 것과 그 값을 지금 값과 견줄 수 있다는 것은 다른 "
                     "이야기입니다. 과거 기록에는 방문할 때마다 달라진 마크 방위각(재방문 "
                     "잔여 RMS 33.7분), 아예 적히지 않은 F 측정시각, 정반대로 적힌 방위 "
                     "방향이 섞여 있습니다. 그래서 편각 D 와, 시간보정이 필요한 F 는 그냥 "
                     "빼서 비교할 수가 없습니다.", warn=True)
    note(ws, 3, "Q", "HOLD 는 범위를 나눈다 — 「현 물리점 식별」이 안 되면 전부 보류지만, "
                     "현 표석이 확실하고 과거 좌표와의 연결만 불명확하면 **신규 측정은 별도 "
                     "Location ID 로 진행**하고 과거 대조만 보류한다.")
    hdr = ["Site ID", "지점명", "지점코드", "종전 성과 연도", "종전 편각(°)", "종전 복각(°)",
           "종전 총자력(nT)", "종전 관측장비",
           "현 물리점 식별", "공식 좌표 확인", "과거점 동일성", "동일 표석", "동일 표지",
           "D 대조가능", "I 대조가능", "F 대조가능", "근거·비고"]
    table(ws, 4, hdr, [11, 12, 11, 13, 12, 12, 13, 13, 13, 13, 13, 11, 11, 11, 11,
                       11, 32], keys={"Site ID"})
    for i, p in enumerate(ex, 1):
        r, s = 4 + i, (p["종전성과"] or {})
        hold = HOLD_SITES.get(p["지점명"])
        for j, v in enumerate([p["_sid"], p["지점명"], p["번호"], s.get("최종관측", ""),
                               s.get("편각"), s.get("복각"), s.get("총자력"),
                               s.get("관측장비", "")], 1):
            cell(ws, r, j, v, fill=FILL_AUTO,
                 fmt={5: "0.0000", 6: "0.0000", 7: "#,##0.0"}.get(j))
        for j in range(9, len(hdr) + 1):
            # ⚠️ 미원은 «과거대조»만 보류다 — 현 물리점 식별·신규 측정은 막지 않는다.
            v = "HOLD" if (hold and 11 <= j <= 16) else None
            cell(ws, r, j, v, fill=FILL_WARN if v else FILL_FIELD,
                 al=AL_L if j == 17 else AL_C)
        if hold:
            ws.cell(row=r, column=17, value=hold[1]).alignment = AL_L
    last = 4 + len(ex)
    dv(ws, f"I5:K{last}", ["확인", "불일치", "PENDING", "HOLD"])
    dv(ws, f"L5:M{last}", ["확인", "불일치", "PENDING", "HOLD"])
    dv(ws, f"N5:P{last}", ["대조가능", "조건부", "불가", "PENDING", "HOLD"])
    return ws


# ═══════════════════════════════════════════════ ⑥ 방위표지
def sheet_mark(wb, pts):
    ws = wb.create_sheet("⑥ 방위표지·참방위각")
    title(ws, "T", "방위표지 · 참방위각 결정 — 편각 30′ 잔차의 주범")
    note(ws, 2, "T", "방위각은 진북을 기준으로 북쪽이 0도, 시계 방향으로 잽니다. 방향은 언제나 "
                     "기준점에서 방위표지를 본 쪽입니다(FROM 기준점 · TO 표지). 현장 원시각은 "
                     "곤(gon)이라 한 바퀴가 400 이고 반대쪽이 200 인데, 여기에 도 단위의 "
                     "180 이 섞여 들어오면 어느 쪽이 반대 방향인지 알 수 없게 되어 대규모 "
                     "오류로 이어집니다.", warn=True)
    note(ws, 3, "T", "재방문 16구간 전부에서 마크 참방위각이 60″ 넘게 변했다. 마크를 바꾸는 것 "
                     "자체는 오류가 아니지만 바뀐 마크의 참방위각이 틀리면 그 방문의 편각이 "
                     "통째로 틀어진다 — 「이전값 대비 차」와 「재결정 사유」가 필수다.")
    hdr = ["Site ID", "Mark ID", "지점명", "표지 구분", "표지 위도", "표지 경도",
           "시준 지점(표지의 어디)", "FROM", "TO",
           "정 시준 원시각(gon)", "반 시준 원시각(gon)", "폐합차(gon)\n계산",
           "참방위각(gon)", "참방위각(°)\n계산", "역방위(°)\n계산",
           "결정방법", "근거좌표 정확도(m)", "이전 확정 방위각(°)",
           "이전값 대비 차(″)\n계산", "재결정 사유"]
    table(ws, 4, hdr, [11, 11, 14, 10, 12, 12, 20, 12, 12, 15, 15, 12, 13, 13, 13,
                       16, 14, 15, 14, 30], keys={"Site ID", "Mark ID"})
    r = 5
    for p in pts:
        for k in (1, 2):
            for j, v in enumerate([p["_sid"], f"{p['_sid']}-M{k}", p["지점명"],
                                   f"방위표지{k}", p.get(f"표지{k}_위도"),
                                   p.get(f"표지{k}_경도")], 1):
                cell(ws, r, j, v, fill=FILL_AUTO,
                     fmt="0.000000" if j in (5, 6) and v is not None else None)
            cell(ws, r, 8, "기준점", fill=FILL_LOCK)
            cell(ws, r, 9, "방위표지", fill=FILL_LOCK)
            for j in (7, 10, 11, 13, 16, 17, 18, 20):
                cell(ws, r, j, fmt={10: "0.0000", 11: "0.0000", 13: "0.0000",
                                    17: "0.00", 18: "0.00"}.get(j))
            cell(ws, r, 12, f'=IF(COUNT(J{r}:K{r})=2,ABS(ABS(K{r}-J{r})-200),"")',
                 fill=FILL_AUTO, fmt="0.0000")
            cell(ws, r, 14, f'=IF(ISNUMBER(M{r}),M{r}*0.9,"")', fill=FILL_AUTO, fmt="0.0000")
            cell(ws, r, 15, f'=IF(ISNUMBER(M{r}),MOD(M{r}*0.9+180,360),"")',
                 fill=FILL_AUTO, fmt="0.0000")
            cell(ws, r, 19, f'=IF(AND(ISNUMBER(M{r}),ISNUMBER(R{r})),'
                            f'(M{r}*0.9-R{r})*3600,"")', fill=FILL_AUTO, fmt="0.0")
            r += 1
    dv(ws, f"P5:P{r-1}", ["천문관측", "자이로", "RTK 장기선", "좌표계산",
                          "기존 성과 인용", "PENDING"])
    dv_warn(ws, f"L5:L{r-1}", 0, 0.02, "마크 폐합차가 0.02 gon(약 65″)을 넘는다. 시준 재확인.")
    return ws


# ═══════════════════════════════════════════════ ⑦ 장비·검교정
def sheet_instrument(wb):
    ws = wb.create_sheet("⑦ 장비·검교정")
    for col, w in zip("ABCDEFGHIJ", [16, 18, 30, 34, 14, 14, 14, 14, 16, 26]):
        ws.column_dimensions[col].width = w
    title(ws, "J", "관측 장비 · 검교정 이력")
    note(ws, 2, "J", "검교정 주기는 KOLAS-G-008 고시(국가표준기본법 §14) 기준이다. 성적서 "
                     "번호와 유효기간을 적지 않으면 그 관측의 신뢰성을 사후에 입증할 수 없다.")
    note(ws, 3, "J", "⚠ 이동 자력계는 GSM-19T **1대**가 전제다. 두 번째 PPM 을 요구하는 절차"
                     "(교차배치 등)는 «2대 확보 시에만» 수행하며 미확보 시 NOT APPLICABLE 로 "
                     "남긴다.", warn=True)
    r = 4
    ws.merge_cells(f"A{r}:J{r}")
    c = ws[f"A{r}"]
    c.value, c.font, c.fill, c.alignment = "■ 보유 장비 (원 자료 31면)", F_SEC, FILL_SEC, AL_L
    r += 1
    table(ws, r, ["구분", "모델", "용도", "유의사항", "기기번호", "검정일", "유효기간",
                  "검정기관", "성적서 번호", "시험 탐사 범위"],
          [16, 18, 30, 34, 14, 14, 14, 14, 16, 26])
    r += 1
    for kind, model, use, caution, scope in SP.INSTRUMENTS:
        for j, v in enumerate([kind, model, use, caution], 1):
            cell(ws, r, j, v, fill=FILL_AUTO, al=AL_L if j in (3, 4) else AL_C)
        for j in range(5, 10):
            cell(ws, r, j, fmt="yyyy-mm-dd" if j in (6, 7) else None)
        # ⚠️ 절대측정 장비는 «범위 밖»으로 표시한다 — 목록에서 지우지 않는 것은
        #    나중 관측 단계에서 쓸 장비이고 검교정 계획이 이어져야 하기 때문이다.
        cell(ws, r, 10, scope, al=AL_L,
             fill=FILL_WARN if "범위 밖" in scope else FILL_AUTO,
             font=F_WARN if "범위 밖" in scope else F_VAL)
        ws.row_dimensions[r].height = 32
        r += 1
    r += 1
    ws.merge_cells(f"A{r}:J{r}")
    c = ws[f"A{r}"]
    c.value = "■ 검교정 주기 — KOLAS-G-008 고시 · 국가표준기본법 §14"
    c.font, c.fill, c.alignment = F_SEC, FILL_SEC, AL_L
    r += 1
    table(ws, r, ["장비 분류", "분류번호", "주기", "적용 기기번호", "직전 검정일",
                  "다음 예정일", "상태", "", "", ""],
          [16, 18, 30, 34, 14, 14, 14, 14, 16, 26])
    r += 1
    cal_start = r
    for nm, code, cyc in SP.CALIB_CYCLE:
        for j, v in enumerate([nm, code, cyc], 1):
            cell(ws, r, j, v, fill=FILL_AUTO)
        for j in range(4, 8):
            cell(ws, r, j, fmt="yyyy-mm-dd" if j in (5, 6) else None)
        r += 1
    dv(ws, f"G{cal_start}:G{r-1}", ["유효", "만료 임박", "만료", "PENDING"])
    return ws


# ═══════════════════════════════════════════════ ⑧ 관측조건·변동관측
def sheet_condition(wb):
    ws = wb.create_sheet("⑧ 관측조건·변동관측")
    title(ws, "P", "관측 조건 · 이동식 Variometer 연속기록")
    note(ws, 2, "P", "P0 를 앞뒤로 재어 그 사이를 직선으로 잇는 방식은, 재는 동안 중심점 "
                     "자기장이 고르게 변했다고 «가정»하는 것입니다. 양 끝 두 점만 가지고는 "
                     "중간에 휘었는지, 순간적으로 튀었는지를 알아낼 방법이 없습니다. "
                     "Variometer 로 연속기록을 남겨야 그 가정이 맞았는지 확인할 수 있어서, "
                     "연속기록이 없는 세트는 ⑬에서 「선형성 미검증」으로 따로 표시합니다.",
         warn=True)
    note(ws, 3, "P", "최근접 상시관측소는 중앙 101 km · 최대 229 km(거제)이고, 그 거리에서 "
                     "「외부장이 전역 균일」 근사가 깨짐을 이 프로젝트가 실증했다. "
                     "관측 조건: " + " / ".join(SP.OBS_CONDITION))
    hdr = ["Site ID", "Visit ID", "관측일자", "야간 관측 여부", "Kp 확인값",
           "우주기상 예보 확인", "Variometer 파일 ID", "위치 ID", "위도", "경도",
           "기록 간격(s)", "기록 시작(KST)", "기록 종료(KST)", "자료공백 구간",
           "시계 동기 방법", "동기 오차(s)"]
    table(ws, 4, hdr, [11, 10, 12, 13, 11, 15, 18, 12, 11, 11, 11, 19, 19, 20,
                       18, 12], keys={"Site ID", "Visit ID", "Variometer 파일 ID"})
    blank_rows(ws, 5, 80, len(hdr),
               fmts={3: "yyyy-mm-dd", 5: "0.0", 9: "0.000000", 10: "0.000000",
                     11: "0", 12: DTFMT, 13: DTFMT, 16: "0.0"})
    dv(ws, "D5:D84", ["야간", "주간", "혼합"])
    dv(ws, "F5:F84", ["확인 — 정온", "확인 — 교란 예보", "미확인"])
    dv(ws, "O5:O84", ["GNSS 시각 동기", "NTP 동기", "수동 대조(오차 기록)",
                      "미동기 — INVALID", "NOT OBSERVED"])
    return ws


# ═══════════════════════════════════════════════ ⑨ 수평 자기구배
def sheet_horizontal(wb, pts):
    ws = wb.create_sheet("⑨ 수평 자기구배")
    title(ws, "AC", f"수평 자기구배 — {SP.IAGA_SOURCE}")
    note(ws, 2, "AC", f"P0 → 방향별 측정 → P0 재측정. 한 방향 "
                      f"{', '.join(str(x) for x in SP.H_OFFSETS_M)} m 를 재고 "
                      f"{'·'.join(SP.H_DIRECTIONS)} 네 방향을 반복한다. 방향마다 "
                      f"P0(전)·P0(후) 행이 따로 있고 그 둘의 «선형보간»으로 각 측정 "
                      f"시각의 중심점을 추정한다.")
    note(ws, 3, "AC", "오른쪽 계산 칸은 「Set 검증」이 OK 로 나올 때만 값을 냅니다. 같은 "
                      "Set ID·같은 차수 안에 유효한 P0 가 앞뒤로 각각 하나씩 있어야 하고 "
                      "측정시각이 그 사이에 들어와야 하는데, 하나라도 어긋나면 숫자 대신 "
                      "무엇이 잘못됐는지가 표시됩니다. ΔF 와 변화율은 부호를 그대로 두었고 "
                      "절댓값은 옆 칸에 따로 두었습니다. 「중심점기준 변화율」은 중심점과 "
                      "그 측점 사이의 평균 기울기이고 「인접구간 변화율」은 이웃한 두 측점 "
                      "사이의 기울기라, 서로 다른 값이니 섞어 쓰지 마세요.", warn=True)
    note(ws, 4, "AC", "현장에서 그 방향을 다시 쟀다면 차수를 2, 3 으로 올려 «행을 새로 "
                      "추가»해 주세요. 1차 행을 고쳐 쓰지 않는 것이 중요합니다. 차수를 "
                      "비워 두고 같은 Set ID 로 행만 늘리면 P0 가 두 벌이 되어 1차 계산까지 "
                      "함께 막힙니다.", warn=True)
    hdr = ["Site ID", "지점명", "Set ID", "차수", "구분", "방향",
           "명목거리(m)", "실측거리(m)", "측정 시각(KST)",
           "F 판독1(nT)", "F 판독2", "F 판독3", "중앙값\n계산",
           "Set 검증\n계산", "P0전 시각\n계산", "P0전 F\n계산",
           "P0후 시각\n계산", "P0후 F\n계산", "P0 추정값\n계산", "ΔF(nT)\n계산",
           "중심점기준\n변화율(nT/m)\n계산", "|변화율|\n계산",
           "인접구간\n변화율(nT/m)\n계산", "Set 지속(분)\n계산",
           "기기번호", "Variometer\n파일 ID", "Set 상태", "중단 시각·사유", "상태·비고"]
    table(ws, 5, hdr, [11, 13, 13, 7, 10, 7, 10, 10, 19, 12, 11, 11, 11, 16, 17,
                       11, 17, 11, 12, 11, 13, 11, 13, 12, 11, 14, 13, 22, 24],
          keys={"Site ID", "Set ID", "차수", "구분"})
    r = 6
    for p in pts:
        for d in SP.H_DIRECTIONS:
            setid = f"{p['_sid']}-H{d}"
            seq = ([("P0(전)", None)] + [("측정", o) for o in SP.H_OFFSETS_M]
                   + [("P0(후)", None)])
            for kind, off in seq:
                isP0 = kind != "측정"
                base = FILL_P0 if isP0 else FILL_AUTO
                for j2, v in enumerate([p["_sid"], p["지점명"], setid, 1, kind, d,
                                        off], 1):
                    cell(ws, r, j2, v, fill=base,
                         fmt="0.0" if j2 == 7 else ("0" if j2 == 4 else None))
                for j2 in (8, 9, 10, 11, 12, 25, 26, 27, 28, 29):
                    cell(ws, r, j2, fmt={8: "0.00", 9: DTFMT, 10: "0.0",
                                         11: "0.0", 12: "0.0"}.get(j2),
                         al=AL_L if j2 in (28, 29) else AL_C)
                cell(ws, r, 13, f'=IF(COUNT(J{r}:L{r})=0,"",MEDIAN(J{r}:L{r}))',
                     fill=FILL_AUTO, fmt="0.0")
                dist = f'IF(ISNUMBER(H{r}),H{r},G{r})'
                if isP0:
                    for j2 in range(14, 25):
                        cell(ws, r, j2, fill=FILL_LOCK)
                else:
                    # ⚠️ 조회 키에 «차수»($D)가 들어간다. 이게 없으면 재측정 행이
                    #    1차 행과 같은 Set 으로 묶여 둘 다 계산이 막힌다.
                    for j2, (col, tag) in enumerate(
                            [("$I", "P0(전)"), ("$M", "P0(전)"),
                             ("$I", "P0(후)"), ("$M", "P0(후)")], 15):
                        f = (f'=IF(COUNTIFS($C:$C,$C{r},$D:$D,$D{r},$E:$E,"{tag}",'
                             f'$M:$M,">0")<>1,"",'
                             f'SUMIFS({col}:{col},$C:$C,$C{r},$D:$D,$D{r},'
                             f'$E:$E,"{tag}"))')
                        cell(ws, r, j2, f, fill=FILL_AUTO,
                             fmt=DTFMT if j2 in (15, 17) else "0.0")
                    cell(ws, r, 14,
                         f'=IF(COUNTIFS($C:$C,$C{r},$D:$D,$D{r},$E:$E,"P0(전)",'
                         f'$M:$M,">0")<>1,"P0(전) 유효 1개 아님",'
                         f'IF(COUNTIFS($C:$C,$C{r},$D:$D,$D{r},$E:$E,"P0(후)",'
                         f'$M:$M,">0")<>1,"P0(후) 유효 1개 아님",'
                         f'IF(NOT(AND(ISNUMBER(O{r}),ISNUMBER(Q{r}),'
                         f'ISNUMBER(I{r}))),"시각 결측",'
                         f'IF(Q{r}<=O{r},"P0 시각 역전",'
                         f'IF(OR(I{r}<O{r},I{r}>Q{r}),"측정시각 구간 밖","OK")))))',
                         fill=FILL_AUTO)
                    cell(ws, r, 19,
                         f'=IF(N{r}<>"OK","",P{r}+(R{r}-P{r})*(I{r}-O{r})/(Q{r}-O{r}))',
                         fill=FILL_AUTO, fmt="0.0")
                    cell(ws, r, 20,
                         f'=IF(OR(N{r}<>"OK",NOT(ISNUMBER(M{r}))),"",M{r}-S{r})',
                         fill=FILL_AUTO, fmt="0.0")
                    cell(ws, r, 21,
                         f'=IF(OR(NOT(ISNUMBER(T{r})),{dist}=0),"",T{r}/{dist})',
                         fill=FILL_AUTO, fmt="0.00")
                    cell(ws, r, 22, f'=IF(ISNUMBER(U{r}),ABS(U{r}),"")',
                         fill=FILL_AUTO, fmt="0.00")
                    prev = r - 1
                    pdist = f'IF(ISNUMBER(H{prev}),H{prev},G{prev})'
                    cell(ws, r, 23,
                         f'=IF(OR($E{prev}<>"측정",$C{prev}<>$C{r},$D{prev}<>$D{r},'
                         f'NOT(ISNUMBER(M{r})),NOT(ISNUMBER(M{prev})),'
                         f'({dist}-{pdist})=0),"",(M{r}-M{prev})/({dist}-{pdist}))',
                         fill=FILL_AUTO, fmt="0.00")
                    cell(ws, r, 24,
                         f'=IF(NOT(AND(ISNUMBER(O{r}),ISNUMBER(Q{r}))),"",'
                         f'(Q{r}-O{r})*1440)', fill=FILL_AUTO, fmt="0.0")
                r += 1
    last = r - 1
    dv(ws, f"D6:D{last}", ["1", "2", "3"])
    dv(ws, f"AA6:AA{last}", ["정상", "중단", "재시작", "무효", "PENDING"])
    dv(ws, f"AC6:AC{last}", ["정상"] + STATUS)
    dv_warn(ws, f"J6:L{last}", 30000, 60000,
            "한반도 총자력 범위 밖이다. 실제 이상일 수 있으니 지우지 말고 사유를 비고에.")
    return ws, last


# ═══════════════════════════════════════════════ ⑩ 수직 자기구배
def sheet_vertical(wb, pts):
    ws = wb.create_sheet("⑩ 수직 자기구배")
    title(ws, "X", "수직 자기구배 — 중심점(P0) 높이별 측정")
    note(ws, 2, "X", f"중심점에서 지상 {SP.V_HEIGHTS_CM[0]}~{SP.V_HEIGHTS_CM[-1]} cm 를 "
                     f"{SP.V_HEIGHTS_CM[1]-SP.V_HEIGHTS_CM[0]} cm 간격으로 순차 측정한다"
                     f"({SP.IAGA_SOURCE}). 기준높이에서 전·후 반복관측해 시간변화를 뺀다.")
    note(ws, 3, "X", "유럽 권고는 수직 자기구배를 점검하라고만 하고 얼마 이하여야 한다는 "
                     "수치는 주지 않았습니다. 그래서 이 시트는 자동으로 합격·불합격을 "
                     "가르지 않고, 보정한 F 프로파일과 전 구간 변화폭, 이웃 높이 사이의 "
                     "변화 같은 기술지표만 냅니다. 높이는 목표값이 아니라 센서 가운데가 "
                     "실제로 어디였는지와 어디를 기준면으로 삼았는지를 적어 주셔야 나중에 "
                     "같은 조건을 재현할 수 있습니다.", warn=True)
    note(ws, 4, "X", "다시 쟀다면 차수를 올려 행을 새로 추가해 주세요. 1차 행을 고쳐 쓰거나 "
                     "차수 없이 행만 늘리면 기준 관측이 두 벌이 되어 1차 계산까지 막힙니다.",
         warn=True)
    hdr = ["Site ID", "지점명", "Set ID", "차수", "구분", "명목 높이(cm)",
           "실제 센서 중심높이(cm)", "측정 기준면", "측정 시각(KST)",
           "F 판독1(nT)", "F 판독2", "F 판독3", "중앙값\n계산", "Set 검증\n계산",
           "기준(전) 시각\n계산", "기준(전) F\n계산", "기준(후) 시각\n계산",
           "기준(후) F\n계산", "기준 추정값\n계산", "ΔF(nT)\n계산",
           "인접높이 변화율\n(nT/m)\n계산", "기기번호", "Set 상태", "상태·비고"]
    table(ws, 5, hdr, [11, 13, 12, 7, 11, 11, 15, 13, 19, 12, 11, 11, 11, 16, 17,
                       12, 17, 12, 12, 11, 14, 11, 12, 24],
          keys={"Site ID", "Set ID", "차수", "구분"})
    r = 6
    for p in pts:
        setid = f"{p['_sid']}-V"
        seq = ([("기준(전)", None)] + [("측정", h) for h in SP.V_HEIGHTS_CM]
               + [("기준(후)", None)])
        for kind, h in seq:
            isR = kind != "측정"
            base = FILL_P0 if isR else FILL_AUTO
            for j2, v in enumerate([p["_sid"], p["지점명"], setid, 1, kind, h], 1):
                cell(ws, r, j2, v, fill=base,
                     fmt="0" if j2 in (4, 6) else None)
            for j2 in (7, 8, 9, 10, 11, 12, 22, 23, 24):
                cell(ws, r, j2, fmt={7: "0.0", 9: DTFMT, 10: "0.0", 11: "0.0",
                                     12: "0.0"}.get(j2),
                     al=AL_L if j2 == 24 else AL_C)
            cell(ws, r, 13, f'=IF(COUNT(J{r}:L{r})=0,"",MEDIAN(J{r}:L{r}))',
                 fill=FILL_AUTO, fmt="0.0")
            hgt = f'IF(ISNUMBER(G{r}),G{r},F{r})'
            if isR:
                for j2 in range(14, 22):
                    cell(ws, r, j2, fill=FILL_LOCK)
            else:
                for j2, (col, tag) in enumerate(
                        [("$I", "기준(전)"), ("$M", "기준(전)"),
                         ("$I", "기준(후)"), ("$M", "기준(후)")], 15):
                    f = (f'=IF(COUNTIFS($C:$C,$C{r},$D:$D,$D{r},$E:$E,"{tag}",'
                         f'$M:$M,">0")<>1,"",'
                         f'SUMIFS({col}:{col},$C:$C,$C{r},$D:$D,$D{r},'
                         f'$E:$E,"{tag}"))')
                    cell(ws, r, j2, f, fill=FILL_AUTO,
                         fmt=DTFMT if j2 in (15, 17) else "0.0")
                cell(ws, r, 14,
                     f'=IF(COUNTIFS($C:$C,$C{r},$D:$D,$D{r},$E:$E,"기준(전)",'
                     f'$M:$M,">0")<>1,"기준(전) 유효 1개 아님",'
                     f'IF(COUNTIFS($C:$C,$C{r},$D:$D,$D{r},$E:$E,"기준(후)",'
                     f'$M:$M,">0")<>1,"기준(후) 유효 1개 아님",'
                     f'IF(NOT(AND(ISNUMBER(O{r}),ISNUMBER(Q{r}),ISNUMBER(I{r}))),'
                     f'"시각 결측",'
                     f'IF(Q{r}<=O{r},"기준 시각 역전",'
                     f'IF(OR(I{r}<O{r},I{r}>Q{r}),"측정시각 구간 밖","OK")))))',
                     fill=FILL_AUTO)
                cell(ws, r, 19,
                     f'=IF(N{r}<>"OK","",P{r}+(R{r}-P{r})*(I{r}-O{r})/(Q{r}-O{r}))',
                     fill=FILL_AUTO, fmt="0.0")
                cell(ws, r, 20,
                     f'=IF(OR(N{r}<>"OK",NOT(ISNUMBER(M{r}))),"",M{r}-S{r})',
                     fill=FILL_AUTO, fmt="0.0")
                prev = r - 1
                phgt = f'IF(ISNUMBER(G{prev}),G{prev},F{prev})'
                cell(ws, r, 21,
                     f'=IF(OR($E{prev}<>"측정",$C{prev}<>$C{r},$D{prev}<>$D{r},'
                     f'NOT(ISNUMBER(M{r})),NOT(ISNUMBER(M{prev})),'
                     f'({hgt}-{phgt})=0),"",(M{r}-M{prev})/(({hgt}-{phgt})/100))',
                     fill=FILL_AUTO, fmt="0.00")
            r += 1
    dv(ws, f"D6:D{r-1}", ["1", "2", "3"])
    dv(ws, f"H6:H{r-1}", ["지면", "표석 상면", "PENDING"])
    dv(ws, f"W6:W{r-1}", ["정상", "중단", "재시작", "무효", "PENDING"])
    dv(ws, f"X6:X{r-1}", ["정상"] + STATUS)
    dv_warn(ws, f"J6:L{r-1}", 30000, 60000, "한반도 총자력 범위 밖이다.")
    return ws


# ═══════════════════════════════════════════════ ⑪ 중심점 최종선정
def sheet_center(wb, pts):
    ws = wb.create_sheet("⑪ 중심점 최종선정")
    title(ws, "P", "중심점(P0) 최종 선정 — 옮겼으면 새 중심에서 다시 재야 한다")
    note(ws, 2, "P", "중심점을 옮기면 거리도 방향도 수직 프로파일도 모두 기준이 달라집니다. "
                     "옮기기 전 자리에서 잰 결과가 새 자리의 10 m 범위나 수직구배를 "
                     "말해 주지는 못합니다. 그래서 처음 후보와 최종 중심에 서로 다른 "
                     "Location ID 를 주고, 새 자리에서 전체를 다시 쟀는지 아니면 "
                     "「잠정 중심 · 확정측정 미실시」로 남겨 두는지를 분명히 적습니다.",
         warn=True)
    note(ws, 3, "P", f"선점 절차: {' → '.join(SP.FIELD_FLOW)} · "
                     f"재관측 주기(안) {SP.REVISIT_YEARS}년. {SP.REVISIT_NOTE}")
    hdr = ["Site ID", "지점명", "최초 후보 Location ID", "사전 위도", "사전 경도",
           "중심점 조정 여부", "최종 중심 Location ID", "최종 위도", "최종 경도",
           "이동거리(m)", "이동 방위(°)", "조정 사유",
           "새 중심 재측정 여부", "표지 설치 여부", "표지 설치일", "설치자"]
    table(ws, 4, hdr, [11, 14, 18, 12, 12, 13, 18, 12, 12, 11, 11, 30, 16, 12,
                       12, 12], keys={"Site ID", "최초 후보 Location ID",
                                      "최종 중심 Location ID"})
    for i, p in enumerate(pts, 1):
        r = 4 + i
        for j, v in enumerate([p["_sid"], p["지점명"], f"{p['_sid']}-P0a",
                               p["위도"], p["경도"]], 1):
            cell(ws, r, j, v, fill=FILL_AUTO, fmt="0.000000" if j in (4, 5) else None)
        for j in range(6, len(hdr) + 1):
            cell(ws, r, j, al=AL_L if j == 12 else AL_C,
                 fmt={8: "0.000000", 9: "0.000000", 10: "0.0", 11: "0.0",
                      15: "yyyy-mm-dd"}.get(j))
    last = 4 + len(pts)
    dv(ws, f"F5:F{last}", ["조정 없음", "조정함", "PENDING"])
    dv(ws, f"M5:M{last}", ["새 중심에서 전체 재측정 완료", "일부 재측정",
                           "잠정 중심 — 확정측정 미실시", "NOT APPLICABLE", "PENDING"])
    dv(ws, f"N5:N{last}", ["설치 완료", "미설치", "PENDING", "NOT APPLICABLE"])
    dv_warn(ws, f"H5:H{last}", 33, 39.5, "한반도 위도 범위 밖이다.")
    dv_warn(ws, f"I5:I{last}", 124, 132, "한반도 경도 범위 밖이다.")
    return ws


# ═══════════════════════════════════════════════ ⑫ 게이트
def sheet_gate(wb, pts):
    ws = wb.create_sheet("⑫ 게이트 점검")
    title(ws, "N", "기록 완전성 점검 — 이 탐사는 법정 선점이 아니다")
    note(ws, 2, "N", "이 탐사는 법정 선점(§14 선점실시)이 아니고 통과 기준도 아직 없습니다"
                     "(2026-09-07 발주자 결정). 아래 칸 이름에 붙은 조항 번호는 어떤 "
                     "지침을 참고했는지를 밝힌 것일 뿐, 그 칸에 「통과」를 적었다고 해서 "
                     "법정 요건을 충족했다는 뜻은 아닙니다. 이 시트가 실제로 묻는 것은 "
                     "하나뿐입니다 — 나중에 분석할 수 있을 만큼 기록이 갖춰졌는가. "
                     "정량 지표는 ⑬에서 봅니다.", warn=True)
    note(ws, 3, "N", "오히려 반대다 — **국내 기준이 없으므로 그 기준을 만들 자료를 "
                     "모으는 것**이 이 탐사의 목적 가운데 하나다. 그래서 참고값"
                     "(IAGA 50 nT · 유럽 3 nT/m)을 넘었다고 탈락시키지 않고 그대로 "
                     "기록한다. 정식 선점·표지 설치로 넘어갈 때 §14·§15·§16·§18 의 "
                     "적용범위를 그때 다시 판단한다.")
    # ⚠️ §19①(1일 6회)·§20(정수차 30′·20분)은 **절대측정의 요건**이라 이 단계에
    #    해당하지 않는다. 시험 탐사는 오버하우저 총자력만 재고, 절대관측은 선점이
    #    끝난 뒤의 별도 단계다(원 계획 31면 선점실시 절차).
    # ⚠️ 조문 번호는 «참고한 지침»의 출처 표시일 뿐 법정 판정이 아니다.
    hdr = ["Site ID", "지점명", "이격 확인\n(§13② 참고)",
           "자기환경 측정 실시\n(§14 참고)", "0.1 nT 판독\n(§17 참고)",
           "자기오염 방지\n(§18③ 참고)", "시간변화 보정자료\n(§21 참고)",
           "수평 4방향 Set 검증 OK", "수직 Set 검증 OK", "P0 전·후 반복관측",
           "참방위각 결정", "좌표 검증", "중심점 확정", "게이트 종합"]
    table(ws, 4, hdr, [11, 15, 11, 17, 13, 14, 14, 15, 13, 13, 12, 11, 12, 13],
          keys={"Site ID"})
    n_gate = len(hdr) - 3
    for i, p in enumerate(pts, 1):
        r = 4 + i
        for j, v in enumerate([p["_sid"], p["지점명"]], 1):
            cell(ws, r, j, v, fill=FILL_AUTO)
        for j in range(3, len(hdr)):
            cell(ws, r, j)
        cell(ws, r, len(hdr),
             f'=IF(COUNTIF(C{r}:M{r},"실패")>0,"실패",'
             f'IF(COUNTIF(C{r}:M{r},"통과")={n_gate},"통과","미완"))', fill=FILL_AUTO)
    last = 4 + len(pts)
    dv(ws, f"C5:M{last}", ["통과", "실패", "PENDING", "NOT APPLICABLE"])
    return ws


# ═══════════════════════════════════════════════ ⑬ 지표·판정
def sheet_result(wb, pts):
    ws = wb.create_sheet("⑬ 지표·판정")
    title(ws, "T", "정량 지표 · 시험 탐사 판정")
    note(ws, 2, "T", f"참고 기준: {SP.CRITERIA_SOURCE['range']} / {SP.CRITERIA_SOURCE['grad']}")
    note(ws, 3, "T", "「관측사실」과 「국제 참고기준 대비 결과」와 「국내 판정상태」를 굳이 "
                     "나눠 놓은 이유가 있습니다. 국내 기준이 없는데 판정 칸 하나만 두면, "
                     "참고기준이 슬그머니 실제 탈락 기준처럼 쓰이게 되고 50점이 저마다 "
                     "다른 잣대로 평가되기 쉽습니다. 어떤 기준의 몇 판을 보고 누가 "
                     "언제 판단했는지까지 남겨 두시면 나중에 그 판단을 되짚을 수 "
                     "있습니다.", warn=True)
    hdr = ["Site ID", "지점명", "구분",
           "4방위 측선 max|ΔF|\n(nT)", "4방위 측선 max−min\n(nT)",
           "1 m 지점 |변화율|\n(nT/m, 시범관측 대조용)",
           "전 거리 최대 |변화율|\n(nT/m)", "전 거리 중앙 |변화율|\n(nT/m)",
           "수직 ΔF 범위\n(nT)", "P0 전후차\n(nT)", "Set 지속 최대\n(분)",
           "시간변화 처리", "IAGA 50 nT\n대비", "유럽 3 nT/m\n대비",
           "예측 |∇ΔT|\n(참고)", "게이트 종합", "중심점 조정",
           "국내 판정상태", "기준명·버전", "판단자·일자", "판단 사유·특이사항"]
    table(ws, 4, hdr, [11, 14, 7, 15, 15, 18, 16, 16, 13, 12, 12, 20, 12, 12, 12,
                       12, 13, 15, 16, 14, 40], keys={"Site ID"})
    for i, p in enumerate(pts, 1):
        r, gap = 4 + i, p["예측구배"] is None
        for j, v in enumerate([p["_sid"], p["지점명"], p["구분"]], 1):
            cell(ws, r, j, v, fill=FILL_AUTO)
        for j in range(4, 12):
            cell(ws, r, j, fmt="0.00" if j in (6, 7, 8) else "0.0")
        for j in (12, 13, 14, 18, 19, 20, 21):
            cell(ws, r, j, al=AL_L if j == 21 else AL_C,
                 fmt="yyyy-mm-dd" if j == 20 else None)
        c = cell(ws, r, 15, "NOT AVAILABLE" if gap else round(p["예측구배"], 1),
                 fill=FILL_AUTO)
        if gap:
            c.font = F_WARN
        cell(ws, r, 16, f"='⑫ 게이트 점검'!N{r}", fill=FILL_AUTO)
        cell(ws, r, 17, f"='⑪ 중심점 최종선정'!F{r}", fill=FILL_AUTO)
    last = 4 + len(pts)
    dv(ws, f"L5:L{last}", ["Variometer 연속기록 + P0 보정 (선형성 검증됨)",
                           "P0 선형보간만 (선형성 미검증)", "미보정 — INVALID"])
    dv(ws, f"M5:N{last}", ["이내", "초과", "판정 불가", "PENDING"])
    # ⚠️ 「적합·부적합」을 쓰지 않는다 — 국내 기준이 없다(발주자 결정).
    dv(ws, f"R5:R{last}", ["연구판단 — 선점 후보 유지", "연구판단 — 조건부",
                           "연구판단 — 중심점 조정 후 재측정", "연구판단 — 재측정 필요",
                           "기준 수립용 자료로만 사용", "HOLD", "PENDING"])
    return ws


# ═══════════════════════════════════════════════ ⑭ 검토 반영 추적
def sheet_trace(wb, pts, counts):
    ws = wb.create_sheet("⑭ 검토 반영·작업량")
    for col, w in zip("ABCDEFG", [10, 34, 26, 34, 20, 20, 24]):
        ws.column_dimensions[col].width = w
    title(ws, "G", "Codex 독립 검토 반영 추적 · 현장 작업량 산정")
    note(ws, 2, "G", "수식에 오류가 없다는 것과 그 수식이 «맞는 행»을 가리킨다는 것은 다른 "
                     "이야기입니다. 엉뚱한 행을 참조해도 계산은 멀쩡히 됩니다. 그래서 "
                     "지적마다 어떤 사례로 확인했고 실제로 무엇이 나왔는지를 함께 "
                     "적어 두었습니다. 아래 ✔ 는 LibreOffice 로 다시 계산해 눈으로 "
                     "확인한 것입니다.", warn=True)
    r = 3
    ws.merge_cells(f"A{r}:G{r}")
    c = ws[f"A{r}"]
    c.value, c.font, c.fill, c.alignment = "■ 지적 반영 추적", F_SEC, FILL_SEC, AL_L
    r += 1
    table(ws, r, ["지적 ID", "지적 요지", "반영 위치(시트·열)", "검증 사례",
                  "기대 결과", "실제 결과", "잔여 위험"],
          [10, 34, 26, 34, 20, 20, 24])
    r += 1
    TRACE = [
        ("2차 C1", "Set 식별 불완전 · SUMIFS 가 누락·중복을 삼킴",
         "⑨·⑩ 「Set 검증」 열", "P0(후) 를 비우고 계산",
         "「P0(후) 유효 1개 아님」", "「P0(후) 유효 1개 아님」 · ΔF 공란 ✔",
         "행 삽입 시 Set ID 수기 입력 의존"),
        ("2차 C1", "측정시각이 P0 구간 밖", "⑨·⑩ 「Set 검증」",
         "측정시각을 P0(후) 이후로", "「측정시각 구간 밖」",
         "「측정시각 구간 밖」 ✔ · 시각 역전도 차단 확인 ✔", ""),
        ("2차 C2", "P0 선형성 미검증", "⑧ Variometer 파일 ID · ⑮ 시간변화 처리",
         "Variometer 미기록 세트", "「선형성 미검증」 선택 가능",
         "⑮ 드롭다운에 존재 ✔ · Set 지속 4분 자동산출 ✔", "임의 시간상한 두지 않음"),
        ("2차 C3", "ΔF/d 와 지점 대표구배 혼동",
         "⑨ 「중심점기준 변화율」·「|변화율|」·「인접구간 변화율」 분리",
         "P0보다 낮은 값(3 m 94 nT)", "ΔF·변화율 음수 보존",
         "ΔF −9.13 · 변화율 −3.04 · |변화율| 3.04 · 인접구간 −12.0 ✔",
         "대표값 확정은 국내기준 대기"),
        ("2차 C4", "4방향으로 반경 10 m 전체 판정 불가",
         "① 기입규약 · ⑨ 머리말 · ⑮ 열 이름", "—",
         "「4개 방위 측선의 관측점」 표기", "", "미관측 방위는 한계로 남음"),
        ("2차 C5", "옮긴 중심점 재검증 절차 없음",
         "⑪ 최초/최종 Location ID · 재측정 여부", "조정함 선택",
         "재측정 여부 필수 선택", "", "재측정 범위는 일정·비용 결정"),
        ("2차 C6", "D·I 와 F 의 연결 불완전", "해당 없음 — 절대측정 제외",
         "—", "D·I 관측행 자체가 없다", "종결 ✔", "절대관측 단계에서 다시 볼 것"),
        ("2차 M1", "수직 기준높이·산출물 미정의",
         "⑩ 실제 센서 중심높이 · 측정 기준면 · 인접높이 변화율", "—",
         "기술지표만 산출", "", "임계 없음 — 자동판정 안 함"),
        ("2차 M2", "작업량이 행 수로 과소 표현", "⑭ 작업량 산정표", "—",
         "판독 기준 재산정", "", "파일럿 실측 필요"),
        ("2차 M3", "1대 조건인데 교차배치 전제", "해당 없음 — 두지점 F 시트 제거",
         "—", "시트 자체가 없다", "종결 ✔", ""),
        ("2차 M4", "§13② 시설별 증빙 없음", "④ 시설 유형마다 한 행", "—",
         "4행/점 · 증빙 10열", "", ""),
        ("2차 M5", "미원 HOLD 범위 불명확", "⑤ 4단 플래그 · ② HOLD 범위", "—",
         "과거대조만 보류", "", "표석 동일성 확인 대기"),
        ("2차 M6", "중단·결측 시 세트 처리 규칙 없음", "⑨·⑩ Set 상태 · 중단 시각·사유",
         "—", "무효 세트 구분", "", ""),
        ("2차 M7", "판정에 근거·기준 버전 없음",
         "⑮ 관측사실/국제기준/국내상태/기준명·버전/판단자 분리", "—",
         "판정 재현 가능", "", "국내기준 미확정"),
        ("2차 m1", "예시에 거리 없어 1.75 재현 불가", "① 기입규약 — 거리 2 m 명시",
         "P0 100@10:00 · 106@10:02 · P0 105@10:04 · 2 m",
         "102.5 / 3.5 / 1.75 nT/m", "102.5 / 3.50 / 1.75 ✔ (원 계획 재현)", ""),
        ("2차 m2", "상태코드에 사유·이력 없음", "⑨⑩ 중단 시각·사유", "—",
         "사유 필수", "", ""),
        ("3차 C1", "카드 1판독 ↔ 사무실 3판독 불일치", "⑭ 작업량 · 현장카드 안내",
         "—", "1 판독 기준 통일 · 64회/점", "반영 ✔",
         "3회 수기 전환은 발주자 확인"),
        ("3차 C2", "Visit·Location·Set·Attempt 식별자 없음",
         "현장카드 머리 · 3·4·5구역", "—", "재방문·중심이동 구분 가능",
         "반영 ✔", "중심 이동 시 새 카드"),
        ("3차 C3", "§18③ 자기오염 방지 없음",
         "현장카드 「0. 측정 전 확인」 · ⑫ 게이트", "—",
         "개인 금속물·차량·전자기기 체크", "반영 ✔", ""),
        ("3차 M1", "절대측정 제외의 잔재", "규격 모듈 · ⑦ 장비 · ⑧ 관측조건",
         "—", "야간 의무 → 안정 구간 · DI-flux 범위 밖", "반영 ✔", ""),
        ("3차 M2", "「법정 요건」 명칭이 과장", "⑫ 제목·안내",
         "—", "「이번 단계 적용요소」 · 부분 적용 명시", "반영 ✔",
         "§15·§16·§18 적용범위는 확인 필요"),
        ("3차 M3", "자동지표만으로 재측정 판단 불가", "현장카드 완전성 점검 줄",
         "—", "필수행·시각·P0·Variometer 자가확인", "반영 ✔", ""),
        ("3차 M4", "방위표지 참방위각 재현 불가", "현장카드 5구역 확장",
         "—", "표지 ID·좌표·취득방법·기존값·결정방법", "반영 ✔", ""),
        ("3차 M5", "파일명만으로 사진 대응 불가", "현장카드 6구역",
         "—", "명명규칙 + 촬영시각", "반영 ✔", "보존·명명 표준 승인 필요"),
        ("3차 M6", "작업량을 확정치로 읽을 위험", "⑭ 작업량",
         "—", "「순수 현장 작업 가정치」 · 미반영 항목 명시", "반영 ✔",
         "파일럿 후 범위값 확정"),
        ("3차 M7", "다중 페이지 지점 식별 소실", "현장카드 인쇄설정",
         "—", "머리 3행 반복 · 머리글 Site ID · 쪽번호", "반영 ✔", ""),
        ("3차 m1", "삭제 전 시트번호 잔재", "① · ⑫ · ⑭",
         "—", "15시트 기준으로 갱신", "반영 ✔", ""),
    ]
    for t in TRACE:
        for j, v in enumerate(t, 1):
            cell(ws, r, j, v, fill=FILL_FIELD if j == 6 else FILL_AUTO,
                 al=AL_L if j in (2, 3, 4, 7) else AL_C)
        ws.row_dimensions[r].height = 30
        r += 1

    r += 1
    ws.merge_cells(f"A{r}:G{r}")
    c = ws[f"A{r}"]
    c.value = "■ 현장 작업량 — 「행 수」가 아니라 「판독 수」로 센다 (2차 M2)"
    c.font, c.fill, c.alignment = F_SEC, FILL_SEC, AL_L
    r += 1
    nH, nV = counts
    # ⚠️ 현장 카드가 위치당 1 판독을 받으므로 여기도 1 을 기준으로 센다
    #    (Codex Delta C1 — 두 산출물이 어긋나 있었다).
    per_read = 1
    rows_pt = nH + nV
    reads_pt = rows_pt * per_read
    table(ws, r, ["구분", "값", "산식", "비고", "", "", ""],
          [10, 34, 26, 34, 20, 20, 24])
    r += 1
    calc = [
        ("수평 관측행/점", nH, f"{len(SP.H_DIRECTIONS)}방향 × "
                             f"({len(SP.H_OFFSETS_M)}측점 + P0 전·후 2)", ""),
        ("수직 관측행/점", nV, f"{len(SP.V_HEIGHTS_CM)}높이 + 기준 전·후 2", ""),
        ("관측행/점", rows_pt, f"{nH} + {nV}", ""),
        ("F 판독/점", reads_pt, f"{rows_pt} × {per_read}회 판독",
         "⚠ 3회 수기로 바꾸면 192회 — 발주자 확인 사항"),
        ("F 판독 · 50점", reads_pt * len(pts), f"{reads_pt} × {len(pts)}", ""),
        ("점당 소요(가정)", "구배 1.6 h + 설치·철수 1 h ≈ 2.6 h",
         "행당 1.5분 가정", "⚠ 가정값 — 파일럿 실측 필요"),
        ("50점 총계(가정)", "약 130 h", "2.6 h × 50점",
         "**순수 현장 작업 가정치 — 이동·기상·재측정 제외**"),
        ("⚠ 반영 안 된 것", "측선 설치·실측거리 · GNSS 메타 · 방위표지 · 사진 · "
                          "Variometer 동기 · 중단·재측정·중심 이동",
         "파일럿에서 따로 재야 한다", "대표 지형별 «범위값»으로 제시할 것"),
        ("절대측정", "이 단계에 없음", "원 계획 31면 선점실시 절차",
         "선점이 끝난 뒤의 별도 관측 단계"),
    ]
    for k, v, f, n in calc:
        for j, val in enumerate([k, v, f, n], 1):
            cell(ws, r, j, val, fill=FILL_AUTO, al=AL_L if j in (3, 4) else AL_C)
        r += 1
    note(ws, r, "G", "⚠ 「점당 소요」와 「50점 총계」는 **가정값이다.** 행당 안정화·3회 판독·"
                     "이동시간을 대표 환경에서 파일럿으로 실측해 확정해야 한다. 또한 5~10 m "
                     "범위를 50점 모두 10 m 로 고정할지는 운영 결정 사항이다.", warn=True)
    return ws


# ═══════════════════════════════════════════════ ⑮ 확인 필요
def sheet_human(wb):
    ws = wb.create_sheet("⑮ 결정·확인 필요")
    for col, w in zip("ABCDE", [6, 34, 56, 20, 16]):
        ws.column_dimensions[col].width = w
    title(ws, "E", "발주자 결정 · 아직 확인이 필요한 것")
    note(ws, 2, "E", "Codex 독립 검토 3회가 「기관 확인 없이는 확정할 수 없다」고 지목한 "
                     "항목이다. **위 4건은 2026-09-07 에 답이 나왔고**, 아래는 아직 열려 "
                     "있다. 결정된 것과 열린 것을 섞지 않는 것이 요점이다.")
    table(ws, 3, ["#", "항목", "내용 · 왜 중요한가", "주체", "상태"],
          [6, 34, 56, 20, 16])
    r = 4
    # ── 결정된 것 ────────────────────────────────────────
    for i2, (k, v, why) in enumerate(SP.DECISIONS, 1):
        for j2, val in enumerate([f"D{i2}", k, f"{v}  {why}", "발주자", "확인 완료"], 1):
            cell(ws, r, j2, val, fill=FILL_AUTO,
                 al=AL_L if j2 in (2, 3) else AL_C,
                 font=F_WARN if j2 == 5 else F_VAL)
        ws.row_dimensions[r].height = 52
        r += 1
    # ── 아직 열린 것 ─────────────────────────────────────
    #
    # ⚠️ §19·§20(절대측정 요건)은 이 목록에서 빠졌다 — 절대측정 자체가 이 단계에
    #    없으므로 물을 것이 없다. 정식 관측 단계에서 다시 물어야 한다.
    OPEN = [
        ("이동식 Variometer 확보와 동시기록 가능 여부",
         "P0 전·후 두 점만으로는 그 사이의 비선형 변화·순간 교란을 «탐지할 수 없다». "
         "연속기록이 없으면 보정식을 적용해도 선형성이 검증되지 않는다.",
         "발주기관 · 예산"),
        ("참고값의 연산 정의 — 3 nT/m 대표지표 · 50 nT 계산법",
         "같은 자료도 최대·중앙·1 m·10 m·부호 처리에 따라 결과가 달라진다. 50 nT 도 "
         "max|ΔF| 인지 max−min 인지 정해야 한다. ⚠ 국내 «판정기준»은 이 탐사 결과로 "
         "만들 것이므로 지금 정하는 것은 «계산 정의»뿐이다.", "연구진"),
        ("최종 중심 이동 후 전체 재측정 여부",
         "옛 중심 자료는 새 중심의 10 m 범위와 수직구배를 입증하지 않는다. 재측정하면 "
         "일정·비용이 늘고, 안 하면 「잠정 중심」으로 남는다.", "발주기관 · 연구진"),
        ("미원의 표석 동일성 · 기존 16점의 성분별 대조 가능성",
         "물리점이 특정되지 않으면 과거 성과와 견줄 수 없다. 현재는 «과거대조만» "
         "보류로 두었고 신규 측정은 진행한다.", "국토지리정보원"),
        ("5~10 m 범위를 50점 모두 10 m 로 적용할지",
         "IAGA 예시가 「5~10 m」라 폭이 있다. 10 m 고정은 작업량을 늘린다. "
         "파일럿에서 함께 정하면 된다.", "연구진 · 운영"),
        ("상시관측소 원자료의 제공조건·시간해상도",
         "자료가 «존재한다»와 «쓸 수 있다»는 다르다. 시간변화 보정의 실현 가능성이 "
         "여기 달렸다.", "청양(INTERMAGNET) · KASA"),
        ("정식 선점·표지 설치로 넘어갈 시점과 그때의 적용 조문",
         "이번은 법정 선점이 아니다. 정식 단계로 갈 때 §14·§15·§16·§18 의 적용범위를 "
         "그때 판단해야 한다.", "국토지리정보원"),
    ]
    for i2, (item, why, who) in enumerate(OPEN, 1):
        for j2, val in enumerate([i2, item, why, who, ""], 1):
            cell(ws, r, j2, val, fill=FILL_FIELD if j2 == 5 else FILL_AUTO,
                 al=AL_L if j2 in (2, 3) else AL_C)
        ws.row_dimensions[r].height = 46
        r += 1
    dv(ws, f"E{4+len(SP.DECISIONS)}:E{r-1}", ["확인 완료", "질의 중", "PENDING"])
    return ws


# ═══════════════════════════════════════════════
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    pts = TP.load_points()
    print(f"대상 {len(pts)}점 로드")
    wb = Workbook()
    wb.remove(wb.active)
    sheet_intro(wb)
    sheet_points(wb, pts)
    sheet_site(wb, pts)
    sheet_clearance(wb, pts)
    sheet_existing(wb, pts)
    sheet_mark(wb, pts)
    sheet_instrument(wb)
    sheet_condition(wb)
    _, h_last = sheet_horizontal(wb, pts)
    sheet_vertical(wb, pts)
    sheet_center(wb, pts)
    sheet_gate(wb, pts)
    sheet_result(wb, pts)
    nH = len(SP.H_DIRECTIONS) * (len(SP.H_OFFSETS_M) + 2)
    nV = len(SP.V_HEIGHTS_CM) + 2
    sheet_trace(wb, pts, (nH, nV))
    sheet_human(wb)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"{dt.datetime.now():%Y%m%d_%H%M%S}_시험탐사_표준야장.xlsx"
    wb.save(path)
    print(f"시트 {len(wb.sheetnames)}장")
    print(f"  수평 {h_last-4}행 · 수직 {len(pts)*nV}행 · §13② {len(pts)*len(FACILITIES)}행")
    print(f"  판독 수(1회 기준): 점당 {nH+nV} · 전체 {(nH+nV)*len(pts):,}")
    print(f"[저장] {path}  ({path.stat().st_size/1e6:.1f} MB)")
    return path


if __name__ == "__main__":
    main()
