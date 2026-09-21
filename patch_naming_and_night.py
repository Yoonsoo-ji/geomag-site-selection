# -*- coding: utf-8 -*-
r"""이름을 값에 맞추고, 야간·자정 시각에 경고를 띄운다 (2026-09-22)

    python patch_naming_and_night.py --file "docs/output/…_오호입력_전체50.xlsx"

**① 이름** — 코덱스 M9 와 발주자 지시. 앞의 두 칸은 «구배»가 아니다.

| 종전 | 바꾼 이름 | 실제로 재는 것 |
|---|---|---|
| P0 전후차 | **P0 후 − 전** | 중심점을 처음과 끝에 두 번 재서 생긴 차이 — 재는 동안 자기장이 흔들린 정도 |
| 관측구간 F 변화폭 | **잰 구간 F 최대 − 최소** | 그 방향에서 나온 값의 폭(범위). 거리로 나누지 않았다 |

**② 야간·자정** — 야간에는 재지 않기로 했다(발주자 결정). 18시 이후나 06시 이전
블록이 있으면 확인 칸과 머리 경고줄에 띄운다. 자정을 넘기면 `hhmmss` 숫자가 작아져
블록 순위가 뒤집히므로, 한 자료 안 시각 폭이 12시간을 넘어도 경고한다.

⚠️ 자정 통과는 «고치지» 않고 «알린다». 파일 안 블록 순서가 시각순이 아닌 경우가
있어(오호 서쪽이 그랬다) 자정 넘김과 순서 뒤섞임을 자료만으로 구분할 수 없다.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

PASTE_TOP, PASTE_N = 9, 700
SUM_TOP, SUM_N = 9, 14
RES_TOP, PRE_TOP = 9, 29
NSET = 5
NOTE_ROW, NOTE_H = 37, 52.0            # ⚠️ 58 이면 1쪽이 넘쳐 인쇄가 5쪽이 된다
NL = chr(10)

HEAD_C = "P0 후 − 전" + NL + "(nT)"
# ⚠️ D 열 폭은 10(약 5 한글자) — 한 줄이 그보다 길면 어절 안에서 갈라진다
HEAD_D = "잰 구간 F" + NL + "최대−최소" + NL + "(nT)"
NOTE = (
    "구간 구배는 바로 앞에 잰 측점과의 F 차이를 거리 차로 나눈 값이고, P0(전)을 0 m 로 "
    "봅니다. 「P0 후 − 전」은 중심점을 처음과 끝에 두 번 재서 생긴 차이입니다 — 자리가 "
    "달라져서가 아니라 재는 동안 자기장이 흔들린 정도이니, 이 값이 크면 다른 칸에도 "
    "그만큼 섞여 있습니다. 「F 최대−최소」도 구배가 아니라 그 방향에서 나온 값의 "
    "폭입니다. 거리로 나눈 것은 구간 구배뿐입니다. 「잰 거리」는 마지막으로 값을 적은 "
    "거리이고, 「P0→끝」 두 칸은 그 거리까지의 차이와 평균입니다. 모두 시간변화를 빼기 "
    "전 값입니다. 국내 기준은 아직 없으니 값이 커도 현장에서 자리를 접지 마시고 나온 "
    "대로 적어 주세요.")


def L(n: int) -> str:
    s = ""
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def cols(d):
    b = 17 + 5 * d
    return dict(zip(("blk", "t", "f", "q", "use"), (L(b + i) for i in range(5))))


def sums(d):
    b = 43 + 9 * d
    k = ("k", "n", "t0", "t1", "med", "used", "avg", "drop", "rank")
    return dict(zip(k, (L(b + i) for i in range(len(k)))))


def t0_range(d: str) -> str:
    return f"${d}${SUM_TOP}:${d}${SUM_TOP + SUM_N - 1}"


def patch(f: Path) -> int:
    import pythoncom
    import win32com.client as win32

    pythoncom.CoInitialize()
    xl = win32.DispatchEx("Excel.Application")
    xl.Visible = False
    xl.DisplayAlerts = False
    xl.AskToUpdateLinks = False
    xl.ScreenUpdating = False
    cards = 0
    try:
        wb = xl.Workbooks.Open(str(f), UpdateLinks=0)
        xl.Calculation = -4135                          # 수동 — 중간 재계산 방지
        g = wb.Worksheets("⑨ 구배 자동계산")
        lo, hi = PASTE_TOP, PASTE_TOP + PASTE_N - 1

        # ── 미리보기 머리글 ────────────────────────────────────
        g.Range("I29").Value = "P0 후 − 전"
        g.Range("J29").Value = "F 최대−최소"

        # ── 확인 칸 — 야간·자정을 앞쪽 순위로 ──────────────────
        for d in range(4):
            r = PRE_TOP + 1 + d
            s, h = sums(d), cols(d)
            t0r = t0_range(s["t0"])
            nb = f"COUNT({t0r})"
            night = (f'COUNTIFS({t0r},">=180000")+'
                     f'COUNTIFS({t0r},">0",{t0r},"<60000")')
            span = f'AND({nb}>1,MAX({t0r})-MIN({t0r})>120000)'
            used = f'${s["used"]}${SUM_TOP}:${s["used"]}${SUM_TOP + SUM_N - 1}'
            drop = f'${s["drop"]}${SUM_TOP}:${s["drop"]}${SUM_TOP + SUM_N - 1}'
            thin = f'COUNTIFS({used},">0",{used},"<10")'
            noisy = f'SUMPRODUCT(({drop}<>"")*({used}>0)*({drop}>{used}))'
            orphan = (f'COUNTIFS(${h["blk"]}${lo}:${h["blk"]}${hi},0,'
                      f'${h["t"]}${lo}:${h["t"]}${hi},">0")')
            g.Cells(r, 15).Formula = (
                f'=IF({orphan}>0,"첫 /time 이 없습니다 — 머리글째 다시",'
                f'IF({nb}=0,"원문을 붙여 넣어 주세요",'
                f'IF({span},"시각이 12시간 넘게 벌어짐 — 자정 확인",'
                f'IF({night}>0,"야간 시각 "&{night}&"개 — 확인",'
                f'IF({nb}>12,"블록 "&{nb}&"개 — 10 m 까지만 표시",'
                f'IF({nb}<3,"블록이 "&{nb}&"개뿐 — P0 전·후가 다 있는지 보세요",'
                f'IF({thin}>0,{thin}&"개 지점 읽음 부족 — 확인",'
                f'IF({noisy}>0,{noisy}&"개 지점 제외 과다 — 확인","정상"))))))))')

        # ── 머리 경고줄 — 수직까지 다섯 자료를 한 번에 ─────────
        orphan_all = "+".join(
            f'COUNTIFS(${cols(d)["blk"]}${lo}:${cols(d)["blk"]}${hi},0,'
            f'${cols(d)["t"]}${lo}:${cols(d)["t"]}${hi},">0")' for d in range(NSET))
        night_all = "+".join(
            f'COUNTIFS({t0_range(sums(d)["t0"])},">=180000")+'
            f'COUNTIFS({t0_range(sums(d)["t0"])},">0",{t0_range(sums(d)["t0"])},"<60000")'
            for d in range(NSET))
        span_all = ",".join(
            f'AND(COUNT({t0_range(sums(d)["t0"])})>1,'
            f'MAX({t0_range(sums(d)["t0"])})-MIN({t0_range(sums(d)["t0"])})>120000)'
            for d in range(NSET))
        g.Range("A7").Formula = (
            f'=IF({orphan_all}>0,'
            '"⚠ 첫 /time 머리글이 빠진 자료가 있습니다 — 그 방향은 첫 측점이 버려지고 '
            '거리가 한 칸씩 당겨집니다. 파일을 머리글째 다시 붙여 넣으세요",'
            f'IF(OR({span_all}),'
            '"⚠ 한 자료 안에서 시각이 12시간 넘게 벌어졌습니다 — 자정을 넘겼다면 블록 '
            '순서가 뒤집힙니다. 시각을 확인하고 필요하면 날짜별로 나눠 붙여 넣으세요",'
            f'IF({night_all}>0,'
            '"⚠ 야간(18시 이후·06시 이전) 시각이 "&' + f'{night_all}' + '&"개 있습니다 — '
            '야간에는 재지 않기로 했으니 기기 시각이 맞는지 확인해 주세요","")))')

        # ── 카드 머리글·각주 ──────────────────────────────────
        for ws in wb.Worksheets:
            if ws.Range("A15").Text != "거리" or not ws.Range("A29").Text.startswith("②"):
                continue
            ws.Range("C30").Value = HEAD_C
            ws.Range("D30").Value = HEAD_D
            ws.Range(f"A{NOTE_ROW}").Value = NOTE
            ws.Rows(NOTE_ROW).RowHeight = NOTE_H
            cards += 1

        xl.Calculation = -4105
        xl.CalculateFull()
        g.Activate()
        g.Range("B4").Select()
        wb.Save()
        wb.Close(True)
    finally:
        xl.ScreenUpdating = True
        xl.Quit()
        pythoncom.CoUninitialize()
    return cards


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True)
    a = ap.parse_args()
    f = Path(a.file).resolve()
    bak = Path(__file__).parent / "_naming_backup.xlsx"
    shutil.copy(f, bak)
    print(f"카드 {patch(f)}장 · 되돌릴 사본 {bak.name}")


if __name__ == "__main__":
    main()
