# -*- coding: utf-8 -*-
r"""코덱스 검토 지적 반영 — 이미 만든 통합문서를 제자리에서 고친다 (2026-09-21)

    python fix_codex_findings.py --file "docs/output/…_오호입력_전체50.xlsx"

| 지적 | 실증 | 고침 |
|---|---|---|
| **C1** 블록 하나가 통째로 무효면 뒤 거리가 한 칸씩 당겨진다 | 재현됨 — 실제 3 m 가 「2 m」로 붙었다 | 블록 시각(`t0`·`t1`)을 «유효 판독»이 아니라 «그 블록의 모든 판독»에서 잡는다. 무효 블록도 제 자리를 지키고 F 만 빈다 |
| **C2** 첫 `/time` 이 없으면 첫 블록이 버려진다 | 재현됨 — 실제 1 m 가 「P0(전)」이 됐다 | 블록 0 에 판독이 있으면 확인 칸과 머리 경고줄에 띄운다 |
| **C3** 13블록 이상이면 11 m 가 조용히 사라진다 | 재현됨 | 발주자가 「11 m 이상은 안 잰다」로 확정 — 경고만 띄운다 |
| **M8** 품질 칸이 비면 유효로 통과 | `RIGHT("",1)<>"0"` 이 참 | 품질이 비면 제외 |
| **M6** 값이 하나뿐이어도 변화폭 0 | — | 두 개부터 계산, 아니면 빈칸 |
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

PASTE_TOP, PASTE_N = 9, 700
SUM_TOP, SUM_N = 9, 14
RES_TOP, PRE_TOP = 9, 29
TRIM = 5
DIRS = ["동", "서", "남", "북"]
NSET = 5


def L(n: int) -> str:
    s = ""
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def cols(d):                    # 줄 단위 숨긴 열
    b = 17 + 5 * d
    return dict(zip(("blk", "t", "f", "q", "use"), (L(b + i) for i in range(5))))


def sums(d):                    # 블록 요약 아홉 칸
    b = 43 + 9 * d
    k = ("k", "n", "t0", "t1", "med", "used", "avg", "drop", "rank")
    return dict(zip(k, (L(b + i) for i in range(len(k)))))


def patch(f: Path) -> None:
    import pythoncom
    import win32com.client as win32

    pythoncom.CoInitialize()
    xl = win32.DispatchEx("Excel.Application")
    xl.Visible = False
    xl.DisplayAlerts = False
    xl.AskToUpdateLinks = False
    xl.ScreenUpdating = False
    try:
        wb = xl.Workbooks.Open(str(f), UpdateLinks=0)
        # ⚠️ Calculation 은 통합문서가 열린 «뒤»에만 설정된다
        xl.Calculation = -4135                  # xlCalculationManual — 중간 재계산 방지
        g = wb.Worksheets("⑨ 구배 자동계산")
        lo, hi = PASTE_TOP, PASTE_TOP + PASTE_N - 1

        for d in range(NSET):
            h, s = cols(d), sums(d)
            B = f"${h['blk']}${lo}:${h['blk']}${hi}"
            T = f"${h['t']}${lo}:${h['t']}${hi}"
            # ① 시각은 블록의 «모든» 판독에서 — 무효 블록도 자리를 지킨다
            for i in range(SUM_N):
                r = SUM_TOP + i
                k = f"{s['k']}{r}"
                g.Range(f"{s['t0']}{r}").FormulaArray = f'=IFERROR(SMALL(IF({B}={k},{T},""),1),"")'
                g.Range(f"{s['t1']}{r}").FormulaArray = f'=IFERROR(LARGE(IF({B}={k},{T},""),1),"")'
            # ② 품질 칸이 비면 유효로 보지 않는다
            g.Range(f"{h['use']}{lo}:{h['use']}{hi}").Formula = (
                f'=IF(OR({h["t"]}{lo}="",{h["f"]}{lo}="",{h["q"]}{lo}=""),0,'
                f'IF(AND({h["f"]}{lo}>1000,RIGHT({h["q"]}{lo},1)<>"0"),1,0))')

        # ③ 미리보기 — 변화폭 하한 · 확인 칸 경고 두 가지
        for d in range(4):
            r = PRE_TOP + 1 + d
            s, h = sums(d), cols(d)
            fc = L(9 + 2 * d)
            nb = f'COUNT(${s["t0"]}${SUM_TOP}:${s["t0"]}${SUM_TOP + SUM_N - 1})'
            g.Cells(r, 10).Formula = (
                f'=IF(COUNT(${fc}${RES_TOP}:${fc}${RES_TOP + 10})<2,"",'
                f'MAX(${fc}${RES_TOP}:${fc}${RES_TOP + 10})'
                f'-MIN(${fc}${RES_TOP}:${fc}${RES_TOP + 10}))')
            used = f'${s["used"]}${SUM_TOP}:${s["used"]}${SUM_TOP + SUM_N - 1}'
            drop = f'${s["drop"]}${SUM_TOP}:${s["drop"]}${SUM_TOP + SUM_N - 1}'
            thin = f'COUNTIFS({used},">0",{used},"<10")'
            noisy = f'SUMPRODUCT(({drop}<>"")*({used}>0)*({drop}>{used}))'
            orphan = (f'COUNTIFS(${h["blk"]}${lo}:${h["blk"]}${hi},0,'
                      f'${h["t"]}${lo}:${h["t"]}${hi},">0")')
            g.Cells(r, 15).Formula = (
                f'=IF({orphan}>0,"첫 /time 이 없습니다 — 머리글째 다시",'
                f'IF({nb}=0,"원문을 붙여 넣어 주세요",'
                f'IF({nb}>12,"블록 "&{nb}&"개 — 10 m 까지만 표시",'
                f'IF({nb}<3,"블록이 "&{nb}&"개뿐 — P0 전·후가 다 있는지 보세요",'
                f'IF({thin}>0,{thin}&"개 지점 읽음 부족 — 확인",'
                f'IF({noisy}>0,{noisy}&"개 지점 제외 과다 — 확인","정상"))))))')

        # ④ 머리 경고줄 — 수직까지 다섯 자료를 한 번에 본다
        terms = "+".join(
            f'COUNTIFS(${cols(d)["blk"]}${lo}:${cols(d)["blk"]}${hi},0,'
            f'${cols(d)["t"]}${lo}:${cols(d)["t"]}${hi},">0")' for d in range(NSET))
        g.Range("A7:O7").Merge()
        g.Range("A7").Formula = (
            f'=IF({terms}=0,"",'
            '"⚠ 첫 /time 머리글이 빠진 자료가 있습니다 — 그 방향은 첫 측점이 버려지고 '
            '거리가 한 칸씩 당겨집니다. 파일을 머리글째 다시 붙여 넣으세요")')
        g.Range("A7").Font.Bold = True
        g.Range("A7").Font.Color = 0x2020C0
        g.Rows(7).RowHeight = 18

        # ⑤ 카드 ② 변화폭 — 값이 하나뿐이면 빈칸
        cards = 0
        for ws in wb.Worksheets:
            if ws.Range("A15").Text != "거리" or not ws.Range("A29").Text.startswith("②"):
                continue
            for r, c in zip(range(31, 35), ("C", "E", "G", "I")):
                ws.Range(f"D{r}").Formula = (
                    f'=IF(COUNT({c}16:{c}26)<2,"",MAX({c}16:{c}26)-MIN({c}16:{c}26))')
            cards += 1
        xl.Calculation = -4105                  # xlCalculationAutomatic
        xl.CalculateFull()
        wb.Worksheets("⑨ 구배 자동계산").Activate()
        g.Range("B4").Select()
        wb.Save()
        wb.Close(True)
        print(f"카드 {cards}장 · 자료 {NSET}종 수식 교체")
    finally:
        xl.ScreenUpdating = True
        xl.Quit()
        pythoncom.CoUninitialize()


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True)
    a = ap.parse_args()
    f = Path(a.file).resolve()
    shutil.copy(f, Path(__file__).parent / "docs" / "output" / ("_codexfix_backup_" + f.name))
    patch(f)


if __name__ == "__main__":
    main()
