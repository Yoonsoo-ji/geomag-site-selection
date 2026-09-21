# -*- coding: utf-8 -*-
r"""⑨ 구배 자동계산 결과를 카드가 «직접» 읽게 한다 (2026-09-21)

    python link_cards_to_gradient.py --into "docs/output/…_구배자동계산_전체50.xlsx"

종전에는 결과표를 복사해 카드에 값으로 붙여 넣어야 했다. 이제 카드의 입력 칸이
⑨ 시트를 가리키는 수식이라, **⑨ 에서 그 지점을 고르는 순간 카드가 채워진다.**

    B16:I27  ← ⑨ H9:O20   (네 방향 시각·평균 F)
    C42:C45  ← ⑨ I23:I26   (수직 봉 4·3·2·1 평균 F)

⚠️ 수식이므로 **다음 지점 자료를 붙여 넣으면 앞 지점 카드는 비워진다.** 한 지점을
끝냈으면 그 카드에서 두 범위를 복사해 «값 붙여넣기»(제자리) 로 고정한다. 손으로 값을
쳐 넣어도 그 자리의 수식이 지워지므로 결과는 같다 — 즉 «고정»과 «수기 입력»이 같은
동작이고, 종이 카드에는 아무 영향이 없다(수식이 빈 문자열을 내므로 인쇄는 그대로 빈 칸).

또 ② 결과표의 방향 이름에 «그 방향에서 마지막으로 잰 거리»를 붙인다 — 동 6 m · 남 4 m
처럼 방향마다 간 거리가 다르므로, 「P0→끝」이 무엇의 끝인지 표에서 바로 읽히게 한다.
"""
from __future__ import annotations

import argparse
import datetime as dt
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "docs" / "output"
G = "'⑨ 구배 자동계산'"
GS = "⑨ 구배 자동계산"
DIRS = {31: "C", 32: "E", 33: "G", 34: "I"}      # 동·서·남·북 결과행 ↔ 평균 F 열
NAMES = {31: "동", 32: "서", 33: "남", 34: "북"}
P0_ROW, FIRST_M, LAST_M = 16, 17, 26

STEP3 = ("③ 위에서 고른 지점 카드에 결과가 «저절로» 들어갑니다 — 옮겨 적지 않으셔도 됩니다.   "
         "④ 다만 다음 지점 자료를 붙여 넣으면 앞 지점 카드는 비워집니다. 한 지점을 끝냈으면 "
         "그 카드의 B16:I27 과 C42:C45 를 복사해 그 자리에 «값 붙여넣기» 로 고정해 주세요.   "
         "⑤ 다음 지점을 하려면 노란 칸을 모두 지우고 다시 붙여 넣으세요.")
GUIDE4 = ('=IF($B$4="","지점을 고르면 그 카드에 결과가 저절로 들어갑니다 — 먼저 골라 주세요",'
          '"「"&$B$4&"」 카드에 결과가 자동으로 들어갑니다  ·  다음 지점 전에 그 카드의 '
          'B16:I27 · C42:C45 를 «값 붙여넣기» 로 고정하세요")')


def card_formulas(name: str) -> dict:
    """카드 한 장에 넣을 수식. 지점이 골라져 있을 때만 값을 낸다."""
    gate = f'{G}!$B$4<>"{name}"'
    out = {
        "B16:I27": f'=IF({gate},"",IF({G}!H9="","",{G}!H9))',
        "C42:C45": f'=IF({gate},"",IF({G}!I23="","",{G}!I23))',
        "A14": ('="① 방향별 관측자료"&IF(' + f'{G}!$B$4="{name}"'
                + ',"    ←  ⑨ 시트 연결 중 · 다음 지점 전에 «값 붙여넣기» 로 고정","")'),
    }
    for r, col in DIRS.items():
        rng = f"{col}{FIRST_M}:{col}{LAST_M}"
        out[f"A{r}"] = (f'=IF(COUNT({rng})=0,"{NAMES[r]}","{NAMES[r]} · 끝 "&'
                        f'LOOKUP(2,1/({rng}<>""),ROW({rng})-{P0_ROW})&" m")')
    return out


def link(target: Path, out: Path) -> list[str]:
    import pythoncom
    import win32com.client as win32

    shutil.copy(target, out)
    pythoncom.CoInitialize()
    xl = win32.DispatchEx("Excel.Application")
    xl.Visible = False
    xl.DisplayAlerts = False
    xl.AskToUpdateLinks = False
    done: list[str] = []
    try:
        wb = xl.Workbooks.Open(str(out), UpdateLinks=0)
        g = wb.Worksheets(GS)
        a2 = str(g.Range("A2").Value or "")
        head = a2.split("③")[0].rstrip()
        g.Range("A2").Value = head + "   " + STEP3
        g.Range("C4").Formula = GUIDE4
        for ws in wb.Worksheets:
            if ws.Range("A15").Text != "거리" or not ws.Range("A29").Text.startswith("②"):
                continue
            for addr, f in card_formulas(ws.Name).items():
                ws.Range(addr).Formula = f
            done.append(ws.Name)
        wb.Worksheets(3).Activate()
        g.Range("B4").Select()
        wb.Save()
        wb.Close(True)
    finally:
        xl.Quit()
        pythoncom.CoUninitialize()
    return done


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--into", required=True)
    ap.add_argument("--tag", default="카드자동연결_전체50")
    a = ap.parse_args()
    out = OUT_DIR / f"{dt.datetime.now():%Y%m%d_%H%M%S}_시험탐사_현장야장_3판_{a.tag}.xlsx"
    done = link(Path(a.into).resolve(), out)
    print(f"연결한 카드 {len(done)}장: {done[:3]} …")
    print(f"[저장] {out}  ({out.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
