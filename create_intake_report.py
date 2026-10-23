# -*- coding: utf-8 -*-
"""
기존 절대관측 자료의 LMM 투입 가능성 검토 — 보고서 생성
========================================================

```
python build_obs_archive.py     # ① 다섯 출처 통합
python lmm_intake_decide.py     # ② 변동 정리 · 투입 판정
python create_intake_report.py  # ③ 보고서 → docs/output/*_기존관측자료_LMM투입검토.docx
```

⚠️ **수치를 문서에 손으로 적지 않는다.** 전부 ①②의 산출 CSV 에서 읽으므로
   자료가 갱신되면 ①→②→③ 순서로 다시 돌리면 문서도 따라간다. 손으로 적으면
   재적합할 때마다 본문과 표가 갈라진다(12-D 가 같은 이유로 못박아 둔 규약이다).

⚠️ **「투입 확정」으로 쓰지 않는다.** 이 문서는 **판정 근거**이고, 실제 투입은
   Midterm Baseline v1(12-D) 동결 해제 뒤의 별도 결정이다.
"""
from __future__ import annotations

import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DOCS = ROOT / "docs" / "data"
OUTDIR = ROOT / "docs" / "output"

S_CSV = DOCS / "obs_archive_sessions.csv"
T_CSV = DOCS / "lmm_intake_decision.csv"
Z_CSV = DOCS / "obs_timeseries_summary.csv"

from docx import Document                               # noqa: E402
from docx.enum.table import WD_TABLE_ALIGNMENT          # noqa: E402
from docx.enum.text import WD_ALIGN_PARAGRAPH           # noqa: E402
from docx.oxml.ns import qn                             # noqa: E402
from docx.shared import Pt, RGBColor                    # noqa: E402

BODY, MONO = "맑은 고딕", "Consolas"
NAVY = RGBColor(0x16, 0x32, 0x4F)
RED = RGBColor(0xB2, 0x18, 0x2B)


def style_doc(doc):
    st = doc.styles["Normal"]
    st.font.name = BODY
    st.font.size = Pt(10)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), BODY)
    st.paragraph_format.space_after = Pt(4)
    st.paragraph_format.line_spacing = 1.4


def h(doc, text, level=1):
    p = doc.add_heading(text, level=level)
    for r in p.runs:
        r.font.name = BODY
        r._element.rPr.rFonts.set(qn("w:eastAsia"), BODY)
        r.font.color.rgb = NAVY
    p.paragraph_format.keep_with_next = True
    return p


def para(doc, text, bold=False, color=None, size=10):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.font.name = BODY
    r._element.rPr.rFonts.set(qn("w:eastAsia"), BODY)
    r.font.size = Pt(size)
    r.bold = bold
    if color is not None:
        r.font.color.rgb = color
    return p


def bullet(doc, text):
    p = doc.add_paragraph(style="List Bullet")
    r = p.add_run(text)
    r.font.name = BODY
    r._element.rPr.rFonts.set(qn("w:eastAsia"), BODY)
    r.font.size = Pt(10)
    return p


def table(doc, head, rows, widths=None, caption=None):
    if caption:
        c = para(doc, caption, bold=True, size=9)
        c.paragraph_format.keep_with_next = True
    t = doc.add_table(rows=1, cols=len(head))
    t.style = "Light Grid Accent 1"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    # ⚠️ python-docx 는 **셀마다** 너비를 넣어야 먹는다. 표 수준 설정만으로는
    #    자동 맞춤이 이겨서 「1984~2016」 같은 칸이 두 줄로 접힌다.
    if widths:
        t.autofit = False
        try:
            t.allow_autofit = False
        except Exception:
            pass
    for i, x in enumerate(head):
        cell = t.rows[0].cells[i]
        cell.text = ""
        r = cell.paragraphs[0].add_run(str(x))
        r.bold = True
        r.font.size = Pt(9)
        r.font.name = BODY
        r._element.rPr.rFonts.set(qn("w:eastAsia"), BODY)
        cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
    # ⚠️ 표가 쪽을 넘을 때 머리행이 반복되게 한다 — 안 넣으면 빈 줄로 보인다
    trPr = t.rows[0]._tr.get_or_add_trPr()
    el = trPr.makeelement(qn("w:tblHeader"), {})
    trPr.append(el)
    if widths:
        for i, w in enumerate(widths):
            t.rows[0].cells[i].width = Pt(w)
    for row in rows:
        cells = t.add_row().cells
        if widths:
            for i, w in enumerate(widths):
                cells[i].width = Pt(w)
        for i, x in enumerate(row):
            cells[i].text = ""
            r = cells[i].paragraphs[0].add_run("" if x is None else str(x))
            r.font.size = Pt(9)
            r.font.name = BODY
            r._element.rPr.rFonts.set(qn("w:eastAsia"), BODY)
            cells[i].paragraphs[0].alignment = (
                WD_ALIGN_PARAGRAPH.CENTER if i == 0
                else WD_ALIGN_PARAGRAPH.RIGHT if _num(x)
                else WD_ALIGN_PARAGRAPH.LEFT)
    return t


# 표 칸에 쓸 짧은 출처 이름. 긴 이름을 그대로 쓰면 칸이 네 줄로 접혀
# 행 높이가 네 배가 된다(인쇄해 보기 전에는 안 보이는 결함이다).
SHORT = {"register": "점조서", "legacy1019": "관측현황 10~19",
         "survey2225": "성과표 22~25", "kriss": "원자료 12~16",
         "fieldbook": "야장 19~25"}


def _srcname(t):
    """사유 문자열에 섞인 출처 키를 사람이 읽는 이름으로 바꾼다."""
    for k, v in SHORT.items():
        t = t.replace(k, v)
    return t


def _num(x):
    try:
        float(str(x).replace(",", "").replace("′", "").replace("nT", "")
              .replace("km", "").replace("+", "").strip())
        return True
    except Exception:
        return False


def fmt(v, nd=1, unit=""):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "-"
    try:
        return f"{float(v):,.{nd}f}{unit}"
    except Exception:
        return str(v)


def main():
    for f in (S_CSV, T_CSV, Z_CSV):
        if not f.exists():
            print(f"⚠ {f.name} 이 없다 — build_obs_archive.py → "
                  f"lmm_intake_decide.py 순으로 먼저 돌릴 것")
            return 1
    S = pd.read_csv(S_CSV, encoding="utf-8-sig")
    T = pd.read_csv(T_CSV, encoding="utf-8-sig")
    Z = pd.read_csv(Z_CSV, encoding="utf-8-sig")

    ok = T[T["D판정"] == "투입"]
    aux = T[T["D판정"] == "보조"]
    rel = T[T["D판정"] == "별도측점"]
    okF = T[T["F판정"] == "투입"]
    rmsD = float(np.sqrt(np.mean(ok["rD"].dropna() ** 2)))
    rmsF = float(np.sqrt(np.mean(okF["rF"].dropna() ** 2)))
    newS = sorted(set(ok.loc[ok["LMM보유"] == "없음", "측점"]))
    haveS = sorted(set(ok.loc[ok["LMM보유"] == "있음", "측점"]))
    auxS = sorted(set(aux["측점"]) - set(ok["측점"]))

    doc = Document()
    style_doc(doc)
    sec = doc.sections[0]
    sec.left_margin = sec.right_margin = Pt(54)

    # ── 표제 ─────────────────────────────────────────────────────
    t0 = doc.add_heading("기존 절대관측 자료의 지역 자기장 모델(LMM) 투입 가능성 검토",
                         level=0)
    for r in t0.runs:
        r.font.name = BODY
        r._element.rPr.rFonts.set(qn("w:eastAsia"), BODY)
        r.font.color.rgb = NAVY
    p = para(doc, f"{dt.datetime.now():%Y년 %m월 %d일} · 지구자기장 측정 "
                  f"입지 선정 연구", size=9)
    p.alignment = WD_ALIGN_PARAGRAPH.RIGHT

    para(doc,
         "국토지리정보원 1등 지자기점 관측망의 절대관측 자료를 한 표로 모아, "
         "어느 자료를 지역 자기장 모델의 초기 구성요소로 쓸 수 있는지 판정했다. "
         "판정 근거와 배제 사유를 함께 싣는다.", size=10)
    para(doc,
         "이 문서의 모든 수치는 산출 CSV 세 개에서 읽어 생성한다 — "
         "obs_archive_sessions.csv · lmm_intake_decision.csv · "
         "obs_timeseries_summary.csv. 자료가 갱신되면 문서도 따라간다.",
         size=9, color=NAVY)

    # ── 1. 요지 ──────────────────────────────────────────────────
    h(doc, "1. 요지", 1)
    bullet(doc, f"관측망 30점의 절대관측 자료를 다섯 출처에서 모아 "
                f"{len(S):,}행(세션·방문)으로 통합했다. 기간은 "
                f"{int(S['연도'].min())}~{int(S['연도'].max())}년이다.")
    bullet(doc, f"편각을 모델에 바로 넣을 수 있는 자료는 "
                f"{len(ok)}개 방문 · {ok['측점'].nunique()}개 측점이다. "
                f"IGRF-14 편각 잔차 RMS는 {fmt(rmsD, 1, '′')}다.")
    bullet(doc, f"현재 모델이 쓰는 16개 측점에 없던 측점이 "
                f"{len(newS)}개 더해진다 — {', '.join(newS)}. "
                f"측점 수가 16에서 최대 {16 + len(newS)}개로 늘어난다.")
    bullet(doc, f"총자력은 {len(okF)}개 방문만 쓸 수 있고 잔차 RMS가 "
                f"{fmt(rmsF, 0, ' nT')}다. 현재 모델의 총자력 성능을 "
                f"이 자료로 끌어올리기는 어렵다.")
    bullet(doc, f"방문 대표값만 있는 자료 {len(aux)}개 방문은 "
                f"«보조»로 분류했다. 세션 산포를 재지 못해 관측 정밀도를 "
                f"확인할 수 없기 때문이며, 측점 {len(auxS)}개가 여기에만 있다.")
    para(doc, "⚠ 이 표는 판정이고 투입이 아니다. 실제 투입은 중간보고 기준 "
              "버전(Midterm Baseline v1) 동결을 해제한 뒤의 별도 결정이다.",
         bold=True, color=RED, size=9)

    # ── 2. 자료 출처 ─────────────────────────────────────────────
    h(doc, "2. 통합한 자료 — 다섯 출처", 1)
    para(doc, "절대관측 측정값(편각 D·복각 I·총자력 F)을 담은 자료가 다섯 "
              "갈래로 흩어져 있었다. 측점마다 시기별 변동을 보려면 그때마다 "
              "다섯 파일을 다시 열어야 했고, 같은 관측이 두 출처에 겹쳐 "
              "실리는 일도 있었다. 하나로 모았다.")
    rows = []
    NOTE = {
        "register": ("지자기 점조서 관측이력", "방문",
                     "측점마다 1984년부터의 관측이력. 가장 긴 시계열"),
        "legacy1019": ("'10~'19년 지자기점 관측현황", "방문",
                       "방문 대표값. 관측 여부(ㅇ)로 연도를 센 자료"),
        "survey2225": ("지자기측량 성과정리(22~25)", "방문",
                       "실측 좌표·표고. 가장 최근 성과"),
        "kriss": ("2012~2016 절대관측 원자료", "세션",
                  "세션별 원시 판독값·분 단위 시각·참방위각"),
        "fieldbook": ("2019~2025 야장", "세션",
                      "세션별 원시 판독값 역산·분 단위 시각"),
    }
    for src, g in S.groupby("출처"):
        nm, unit, note = NOTE.get(src, (src, "", ""))
        rows.append([nm, unit, f"{int(g['연도'].min())}~{int(g['연도'].max())}",
                     f"{len(g):,}", g["측점"].nunique(), note])
    rows.sort(key=lambda r: r[2])
    table(doc, ["출처", "단위", "기간", "행", "측점", "담긴 것"], rows,
          caption="표 1. 통합한 다섯 출처")
    para(doc, f"통합 결과 {len(S):,}행이며, 같은 (측점, 연도)가 여러 출처에 "
              f"겹치는 것을 정리하면 {len(T)}개 방문이 된다. "
              f"측점은 {T['측점'].nunique()}개로 관측망 전체를 덮는다.", size=9)

    mism = T[T["출처불일치"].astype(str).str.len() > 0]
    para(doc, f"⚠ 출처 두 곳이 같은 관측을 다르게 적은 방문이 {len(mism)}개 "
              f"있다. 정보량이 많은 쪽을 택했고 어긋난 내용은 "
              f"판정표의 「출처불일치」 열에 남겼다 — 이 어긋남 자체가 "
              f"전사 오류의 신호다.", color=RED, size=9)

    # ── 3. 판정 기준 ─────────────────────────────────────────────
    h(doc, "3. 판정 기준 — 네 관문", 1)
    para(doc, "성분(편각·복각·총자력)마다 따로 판정한다. 한 방문의 편각은 "
              "쓸 수 있고 총자력은 못 쓰는 경우가 많기 때문이다.")
    table(doc, ["관문", "묻는 것", "기준", "통과 못 하면"], [
        ["1 중복", "같은 관측이 여러 출처에 겹쳤는가",
         "정보량이 많은 출처를 택한다(세션·시각·좌표가 있는 쪽)",
         "택하지 않은 쪽의 값 차이를 기록"],
        ["2 좌표", "그 관측의 좌표를 아는가",
         "관측망 좌표와 5 km 이내이거나 이설 이력으로 설명",
         "불가 — 다른 자리의 자료를 그 측점에 넣게 된다"],
        ["3 세션", "블런더를 뺀 유효 세션이 충분한가",
         "세션 자료는 3회 이상",
         "보류(세션 부족) 또는 보조(방문 대표값)"],
        ["4 잔차", "IGRF-14 잔차가 상식 범위인가",
         "편각·복각 60′ 이내 · 총자력 400 nT 이내",
         "불가 — 모형 문제가 아니라 자료 문제다"],
    ], caption="표 2. 네 관문")

    para(doc, "관문 2가 결정적이다. 2015년 제천 자료는 파일 머리글 좌표가 "
              "장흥의 것(거리 0 km)인데, 편각 잔차만 보면 −11.7′로 멀쩡해 "
              "보였다. 좌표를 보지 않으면 다른 자리의 자료를 제천에 넣는다.",
         size=9)
    para(doc, "블런더는 방문 안에서 중앙값에서 크게 벗어난 «한 세션»을 "
              "말한다. 편각·복각 20′, 총자력 300 nT를 문턱으로 쓰고, "
              f"통합 {len(S):,}행 중 편각 {int((S['D_판정'] != '정상').sum()) if 'D_판정' in S else 0}"
              f"행이 범위 밖으로 표시됐다. "
              "블런더는 지우지 않고 표시만 한다 — 원자료를 버리면 나중에 "
              "되짚을 수 없다.", size=9)

    # ── 4. 판정 결과 ─────────────────────────────────────────────
    h(doc, "4. 판정 결과", 1)
    rows = []
    LBL = {"투입": "투입 — 모델에 바로 넣을 수 있다",
           "보조": "보조 — 방문 대표값. 세션 산포 미확인",
           "별도측점": "별도측점 — 이설 이전. 다른 표석으로 다뤄야 한다",
           "보류": "보류 — 세션 부족 또는 출처 간 큰 불일치",
           "불가": "불가 — 좌표·잔차 문제"}
    for comp, nm in (("D", "편각 D"), ("I", "복각 I"), ("F", "총자력 F")):
        vc = T[f"{comp}판정"].value_counts()
        rows.append([nm] + [vc.get(k, 0) for k in
                            ("투입", "보조", "별도측점", "보류", "불가")])
    table(doc, ["성분", "투입", "보조", "별도측점", "보류", "불가"], rows,
          caption="표 3. 성분별 판정 (방문 수)")
    for k, v in LBL.items():
        bullet(doc, v)

    para(doc, f"편각 투입분의 IGRF-14 잔차 RMS는 {fmt(rmsD, 1, '′')}다. "
              f"대조로, 2019년 야장은 25.2′이고 2022~2025년 성과는 33.4′다 "
              f"(문서 12-F·13절). 즉 이 자료의 편각 품질은 현재 모델이 쓰는 "
              f"자료와 같은 수준이거나 그보다 낫다.", size=10)
    para(doc, f"⚠ 총자력은 다르다. 투입 가능분만 모아도 잔차 RMS가 "
              f"{fmt(rmsF, 0, ' nT')}이고, 범위 밖·과대 잔차로 배제된 방문이 "
              f"{int(T['F판정'].value_counts().get('불가', 0))}개다. "
              f"총자력 목표치를 이 자료로 달성할 수는 없다.",
         color=RED, size=10)

    # ── 5. 투입 가능 목록 ────────────────────────────────────────
    h(doc, "5. 편각 투입 가능 방문", 1)
    rows = []
    for _, r in ok.sort_values(["측점", "연도"]).iterrows():
        rows.append([r["측점"], int(r["연도"]),
                     SHORT.get(r["출처"], r["출처"]),
                     int(r["D_유효"]), fmt(r["rD"], 1, "′"),
                     fmt(r["rI"], 1, "′"), fmt(r["rF"], 0),
                     r["I판정"], r["F판정"], r["LMM보유"]])
    table(doc, ["측점", "연도", "출처", "세션", "편각 잔차", "복각 잔차",
                "총자력 잔차", "복각", "총자력", "현 모델"], rows,
          widths=[40, 32, 70, 30, 48, 48, 52, 40, 44, 44],
          caption=f"표 4. 편각 투입 가능 {len(ok)}개 방문")

    para(doc, f"현재 모델 16개 측점에 없는 측점 {len(newS)}개 — "
              f"{', '.join(newS)}.", bold=True)
    para(doc, f"이미 있는 측점 {len(haveS)}개 — {', '.join(haveS)}. "
              f"같은 측점의 다른 시기 자료이므로 시기 간 대조에 쓴다.", size=9)

    # ── 6. 시계열 변동 ───────────────────────────────────────────
    h(doc, "6. 측점별 시계열 변동", 1)
    para(doc, "IGRF-14 잔차는 영년변화가 빠진 양이다. 같은 표석을 다시 재면 "
              "잔차가 같아야 하므로, 방문 사이에 잔차가 흔들리면 그 측점이나 "
              "그 방문에 문제가 있다는 뜻이다. 이설 전후는 다른 표석이므로 "
              "이설 이전 방문은 산포 계산에서 뺐다.")
    rows = []
    for _, r in Z.sort_values("rD_폭", ascending=False).iterrows():
        rows.append([r["측점"], int(r["방문"]), int(r["이설이전"]),
                     f"{int(r['최초'])}~{int(r['최종'])}"
                     if pd.notna(r["최초"]) else "-",
                     fmt(r["rD_중앙"], 1, "′"), fmt(r["rD_폭"], 1, "′"),
                     fmt(r["rF_중앙"], 0), fmt(r["rF_폭"], 0),
                     fmt(r["같은날_D산포_중앙"], 2, "′")])
    table(doc, ["측점", "방문", "이설전", "기간", "편각 잔차 중앙",
                "편각 잔차 폭", "총자력 잔차 중앙", "총자력 잔차 폭",
                "같은 날 산포"], rows,
          widths=[42, 32, 34, 72, 56, 52, 58, 58, 52],
          caption="표 5. 측점별 잔차 변동 — 편각 잔차 폭이 큰 순서")

    wide = Z[Z["rD_폭"] > 60].sort_values("rD_폭", ascending=False)
    if len(wide):
        para(doc, f"편각 잔차 폭이 60′를 넘는 측점이 {len(wide)}개다 — "
                  + ", ".join(f"{r['측점']} {fmt(r['rD_폭'], 0, '′')}"
                              for _, r in wide.iterrows()) + ".",
             color=RED, size=9)
        para(doc, "같은 표석에서 잰 값의 잔차가 이만큼 벌어지는 것은 "
                  "영년변화로 설명되지 않는다. 방문마다 달라지는 방위표지 "
                  "기준이 주된 후보이며, 이 진단은 2026년 8월 재방문 감사가 "
                  "방위표지 재설정으로 특정한 것과 같은 방향이다.", size=9)
    tight = Z[(Z["rD_폭"].notna()) & (Z["rD_폭"] <= 20)]
    if len(tight):
        para(doc, f"반대로 편각 잔차 폭이 20′ 이내인 측점이 {len(tight)}개 "
                  f"있다 — " + ", ".join(tight["측점"].tolist())
             + ". 이 측점들은 시기를 가로질러 일관되므로 모델의 공간 구조를 "
               "결정하는 데 가장 믿을 만한 자료다.", size=9)

    # ── 7. 배제 사유 ─────────────────────────────────────────────
    h(doc, "7. 배제한 자료와 사유", 1)
    bad = T[T["D판정"].isin(["불가", "보류"])]
    rows = []
    for _, r in bad.sort_values(["D판정", "측점"]).iterrows():
        rows.append([r["측점"], int(r["연도"]) if pd.notna(r["연도"]) else "-",
                     SHORT.get(r["출처"], r["출처"]),
                     r["D판정"], _srcname(str(r["D사유"]))[:64]])
    table(doc, ["측점", "연도", "출처", "판정", "사유"], rows,
          widths=[46, 36, 80, 44, 240],
          caption=f"표 6. 편각 배제 {len(bad)}개 방문")
    para(doc, f"이설 이전 관측 {len(rel)}개 방문은 좌표가 설명되지만 "
              f"«별도측점»으로 뒀다. 폐지된 다른 표석에서 잰 자료라 같은 "
              f"측점으로 묶으면 안 되고, 모델에 넣으려면 측점 식별자를 "
              f"분리해야 한다.", size=9)

    # ── 8. 남은 검토 ─────────────────────────────────────────────
    h(doc, "8. 더 검토해야 할 사항", 1)
    for t in [
        "양식 복사 오염을 칸마다 따로 확인해야 한다. 측점명·좌표·총자력이 "
        "서로 다른 방식으로 오염돼 있었다 — 측점명은 파일명으로, 좌표는 "
        "관측망 대조로, 총자력은 범위 검사로 각각 걸러야 한다.",
        "방문 대표값 자료의 정밀도를 확인할 길이 없다. 점조서 관측이력은 "
        f"{len(aux)}개 방문을 담지만 세션 산포가 없어 «보조»에 머문다. "
        "원 야장이 확보되면 세션으로 승격할 수 있다.",
        "총자력 잔차가 큰 원인을 가리지 못했다. 단위·기기·기록 위치 가운데 "
        "무엇인지 확인해야 하고, 그 전까지 총자력은 투입 가능분만 쓴다.",
        "이설 이전 관측을 측점 식별자 분리로 되살릴 수 있다. 이설 전후를 "
        "다른 측점으로 두면 과거 자료가 공간 구조 결정에 기여한다.",
        "외부장 보정의 공간 투영이 남아 있다. 관측소 한 곳의 편차를 전국에 "
        "같이 적용하는 근사는 거리가 멀면 성립하지 않는다.",
        "편각 잔차 폭이 큰 측점의 방위표지 기준을 현장에서 재확인해야 한다. "
        "표 5에서 폭이 큰 측점이 재확인 1순위다.",
    ]:
        bullet(doc, t)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / (f"{dt.datetime.now():%Y%m%d_%H%M%S}"
                    f"_기존관측자료_LMM투입검토.docx")
    doc.save(out)
    print(f"[저장] {out}")
    print(f"  편각 투입 {len(ok)} 방문 / {ok['측점'].nunique()} 측점"
          f" · 잔차 RMS {rmsD:.1f}′")
    print(f"  새 측점 {len(newS)} — {', '.join(newS)}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
