# -*- coding: utf-8 -*-
"""
LMM 4개 층 자료 확보 현황 — 한 장으로 보는 정리
================================================

Core / Regional / Crustal / External 네 층의 자료를 실제 파일에서 세어 HTML 한 장으로
낸다. 수치를 손으로 적지 않고 **파일을 직접 읽어** 세므로, 자료가 늘면 재실행만으로
갱신된다.

    python make_data_inventory.py      ->  docs/data_inventory.html
"""
import datetime as dt
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import pandas as pd

ROOT = Path(__file__).parent
DATA = ROOT / "data"
DOCS = ROOT / "docs"
DDATA = DOCS / "data"
OUT = DOCS / "data_inventory.html"

KIGAM_CITE = ("박영수, 임형래, 임무택, 신영홍, 2019, 한국의 자력이상도, "
              "지구물리와 물리탐사, 22(1), 29–36.")


def human(n):
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024 or u == "GB":
            return f"{n:,.0f} {u}" if u == "B" else f"{n:.1f} {u}"
        n /= 1024


def size_of(p):
    try:
        return p.stat().st_size
    except OSError:
        return 0


# ── 층별 실측 ───────────────────────────────────────────────────────────────
def core_layer():
    rows = []
    try:
        import ppigrf
        shc = Path(ppigrf.__file__).parent / "IGRF14.shc"
        rows.append(("IGRF-14 계수", "ppigrf/IGRF14.shc",
                     "차수 13 · 5년 주기", human(size_of(shc)) if shc.exists() else "패키지 내장"))
    except Exception:
        rows.append(("IGRF-14 계수", "ppigrf", "차수 13", "미확인"))
    js = DOCS / "igrf14.js"
    if js.exists():
        rows.append(("웹 IGRF 엔진", "docs/igrf14.js", "lmm.html 엔진 추출본",
                     human(size_of(js))))
    return rows, "완비", "ppigrf 와 표시 자릿수 내 일치 검증 통과"


def regional_layer():
    rows = []
    sp = DATA / "지자기측량 성과정리(22_25).xlsx"
    n_sp = None
    if sp.exists():
        try:
            from lmm_build import load_survey_points
            n_sp = len(load_survey_points())
        except Exception:
            pass
        rows.append(("성과표", "data/지자기측량 성과정리(22_25).xlsx",
                     f"{n_sp}행" if n_sp else "-", human(size_of(sp))))
    fb = DDATA / "fieldbook_sessions.csv"
    n_fb = n_site = n_time = 0
    yrs = ""
    if fb.exists():
        d = pd.read_csv(fb, encoding="utf-8-sig")
        n_fb, n_site = len(d), d["측점"].nunique()
        n_time = d["시작"].notna().sum()
        y = pd.to_datetime(d["날짜"]).dt.year
        yrs = f"{y.min()}~{y.max()}"
        rows.append(("야장 세션표", "docs/data/fieldbook_sessions.csv",
                     f"{n_fb}세션 · {n_site}측점 · {yrs}", human(size_of(fb))))
    rows.append(("야장 원본", "D:\\...\\20260811_지리원 야장 자료\\지구물리측량\\지자기측량",
                 "연도·측점별 폴더 (2020~2025)", "약 4 GB"))
    note = (f"관측시각 {n_time}/{n_fb}세션 확보. 성과표 25행 중 18행에 실일시 부착. "
            "권고 30측점 대비 미달")
    return rows, "부족", note


def crustal_layer():
    rows = []
    dat = DATA / "mag_1982-2018_1.5min_ed.dat"
    n = 0
    if dat.exists():
        n = sum(1 for _ in open(dat, encoding="utf-8", errors="ignore"))
        rows.append(("KIGAM 수치자력이상도", "data/mag_1982-2018_1.5min_ed.dat",
                     f"1.5분 등격자 · {n:,}행 · WGS84", human(size_of(dat))))
    for f in ("수치자력이상도를+사용할+때+주의할+사항.pdf",
              "수치자력이상도를+사용할+때+주의할+사항.hwp"):
        p = DATA / f
        if p.exists():
            rows.append(("사용 주의사항", f"data/{f}", "이용조건·인용 규정",
                         human(size_of(p))))
    note = ("2.8 km 격자만 보유. 원측선 자료는 존재하지 않으므로(2026-08 확인) "
            "현 격자가 지각층 해상도의 상한이다")
    return rows, "해상도 부족", note


def external_layer():
    rows = []
    cyg = DATA / "cyg"
    if cyg.exists():
        for grade in ("best-avail", "quasi-def"):
            g = cyg / grade
            if not g.exists():
                continue
            years = sorted(p.name for p in g.iterdir() if p.is_dir())
            nday = sum(len(list((g / y).glob('*.csv'))) for y in years)
            rows.append((f"CYG 청양 ({grade})", f"data/cyg/{grade}/",
                         f"{years[0]}~{years[-1]} · {nday:,}일", "1분 · UTC"))
    kasa = DATA / "kasa"
    if kasa.exists():
        by = {}
        for p in kasa.glob("kasa_*.csv"):
            st = p.stem.split("_")[1]
            by.setdefault(st, []).append(p.stem.split("_")[2])
        names = {"GN": "강릉", "ICH": "이천", "JJ": "제주", "BOH": "보현산"}
        for st, ys in sorted(by.items()):
            ys = sorted(ys)
            rows.append((f"KASA {names.get(st, st)}", f"data/kasa/kasa_{st}_*.csv",
                         f"{ys[0]}~{ys[-1]} · {len(ys)}개 파일", "1분 · KST"))
    kp = DATA / "cyg" / "Kp_ap_since_1932.txt"
    if kp.exists():
        n = sum(1 for _ in open(kp, encoding="utf-8", errors="ignore"))
        rows.append(("Kp 지수", "data/cyg/Kp_ap_since_1932.txt",
                     f"{n:,}행 · 1932~현재", "3시간 · UT"))
    ec = DDATA / "external_corrections.csv"
    if ec.exists():
        d = pd.read_csv(ec, encoding="utf-8-sig")
        ok = (d["상태"] == "정상").sum()
        rows.append(("산출 보정량", "docs/data/external_corrections.csv",
                     f"{len(d)}세션 중 정상 {ok}", human(size_of(ec))))
    note = ("자료는 모두 확보. 그러나 균일 V 근사로 차감하면 LOO 가 악화되어 "
            "미적용(EXTERNAL_MODE=none). NOC 공간투영 구현이 다음 요건")
    return rows, "자료 확보 · 적용 보류", note


def model_numbers():
    p = DDATA / "lmm_model.json"
    if not p.exists():
        return {}
    m = json.load(open(p, encoding="utf-8"))
    return {"D": m["loo_cv"]["D"], "I": m["loo_cv"]["I"], "F": m["loo_cv"]["F"],
            "n": len(m.get("sites", [])), "epoch": m.get("epoch")}


STATUS_COLOR = {"완비": "#1B7F5A", "부족": "#C2410C", "해상도 부족": "#C2410C",
                "자료 확보 · 적용 보류": "#B45309"}


def build_html():
    layers = [("① Core", "지구 내부 장주기장", *core_layer()),
              ("② Regional", "한반도 지역장 · 지표 절대측정", *regional_layer()),
              ("③ Crustal", "지각 자기이상", *crustal_layer()),
              ("④ External", "외부장 시간변동", *external_layer())]
    mm = model_numbers()

    def rowhtml(r):
        return ("<tr><td class='k'>" + r[0] + "</td><td class='p'>" + r[1] +
                "</td><td>" + r[2] + "</td><td class='s'>" + r[3] + "</td></tr>")

    secs = []
    for name, sub, rows, status, note in layers:
        col = STATUS_COLOR.get(status, "#475569")
        secs.append(f"""
<section class="layer">
  <div class="lhead">
    <div><h2>{name}</h2><p class="sub">{sub}</p></div>
    <span class="badge" style="background:{col}">{status}</span>
  </div>
  <table><thead><tr><th>자료</th><th>경로</th><th>내용</th><th>비고</th></tr></thead>
  <tbody>{''.join(rowhtml(r) for r in rows)}</tbody></table>
  <p class="note">{note}</p>
</section>""")

    kpi = ""
    if mm:
        def cell(lab, val, target, ok):
            c = "#1B7F5A" if ok else "#C2410C"
            return (f"<div class='kpi'><div class='v' style='color:{c}'>{val}</div>"
                    f"<div class='l'>{lab}</div><div class='t'>{target}</div></div>")
        kpi = ("<div class='kpis'>"
               + cell("LOO 편각 D", f"{mm['D']:.3f}°", "목표 &lt; 0.1°", False)
               + cell("LOO 총자력 F", f"{mm['F']:.1f} nT", "목표 &lt; 50 nT", False)
               + cell("LOO 복각 I", f"{mm['I']:.3f}°", "별도 목표 없음", True)
               + cell("Regional 측점", f"{mm['n']}점", "권고 30점 이상", False)
               + "</div>")

    return f"""<!doctype html>
<html lang="ko"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>LMM 자료 확보 현황</title>
<style>
:root{{--bg:#fff;--fg:#182638;--mut:#607184;--line:#E3E9F0;--card:#F7F9FC}}
@media (prefers-color-scheme:dark){{:root{{--bg:#0F172A;--fg:#E8EEF6;--mut:#94A8BD;
  --line:#24344A;--card:#16233A}}}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--fg);
  font-family:Pretendard,"맑은 고딕","Malgun Gothic",system-ui,sans-serif;
  line-height:1.55;padding:28px 20px 60px}}
.wrap{{max-width:1080px;margin:0 auto}}
h1{{font-size:26px;margin:0 0 4px}}
.lead{{color:var(--mut);margin:0 0 22px;font-size:14px}}
.kpis{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));
  gap:12px;margin:0 0 26px}}
.kpi{{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}}
.kpi .v{{font-size:24px;font-weight:700}}
.kpi .l{{font-size:12.5px;margin-top:2px}}
.kpi .t{{font-size:11.5px;color:var(--mut)}}
.layer{{border:1px solid var(--line);border-radius:12px;padding:16px 18px;margin-bottom:16px}}
.lhead{{display:flex;justify-content:space-between;align-items:flex-start;gap:12px}}
h2{{font-size:17px;margin:0}}
.sub{{color:var(--mut);font-size:12.5px;margin:2px 0 10px}}
.badge{{color:#fff;font-size:12px;font-weight:600;padding:4px 10px;border-radius:999px;
  white-space:nowrap}}
table{{width:100%;border-collapse:collapse;font-size:13px;margin-top:6px;display:block;
  overflow-x:auto}}
th,td{{text-align:left;padding:7px 10px;border-bottom:1px solid var(--line);
  vertical-align:top;white-space:nowrap}}
th{{color:var(--mut);font-weight:600;font-size:12px}}
td.k{{font-weight:600}}
td.p{{font-family:Consolas,monospace;font-size:12px;color:var(--mut)}}
td.s{{color:var(--mut)}}
.note{{background:var(--card);border-radius:8px;padding:9px 12px;margin:10px 0 0;
  font-size:12.5px;color:var(--mut)}}
footer{{margin-top:26px;font-size:12px;color:var(--mut);border-top:1px solid var(--line);
  padding-top:14px}}
code{{font-family:Consolas,monospace;font-size:12px}}
</style></head><body><div class="wrap">
<h1>LMM 4개 층 자료 확보 현황</h1>
<p class="lead">파일을 직접 읽어 센 값이다. 자료가 늘면
<code>python make_data_inventory.py</code> 재실행으로 갱신된다.
&nbsp;·&nbsp; 생성 {dt.datetime.now():%Y-%m-%d %H:%M}</p>
{kpi}
{''.join(secs)}
<footer>
<b>시간대 주의</b> — 야장 KST · CYG UTC · Kp UT · KASA KST. 결합 지점마다 변환 확인.<br>
<b>KIGAM 자료 이용조건</b> — 한국지질자원연구원장의 사전승인 없이 무단 복제·대여·양도·
판매·국외반출 금지. 인용: {KIGAM_CITE}
</footer>
</div></body></html>"""


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(build_html(), encoding="utf-8")
    print(f"[저장] {OUT}")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
