# -*- coding: utf-8 -*-
"""
한반도 LMM 엑셀 계산기 생성기
==============================

수식만으로(매크로·VBA 없이) LMM 을 계산하는 xlsx 를 생성한다.

구조
----
IGRF 층은 13차 구면조화 합성이 필요하지만, 엑셀 수식으로 직접 구현하면
105개 (n,m) 항의 르장드르 재귀를 셀로 풀어야 해 취약하다.
IGRF 는 공간적으로 매우 매끄러운 장이므로, 대신 0.1° 격자로 미리 계산해
쌍선형 보간한다(검증: 최대오차 D 0.09", F 0.017 nT — 모델 자체 불확도의 1/5000).

    시트          내용
    계산          입력(위경도·표고·연도)과 결과. 사용자가 쓰는 유일한 시트
    IGRF_D/I/F    0.1° 격자의 2025.0 기준값 + 연변화율(SV)
    CRUSTAL       KIGAM 자력이상 격자 (0.025°)
    REGIONAL      지표 절대측정 적합 다항식 계수
    검증          Python 모델과 엑셀 수식 결과 대조표
    설명          방법론·정확도 고지

실행:
    python lmm_excel.py
"""

import json
import datetime as dt
from pathlib import Path

import numpy as np
import xlsxwriter

from lmm_build import igrf_dif, load_kigam_grid, CrustalGrid, poly_terms

BASE = Path(__file__).parent
MODEL_JSON = BASE / "docs" / "data" / "lmm_model.json"
DOCS_OUT = BASE / "docs" / "output"

# IGRF 사전격자 사양
IG_STEP = 0.1
IG_LAT = np.round(np.arange(33.0, 38.60001, IG_STEP), 4)
IG_LON = np.round(np.arange(125.5, 129.90001, IG_STEP), 4)
IG_EPOCH = 2025.0
IG_SV_YEARS = 5.0  # 2025 -> 2030 차분으로 연변화율 산정

# 표고 보정 상수 (한반도 중심 기준 선형 기울기)
GRAD_REF = (36.0, 128.0)

# 계산 시트 기본 입력값 (캐시값 산정 기준)
DEFAULT_INPUT = dict(lat=36.0, lon=128.0, elev=0.0, year=2027.0)


def bilinear_np(lats, lons, grid, lat, lon):
    """
    엑셀 수식(MATCH+INDEX 쌍선형)과 동일한 보간을 파이썬으로 수행.

    xlsxwriter 는 수식의 계산 결과를 알지 못해 캐시값으로 0 을 기록한다.
    Excel 은 fullCalcOnLoad 플래그로 열 때 재계산하지만, 미리보기·
    LibreOffice 변환·일부 뷰어는 캐시값을 그대로 보여 0 으로 표시된다.
    그래서 여기서 계산한 값을 write_formula(value=...) 로 함께 넣는다.
    """
    lats = np.asarray(lats, float)
    lons = np.asarray(lons, float)
    j = int(np.clip(np.searchsorted(lats, lat, side="right") - 1, 0, lats.size - 2))
    i = int(np.clip(np.searchsorted(lons, lon, side="right") - 1, 0, lons.size - 2))
    ty = (lat - lats[j]) / (lats[j + 1] - lats[j])
    tx = (lon - lons[i]) / (lons[i + 1] - lons[i])

    v = [grid[j, i], grid[j, i + 1], grid[j + 1, i], grid[j + 1, i + 1]]
    w = [(1 - tx) * (1 - ty), tx * (1 - ty), (1 - tx) * ty, tx * ty]
    num = sum(wi * vi for wi, vi in zip(w, v) if np.isfinite(vi))
    den = sum(wi for wi, vi in zip(w, v) if np.isfinite(vi))
    return num / den if den > 0 else 0.0


def build_igrf_grids():
    """0.1° 격자에서 D·I·F 의 2025.0 값과 연변화율을 계산."""
    LA, LO = np.meshgrid(IG_LAT, IG_LON, indexing="ij")
    flat_la, flat_lo = LA.ravel(), LO.ravel()
    z = np.zeros(flat_la.size)

    d0 = dt.datetime(2025, 1, 1)
    d1 = dt.datetime(2030, 1, 1)

    D0, I0, F0, *_ = igrf_dif(flat_la, flat_lo, z, d0)
    D1, I1, F1, *_ = igrf_dif(flat_la, flat_lo, z, d1)

    shape = LA.shape
    out = {}
    for key, a, b in [("D", D0, D1), ("I", I0, I1), ("F", F0, F1)]:
        out[key] = {
            "val": a.reshape(shape),
            "sv": ((b - a) / IG_SV_YEARS).reshape(shape),
        }
    return out


def height_gradients():
    """표고 1 m 당 D·I·F 변화율 (한반도 중심 기준)."""
    la, lo = GRAD_REF
    d = dt.datetime(2027, 1, 1)
    lo_h, hi_h = 0.0, 1000.0
    D0, I0, F0, *_ = igrf_dif([la], [lo], [lo_h], d)
    D1, I1, F1, *_ = igrf_dif([la], [lo], [hi_h], d)
    dh = hi_h - lo_h
    return {
        "D": float((D1[0] - D0[0]) / dh),
        "I": float((I1[0] - I0[0]) / dh),
        "F": float((F1[0] - F0[0]) / dh),
    }


# --------------------------------------------------------------- 수식 조립
def bilinear_formula(lat_cell, lon_cell, lat_name, lon_name, val_name,
                     handle_gaps=False):
    """
    명명범위 기반 쌍선형 보간 수식 문자열을 만든다.

    handle_gaps=True 이면 빈 셀(자료 공백)을 제외하고 유효 이웃만
    가중평균한다 — CRUSTAL 격자용.
    """
    j = f"MATCH({lat_cell},{lat_name},1)"
    i = f"MATCH({lon_cell},{lon_name},1)"

    la0 = f"INDEX({lat_name},{j})"
    la1 = f"INDEX({lat_name},{j}+1)"
    lo0 = f"INDEX({lon_name},{i})"
    lo1 = f"INDEX({lon_name},{i}+1)"

    ty = f"(({lat_cell}-{la0})/({la1}-{la0}))"
    tx = f"(({lon_cell}-{lo0})/({lo1}-{lo0}))"

    v00 = f"INDEX({val_name},{j},{i})"
    v10 = f"INDEX({val_name},{j},{i}+1)"
    v01 = f"INDEX({val_name},{j}+1,{i})"
    v11 = f"INDEX({val_name},{j}+1,{i}+1)"

    w00 = f"(1-{tx})*(1-{ty})"
    w10 = f"{tx}*(1-{ty})"
    w01 = f"(1-{tx})*{ty}"
    w11 = f"{tx}*{ty}"

    if not handle_gaps:
        return f"={v00}*{w00}+{v10}*{w10}+{v01}*{w01}+{v11}*{w11}"

    num = "+".join(
        f"IF(ISNUMBER({v}),{v}*({w}),0)"
        for v, w in [(v00, w00), (v10, w10), (v01, w01), (v11, w11)]
    )
    den = "+".join(
        f"IF(ISNUMBER({v}),({w}),0)"
        for v, w in [(v00, w00), (v10, w10), (v01, w01), (v11, w11)]
    )
    return f"=IF(({den})=0,0,({num})/({den}))"


def write_grid_sheet(wb, name, lats, lons, blocks, fmt_hdr, fmt_num, prefix):
    """
    격자 시트를 쓰고 명명범위를 등록한다.

    blocks: [(제목, 2차원배열, 이름접미사), ...] — 세로로 이어 붙인다.
    """
    ws = wb.add_worksheet(name)
    ws.freeze_panes(0, 1)
    row = 0
    for title, arr, suffix in blocks:
        ws.write(row, 0, title, fmt_hdr)
        row += 1
        hdr_row = row
        for c, lo in enumerate(lons):
            ws.write_number(hdr_row, c + 1, float(lo), fmt_hdr)
        row += 1
        first_data = row
        for r, la in enumerate(lats):
            ws.write_number(row, 0, float(la), fmt_hdr)
            for c in range(len(lons)):
                v = arr[r, c]
                if v is None or (isinstance(v, float) and not np.isfinite(v)):
                    continue
                ws.write_number(row, c + 1, float(v), fmt_num)
            row += 1
        last_data = row - 1

        def a1(r0, c0, r1, c1):
            return (f"'{name}'!${xlsxwriter.utility.xl_col_to_name(c0)}${r0+1}"
                    f":${xlsxwriter.utility.xl_col_to_name(c1)}${r1+1}")

        wb.define_name(f"{prefix}_LAT{suffix}", a1(first_data, 0, last_data, 0))
        wb.define_name(f"{prefix}_LON{suffix}", a1(hdr_row, 1, hdr_row, len(lons)))
        wb.define_name(f"{prefix}_VAL{suffix}",
                       a1(first_data, 1, last_data, len(lons)))
        row += 2
    ws.set_column(0, 0, 9)
    ws.set_column(1, len(lons), 9)
    return ws


def main():
    model = json.loads(MODEL_JSON.read_text(encoding="utf-8"))
    reg = model["regional"]
    deg = reg["degree"]

    print("[1/4] IGRF 사전격자 계산 중...")
    ig = build_igrf_grids()
    grad = height_gradients()
    print(f"      표고 기울기 F={grad['F']*1000:.2f} nT/km  "
          f"D={grad['D']*1000*3600:.2f}\"/km")

    print("[2/4] KIGAM 격자 적재 중...")
    lons, lats, cgrid = load_kigam_grid()

    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out = DOCS_OUT / f"{stamp}_LMM_계산기.xlsx"
    DOCS_OUT.mkdir(parents=True, exist_ok=True)

    wb = xlsxwriter.Workbook(str(out), {"nan_inf_to_errors": True})

    f_title = wb.add_format({"bold": True, "font_size": 14})
    f_hdr = wb.add_format({"bold": True, "bg_color": "#EDF1F5", "border": 1,
                           "align": "center", "num_format": "0.###"})
    f_lbl = wb.add_format({"bold": True, "bg_color": "#F7F9FB", "border": 1})
    f_in = wb.add_format({"bg_color": "#FFF8E1", "border": 1,
                          "num_format": "0.######"})
    f_out = wb.add_format({"border": 1, "num_format": "0.0000"})
    f_outn = wb.add_format({"border": 1, "num_format": "0.0"})
    f_num = wb.add_format({"num_format": "0.###"})
    f_note = wb.add_format({"text_wrap": True, "valign": "top"})
    f_warn = wb.add_format({"text_wrap": True, "valign": "top",
                            "bg_color": "#FFF4E5", "border": 1})

    # ---------------------------------------------------------- 격자 시트
    print("[3/4] 시트 작성 중...")
    for comp in ("D", "I", "F"):
        write_grid_sheet(
            wb, f"IGRF_{comp}", IG_LAT, IG_LON,
            [(f"IGRF-14 {comp} @ Epoch 2025.0 (표고 0 m)", ig[comp]["val"], ""),
             (f"연변화율 SV ({comp}/year)", ig[comp]["sv"], "SV")],
            f_hdr, f_num, f"IG{comp}",
        )

    write_grid_sheet(
        wb, "CRUSTAL", lats, lons,
        [("KIGAM 자력이상 (nT) — 빈 셀은 항공자력 미측선 구역", cgrid, "")],
        f_hdr, f_num, "CR",
    )

    # 지각 이상벡터 — 스칼라 이상에서 파수영역 역산으로 복원한 북·동 성분.
    # 하향 성분은 스칼라에서 유도하므로 시트를 따로 두지 않는다:
    #   b_down = (dF + dc - l*b_north - m*b_east) / n
    vec = model["crustal"].get("vector")
    if vec:
        bn = np.array(vec["b_north"], float).reshape(len(lats), len(lons))
        be = np.array(vec["b_east"], float).reshape(len(lats), len(lons))
        write_grid_sheet(
            wb, "CRUSTAL_VEC", lats, lons,
            [("지각 이상벡터 북성분 b_north (nT)", bn, ""),
             ("지각 이상벡터 동성분 b_east (nT)", be, "E")],
            f_hdr, f_num, "CV",
        )

    # ------------------------------------------------------- REGIONAL 시트
    wsr = wb.add_worksheet("REGIONAL")
    wsr.write(0, 0, "Regional 층 다항식 계수", f_title)
    wsr.write(1, 0, f"기준점 lat0={reg['lat0']}, lon0={reg['lon0']}, 차수={deg}")
    wsr.write_row(3, 0, ["항", "D (deg)", "I (deg)", "F (nT)"], f_hdr)
    for k, term in enumerate(reg["terms"]):
        wsr.write_row(4 + k, 0, [term, reg["D"][k], reg["I"][k], reg["F"][k]])
    wsr.set_column(0, 0, 12)
    wsr.set_column(1, 3, 16)
    n_terms = len(reg["terms"])
    wb.define_name("REG_D", f"'REGIONAL'!$B$5:$B${4+n_terms}")
    wb.define_name("REG_I", f"'REGIONAL'!$C$5:$C${4+n_terms}")
    wb.define_name("REG_F", f"'REGIONAL'!$D$5:$D${4+n_terms}")

    # --------------------------------------------------------- 계산 시트
    ws = wb.add_worksheet("계산")
    wb.worksheets_objs.insert(0, wb.worksheets_objs.pop())  # 첫 시트로 이동
    ws.set_column(0, 0, 26)
    ws.set_column(1, 4, 16)
    # 인쇄·PDF 변환 시 한 장에 담기도록 (A~E 합계가 세로 용지 폭을 넘고,
    # 정확도 고지가 페이지 경계에 걸려 잘리는 것을 막는다)
    ws.fit_to_pages(1, 1)
    ws.set_margins(0.4, 0.4, 0.4, 0.4)
    ws.write(0, 0, "한반도 LMM 계산기", f_title)
    ws.write(1, 0, "노란 칸에 값을 입력하면 아래 결과가 자동 계산됩니다. "
                   "매크로 없이 수식만 사용합니다.")

    ws.write(3, 0, "위도 (°N)", f_lbl);  ws.write_number(3, 1, 36.0, f_in)
    ws.write(4, 0, "경도 (°E)", f_lbl);  ws.write_number(4, 1, 128.0, f_in)
    ws.write(5, 0, "표고 (m)", f_lbl);   ws.write_number(5, 1, 0.0, f_in)
    ws.write(6, 0, "연도 (예: 2027.5)", f_lbl); ws.write_number(6, 1, 2027.0, f_in)

    LAT, LON, ELV, YR = "$B$4", "$B$5", "$B$6", "$B$7"

    # 중간 계산(감사 가능하도록 노출)
    ws.write(8, 0, "── 층별 계산 ──", f_lbl)

    # 기본 입력값에 대한 캐시값 (뷰어가 재계산하지 않아도 값이 보이도록)
    d_in = DEFAULT_INPUT
    cache = {}
    for comp in ("D", "I", "F"):
        base = bilinear_np(IG_LAT, IG_LON, ig[comp]["val"], d_in["lat"], d_in["lon"])
        sv = bilinear_np(IG_LAT, IG_LON, ig[comp]["sv"], d_in["lat"], d_in["lon"])
        h = d_in["elev"] * grad[comp]
        cache[comp] = dict(base=base, sv=sv, h=h,
                           fin=base + sv * (d_in["year"] - IG_EPOCH) + h)

    r = 9
    rows_layer = []
    for comp, unit in [("D", "deg"), ("I", "deg"), ("F", "nT")]:
        c = cache[comp]
        ws.write(r, 0, f"① IGRF {comp} @2025 ({unit})", f_lbl)
        ws.write_formula(
            r, 1,
            bilinear_formula(LAT, LON, f"IG{comp}_LAT", f"IG{comp}_LON",
                             f"IG{comp}_VAL"),
            f_out, c["base"])
        ws.write(r + 1, 0, f"   SV {comp} (/year)", f_lbl)
        ws.write_formula(
            r + 1, 1,
            bilinear_formula(LAT, LON, f"IG{comp}_LATSV", f"IG{comp}_LONSV",
                             f"IG{comp}_VALSV"),
            f_out, c["sv"])
        ws.write(r + 2, 0, f"   표고보정 {comp}", f_lbl)
        ws.write_formula(r + 2, 1, f"={ELV}*{grad[comp]!r}", f_out, c["h"])
        ws.write(r + 3, 0, f"   IGRF {comp} 최종", f_lbl)
        ws.write_formula(
            r + 3, 1,
            f"=B{r+1}+B{r+2}*({YR}-{IG_EPOCH})+B{r+3}", f_out, c["fin"])
        rows_layer.append((comp, r + 3))  # 0-based row of 최종
        r += 5

    # Regional 다항식
    ws.write(r, 0, "② Regional 다항식", f_lbl)
    terms = []
    for total in range(deg + 1):
        for p in range(total + 1):
            terms.append(f"(({LAT}-{reg['lat0']})^{total-p})*(({LON}-{reg['lon0']})^{p})")
    reg_rows = {}
    reg_val = {}
    A_def = poly_terms(np.array([d_in["lat"]]), np.array([d_in["lon"]]), deg)
    for comp in ("D", "I", "F"):
        expr = "+".join(
            f"INDEX(REG_{comp},{k+1})*{t}" for k, t in enumerate(terms)
        )
        reg_val[comp] = float((A_def @ np.array(reg[comp]))[0])
        ws.write(r + 1 + ("DIF".index(comp)), 0, f"   Regional {comp}", f_lbl)
        ws.write_formula(r + 1 + ("DIF".index(comp)), 1, f"={expr}", f_out,
                         reg_val[comp])
        reg_rows[comp] = r + 1 + ("DIF".index(comp))
    r += 5

    # Crustal
    crust_val = bilinear_np(lats, lons, cgrid, d_in["lat"], d_in["lon"])
    ws.write(r, 0, "③ Crustal F (nT)", f_lbl)
    ws.write_formula(
        r, 1,
        bilinear_formula(LAT, LON, "CR_LAT", "CR_LON", "CR_VAL",
                         handle_gaps=True),
        f_outn, crust_val)
    crust_row = r
    ws.write(r + 1, 0, "   자료 유무", f_lbl)
    ws.write_formula(
        r + 1, 1,
        f'=IF(OR({LAT}<MIN(CR_LAT),{LAT}>MAX(CR_LAT),{LON}<MIN(CR_LON),'
        f'{LON}>MAX(CR_LON)),"격자 범위 밖",IF(B{r+1}=0,"자료 공백(0 적용)","정상"))',
        None, "정상" if crust_val != 0 else "자료 공백(0 적용)")
    r += 3

    # ---- 지각 이상벡터 -> 편각·복각 기여 -------------------------------
    crD_row = crI_row = None
    if vec:
        l_, m_, n_, dc_ = vec["l"], vec["m"], vec["n"], vec["dc_nT"]
        bn_v = bilinear_np(lats, lons, bn, d_in["lat"], d_in["lon"])
        be_v = bilinear_np(lats, lons, be, d_in["lat"], d_in["lon"])
        bd_v = (crust_val + dc_ - l_ * bn_v - m_ * be_v) / n_
        igD_c, igI_c, igF_c = (cache["D"]["fin"], cache["I"]["fin"],
                               cache["F"]["fin"])
        Xi = igF_c * np.cos(np.radians(igI_c)) * np.cos(np.radians(igD_c))
        Yi = igF_c * np.cos(np.radians(igI_c)) * np.sin(np.radians(igD_c))
        Zi = igF_c * np.sin(np.radians(igI_c))
        crD_v = float(np.degrees(np.arctan2(Yi + be_v, Xi + bn_v)) - igD_c)
        crI_v = float(np.degrees(np.arctan2(Zi + bd_v,
                                            np.hypot(Xi + bn_v, Yi + be_v)))
                      - igI_c)

        # Excel 행 참조 (1-based)
        eD, eI, eF = (rows_layer[0][1] + 1, rows_layer[1][1] + 1,
                      rows_layer[2][1] + 1)
        ws.write(r, 0, "③-b 지각벡터 b_north (nT)", f_lbl)
        ws.write_formula(r, 1, bilinear_formula(LAT, LON, "CV_LAT", "CV_LON",
                                                "CV_VAL"), f_outn, bn_v)
        ws.write(r + 1, 0, "     b_east (nT)", f_lbl)
        ws.write_formula(r + 1, 1, bilinear_formula(LAT, LON, "CV_LATE",
                                                    "CV_LONE", "CV_VALE"),
                         f_outn, be_v)
        ws.write(r + 2, 0, "     b_down (유도, nT)", f_lbl)
        ws.write_formula(
            r + 2, 1,
            f"=(B{crust_row+1}+{dc_!r}-{l_!r}*B{r+1}-{m_!r}*B{r+2})/{n_!r}",
            f_outn, bd_v)
        # IGRF 벡터 성분 (표시하지 않고 식 안에서 조립)
        Xe = f"(B{eF}*COS(RADIANS(B{eI}))*COS(RADIANS(B{eD}))+B{r+1})"
        Ye = f"(B{eF}*COS(RADIANS(B{eI}))*SIN(RADIANS(B{eD}))+B{r+2})"
        Ze = f"(B{eF}*SIN(RADIANS(B{eI}))+B{r+3})"
        ws.write(r + 3, 0, "     Crustal D 기여 (°)", f_lbl)
        ws.write_formula(r + 3, 1,
                         f"=DEGREES(ATAN2({Xe},{Ye}))-B{eD}", f_out, crD_v)
        ws.write(r + 4, 0, "     Crustal I 기여 (°)", f_lbl)
        ws.write_formula(
            r + 4, 1,
            f"=DEGREES(ATAN2(SQRT({Xe}^2+{Ye}^2),{Ze}))-B{eI}", f_out, crI_v)
        crD_row, crI_row = r + 3, r + 4
        r += 6

    ws.write(r, 0, "④ External (CYG)", f_lbl)
    ws.merge_range(r, 1, r, 4,
                   "예측 시점에는 미적용 — LMM 은 정온시 기준값 모형입니다. "
                   "외부장은 관측 성과 정리 단계에서 관측소 4소(청양·제주·강릉·이천) "
                   "공간보간으로 편각에 반영했습니다 ('설명' 시트 참조).")
    r += 2

    # 최종 결과
    igD, igI, igF = (rows_layer[0][1], rows_layer[1][1], rows_layer[2][1])
    ws.write(r, 0, "── LMM 최종 결과 ──", f_lbl)
    r += 1
    res_start = r
    # 최종 결과 캐시값
    add_D = crD_v if crD_row is not None else 0.0
    add_I = crI_v if crI_row is not None else 0.0
    vD = cache["D"]["fin"] + add_D + reg_val["D"]
    vI = cache["I"]["fin"] + add_I + reg_val["I"]
    vF = cache["F"]["fin"] + crust_val + reg_val["F"]
    vH = vF * np.cos(np.radians(vI))
    vX, vY = vH * np.cos(np.radians(vD)), vH * np.sin(np.radians(vD))
    vZ = vF * np.sin(np.radians(vI))

    exD = f"+B{crD_row+1}" if crD_row is not None else ""
    exI = f"+B{crI_row+1}" if crI_row is not None else ""
    ws.write(r, 0, "편각 D (°)", f_lbl)
    ws.write_formula(r, 1, f"=B{igD+1}{exD}+B{reg_rows['D']+1}", f_out, vD)
    ws.write(r + 1, 0, "복각 I (°)", f_lbl)
    ws.write_formula(r + 1, 1, f"=B{igI+1}{exI}+B{reg_rows['I']+1}", f_out, vI)
    ws.write(r + 2, 0, "총자력 F (nT)", f_lbl)
    ws.write_formula(r + 2, 1,
                     f"=B{igF+1}+B{crust_row+1}+B{reg_rows['F']+1}", f_outn, vF)
    ws.write(r + 3, 0, "수평 H (nT)", f_lbl)
    ws.write_formula(r + 3, 1, f"=B{r+3}*COS(RADIANS(B{r+2}))", f_outn, vH)
    ws.write(r + 4, 0, "북 X (nT)", f_lbl)
    ws.write_formula(r + 4, 1, f"=B{r+4}*COS(RADIANS(B{r+1}))", f_outn, vX)
    ws.write(r + 5, 0, "동 Y (nT)", f_lbl)
    ws.write_formula(r + 5, 1, f"=B{r+4}*SIN(RADIANS(B{r+1}))", f_outn, vY)
    ws.write(r + 6, 0, "연직 Z (nT)", f_lbl)
    ws.write_formula(r + 6, 1, f"=B{r+3}*SIN(RADIANS(B{r+2}))", f_outn, vZ)

    # 편각 도분초 표기
    sign = "-" if vD < 0 else ""
    a = abs(vD)
    dd, mm = int(a), int((a - int(a)) * 60)
    ss = ((a - int(a)) * 60 - mm) * 60
    dms_txt = f'{sign}{dd}° {mm:02d}\' {ss:04.1f}"'

    ws.write(r + 8, 0, "편각 D (도 분 초)", f_lbl)
    ws.write_formula(
        r + 8, 1,
        f'=TEXT(TRUNC(B{res_start+1}),"0")&"° "&'
        f'TEXT(TRUNC(ABS(B{res_start+1}-TRUNC(B{res_start+1}))*60),"00")&"\' "&'
        f'TEXT(ROUND((ABS(B{res_start+1}-TRUNC(B{res_start+1}))*60-'
        f'TRUNC(ABS(B{res_start+1}-TRUNC(B{res_start+1}))*60))*60,1),"00.0")&""""',
        None, dms_txt)

    ws.write(r + 10, 0, "정확도 고지", f_lbl)
    # 병합 셀은 자동 행높이 조정이 되지 않으므로 명시적으로 지정한다
    for k in range(r + 10, r + 15):
        ws.set_row(k, 22)
    cv = model["loo_cv"]
    ws.merge_range(
        r + 10, 1, r + 14, 4,
        "본 모델은 4개 층 중 3개(IGRF·Regional·Crustal)만으로 구축된 잠정판입니다. "
        f"교차검증 실측오차 D {cv['D']:.2f}°, I {cv['I']:.2f}°, F {cv['F']:.0f} nT 로 "
        "목표 KPI(D<0.1°, F<50 nT)에 미달합니다. "
        "IGRF-14 대비 개선된 참고값으로만 사용하고, "
        "지형도 자침편차 등 정식 편각 산출에는 사용하지 마십시오. "
        "※ 위 KPI 는 설명자료가 제시한 공학적 목표치이며 법정 기준이 아닙니다. "
        "법정 측정오차 한계는 편각 정수차 30'(지구물리측량 작업규정 제20조)입니다. "
        "자세한 내용은 '설명' 시트를 참조하십시오.",
        f_warn)

    # ---------------------------------------------------------- 검증 시트
    wsv = wb.add_worksheet("검증")
    wsv.write(0, 0, "엑셀 수식 vs Python 모델 대조", f_title)
    wsv.write(1, 0, "동일 좌표에 대해 Python 이 계산한 값을 함께 기록했습니다. "
                    "엑셀 수식 결과와 일치해야 합니다.")
    wsv.write_row(3, 0,
                  ["위도", "경도", "표고", "연도",
                   "D(Python)", "I(Python)", "F(Python)"], f_hdr)

    cg = CrustalGrid(lons, lats, cgrid)
    tests = [(36.0, 128.0, 0, 2027.0), (37.5665, 126.978, 50, 2026.5),
             (34.870311, 128.592862, 91.0, 2027.0),
             (35.800783, 128.187658, 190.9, 2027.0),
             (33.524, 126.894, 30, 2028.0)]
    for k, (la, lo, el, yr) in enumerate(tests):
        y = int(yr)
        d = dt.datetime(y, 1, 1) + dt.timedelta(days=(yr - y) * 365.25)
        D0, I0, F0, *_ = igrf_dif([la], [lo], [el], d)
        A = poly_terms(np.array([la]), np.array([lo]), deg)
        cr = cg(np.array([la]), np.array([lo]))[0]
        crv = 0.0 if not np.isfinite(cr) else cr
        wsv.write_row(4 + k, 0, [
            la, lo, el, yr,
            float(D0[0] + (A @ np.array(reg["D"]))[0]),
            float(I0[0] + (A @ np.array(reg["I"]))[0]),
            float(F0[0] + crv + (A @ np.array(reg["F"]))[0]),
        ])
    wsv.set_column(0, 6, 15)
    wsv.fit_to_pages(1, 0)

    # ---------------------------------------------------------- 설명 시트
    wsd = wb.add_worksheet("설명")
    wsd.set_column(0, 0, 100)
    wsd.fit_to_pages(1, 0)
    lines = [
        ("한반도 LMM 계산기 — 방법론 및 한계", f_title),
        ("", None),
        ("1. 모델 구조", f_lbl),
        ("   B_LMM = B_IGRF + B_Regional + B_Crustal + B_External", None),
        ("   ① IGRF-14 (13차 구면조화) — 0.1° 격자 사전계산 후 쌍선형 보간.", None),
        ("      사전격자 오차: D 0.09\", F 0.017 nT (직접계산 대비)", None),
        ("   ② Regional — 지표 절대측정 성과에 적합한 다항식", None),
        ("   ③ Crustal — KIGAM 자력이상도 (1.5분 ≒ 2.8 km) 쌍선형 보간", None),
        ("      총자력 F 뿐 아니라 **편각·복각에도** 기여합니다. 항공자력이 주는", None),
        ("      스칼라 이상 ΔF 는 이상벡터를 주자기장 방향으로 투영한 값이므로,", None),
        ("      파수영역 역산으로 벡터 3성분(b_north·b_east·b_down)을 복원해", None),
        ("      IGRF 벡터에 더한 뒤 각을 다시 잽니다(CRUSTAL_VEC 시트).", None),
        ("   ④ External — 예측 시점 미적용 (LMM 은 정온시 기준값 모형)", None),
        ("      외부장이 들어갈 자리는 예측 시점이 아니라 관측 성과 정리 단계입니다.", None),
        ("      관측소 4소(청양·제주·강릉·이천)의 정온야간 기준선 대비 세션 편차를", None),
        ("      1차 평면으로 공간보간해 측점 위치의 편차를 추정하고, 성과 편각에서", None),
        ("      뺐습니다. 교란시(Kp>2) 세션은 공간보간이 먼저 깨지므로 보정하지 않습니다.", None),
        ("      총자력에는 적용하지 않습니다 — F 잔차는 외부장(약 38 nT)이 아니라", None),
        ("      지각장 불일치(약 190 nT)가 지배해 잡음만 늘기 때문입니다.", None),
        ("", None),
        ("2. 표고 보정", f_lbl),
        (f"   F {grad['F']*1000:.2f} nT/km, D {grad['D']*1000*3600:.2f} arcsec/km "
         "의 선형 기울기를 적용합니다(한반도 중심 기준 상수).", None),
        ("", None),
        ("3. 한계 — 반드시 읽어주십시오", f_lbl),
        ("   · External 층은 성과 정리 단계에만 적용: 외부장 일변동은 F 중앙값 37.5 nT,", None),
        ("     D 중앙값 0.112° 로 편각은 전형적인 하루 변동만으로도 KPI(0.1°)를 초과한다.", None),
        ("     관측소 4소 공간보간으로 편각을 보정했으나, 4소는 NOC 모드 분해에", None),
        ("     필요한 밀도에 못 미쳐 근사에 머문다.", None),
        ("   · Crustal 해상도 부족: 보유 격자는 2.8 km 간격 전국 컴필레이션으로,", None),
        ("     지표 점 잔차와의 상관이 낮습니다. 측선간격 250 m 원측선 자료는 존재하지", None),
        ("     않으므로 현 격자가 해상도 상한입니다(신규 자력측량 여부가 갈림길).", None),
        ("   · 지표 측점 부족 및 품질: 유효 측점이 권고치(30점)에 미달하며,", None),
        ("     일부 측점은 재방문 편각 산포가 영년변화로 설명 불가한 수준입니다.", None),
        ("", None),
        ("4. 좌표계와 높이 기준", f_lbl),
        ("   위경도는 WGS84 측지좌표를 입력합니다.", None),
        ("   높이는 '표고'(정표고, 평균해면 기준)를 입력합니다. 본 모델이 표고로 적합되었기 때문입니다.", None),
        ("   GNSS·RTK 수신값은 '타원체고(h)'이므로 지오이드고 N 을 빼서 표고로 바꿔야 합니다.", None),
        ("     H(표고) = h(타원체고) − N(지오이드고),  한반도 N ≈ 20~30 m (KNGeoid)", None),
        ("", None),
        ("   참고 — IGRF 자체는 타원체고를 요구하지만(ppigrf: height above ellipsoid),", None),
        ("   지오이드 오프셋은 공간적으로 완만해 Regional 다항식의 상수항에 흡수됩니다.", None),
        ("   실측: N 을 20~30 m 부여해도 교차검증 오차는 58.64 nT 로 변하지 않고", None),
        ("   Regional F 상수항만 약 0.7 nT 이동합니다. 높이 감도는 F 기준 −0.026 nT/m 이므로", None),
        ("   표고/타원체고 혼동(약 25 m)은 0.7 nT 수준으로 현재 모델 불확도(59 nT)에 비해 미미합니다.", None),
        ("", None),
        ("5. KPI 의 성격 — 법정 기준이 아님", f_lbl),
        ("   D<0.1°, F<50 nT 는 「한반도 LMM」 설명자료(오석훈, 강원대)가 제시한 공학적 목표치입니다.", None),
        ("", None),
        ("   법정 기준은 「지구물리측량 작업규정」(고시 제2021-2985호)에 있습니다.", None),
        ("     · 제20조 측정오차의 한계: 편각·복각 정수차 30'(=0.5°), 관측시간 20분 이내", None),
        ("     · 제17조: 총자기장 측정기 0.1 nT 단위 판독", None),
        ("     · 제19조: 1등 1일 6회 이상 측정, 편각·복각 측정 시 총자기장과 시간 동시 측정", None),
        ("     · 제21조: 상시 기준점과 비교하여 일변화·영년변화 보정 후 기준년 환산", None),
        ("   법정 한계(30')는 측정 품질관리 기준이고 KPI(0.1°)는 모델 예측정확도 목표로,", None),
        ("   성격이 다르며 KPI 가 5배 엄격합니다. 따라서 KPI 미달은 법령 위반이 아닙니다.", None),
        ("", None),
        ("   자편각 표기의 근거: 같은 규정 제12조제1항은 1등 지자기점의 설치 목적에", None),
        ("   '국가기본도의 자편각 표기'를 명시합니다. 다만 지형도 도식 규정상 자침편차 도표는", None),
        ("   난외사항 필수목록에 없어 표기 자체는 임의사항이며, 적용 축척도 특정되어 있지 않습니다.", None),
        ("", None),
        ("   다만 도폭 규격은 법정 사항이며, 1:25,000 을 기준으로 삼는 근거는 다음과 같습니다.", None),
        ("     · 도식적용규정 제4조: 1:25,000 도곽 = 경위도 7'30\" (0.125°)", None),
        ("     · 도폭당 자침편차를 단일 값으로 표기하려면 도폭 내 편각 변화가 표기 정밀도보다 작아야 함", None),
        ("     · 실측 계산: 도폭 내 D 변화가 1:50,000 은 0.091°(0.1° 예산의 91%),", None),
        ("       1:25,000 은 0.045°(45%). 모델 오차분을 남길 수 있는 가장 성긴 축척이 1:25,000", None),
    ]
    for k, (txt, fmt) in enumerate(lines):
        wsd.write(k, 0, txt, fmt if fmt else f_note)

    wb.close()
    print(f"[4/4] [저장] {out}  ({out.stat().st_size/1024/1024:.1f} MB)")
    return out


if __name__ == "__main__":
    main()
