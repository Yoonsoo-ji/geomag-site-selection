# -*- coding: utf-8 -*-
"""
한반도 LMM 웹 계산기 생성기
============================

lmm_build.py 가 만든 docs/data/lmm_model.json 과 IGRF14.shc 를 읽어
**단일 파일 오프라인 HTML 계산기**(docs/lmm.html)를 생성한다.

외부 네트워크·CDN 의존성이 전혀 없으므로 현장에서 파일만 열면 동작한다.

실행:
    python lmm_web.py
"""

import json
import datetime as dt
from pathlib import Path

import numpy as np
import ppigrf

BASE = Path(__file__).parent
MODEL_JSON = BASE / "docs" / "data" / "lmm_model.json"
COAST_JSON = BASE / "data" / "korea_boundary.geojson"
SHC = Path(ppigrf.__file__).parent / "IGRF14.shc"
OUT_HTML = BASE / "docs" / "lmm.html"

# 계산기가 지원할 Epoch 범위. 앞뒤 epoch 를 포함해 선형보간/외삽한다.
EPOCHS_WANTED = [2015.0, 2020.0, 2025.0, 2030.0]


def load_shc(path: Path, wanted):
    """IGRF .shc 파일에서 필요한 epoch 의 Gauss 계수만 추출."""
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    header = lines[0].split()
    nmin, nmax = int(header[0]), int(header[1])
    epochs = [float(x) for x in lines[1].split()]

    idx = [epochs.index(e) for e in wanted]

    coeffs = {}  # (n, m) -> [values per wanted epoch]
    for ln in lines[2:]:
        parts = ln.split()
        if len(parts) < 3:
            continue
        n, m = int(parts[0]), int(parts[1])
        vals = [float(parts[2 + i]) for i in idx]
        coeffs[(n, m)] = vals

    # g[n][m], h[n][m] 형태의 조밀 배열로 변환 (epoch 별)
    g = [[[0.0] * (nmax + 1) for _ in range(nmax + 1)] for _ in wanted]
    h = [[[0.0] * (nmax + 1) for _ in range(nmax + 1)] for _ in wanted]
    for (n, m), vals in coeffs.items():
        for k, v in enumerate(vals):
            if m >= 0:
                g[k][n][m] = v
            else:
                h[k][n][-m] = v

    return nmax, list(wanted), g, h


def load_coastline(path: Path, ndigits=3):
    """
    해안선 폴리곤을 [[lon,lat,...평탄화], ...] 형태의 링 목록으로 변환.

    좌표를 소수 3자리(약 100 m)로 반올림하고 중복 정점을 제거해
    인라인 용량을 줄인다. 지도 표시용이므로 이 정도 해상도로 충분하다.
    """
    gj = json.loads(path.read_text(encoding="utf-8"))
    rings = []

    def walk(coords, depth):
        # MultiPolygon -> Polygon -> ring -> [lon,lat]
        if depth == 0:
            flat, prev = [], None
            for lon, lat in coords:
                p = (round(lon, ndigits), round(lat, ndigits))
                if p != prev:
                    flat.extend(p)
                    prev = p
            if len(flat) >= 8:
                rings.append(flat)
            return
        for c in coords:
            walk(c, depth - 1)

    for feat in gj["features"]:
        geom = feat["geometry"]
        if geom["type"] == "Polygon":
            walk(geom["coordinates"], 1)
        elif geom["type"] == "MultiPolygon":
            walk(geom["coordinates"], 2)

    rings.sort(key=len, reverse=True)
    return rings


def build_html(model, nmax, epochs, g, h) -> str:
    igrf_data = {
        "nmax": nmax,
        "epochs": epochs,
        "g": g,
        "h": h,
    }

    coast = load_coastline(COAST_JSON)
    payload = json.dumps(
        {"igrf": igrf_data, "lmm": model, "coast": coast},
        ensure_ascii=False, separators=(",", ":"),
    )

    val = {(r["성분"], r["단계"]): r["RMS"] for r in model["validation"]}
    cv = model["loo_cv"]
    cd = model["crustal_diagnostics"]

    return _TEMPLATE.replace("/*__PAYLOAD__*/", payload).replace(
        "__SUMMARY__",
        json.dumps(
            {
                "generated": model["generated"],
                "d_loo": cv["D"],
                "i_loo": cv["I"],
                "f_loo": cv["F"],
                "crust_r": cd["corr"],
                "crust_n": cd["n"],
                "n_sites": len(model["sites"]),
                "epoch": model.get("epoch_label") or "—",
            },
            ensure_ascii=False,
        ),
    )


_TEMPLATE = r"""<!doctype html>
<html lang="ko">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>한반도 LMM 계산기 — Korea Local Magnetic Model</title>
<style>
  :root{
    --bg:#f6f7f9; --panel:#fff; --ink:#1b1f24; --muted:#5b6672;
    --line:#dfe3e8; --accent:#1f6feb; --warn:#b45309; --warnbg:#fff7ed;
    --ok:#0f7b3f;
  }
  @media (prefers-color-scheme:dark){
    :root{ --bg:#0f1216; --panel:#171b21; --ink:#e6e9ee; --muted:#98a2b0;
           --line:#2a313a; --accent:#4c8dff; --warn:#f0a355; --warnbg:#2a1f12; }
  }
  *{box-sizing:border-box}
  body{margin:0;background:var(--bg);color:var(--ink);
       font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI","Malgun Gothic",sans-serif}
  .wrap{max-width:1080px;margin:0 auto;padding:24px 18px 60px}
  h1{font-size:22px;margin:0 0 4px}
  .sub{color:var(--muted);font-size:13px;margin-bottom:20px}
  .card{background:var(--panel);border:1px solid var(--line);border-radius:10px;
        padding:18px;margin-bottom:16px}
  .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px}
  label{display:block;font-size:12px;color:var(--muted);margin-bottom:4px}
  .hint{display:inline-block;width:14px;height:14px;line-height:14px;text-align:center;
        border-radius:50%;background:var(--line);color:var(--ink);font-size:10px;
        cursor:help;vertical-align:1px}
  input,select{width:100%;padding:8px 10px;border:1px solid var(--line);
        border-radius:6px;background:var(--bg);color:var(--ink);font-size:14px}
  button{background:var(--accent);color:#fff;border:0;border-radius:6px;
        padding:10px 18px;font-size:14px;cursor:pointer;margin-top:12px}
  button.sec{background:transparent;color:var(--accent);border:1px solid var(--accent)}
  table{width:100%;border-collapse:collapse;font-size:14px}
  th,td{text-align:left;padding:7px 8px;border-bottom:1px solid var(--line)}
  th{color:var(--muted);font-weight:600;font-size:12px}
  td.num{text-align:right;font-variant-numeric:tabular-nums;font-weight:600}
  .layers td.num{font-weight:400;color:var(--muted)}
  .big{font-size:26px;font-weight:700}
  .note{background:var(--warnbg);border-left:3px solid var(--warn);
        padding:12px 14px;border-radius:0 6px 6px 0;font-size:13px;color:var(--ink)}
  .note b{color:var(--warn)}
  .scroll{overflow-x:auto}
  code{background:var(--bg);padding:1px 5px;border-radius:4px;font-size:12px}
  .foot{color:var(--muted);font-size:12px;margin-top:24px}

  /* ---- 분포도 ---- */
  .maphead{display:flex;justify-content:space-between;align-items:flex-start;
           gap:12px;flex-wrap:wrap;margin-bottom:14px}
  .mtitle{font-weight:700;font-size:15px}
  .msub{color:var(--muted);font-size:12px;margin-top:2px}
  .seg{display:flex;flex-wrap:wrap;gap:4px;background:var(--bg);
       border:1px solid var(--line);border-radius:8px;padding:3px}
  .seg button{margin:0;padding:6px 11px;font-size:12px;border-radius:6px;
       background:transparent;color:var(--muted);font-weight:600}
  .seg button[aria-pressed="true"]{background:var(--panel);color:var(--ink);
       box-shadow:0 1px 3px rgba(0,0,0,.13)}
  .maprow{display:flex;gap:20px;align-items:flex-start;flex-wrap:wrap}
  .mapwrap{position:relative;flex:1 1 340px;min-width:280px;max-width:520px}
  #map{width:100%;height:auto;display:block;border-radius:8px;cursor:crosshair;
       background:var(--bg)}
  .tip{position:absolute;pointer-events:none;opacity:0;transition:opacity .1s;
       background:var(--ink);color:var(--panel);font-size:11.5px;line-height:1.45;
       padding:6px 9px;border-radius:6px;white-space:nowrap;z-index:5;
       font-variant-numeric:tabular-nums}
  .legend{flex:0 1 210px;min-width:180px;font-size:12px}
  .lgt{font-weight:600;margin-bottom:8px}
  .ramp{margin-bottom:4px}
  #ramp{width:100%;height:14px;display:block;border-radius:3px}
  .ticks{display:flex;justify-content:space-between;color:var(--muted);
         font-size:11px;margin-top:3px;font-variant-numeric:tabular-nums}
  .lgkey{margin-top:14px;display:flex;flex-direction:column;gap:6px;
         color:var(--muted)}
  .lgkey .sw{display:inline-block;width:11px;height:11px;border-radius:50%;
       margin-right:7px;vertical-align:-1px}
  .sw.site{background:#eda100;border:1.5px solid var(--panel);
       box-shadow:0 0 0 1px #eda100}
  .sw.gap{border-radius:2px;background:repeating-linear-gradient(45deg,
       var(--muted) 0 1.5px,transparent 1.5px 4px)}
  .sw.iso{border-radius:2px;height:2px;background:var(--ink);opacity:.55;
       margin-bottom:3px}
  .lgnote{margin-top:14px;color:var(--muted);font-size:11.5px;line-height:1.5}
</style>
</head>
<body>
<div class="wrap">
  <h1>한반도 LMM 계산기</h1>
  <div class="sub">Korea Local Magnetic Model v1 — IGRF-14 + Regional + Crustal 결합 · 오프라인 단일 파일</div>

  <div class="card">
    <div class="grid">
      <div><label>위도 (°N)</label><input id="lat" value="37.5665"></div>
      <div><label>경도 (°E)</label><input id="lon" value="126.9780"></div>
      <div><label>표고 (m) <span class="hint" title="정표고(평균해면 기준). GNSS·RTK 의 타원체고를 쓰려면 지오이드고 N 을 빼십시오 (H = h − N, 한반도 N ≈ 20~30 m). 높이 감도는 F 기준 −0.026 nT/m 로, 25 m 혼동은 약 0.7 nT 입니다.">ⓘ</span></label><input id="elev" value="50"></div>
      <div><label>날짜</label><input id="date" type="date"></div>
    </div>
    <button onclick="run()">계산</button>
    <button class="sec" onclick="batch()">CSV 일괄계산</button>
    <input id="file" type="file" accept=".csv" style="display:none">
  </div>

  <div class="card" id="out" style="display:none">
    <div class="scroll"><table>
      <thead><tr><th>성분</th><th style="text-align:right">LMM 값</th><th style="text-align:right">IGRF-14 단독</th><th style="text-align:right">차이</th></tr></thead>
      <tbody id="res"></tbody>
    </table></div>
    <div style="margin-top:16px"><label>층별 기여 (총자력 F)</label>
    <div class="scroll"><table class="layers"><tbody id="lay"></tbody></table></div></div>
  </div>

  <div class="card">
    <div class="maphead">
      <div>
        <div class="mtitle">모델 분포도</div>
        <div class="msub" id="mapsub"></div>
      </div>
      <div class="seg" id="fieldsel"></div>
    </div>
    <div class="maprow">
      <div class="mapwrap">
        <canvas id="map"></canvas>
        <div id="tip" class="tip"></div>
      </div>
      <div class="legend">
        <div class="lgt" id="lgtitle"></div>
        <div class="ramp"><canvas id="ramp"></canvas>
          <div class="ticks" id="ticks"></div></div>
        <div class="lgkey">
          <div><span class="sw site"></span>지표 절대측정 측점</div>
          <div><span class="sw gap"></span>항공자력 자료 공백</div>
          <div><span class="sw iso"></span>등치선</div>
        </div>
        <div class="lgnote" id="lgnote"></div>
      </div>
    </div>
  </div>

  <div class="card">
    <div class="note" id="acc"></div>
  </div>

  <div class="foot" id="foot"></div>
</div>

<script>
const DATA = /*__PAYLOAD__*/;
const SUM  = __SUMMARY__;
const IG = DATA.igrf, LM = DATA.lmm;

// ---- WGS84 ----
const A_WGS = 6378.137, F_WGS = 1/298.257223563;
const B_WGS = A_WGS*(1-F_WGS);
const RE = 6371.2;                      // IGRF 기준 반경 (km)

/* Schmidt 준정규화 결합 르장드르 함수와 그 theta 미분.
   표준 IGRF/WMM 재귀식을 사용한다. */
function legendre(nmax, theta){
  const P = [], dP = [];
  for(let n=0;n<=nmax;n++){ P.push(new Float64Array(nmax+1)); dP.push(new Float64Array(nmax+1)); }
  const c = Math.cos(theta), s = Math.sin(theta);
  P[0][0] = 1; dP[0][0] = 0;
  for(let n=1;n<=nmax;n++){
    for(let m=0;m<=n;m++){
      if(n === m){
        const k = (n === 1) ? 1 : Math.sqrt((2*n-1)/(2*n));
        const Pp = P[n-1][n-1], dPp = dP[n-1][n-1];
        P[n][n]  = s*Pp*k;
        dP[n][n] = (c*Pp + s*dPp)*k;
      } else {
        const k1 = (2*n-1)/Math.sqrt(n*n - m*m);
        const k2 = (n-1 >= m) ? Math.sqrt((n-1)*(n-1) - m*m)/Math.sqrt(n*n - m*m) : 0;
        const Pm1 = P[n-1][m], dPm1 = dP[n-1][m];
        const Pm2 = (n-2 >= m) ? P[n-2][m] : 0;
        const dPm2 = (n-2 >= m) ? dP[n-2][m] : 0;
        P[n][m]  = k1*c*Pm1 - k2*Pm2;
        dP[n][m] = k1*(c*dPm1 - s*Pm1) - k2*dPm2;
      }
    }
  }
  return {P, dP};
}

/* epoch 간 Gauss 계수 선형보간 */
function interpCoef(year){
  const E = IG.epochs;
  let i = 0;
  while(i < E.length-2 && year > E[i+1]) i++;
  const t = (year - E[i])/(E[i+1]-E[i]);
  const g = [], h = [];
  for(let n=0;n<=IG.nmax;n++){
    g.push(new Float64Array(IG.nmax+1));
    h.push(new Float64Array(IG.nmax+1));
    for(let m=0;m<=n;m++){
      g[n][m] = IG.g[i][n][m] + t*(IG.g[i+1][n][m]-IG.g[i][n][m]);
      h[n][m] = IG.h[i][n][m] + t*(IG.h[i+1][n][m]-IG.h[i][n][m]);
    }
  }
  return {g, h};
}

/* IGRF-14 구면조화 합성 -> 측지좌표계 X(북) Y(동) Z(하) */
function igrf(latDeg, lonDeg, elevM, year){
  const {g, h} = interpCoef(year);
  const lat = latDeg*Math.PI/180, lon = lonDeg*Math.PI/180;
  const hkm = elevM/1000;

  // 측지 -> 지심 변환
  const sl = Math.sin(lat), cl = Math.cos(lat);
  const a2=A_WGS*A_WGS, b2=B_WGS*B_WGS;
  const rho = Math.sqrt(a2*cl*cl + b2*sl*sl);
  const r = Math.sqrt(hkm*hkm + 2*hkm*rho + (a2*a2*cl*cl + b2*b2*sl*sl)/(rho*rho));
  const cd = (hkm + rho)/r;                        // cos(psi)
  const sd = (a2 - b2)/rho * cl*sl / r;            // sin(psi)
  const slc = sl*cd - cl*sd;                       // 지심 sin(lat)
  const clc = cl*cd + sl*sd;
  const theta = Math.atan2(clc, slc);              // 여위도

  const {P, dP} = legendre(IG.nmax, theta);
  const st = Math.sin(theta);

  let Br=0, Bt=0, Bp=0;
  const ratio = RE/r;
  const cosm = new Float64Array(IG.nmax+1), sinm = new Float64Array(IG.nmax+1);
  for(let m=0;m<=IG.nmax;m++){ cosm[m]=Math.cos(m*lon); sinm[m]=Math.sin(m*lon); }

  for(let n=1;n<=IG.nmax;n++){
    const rr = Math.pow(ratio, n+2);
    for(let m=0;m<=n;m++){
      const gh = g[n][m]*cosm[m] + h[n][m]*sinm[m];
      const dgh = -g[n][m]*sinm[m] + h[n][m]*cosm[m];
      Br += rr*(n+1)*gh*P[n][m];
      Bt -= rr*gh*dP[n][m];
      if(st !== 0) Bp -= rr*m*dgh*P[n][m]/st;
      else if(m === 1) Bp -= rr*dgh*dP[n][m];
    }
  }

  // 지심 -> 측지 회전 (psi = 측지위도 - 지심위도)
  const X = -Bt, Y = Bp, Z = -Br;
  const Xg = X*cd + Z*sd;
  const Zg = Z*cd - X*sd;
  return {X:Xg, Y:Y, Z:Zg};
}

/* Regional 층: 정규화 좌표 다항식 */
function polyTerms(lat, lon, deg){
  const u = lat - LM.regional.lat0, v = lon - LM.regional.lon0;
  const t = [];
  for(let total=0; total<=deg; total++)
    for(let p=0; p<=total; p++)
      t.push(Math.pow(u, total-p)*Math.pow(v, p));
  return t;
}
function regional(lat, lon, key){
  const c = LM.regional[key], deg = LM.regional.degree;
  const t = polyTerms(lat, lon, deg);
  let s = 0;
  for(let i=0;i<c.length;i++) s += c[i]*t[i];
  return s;
}

/* Crustal 층: KIGAM 격자 쌍선형 보간.
   반환 status — 'ok' | 'outside'(격자 경계 밖) | 'gap'(격자 내부이나 자료 공백) */
function crustal(lat, lon){
  const C = LM.crustal;
  const fx = (lon - C.lon0)/C.dlon, fy = (lat - C.lat0)/C.dlat;
  if(fx < -0.5 || fy < -0.5 || fx > C.nlon-0.5 || fy > C.nlat-0.5)
    return {value:null, status:'outside'};

  const i0 = Math.min(Math.max(Math.floor(fx),0), C.nlon-2);
  const j0 = Math.min(Math.max(Math.floor(fy),0), C.nlat-2);
  const tx = Math.min(Math.max(fx-i0,0),1), ty = Math.min(Math.max(fy-j0,0),1);
  const at = (j,i) => C.values[j*C.nlon + i];
  const v = [at(j0,i0), at(j0,i0+1), at(j0+1,i0), at(j0+1,i0+1)];
  const w = [(1-tx)*(1-ty), tx*(1-ty), (1-tx)*ty, tx*ty];
  let sw=0, sv=0;
  for(let k=0;k<4;k++) if(v[k] !== null){ sw += w[k]; sv += w[k]*v[k]; }
  return sw > 0 ? {value:sv/sw, status:'ok'} : {value:null, status:'gap'};
}

/* Crustal 벡터 성분 (b_north, b_east, b_down).
   항공자력은 총자력 이상 dF 만 주지만, 이는 이상벡터를 주자기장 방향
   (l,m,n) 으로 투영한 값이라 파수영역 역산으로 3성분을 되찾을 수 있다.
   b_down 은 용량 절감을 위해 여기서 유도한다 (스칼라 격자와 정의상 정합):
     b_down = (dF + dc - l*b_north - m*b_east) / n
   자료 공백 화소는 역산에서 0 으로 채웠으므로 여기서도 0 으로 읽는다. */
function crustalVec(lat, lon){
  const C = LM.crustal, V = C && C.vector;
  if(!V) return null;
  const fx = (lon - C.lon0)/C.dlon, fy = (lat - C.lat0)/C.dlat;
  if(fx < -0.5 || fy < -0.5 || fx > C.nlon-0.5 || fy > C.nlat-0.5) return null;
  const i0 = Math.min(Math.max(Math.floor(fx),0), C.nlon-2);
  const j0 = Math.min(Math.max(Math.floor(fy),0), C.nlat-2);
  const tx = Math.min(Math.max(fx-i0,0),1), ty = Math.min(Math.max(fy-j0,0),1);
  const w = [(1-tx)*(1-ty), tx*(1-ty), (1-tx)*ty, tx*ty];
  const bil = a => {
    const v = [a[j0*C.nlon+i0], a[j0*C.nlon+i0+1],
               a[(j0+1)*C.nlon+i0], a[(j0+1)*C.nlon+i0+1]];
    let s = 0;
    for(let k=0;k<4;k++) s += w[k]*(v[k]===null?0:v[k]);
    return s;
  };
  const bn = bil(V.b_north), be = bil(V.b_east);
  const bd = (bil(C.values) + V.dc_nT - V.l*bn - V.m*be)/V.n;
  return {n:bn, e:be, d:bd};
}

/* 4-층 결합 */
function lmm(lat, lon, elev, year){
  const b = igrf(lat, lon, elev, year);
  const H0 = Math.hypot(b.X, b.Y);
  const F0 = Math.hypot(H0, b.Z);
  const D0 = Math.atan2(b.Y, b.X)*180/Math.PI;
  const I0 = Math.atan2(b.Z, H0)*180/Math.PI;

  const cr = crustal(lat, lon);
  const cv = crustalVec(lat, lon);
  const dD = regional(lat, lon, 'D');
  const dI = regional(lat, lon, 'I');
  const dF = regional(lat, lon, 'F');

  // 지각 벡터가 만드는 편각·복각 기여 (소각 근사 없이 벡터합 후 재측정)
  let crD = 0, crI = 0;
  if(cv){
    const Xc = b.X + cv.n, Yc = b.Y + cv.e, Zc = b.Z + cv.d;
    crD = Math.atan2(Yc, Xc)*180/Math.PI - D0;
    crI = Math.atan2(Zc, Math.hypot(Xc, Yc))*180/Math.PI - I0;
  }

  const D = D0 + crD + dD, I = I0 + crI + dI,
        F = F0 + (cr.value === null ? 0 : cr.value) + dF;
  const H = F*Math.cos(I*Math.PI/180);
  return {
    igrf:{D:D0,I:I0,F:F0,H:H0,X:b.X,Y:b.Y,Z:b.Z},
    lmm:{D,I,F,H,
         X:H*Math.cos(D*Math.PI/180), Y:H*Math.sin(D*Math.PI/180),
         Z:F*Math.sin(I*Math.PI/180)},
    layers:{core:F0, crustal:cr.value, crustalStatus:cr.status, regional:dF,
            crustalD:crD, crustalI:crI, regionalD:dD, regionalI:dI, external:null}
  };
}

// ---- UI ----
function dms(d){
  const s = d<0?'-':'', a=Math.abs(d);
  const deg=Math.floor(a), mn=Math.floor((a-deg)*60), sc=((a-deg)*60-mn)*60;
  return `${s}${deg}° ${mn}' ${sc.toFixed(1)}"`;
}
function yearFrac(dstr){
  const dt = new Date(dstr);
  const y = dt.getFullYear();
  const s = new Date(y,0,1), e = new Date(y+1,0,1);
  return y + (dt-s)/(e-s);
}

function run(){
  const lat=+document.getElementById('lat').value;
  const lon=+document.getElementById('lon').value;
  const el =+document.getElementById('elev').value;
  const yr = yearFrac(document.getElementById('date').value);
  const r = lmm(lat, lon, el, yr);

  const rows = [
    ['편각 D', dms(r.lmm.D), dms(r.igrf.D), (r.lmm.D-r.igrf.D).toFixed(3)+'°'],
    ['복각 I', dms(r.lmm.I), dms(r.igrf.I), (r.lmm.I-r.igrf.I).toFixed(3)+'°'],
    ['총자력 F', r.lmm.F.toFixed(1)+' nT', r.igrf.F.toFixed(1)+' nT', (r.lmm.F-r.igrf.F).toFixed(1)+' nT'],
    ['수평 H', r.lmm.H.toFixed(1)+' nT', r.igrf.H.toFixed(1)+' nT', (r.lmm.H-r.igrf.H).toFixed(1)+' nT'],
    ['북 X', r.lmm.X.toFixed(1)+' nT', r.igrf.X.toFixed(1)+' nT', (r.lmm.X-r.igrf.X).toFixed(1)+' nT'],
    ['동 Y', r.lmm.Y.toFixed(1)+' nT', r.igrf.Y.toFixed(1)+' nT', (r.lmm.Y-r.igrf.Y).toFixed(1)+' nT'],
    ['연직 Z', r.lmm.Z.toFixed(1)+' nT', r.igrf.Z.toFixed(1)+' nT', (r.lmm.Z-r.igrf.Z).toFixed(1)+' nT'],
  ];
  document.getElementById('res').innerHTML = rows.map(x=>
    `<tr><td>${x[0]}</td><td class="num">${x[1]}</td><td class="num">${x[2]}</td><td class="num">${x[3]}</td></tr>`).join('');

  const L = r.layers;
  document.getElementById('lay').innerHTML = [
    ['① Core (IGRF-14)', L.core.toFixed(1)+' nT'],
    ['② Regional (지표 절대측정 '+SUM.n_sites+'점)', (L.regional>=0?'+':'')+L.regional.toFixed(1)+' nT'],
    ['③ Crustal (KIGAM 항공자력)',
       (L.crustalStatus==='ok' ? (L.crustal>=0?'+':'')+L.crustal.toFixed(1)+' nT'
     : L.crustalStatus==='gap' ? '자료 공백 구역 — 0 적용 (항공자력 미측선)'
     : '격자 범위 밖 — 0 적용')
     + (L.crustalD ? `　D ${(L.crustalD*60>=0?'+':'')+(L.crustalD*60).toFixed(2)}′ ·`
                   + ` I ${(L.crustalI*60>=0?'+':'')+(L.crustalI*60).toFixed(2)}′` : '')],
    ['④ External (관측소 4소)',
     '예측 시점 미적용 — 정온시 기준값 모형 (성과 정리 단계에서 편각·복각에 반영)'],
  ].map(x=>`<tr><td>${x[0]}</td><td class="num">${x[1]}</td></tr>`).join('');

  document.getElementById('out').style.display='block';
}

function batch(){
  const f = document.getElementById('file');
  f.onchange = e => {
    const rd = new FileReader();
    rd.onload = ev => {
      const lines = ev.target.result.trim().split(/\r?\n/);
      const out = ['lat,lon,elev_m,date,D_deg,I_deg,F_nT,H_nT,X_nT,Y_nT,Z_nT'];
      for(let i=1;i<lines.length;i++){
        const p = lines[i].split(',');
        if(p.length<4) continue;
        const r = lmm(+p[0], +p[1], +p[2], yearFrac(p[3]));
        out.push([p[0],p[1],p[2],p[3],
          r.lmm.D.toFixed(4), r.lmm.I.toFixed(4), r.lmm.F.toFixed(1),
          r.lmm.H.toFixed(1), r.lmm.X.toFixed(1), r.lmm.Y.toFixed(1), r.lmm.Z.toFixed(1)].join(','));
      }
      const blob = new Blob([out.join('\n')], {type:'text/csv'});
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = 'lmm_result.csv';
      a.click();
    };
    rd.readAsText(f.files[0]);
  };
  f.click();
}

/* ==================== 모델 분포도 ====================
   IGRF 는 공간적으로 매끄러우므로 0.25° 조밀격자에서 한 번만 구면조화
   합성을 수행한 뒤 화소별로는 쌍선형 보간한다. Crustal 층만 원해상도
   (0.025°)로 조회하므로 실시간 렌더링이 가능하다. */

const MAP = {
  lon0:125.4, lon1:130.3, lat0:32.9, lat1:38.8,
  step:0.25, W:0, H:0, cache:null, cacheYear:null,
  field:'D', sites:LM.sites || []
};

/* 색상 — 순차형(크기)은 단일 색조 blue 램프, 발산형(부호)은 blue↔red + 회색 중간 */
const SEQ = ['#cde2fb','#b7d3f6','#9ec5f4','#86b6ef','#6da7ec','#5598e7',
             '#3987e5','#2a78d6','#256abf','#1c5cab','#184f95','#104281','#0d366b'];
const DIV_LO = ['#0d366b','#184f95','#256abf','#3987e5','#6da7ec','#9ec5f4'];
const DIV_HI = ['#f5b3b2','#ef8d8c','#e96b6a','#e34948','#c73433','#a52625'];
const isDark = () => matchMedia('(prefers-color-scheme:dark)').matches
                     && document.documentElement.dataset.theme !== 'light';
const DIV_MID = () => isDark() ? '#383835' : '#f0efec';

const FIELDS = {
  D:{label:'편각 D', unit:'°', kind:'seq', dec:2,
     note:'진북과 자북의 각도 차. 음수는 자북이 진북보다 서쪽에 있음을 뜻합니다.'},
  I:{label:'복각 I', unit:'°', kind:'seq', dec:2,
     note:'자기력선이 수평면과 이루는 각도.'},
  F:{label:'총자력 F', unit:'nT', kind:'seq', dec:0,
     note:'자기장 벡터의 크기.'},
  C:{label:'지각 자기이상', unit:'nT', kind:'div', dec:0,
     note:'KIGAM 항공자력 이상값(③ Crustal 층). 빗금은 항공자력 미측선 구역으로 0이 적용됩니다.'},
  R:{label:'LMM − IGRF', unit:'nT', kind:'div', dec:0,
     note:'IGRF-14 단독 대비 본 모델의 총자력 보정량. Regional + Crustal 층의 합입니다.'}
};

function hex2rgb(h){return [parseInt(h.slice(1,3),16),parseInt(h.slice(3,5),16),parseInt(h.slice(5,7),16)];}
function rampColor(t, kind){
  t = Math.max(0, Math.min(1, t));
  if(kind === 'seq'){
    const x = t*(SEQ.length-1), i = Math.min(Math.floor(x), SEQ.length-2);
    return mix(hex2rgb(SEQ[i]), hex2rgb(SEQ[i+1]), x-i);
  }
  // 발산형: 0.5 를 중립 회색으로 두고 양쪽 팔을 동일 단계로
  const mid = hex2rgb(DIV_MID());
  if(t < 0.5){
    const x = (0.5-t)/0.5*(DIV_LO.length-1);
    const i = Math.min(Math.floor(x), DIV_LO.length-2);
    const arm = mix(hex2rgb(DIV_LO[DIV_LO.length-1-i]),
                    hex2rgb(DIV_LO[Math.max(DIV_LO.length-2-i,0)]), x-i);
    return mix(mid, arm, Math.min((0.5-t)/0.5*1.15, 1));
  }
  const x = (t-0.5)/0.5*(DIV_HI.length-1);
  const i = Math.min(Math.floor(x), DIV_HI.length-2);
  const arm = mix(hex2rgb(DIV_HI[i]), hex2rgb(DIV_HI[i+1]), x-i);
  return mix(mid, arm, Math.min((t-0.5)/0.5*1.15, 1));
}
function mix(a,b,t){return [a[0]+(b[0]-a[0])*t, a[1]+(b[1]-a[1])*t, a[2]+(b[2]-a[2])*t];}

/* 0.25° 격자에서 IGRF D·I·F 를 미리 계산 (연도 변경 시에만 재계산) */
function igrfCache(year){
  if(MAP.cache && MAP.cacheYear === year) return MAP.cache;
  const nx = Math.round((MAP.lon1-MAP.lon0)/MAP.step)+1;
  const ny = Math.round((MAP.lat1-MAP.lat0)/MAP.step)+1;
  const D=new Float64Array(nx*ny), I=new Float64Array(nx*ny), F=new Float64Array(nx*ny);
  for(let j=0;j<ny;j++) for(let i=0;i<nx;i++){
    const b = igrf(MAP.lat0+j*MAP.step, MAP.lon0+i*MAP.step, 0, year);
    const H = Math.hypot(b.X,b.Y);
    const k = j*nx+i;
    D[k]=Math.atan2(b.Y,b.X)*180/Math.PI;
    I[k]=Math.atan2(b.Z,H)*180/Math.PI;
    F[k]=Math.hypot(H,b.Z);
  }
  MAP.cache={nx,ny,D,I,F}; MAP.cacheYear=year;
  return MAP.cache;
}
function sampleCache(c, arr, lat, lon){
  const fx=(lon-MAP.lon0)/MAP.step, fy=(lat-MAP.lat0)/MAP.step;
  const i0=Math.min(Math.max(Math.floor(fx),0),c.nx-2);
  const j0=Math.min(Math.max(Math.floor(fy),0),c.ny-2);
  const tx=fx-i0, ty=fy-j0;
  const a=arr[j0*c.nx+i0], b=arr[j0*c.nx+i0+1];
  const d=arr[(j0+1)*c.nx+i0], e=arr[(j0+1)*c.nx+i0+1];
  return a*(1-tx)*(1-ty)+b*tx*(1-ty)+d*(1-tx)*ty+e*tx*ty;
}

/* 화면에 표시할 스칼라장 값. null 이면 자료 공백 */
function fieldValue(f, c, lat, lon){
  if(f === 'C'){ const r=crustal(lat,lon); return r.status==='ok'?r.value:null; }
  if(f === 'R'){
    const r=crustal(lat,lon);
    return (r.status==='ok'?r.value:0) + regional(lat,lon,'F');
  }
  if(f === 'F'){
    const r=crustal(lat,lon);
    return sampleCache(c,c.F,lat,lon)+(r.status==='ok'?r.value:0)+regional(lat,lon,'F');
  }
  return sampleCache(c, f==='D'?c.D:c.I, lat, lon) + regional(lat, lon, f);
}

function projX(lon,W){return (lon-MAP.lon0)/(MAP.lon1-MAP.lon0)*W;}
function projY(lat,H){return (1-(lat-MAP.lat0)/(MAP.lat1-MAP.lat0))*H;}
function invLon(x,W){return MAP.lon0+x/W*(MAP.lon1-MAP.lon0);}
function invLat(y,H){return MAP.lat0+(1-y/H)*(MAP.lat1-MAP.lat0);}

function coastPath(ctx,W,H){
  ctx.beginPath();
  for(const ring of DATA.coast){
    for(let k=0;k<ring.length;k+=2){
      const x=projX(ring[k],W), y=projY(ring[k+1],H);
      if(k===0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
    }
    ctx.closePath();
  }
}

/* 육지 마스크.
   putImageData 는 clip 경로를 무시하므로(캔버스 사양) 화소 단위 마스크를
   따로 만들어 색범위 산정과 알파 처리에 사용한다. */
function landMask(W,H){
  if(MAP.mask && MAP.maskW===W && MAP.maskH===H) return MAP.mask;
  const off=document.createElement('canvas');
  off.width=W; off.height=H;
  const c=off.getContext('2d');
  c.fillStyle='#fff'; coastPath(c,W,H); c.fill('evenodd');
  const d=c.getImageData(0,0,W,H).data;
  const m=new Uint8Array(W*H);
  for(let i=0,k=0;i<d.length;i+=4,k++) m[k]=d[i+3]>128?1:0;
  MAP.mask=m; MAP.maskW=W; MAP.maskH=H;
  return m;
}

/* 극단값 몇 개가 색범위를 독점하지 않도록 분위수로 자른다 */
function robustRange(arr, lo=0.02, hi=0.98){
  const v=arr.filter(x=>isFinite(x)).sort((a,b)=>a-b);
  if(!v.length) return [0,1];
  return [v[Math.floor(v.length*lo)], v[Math.min(v.length-1, Math.floor(v.length*hi))]];
}

function drawMap(){
  const cv=document.getElementById('map'), ctx=cv.getContext('2d');
  const cssW=cv.clientWidth||420;
  const midLat=(MAP.lat0+MAP.lat1)/2;
  const aspect=((MAP.lon1-MAP.lon0)*Math.cos(midLat*Math.PI/180))/(MAP.lat1-MAP.lat0);
  const cssH=Math.round(cssW/aspect);
  const dpr=Math.min(devicePixelRatio||1, 2);
  cv.width=Math.round(cssW*dpr); cv.height=Math.round(cssH*dpr);
  cv.style.height=cssH+'px';
  const W=cv.width, H=cv.height;
  MAP.W=W; MAP.H=H;

  const yr=yearFrac(document.getElementById('date').value);
  const c=igrfCache(yr);
  const f=MAP.field, spec=FIELDS[f];

  // 1) 육지 화소에서만 값을 구해 색범위를 정한다
  ctx.clearRect(0,0,W,H);
  const mask=landMask(W,H);
  const img=ctx.createImageData(W,H);
  const vals=new Float64Array(W*H).fill(NaN);
  const land=[];
  for(let y=0;y<H;y++) for(let x=0;x<W;x++){
    const k=y*W+x;
    if(!mask[k]) continue;
    const v=fieldValue(f,c,invLat(y+0.5,H),invLon(x+0.5,W));
    vals[k]=(v===null?NaN:v);
    if(v!==null && isFinite(v)) land.push(v);
  }
  let [vmin,vmax]=robustRange(land);
  if(spec.kind==='div'){ const m=Math.max(Math.abs(vmin),Math.abs(vmax))||1; vmin=-m; vmax=m; }
  if(!(vmax>vmin)) vmax=vmin+1;
  MAP.vmin=vmin; MAP.vmax=vmax; MAP.vals=vals;

  const gapRGB=hex2rgb(isDark()?'#4a4a46':'#b9b7ae');
  const period=Math.max(4,Math.round(7*dpr)), on=Math.max(1,Math.round(2*dpr));
  for(let y=0;y<H;y++) for(let x=0;x<W;x++){
    const k=y*W+x, o=k*4;
    if(!mask[k]){ img.data[o+3]=0; continue; }   // 해상 = 투명
    const v=vals[k];
    if(isNaN(v)){
      // 자료 공백: 45° 빗금
      const hit=((x+y)%period)<on;
      img.data[o]=gapRGB[0]; img.data[o+1]=gapRGB[1]; img.data[o+2]=gapRGB[2];
      img.data[o+3]=hit?205:60;
    } else {
      const rgb=rampColor((v-vmin)/(vmax-vmin), spec.kind);
      img.data[o]=rgb[0]; img.data[o+1]=rgb[1]; img.data[o+2]=rgb[2]; img.data[o+3]=255;
    }
  }
  ctx.putImageData(img,0,0);

  // 2) 등치선
  drawContours(ctx,W,H,c,f,spec,dpr);

  // 3) 해안선
  ctx.save();
  ctx.strokeStyle=isDark()?'rgba(255,255,255,.5)':'rgba(11,11,11,.45)';
  ctx.lineWidth=1.1*dpr; coastPath(ctx,W,H); ctx.stroke(); ctx.restore();

  // 4) 측점
  ctx.save();
  for(const s of MAP.sites){
    const x=projX(s.lon,W), y=projY(s.lat,H);
    ctx.beginPath(); ctx.arc(x,y,3.4*dpr,0,7);
    ctx.fillStyle='#eda100'; ctx.fill();
    ctx.lineWidth=1.4*dpr;
    ctx.strokeStyle=isDark()?'#1a1a19':'#fcfcfb'; ctx.stroke();
  }
  ctx.restore();

  drawLegend(spec);
  document.getElementById('mapsub').textContent =
    `${document.getElementById('date').value} 기준 · 표고 0 m · 클릭하면 해당 지점이 계산됩니다`;
  document.getElementById('lgnote').textContent = spec.note;
}

/* 마칭스퀘어 등치선 */
function drawContours(ctx,W,H,c,f,spec,dpr){
  const nx=90, ny=Math.round(90/((MAP.lon1-MAP.lon0)/(MAP.lat1-MAP.lat0)));
  const gx=(MAP.lon1-MAP.lon0)/(nx-1), gy=(MAP.lat1-MAP.lat0)/(ny-1);
  const mask=landMask(W,H);
  const G=new Float64Array(nx*ny);
  for(let j=0;j<ny;j++) for(let i=0;i<nx;i++){
    const lat=MAP.lat0+j*gy, lon=MAP.lon0+i*gx;
    // 해상 격자점은 등치선 계산에서 제외 (외삽 등치선 방지)
    const px=Math.min(W-1,Math.max(0,projX(lon,W)|0));
    const py=Math.min(H-1,Math.max(0,projY(lat,H)|0));
    if(!mask[py*W+px]){ G[j*nx+i]=NaN; continue; }
    const v=fieldValue(f,c,lat,lon);
    G[j*nx+i]=(v===null?NaN:v);
  }
  const span=MAP.vmax-MAP.vmin;
  if(!(span>0)) return;
  const raw=span/7;
  const mag=Math.pow(10,Math.floor(Math.log10(raw)));
  const step=[1,2,2.5,5,10].map(m=>m*mag).find(s=>span/s<=9) || mag*10;
  const start=Math.ceil(MAP.vmin/step)*step;

  ctx.save(); coastPath(ctx,W,H); ctx.clip();
  ctx.lineWidth=1*dpr;
  ctx.strokeStyle=isDark()?'rgba(255,255,255,.55)':'rgba(11,11,11,.5)';
  ctx.font=`${11*dpr}px system-ui,sans-serif`;
  ctx.textAlign='center'; ctx.textBaseline='middle';

  for(let lv=start; lv<=MAP.vmax; lv+=step){
    const segs=[];
    for(let j=0;j<ny-1;j++) for(let i=0;i<nx-1;i++){
      const v=[G[j*nx+i],G[j*nx+i+1],G[(j+1)*nx+i+1],G[(j+1)*nx+i]];
      if(v.some(isNaN)) continue;
      const p=[[i,j],[i+1,j],[i+1,j+1],[i,j+1]];
      const pts=[];
      for(let e=0;e<4;e++){
        const a=v[e], b=v[(e+1)%4];
        if((a<lv)!==(b<lv)){
          const t=(lv-a)/(b-a);
          const pa=p[e], pb=p[(e+1)%4];
          pts.push([pa[0]+(pb[0]-pa[0])*t, pa[1]+(pb[1]-pa[1])*t]);
        }
      }
      if(pts.length===2) segs.push(pts);
    }
    ctx.beginPath();
    for(const s of segs){
      ctx.moveTo(projX(MAP.lon0+s[0][0]*gx,W), projY(MAP.lat0+s[0][1]*gy,H));
      ctx.lineTo(projX(MAP.lon0+s[1][0]*gx,W), projY(MAP.lat0+s[1][1]*gy,H));
    }
    ctx.stroke();

    // 대표 위치 1곳에만 라벨 (숫자 과밀 방지)
    if(segs.length){
      const s=segs[Math.floor(segs.length*0.5)];
      const x=projX(MAP.lon0+s[0][0]*gx,W), y=projY(MAP.lat0+s[0][1]*gy,H);
      const txt=lv.toFixed(spec.dec===0?0:1);
      const w=ctx.measureText(txt).width;
      ctx.fillStyle=isDark()?'rgba(26,26,25,.85)':'rgba(252,252,251,.85)';
      ctx.fillRect(x-w/2-3*dpr, y-7*dpr, w+6*dpr, 14*dpr);
      ctx.fillStyle=isDark()?'#fff':'#0b0b0b';
      ctx.fillText(txt,x,y);
    }
  }
  ctx.restore();
}

function drawLegend(spec){
  const cv=document.getElementById('ramp'), ctx=cv.getContext('2d');
  const dpr=Math.min(devicePixelRatio||1,2);
  cv.width=Math.round((cv.clientWidth||190)*dpr); cv.height=Math.round(14*dpr);
  for(let x=0;x<cv.width;x++){
    const rgb=rampColor(x/(cv.width-1), spec.kind);
    ctx.fillStyle=`rgb(${rgb[0]|0},${rgb[1]|0},${rgb[2]|0})`;
    ctx.fillRect(x,0,1,cv.height);
  }
  document.getElementById('lgtitle').textContent=`${spec.label} (${spec.unit})`;
  const fmt=v=>spec.dec===0?Math.round(v).toLocaleString():v.toFixed(spec.dec);
  document.getElementById('ticks').innerHTML=
    `<span>${fmt(MAP.vmin)}</span><span>${fmt((MAP.vmin+MAP.vmax)/2)}</span><span>${fmt(MAP.vmax)}</span>`;
}

// 필드 선택 버튼
const fs=document.getElementById('fieldsel');
fs.innerHTML=Object.entries(FIELDS).map(([k,v])=>
  `<button data-f="${k}" aria-pressed="${k==='D'}">${v.label}</button>`).join('');
fs.addEventListener('click', e=>{
  const b=e.target.closest('button'); if(!b) return;
  MAP.field=b.dataset.f;
  [...fs.querySelectorAll('button')].forEach(x=>
    x.setAttribute('aria-pressed', x===b));
  drawMap();
});

// 클릭 → 해당 지점 계산
document.getElementById('map').addEventListener('click', e=>{
  const cv=e.currentTarget, r=cv.getBoundingClientRect();
  const lon=invLon((e.clientX-r.left)/r.width*MAP.W, MAP.W);
  const lat=invLat((e.clientY-r.top)/r.height*MAP.H, MAP.H);
  document.getElementById('lat').value=lat.toFixed(4);
  document.getElementById('lon').value=lon.toFixed(4);
  run();
  document.getElementById('out').scrollIntoView({behavior:'smooth',block:'nearest'});
});

// 호버 → 툴팁
const tip=document.getElementById('tip');
document.getElementById('map').addEventListener('mousemove', e=>{
  const cv=e.currentTarget, r=cv.getBoundingClientRect();
  const px=(e.clientX-r.left)/r.width*MAP.W, py=(e.clientY-r.top)/r.height*MAP.H;
  const ix=Math.floor(px), iy=Math.floor(py);
  if(ix<0||iy<0||ix>=MAP.W||iy>=MAP.H){ tip.style.opacity=0; return; }

  const lon=invLon(px,MAP.W), lat=invLat(py,MAP.H);
  const k=iy*MAP.W+ix;
  const spec=FIELDS[MAP.field];
  const head=`${lat.toFixed(3)}°N ${lon.toFixed(3)}°E`;

  /* NaN 은 두 가지 원인이 있으므로 구분한다.
     · 해상 — 육지 마스크 밖. Regional 층이 육상 측점으로만 적합되어 적용 범위 밖
     · 자료 공백 — 육지이나 항공자력 미측선 구간 */
  const onLand=MAP.mask && MAP.maskW===MAP.W && MAP.maskH===MAP.H
               ? MAP.mask[k] : 1;
  const v=MAP.vals?MAP.vals[k]:NaN;

  tip.innerHTML = !onLand
    ? `${head}<br>해상 — 모델 적용 범위 밖`
    : isNaN(v)
      ? `${head}<br>항공자력 자료 공백`
      : `${head}<br>${spec.label} ${
          spec.dec===0?Math.round(v).toLocaleString():v.toFixed(spec.dec)} ${spec.unit}`;
  tip.style.opacity=1;
  const tw=tip.offsetWidth, th=tip.offsetHeight;
  tip.style.left=Math.min(Math.max(e.clientX-r.left+12,0), r.width-tw)+'px';
  tip.style.top=Math.max(e.clientY-r.top-th-8,0)+'px';
});
document.getElementById('map').addEventListener('mouseleave',()=>tip.style.opacity=0);

// 리사이즈·화면회전 시 툴팁은 즉시 감춘다. 이전 폭 기준으로 잡힌 left 가
// 남아 있으면 좁은 화면에서 문서 가로스크롤을 만든다.
let rz; addEventListener('resize',()=>{
  tip.style.opacity=0;
  clearTimeout(rz); rz=setTimeout(drawMap,180);
});
matchMedia('(prefers-color-scheme:dark)').addEventListener('change',drawMap);
document.getElementById('date').addEventListener('change',()=>{run();drawMap();});

document.getElementById('date').value = new Date().toISOString().slice(0,10);
document.getElementById('acc').innerHTML =
  `<b>정확도 고지</b> — 본 모델은 <b>4개 층 중 3개</b>만으로 구축된 잠정판(v1)입니다.<br>
   교차검증(LOO) 실측 오차: <code>D ${SUM.d_loo.toFixed(2)}°</code>
   <code>I ${SUM.i_loo.toFixed(2)}°</code> <code>F ${SUM.f_loo.toFixed(0)} nT</code><br>
   목표 KPI(D&lt;0.1°, F&lt;50 nT)에 아직 도달하지 못했습니다. 주요 원인:
   <ul style="margin:6px 0 0 18px;padding:0">
     <li><b>External 층 미적용</b> — CYG 1분 자료 자체는 INTERMAGNET HAPI 로
         확보 가능하다. 다만 ① 성과표에 관측<b>시각</b>이 없어 보정할 시점을 알 수 없고,
         ② 2019 야장으로 시험 적용한 결과 「변동이 전역 균일」이라는 단순 차감은
         정온시 실익이 없고 자기폭풍 시 오히려 악화되었다. NOC 공간투영이 필요하다</li>
     <li><b>Crustal 층 해상도 부족</b> — 보유 격자는 1.5분(약 2.8 km) 전국 컴필레이션.
         지표 점 잔차와의 상관 r=${SUM.crust_r.toFixed(2)} 이나
         표본이 ${SUM.crust_n}개뿐이라 측점 하나만 바뀌어도 크게 흔들린다.
         측선간격 250 m 원측선 자료는 존재하지 않아 현 격자가 해상도 상한</li>
     <li><b>지표 측점 부족</b> — 유효 ${SUM.n_sites}개 측점 (설명자료 권고 30점)</li>
   </ul>
   <div style="margin-top:8px">따라서 본 계산기는 <b>IGRF-14 대비 개선된 참고값</b> 및
   계산 체계 검증용으로 사용하고, <b>지형도 자침편차 등 정식 편각 산출에는 사용하지 마십시오.</b></div>
   <div style="margin-top:8px;font-size:12px;opacity:.85">
   ※ 위 KPI(D&lt;0.1°, F&lt;50 nT)는 「한반도 LMM」 설명자료(오석훈, 강원대)가 제시한
   <b>공학적 목표치</b>이며 법정 기준이 아닙니다.
   「지구물리측량 작업규정」(국토지리정보원고시 제2021-2985호) 제20조의 법정 측정오차 한계는
   <b>편각·복각 정수차 30′(=0.5°), 관측시간 20분 이내</b>로, KPI 보다 5배 느슨합니다.
   다만 그것은 측정 품질관리 기준이고 KPI 는 모델 예측정확도 목표이므로 성격이 다릅니다.
   같은 규정 제12조제1항은 「국가기본도의 자편각 표기」를 1등 지자기점의 법정 설치목적으로
   정하고 있습니다. 상세는 구축현황보고서 7장 참조.</div>`;
document.getElementById('foot').textContent =
  `모델 생성 ${SUM.generated} · IGRF-14 (deg 13) · 지표 절대측정 ${SUM.epoch} · KIGAM 자력이상도 1982–2018`;

run();
drawMap();
</script>
</body>
</html>
"""


def main():
    model = json.loads(MODEL_JSON.read_text(encoding="utf-8"))
    nmax, epochs, g, h = load_shc(SHC, EPOCHS_WANTED)
    print(f"[IGRF] degree {nmax}, epochs {epochs}")

    html = build_html(model, nmax, epochs, g, h)
    OUT_HTML.write_text(html, encoding="utf-8")
    print(f"[저장] {OUT_HTML}  ({OUT_HTML.stat().st_size/1024:.0f} KB, 외부 의존성 없음)")


if __name__ == "__main__":
    main()
