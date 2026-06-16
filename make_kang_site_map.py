"""
20260616_kang_site_v1.kml → docs/kang_site_map.html
도상 선점 후보지 + 기존 측정점 + P1 후보 + 제외구역 통합 지도
2D(Leaflet/Folium) + 3D(MapLibre) 토글 지원
"""
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import folium

KML_PATH = Path(__file__).parent / "docs/output/20260616_kang_site_v1.kml"
OUT_PATH  = Path(__file__).parent / "docs/kang_site_map.html"

# AWS 공개 지형 타일 — API 키 불필요, GitHub Pages에서 바로 동작
# Terrarium 인코딩: R*256 + G + B/256 - 32768 = 고도(m)
TERRAIN_URL      = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"
TERRAIN_ENCODING = "terrarium"


def strip_cdata_html(raw: str) -> str:
    if not raw:
        return ""
    s = re.sub(r"<br\s*/?>", "\n", raw, flags=re.I)
    s = re.sub(r"</div>",    "\n", s,   flags=re.I)
    s = re.sub(r"<[^>]+>",   "",   s)
    s = s.replace("&nbsp;", " ").replace("&amp;", "&")
    return "\n".join(ln.strip() for ln in s.splitlines() if ln.strip())


def classify(name: str, desc: str):
    combined = (name + " " + desc).lower()
    if "기존점"  in combined:                          return "orange",    "star",   "기존점"
    if "간선임도" in combined:                          return "green",     "tree",   "간선임도"
    if "작업임도" in combined:                          return "darkgreen", "tree",   "작업임도"
    if "임도 아님" in combined or "임도아님" in combined: return "blue",      "info",   "임도 외"
    return "gray", "circle", "기타"


def parse_kml():
    tree = ET.parse(KML_PATH)
    root = tree.getroot()
    points = []
    for pm in root.iter("{http://www.opengis.net/kml/2.2}Placemark"):
        name_el  = pm.find("{http://www.opengis.net/kml/2.2}name")
        desc_el  = pm.find("{http://www.opengis.net/kml/2.2}description")
        coord_el = pm.find(".//{http://www.opengis.net/kml/2.2}coordinates")
        if coord_el is None:
            continue
        coords = coord_el.text.strip().split(",")
        if len(coords) < 2:
            continue
        lon, lat = float(coords[0]), float(coords[1])
        name     = name_el.text.strip()  if name_el  is not None and name_el.text  else "이름 없음"
        desc_raw = desc_el.text.strip()  if desc_el  is not None and desc_el.text  else ""
        desc     = strip_cdata_html(desc_raw)
        points.append({"name": name, "desc": desc, "lat": lat, "lon": lon})
    return points


def popup_html(p):
    BADGE = {"기존점": "#FF8C00", "간선임도": "#2E8B57",
             "작업임도": "#5B8C5A", "임도 외": "#4488CC", "기타": "#888"}
    _, _, label = classify(p["name"], p["desc"])
    color   = BADGE[label]
    desc_h  = p["desc"].replace("\n", "<br>") if p["desc"] else "—"
    lat, lon = p["lat"], p["lon"]
    enc      = p["name"].replace(" ", "%20")
    kakao    = f"https://map.kakao.com/link/map/{enc},{lat:.6f},{lon:.6f}"
    naver    = f"https://map.naver.com/v5/search/{lat:.6f},{lon:.6f}"
    return f"""<div style="font-family:'Malgun Gothic',sans-serif;font-size:13px;min-width:240px;max-width:320px;">
<b style="color:{color};">{p['name']}</b>
&nbsp;<span style="background:{color};color:white;padding:1px 6px;border-radius:3px;font-size:11px;">{label}</span>
<hr style="margin:5px 0;">
<b>위도:</b> {lat:.5f}°&thinsp;N &nbsp; <b>경도:</b> {lon:.5f}°&thinsp;E<br>
<hr style="margin:5px 0;border-color:#ddd;">{desc_h}
<hr style="margin:5px 0;border-color:#ddd;">
<div style="display:flex;gap:6px;margin-top:2px;">
  <a href="{kakao}" target="_blank"
     style="flex:1;text-align:center;background:#FFCD00;color:#3A1D1D;text-decoration:none;
            font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 카카오맵</a>
  <a href="{naver}" target="_blank"
     style="flex:1;text-align:center;background:#03C75A;color:white;text-decoration:none;
            font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 네이버지도</a>
</div></div>"""


# ─────────────────────────────────────────────────────────────────────────────
# JS 블록 — 제외구역 + P1 + 기존 측정점 (FeatureGroup 변수명은 나중에 치환)
# ─────────────────────────────────────────────────────────────────────────────
JS_LAYERS = """\
<script>
(function(){
  /* 고압철탑·송전탑 */
  fetch('data/zone_power.geojson').then(r=>r.json()).then(d=>{
    L.geoJSON(d,{
      style:()=>({fillColor:'#FF3300',color:'#CC0000',weight:1,fillOpacity:0.28}),
      onEachFeature:(f,l)=>l.bindTooltip('⚡ 고압철탑·송전탑 제외 (1.0 km)')
    }).addTo(%%FG_POWER%%);
  }).catch(e=>console.warn('zone_power:',e));

  /* 철도 */
  fetch('data/zone_railway.geojson').then(r=>r.json()).then(d=>{
    L.geoJSON(d,{
      style:()=>({fillColor:'#FF7700',color:'#CC4400',weight:1,fillOpacity:0.28}),
      onEachFeature:(f,l)=>l.bindTooltip('🚆 철도 제외 (5.0 km)')
    }).addTo(%%FG_RAIL%%);
  }).catch(e=>console.warn('zone_railway:',e));

  /* P1 후보 */
  fetch('data/candidates_p1.geojson').then(r=>r.json()).then(d=>{
    L.geoJSON(d,{
      pointToLayer:(f,ll)=>L.circleMarker(ll,{radius:7,color:'#CC0000',fillColor:'#FF2200',fillOpacity:0.85,weight:1.5}),
      onEachFeature(f,l){
        const p=f.properties;
        const vn=x=>(x==null||isNaN(+x))?'-':(+x).toFixed(1);
        const lat=(p.lat!=null)?p.lat:f.geometry.coordinates[1];
        const lon=(p.lon!=null)?p.lon:f.geometry.coordinates[0];
        const kakao='https://map.kakao.com/link/map/P1%ED%9B%84%EB%B3%B4'+p.idx+','+lat.toFixed(6)+','+lon.toFixed(6);
        const naver='https://map.naver.com/v5/search/'+lat.toFixed(6)+','+lon.toFixed(6);
        l.bindPopup(`<div style="font-family:Malgun Gothic,sans-serif;font-size:12.5px;min-width:240px;">
          <b style="color:#CC0000;">후보지 #${p.idx}</b>
          <span style="background:#FF2200;color:white;padding:1px 5px;border-radius:3px;font-size:11px;margin-left:4px;">P1 최우선</span><br>
          <hr style="margin:4px 0;">
          <b>위도:</b> ${lat.toFixed(5)}° N &nbsp; <b>경도:</b> ${lon.toFixed(5)}° E<br>
          <b>입지점수:</b> ${vn(p.score)}/100<br>
          <hr style="margin:4px 0;border-color:#ddd;">
          <div style="display:flex;gap:6px;margin-top:2px;">
            <a href="${kakao}" target="_blank" style="flex:1;text-align:center;background:#FFCD00;color:#3A1D1D;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 카카오맵</a>
            <a href="${naver}" target="_blank" style="flex:1;text-align:center;background:#03C75A;color:white;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 네이버지도</a>
          </div></div>`,{maxWidth:300});
        l.bindTooltip('P1 #'+p.idx+' ('+lat.toFixed(4)+'°N)');
      }
    }).addTo(%%FG_P1%%);
  }).catch(e=>console.warn('candidates_p1:',e));

  /* 기존 측정점 */
  fetch('data/existing_sites.geojson').then(r=>r.json()).then(d=>{
    L.geoJSON(d,{
      pointToLayer:(f,ll)=>{
        const icon=L.divIcon({html:'<div style="font-size:22px;line-height:1;text-shadow:1px 1px 2px rgba(0,0,0,0.5);">⭐</div>',className:'',iconAnchor:[11,11]});
        return L.marker(ll,{icon});
      },
      onEachFeature(f,l){
        const p=f.properties;
        const vf=(x,fmt)=>(x==null)?'-':fmt(x);
        const lat=p.lat, lon=p.lon;
        const kakao='https://map.kakao.com/link/map/'+encodeURIComponent(p.name)+','+lat.toFixed(6)+','+lon.toFixed(6);
        const naver='https://map.naver.com/v5/search/'+lat.toFixed(6)+','+lon.toFixed(6);
        l.bindPopup(`<div style="font-family:Malgun Gothic,sans-serif;font-size:12.5px;min-width:240px;">
          <b style="color:#8B4513;">⭐ ${p.name}</b><br>
          <hr style="margin:4px 0;">
          <b>위도:</b> ${lat.toFixed(5)}° N &nbsp; <b>경도:</b> ${lon.toFixed(5)}° E<br>
          <b>표고:</b> ${vf(p.elev,x=>x.toFixed(1)+' m')} &nbsp; <b>최신관측:</b> ${p.obs_year}년<br>
          <b>총자력:</b> ${vf(p.total,x=>x.toLocaleString('ko-KR',{minimumFractionDigits:1,maximumFractionDigits:1})+' nT')}<br>
          <b>편각:</b> ${vf(p.decl,x=>x.toFixed(2)+'°')} &nbsp; <b>복각:</b> ${vf(p.incl,x=>x.toFixed(2)+'°')}<br>
          <hr style="margin:4px 0;border-color:#ddd;">
          <div style="display:flex;gap:6px;margin-top:2px;">
            <a href="${kakao}" target="_blank" style="flex:1;text-align:center;background:#FFCD00;color:#3A1D1D;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 카카오맵</a>
            <a href="${naver}" target="_blank" style="flex:1;text-align:center;background:#03C75A;color:white;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 네이버지도</a>
          </div></div>`,{maxWidth:300});
        l.bindTooltip('⭐ '+p.name+' ('+p.obs_year+'년)');
      }
    }).addTo(%%FG_EXISTING%%);
  }).catch(e=>console.warn('existing_sites:',e));
})();
</script>
"""

# ─────────────────────────────────────────────────────────────────────────────
# 3D 토글 + MapLibre (Leaflet map_id는 나중에 치환)
# ─────────────────────────────────────────────────────────────────────────────
JS_3D = """\
<!-- MapLibre GL CSS/JS -->
<link href="https://unpkg.com/maplibre-gl@3/dist/maplibre-gl.css" rel="stylesheet">
<script src="https://unpkg.com/maplibre-gl@3/dist/maplibre-gl.js"></script>

<!-- 3D 뷰 컨테이너 -->
<div id="maplibre-container" style="display:none;position:fixed;top:0;left:0;width:100%;height:100%;z-index:500;"></div>

<!-- 2D/3D 토글 버튼 -->
<div id="toggle3d-btn" style="
    position:fixed; top:12px; left:50%; transform:translateX(-50%);
    z-index:10000; background:#1a73e8; color:white;
    padding:8px 22px; border-radius:20px; border:none;
    font-family:'Malgun Gothic',sans-serif; font-size:14px; font-weight:bold;
    cursor:pointer; box-shadow:0 2px 8px rgba(0,0,0,0.35);
    user-select:none;" onclick="toggle3D()">
  🏔 3D 지형 보기
</div>

<!-- 3D 패널 (3D 모드일 때만 표시) -->
<div id="panel3d" style="display:none;position:fixed;top:60px;left:12px;z-index:10001;
    background:rgba(10,10,30,0.88);color:#e0e0f0;padding:13px 15px;border-radius:10px;
    font-family:'Malgun Gothic',sans-serif;font-size:12.5px;line-height:1.75;
    border:1px solid rgba(255,255,255,0.15);min-width:210px;">
  <b style="font-size:13px;color:#a8e6a3;">🏔 3D 뷰 설정</b><br>
  <hr style="margin:7px 0;border-color:rgba(255,255,255,0.15);">
  <label style="display:flex;align-items:center;gap:8px;margin-bottom:6px;">
    지형 과장
    <input type="range" id="exag3d" min="1" max="20" step="0.5" value="4"
           style="width:100px;accent-color:#a8e6a3;"
           oninput="setExag(this.value)">
    <span id="exagVal" style="min-width:28px;">4x</span>
  </label>
  <label style="display:flex;align-items:center;gap:8px;margin-bottom:6px;">
    <input type="checkbox" id="hillshade3d" onchange="toggleHillshade(this.checked)"> 힐쉐이드
  </label>
  <label style="display:flex;align-items:center;gap:8px;">
    배경지도
    <select id="basemap3d" onchange="changeBasemap(this.value)"
            style="background:#1a1a2e;color:#e0e0f0;border:1px solid #444;border-radius:4px;padding:2px 4px;">
      <option value="satellite">🛰 위성 (ESRI)</option>
      <option value="osm">🌍 OpenStreetMap</option>
      <option value="dark">🌑 다크</option>
    </select>
  </label>
  <hr style="margin:7px 0;border-color:rgba(255,255,255,0.15);">
  <div id="terrain-status" style="font-size:11px;color:#a8e6a3;">
    ⏳ 3D 지형 로딩 중…
  </div>
</div>

<script>
const BASEMAPS_3D = {
  satellite: 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
  osm:       'https://tile.openstreetmap.org/{z}/{x}/{y}.png',
  dark:      'https://tiles.stadiamaps.com/tiles/alidade_smooth_dark/{z}/{x}/{y}.png',
};
const TERRAIN_URL_3D = '%%TERRAIN_URL%%';

let mlMap = null;
let is3D  = false;
let exag3d = 4;

function makeStyle3D(basemap) {
  return {
    version: 8,
    sources: {
      basemap: { type:'raster', tiles:[BASEMAPS_3D[basemap]], tileSize:256 },
      terrain: { type:'raster-dem', tiles:[TERRAIN_URL_3D], tileSize:256, encoding:'terrarium' },
    },
    layers:[{id:'basemap',type:'raster',source:'basemap'}],
    terrain:{ source:'terrain', exaggeration: exag3d },
    sky:{
      'sky-color':'#199EF3','sky-horizon-blend':0.5,
      'horizon-color':'#ffffff','horizon-fog-blend':0.5,
    },
  };
}

function getLeafletCenter() {
  /* Folium이 생성한 map 객체 찾기 */
  const mapObj = window.%%LEAFLET_MAP_VAR%%;
  if (!mapObj) return { lng:127.5, lat:36.5, zoom:7 };
  const c = mapObj.getCenter();
  return { lng: c.lng, lat: c.lat, zoom: Math.min(mapObj.getZoom() + 3, 16) };
}

function toggle3D() {
  is3D = !is3D;
  const btn = document.getElementById('toggle3d-btn');
  const container = document.getElementById('maplibre-container');
  const panel = document.getElementById('panel3d');

  if (is3D) {
    // Leaflet 지도 숨기기
    document.getElementById('%%LEAFLET_DIV_ID%%').style.visibility = 'hidden';
    container.style.display = 'block';
    panel.style.display = 'block';
    btn.textContent = '🗺 2D 지도로 돌아가기';
    btn.style.background = '#444';

    const { lng, lat, zoom } = getLeafletCenter();

    if (!mlMap) {
      mlMap = new maplibregl.Map({
        container: 'maplibre-container',
        style: makeStyle3D('satellite'),
        center: [lng, lat],
        zoom: zoom,
        pitch: 55,
        bearing: -15,
        antialias: true,
        maxZoom: 18,
      });
      mlMap.addControl(new maplibregl.NavigationControl());
      mlMap.addControl(new maplibregl.ScaleControl({ unit:'metric' }));
      mlMap.on('error', e => {
        console.warn('MapLibre error:', e);
      });
      mlMap.on('load', () => {
        document.getElementById('terrain-status').innerHTML =
          '✅ 3D 지형 로드 완료';
        setTimeout(() => {
          document.getElementById('terrain-status').style.display = 'none';
        }, 3000);
      });
    } else {
      mlMap.jumpTo({ center:[lng, lat], zoom, pitch:55, bearing:-15 });
    }
  } else {
    // 2D 복귀
    document.getElementById('%%LEAFLET_DIV_ID%%').style.visibility = 'visible';
    container.style.display = 'none';
    panel.style.display = 'none';
    btn.textContent = '🏔 3D 지형 보기';
    btn.style.background = '#1a73e8';
  }
}

function setExag(v) {
  exag3d = parseFloat(v);
  document.getElementById('exagVal').textContent = exag3d + 'x';
  if (mlMap) mlMap.setTerrain({ source:'terrain', exaggeration: exag3d });
}

function toggleHillshade(on) {
  if (!mlMap) return;
  if (on) {
    if (!mlMap.getSource('terrain-hs'))
      mlMap.addSource('terrain-hs', { type:'raster-dem', tiles:[TERRAIN_URL_3D], tileSize:256, encoding:'terrarium' });
    if (!mlMap.getLayer('hillshade'))
      mlMap.addLayer({ id:'hillshade', type:'hillshade', source:'terrain-hs',
        paint:{ 'hillshade-exaggeration':0.6, 'hillshade-shadow-color':'#222' } });
  } else {
    if (mlMap.getLayer('hillshade')) mlMap.removeLayer('hillshade');
  }
}

function changeBasemap(v) {
  if (!mlMap) return;
  mlMap.setStyle(makeStyle3D(v));
  mlMap.once('style.load', () => {
    mlMap.setTerrain({ source:'terrain', exaggeration: exag3d });
  });
}
</script>
"""

# ─────────────────────────────────────────────────────────────────────────────
# 범례 패널 (좌측 하단 — 레이어 컨트롤과 겹치지 않도록)
# ─────────────────────────────────────────────────────────────────────────────
LEGEND_HTML = """\
<div style="
    position:fixed; bottom:24px; left:14px; width:230px;
    background:rgba(255,255,255,0.96); border:2px solid #555;
    z-index:9999; padding:11px 13px; border-radius:8px;
    font-family:'Malgun Gothic',sans-serif; font-size:12.5px; line-height:1.75;
    box-shadow:2px 2px 6px rgba(0,0,0,0.28);">
  <b style="font-size:13.5px;">🗺 도상 선점 후보지</b><br>
  <span style="font-size:11px;color:#666;">총 %%TOTAL%%개 지점 | 2026.06.16</span>
  <hr style="margin:6px 0;">
  <span style="color:#FF8C00;">⭐</span>&thinsp;기존점 &nbsp;<b>%%C_기존점%%</b>개<br>
  <span style="color:#2E8B57;">🌲</span>&thinsp;간선임도 &nbsp;<b>%%C_간선임도%%</b>개<br>
  <span style="color:#5B8C5A;">🌿</span>&thinsp;작업임도 &nbsp;<b>%%C_작업임도%%</b>개<br>
  <span style="color:#4488CC;">🔵</span>&thinsp;임도 외 &nbsp;<b>%%C_임도외%%</b>개<br>
  <span style="color:#888;">⚫</span>&thinsp;기타 &nbsp;<b>%%C_기타%%</b>개<br>
  <hr style="margin:6px 0;border-color:#ccc;">
  <a href="index.html" style="color:#1a73e8;text-decoration:none;font-weight:bold;">
    ← 입지 선정 분석 지도</a>
</div>
"""


# ─────────────────────────────────────────────────────────────────────────────
def make_map():
    points = parse_kml()
    print(f"파싱된 지점 수: {len(points)}")

    lats = [p["lat"] for p in points]
    lons = [p["lon"] for p in points]
    center = [sum(lats)/len(lats), sum(lons)/len(lons)]

    m = folium.Map(location=center, zoom_start=7, tiles=None)

    folium.TileLayer("https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}",
                     attr="Google Hybrid", name="🛰 구글 위성+도로", max_zoom=20).add_to(m)
    folium.TileLayer("https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}",
                     attr="Google Maps",   name="🗺 구글 지도",       max_zoom=20).add_to(m)
    folium.TileLayer("OpenStreetMap", name="🌍 OpenStreetMap").add_to(m)

    # FeatureGroup — 순서가 overlays 순서와 일치해야 함
    fg_power    = folium.FeatureGroup(name="⚡ 고압철탑·송전탑 제외 (1km)", show=False)
    fg_rail     = folium.FeatureGroup(name="🚆 철도 제외 (5km)",             show=False)
    fg_p1       = folium.FeatureGroup(name="🔴 P1 우선순위 1등급 후보지",   show=True)
    fg_existing = folium.FeatureGroup(name="⭐ 기존 측정점",                 show=True)
    fg_cat = {
        "기존점":   folium.FeatureGroup(name="🟠 도상선점 — 기존점",   show=True),
        "간선임도": folium.FeatureGroup(name="🟢 도상선점 — 간선임도", show=True),
        "작업임도": folium.FeatureGroup(name="🟩 도상선점 — 작업임도", show=True),
        "임도 외":  folium.FeatureGroup(name="🔵 도상선점 — 임도 외",  show=True),
        "기타":     folium.FeatureGroup(name="⚫ 도상선점 — 기타",     show=True),
    }
    for fg in [fg_power, fg_rail, fg_p1, fg_existing] + list(fg_cat.values()):
        m.add_child(fg)

    folium.LayerControl(collapsed=False).add_to(m)

    # 도상 선점 마커
    for p in points:
        mc, fi, label = classify(p["name"], p["desc"])
        folium.Marker(
            location=[p["lat"], p["lon"]],
            popup=folium.Popup(popup_html(p), max_width=340),
            tooltip=p["name"],
            icon=folium.Icon(color=mc, icon=fi, prefix="fa"),
        ).add_to(fg_cat[label])

    m.save(str(OUT_PATH))
    html = OUT_PATH.read_text(encoding="utf-8")

    # ── FeatureGroup 변수명 추출 (overlays 순서 그대로) ──────────────────────
    ov_m = re.search(r'overlays\s*:\s*\{([^}]+)\}', html, re.DOTALL)
    if not ov_m:
        print("ERROR: overlays 섹션 없음"); return len(points)
    fg_vars = re.findall(r'(feature_group_[a-f0-9]+)', ov_m.group(1))
    print("overlays 변수:", fg_vars)
    assert len(fg_vars) >= 4, "FeatureGroup 수 부족"
    fgv_power, fgv_rail, fgv_p1, fgv_existing = fg_vars[:4]

    # ── Leaflet map div ID 추출 ──────────────────────────────────────────────
    div_m = re.search(r'id="(map_[a-f0-9]+)"', html)
    map_div_id = div_m.group(1) if div_m else "map"

    # ── Leaflet map JS 변수명 추출 ───────────────────────────────────────────
    var_m = re.search(r'var\s+(map_[a-f0-9]+)\s*=\s*L\.map\(', html)
    map_var = var_m.group(1) if var_m else map_div_id

    print(f"  map div: {map_div_id},  map var: {map_var}")

    # ── JS 치환 ──────────────────────────────────────────────────────────────
    js_layers = (JS_LAYERS
                 .replace("%%FG_POWER%%",    fgv_power)
                 .replace("%%FG_RAIL%%",     fgv_rail)
                 .replace("%%FG_P1%%",       fgv_p1)
                 .replace("%%FG_EXISTING%%", fgv_existing))

    js_3d = (JS_3D
             .replace("%%TERRAIN_URL%%",    TERRAIN_URL)
             .replace("%%LEAFLET_MAP_VAR%%", map_var)
             .replace("%%LEAFLET_DIV_ID%%",  map_div_id))

    # ── 범례 치환 ────────────────────────────────────────────────────────────
    cat_counts = {}
    for p in points:
        _, _, label = classify(p["name"], p["desc"])
        cat_counts[label] = cat_counts.get(label, 0) + 1

    legend = (LEGEND_HTML
              .replace("%%TOTAL%%",      str(len(points)))
              .replace("%%C_기존점%%",   str(cat_counts.get("기존점", 0)))
              .replace("%%C_간선임도%%", str(cat_counts.get("간선임도", 0)))
              .replace("%%C_작업임도%%", str(cat_counts.get("작업임도", 0)))
              .replace("%%C_임도외%%",   str(cat_counts.get("임도 외", 0)))
              .replace("%%C_기타%%",     str(cat_counts.get("기타", 0))))

    # ── </body> 앞에 주입 ────────────────────────────────────────────────────
    html = html.replace("</body>", js_layers + js_3d + legend + "\n</body>")

    OUT_PATH.write_text(html, encoding="utf-8")
    print(f"저장 완료: {OUT_PATH}")
    return len(points)


if __name__ == "__main__":
    n = make_map()
    print(f"완료: {n}개 지점")
