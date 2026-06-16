"""
20260616_kang_site_v1.kml → docs/kang_site_map.html
2D(Leaflet) + 3D(MapLibre) 통합 지도
3D 모드에서도 도상선점·P1·기존측정점·제외구역 모두 표시
"""
import json, re, xml.etree.ElementTree as ET
from pathlib import Path
import folium

KML_PATH = Path(__file__).parent / "docs/output/20260616_kang_site_v1.kml"
OUT_PATH  = Path(__file__).parent / "docs/kang_site_map.html"

# AWS 공개 지형 타일 (API 키 불필요, GitHub Pages 동작)
TERRAIN_URL      = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"
TERRAIN_ENCODING = "terrarium"

# ── 헬퍼 ──────────────────────────────────────────────────────────────────
def strip_html(raw):
    if not raw: return ""
    s = re.sub(r"<br\s*/?>", "\n", raw, flags=re.I)
    s = re.sub(r"</div>",    "\n", s,   flags=re.I)
    s = re.sub(r"<[^>]+>",   "",   s)
    return "\n".join(l.strip() for l in s.replace("&nbsp;"," ").replace("&amp;","&").splitlines() if l.strip())

def classify(name, desc):
    c = (name + " " + desc).lower()
    if "기존점"  in c:                         return "orange",    "star",   "기존점",   "#FF8C00"
    if "간선임도" in c:                         return "green",     "tree",   "간선임도", "#2E8B57"
    if "작업임도" in c:                         return "darkgreen", "tree",   "작업임도", "#5B8C5A"
    if "임도 아님" in c or "임도아님" in c:    return "blue",      "info",   "임도 외",  "#4488CC"
    return "gray", "circle", "기타", "#888888"

def parse_kml():
    root = ET.parse(KML_PATH).getroot()
    pts  = []
    for pm in root.iter("{http://www.opengis.net/kml/2.2}Placemark"):
        ne = pm.find("{http://www.opengis.net/kml/2.2}name")
        de = pm.find("{http://www.opengis.net/kml/2.2}description")
        ce = pm.find(".//{http://www.opengis.net/kml/2.2}coordinates")
        if ce is None: continue
        coords = ce.text.strip().split(",")
        if len(coords) < 2: continue
        lon, lat = float(coords[0]), float(coords[1])
        name = ne.text.strip() if ne is not None and ne.text else "이름 없음"
        desc = strip_html(de.text.strip() if de is not None and de.text else "")
        _, _, label, color = classify(name, desc)
        pts.append({"name": name, "desc": desc, "lat": lat, "lon": lon,
                    "label": label, "color": color})
    return pts

def map_links(lat, lon, name):
    enc   = name.replace(" ", "%20")
    kakao = f"https://map.kakao.com/link/map/{enc},{lat:.6f},{lon:.6f}"
    naver = f"https://map.naver.com/v5/search/{lat:.6f},{lon:.6f}"
    return kakao, naver

def link_btns(kakao, naver):
    return (f'<div style="display:flex;gap:6px;margin-top:4px;">'
            f'<a href="{kakao}" target="_blank" style="flex:1;text-align:center;'
            f'background:#FFCD00;color:#3A1D1D;text-decoration:none;font-weight:bold;'
            f'font-size:12px;padding:5px 0;border-radius:4px;">🗺 카카오맵</a>'
            f'<a href="{naver}" target="_blank" style="flex:1;text-align:center;'
            f'background:#03C75A;color:white;text-decoration:none;font-weight:bold;'
            f'font-size:12px;padding:5px 0;border-radius:4px;">🗺 네이버지도</a>'
            f'</div>')

def popup_html(p):
    kakao, naver = map_links(p["lat"], p["lon"], p["name"])
    dh = p["desc"].replace("\n","<br>") if p["desc"] else "—"
    return (f'<div style="font-family:\'Malgun Gothic\',sans-serif;font-size:13px;'
            f'min-width:240px;max-width:320px;">'
            f'<b style="color:{p["color"]};">{p["name"]}</b>'
            f'&nbsp;<span style="background:{p["color"]};color:white;padding:1px 6px;'
            f'border-radius:3px;font-size:11px;">{p["label"]}</span>'
            f'<hr style="margin:5px 0;">'
            f'<b>위도:</b> {p["lat"]:.5f}°&thinsp;N &nbsp; <b>경도:</b> {p["lon"]:.5f}°&thinsp;E<br>'
            f'<hr style="margin:5px 0;border-color:#ddd;">{dh}'
            f'<hr style="margin:5px 0;border-color:#ddd;">'
            + link_btns(kakao, naver) + '</div>')


# ── Leaflet용 GeoJSON fetch JS ──────────────────────────────────────────────
LEAFLET_JS = """\
<script>
(function(){
  fetch('data/zone_power.geojson').then(r=>r.json()).then(d=>{
    L.geoJSON(d,{style:()=>({fillColor:'#FF3300',color:'#CC0000',weight:1,fillOpacity:0.28}),
      onEachFeature:(f,l)=>l.bindTooltip('⚡ 고압철탑·송전탑 제외 (1.0 km)')
    }).addTo(%%FG_POWER%%);
  }).catch(e=>console.warn(e));

  fetch('data/zone_railway.geojson').then(r=>r.json()).then(d=>{
    L.geoJSON(d,{style:()=>({fillColor:'#FF7700',color:'#CC4400',weight:1,fillOpacity:0.28}),
      onEachFeature:(f,l)=>l.bindTooltip('🚆 철도 제외 (5.0 km)')
    }).addTo(%%FG_RAIL%%);
  }).catch(e=>console.warn(e));

  fetch('data/candidates_p1.geojson').then(r=>r.json()).then(d=>{
    L.geoJSON(d,{
      pointToLayer:(f,ll)=>L.circleMarker(ll,{radius:7,color:'#CC0000',fillColor:'#FF2200',fillOpacity:0.85,weight:1.5}),
      onEachFeature(f,l){
        const p=f.properties, vn=x=>(x==null||isNaN(+x))?'-':(+x).toFixed(1);
        const lat=(p.lat!=null)?p.lat:f.geometry.coordinates[1];
        const lon=(p.lon!=null)?p.lon:f.geometry.coordinates[0];
        const kakao='https://map.kakao.com/link/map/P1%ED%9B%84%EB%B3%B4'+p.idx+','+lat.toFixed(6)+','+lon.toFixed(6);
        const naver='https://map.naver.com/v5/search/'+lat.toFixed(6)+','+lon.toFixed(6);
        l.bindPopup(`<div style="font-family:Malgun Gothic,sans-serif;font-size:12.5px;min-width:240px;">
          <b style="color:#CC0000;">후보지 #${p.idx}</b>
          <span style="background:#FF2200;color:white;padding:1px 5px;border-radius:3px;font-size:11px;margin-left:4px;">P1 최우선</span><br>
          <hr style="margin:4px 0;"><b>위도:</b> ${lat.toFixed(5)}° N &nbsp; <b>경도:</b> ${lon.toFixed(5)}° E<br>
          <b>입지점수:</b> ${vn(p.score)}/100
          <hr style="margin:4px 0;border-color:#ddd;">
          <div style="display:flex;gap:6px;">
            <a href="${kakao}" target="_blank" style="flex:1;text-align:center;background:#FFCD00;color:#3A1D1D;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 카카오맵</a>
            <a href="${naver}" target="_blank" style="flex:1;text-align:center;background:#03C75A;color:white;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 네이버지도</a>
          </div></div>`,{maxWidth:300});
        l.bindTooltip('P1 #'+p.idx);
      }
    }).addTo(%%FG_P1%%);
  }).catch(e=>console.warn(e));

  fetch('data/existing_sites.geojson').then(r=>r.json()).then(d=>{
    L.geoJSON(d,{
      pointToLayer:(f,ll)=>{
        const icon=L.divIcon({html:'<div style="font-size:22px;line-height:1;text-shadow:1px 1px 2px rgba(0,0,0,0.5);">⭐</div>',className:'',iconAnchor:[11,11]});
        return L.marker(ll,{icon});
      },
      onEachFeature(f,l){
        const p=f.properties, vf=(x,fmt)=>(x==null)?'-':fmt(x);
        const lat=p.lat, lon=p.lon;
        const kakao='https://map.kakao.com/link/map/'+encodeURIComponent(p.name)+','+lat.toFixed(6)+','+lon.toFixed(6);
        const naver='https://map.naver.com/v5/search/'+lat.toFixed(6)+','+lon.toFixed(6);
        l.bindPopup(`<div style="font-family:Malgun Gothic,sans-serif;font-size:12.5px;min-width:240px;">
          <b style="color:#8B4513;">⭐ ${p.name}</b><br>
          <hr style="margin:4px 0;"><b>위도:</b> ${lat.toFixed(5)}° N &nbsp; <b>경도:</b> ${lon.toFixed(5)}° E<br>
          <b>표고:</b> ${vf(p.elev,x=>x.toFixed(1)+' m')} &nbsp; <b>최신관측:</b> ${p.obs_year}년<br>
          <b>총자력:</b> ${vf(p.total,x=>x.toLocaleString('ko-KR',{minimumFractionDigits:1,maximumFractionDigits:1})+' nT')}<br>
          <b>편각:</b> ${vf(p.decl,x=>x.toFixed(2)+'°')} &nbsp; <b>복각:</b> ${vf(p.incl,x=>x.toFixed(2)+'°')}
          <hr style="margin:4px 0;border-color:#ddd;">
          <div style="display:flex;gap:6px;">
            <a href="${kakao}" target="_blank" style="flex:1;text-align:center;background:#FFCD00;color:#3A1D1D;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 카카오맵</a>
            <a href="${naver}" target="_blank" style="flex:1;text-align:center;background:#03C75A;color:white;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 네이버지도</a>
          </div></div>`,{maxWidth:300});
        l.bindTooltip('⭐ '+p.name);
      }
    }).addTo(%%FG_EXISTING%%);
  }).catch(e=>console.warn(e));
})();
</script>
"""


# ── 3D MapLibre 블록 ────────────────────────────────────────────────────────
# %%POINTS_GEOJSON%% = 도상선점 GeoJSON, %%TERRAIN_URL%% = DEM 타일 URL
MAPLIBRE_BLOCK = """\
<link href="https://unpkg.com/maplibre-gl@3/dist/maplibre-gl.css" rel="stylesheet">
<script src="https://unpkg.com/maplibre-gl@3/dist/maplibre-gl.js"></script>

<div id="ml-wrap" style="display:none;position:fixed;top:0;left:0;width:100%;height:100%;z-index:500;">
  <div id="ml-map" style="width:100%;height:100%;"></div>
</div>

<!-- 토글 버튼 -->
<div id="toggle3d-btn" onclick="toggle3D()" style="
    position:fixed;top:12px;left:50%;transform:translateX(-50%);
    z-index:10000;background:#1a73e8;color:white;
    padding:8px 24px;border-radius:20px;border:none;
    font-family:'Malgun Gothic',sans-serif;font-size:14px;font-weight:bold;
    cursor:pointer;box-shadow:0 2px 8px rgba(0,0,0,0.35);user-select:none;">
  🏔 3D 지형 보기
</div>

<!-- 3D 패널 -->
<div id="panel3d" style="display:none;position:fixed;top:56px;left:12px;z-index:10001;
    background:rgba(10,10,30,0.88);color:#e0e0f0;padding:13px 15px;border-radius:10px;
    font-family:'Malgun Gothic',sans-serif;font-size:12.5px;line-height:1.8;
    border:1px solid rgba(255,255,255,0.15);min-width:220px;">
  <b style="font-size:13px;color:#a8e6a3;">🏔 3D 뷰 설정</b>
  <hr style="margin:7px 0;border-color:rgba(255,255,255,0.15);">
  <label style="display:flex;align-items:center;gap:8px;margin-bottom:5px;">
    지형 과장&thinsp;
    <input type="range" id="exag3d" min="1" max="20" step="0.5" value="4"
      style="width:90px;accent-color:#a8e6a3;" oninput="setExag(this.value)">
    <span id="exagVal">4x</span>
  </label>
  <label style="display:flex;align-items:center;gap:8px;margin-bottom:5px;">
    <input type="checkbox" id="hs3d" onchange="toggleHillshade(this.checked)"> 힐쉐이드
  </label>
  <label style="display:flex;align-items:center;gap:8px;margin-bottom:8px;">배경지도&thinsp;
    <select id="bm3d" onchange="changeBasemap(this.value)"
      style="background:#1a1a2e;color:#e0e0f0;border:1px solid #444;border-radius:4px;padding:2px 4px;">
      <option value="satellite">🛰 위성 (ESRI)</option>
      <option value="osm">🌍 OpenStreetMap</option>
      <option value="dark">🌑 다크</option>
    </select>
  </label>
  <hr style="margin:6px 0;border-color:rgba(255,255,255,0.15);">
  <b style="font-size:12px;color:#ccc;">레이어</b><br>
  <label style="display:flex;align-items:center;gap:6px;margin-top:4px;">
    <input type="checkbox" id="l3d-dosung"   checked onchange="set3DLayer('dosung',   this.checked)">
    <span style="color:#4ECDC4;">●</span> 도상 선점
  </label>
  <label style="display:flex;align-items:center;gap:6px;">
    <input type="checkbox" id="l3d-p1"       checked onchange="set3DLayer('p1',       this.checked)">
    <span style="color:#FF2200;">●</span> P1 우선순위 후보
  </label>
  <label style="display:flex;align-items:center;gap:6px;">
    <input type="checkbox" id="l3d-existing" checked onchange="set3DLayer('existing', this.checked)">
    <span style="color:#FFD700;">★</span> 기존 측정점
  </label>
  <label style="display:flex;align-items:center;gap:6px;">
    <input type="checkbox" id="l3d-power"          onchange="set3DLayer('power',    this.checked)">
    <span style="color:#FF3300;">■</span> 고압철탑 제외
  </label>
  <label style="display:flex;align-items:center;gap:6px;">
    <input type="checkbox" id="l3d-rail"           onchange="set3DLayer('rail',     this.checked)">
    <span style="color:#FF7700;">■</span> 철도 제외
  </label>
  <div id="ml-status" style="margin-top:7px;font-size:11px;color:#a8e6a3;">⏳ 지형 로딩 중…</div>
</div>

<script>
// ── 도상 선점 데이터 (인라인 GeoJSON) ──────────────────────────────────────
const DOSUNG_GJ = %%POINTS_GEOJSON%%;

const BASEMAPS_3D = {
  satellite:'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
  osm:      'https://tile.openstreetmap.org/{z}/{x}/{y}.png',
  dark:     'https://tiles.stadiamaps.com/tiles/alidade_smooth_dark/{z}/{x}/{y}.png',
};
const TERRAIN_URL_3D = '%%TERRAIN_URL%%';

let mlMap = null, is3D = false, exag3d = 4;
// 로드 완료된 GeoJSON 캐시
const gjCache = {};

function makeStyle3D(bm){
  return {
    version:8,
    sources:{
      basemap:{ type:'raster', tiles:[BASEMAPS_3D[bm]], tileSize:256 },
      terrain:{ type:'raster-dem', tiles:[TERRAIN_URL_3D], tileSize:256, encoding:'terrarium' },
    },
    layers:[{id:'basemap',type:'raster',source:'basemap'}],
    terrain:{ source:'terrain', exaggeration:exag3d },
    sky:{'sky-color':'#199EF3','sky-horizon-blend':0.5,'horizon-color':'#ffffff','horizon-fog-blend':0.5},
  };
}

function getLeafletCenter(){
  const mo = window.%%MAP_VAR%%;
  if(!mo) return {lng:127.5,lat:36.5,zoom:10};
  const c = mo.getCenter();
  return {lng:c.lng, lat:c.lat, zoom:Math.min(mo.getZoom()+3, 17)};
}

function toggle3D(){
  is3D = !is3D;
  const btn = document.getElementById('toggle3d-btn');
  const wrap = document.getElementById('ml-wrap');
  const panel= document.getElementById('panel3d');
  document.getElementById('%%MAP_DIV%%').style.visibility = is3D ? 'hidden' : 'visible';
  wrap.style.display  = is3D ? 'block' : 'none';
  panel.style.display = is3D ? 'block' : 'none';
  btn.textContent     = is3D ? '🗺 2D 지도로 돌아가기' : '🏔 3D 지형 보기';
  btn.style.background= is3D ? '#555' : '#1a73e8';

  if(!is3D) return;
  const {lng,lat,zoom} = getLeafletCenter();

  if(!mlMap){
    mlMap = new maplibregl.Map({
      container:'ml-map', style:makeStyle3D('satellite'),
      center:[lng,lat], zoom:zoom, pitch:55, bearing:-15,
      antialias:true, maxZoom:18,
    });
    mlMap.addControl(new maplibregl.NavigationControl());
    mlMap.addControl(new maplibregl.ScaleControl({unit:'metric'}));
    mlMap.on('load', onMlLoad);
    mlMap.on('error', e=>console.warn('ML error:',e));
  } else {
    mlMap.jumpTo({center:[lng,lat], zoom, pitch:55, bearing:-15});
  }
}

// ── 마커 풀 (MapLibre Marker 재사용) ────────────────────────────────────────
const markerPool = {dosung:[], p1:[], existing:[]};

function clearMarkers(key){ markerPool[key].forEach(m=>m.remove()); markerPool[key]=[]; }

// ── 팝업 HTML ──────────────────────────────────────────────────────────────
function mkPopup(html){ return new maplibregl.Popup({maxWidth:'330px',closeButton:true}).setHTML(html); }

function dosungPopupHtml(f){
  const p=f.properties;
  const lat=f.geometry.coordinates[1], lon=f.geometry.coordinates[0];
  const kakao='https://map.kakao.com/link/map/'+encodeURIComponent(p.name)+','+lat.toFixed(6)+','+lon.toFixed(6);
  const naver='https://map.naver.com/v5/search/'+lat.toFixed(6)+','+lon.toFixed(6);
  const dh=(p.desc||'—').replace(/\\n/g,'<br>');
  return `<div style="font-family:Malgun Gothic,sans-serif;font-size:12.5px;min-width:230px;max-width:310px;">
    <b style="color:${p.color};">${p.name}</b>
    <span style="background:${p.color};color:white;padding:1px 5px;border-radius:3px;font-size:11px;margin-left:4px;">${p.label}</span>
    <hr style="margin:4px 0;"><b>위도:</b> ${lat.toFixed(5)}°N &nbsp;<b>경도:</b> ${lon.toFixed(5)}°E
    <hr style="margin:4px 0;border-color:#ddd;">${dh}
    <hr style="margin:4px 0;border-color:#ddd;">
    <div style="display:flex;gap:5px;">
      <a href="${kakao}" target="_blank" style="flex:1;text-align:center;background:#FFCD00;color:#3A1D1D;text-decoration:none;font-weight:bold;font-size:12px;padding:4px 0;border-radius:4px;">🗺 카카오맵</a>
      <a href="${naver}" target="_blank" style="flex:1;text-align:center;background:#03C75A;color:white;text-decoration:none;font-weight:bold;font-size:12px;padding:4px 0;border-radius:4px;">🗺 네이버지도</a>
    </div></div>`;
}

function p1PopupHtml(f){
  const p=f.properties, vn=x=>(x==null||isNaN(+x))?'-':(+x).toFixed(1);
  const lat=(p.lat!=null)?p.lat:f.geometry.coordinates[1];
  const lon=(p.lon!=null)?p.lon:f.geometry.coordinates[0];
  const kakao='https://map.kakao.com/link/map/P1%ED%9B%84%EB%B3%B4'+p.idx+','+lat.toFixed(6)+','+lon.toFixed(6);
  const naver='https://map.naver.com/v5/search/'+lat.toFixed(6)+','+lon.toFixed(6);
  return `<div style="font-family:Malgun Gothic,sans-serif;font-size:12.5px;min-width:230px;">
    <b style="color:#CC0000;">후보지 #${p.idx}</b>
    <span style="background:#FF2200;color:white;padding:1px 5px;border-radius:3px;font-size:11px;margin-left:4px;">P1 최우선</span>
    <hr style="margin:4px 0;"><b>위도:</b> ${lat.toFixed(5)}°N &nbsp;<b>경도:</b> ${lon.toFixed(5)}°E<br>
    <b>입지점수:</b> ${vn(p.score)}/100
    <hr style="margin:4px 0;border-color:#ddd;">
    <div style="display:flex;gap:5px;">
      <a href="${kakao}" target="_blank" style="flex:1;text-align:center;background:#FFCD00;color:#3A1D1D;text-decoration:none;font-weight:bold;font-size:12px;padding:4px 0;border-radius:4px;">🗺 카카오맵</a>
      <a href="${naver}" target="_blank" style="flex:1;text-align:center;background:#03C75A;color:white;text-decoration:none;font-weight:bold;font-size:12px;padding:4px 0;border-radius:4px;">🗺 네이버지도</a>
    </div></div>`;
}

function existingPopupHtml(p, lat, lon){
  const vf=(x,fmt)=>(x==null)?'-':fmt(x);
  const kakao='https://map.kakao.com/link/map/'+encodeURIComponent(p.name)+','+lat.toFixed(6)+','+lon.toFixed(6);
  const naver='https://map.naver.com/v5/search/'+lat.toFixed(6)+','+lon.toFixed(6);
  return `<div style="font-family:Malgun Gothic,sans-serif;font-size:12.5px;min-width:230px;">
    <b style="color:#8B4513;">⭐ ${p.name}</b>
    <hr style="margin:4px 0;"><b>위도:</b> ${lat.toFixed(5)}°N &nbsp;<b>경도:</b> ${lon.toFixed(5)}°E<br>
    <b>표고:</b> ${vf(p.elev,x=>x.toFixed(1)+' m')} &nbsp;<b>최신관측:</b> ${p.obs_year}년<br>
    <b>총자력:</b> ${vf(p.total,x=>x.toLocaleString('ko-KR',{minimumFractionDigits:1,maximumFractionDigits:1})+' nT')}
    <hr style="margin:4px 0;border-color:#ddd;">
    <div style="display:flex;gap:5px;">
      <a href="${kakao}" target="_blank" style="flex:1;text-align:center;background:#FFCD00;color:#3A1D1D;text-decoration:none;font-weight:bold;font-size:12px;padding:4px 0;border-radius:4px;">🗺 카카오맵</a>
      <a href="${naver}" target="_blank" style="flex:1;text-align:center;background:#03C75A;color:white;text-decoration:none;font-weight:bold;font-size:12px;padding:4px 0;border-radius:4px;">🗺 네이버지도</a>
    </div></div>`;
}

// ── MapLibre load 후 레이어 추가 ────────────────────────────────────────────
function onMlLoad(){
  document.getElementById('ml-status').textContent = '✅ 3D 지형 로드 완료';
  setTimeout(()=>{ document.getElementById('ml-status').style.display='none'; }, 3000);

  // 도상 선점 마커
  addDosungMarkers();

  // P1 후보 — GeoJSON 소스+레이어
  fetchGj('data/candidates_p1.geojson','p1-src', gj=>{
    mlMap.addSource('p1-src',{type:'geojson',data:gj});
    mlMap.addLayer({id:'p1-circles',type:'circle',source:'p1-src',paint:{
      'circle-radius':7,'circle-color':'#FF2200','circle-stroke-color':'#CC0000','circle-stroke-width':1.5,'circle-opacity':0.9
    }});
    mlMap.on('click','p1-circles', e=>{
      const f=e.features[0];
      mkPopup(p1PopupHtml(f)).setLngLat(f.geometry.coordinates.slice()).addTo(mlMap);
    });
    mlMap.on('mouseenter','p1-circles',()=>mlMap.getCanvas().style.cursor='pointer');
    mlMap.on('mouseleave','p1-circles',()=>mlMap.getCanvas().style.cursor='');
  });

  // 기존 측정점 — HTML 마커
  fetchGj('data/existing_sites.geojson','existing', gj=>{
    gj.features.forEach(f=>{
      const p=f.properties;
      const lat=p.lat, lon=p.lon;
      const el=document.createElement('div');
      el.textContent='⭐';
      el.style.cssText='font-size:20px;line-height:1;cursor:pointer;text-shadow:0 0 3px rgba(0,0,0,0.8);';
      const m=new maplibregl.Marker({element:el})
        .setLngLat([lon,lat])
        .setPopup(mkPopup(existingPopupHtml(p,lat,lon)))
        .addTo(mlMap);
      markerPool.existing.push(m);
    });
  });

  // 제외구역 소스만 미리 로드 (레이어는 체크박스로 제어)
  fetchGj('data/zone_power.geojson','power-src', gj=>{
    mlMap.addSource('power-src',{type:'geojson',data:gj});
    // 기본 비표시
  });
  fetchGj('data/zone_railway.geojson','rail-src', gj=>{
    mlMap.addSource('rail-src',{type:'geojson',data:gj});
  });
}

function fetchGj(url, key, cb){
  if(gjCache[key]){ cb(gjCache[key]); return; }
  fetch(url).then(r=>r.json()).then(d=>{ gjCache[key]=d; cb(d); }).catch(e=>console.warn(url,e));
}

function addDosungMarkers(){
  clearMarkers('dosung');
  DOSUNG_GJ.features.forEach(f=>{
    const p=f.properties;
    const el=document.createElement('div');
    el.style.cssText=`width:12px;height:12px;border-radius:50%;background:${p.color};
      border:2px solid rgba(255,255,255,0.7);cursor:pointer;
      box-shadow:0 0 4px rgba(0,0,0,0.6);`;
    const m=new maplibregl.Marker({element:el, anchor:'center'})
      .setLngLat(f.geometry.coordinates)
      .setPopup(mkPopup(dosungPopupHtml(f)))
      .addTo(mlMap);
    markerPool.dosung.push(m);
  });
}

// ── 레이어 토글 ─────────────────────────────────────────────────────────────
const ZONE_LAYERS = {
  power:{ src:'power-src', lid:'power-fill', color:'#FF3300', stroke:'#CC0000', label:'고압철탑·송전탑 제외 (1.0 km)' },
  rail: { src:'rail-src',  lid:'rail-fill',  color:'#FF7700', stroke:'#CC4400', label:'철도 제외 (5.0 km)' },
};

function set3DLayer(key, on){
  if(!mlMap || !mlMap.isStyleLoaded()) return;

  if(key==='dosung'){
    markerPool.dosung.forEach(m=>{ m.getElement().style.display = on?'':'none'; });
    return;
  }
  if(key==='p1'){
    if(mlMap.getLayer('p1-circles')) mlMap.setLayoutProperty('p1-circles','visibility',on?'visible':'none');
    return;
  }
  if(key==='existing'){
    markerPool.existing.forEach(m=>{ m.getElement().style.display = on?'':'none'; });
    return;
  }
  // power / rail
  const cfg=ZONE_LAYERS[key]; if(!cfg) return;
  if(on){
    if(!mlMap.getSource(cfg.src)) return; // 아직 로드 안 됨
    if(!mlMap.getLayer(cfg.lid)){
      mlMap.addLayer({
        id:cfg.lid, type:'fill', source:cfg.src,
        paint:{'fill-color':cfg.color,'fill-opacity':0.28,'fill-outline-color':cfg.stroke}
      }, 'p1-circles'); // P1 레이어 아래에 삽입
    } else {
      mlMap.setLayoutProperty(cfg.lid,'visibility','visible');
    }
  } else {
    if(mlMap.getLayer(cfg.lid)) mlMap.setLayoutProperty(cfg.lid,'visibility','none');
  }
}

function setExag(v){
  exag3d=parseFloat(v);
  document.getElementById('exagVal').textContent=exag3d+'x';
  if(mlMap) mlMap.setTerrain({source:'terrain',exaggeration:exag3d});
}

function toggleHillshade(on){
  if(!mlMap) return;
  if(on){
    if(!mlMap.getSource('hs-src'))
      mlMap.addSource('hs-src',{type:'raster-dem',tiles:[TERRAIN_URL_3D],tileSize:256,encoding:'terrarium'});
    if(!mlMap.getLayer('hillshade'))
      mlMap.addLayer({id:'hillshade',type:'hillshade',source:'hs-src',
        paint:{'hillshade-exaggeration':0.6,'hillshade-shadow-color':'#222'}});
  } else {
    if(mlMap.getLayer('hillshade')) mlMap.removeLayer('hillshade');
  }
}

function changeBasemap(v){
  if(!mlMap) return;
  mlMap.setStyle(makeStyle3D(v));
  mlMap.once('style.load',()=>{
    mlMap.setTerrain({source:'terrain',exaggeration:exag3d});
    onMlLoad(); // 레이어 재추가
  });
}
</script>
"""

# ── 범례 (좌측 하단) ──────────────────────────────────────────────────────
LEGEND_HTML = """\
<div style="position:fixed;bottom:24px;left:14px;width:230px;
    background:rgba(255,255,255,0.96);border:2px solid #555;z-index:9999;
    padding:11px 13px;border-radius:8px;
    font-family:'Malgun Gothic',sans-serif;font-size:12.5px;line-height:1.75;
    box-shadow:2px 2px 6px rgba(0,0,0,0.28);">
  <b style="font-size:13.5px;">🗺 도상 선점 후보지</b><br>
  <span style="font-size:11px;color:#666;">총 %%TOTAL%%개 지점 | 2026.06.16</span>
  <hr style="margin:6px 0;">
  <span style="color:#FF8C00;">⭐</span>&thinsp;기존점 &nbsp;<b>%%C0%%</b>개<br>
  <span style="color:#2E8B57;">🌲</span>&thinsp;간선임도 &nbsp;<b>%%C1%%</b>개<br>
  <span style="color:#5B8C5A;">🌿</span>&thinsp;작업임도 &nbsp;<b>%%C2%%</b>개<br>
  <span style="color:#4488CC;">🔵</span>&thinsp;임도 외 &nbsp;<b>%%C3%%</b>개<br>
  <span style="color:#888;">⚫</span>&thinsp;기타 &nbsp;<b>%%C4%%</b>개<br>
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

    # ── Folium 지도 ──────────────────────────────────────────────────────────
    m = folium.Map(location=center, zoom_start=7, tiles=None)
    folium.TileLayer("https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}",
                     attr="Google Hybrid", name="🛰 구글 위성+도로", max_zoom=20).add_to(m)
    folium.TileLayer("https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}",
                     attr="Google Maps",   name="🗺 구글 지도",       max_zoom=20).add_to(m)
    folium.TileLayer("OpenStreetMap", name="🌍 OpenStreetMap").add_to(m)

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

    for p in points:
        mc, fi, label, _ = classify(p["name"], p["desc"])
        folium.Marker(
            location=[p["lat"], p["lon"]],
            popup=folium.Popup(popup_html(p), max_width=340),
            tooltip=p["name"],
            icon=folium.Icon(color=mc, icon=fi, prefix="fa"),
        ).add_to(fg_cat[label])

    m.save(str(OUT_PATH))
    html = OUT_PATH.read_text(encoding="utf-8")

    # ── 변수명 추출 ──────────────────────────────────────────────────────────
    ov_m   = re.search(r'overlays\s*:\s*\{([^}]+)\}', html, re.DOTALL)
    fg_vars= re.findall(r'(feature_group_[a-f0-9]+)', ov_m.group(1))
    fgv_power, fgv_rail, fgv_p1, fgv_existing = fg_vars[:4]

    div_m   = re.search(r'id="(map_[a-f0-9]+)"', html)
    map_div = div_m.group(1) if div_m else "map"
    var_m   = re.search(r'var\s+(map_[a-f0-9]+)\s*=\s*L\.map\(', html)
    map_var = var_m.group(1) if var_m else map_div
    print(f"  map_var={map_var}  map_div={map_div}")

    # ── 도상 선점 인라인 GeoJSON 생성 ────────────────────────────────────────
    features = []
    for p in points:
        kakao, naver = map_links(p["lat"], p["lon"], p["name"])
        features.append({
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [p["lon"], p["lat"]]},
            "properties": {
                "name":  p["name"],
                "desc":  p["desc"],
                "label": p["label"],
                "color": p["color"],
                "kakao": kakao,
                "naver": naver,
            }
        })
    dosung_gj_str = json.dumps({"type":"FeatureCollection","features":features},
                               ensure_ascii=False)

    # ── Leaflet JS 치환 ──────────────────────────────────────────────────────
    leaflet_js = (LEAFLET_JS
                  .replace("%%FG_POWER%%",    fgv_power)
                  .replace("%%FG_RAIL%%",     fgv_rail)
                  .replace("%%FG_P1%%",       fgv_p1)
                  .replace("%%FG_EXISTING%%", fgv_existing))

    # ── MapLibre JS 치환 ─────────────────────────────────────────────────────
    ml_block = (MAPLIBRE_BLOCK
                .replace("%%POINTS_GEOJSON%%", dosung_gj_str)
                .replace("%%TERRAIN_URL%%",    TERRAIN_URL)
                .replace("%%MAP_VAR%%",        map_var)
                .replace("%%MAP_DIV%%",        map_div))

    # ── 범례 치환 ────────────────────────────────────────────────────────────
    cc = {}
    for p in points:
        cc[p["label"]] = cc.get(p["label"], 0) + 1
    legend = (LEGEND_HTML
              .replace("%%TOTAL%%", str(len(points)))
              .replace("%%C0%%",    str(cc.get("기존점",  0)))
              .replace("%%C1%%",    str(cc.get("간선임도",0)))
              .replace("%%C2%%",    str(cc.get("작업임도",0)))
              .replace("%%C3%%",    str(cc.get("임도 외", 0)))
              .replace("%%C4%%",    str(cc.get("기타",    0))))

    html = html.replace("</body>", leaflet_js + ml_block + legend + "\n</body>")
    OUT_PATH.write_text(html, encoding="utf-8")
    print(f"저장 완료: {OUT_PATH}")
    return len(points)


if __name__ == "__main__":
    n = make_map()
    print(f"완료: {n}개 지점")
