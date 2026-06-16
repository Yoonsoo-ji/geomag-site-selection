"""
20260616_kang_site_v1.kml → docs/kang_site_map.html
도상 선점 후보지 + 기존 측정점 + P1 후보 + 제외구역 통합 지도
"""
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import folium

KML_PATH = Path(__file__).parent / "docs/output/20260616_kang_site_v1.kml"
OUT_PATH = Path(__file__).parent / "docs/kang_site_map.html"

# FeatureGroup 추가 순서 (overlays에서 이 순서대로 등장)
FG_ORDER = ["fg_power", "fg_rail", "fg_p1", "fg_existing",
            "fg_기존점", "fg_간선임도", "fg_작업임도", "fg_임도외", "fg_기타"]


def strip_cdata_html(raw: str) -> str:
    if not raw:
        return ""
    s = re.sub(r"<br\s*/?>", "\n", raw, flags=re.I)
    s = re.sub(r"</div>", "\n", s, flags=re.I)
    s = re.sub(r"<[^>]+>", "", s)
    s = s.replace("&nbsp;", " ").replace("&amp;", "&")
    return "\n".join(ln.strip() for ln in s.splitlines() if ln.strip())


def classify(name: str, desc: str):
    """(marker_color, fa_icon, label)"""
    combined = (name + " " + desc).lower()
    if "기존점" in combined:
        return "orange", "star", "기존점"
    if "간선임도" in combined:
        return "green", "tree", "간선임도"
    if "작업임도" in combined:
        return "darkgreen", "tree", "작업임도"
    if "임도 아님" in combined or "임도아님" in combined:
        return "blue", "info", "임도 외"
    return "gray", "circle", "기타"


def parse_kml():
    tree = ET.parse(KML_PATH)
    root = tree.getroot()
    points = []
    for pm in root.iter("{http://www.opengis.net/kml/2.2}Placemark"):
        name_el = pm.find("{http://www.opengis.net/kml/2.2}name")
        desc_el = pm.find("{http://www.opengis.net/kml/2.2}description")
        coord_el = pm.find(".//{http://www.opengis.net/kml/2.2}coordinates")
        if coord_el is None:
            continue
        coords = coord_el.text.strip().split(",")
        if len(coords) < 2:
            continue
        lon, lat = float(coords[0]), float(coords[1])
        name = name_el.text.strip() if name_el is not None and name_el.text else "이름 없음"
        desc_raw = desc_el.text.strip() if desc_el is not None and desc_el.text else ""
        desc = strip_cdata_html(desc_raw)
        points.append({"name": name, "desc": desc, "lat": lat, "lon": lon})
    return points


def make_popup_html(p):
    color_map = {"기존점": "#FF8C00", "간선임도": "#2E8B57",
                 "작업임도": "#5B8C5A", "임도 외": "#4488CC", "기타": "#888"}
    _, _, label = classify(p["name"], p["desc"])
    badge_color = color_map[label]
    desc_html = p["desc"].replace("\n", "<br>") if p["desc"] else "—"
    lat, lon = p["lat"], p["lon"]
    name_enc = p["name"].replace(" ", "%20")
    kakao_url = f"https://map.kakao.com/link/map/{name_enc},{lat:.6f},{lon:.6f}"
    naver_url = f"https://map.naver.com/v5/?c={lon:.6f},{lat:.6f},15,0,0,0,dh"

    return f"""<div style="font-family:'Malgun Gothic',sans-serif;font-size:13px;min-width:240px;max-width:320px;">
<b style="color:{badge_color};">{p['name']}</b>
&nbsp;<span style="background:{badge_color};color:white;padding:1px 6px;border-radius:3px;font-size:11px;">{label}</span>
<hr style="margin:5px 0;">
<b>위도:</b> {lat:.5f}°&thinsp;N &nbsp; <b>경도:</b> {lon:.5f}°&thinsp;E<br>
<hr style="margin:5px 0;border-color:#ddd;">
{desc_html}
<hr style="margin:5px 0;border-color:#ddd;">
<div style="display:flex;gap:6px;margin-top:2px;">
  <a href="{kakao_url}" target="_blank"
     style="flex:1;text-align:center;background:#FFCD00;color:#3A1D1D;text-decoration:none;
            font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">
    🗺 카카오맵
  </a>
  <a href="{naver_url}" target="_blank"
     style="flex:1;text-align:center;background:#03C75A;color:white;text-decoration:none;
            font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">
    🗺 네이버지도
  </a>
</div>
</div>"""


def make_map():
    points = parse_kml()
    print(f"파싱된 지점 수: {len(points)}")

    lats = [p["lat"] for p in points]
    lons = [p["lon"] for p in points]
    center = [sum(lats) / len(lats), sum(lons) / len(lons)]

    m = folium.Map(location=center, zoom_start=7, tiles=None)

    # 베이스 타일
    folium.TileLayer(
        "https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}",
        attr="Google Hybrid", name="🛰 구글 위성+도로", max_zoom=20,
    ).add_to(m)
    folium.TileLayer(
        "https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}",
        attr="Google Maps", name="🗺 구글 지도", max_zoom=20,
    ).add_to(m)
    folium.TileLayer("OpenStreetMap", name="🌍 OpenStreetMap").add_to(m)

    # FeatureGroup 순서 엄격히 유지 (FG_ORDER와 동일)
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

    # 지도에 추가 (이 순서 = overlays 순서)
    for fg in [fg_power, fg_rail, fg_p1, fg_existing] + list(fg_cat.values()):
        m.add_child(fg)

    folium.LayerControl(collapsed=False).add_to(m)

    # 도상 선점 마커
    for p in points:
        mc, fi, label = classify(p["name"], p["desc"])
        folium.Marker(
            location=[p["lat"], p["lon"]],
            popup=folium.Popup(make_popup_html(p), max_width=340),
            tooltip=p["name"],
            icon=folium.Icon(color=mc, icon=fi, prefix="fa"),
        ).add_to(fg_cat[label])

    # HTML 저장
    m.save(str(OUT_PATH))
    html_src = OUT_PATH.read_text(encoding="utf-8")

    # overlays 섹션에서 변수명 순서대로 추출
    overlays_m = re.search(r'overlays\s*:\s*\{([^}]+)\}', html_src, re.DOTALL)
    if not overlays_m:
        print("ERROR: overlays 섹션을 찾지 못했습니다")
        return len(points)

    fg_vars = re.findall(r'(feature_group_[a-f0-9]+)', overlays_m.group(1))
    print("overlays 순서:", fg_vars)

    if len(fg_vars) < 4:
        print("ERROR: FeatureGroup 변수 수가 부족합니다")
        return len(points)

    fgv_power    = fg_vars[0]  # fg_power
    fgv_rail     = fg_vars[1]  # fg_rail
    fgv_p1       = fg_vars[2]  # fg_p1
    fgv_existing = fg_vars[3]  # fg_existing

    # JS 주입 블록
    js_inject = f"""
<script>
(function(){{
  // ── 고압철탑·송전탑 제외구역 ──────────────────────────────────────
  fetch('data/zone_power.geojson')
    .then(function(r){{return r.json();}})
    .then(function(d){{
      L.geoJSON(d,{{
        style:function(){{return{{fillColor:'#FF3300',color:'#CC0000',weight:1,fillOpacity:0.28}};}},
        onEachFeature:function(f,l){{l.bindTooltip('⚡ 고압철탑·송전탑 제외 (1.0 km)');}}
      }}).addTo({fgv_power});
    }}).catch(function(e){{console.warn('zone_power:',e);}});

  // ── 철도 제외구역 ─────────────────────────────────────────────────
  fetch('data/zone_railway.geojson')
    .then(function(r){{return r.json();}})
    .then(function(d){{
      L.geoJSON(d,{{
        style:function(){{return{{fillColor:'#FF7700',color:'#CC4400',weight:1,fillOpacity:0.28}};}},
        onEachFeature:function(f,l){{l.bindTooltip('🚆 철도 제외 (5.0 km)');}}
      }}).addTo({fgv_rail});
    }}).catch(function(e){{console.warn('zone_railway:',e);}});

  // ── P1 우선순위 1등급 후보지 ─────────────────────────────────────
  fetch('data/candidates_p1.geojson')
    .then(function(r){{return r.json();}})
    .then(function(d){{
      L.geoJSON(d,{{
        pointToLayer:function(f,ll){{
          return L.circleMarker(ll,{{radius:7,color:'#CC0000',fillColor:'#FF2200',fillOpacity:0.85,weight:1.5}});
        }},
        onEachFeature:function(f,l){{
          var p=f.properties;
          var vn=function(x){{return(x===null||x===undefined||isNaN(+x))?'-':(+x).toFixed(1);}};
          var lat=(p.lat!=null)?p.lat:f.geometry.coordinates[1];
          var lon=(p.lon!=null)?p.lon:f.geometry.coordinates[0];
          var kakao='https://map.kakao.com/link/map/P1%ED%9B%84%EB%B3%B4'+p.idx+','+lat.toFixed(6)+','+lon.toFixed(6);
          var naver='https://map.naver.com/v5/?c='+lon.toFixed(6)+','+lat.toFixed(6)+',15,0,0,0,dh';
          var html='<div style="font-family:Malgun Gothic,sans-serif;font-size:12.5px;min-width:240px;">'
            +'<b style="color:#CC0000;">후보지 #'+p.idx+'</b>'
            +'&nbsp;<span style="background:#FF2200;color:white;padding:1px 5px;border-radius:3px;font-size:11px;">P1 최우선</span><br>'
            +'<hr style="margin:4px 0;">'
            +'<b>위도:</b> '+lat.toFixed(5)+'° N　<b>경도:</b> '+lon.toFixed(5)+'° E<br>'
            +'<b>입지점수:</b> '+vn(p.score)+'/100<br>'
            +'<hr style="margin:4px 0;border-color:#ddd;">'
            +'<div style="display:flex;gap:6px;margin-top:2px;">'
            +'<a href="'+kakao+'" target="_blank" style="flex:1;text-align:center;background:#FFCD00;color:#3A1D1D;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 카카오맵</a>'
            +'<a href="'+naver+'" target="_blank" style="flex:1;text-align:center;background:#03C75A;color:white;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 네이버지도</a>'
            +'</div></div>';
          l.bindPopup(html,{{maxWidth:300}});
          l.bindTooltip('P1 #'+p.idx+' ('+lat.toFixed(4)+'°N)');
        }}
      }}).addTo({fgv_p1});
    }}).catch(function(e){{console.warn('candidates_p1:',e);}});

  // ── 기존 측정점 ───────────────────────────────────────────────────
  fetch('data/existing_sites.geojson')
    .then(function(r){{return r.json();}})
    .then(function(d){{
      L.geoJSON(d,{{
        pointToLayer:function(f,ll){{
          var icon=L.divIcon({{html:'<div style="font-size:22px;line-height:1;text-shadow:1px 1px 2px rgba(0,0,0,0.5);">⭐</div>',className:'',iconAnchor:[11,11]}});
          return L.marker(ll,{{icon:icon}});
        }},
        onEachFeature:function(f,l){{
          var p=f.properties;
          var vf=function(x,fmt){{return(x===null||x===undefined)?'-':fmt(x);}};
          var lat=p.lat, lon=p.lon;
          var kakao='https://map.kakao.com/link/map/'+encodeURIComponent(p.name)+','+lat.toFixed(6)+','+lon.toFixed(6);
          var naver='https://map.naver.com/v5/?c='+lon.toFixed(6)+','+lat.toFixed(6)+',15,0,0,0,dh';
          var html='<div style="font-family:Malgun Gothic,sans-serif;font-size:12.5px;min-width:240px;">'
            +'<b style="color:#8B4513;">⭐ '+p.name+'</b><br>'
            +'<hr style="margin:4px 0;">'
            +'<b>위도:</b> '+lat.toFixed(5)+'° N　<b>경도:</b> '+lon.toFixed(5)+'° E<br>'
            +'<b>표고:</b> '+vf(p.elev,function(x){{return x.toFixed(1)+' m'}})+' &nbsp; '
            +'<b>최신관측:</b> '+p.obs_year+'년<br>'
            +'<b>총자력:</b> '+vf(p.total,function(x){{return x.toLocaleString('ko-KR',{{minimumFractionDigits:1,maximumFractionDigits:1}})+' nT'}})+' &nbsp; '
            +'<b>편각:</b> '+vf(p.decl,function(x){{return x.toFixed(2)+'°'}})+' &nbsp; '
            +'<b>복각:</b> '+vf(p.incl,function(x){{return x.toFixed(2)+'°'}})+' <br>'
            +'<hr style="margin:4px 0;border-color:#ddd;">'
            +'<div style="display:flex;gap:6px;margin-top:2px;">'
            +'<a href="'+kakao+'" target="_blank" style="flex:1;text-align:center;background:#FFCD00;color:#3A1D1D;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 카카오맵</a>'
            +'<a href="'+naver+'" target="_blank" style="flex:1;text-align:center;background:#03C75A;color:white;text-decoration:none;font-weight:bold;font-size:12px;padding:5px 0;border-radius:4px;">🗺 네이버지도</a>'
            +'</div></div>';
          l.bindPopup(html,{{maxWidth:300}});
          l.bindTooltip('⭐ '+p.name+' ('+p.obs_year+'년)');
        }}
      }}).addTo({fgv_existing});
    }}).catch(function(e){{console.warn('existing_sites:',e);}});
}})();
</script>"""

    # 우측 상단 범례 패널
    total = len(points)
    cat_counts = {}
    for p in points:
        _, _, label = classify(p["name"], p["desc"])
        cat_counts[label] = cat_counts.get(label, 0) + 1

    legend_html = f"""
<div style="
    position:fixed; top:10px; right:10px; width:250px;
    background:rgba(255,255,255,0.96); border:2px solid #555;
    z-index:9999; padding:12px 14px; border-radius:8px;
    font-family:'Malgun Gothic',sans-serif; font-size:13px;
    box-shadow:2px 2px 6px rgba(0,0,0,0.3);">
  <b style="font-size:14px;">🗺 도상 선점 후보지</b><br>
  <span style="font-size:11px;color:#555;">총 {total}개 지점 | 2026.06.16</span>
  <hr style="margin:7px 0;">
  <span style="color:#FF8C00;">⭐</span> 기존점 &nbsp; <b>{cat_counts.get('기존점',0)}</b>개<br>
  <span style="color:#2E8B57;">🌲</span> 간선임도 &nbsp; <b>{cat_counts.get('간선임도',0)}</b>개<br>
  <span style="color:#5B8C5A;">🌿</span> 작업임도 &nbsp; <b>{cat_counts.get('작업임도',0)}</b>개<br>
  <span style="color:#4488CC;">🔵</span> 임도 외 &nbsp; <b>{cat_counts.get('임도 외',0)}</b>개<br>
  <span style="color:#888;">⚫</span> 기타 &nbsp; <b>{cat_counts.get('기타',0)}</b>개<br>
  <hr style="margin:7px 0;border-color:#ccc;">
  <a href="index.html" style="color:#1a73e8;text-decoration:none;font-weight:bold;">
    ← 입지 선정 분석 지도</a>
</div>"""

    # </body> 앞에 JS + 범례 주입
    html_src = html_src.replace("</body>", js_inject + "\n" + legend_html + "\n</body>")
    OUT_PATH.write_text(html_src, encoding="utf-8")
    print(f"저장 완료: {OUT_PATH}")
    return len(points)


if __name__ == "__main__":
    n = make_map()
    print(f"완료: {n}개 지점")
