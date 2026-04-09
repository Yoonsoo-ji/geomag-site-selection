"""
docs/ 폴더 빌드 스크립트 — GitHub Pages 배포용
- 대형 GeoJSON 파일은 shapely로 단순화 (좌표 정밀도 축소 포함)
- HTML을 index.html로 복사
- output/data/ → docs/data/
"""
import json
import shutil
from pathlib import Path

import geopandas as gpd
import shapely

SRC_DATA = Path("output/data")
SRC_HTML = Path("output/geomag_site_selection.html")
DST_DIR  = Path("docs")
DST_DATA = DST_DIR / "data"

# 단순화 허용 오차 (도, WGS84 기준)
# 0.001° ≈ 100 m,  0.0002° ≈ 20 m
SIMPLIFY = {
    "zone_urban_resid.geojson": 0.001,   # 233 MB → 목표 3 MB 이하
    "zone_urban_dense.geojson": 0.001,   # 46 MB  → 목표 2 MB 이하
    "zone_power.geojson":       0.0005,  # 17 MB  → 목표 2 MB 이하
    "zone_railway.geojson":     0.0002,  # 1.8 MB → 이미 작지만 추가 압축
}

DST_DIR.mkdir(exist_ok=True)
DST_DATA.mkdir(exist_ok=True)

# ── HTML 복사 ──────────────────────────────────────────────────────────
shutil.copy(SRC_HTML, DST_DIR / "index.html")
print(f"[OK] index.html 복사")

# ── GeoJSON 처리 ──────────────────────────────────────────────────────
for src in sorted(SRC_DATA.glob("*.geojson")):
    dst = DST_DATA / src.name
    tol = SIMPLIFY.get(src.name)

    if tol is None:
        # 소형 파일 → 그대로 복사
        shutil.copy(src, dst)
        mb = src.stat().st_size / 1e6
        print(f"[CP] {src.name:<35s}  {mb:6.1f} MB")
        continue

    # 대형 파일 → GeoDataFrame으로 읽어 단순화 후 저장
    src_mb = src.stat().st_size / 1e6
    print(f"[..] {src.name:<35s}  {src_mb:6.1f} MB  단순화 중 (tol={tol})...", end="", flush=True)

    try:
        gdf = gpd.read_file(src)
        gdf["geometry"] = gdf.geometry.simplify(tolerance=tol, preserve_topology=True)
        gdf = gdf[~gdf.geometry.is_empty]
        gdf["geometry"] = gdf.geometry.apply(
            lambda g: shapely.set_precision(g, grid_size=10**(-5))
        )
        gdf.to_file(dst, driver="GeoJSON")
    except Exception:
        # geopandas 실패 시 json + shapely 직접 처리 (초대형 단일 피처용)
        import json
        from shapely.geometry import shape, mapping

        with open(src, encoding="utf-8") as f:
            fc = json.load(f)

        out_features = []
        for feat in fc.get("features", []):
            geom = shape(feat["geometry"])
            geom = geom.simplify(tol, preserve_topology=True)
            geom = shapely.set_precision(geom, grid_size=1e-5)
            if not geom.is_empty:
                out_features.append({
                    "type": "Feature",
                    "properties": feat.get("properties", {}),
                    "geometry": mapping(geom),
                })

        with open(dst, "w", encoding="utf-8") as f:
            json.dump({"type": "FeatureCollection", "features": out_features}, f)
    dst_mb = dst.stat().st_size / 1e6
    print(f" → {dst_mb:5.1f} MB  ({dst_mb/src_mb*100:.0f}%)")

print("\n배포 파일 준비 완료: docs/")
print(f"  HTML : docs/index.html")
print(f"  Data : docs/data/ ({len(list(DST_DATA.glob('*.geojson')))} 파일)")
print()

# 파일별 최종 크기 요약
total = 0
for f in sorted(DST_DATA.glob("*.geojson")):
    mb = f.stat().st_size / 1e6
    total += mb
    flag = " ⚠ 100MB 초과!" if mb > 100 else (" ⚠ 50MB 초과" if mb > 50 else "")
    print(f"  {f.name:<35s}  {mb:6.1f} MB{flag}")
print(f"  {'합계':<35s}  {total:6.1f} MB")
