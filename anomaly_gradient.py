# -*- coding: utf-8 -*-
"""
항공자력 공간구배 — 입지 대표성과 관측망 밀도 설계의 단일 출처
================================================================

KIGAM 자력이상도(1.5분 ≈ 2.8 km)에서 **두 개의 서로 다른 축**을 산출한다.
둘을 섞으면 결론이 뒤집히므로 코드에서 분리해 둔다.

  ① 점(site) 축 — 국소 공간구배 |∇ΔT| (nT/km)
     한 측점이 «자기 주변을 대표하는가». 구배가 클수록 그 지점의 값은
     주변을 대표하지 못하고, 2.8 km 격자로 만든 Crustal 보정도 그 지점에서
     빗나간다. **저구배가 좋다.** (IAGA 반복관측점 입지 지침과 같은 방향)

  ② 권역(region) 축 — 자기복잡도 σ(R) (nT) → 권장 측점 밀도
     한 권역이 «몇 점을 필요로 하는가». 지각장이 복잡할수록 한 점의
     대표성이 떨어지므로 더 많은 점으로 평균해야 한다. **고복잡도일수록
     조밀해야 한다.**

⚠️ 두 축은 방향이 반대로 보이지만 모순이 아니다 — **복잡한 권역에 점을 더
   많이 두되, 각 점은 그 권역 안의 조용한 자리에 놓는다**는 하나의 설계다.
   종전 `geomag_site_selection.anomaly_gradient_score()` 는 「변동폭이 클수록
   모델 기여가 크다」며 중간 변동폭(250 nT)에 최고점을 줬는데, 그것은 ②의
   논리를 ①의 자리에 넣은 것이다. 이 모듈이 그 둘을 갈라 놓는다.

근거
----
* 일본 국토지리원 — 2010.0 자기도 검증에서 야쓰가타케·홋카이도 동부의 큰
  오차를 **「자기구배의 크기에 비해 1·2등 자기점의 배점밀도가 낮은 것」**으로
  판단했다. 지상망(평균 약 20 km 간격)이 표현하지 못하는 단파장 자기이상을
  전국 항공자력측량으로 보완한다고 밝히고 있다.
* 국토지리정보원 2010 「지구물리측량 연구」 — 신규 지자기점 설치 시 전국
  지자기이상도와 현지 자력탐사로 큰 이상체를 사전 확인할 것을 요구했고,
  제주는 **「지자기이상이 복잡하게 나타나 현재까지의 3개 지점 측정으로는
  지역을 대표할 만한 지자기점을 매설하기 매우 어렵다」**고 적었다.
* 영국 BGS — 반복관측점은 core/regional, 항공자력은 crustal 로 **역할을
  분리**한다. 우리 LMM 의 Regional / Crustal 층 분리와 같은 구조다.

⚠️ 위 어느 사례도 **「구배가 X 이상이면 밀도를 Y 배」 같은 정량 규칙까지
   가지는 않았다.** 이 모듈의 Neyman 배분(아래)은 그 빈자리를 메우려는
   이 연구의 제안이지 선례의 인용이 아니다. 산출물에 선례로 쓰지 말 것.

⚠️ **우리 16측점으로는 ①이 검증되지 않는다** (`validate_against_lmm()`, 2026-08-27)

    |잔차| vs 국소구배 |∇ΔT|      n=14   rD −0.22 (p 0.45) · rI −0.40 (p 0.16)
                                        · rF −0.18 (p 0.54)
    |잔차| vs 권역복잡도 σ(25km)  n=15   rD **+0.47 (p 0.080)** · rI −0.06 · rF −0.03

  즉 **②(권역 축)는 약하게나마 예측 방향으로 나오는데, ①(점 축)은 나오지
  않는다.** 그렇다고 ①이 틀렸다고 읽어서는 안 된다 — 이 표본의 편각 잔차는
  **방문마다 달라지는 마크 방위각 오차가 지배**하고(재방문 잔여 RMS 33.7′,
  CLAUDE.md §13), 14점은 상관을 가릴 검정력이 없다.

  ⇒ ①의 근거는 **물리 추론과 IAGA 관행**이지 이 연구의 실측이 아니다.
     그렇게 표기할 것. 방위각 계통오차를 걷어낸 뒤 다시 검정해야 판정된다.
     되돌리려면 `geomag_site_selection.USE_GRADIENT_SCORE = False` 하나면 된다.

권장 밀도 — Neyman 배분
------------------------
층화표본에서 층별 평균을 같은 정밀도로 추정하려면 표본을 **층 크기 × 층
표준편차**에 비례해 배분한다(Neyman allocation). 권역을 층, 지각장 산포
σ(R) 을 층 표준편차로 놓으면

    n_h ∝ A_h · σ_h          →      밀도 ρ(x) ∝ σ(x)
    권장 측점 간격 L(x) = 1 / sqrt(ρ(x))

가 된다. 「구배가 크면 조밀하게」라는 정성 논리를 예산(총 측점 수) 제약
아래에서 재현 가능한 수치로 바꾼 것이다.

    python anomaly_gradient.py        # 자기점검 — 분포·권역 요약 출력
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent
KIGAM_DAT = ROOT / "data" / "mag_1982-2018_1.5min_ed.dat"

EARTH_R_KM = 6371.0
GRID_STEP_DEG = 0.025          # KIGAM 1.5분

# ── ① 점 축 — 대표성 점수 곡선 ────────────────────────────────────────
# 로그 구배축의 반쪽 정규곡선. G0 이하는 평탄으로 보고 만점.
#   5 nT/km→10.0 · 10→8.9 · 15→7.4 · 25→5.2 · 50→2.6 · 100→1.0 · 180→0.4
# 후보 178점 실측 분포: P10 3.4 · 중앙 14.3 · P75 24.6 · P90 45.1 · 최대 180.6
GRAD_FLAT_NT_KM = 5.0          # 이하 만점 — 자기적으로 평탄
GRAD_SIGMA_LOG = 1.4           # 로그축 감쇠 폭
GRAD_NODATA_SCORE = 6.0        # 자료 미취득 — 검증되지 않았으므로 중앙값보다 보수적

# ── 구배 등급 구간 (지도 레이어 · nT/km) ──────────────────────────────
#   국토 격자 실측: P10 4.2 · P25 7.8 · 중앙 16.2 · P75 35 · P90 65.9 · P99 142
GRAD_TIERS = [
    ("flat",     0.0,    5.0, "매우 평탄"),
    ("gentle",   5.0,   15.0, "양호"),
    ("moderate", 15.0,  40.0, "보통"),
    ("steep",    40.0, 100.0, "급변 — 현장 검증"),
    ("extreme", 100.0,  np.inf, "극심 — 현장 정밀 자력측량"),
]

# ── ② 권역 축 ─────────────────────────────────────────────────────────
SIGMA_RADIUS_KM = 25.0         # 복잡도 평가 반경 — 측점 간격 규모
SIGMA_MIN_PTS = 8              # 반경 내 최소 유효점

# 배분 대상 총 측점 수. 현 LMM 투입 16점 + 후반기 Track C 신규 50점.
# ⚠ 이것은 설계 가정이다 — 사업 계획이 바뀌면 이 상수 하나만 고친다.
N_TARGET_SITES = 66

# 권장 측점 간격 등급 (km)
SPACING_TIERS = [
    ("dense",  0.0,  35.0, "조밀 배치 필요"),
    ("mid",   35.0,  50.0, "중간"),
    ("coarse", 50.0, 65.0, "보통"),
    ("sparse", 65.0, np.inf, "성김 허용"),
]

# 밀도 충족도 = 현 측점 간격 / 권장 간격 (1 미만이면 충족)
DEFICIT_TIERS = [
    ("ok",       0.0, 1.0, "충족"),
    ("mild",     1.0, 1.5, "다소 부족"),
    ("short",    1.5, 2.0, "부족"),
    ("critical", 2.0, np.inf, "크게 부족"),
]


# ======================================================================
# 격자 적재
# ======================================================================

def load_points(path: Path | None = None):
    """KIGAM 원자료를 (lon, lat, anomaly) 배열로. 격자화하지 않은 원점."""
    path = path or KIGAM_DAT
    a = np.loadtxt(path, skiprows=9)
    return a[:, 0], a[:, 1], a[:, 2]


def load_grid(path: Path | None = None):
    """
    KIGAM 자력이상도를 **균일 간격** 규칙격자로 반환.

    ⚠️ 원자료 위도축에는 제주해협(33.55°~34.1°)에 0.55° 공백이 있어
    `np.unique()` 결과를 그대로 축으로 쓰면 (값−원점)/간격 색인이 공백 위쪽
    전 구간에서 21행 어긋난다. 결측행을 NaN 으로 채운 완전 균일축을 만든다.
    (`lmm_build.load_kigam_grid` 와 같은 처리 — 사본이 아니라 같은 규약이다)
    """
    lon, lat, anom = load_points(path)
    step = GRID_STEP_DEG
    lons = np.round(np.arange(lon.min(), lon.max() + step / 2, step), 4)
    lats = np.round(np.arange(lat.min(), lat.max() + step / 2, step), 4)
    grid = np.full((lats.size, lons.size), np.nan)
    grid[np.rint((lat - lats[0]) / step).astype(int),
         np.rint((lon - lons[0]) / step).astype(int)] = anom
    return lons, lats, grid


def gradient_field(lons, lats, grid):
    """
    중앙차분으로 수평 공간구배를 계산해 (|∇ΔT|, ∂/∂동, ∂/∂북) 를 nT/km 로 반환.

    경도 방향 격자 폭은 위도에 따라 달라지므로(cos φ) 행마다 다른 간격을 쓴다.
    NaN(미측선) 이웃이 하나라도 있으면 그 셀의 구배도 NaN 이다 — 자료 공백을
    0 으로 메우면 인접 셀에 가짜 급변이 생기기 때문이다.
    """
    step = GRID_STEP_DEG
    dy_km = step * np.pi / 180.0 * EARTH_R_KM
    dx_km = step * np.pi / 180.0 * EARTH_R_KM * np.cos(np.radians(lats))

    gy = np.full_like(grid, np.nan)
    gx = np.full_like(grid, np.nan)
    gy[1:-1, :] = (grid[2:, :] - grid[:-2, :]) / (2.0 * dy_km)
    gx[:, 1:-1] = (grid[:, 2:] - grid[:, :-2]) / (2.0 * dx_km[:, None])
    return np.hypot(gx, gy), gx, gy


# ======================================================================
# 조회기
# ======================================================================

class AnomalyGradient:
    """항공자력 격자에서 구배·복잡도를 조회한다."""

    def __init__(self, path: Path | None = None):
        self.lons, self.lats, self.grid = load_grid(path)
        self.G, self.gx, self.gy = gradient_field(self.lons, self.lats, self.grid)
        self.step = GRID_STEP_DEG
        # 반경 질의용 원점 배열 (격자화 전 — NaN 이 섞이지 않는다)
        self.plon, self.plat, self.panom = load_points(path)

    # ── ① 국소 구배 ──────────────────────────────────────────────
    def grad(self, lat, lon):
        """|∇ΔT| (nT/km). 격자 밖·자료 공백은 NaN, 인접 3×3 에 값이 있으면 평균."""
        lat = np.atleast_1d(np.asarray(lat, float))
        lon = np.atleast_1d(np.asarray(lon, float))
        out = np.full(lat.shape, np.nan)

        j = np.rint((lat - self.lats[0]) / self.step).astype(int)
        i = np.rint((lon - self.lons[0]) / self.step).astype(int)
        inside = ((j >= 0) & (j < self.lats.size) &
                  (i >= 0) & (i < self.lons.size))
        out[inside] = self.G[j[inside], i[inside]]

        # 자료 공백 경계 — 3×3 이웃 평균으로 한 겹 메운다
        for k in np.where(inside & ~np.isfinite(out))[0]:
            w = self.G[max(0, j[k] - 1):j[k] + 2, max(0, i[k] - 1):i[k] + 2]
            if np.isfinite(w).any():
                out[k] = np.nanmean(w)
        return out

    # ── ② 권역 복잡도 ────────────────────────────────────────────
    def sigma(self, lat, lon, r_km: float = SIGMA_RADIUS_KM, detrend: bool = True):
        """
        반경 r_km 내 자기이상의 표준편차 (nT) — 권역 자기복잡도.

        `detrend=True` 면 국소 1차 평면을 뺀다. 광역 추세는 Regional 층이
        담당하므로 밀도 설계의 근거에서 제외하는 것이 옳다(실측상 차이는
        3% 로 작지만 정의가 명확해진다).
        """
        lat = np.atleast_1d(np.asarray(lat, float))
        lon = np.atleast_1d(np.asarray(lon, float))
        out = np.full(lat.shape, np.nan)

        for k, (la, lo) in enumerate(zip(lat, lon)):
            dlat = r_km / 110.574
            dlon = r_km / (111.320 * np.cos(np.radians(la)))
            msk = ((np.abs(self.plat - la) <= dlat) &
                   (np.abs(self.plon - lo) <= dlon))
            v = self.panom[msk]
            if v.size < SIGMA_MIN_PTS:
                continue
            if not detrend:
                out[k] = float(np.std(v))
                continue
            x = (self.plon[msk] - lo) * 111.320 * np.cos(np.radians(la))
            y = (self.plat[msk] - la) * 110.574
            A = np.column_stack([np.ones_like(x), x, y])
            try:
                coef, *_ = np.linalg.lstsq(A, v, rcond=None)
                out[k] = float(np.std(v - A @ coef))
            except np.linalg.LinAlgError:
                out[k] = float(np.std(v))
        return out


# ======================================================================
# 점수·등급
# ======================================================================

def representativeness_score(g_nt_km, max_score: float = 10.0):
    """
    국소 구배 → 자기환경 대표성 점수 (0~max_score). 저구배가 높은 점수.

    ⚠️ 이 함수는 종전 `anomaly_gradient_score()`(변동폭 250 nT 에서 최고점)와
    **방향이 반대**다. 종전 곡선은 「변동이 큰 곳이 모델에 기여한다」는 권역
    논리를 점 점수에 넣은 것이었다. 권역 논리는 `neyman_density()` 로 옮겼다.
    """
    g = np.asarray(g_nt_km, float)
    out = np.full(g.shape, np.nan)
    ok = np.isfinite(g) & (g >= 0)
    x = np.log(np.maximum(np.where(ok, g, 1.0), 1e-6) / GRAD_FLAT_NT_KM)
    x = np.maximum(x, 0.0)                       # 평탄 구간은 감점 없음
    out[ok] = max_score * np.exp(-0.5 * (x[ok] / GRAD_SIGMA_LOG) ** 2)
    return out


def _tier_of(value, tiers):
    if not np.isfinite(value):
        return None
    for key, lo, hi, _label in tiers:
        if lo <= value < hi:
            return key
    return tiers[-1][0]


def grad_tier(value):
    """|∇ΔT| → 등급 키."""
    return _tier_of(value, GRAD_TIERS)


def spacing_tier(value):
    """권장 측점 간격(km) → 등급 키."""
    return _tier_of(value, SPACING_TIERS)


def deficit_tier(value):
    """밀도 충족도(현 간격 / 권장 간격) → 등급 키."""
    return _tier_of(value, DEFICIT_TIERS)


TIER_LABEL = {k: lab for k, _, _, lab in GRAD_TIERS}
TIER_LABEL.update({k: lab for k, _, _, lab in SPACING_TIERS})
TIER_LABEL.update({k: lab for k, _, _, lab in DEFICIT_TIERS})


# ======================================================================
# Neyman 배분 — 권장 측점 밀도·간격
# ======================================================================

def neyman_density(sigma_nt, area_km2, n_target: int = N_TARGET_SITES):
    """
    권역 자기복잡도 σ 에서 권장 측점 밀도·간격을 산출.

        n_h ∝ A_h · σ_h        (Neyman allocation)
        ρ(x) = n_target · σ(x) / Σ(A·σ)          [측점 / km²]
        L(x) = 1 / sqrt(ρ(x))                     [km]

    Parameters
    ----------
    sigma_nt  : 격자별 σ(R) (nT). NaN 은 유한값 중앙값으로 대체된다.
    area_km2  : 격자별 면적 (km²)
    n_target  : 배분할 총 측점 수

    Returns
    -------
    (density, spacing, sigma_filled) — 측점/km², km, 대체 후 σ
    """
    sigma_nt = np.asarray(sigma_nt, float)
    area_km2 = np.asarray(area_km2, float)
    finite = np.isfinite(sigma_nt)
    if not finite.any():
        raise ValueError("σ 가 전부 결측 — 항공자력 자료 범위를 확인할 것")
    filled = np.where(finite, sigma_nt, np.median(sigma_nt[finite]))

    w = area_km2 * filled
    density = n_target * filled / w.sum()        # = n·A·σ/Σ(Aσ) / A
    spacing = 1.0 / np.sqrt(density)
    return density, spacing, filled


def uniform_spacing(area_km2_total, n_target: int = N_TARGET_SITES):
    """비교 기준 — 같은 측점 수를 균등 배치했을 때의 간격 (km)."""
    return float(np.sqrt(area_km2_total / n_target))


def cell_area_km2(lat, dlat_deg, dlon_deg):
    """경위도 셀의 근사 면적 (km²)."""
    lat = np.asarray(lat, float)
    return (dlat_deg * 110.574) * (dlon_deg * 111.320 * np.cos(np.radians(lat)))


# ======================================================================
# 검정 — 우리 측점 잔차로 두 축을 확인한다
# ======================================================================

def validate_against_lmm(ag: "AnomalyGradient | None" = None, verbose: bool = True):
    """
    LMM 측점의 IGRF 잔차가 ① 국소구배 ② 권역복잡도 와 어떤 관계인지 잰다.

    ①이 옳다면 |잔차| 가 |∇ΔT| 와 **양의 상관**이어야 한다(구배가 크면 그
    지점의 값이 주변을 대표하지 못하므로). ②가 옳다면 σ 와 양의 상관이다.

    ⚠️ 이 표본의 편각 잔차는 방위 기준 계통오차가 지배하므로 검정력이 낮다.
    귀무가설을 기각하지 못한 것을 「관계가 없다」로 읽지 말 것.

    Returns
    -------
    dict — 성분별 {n, r, p} (grad / sigma 각각)
    """
    import json as _json

    from scipy import stats

    ag = ag or AnomalyGradient()
    path = ROOT / "docs" / "data" / "lmm_diagnosis.json"
    if not path.exists():
        if verbose:
            print(f"[건너뜀] {path.name} 없음 — lmm_diagnose.py 먼저 실행")
        return {}
    tbl = _json.load(open(path, encoding="utf-8")).get("site_table", [])
    if not tbl:
        return {}

    lat = np.array([r["lat"] for r in tbl], float)
    lon = np.array([r["lon"] for r in tbl], float)
    g = ag.grad(lat, lon)
    s = ag.sigma(lat, lon)

    out = {}
    if verbose:
        print(f"\n■ 측점 잔차 검정 (n={len(tbl)})  — 양의 상관이면 가설 지지")
        print(f"   {'':10s}{'|∇ΔT| r':>10}{'p':>8}   {'σ(25km) r':>11}{'p':>8}")
    for comp in ("rD", "rI", "rF"):
        y = np.abs(np.array([r.get(comp, np.nan) for r in tbl], float))
        row = {}
        for key, x in (("grad", g), ("sigma", s)):
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 4:
                row[key] = {"n": int(ok.sum()), "r": None, "p": None}
                continue
            r, p = stats.pearsonr(x[ok], y[ok])
            row[key] = {"n": int(ok.sum()), "r": float(r), "p": float(p)}
        out[comp] = row
        if verbose:
            a, b = row["grad"], row["sigma"]
            print(f"   |{comp}|{'':5s}{a['r']:>+10.3f}{a['p']:>8.3f}   "
                  f"{b['r']:>+11.3f}{b['p']:>8.3f}   (n {a['n']}/{b['n']})")
    if verbose:
        print("   ⚠ 점 축(①)은 이 표본에서 지지되지 않는다 — 편각 잔차를 방위")
        print("     기준 계통오차가 덮고 있고 14점은 검정력이 없다. ①의 근거는")
        print("     물리 추론과 IAGA 관행이며 이 연구의 실측이 아니다.")
    return out


# ======================================================================
# 자기점검
# ======================================================================

def _self_check():
    import json

    ag = AnomalyGradient()
    G = ag.G
    ok = np.isfinite(G)
    print(f"KIGAM 격자 {ag.grid.shape}  유효 {np.isfinite(ag.grid).sum():,} "
          f"({np.isfinite(ag.grid).mean()*100:.0f}%)")
    print(f"구배 유효셀 {ok.sum():,}")
    q = np.percentile(G[ok], [10, 25, 50, 75, 90, 99])
    print("|∇ΔT| nT/km  P10/25/50/75/90/99:", np.round(q, 1),
          f" 최대 {G[ok].max():.0f}")

    print("\n대표성 점수 곡선")
    for g in (2, 5, 10, 15, 25, 40, 60, 100, 180):
        print(f"   {g:>4} nT/km → {float(representativeness_score(g)):5.2f} / 10"
              f"   [{TIER_LABEL[grad_tier(g)]}]")

    # 국토 격자에서 Neyman 배분
    data = ROOT / "docs" / "data"
    gj = json.load(open(data / "korea_boundary.geojson", encoding="utf-8"))
    from matplotlib.path import Path as MPath
    polys = []
    for f in gj["features"]:
        g = f["geometry"]
        for poly in ([g["coordinates"]] if g["type"] == "Polygon"
                     else g["coordinates"]):
            polys.append(MPath(np.asarray(poly[0])))
    la = np.arange(33.1, 38.7, 0.1)
    lo = np.arange(125.0, 129.7, 0.1)
    LO, LA = np.meshgrid(lo, la)
    pts = np.column_stack([LO.ravel(), LA.ravel()])
    m = np.zeros(len(pts), bool)
    for pp in polys:
        m |= pp.contains_points(pts)
    m = m.reshape(LA.shape)
    gn, gl = LA[m], LO[m]

    sg = ag.sigma(gn, gl)
    A = cell_area_km2(gn, 0.1, 0.1)
    dens, spac, _ = neyman_density(sg, A)
    print(f"\n국토 격자 {gn.size}점 · σ(R={SIGMA_RADIUS_KM:.0f}km) "
          f"결측 {int(np.isnan(sg).sum())}")
    print("σ nT  P10/50/90:", np.round(np.nanpercentile(sg, [10, 50, 90]), 0))
    print(f"권장 간격 km  P10/50/90: {np.round(np.percentile(spac,[10,50,90]),0)}"
          f"   (균등 배치 기준 {uniform_spacing(A.sum()):.0f} km)")

    def region(a, b):
        if a >= 37.4:
            return "경기·강원북"
        if a >= 36.6:
            return "충청·강원남"
        return "영남" if b >= 127.9 else "호남"

    regs = np.array([region(a, b) for a, b in zip(gn, gl)])
    tot = A.sum()
    print(f"\n{'권역':<12}{'면적%':>7}{'σ중앙':>8}{'권장간격':>9}{'배분점수':>9}{'균등대비':>9}")
    for rg in ["경기·강원북", "충청·강원남", "영남", "호남"]:
        s = regs == rg
        n_a = float((dens[s] * A[s]).sum())
        n_u = N_TARGET_SITES * A[s].sum() / tot
        print(f"{rg:<12}{100*A[s].sum()/tot:>7.1f}{np.nanmedian(sg[s]):>8.0f}"
              f"{np.median(spac[s]):>9.0f}{n_a:>9.1f}{n_a/n_u:>8.2f}x")

    validate_against_lmm(ag)


if __name__ == "__main__":
    import sys

    sys.stdout.reconfigure(encoding="utf-8")
    _self_check()
