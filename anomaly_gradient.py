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

# ── 배분 대상 총 측점 수 ──────────────────────────────────────────────
#
# ⚠️ **분모와 분자는 같은 모집단이어야 한다** (2026-08-27 정정).
# 종전에는 목표를 66 = 「현 LMM 투입 16 + 신규 50」으로 잡아 놓고, 충족도의
# 「현 간격」은 **1등 지자기점 관측망 30점**으로 쟀다. 목표 모집단(LMM 투입망)과
# 현황 모집단(관측망)이 달라 비율이 뜻을 갖지 못했다.
#
# 모집단을 **1등 지자기점 반복관측망** 하나로 통일한다:
#     목표 = 현 관측망 + Track C 신규 50점
# LMM 투입 16점은 「그중 현재 모델에 들어간 자료」이지 관측망이 아니다.
#
# ⚠️ 「현 관측망」은 **좌표를 실제로 가진 측점 수**다. 2026-08-27 이전에는
#    남양·서산이 같은 좌표라 공간적으로 29점이었으나, 지리원 원장 좌표를
#    적용해 중복이 해소돼 **30점**이 됐다.
#
# **권역 간 비(比)는 N 과 무관하다.** ρ(x) ∝ σ(x) 이고 L = 1/√ρ 이므로
#    L ∝ 1/√N — N 을 바꾸면 전국이 같은 배율로 늘거나 줄 뿐 어느 권역이 몇 배
#    조밀해야 하는가는 바뀌지 않는다.
#
# ⚠️ **그러나 그것을 「확정 사실」로 쓰지 말 것.** N 에 무관하다는 것은 계산이
#    N 에 딸리지 않는다는 뜻이지 결론이 검증됐다는 뜻이 아니다 — σ 결측
#    (항공자력 공백) 위에 서 있고 Neyman 배분 가정 자체가 미검증이다.
#    정확한 표현은 **「N 가정에 독립적인 잠정 설계비」**이며, 절대 km 값과
#    부족도는 그보다도 더 잠정적이다.
N_NEW_SITES = 50               # Track C 신규 (35~40 구축용 + 10~15 검증 전용)


def _network_size() -> int:
    """
    **공간적으로 쓸 수 있는** 현 관측망 크기.

    ⚠️ 이것은 **기본값일 뿐**이다. 실제 계산에서는 좌표를 실제로 가진 측점
    수를 세어 목표 N 을 잡는다(`compute_density_design` 이 `n_sites` 로
    유도). 좌표가 겹치는 측점은 담당면적을 하나도 못 받아 공간 계산에서
    사라지므로, 행 수를 그대로 분모로 쓰면 모집단이 갈라진다.
    2026-08-27 에 지리원 원장 좌표를 적용해 중복이 해소돼 30 이 됐다.
    """
    try:
        from existing_network import EXISTING_NETWORK
        return len(EXISTING_NETWORK)
    except Exception:
        return 30


N_TARGET_SITES = _network_size() + N_NEW_SITES      # = 80

# 권장 측점 간격 등급 (km)
SPACING_TIERS = [
    ("dense",  0.0,  35.0, "조밀 배치 필요"),
    ("mid",   35.0,  50.0, "중간"),
    ("coarse", 50.0, 65.0, "보통"),
    ("sparse", 65.0, np.inf, "성김 허용"),
]

# ── 밀도 부족도 ───────────────────────────────────────────────────────
#
#   충족도(raw)   = 현 관측망 등가 간격 / 권장 간격
#   상대 부족도    = 충족도 / sqrt(N_target / N_now)
#
# ⚠️ **raw 충족도의 «수준»은 산술이 정한다** (2026-08-27 정정). 현 30점을
# 목표 80점 설계와 견주면 전국이 평균 sqrt(80/30) = 1.63 배 성긴 것이 당연하고,
# 「93 %가 부족」은 발견이 아니라 30→80 이라는 뺄셈이다. 지도가 실제로 말하는
# 것은 **어디가 상대적으로 더 부족한가**이므로 등급은 그 균일부족선으로
# 나눈 **상대 부족도**에 매긴다. 1.0 이 「전국 평균만큼 부족」이다.
#
# 이렇게 하면 등급이 N 가정에 딸리지 않는다 — N 을 바꿔도 raw 충족도와
# 균일부족선이 같은 배율로 움직여 상대 부족도는 그대로다.
DEFICIT_TIERS = [
    ("ok",       0.0, 0.7, "충분 — 평균보다 크게 양호"),
    ("mild",     0.7, 1.0, "평균보다 양호"),
    ("short",    1.0, 1.4, "평균보다 부족"),
    ("critical", 1.4, np.inf, "크게 부족 — 보강 1순위"),
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
        **반경 r_km 원 안**의 자기이상 표준편차 (nT) — 권역 자기복잡도.

        `detrend=True` 면 국소 1차 평면을 뺀다. 광역 추세는 Regional 층이
        담당하므로 밀도 설계의 근거에서 제외하는 것이 옳다(실측상 차이는
        3% 로 작지만 정의가 명확해진다).

        ⚠️ **경위도 상자로 자르면 반경이 아니다** (2026-08-27 정정). 종전에는
        |Δlat|·|Δlon| 을 각각 비교해 약 2r × 2r **사각창**을 썼고, 모서리가
        중심에서 r·√2 = 35.4 km 나 떨어졌다. 실제 원형 창으로 바꾸면
        σ 중앙값 88.2 → 84.5 nT, 결측 44 → 57셀, 권장간격 등급이 1,090셀 중
        125셀(11.5 %) 바뀐다. 상자는 먼저 걸러 내는 용도로만 쓰고 반드시
        원형 거리로 다시 자를 것.
        """
        lat = np.atleast_1d(np.asarray(lat, float))
        lon = np.atleast_1d(np.asarray(lon, float))
        out = np.full(lat.shape, np.nan)

        for k, (la, lo) in enumerate(zip(lat, lon)):
            coslat = np.cos(np.radians(la))
            dlat = r_km / 110.574
            dlon = r_km / (111.320 * coslat)
            # 1차: 상자로 후보를 줄인다(전수 거리계산 회피)
            box = ((np.abs(self.plat - la) <= dlat) &
                   (np.abs(self.plon - lo) <= dlon))
            if not box.any():
                continue
            x = (self.plon[box] - lo) * 111.320 * coslat
            y = (self.plat[box] - la) * 110.574
            # 2차: 실제 원형 반경으로 자른다
            inr = (x * x + y * y) <= (r_km * r_km)
            if inr.sum() < SIGMA_MIN_PTS:
                continue
            v = self.panom[box][inr]
            if not detrend:
                out[k] = float(np.std(v))
                continue
            A = np.column_stack([np.ones(inr.sum()), x[inr], y[inr]])
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


def catchment_spacing(qlat, qlon, area_km2, slat, slon):
    """
    기존 측점망의 **등가 측점 간격** L_now (km) — 권장 간격과 같은 정의.

    ⚠️ **「최근접 측점 거리 × 2」를 쓰면 안 된다** (2026-08-27 정정).
    그 값은 측점 위에서 0 으로 무너지고 측점 사이·외곽에서 성질이 달라져,
    권장 간격과 나눠 「충족/부족」을 판정할 수 있는 양이 아니었다. 실제로
    측점 근처 셀이 자동으로 「충족」으로 찍혔다.
    ⚠️ k-최근접 밀도 추정도 쓰지 않는다 — 20 km 정사각 격자로 검산하면 셀
    중심에서 14.4 km, 측점 위에서 20.5 km 로 **위치에 따라 ±28 % 흔들린다.**

    권장 간격은 밀도의 함수 L_req = 1/√ρ_req 이므로, 현황도 **밀도**에서
    같은 방식으로 유도해야 비율이 뜻을 갖는다. 여기서는 **담당 면적**을 쓴다
    (국토 격자 위의 이산 보로노이):

        각 격자셀을 가장 가까운 측점에 배정
        → 측점 i 의 담당 면적 A_i = Σ(배정된 셀 면적)
        → 그 셀의 등가 간격 L_now = √(A_i)

    정규 격자 배치에서 정확히 L 을 돌려주고, 국토 경계로 자동으로 잘리며
    (해상은 격자에 없다), 「이 측점이 실제로 담당하는 면적 대 담당해야 할
    면적」이라는 Neyman 배분과 같은 물음이 된다.

    ⚠️ 값은 **측점 담당구역 단위로 계단**이다(보로노이 경계에서 불연속).
    격자 해상도만큼 면적이 양자화되므로 0.1° 격자에서 간격 오차는 약 1 %.

    Returns
    -------
    (spacing_km, cover_km, area_of_cell, area_by_site)
      spacing_km   : 셀별 등가 측점 간격
      cover_km     : **최근접 측점 거리**(coverage radius) — 「가장 가까운 측점이
                     얼마나 먼가」라는 별개의 물음. 이름을 나눠 함께 돌려준다.
      area_of_cell : 셀이 속한 측점의 담당 면적 (길이 = 격자 수)
      area_by_site : **측점별** 담당 면적 (길이 = 측점 수)

    ⚠️ **`area_of_cell` 로 평균을 내지 말 것** — 셀마다 소유 측점의 전체
    담당면적이 반복돼 있어 `mean()` 이 **셀 수로 가중**된다(큰 담당구역이
    과대 반영). 측점당 통계는 `area_by_site` 를 쓴다. 담당 셀을 하나도 못 받은
    측점은 0 이므로 `area_by_site > 0` 으로 걸러야 한다 — 좌표가 겹치는 측점이
    실제로 그렇게 된다(`existing_network.drop_duplicate_coords` 로 사전 제거).
    """
    qlat = np.atleast_1d(np.asarray(qlat, float))
    qlon = np.atleast_1d(np.asarray(qlon, float))
    area = np.atleast_1d(np.asarray(area_km2, float))
    slat = np.asarray(slat, float)
    slon = np.asarray(slon, float)
    n = slat.size
    if n == 0:
        nan = np.full(qlat.shape, np.nan)
        return nan, nan.copy(), nan.copy(), np.zeros(0)

    d = np.empty((qlat.size, n))
    for j in range(n):
        m = np.radians((qlat + slat[j]) / 2.0)
        d[:, j] = np.hypot((qlon - slon[j]) * 111.320 * np.cos(m),
                           (qlat - slat[j]) * 110.574)
    owner = np.argmin(d, axis=1)
    cover = d[np.arange(qlat.size), owner]

    site_area = np.zeros(n)
    np.add.at(site_area, owner, area)
    a_of_cell = site_area[owner]
    return np.sqrt(a_of_cell), cover, a_of_cell, site_area


def uniform_deficit(n_now: int, n_target: int = N_TARGET_SITES) -> float:
    """
    **균일부족선** — 현 n_now 점이 목표 n_target 점 설계에 견줘 전국 평균으로
    몇 배나 성긴가. 밀도가 N 에 비례하고 간격이 1/sqrt(밀도) 이므로
    sqrt(n_target / n_now) 다.

    상대 부족도 = 충족도 / 이 값. 1.0 이 「전국 평균만큼 부족」이다.
    """
    if n_now <= 0:
        return float("nan")
    return float(np.sqrt(n_target / n_now))


def uniform_spacing(area_km2_total, n_target: int = N_TARGET_SITES):
    """비교 기준 — 같은 측점 수를 균등 배치했을 때의 간격 (km)."""
    return float(np.sqrt(area_km2_total / n_target))


DENSITY_CELL_DEG = 0.10        # 권역 밀도 격자 (약 11 km)


def land_grid(korea_geojson=None, cell_deg: float = DENSITY_CELL_DEG):
    """
    국토 내부 평가 격자 (clat, clon, area_km2) — **밀도 설계의 단일 격자**.

    ⚠️ 자기점검과 본 분석이 서로 다른 격자를 쓰면 셀 수·등급 수가 갈라져
    「같은 정의를 쓴다」는 말이 성립하지 않는다(2026-08-27 검토 지적:
    자기점검 1,088셀 vs 배포 1,090셀). 두 곳이 이 함수를 함께 쓴다.

    격자 원점·범위는 **경계 bbox 를 cell_deg 로 내림**해 잡는다 —
    `np.arange(33.1, ...)` 처럼 손으로 적으면 경계가 바뀔 때 갈라진다.
    """
    import json as _json

    from matplotlib.path import Path as _MPath

    # 원본(data/)을 먼저 본다 — docs/ 는 배포 사본이라 빌드 전후로 갈릴 수 있다
    path = korea_geojson
    if path is None:
        for cand in (ROOT / "data" / "korea_boundary.geojson",
                     ROOT / "docs" / "data" / "korea_boundary.geojson"):
            if cand.exists():
                path = cand
                break
    gj = _json.load(open(path, encoding="utf-8"))
    polys = []
    for f in gj["features"]:
        g = f["geometry"]
        for poly in ([g["coordinates"]] if g["type"] == "Polygon"
                     else g["coordinates"]):
            polys.append(_MPath(np.asarray(poly[0])))

    xs = np.concatenate([p.vertices[:, 0] for p in polys])
    ys = np.concatenate([p.vertices[:, 1] for p in polys])
    glon = np.arange(np.floor(xs.min() / cell_deg) * cell_deg,
                     xs.max() + cell_deg, cell_deg)
    glat = np.arange(np.floor(ys.min() / cell_deg) * cell_deg,
                     ys.max() + cell_deg, cell_deg)
    LO, LA = np.meshgrid(glon, glat)
    pts = np.column_stack([LO.ravel(), LA.ravel()])
    m = np.zeros(len(pts), bool)
    for pp in polys:
        m |= pp.contains_points(pts)
    clat, clon = LA.ravel()[m], LO.ravel()[m]
    return clat, clon, cell_area_km2(clat, cell_deg, cell_deg)


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

    # 국토 격자에서 Neyman 배분 — **본 분석과 같은 격자**를 쓴다
    data = ROOT / "docs" / "data"
    gn, gl, A = land_grid()

    sg = ag.sigma(gn, gl)
    dens, spac, _ = neyman_density(sg, A)
    print(f"\n국토 격자 {gn.size}점 · σ(R={SIGMA_RADIUS_KM:.0f}km) "
          f"결측 {int(np.isnan(sg).sum())}")
    print("σ nT  P10/50/90:", np.round(np.nanpercentile(sg, [10, 50, 90]), 0))
    print(f"권장 간격 km  P10/50/90: {np.round(np.percentile(spac,[10,50,90]),0)}"
          f"   (균등 배치 기준 {uniform_spacing(A.sum()):.0f} km)")
    print(f"배분 총점수 N = {N_TARGET_SITES} "
          f"(공간 관측망 {_network_size()} + 신규 {N_NEW_SITES})")

    # 현 관측망 등가 간격 — 담당면적(이산 보로노이) 방식
    try:
        import pandas as _pd

        from existing_network import select_rows as _sel
        ex = _pd.DataFrame(
            [{"도엽명": f["properties"]["name"],
              "위도": f["properties"]["lat"], "경도": f["properties"]["lon"]}
             for f in json.load(open(data / "existing_sites.geojson",
                                     encoding="utf-8"))["features"]])
        from existing_network import drop_duplicate_coords as _dedup
        net_all = _sel(ex)
        net, dropped = _dedup(net_all)
        if dropped:
            print(f"  ⚠ 좌표 중복으로 제외: {' · '.join(dropped)} "
                  f"(명목 {len(net_all)} → 공간 {len(net)}점)")
        sp_now, cover, a_cell, a_site = catchment_spacing(
            gn, gl, A, net["위도"].values, net["경도"].values)
        n_now = len(net)
        n_tgt = n_now + N_NEW_SITES
        # ⚠️ 등급은 **상대 부족도**에 매긴다 — raw 를 그대로 넣으면 「전국이
        #    보강 1순위」라는 뻔한 결과가 나온다(수준은 n_now→n_tgt 산술이 정한다).
        dfc = sp_now / spac
        dfr = dfc / uniform_deficit(n_now, n_tgt)
        used = a_site > 0
        print(f"\n현 관측망 {n_now}점 · 담당면적 등가 간격 중앙 "
              f"{np.median(sp_now):.0f} km "
              f"(측점당 {a_site[used].mean():.0f} km2 · 국토 {A.sum():.0f} km2)")
        print(f"  [별도] 최근접 측점 거리 중앙 {np.median(cover):.0f} km "
              f"· 최대 {cover.max():.0f} km")
        print(f"  균일부족선 sqrt({n_tgt}/{n_now}) = "
              f"{uniform_deficit(n_now, n_tgt):.2f} · raw 중앙 "
              f"{np.median(dfc):.2f} → 상대 부족도 중앙 {np.median(dfr):.2f}")
        print("  상대 부족도: " + " · ".join(
            f"{TIER_LABEL[k]} {int(((dfr >= lo) & (dfr < hi)).sum())}"
            for k, lo, hi, _lab in DEFICIT_TIERS))
    except Exception as exc:
        print(f"\n[건너뜀] 현 관측망 등가 간격: {exc}")


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
