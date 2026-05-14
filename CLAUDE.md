# 지구자기장 측정 입지 선정 프로젝트 — Claude 작업 규칙

## 출력 파일 저장 경로

**모든 작업 산출물은 두 곳에 동시 저장한다.**

| 경로 | 용도 |
|---|---|
| `docs/output/` (스크립트 기준 상대경로) | git 커밋 대상, GitHub Pages 배포 |
| `C:/LG_gram_backup_users/LX/2026_geomag/docs/output/` | 로컬 백업 |

- 스크립트 내 저장 순서: ① `docs/output/` 에 `doc.save()` → ② `shutil.copy2()` 로 로컬 경로 복사
- `docs/output/` 경로는 `Path(__file__).parent / "docs" / "output"` 패턴 사용 (절대경로 하드코딩 금지)
- 루트의 `output/` 단독 저장 금지

## 파일명 규칙

- 날짜/시간이 포함된 파일을 저장할 때는 반드시 `YYYYMMDD_HHMMSS` 접두사 사용
  - 예: `20260514_102055_site_selection_methodology_geology.docx`

## Git 커밋 규칙

- 커밋 대상: `docs/` 폴더 전체 (HTML, GeoJSON, output 산출물 포함)
- `output/` 루트 폴더는 커밋하지 않는다
- push는 `git push origin HEAD:main` 사용 (브랜치명 불일치 방지)

## 주요 스크립트

| 스크립트 | 역할 |
|---|---|
| `geomag_site_selection.py` | 메인 분석 + Folium 지도 생성 |
| `build_docs.py` | `output/` → `docs/` 빌드 (GeoJSON 단순화 포함) |
| `create_methodology_doc.py` | 기술 방법론 Word 문서 생성 → `docs/output/` 저장 |

## Python 환경

- 인터프리터: `C:\Users\YOONS\anaconda3\python.exe`
- 실행 예시: `& "C:\Users\YOONS\anaconda3\python.exe" create_methodology_doc.py`
