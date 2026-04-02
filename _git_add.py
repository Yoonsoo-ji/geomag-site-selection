"""Stage all untracked/modified project files for commit."""
import subprocess, sys, os

os.chdir(r"C:\LG_gram_backup_users\LX\2026_geomag")

files_to_add = [
    "data/국가기본도_도엽인덱스50K/",
    "data/국가기본도_도엽인덱스50K.zip",
    "data/국가기본도_도엽인덱스50K_테이블정의서.xlsx",
    "data/국가기본도_도엽인덱스25K_테이블정의서.xlsx",
    "data/수치자력이상도를+사용할+때+주의할+사항.pdf",
]

for f in files_to_add:
    r = subprocess.run(["git", "add", f], capture_output=True, text=True, encoding="utf-8")
    if r.returncode == 0:
        print(f"  OK staged: {f}")
    else:
        print(f"  FAIL: {f}\n    {r.stderr.strip()}")

# Final status
r = subprocess.run(["git", "status", "--short"], capture_output=True, text=True, encoding="utf-8")
print("\nGit status:")
print(r.stdout)
