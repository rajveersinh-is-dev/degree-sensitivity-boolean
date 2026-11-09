import subprocess
import sys
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable

def main():
    cmd = [PY, "-u", os.path.join(ROOT, "experiments", "scan_level_d.py"), "16"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)

    out = r.stdout + r.stderr
    print(out)

    if "feasible=True" in out and "d=16" in out:
        print("PASS d=16 test")
        sys.exit(0)
    else:
        print("FAIL d=16 test")
        sys.exit(1)

if __name__ == "__main__":
    main()
