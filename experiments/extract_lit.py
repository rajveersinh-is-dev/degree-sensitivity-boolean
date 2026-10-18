"""Extract relevant passages from the literature texts (Proskurin 2021, lowdeg 2021).
Writes UTF-8 output file only (console is cp1252).
"""
import re

LIT = r"C:\Users\davea\AppData\Local\Temp\opencode\lit"

out = []


def pr(*a):
    """Pr.
    
    """
    out.append(" ".join(str(x) for x in a))


t = open(LIT + r"\proskurin.txt", encoding="utf-8").read()
i = t.find("Fully Sensitive Function of Small Degree")
if i < 0:
    i = t.find("fully sensitive function of small degree")
pr("=" * 70)
pr(f"PROSKURIN SECTION 5 (offset {i}):")
pr(t[i:i + 5000] if i >= 0 else "NOT FOUND")

pr("=" * 70)
pr("proskurin: occurrences of the word '14':")
for m in re.finditer(r"\b14\b", t):
    s = max(0, m.start() - 150)
    e = min(len(t), m.end() + 150)
    pr(repr(t[s:e])[:400])
    pr("---")

t2 = open(LIT + r"\lowdeg.txt", encoding="utf-8").read()
pr("=" * 70)
pr("lowdeg (arXiv 2107.11205): occurrences of '14' and key claims:")
for pat in [r"\b14\b", r"[7-9] variables", r"sensitivity 7", r"seven"]:
    pr(f"-- pattern {pat}:")
    cnt = 0
    for m in re.finditer(pat, t2):
        s = max(0, m.start() - 120)
        e = min(len(t2), m.end() + 120)
        pr(repr(t2[s:e])[:300])
        cnt += 1
        if cnt >= 6:
            break

with open(r"C:\Users\davea\Downloads\newp\lit_extract_out.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(out) + "\n")
print("done", len(out), "lines")
