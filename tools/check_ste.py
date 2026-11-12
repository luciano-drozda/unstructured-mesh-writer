#!/usr/bin/env python3
"""Scan the comments of a Julia file for common Simplified Technical English
(ASD-STE100) problems. Usage: check_ste.py FILE... Exit code 1 on findings.
This is a partial check. It tests what a script can test:
  - a sentence of more than 25 words
  - a contraction (don't, it's, ...)
  - a word or phrase that STE does not allow in plain form: etc, e.g., i.e.,
    via, and/or, should, approximately, utilize
It does not test the dictionary, the grammar, or the passive voice."""
import re, sys
BANNED = [r"\betc\b\.?", r"\be\.g\.", r"\bi\.e\.", r"\bvia\b", r"and/or", r"\bshould\b", r"\bapproximately\b", r"\butili[sz]e\b"]
CONTR = r"\b\w+'(t|s|re|ve|ll|d|m)\b"

def comment_text(text):
    """Return a list of (first_line_no, text) blocks: each block is one run of
    comment lines, or the text of a #= ... =# block, with code-like lines removed."""
    out = []
    lines = text.split("\n"); i = 0
    while i < len(lines):
        l = lines[i]
        if l.startswith("#="):
            j = i
            while not lines[j].startswith("=#"): j += 1
            out.append((i + 1, "\n".join(lines[i+1:j]))); i = j + 1; continue
        if l.lstrip().startswith("#") and not l.startswith("#!") :
            j = i; blk = []
            while j < len(lines) and lines[j].lstrip().startswith("#") and not re.match(r"^\s*#\s*[=─]{5,}", lines[j]):
                blk.append(re.sub(r"^\s*#\s?", "", lines[j])); j += 1
            if blk: out.append((i + 1, "\n".join(blk)))
            i = max(j, i + 1); continue
        # trailing comment on a code line
        m = re.search(r"\s#\s(.*)$", l)
        if m and '"' not in l[:m.start()].split("#")[0][-1:]:
            out.append((i + 1, m.group(1)))
        i += 1
    return out

def sentences(block):
    """Split one comment block into sentences. An argument line 'name: text',
    a blank line, a table line, and a function-signature line start a new item."""
    items = []; cur = []
    def flush():
        if cur: items.append(" ".join(cur)); cur.clear()
    for ln in block.split("\n"):
        raw = ln.rstrip()
        if not raw.strip():
            flush(); continue
        if re.match(r"^\s*[-=─]{3,}", raw): flush(); continue
        if re.match(r"^\s{4,}\S", raw) and not re.match(r"^\s{4,}[A-Za-z].*[a-z]\.$", raw):
            flush(); continue                                           # usage and table lines
        if re.match(r"^\s{0,4}\d+, \w+\s{2,}\S", raw):
            flush(); continue                                           # table row of a return code
        if re.match(r"^[A-Za-z_][\w\[\], ]*:\s", raw) or re.match(r"^[a-z_0-9]+\(.*\)$", raw.strip()) \
           or re.match(r"^\s{0,3}[A-Za-z_][\w\[\], ]*\s{2,}\S", raw):
            flush()
        cur.append(raw.strip())
    flush()
    out = []
    for it in items:
        out += [p for p in re.split(r"(?<=[.!?])\s+(?=[A-Z])", it) if p]
    return out

nf = 0
for path in sys.argv[1:]:
    text = open(path, encoding="utf-8").read()
    n = 0
    for ln, blk in comment_text(text):
        for s in sentences(blk):
            words = re.findall(r"[A-Za-z0-9_'\-]+", s)
            if len(words) > 25:
                n += 1; print(f"{path.split('/')[-1]}:{ln}: long sentence ({len(words)} words): {s[:90]}...")
            if re.search(CONTR, s) and not re.search(r"\w+'s\b(?=\s+(own|value))", s):
                m = re.search(CONTR, s)
                if m and m.group(1) != "s":
                    n += 1; print(f"{path.split('/')[-1]}:{ln}: contraction '{m.group(0)}'")
            for pat in BANNED:
                if re.search(pat, s, re.I):
                    n += 1; print(f"{path.split('/')[-1]}:{ln}: banned word/phrase {pat}: {s[:80]}")
    print(f"{path.split('/')[-1]}: {n} findings"); nf += n
sys.exit(1 if nf else 0)
