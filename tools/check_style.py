#!/usr/bin/env python3
"""Scan the computing-kernel block of a Julia file for skill-stade violations.
Usage: check_style.py FILE...   Exit code 1 if any violation is found.
The block lies between the 'computing kernels' banner and the
'end of computing kernels' banner. Heuristic checks on text, not a parser."""
import re, sys
ALLOWED = {"abs","sqrt","exp","log","log10","sin","cos","tan","asin","acos","atan","sinh","cosh","tanh",
           "mod","div","max","min","sign","floor","ceil","trunc"}
KEYWORDS = {"if","elseif","while","for","function","return","end","else"}

def kernel_block(text):
    a = text.index("# computing kernels")
    b = text.index("# end of computing kernels")
    return text[a:b]

def split_functions(block):
    L = block.split("\n"); out = []
    i = 0
    while i < len(L):
        m = re.match(r"^function\s+([A-Za-z_][A-Za-z0-9_!]*)\s*\(", L[i])
        if m:
            j = i
            while L[j] != "end": j += 1
            h = i
            while h > 0 and L[h-1].startswith("#"): h -= 1
            out.append((m.group(1), L[h:i], L[i:j+1], i + 1)); i = j
        i += 1
    return out

def strip_comment(s):
    return re.sub(r"#.*$", "", s)

def check_function(name, header, body):
    errs = []
    sig_lines = []
    k = 0
    depth = 0
    while True:
        sig_lines.append(body[k]); depth += body[k].count("(") - body[k].count(")")
        if depth <= 0: break
        k += 1
    sig = " ".join(sig_lines)
    code = [strip_comment(l) for l in body[k+1:-1]]
    ctext = "\n".join(code)
    # 5 keyword arguments
    inner = sig[sig.index("(")+1: sig.rindex(")")]
    if ";" in inner: errs.append("keyword-argument section in the signature")
    # 12 header
    if not header or not header[0].startswith("# " + name + "("):
        errs.append("comment header does not start with '# name(args)'")
    if '"""' in "\n".join(body) or any('"""' in h for h in header): errs.append("docstring")
    # 14 return nothing
    if body[-2].strip() != "return nothing": errs.append("last statement is not 'return nothing'")
    # 7 allocation
    for pat in (r"\bzeros\(", r"\bones\(", r"\bsimilar\(", r"\bArray\{", r"\bVector\{", r"\bMatrix\{", r"\bcollect\(", r"\bfill\(", r"\bcopy\("):
        if re.search(pat, ctext): errs.append(f"allocation ({pat})")
    if re.search(r"\[[^\]\n]*\bfor\b[^\]\n]*\]", ctext): errs.append("array comprehension")
    # 8 inline indexing a[b[i]]
    if re.search(r"[A-Za-z_]\w*\[[^\[\]\n]*[A-Za-z_]\w*\[", ctext): errs.append("inline indexing a[b[i]]")
    # 9 broadcast
    if re.search(r"(?<![0-9])\.[+\-*/^=]", ctext) or re.search(r"[A-Za-z_]\w*\.\(", ctext): errs.append("broadcast operator")
    # 10 allowed functions
    for m in re.finditer(r"(?<![A-Za-z0-9_.])([A-Za-z_][A-Za-z0-9_!]*)\(", ctext):
        n = m.group(1)
        if n in KEYWORDS or n in ALLOWED: continue
        errs.append(f"call to non-whitelisted function '{n}'")
    if re.search(r"\bifelse\b", ctext): errs.append("ifelse")
    # 11 integer division
    if re.search(r"÷|\bfld\(|Int\(floor|Int\(\s*floor", ctext): errs.append("integer division not via div")
    if re.search(r"÷=|%=|\\=", ctext): errs.append("forbidden compound assignment")
    # 6 stored comparison
    for l in code:
        if re.match(r"^\s*(if|elseif|while|for|return)\b", l): continue
        if re.match(r"^\s*[A-Za-z_]\w*(\[[^\]]*\])?\s*=[^=]", l) and re.search(r"==|!=|<=|>=|(?<![<-])<(?![<=])|(?<![>-])>(?![>=])", l.split("=", 1)[1]):
            errs.append(f"stored comparison: {l.strip()}")
        if re.search(r"\btrue\b|\bfalse\b", l): errs.append(f"Bool literal: {l.strip()}")
    # 4 conditionals only where needed is a judgement; count them
    # 2 unfused rectangular nest of two independent loops: inner for is the only statement
    ind = lambda s: len(s) - len(s.lstrip())
    stmts = [(i, l) for i, l in enumerate(code) if l.strip()]
    for pos, (i, l) in enumerate(stmts):
        if re.match(r"^\s*for\b", l) and pos + 1 < len(stmts) and re.match(r"^\s*for\b", stmts[pos+1][1]) and ind(stmts[pos+1][1]) > ind(l):
            # find the end of the outer loop and the end of the inner loop
            def end_of(idx):
                d = ind(code[idx])
                for jj in range(idx + 1, len(code)):
                    if code[jj].strip() == "end" and ind(code[jj]) == d: return jj
            eo, ei = end_of(i), end_of(stmts[pos+1][0])
            rest = [s for s in code[ei+1:eo] if s.strip()]
            if not rest:
                errs.append(f"unfused rectangular nest at: {l.strip()}")
    # a kernel must not call another kernel
    for m in re.finditer(r"\b([A-Za-z_]\w*_kernel)\(", ctext):
        errs.append(f"calls kernel {m.group(1)}")
    return errs

bad = 0
for path in sys.argv[1:]:
    text = open(path, encoding="utf-8").read()
    fs = split_functions(kernel_block(text))
    nb = 0
    for name, header, body, ln in fs:
        for e in check_function(name, header, body):
            nb += 1; print(f"{path.split('/')[-1]}: {name}: {e}")
    print(f"{path.split('/')[-1]}: {len(fs)} kernels scanned, {nb} violations")
    bad += nb
sys.exit(1 if bad else 0)
