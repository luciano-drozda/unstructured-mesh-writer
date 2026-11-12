#!/usr/bin/env python3
"""List every top-level function (name, argument count) and top-level
constant in the given Julia files. Report a name that appears twice with the
same argument count (a collision). Exit code 1 if any collision exists.
With one file: duplicates inside that file. With several files: duplicates
across the files and inside each file."""
import re, sys, collections

def count_args(sig):
    # sig: text between the outer parentheses of a function header
    depth = 0; n = 0; seen = False; stop = False
    for ch in sig:
        if ch in "([{": depth += 1
        elif ch in ")]}": depth -= 1
        elif ch == ";" and depth == 0: break
        elif ch == "," and depth == 0: n += 1
        if not ch.isspace(): seen = True
    return (n + 1) if seen else 0

def scan(path):
    funcs, consts = [], []
    lines = open(path, encoding="utf-8").read().split("\n")
    i = 0
    while i < len(lines):
        ln = lines[i]
        m = re.match(r"^function\s+([A-Za-z_][A-Za-z0-9_!]*)\s*\(", ln)
        if m:
            text = ln; j = i
            depth = text.count("(") - text.count(")")
            while depth > 0 and j + 1 < len(lines):
                j += 1; text += "\n" + lines[j]
                depth = text.count("(") - text.count(")")
            a = text.index("(")
            d = 0
            for k in range(a, len(text)):
                if text[k] == "(": d += 1
                elif text[k] == ")":
                    d -= 1
                    if d == 0: break
            funcs.append((m.group(1), count_args(text[a + 1:k]), i + 1))
        m = re.match(r"^(?:const\s+)?([A-Za-z_][A-Za-z0-9_]*)\s*=[^=]", ln)
        if m and not ln.startswith("function"):
            consts.append((m.group(1), i + 1))
        m = re.match(r"^(?:mutable\s+)?struct\s+([A-Za-z_][A-Za-z0-9_]*)", ln)
        if m: consts.append((m.group(1), i + 1))
        i += 1
    return funcs, consts

paths = sys.argv[1:]
fseen = collections.defaultdict(list)
cseen = collections.defaultdict(list)
for p in paths:
    f, c = scan(p)
    for name, nargs, ln in f: fseen[(name, nargs)].append((p.split("/")[-1], ln))
    for name, ln in c: cseen[name].append((p.split("/")[-1], ln))
bad = 0
for (name, nargs), where in sorted(fseen.items()):
    if len(where) > 1:
        bad += 1
        print(f"COLLISION function {name}/{nargs}: " + ", ".join(f"{f}:{l}" for f, l in where))
for name, where in sorted(cseen.items()):
    if len(where) > 1:
        bad += 1
        print(f"COLLISION constant {name}: " + ", ".join(f"{f}:{l}" for f, l in where))
byname = collections.defaultdict(set)
for (name, nargs) in fseen: byname[name].add(nargs)
for name, s in sorted(byname.items()):
    if len(s) > 1: print(f"note: {name} has several arities {sorted(s)}")
nf = sum(len(v) for v in fseen.values()); print(f"{nf} functions, {len(cseen)} constants, {bad} collisions")
sys.exit(1 if bad else 0)
