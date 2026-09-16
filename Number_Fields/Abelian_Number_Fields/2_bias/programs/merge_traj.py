#!/usr/bin/env python3
"""merge_traj.py — merge per-node outputs of drh_scan_mp into one global trajectory.

usage:  python3 merge_traj.py part0.csv part1.csv ... [-o merged.csv] [--ref serial.csv]

The node files carry node-relative cumulative sums; this script prefix-sums the
node totals (math.fsum) and shifts each node's rows, producing a CSV identical
in format to the serial drh_scan.c output. pi(10^k) is checked against known
values where available.
"""
import sys, math, csv

KNOWN_PI = {10**5: 9592, 10**6: 78498, 10**7: 664579, 10**8: 5761455,
            10**9: 50847534, 10**10: 455052511, 10**11: 4118054813,
            10**12: 37607912018, 10**13: 346065536839,
            10**14: 3204941750802, 10**15: 29844570422669,
            10**16: 279238341033925}

def read_part(path):
    meta, rows, total, header = {}, [], None, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith("#NODE"):
                _, node, of, N, per = line.split(",")
                meta.update(node=int(node), of=int(of), N=int(N), per=int(per))
            elif line.startswith("#RANGE"):
                _, lo, hi = line.split(",")
                meta.update(lo=int(lo), hi=int(hi))
            elif line.startswith("#TOTAL"):
                p = line.split(",")
                total = (int(p[1]), [float(v) for v in p[2:]])
            elif line.startswith("x,"):
                header = line.split(",")
                meta["header"] = header
                continue
            else:
                p = line.split(",")
                rows.append((int(p[0]), int(p[1]), [float(v) for v in p[2:]]))
    if total is None:
        sys.exit(f"{path}: missing #TOTAL line (run incomplete?)")
    return meta, rows, total

def main():
    args = sys.argv[1:]
    out, ref = "merged.csv", None
    if "-o" in args:
        i = args.index("-o"); out = args[i+1]; del args[i:i+2]
    if "--ref" in args:
        i = args.index("--ref"); ref = args[i+1]; del args[i:i+2]

    parts = [read_part(p) for p in args]
    parts.sort(key=lambda t: t[0]["node"])
    of, N = parts[0][0]["of"], parts[0][0]["N"]
    if len(parts) != of:
        sys.exit(f"expected {of} node files, got {len(parts)}")
    for k, (m, _, _) in enumerate(parts):
        if m["node"] != k or m["of"] != of or m["N"] != N:
            sys.exit(f"inconsistent metadata in node {k}")
    for k in range(1, of):
        if parts[k][0]["lo"] != parts[k-1][0]["hi"]:
            sys.exit(f"ranges of nodes {k-1},{k} do not tile")

    # prefix over node totals, shift each node's rows
    header = parts[0][0].get("header") or ["x","pi_x","S4","S3","S5","S5c_re","S5c_im","S7"]
    nv = len(header) - 2
    off_cnt, off_s = 0, [0.0]*nv
    acc = [[] for _ in range(nv)]
    merged = []
    for m, rows, (tcnt, tsum) in parts:
        for x, c, s in rows:
            merged.append((x, c + off_cnt, [s[j] + off_s[j] for j in range(nv)]))
        off_cnt += tcnt
        for j in range(nv):
            acc[j].append(tsum[j])
            off_s[j] = math.fsum(acc[j])

    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for x, c, s in merged:
            w.writerow([x, c] + [f"{v:.15g}" for v in s])
    print(f"wrote {out}: {len(merged)} checkpoints, x up to {merged[-1][0]:,}")
    print(f"pi({N:,}) = {off_cnt:,}", end="")
    if N in KNOWN_PI:
        print("  [MATCHES known value]" if off_cnt == KNOWN_PI[N]
              else f"  [MISMATCH! expected {KNOWN_PI[N]:,}]")
    else:
        print()
    for x, c, _ in merged:
        if x in KNOWN_PI and KNOWN_PI[x] != c:
            print(f"  WARNING pi({x:,}) = {c:,} != known {KNOWN_PI[x]:,}")

    if ref:
        rrows = {}
        with open(ref) as f:
            for row in csv.DictReader(f):
                rrows[int(row["x"])] = (int(row["pi_x"]),
                    [float(row[k]) for k in header[2:]])
        dmax, cbad = 0.0, 0
        for x, c, s in merged:
            if x not in rrows: continue
            rc, rs = rrows[x]
            if rc != c: cbad += 1
            dmax = max(dmax, max(abs(s[j]-rs[j]) for j in range(nv)))
        print(f"vs {ref}: pi mismatches = {cbad},  max |Delta S| = {dmax:.3g}")

if __name__ == "__main__":
    main()
