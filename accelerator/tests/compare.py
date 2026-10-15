#!/usr/bin/env python3
"""Compare the H that h_construct emits against the reference H from build_ball_graph.

Usage:  compare.py <vector-dir> <num-shots>

For shot k it reads ref_<k>.txt (software) and dut_<k>.txt (RTL), canonicalises both and
diffs them. Canonical form is a sorted list of records, so a duplicate record on one side is
a mismatch rather than being silently absorbed:

    ("E", min(u, v), max(u, v), w)      a defect-defect edge
    ("B", u, w)                        a boundary edge

The hardware emits edges in whatever order the scan reaches them, so order is not compared;
multiplicity is. Exits non-zero if any shot mismatches.

Needs nothing outside the standard library.
"""

import os
import sys

# Differing records printed per failing shot before the listing is truncated.
MAX_SHOWN = 6


def parse(path):
    """Canonical sorted record list for one H dump file."""
    records = []
    with open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if fields[0] == "E" and len(fields) == 4:
                u, v, w = (int(x) for x in fields[1:])
                records.append(("E", min(u, v), max(u, v), w))
            elif fields[0] == "B" and len(fields) == 3:
                u, w = (int(x) for x in fields[1:])
                records.append(("B", u, w))
            else:
                raise ValueError("%s:%d: cannot parse %r" % (path, lineno, line))
    records.sort()
    return records


def show(record):
    return " ".join(str(f) for f in record)


def diff(ref, dut):
    """Records present in exactly one side, as (only_in_ref, only_in_dut).

    Both inputs are sorted, so this is a merge walk and it keeps multiplicity: a record that
    appears twice on one side and once on the other leaves one copy in the result.
    """
    only_ref, only_dut = [], []
    i = j = 0
    while i < len(ref) and j < len(dut):
        if ref[i] == dut[j]:
            i += 1
            j += 1
        elif ref[i] < dut[j]:
            only_ref.append(ref[i])
            i += 1
        else:
            only_dut.append(dut[j])
            j += 1
    only_ref.extend(ref[i:])
    only_dut.extend(dut[j:])
    return only_ref, only_dut


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    vec_dir, num_shots = argv[1], int(argv[2])

    passed = 0
    failures = []
    for k in range(num_shots):
        ref_path = os.path.join(vec_dir, "ref_%d.txt" % k)
        dut_path = os.path.join(vec_dir, "dut_%d.txt" % k)
        if not os.path.exists(dut_path):
            print("shot %-3d FAIL  no DUT output (%s missing)" % (k, dut_path))
            failures.append(k)
            continue

        ref, dut = parse(ref_path), parse(dut_path)
        only_ref, only_dut = diff(ref, dut)
        if not only_ref and not only_dut:
            n_edge = sum(1 for r in ref if r[0] == "E")
            print(
                "shot %-3d PASS  %d edges, %d boundary" % (k, n_edge, len(ref) - n_edge)
            )
            passed += 1
            continue

        failures.append(k)
        print(
            "shot %-3d FAIL  ref %d records, dut %d records, %d only in ref, %d only in dut"
            % (k, len(ref), len(dut), len(only_ref), len(only_dut))
        )
        for record in only_ref[:MAX_SHOWN]:
            print("           ref only: %s" % show(record))
        if len(only_ref) > MAX_SHOWN:
            print("           ... %d more only in ref" % (len(only_ref) - MAX_SHOWN))
        for record in only_dut[:MAX_SHOWN]:
            print("           dut only: %s" % show(record))
        if len(only_dut) > MAX_SHOWN:
            print("           ... %d more only in dut" % (len(only_dut) - MAX_SHOWN))

    if failures:
        print(
            "CROSSCHECK FAIL %d/%d shots; mismatching shots: %s"
            % (passed, num_shots, ", ".join(str(k) for k in failures))
        )
        return 1
    print("CROSSCHECK PASS %d/%d shots" % (passed, num_shots))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
