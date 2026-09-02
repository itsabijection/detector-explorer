#!/usr/bin/env python3
"""Generate an interactive HTML explorer for detector bases of a stim circuit.

Reads a stim circuit, finds its missing detectors (deterministic measurement
products independent of the annotated detectors/observables), computes the
detecting region of every detector/observable/missing degree of freedom, and
emits a single self-contained HTML file where you can pick a target and
multiply other detectors into it by clicking their polygons, detslice-style.

Usage:
    python explore.py circuit.stim -o out.html
    python explore.py circuit.stim --unknown-input --open
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import webbrowser

import stim


def measurement_sets(circuit: stim.Circuit) -> tuple[list[set[int]], dict[int, set[int]]]:
    """Absolute measurement indices for each detector and observable.

    Returns (detector_rec_sets, observable_rec_sets). Observables accumulate
    across multiple OBSERVABLE_INCLUDE instructions via symmetric difference.
    """
    det_recs: list[set[int]] = []
    obs_recs: dict[int, set[int]] = {}
    n_meas = 0
    for inst in circuit.flattened():
        if inst.name == "DETECTOR":
            recs = {n_meas + t.value for t in inst.targets_copy()}
            det_recs.append(recs)
        elif inst.name == "OBSERVABLE_INCLUDE":
            idx = int(inst.gate_args_copy()[0])
            recs = {n_meas + t.value for t in inst.targets_copy() if t.is_measurement_record_target}
            obs_recs[idx] = obs_recs.get(idx, set()) ^ recs
        else:
            n_meas += inst.num_measurements
    return det_recs, obs_recs


# Pauli algebra on codes 0=I 1=X 2=Y 3=Z via (x, z) bits.
_XB = {0: 0, 1: 1, 2: 1, 3: 0}
_ZB = {0: 0, 1: 0, 2: 1, 3: 1}
_P_OF = {(0, 0): 0, (1, 0): 1, (1, 1): 2, (0, 1): 3}


def region_locs(tick_regions: dict[int, stim.PauliString]) -> dict[tuple[int, int], int]:
    """Detecting region as a {(tick, qubit): pauli} dict."""
    locs = {}
    for tick, ps in tick_regions.items():
        for q in ps.pauli_indices():
            locs[(tick, q)] = ps[q]
    return locs


def xor_into(m: dict[tuple[int, int], int], c: dict[tuple[int, int], int]) -> None:
    """Multiply region c into region m, in place."""
    for loc, p in c.items():
        prev = m.get(loc, 0)
        np = _P_OF[(_XB[prev] ^ _XB[p], _ZB[prev] ^ _ZB[p])]
        if np == 0:
            m.pop(loc, None)
        else:
            m[loc] = np


def greedy_reduce(
    m: dict[tuple[int, int], int],
    candidates: list[dict[tuple[int, int], int]],
) -> tuple[dict[tuple[int, int], int], set[int]]:
    """Greedily multiply candidate regions into m to minimize m's location count.

    Each step picks the candidate with the best (most negative) change in total
    number of sensitive (tick, qubit) locations, until no candidate helps.
    Candidates must not include m itself (or anything spanning the same missing
    degree of freedom), or the reduction could cancel m to the identity.
    Returns (reduced region, set of candidate indices folded in).
    """
    loc_index: dict[tuple[int, int], list[int]] = {}
    for ci, c in enumerate(candidates):
        for loc in c:
            loc_index.setdefault(loc, []).append(ci)

    factors: set[int] = set()
    while True:
        overlapping = {ci for loc in m for ci in loc_index.get(loc, ())}
        best, best_delta = None, 0
        for ci in overlapping:
            delta = 0
            for loc, p in candidates[ci].items():
                mp = m.get(loc, 0)
                if mp == 0:
                    delta += 1      # new location appears
                elif mp == p:
                    delta -= 1      # same pauli cancels
                # different pauli: location stays occupied, delta 0
            if delta < best_delta:
                best, best_delta = ci, delta
        if best is None:
            break
        xor_into(m, candidates[best])
        factors ^= {best}
    return m, factors


NOISE_PREFIXES = ("DEPOLARIZE", "X_ERROR", "Y_ERROR", "Z_ERROR", "PAULI_CHANNEL",
                  "CORRELATED_ERROR", "ELSE_CORRELATED_ERROR", "E", "HERALDED", "II_ERROR")
ANNOTATIONS = {"DETECTOR", "OBSERVABLE_INCLUDE", "QUBIT_COORDS", "SHIFT_COORDS", "TICK", "MPAD", "I", "II"}


def ops_by_tick(circuit: stim.Circuit) -> dict[int, list]:
    """Group circuit operations by the tick slice they lead into.

    Slice t of detecting_regions reports sensitivities at the t-th TICK, i.e.
    after the ops between TICK t-1 and TICK t. Those ops get key t. Ops after
    the final TICK (e.g. final data measurements) get key num_ticks.
    """
    out: dict[int, list] = {}
    tick = 0
    for inst in circuit.flattened():
        if inst.name == "TICK":
            tick += 1
            continue
        if inst.name in ANNOTATIONS or inst.name.startswith(NOISE_PREFIXES):
            continue
        groups = []
        for grp in inst.target_groups():
            entry = []
            for t in grp:
                if t.qubit_value is None:
                    continue
                letter = "X" if t.is_x_target else "Y" if t.is_y_target else "Z" if t.is_z_target else None
                entry.append([t.qubit_value, letter])
            if entry:
                groups.append(entry)
        if groups:
            out.setdefault(tick, []).append({"g": inst.name, "t": groups})
    return out


def extract_data(circuit: stim.Circuit, *, unknown_input: bool, ignore_anticommutation_errors: bool,
                 reduce_missing: bool = True) -> dict:
    missing = circuit.missing_detectors(unknown_input=unknown_input)
    full = circuit + missing
    n_orig_dets = circuit.num_detectors

    regions = full.detecting_regions(ignore_anticommutation_errors=ignore_anticommutation_errors)
    det_recs, obs_recs = measurement_sets(full)
    det_coords = full.get_detector_coordinates()

    coords = full.get_final_qubit_coordinates()
    all_qubits = set(coords)
    for tick_regions in regions.values():
        for ps in tick_regions.values():
            all_qubits.update(ps.pauli_indices())
    # Fallback layout for qubits without coordinates: a row below everything.
    max_y = max((c[1] for c in coords.values() if len(c) >= 2), default=0)
    uncoordinated = sorted(q for q in all_qubits if q not in coords or len(coords[q]) < 2)
    qubits = {q: [float(c[0]), float(c[1])] for q, c in coords.items() if len(c) >= 2 and q in all_qubits}
    for i, q in enumerate(uncoordinated):
        qubits[q] = [float(i), max_y + 2.0]

    # Enumerate all targets explicitly: a target whose detecting region is empty
    # at every TICK (e.g. a measurement product that multiplies to identity) has
    # no detecting_regions entry but should still be listed.
    all_dem_targets = [stim.DemTarget.relative_detector_id(i) for i in range(full.num_detectors)]
    all_dem_targets += [stim.DemTarget.logical_observable_id(i) for i in range(full.num_observables)]

    entries = []
    for dem_target in all_dem_targets:
        tick_regions = regions.get(dem_target, {})
        idx = dem_target.val
        if dem_target.is_logical_observable_id():
            kind, recs, tcoords = "observable", obs_recs.get(idx, set()), None
        elif idx >= n_orig_dets:
            kind, recs, tcoords = "missing", det_recs[idx], det_coords.get(idx)
        else:
            kind, recs, tcoords = "detector", det_recs[idx], det_coords.get(idx)
        entries.append({
            "id": f"L{idx}" if kind == "observable" else f"D{idx}",
            "kind": kind, "index": idx, "coords": tcoords,
            "recs": set(recs), "locs": region_locs(tick_regions),
        })

    if reduce_missing:
        # Shrink each missing DOF's total sensitivity by greedily folding in
        # annotated detectors/observables. Never other missing DOFs (and never
        # itself): only annotated candidates, so the product stays a valid,
        # independent representative of the same missing degree of freedom.
        cand = [e for e in entries if e["kind"] != "missing"]
        cand_locs = [e["locs"] for e in cand]
        for e in entries:
            if e["kind"] != "missing":
                continue
            before_locs, before_recs = len(e["locs"]), len(e["recs"])
            e["locs"], factors = greedy_reduce(dict(e["locs"]), cand_locs)
            for fi in factors:
                e["recs"] ^= cand[fi]["recs"]
            print(f"  reduced {e['id']}: {before_locs} -> {len(e['locs'])} sensitive locations, "
                  f"{before_recs} -> {len(e['recs'])} recs ({len(factors)} factors folded in)")

    targets = []
    for e in entries:
        regs: dict[str, list] = {}
        for (tick, q), p in sorted(e["locs"].items()):
            regs.setdefault(str(tick), []).append([q, p])
        targets.append({
            "id": e["id"], "kind": e["kind"], "index": e["index"], "coords": e["coords"],
            "recs": sorted(e["recs"]), "regions": regs,
        })

    return {
        "num_ticks": full.num_ticks,
        "num_measurements": full.num_measurements,
        "num_missing": full.num_detectors - n_orig_dets,
        "qubits": {str(q): xy for q, xy in qubits.items()},
        "ops": {str(t): ops for t, ops in ops_by_tick(full).items()},
        "targets": targets,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("circuit", type=pathlib.Path, help="Path to a .stim circuit file.")
    parser.add_argument("-o", "--out", type=pathlib.Path, default=None,
                        help="Output HTML path (default: <circuit>.html).")
    parser.add_argument("--unknown-input", action="store_true",
                        help="Treat circuit inputs as unknown random states when finding missing detectors.")
    parser.add_argument("--ignore-anticommutation-errors", action="store_true",
                        help="Silently drop detecting-region components that anticommute with a reset.")
    parser.add_argument("--no-reduce", action="store_true",
                        help="Skip the greedy pass that shrinks each missing detector's total "
                             "sensitivity by folding in annotated detectors/observables.")
    parser.add_argument("--title", default=None, help="Title shown in the viewer (default: circuit filename).")
    parser.add_argument("--open", action="store_true", help="Open the generated HTML in a browser.")
    args = parser.parse_args()

    circuit = stim.Circuit(args.circuit.read_text())
    data = extract_data(
        circuit,
        unknown_input=args.unknown_input,
        ignore_anticommutation_errors=args.ignore_anticommutation_errors,
        reduce_missing=not args.no_reduce,
    )
    data["title"] = args.title or args.circuit.name

    template = (pathlib.Path(__file__).parent / "viewer_template.html").read_text()
    html = template.replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")))

    out = args.out or args.circuit.with_suffix(".html")
    out.write_text(html)
    n_det = sum(1 for t in data["targets"] if t["kind"] == "detector")
    n_obs = sum(1 for t in data["targets"] if t["kind"] == "observable")
    print(f"{args.circuit}: {n_det} detectors, {n_obs} observables, "
          f"{data['num_missing']} missing detectors, {data['num_ticks']} ticks")
    print(f"wrote {out}")
    if args.open:
        webbrowser.open(out.resolve().as_uri())


if __name__ == "__main__":
    main()
