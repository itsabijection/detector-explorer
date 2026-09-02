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


def extract_data(circuit: stim.Circuit, *, unknown_input: bool, ignore_anticommutation_errors: bool) -> dict:
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

    targets = []
    for dem_target, tick_regions in sorted(
        regions.items(), key=lambda kv: (kv[0].is_logical_observable_id(), kv[0].val)
    ):
        idx = dem_target.val
        if dem_target.is_logical_observable_id():
            kind, recs, tcoords = "observable", obs_recs.get(idx, set()), None
        elif idx >= n_orig_dets:
            kind, recs, tcoords = "missing", det_recs[idx], det_coords.get(idx)
        else:
            kind, recs, tcoords = "detector", det_recs[idx], det_coords.get(idx)
        targets.append({
            "id": f"L{idx}" if kind == "observable" else f"D{idx}",
            "kind": kind,
            "index": idx,
            "coords": tcoords,
            "recs": sorted(recs),
            "regions": {
                str(tick): [[q, ps[q]] for q in ps.pauli_indices()]
                for tick, ps in tick_regions.items()
            },
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
    parser.add_argument("--title", default=None, help="Title shown in the viewer (default: circuit filename).")
    parser.add_argument("--open", action="store_true", help="Open the generated HTML in a browser.")
    args = parser.parse_args()

    circuit = stim.Circuit(args.circuit.read_text())
    data = extract_data(
        circuit,
        unknown_input=args.unknown_input,
        ignore_anticommutation_errors=args.ignore_anticommutation_errors,
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
