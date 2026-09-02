#!/usr/bin/env python3
"""Generate demo circuits with deliberately missing detectors."""

import pathlib
import random

import stim

HERE = pathlib.Path(__file__).parent


def remove_random_detectors(circuit: stim.Circuit, n: int, seed: int) -> stim.Circuit:
    """Delete n randomly chosen DETECTOR instructions (seeded for reproducibility)."""
    total = circuit.num_detectors
    doomed = set(random.Random(seed).sample(range(total), n))
    out = stim.Circuit()
    det_index = 0
    for inst in circuit.flattened():
        if inst.name == "DETECTOR":
            det_index += 1
            if det_index - 1 in doomed:
                continue
        out.append(inst)
    return out


def surface_code_with_missing() -> stim.Circuit:
    """d=3 rotated surface code memory with a few detectors deleted."""
    base = stim.Circuit.generated("surface_code:rotated_memory_z", distance=3, rounds=3)
    out = stim.Circuit()
    removed = 0
    det_index = 0
    for inst in base.flattened():
        if inst.name == "DETECTOR":
            det_index += 1
            # drop a couple of mid-circuit detectors and one final one
            if det_index in (6, 11, 19):
                removed += 1
                continue
        out.append(inst)
    assert removed == 3
    return out


def three_basis_measurements() -> stim.Circuit:
    """Repeated MZZ/MYY/MXX pairs: the missing DOF is a *product* of detectors."""
    return stim.Circuit("""
        QUBIT_COORDS(0, 0) 0
        QUBIT_COORDS(2, 0) 1
        MZZ 0 1
        MYY 0 1
        MXX 0 1
        TICK
        DEPOLARIZE1(0.001) 0 1
        MZZ 0 1
        MYY 0 1
        MXX 0 1
        DETECTOR rec[-1] rec[-4]
        DETECTOR rec[-2] rec[-5]
        DETECTOR rec[-3] rec[-6]
        TICK
    """)


def repetition_code_with_missing() -> stim.Circuit:
    base = stim.Circuit.generated("repetition_code:memory", distance=5, rounds=4)
    out = stim.Circuit()
    det_index = 0
    for inst in base.flattened():
        if inst.name == "DETECTOR":
            det_index += 1
            if det_index in (5, 9):
                continue
        out.append(inst)
    return out


def main() -> None:
    for name, circuit in [
        ("surface_d3_missing", surface_code_with_missing()),
        ("three_basis_measure", three_basis_measurements()),
        ("repetition_d5_missing", repetition_code_with_missing()),
        ("surface_d5_r12_z_missing", remove_random_detectors(
            stim.Circuit.generated("surface_code:rotated_memory_z", distance=5, rounds=12), 3, seed=5)),
        ("surface_d5_r12_x_missing", remove_random_detectors(
            stim.Circuit.generated("surface_code:rotated_memory_x", distance=5, rounds=12), 3, seed=7)),
    ]:
        path = HERE / f"{name}.stim"
        path.write_text(str(circuit))
        n_missing = circuit.missing_detectors().num_detectors
        print(f"{path.name}: {circuit.num_detectors} detectors, {n_missing} missing")


if __name__ == "__main__":
    main()
