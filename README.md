# Missing Detector Explorer

An interactive tool for finding *explanations* of missing degrees of freedom in
stim circuits by exploring alternative detector bases visually.

`stim.Circuit.missing_detectors()` finds products of measurements that are
deterministic but independent of the annotated detectors/observables — but the
products it returns are usually not a good explanation, just a pile of
seemingly unrelated measurements. This tool renders the detecting region of
every detector (annotated and missing) as time slices, detslice-style, and lets
you multiply annotated detectors into a missing one by clicking their polygons,
until the product becomes something you recognize.

## Usage

```bash
python3 -m venv .venv && .venv/bin/pip install stim
.venv/bin/python explore.py circuit.stim -o out.html
open out.html   # self-contained, no server needed
```

By default each missing degree of freedom is preprocessed by a greedy pass that
repeatedly multiplies in whichever annotated detector/observable most reduces
the total number of sensitive (tick, qubit) locations in its detecting region,
stopping when no single multiplication helps. This usually turns a sprawling
`missing_detectors()` product into a compact region so fewer manual clicks are
needed. Only annotated detectors and observables are used as factors — never
the missing DOFs themselves — so each reduced product remains a valid,
independent representative of its missing degree of freedom (verified: the
reduced products are still deterministic, and appending them leaves no
remaining missing detectors). It is a heuristic with no optimality guarantee;
disable it with `--no-reduce` to see the raw products.

Options:

- `--no-reduce` — skip the greedy sensitivity-reduction pass.
- `--unknown-input` — treat circuit inputs as unknown random states when
  finding missing detectors (passed to `missing_detectors`).
- `--ignore-anticommutation-errors` — silently drop detecting-region components
  that anticommute with a reset instead of raising.
- `--title NAME`, `--open`

## The viewer

- **Sidebar** lists missing degrees of freedom first, then detectors, then
  observables, each with its weight (number of measurements) and coordinates.
  Click a row to make it the **target**; click its `⊗` button to multiply it
  into the current product.
- **Slices** show every tick as a film reel. Each panel draws the qubits, the
  gates applied in that time step (toggle with the *gates* checkbox), faint
  dashed polygons for every other detector's region (clickable), and the
  current product's region filled in X/Y/Z colors, matching
  `diagram("detslice-svg")` conventions.
- **Click a faint polygon** to multiply that detector into the target. If
  several regions overlap at the click point, a menu pops up to disambiguate.
  Clicking a factor again (or its chip's ×) removes it. `⌘Z`/`ctrl+Z` undoes.
- **Product bar** shows the current factor chips, the product's weight vs. the
  target's original weight (green when you've simplified it), and a ready-to-
  paste `DETECTOR rec[...]` line with a copy button.
- **Timeline strip** shows the product's support weight per tick; click to jump
  to that tick's panel. *Only ticks with support* filters the reel down.
- **Zoom** with the slider, `+`/`-`, or `ctrl`/`⌘` + scroll. Panels drop the
  Pauli letter badges below ~170px and use colored dots instead.

The generated HTML is fully self-contained (data embedded as JSON, no network
access) and can be shared or archived.

## How it works

1. `missing = circuit.missing_detectors()` returns a circuit of `DETECTOR`
   instructions; appending it (`circuit + missing`) makes the missing degrees
   of freedom ordinary detectors with indices `>= circuit.num_detectors`.
2. `full.detecting_regions()` gives each detector/observable's Pauli
   sensitivity per qubit per tick.
3. Multiplying detectors is done client-side: XOR of (x,z) bits per qubit per
   tick for the regions, symmetric difference of measurement-record sets for
   the exported `DETECTOR` instruction.

Detecting regions are only reported at `TICK` instructions — add more `TICK`s
to a circuit for finer time resolution.

## Demos

```bash
.venv/bin/python demos/make_demos.py                                  # write .stim files
.venv/bin/python explore.py demos/surface_d3_missing.stim             # d=3 surface code, 3 detectors deleted
.venv/bin/python explore.py demos/three_basis_measure.stim --unknown-input   # missing DOF is a *product* of detectors
.venv/bin/python explore.py demos/repetition_d5_missing.stim
```
