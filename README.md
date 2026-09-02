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

By default each missing degree of freedom is preprocessed to minimize the total
number of sensitive (tick, qubit) locations in its detecting region, so the
viewer starts from a compact product instead of a sprawling
`missing_detectors()` one. The search is randomized greedy with restarts:
repeatedly fold in the factor with the most negative location-count delta
(ties broken randomly; later restarts sometimes explore the second-best tier),
and when no single factor helps, escape the plateau with *chain moves* — start
from a low-delta seed, then keep folding candidates that share a location with
the previous factor until the cumulative delta turns negative (round-to-round
detector chains often need 100+ links before paying off, so chains are capped
only by the candidate count). Factors may be annotated detectors, observables,
and the *other* missing DOFs — each step is an elementary row operation on the
missing subspace, so the reduced set stays an independent basis; only
multiplying a DOF into itself is forbidden (it would cancel to identity).
Verified after reduction: the products remain deterministic and appending them
leaves no remaining missing detectors. No optimality guarantee; disable with
`--no-reduce` to see the raw products.

Options:

- `--no-reduce` — skip the sensitivity-reduction pass.
- `--restarts N` — randomized-greedy restarts per missing detector (default 8).
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
  paste `DETECTOR rec[...]` line with copy and **save** buttons.
- **Save** stores the current product under *Rewritten missing DOFs* in the
  sidebar (labeled `R0`, `R1`, … with their composition, deletable with ✕).
  Rewritten DOFs are full targets: click to select, `⊗` to multiply in.
  The group's **check independence** button verifies, over GF(2) modulo the
  annotated detectors/observables, that each rewrite represents a genuinely
  distinct missing DOF (flagging duplicates and rewrites that collapsed into
  the annotated span), shows each rewrite's decomposition over the original
  missing DOFs (a rewrite can represent a *product* like `D22⊗D23`, which is
  then covered as a whole even though neither factor is individually), and
  reports which original missing DOFs are covered — all green means the saved
  set is a valid replacement basis ready to paste into the circuit.
- **Paste a detector** as raw measurement records (e.g.
  `DETECTOR rec[-1224] rec[-669]`; negative or absolute indices) into the box
  under the export line. Its detecting region is derived client-side: the
  rec-vector is decomposed over span(annotated detectors ∪ observables ∪
  missing DOFs) by Gaussian elimination over GF(2) — that span is exactly the
  deterministic measurement products, so anything that doesn't reduce to zero
  is rejected as not-a-detector — and, since Pauli webs compose linearly, the
  region is the XOR of the decomposition terms' regions. Accepted pastes are
  saved under *Rewritten missing DOFs* with their decomposition reported.
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
