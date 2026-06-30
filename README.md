# Meshless Boundary Collocation for Ultrasound Field Modelling

This repository solves 2D ultrasound field modelling (Helmholtz PDE) using meshless/boundary-collocation ideas inspired by the **Method of Fundamental Solutions (MFS)**. The propagating, reflecting and transmitting fields are assembled from real (Huygens) and virtual (MFS) sources and validated over an **infinite plane fluid–fluid interface** (current scope).

## Reference

If you use this code, please cite:

> Agesinaldo M. Silva, Naser Tanabi, Luiz Octavio V. Pereira, Marcelo Y. Matuda, Flávio Buiochi, Nicolas Perez, Marcos S. G. Tsuzuki; Through-transmission ultrasound detection across liquid–liquid interfaces. Proc. A 1 June 2026; 482 (2340): 20250756. https://doi.org/10.1098/rspa.2025.0756

The following open-access dataset features paired experiment–simulation records for through-transmission ultrasound at extended liquid–liquid interfaces.

> Silva Jr AM 2026 Through-transmission ultrasound detection across liquid–liquid interfaces https://doi.org/10.5281/zenodo.17155945

---

## What’s here

- **mfsolution** – discretises the scene, assembles/solves the MFS system, and writes per-sweep HDF5 results (or to an in‑memory `Store`).
- **analyse_concurrent** – post-processes the MFS results by convolving with an **apodization window** (across piston taps per surface) using multi‑process workers; appends results into HDF5 or `Store`.
- **tseries_parallel** – builds time‑series synthesis matrices from the analysed spectra; writes per‑dataset outputs.
- **config.build_geometry** – thin wrapper that hides the legacy `create_configfile(parse_config, ...)` signature. Use this for new code.

Refactored drop‑in modules (kept behaviour, improved readability):

- `analyse_concurrent.py`
- `run_tseries_parallel.py`

> You can continue to use the legacy modules; the refactored ones provide the same I/O with clearer APIs and named arguments.

---

## Installation

```bash
python -m venv .venv && source .venv/bin/activate  # (Windows: .venv\Scripts\activate)
pip install -U pip
pip install -U -r requirements.txt
```

**Tested with** Python 3.10+ and NumPy/SciPy recent releases. OpenCL/CPU is required GPU is optional but recommended for `mfsolution` (via `Compute.InitCL("GPU")`).

---

## Quickstart (CLI)

### 1) MFS only

```bash
python main.py INPUT.yaml
# or
python mfsolution.py INPUT.yaml
```

### 2) MFS + Analysis (requires output directory)

```bash
python main.py INPUT.yaml /path/to/output
# or
python analyse.py INPUT.yaml /path/to/output
```

### 3) MFS + Analysis + Time‑series (requires tseries config)

```bash
python main.py INPUT.yaml /path/to/output TSERIES.yaml
# or
python tseries.py /path/to/output TSERIES.yaml
```

---

## Configuration & Data Flow

### Geometry/Config
- **New code**: call `config.build_geometry(INPUT.yaml, /path/to/output)` which internally runs `parse_config` and `create_configfile`, returning the standard tuple:
  ```text
  T, M, S, D, R, Neltoverlambda, nRD, g
  ```
- **Legacy code**: still supported: `create_configfile(parse_config, INPUT.yaml, /path/to/output)`.

### MFS outputs
For each frequency index `elt` and sweep pair `(skr, sdr)`, MFS writes/produces an HDF5 file:
```
{dataroot}_{elt+1}_{skr}_{sdr}.h5
```
containing at least the datasets:
- `domain` (MH) – propagation in the domain
- `receiver` (MR) – pressure at the receiver grid

`analyse_concurrent` reads these back (or from `Store`), applies apodization, and appends to per‑sweep files:
```
doma_{dataroot}_{skr}_{sdr}.h5
resp_{dataroot}_{skr}_{sdr}.h5
```
with shapes:
- `doma`: (|sfr|, |D.c|, |T.c| − PPT + 1)
- `resp`: (|sfr|, |R.c| − PPT + 1, |T.c| − PPT + 1)

`tseries_parallel` consumes `doma_*`/`resp_*` plus its own `config.yaml` to build synthesis matrices and tables.

---

## Notes on the refactors

- **No algorithm changes**, only readability:
  - Named arguments at all helper callsites (no ambiguous tuples)
  - Small dataclasses (`DatasetBundle`, `Limits`, `RunConfig`, `CreateMatrixResult`) in the refactored modules
  - Clear separation of concerns: load→map→compute→save
  - Workers take only what they need (e.g., `compute_analyse(elt, apod, MH, MR)`)

You can migrate gradually: call the refactored modules from your existing scripts, or replace imports one by one.

---

## Citation

If you use this code in research, please cite this repository and standard references on MFS relevant to your work.

---

## Troubleshooting

- *OpenCL errors*: confirm your platform is detected; fall back to CPU by switching the device in `Compute.InitCL(...)` if needed.
- *Shape mismatches in analysis*: ensure the same `PPT` (piston taps per surface) is used when creating and consuming datasets.
- *NaNs in tseries*: verify that height refs are interpolated and that full profiles are exported independently (per earlier guidance).
