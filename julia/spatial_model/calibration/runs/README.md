# Saved calibration evidence

This directory retains selected outputs cited in the [calibration note](../spatial_calibration.md)
and [calibration log](../calibration_log.md). The local `.gitignore` excludes other
run artifacts by default; add explicit exceptions when retaining new evidence.

- `screen-1-512/screen.jsonl` and `screen/screen.jsonl`: both Sobol screen batches.
  `screen/plots/` retains the figures and combined `screen_points.csv`.
- `exact-polish-base/`: current exactly identified baseline, including fitted
  parameters, moment diagnostics, Jacobian and first-stage consistency checks.
- `exact-grid-*/`: finer-grid checks at that baseline.
- `exact-ms/*/`: the 15 multistart fits used to check for other roots.
- `polish-noincome/` and `grid-*/`: the earlier overidentified baseline and its
  grid checks, retained as historical evidence.
- `panel/*/`: sensitivity reports and fitted parameters around that earlier
  baseline, including preliminary narrow-box cases. These are not sensitivity
  results for the current exactly identified baseline.
- `t7-2016-19/`: the historical alternative occupational-sample fit.

Only `report.txt` and `theta.toml` are retained from the diagnostic runs.
Optimization histories, locks, scheduler logs, scratch scripts and serialized
solutions remain ignored. These outputs preserve evidence; they are not a full
checkpoint archive for resuming optimization.

To regenerate the screen plots from the repository root:

```bash
python3 julia/spatial_model/calibration/screen_plots.py \
    julia/spatial_model/calibration/runs/screen-1-512/screen.jsonl \
    julia/spatial_model/calibration/runs/screen/screen.jsonl \
    --out julia/spatial_model/calibration/runs/screen/plots
```
