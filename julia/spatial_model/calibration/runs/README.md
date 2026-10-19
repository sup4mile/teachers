# Saved calibration evidence

This directory retains selected outputs cited in the [calibration note](../spatial_calibration.md)
and [calibration log](../calibration_log.md). The local `.gitignore` excludes other
run artifacts by default; add explicit exceptions when retaining new evidence.

## Preservation for the accepted next phase (2026-10-15)

**Update (2026-10-16).** The user provisionally accepts the next measurement and
mechanism diagnostics, conditional on preserving previous work, and prefers
simplicity. College choice and licensure are optional experimental benchmarks;
new comparisons should include simpler alternatives close to the earlier model.
The [current preferences and sequence](../calibration_log.md#accepted-preferences-and-next-sequence-2026-10-16)
qualify the earlier plan below. Historical candidate designations are unchanged.

[`_code/accepted-preferences-20261016/`](_code/accepted-preferences-20261016/README.md)
holds the three pre-edit documentation files in a hash-verified local archive.
It supplements the earlier source/evidence checkpoint; it is not a new full
source or run backup. This documentation update changes no source or fit,
launches no jobs, and leaves B1 sealed.

The user accepted the [development plan](../calibration_log.md#accepted-development-plan-2026-10-15)
conditional on preserving previous work. Existing runs, including unsuccessful
experiments, preregistrations, working logs and code snapshots must remain intact.
New experiments use new named output directories; do not reuse historical paths.
The September exact fit remains the reference, C4′ the designated development
candidate and S1·R·E a diagnostic benchmark.

[`_code/accepted-direction-20261015/`](_code/accepted-direction-20261015/README.md)
contains a verified local archive of 505 files: current source/configuration,
pre-update spatial documentation and selected fit evidence, including uncommitted
source work. Its manifest records SHA-256 hashes; Git HEAD and status are recorded
separately. The archive is local and Git-ignored, not an off-site backup or a full
optimizer/raw-data checkpoint. Existing empirical inputs and other run artifacts
remain in their original locations. Sealed data were not read or archived; B1
still requires explicit user instruction before opening.

## Evidence inventory

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
- `exp-c3/c4p/`: C4′, the development candidate until round 2026-10-16 (`report.txt`, `theta.toml`, `moments.toml`), and
  `exp-1015/reeval/c4p/`, its re-solve on the 2026-10-15 code.
- `exp-1015b/s1_r_e_fit/`: S1·R·E's exact fit, the diagnostic benchmark.
- `exp-1015/audit/README.md`: the teacher-supply audit table.

## Local experiment runs (not tracked)

The October 2026 experiment rounds stay on disk for reference but are not tracked; their
results are summarized in the [calibration log](../calibration_log.md). Only the
candidates listed above are tracked.

- `exp-final/*/`: the October 2026 calibration experiments, every case re-solved at its
  estimate on the final code (`../slurm/cases_final.tsv`; `moments.toml` holds every model
  moment and the case's targets). The prefix names the experiment: `p0_` baseline, `p1_`
  revenue, `p2_` target tradeoffs, `p3_` focused panel, `p4_` teacher pay, `p5_`/`pc_`/`pr3_`
  sorting and combined specifications (`pr3_teacher_sel` is the teacher block with ability
  selection), `p6_` hold-out half-A refits.
- `panel-exact/*/`: the focused sensitivity panel's Levenberg–Marquardt runs.
- `exp-counterfactual/*/`: validation outcomes with teacher composition held fixed
  (`spatial_counterfactual.jl`).
- `exp-roy/*/` and `exp-peer/*/`: the 2026-10-07 experiments, the regular Roy model (γ = 1) and
  local peer effects (ψ_a), as Levenberg–Marquardt refits (`report.txt`, `theta.toml`, `moments.toml`);
  case lists and the working log are in `experiments-2026-10c/` (untracked).
- `exp-1015/`: the pre-registered 2026-10-15 round (`report.txt`, `theta.toml`, `moments.toml` per case): the
  supply repair `s1` and its combinations with the utility moving cost (`_u`), the BA suburban amenity (`_e`) and
  the staffing closure (`_sc`); `u_c4p`, `e_c4p`, `u_e_c4p`, `theta_c4p` on C4′; the θν sensitivity `c4p_tnu*`.
  `reeval/` re-solves fits on the final code (with the own-ability mobility moments), `explore/` holds the
  persistence exploration (ρz raised at fixed ϑ), and `audit/` the teacher-supply audit
  (`spatial_supply_audit.jl`; table in `audit/README.md`). `exp-1015b/` is the follow-up round with the
  non-pecuniary staffing closure (b_T):
  - `s1_r` is S1·R.
  - `s1_r_e_fit` is S1·R·E's exact fit, with its Nz = 9 re-solve `_nz9` and identification Jacobian `_jac`;
    `s1_r_e_root` re-fits it from the same start.
  - The fit was found through the ζ_BA profile `s1_r_e_pz*` (diagnostics with ζ_BA fixed; `_nz9` at Nz = 9).
  - `s1_r_e_pzroot_ms*` is the multistart at the root.
  - `c4p_nz9` re-solves C4′ at Nz = 9.

  Pre-registrations, case lists, scores and working logs are in `experiments-2026-10d/` and
  `experiments-2026-10e/` (untracked); frozen code in `_code/` (untracked).
- `exp-1015c/`: round 2026-10-15c, income-graded sorting (Y, `ζ_y`) and teacher demand (D, `fte_rule`):
  - the GE fits at Nz = 9 (`s1_d0`, `s1_y`, `s1_y_d0`, `s1_y_r`, `s1_y_d00`) and their Nz = 13 re-solves
    (`*_nz13`);
  - `*_c2eval`: the r1015c1 fits of S1·D0 and S1·Y·D0 re-evaluated with r1015c2's closure update;
  - `bitwise_pre/`, `bitwise_dev/` and `tests/`, `tests_c2/`, `tests_c3/`: the switch-off checks and the closure tests
    for each code snapshot;
  - `bounded/`: the bounded test at S1·R's aggregates (`{y,z,yz}_root`).

  No candidate met every criterion, so C4′ stays. S1·Y·D0 is the best without a severe failure.

  The round's files are in `experiments-2026-10f/` (untracked): the pre-registration, case lists, scripts
  (`spatial_bounded.jl`, `score_1015c.py`), scores (`scores_final.md`), working log and the race-composition
  diagnostic (`race_concept/`).

  Code:
  - The fits ran on `_code/r1015c1`.
  - S1·Y·D00's fit, the Nz = 13 re-solves and the `*_c2eval` re-evaluations ran on `r1015c2`, which changes only the
    closure update.
  - `tests_c3/`, the C4′ closure check and the smoke test ran on `r1015c3` (r1015c2 plus review fixes). The main
    tree equals `r1015c3`.
- `exp-1016/`: round 2026-10-16, which covers measurement concepts, the teacher-quality route and the simpler
  (no-college) comparison:
  - fits at Nz = 9:
    - `p_d0`, `p_y_d0`, `p_y_d0_rp04`, `p_tb_y_d0`, `m_y_d0_rp0`;
    - the CFR-ρ variants `m_y_d0_rho`, `m_y_r_rho` and `p_tb_y_d0_rho` (the selected development candidate);
    - their Nz = 13 re-solves (`*_nz13`);
  - `tx_m_y_d0`, `tx_p_tb_y_d0`: diagnostic refits with teachers exempt from the income taste (never selectable),
    with their Nz = 13 re-solves (`*_nz13`, run 2026-10-17);
  - `t/`: the fixed-ϑ route diagnostics (`*_Ta`, `*_Tb`);
  - `jac/`: the identification check of the selection and of its ρ = 1 version;
  - `ms/`: the 8-start multistart of the selection;
  - `reeval/`: the ten inventory fits and the five S fits re-evaluated on r1016a, which adds the teacher-sorting moments
    (the ρ and TX fits ran on r1016a directly). It also holds `s1`, the clean S1 fit of exp-1015 used for the S1
    removal row;
  - `reeval_pre/`, `reeval_pre_pair/`: the same on r1015c3, alone and in same-node pairs with r1016a (the bitwise check).

  The selection is P·TB·Y·D0·ρ, the first candidate to meet every criterion (under the round's pre-registered A4 and
  A5 rules for fits without college and for overidentified fits), and it replaces C4′ as the development candidate. The September fit stays the reference.
  The new development candidate's files (`exp-1016/p_tb_y_d0_rho/report.txt`, `theta.toml`, `moments.toml`) are local
  and untracked. Tracking them needs a `.gitignore` exception like C4′'s, which is the user's call.

  The round's files are in `experiments-2026-10g/` (untracked):
  - the pre-registration and its amendments;
  - case lists and the working log;
  - measurement scripts and outputs (`measure/`);
  - the scorer (`score_1016.py`, `score_final.sh`) and the final scores (`scores_1016.md`, plus the
    partition-share sensitivity);
  - the route rule and multistart scripts.

  Code: the five S fits (`p_d0`, `p_y_d0`, `p_y_d0_rp04`, `p_tb_y_d0`, `m_y_d0_rp0`), their `*_nz13` re-solves,
  `reeval_pre/` and `reeval_pre_pair/` ran on `_code/r1015c3`. Everything else, including the selected `p_tb_y_d0_rho`
  and both `tx_*` refits, ran on `_code/r1016a`, which the main tree equals.
- The other experiment directories (`exp-teacher/`, `exp-round3/`, …) hold the original
  Levenberg–Marquardt runs; `exp-final/` reproduces their estimates.

Only `report.txt` and `theta.toml` (and, for the 2026-10 candidates, `moments.toml`) are
tracked from the diagnostic runs.
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
