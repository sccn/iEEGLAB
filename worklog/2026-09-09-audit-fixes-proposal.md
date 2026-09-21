# 2026-09-09 → 09-16 · Audit, fixes, contract proposal

## Set out to do
- "Catch up on what this project is about, where it's left at", check open GitHub
  issues, full audit and smoke tests; then answer Arno's contract request
  (timeline, deliverables, milestones, cost; MATLAB tool + tutorial, 10 public
  datasets, Python port EEGPrep, STUDY, tests/CI, comparison with another tool, paper).
- Fix existing issues including CARLA; check Dora Hermes' lab code and literature.
- Fill the CCEP gap (connectivity matrix, export), bad channels from clinician
  labels, everything scriptable; adversarial audit; budget.
- 09-16: email to Dora (Arno asked to meet with her first), with extra figures.

## What changed (branch `fix/carla-and-open-issues`, pushed)
- `dde30fe` real CARLA, stim-contact leak fixed, headless API, tests.
- `f4e669a` CARLA validated vs reference; CCEP-only; 1/f baseline removed; N1 + blanking.
- `5c0fc03`, `c225145` baseline correction guards.
- `c8e7760` connectivity matrix, export, clinician bad channels, headless load, electrode-value plots.
- `f8bd98d` automatic trial rejection, tsv_row traceability, CI workflow.
- `4a215b1` all 48 adversarial-audit findings fixed; `tests/test_ieeglab_audit_regressions.m`.
- `603a816` CI posts failing tests as annotations; `4ee761f` single-precision tolerance fix.

## Measured
- Tests: 79/79 locally (`run_all2.m` runner) and on GitHub Actions run 34578199143 (Ubuntu).
- CARLA vs published code (`tests/test_carla_vs_reference.m`): identical to 1e-9 on single
  trials, ranking to 1e-10, multi-trial agreement distributional (bootstrap).
- Old CAR: stimulated contact inside the reference on 22/30 checked trials; now 0/14 sites.
- Adversarial audit: 48 confirmed, 0 refuted (2 blockers: N1 and CRP significance).
- Tutorial step 9 (128 Hz sEEG): 168 trials, CARLA; CRP 134/182 pairs; N1 74/182;
  matrix 13 sites x 16 contacts, density 0.41; 2 N1 values >1000 uV flagged as artifact.
- Noise calibration (regression tests): N1 and CRP permutation tests flag <5% of noise pairs.

## Decided, and why
- CARLA only on CCEP data (it ranks against a known stimulated pair); others fall back to CAR.
- Significance by permutation tests that repeat the peak / tau_R selection (the old
  z-test and run_CRP's `p_value_tR` ignore the selection and inflate false positives).
- Re-running preprocessing never re-applies stored step switches; history lines replay.
- Electrode localisation stays with external tools; the plugin consumes `electrodes.tsv`.
- Figures are drawn in Python from MATLAB JSON exports (MATLAB graphics hang in -batch here).
- Proposal (artifact `bc1248ae-...`): 3 milestones $51,500 at $80/h, 15% contingency,
  4 days/week; 7 optional additions at cost ($6,720). Work already done is not billed.

## Unfinished / next session
- Email to Dora (drafted 09-16); meeting with Dora and Arno to settle scope.
- Open a PR from `fix/carla-and-open-issues` to `main` once Arno/Dora agree.
- Everything in `TODO.md` (backlog, methods, surfaces, STUDY, Python arm).
- Dialogs were not clicked through after the audit fixes (graphics hang in -batch):
  manual GUI check of Preprocess, CCEP analysis and Export still needed.
