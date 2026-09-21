# Changelog

## Unreleased (2026-09-21)

- **N1 detection gains erdetect's rule as options.** `ieeglab_detect_n1` now
  takes `polarity` ('abs' default, 'negative' as erdetect, 'positive') and
  `min_baseline_sd` (default 50 uV), a floor on the baseline SD, so threshold
  3.4 gives erdetect's effective 170 uV criterion. erdetect-equivalent settings
  are documented in the header. Not yet covered by tests.
- **A baseline that does not fit the epoch is clamped, not fatal.** It is
  reduced to the available pre-stimulus period with warning
  `ieeglab_detect_n1:baselineClamped`. erdetect with its default -1000 to
  -100 ms baseline returns zero detections silently on shorter epochs.
- **Tutorial: real `channels.tsv` added** for sub-02 (from ds004696), so the
  clinician bad-channel path is exercised; ROP3, ROP5, ROP6 are marked bad.
- **Test: `tests/test_crp_vs_python.m`** checks `run_CRP.m` against an
  independent implementation on fixed synthetic data
  (`tests/reference/crp_test_*.csv`, `crp_python_reference.json`).
- **Tests follow the tutorial's clinician labels.** Four tests assumed the tutorial had
  no `channels.tsv`; they now expect ROP3, ROP5 and ROP6 to be bad from load, and
  the matrix in-degree to be NaN for contacts never tested. Suite: 80 tests.

## 1.1.0 — 2026-09

### Scientific fixes — re-run any CCEP analysis made with 1.0

- **Re-referencing was not CARLA and leaked the stimulation artifact.** The step
  labelled "aCAR" and cited as Huang et al. (2024) implemented the older fixed-
  fraction variance CAR (Valencia et al. 2023), and its stimulated-contact
  exclusion silently never ran: on the tutorial dataset 22 of 30 checked trials
  had a stimulated contact inside the common average. CARLA is now implemented
  from the published code and **validated against it** (identical to 1e-9 on
  the deterministic path), and 0 of 14 stimulation sites leak.
- **CARLA is CCEP-only.** On non-CCEP data it falls back to a plain CAR with an
  explanation.
- **The "1/f" baseline was removed.** It performed spectral whitening, not
  baseline correction, and broke Hermitian symmetry so part of the signal was
  discarded.
- **Baseline correction is guarded** against windows outside the epoch, windows
  crossing t = 0, too few samples, unsafe divisive baselining, and overlap with
  the blanked stimulation window.
- **Removing a contact no longer drops unrelated trials.** Several steps matched
  stimulation sites with a substring test, so removing `RA1` also removed every
  `RA10-RA9` trial.
- **Electrode coordinates are never assigned by row order silently.** When no
  label matched, the loader relabelled every channel by TSV row and reported a
  100% match.

### Adversarial audit — 48 confirmed findings, all fixed

- **N1 significance ignored that the peak was searched for.** The z-score of
  the largest deflection in the window was tested as if the latency were fixed.
  p-values now come from a sign-flip permutation test whose statistic is the
  same maximum over the window (exact for few trials), FDR across contacts.
  The z >= 3.4 rule is still available as `method = 'sd'`.
- **CRP significance was biased by the choice of response duration.** tau_R
  is the argmax of the projection profile; a t-test at that duration rejected
  far above its nominal rate on noise. Replaced by a permutation null that
  repeats the selection.
- **Bad channels:** indices refer to the dataset passed in (they were applied
  after other channels had been removed, hitting the wrong contact); trials
  that stimulate a marked-bad contact are always dropped, not only when an
  unrelated option is on; `exclude_soz` removes the contacts as its label says.
- **Event selection:** `boundary` markers are never deleted (epochs spanning
  removed data were being kept); a filter given as text no longer matches
  single characters and silently empties the dataset; rare conditions are
  counted per site, so `ROP2-ROP4` and `ROP4-ROP2` are one site.
- **Re-running preprocessing** no longer re-applies steps stored by an earlier
  run (a second high-pass, a second CARLA), and pre-1.0 option names stored on
  a dataset no longer override explicit options.
- **History lines replay.** Preprocessing, CCEP analysis, the matrix, export and
  the electrode plot now return commands that reproduce the run without a
  dialog (preprocessing previously returned a constant string).
- **Blanking** works on every epoch of epoched data (it blanked only the first)
  and `blank_method = 'nan'` survives filtering (the whole recording became NaN).
- **Re-referencing:** the legacy variance-subset CAR reproduces the HAPwave code
  exactly (block-wise means); non-CCEP data form one reference group instead of
  one per condition; `EEG.ref` is set only when a reference was applied; CARLA
  ranks on finite samples instead of collapsing to two channels.
- **Connectivity matrix:** contacts never tested have NaN degree, not 0; a
  matrix that no longer matches recomputed N1/CRP results is detected, and
  export and plots refuse it.
- **Export** works on datasets without clinician annotations, checks every
  target before writing when `overwrite = false`, rejects unknown formats, and
  writes amplitude/latency matrices for significant responses only.
- **Electrode maps** show only significant responses for per-site metrics, and
  latency maps leave out each contact's own stimulation trials (its artifact
  dominated the map) and honour `'site'`.
- **Loader:** BIDS sidecars are matched by entities, not by the first file in
  the folder; events with `n/a`, negative or out-of-range onsets are dropped
  explicitly; `boundary` events are kept; coordinates with comma decimals or
  duplicate rows are handled.

### New

- CCEP connectivity matrix (stimulation sites × contacts) with in/out degree,
  and a matrix plot.
- N1 detection (amplitude, latency, FDR-corrected significance).
- Stimulation-artifact blanking before filtering (issue #10).
- Bad channels from clinician labels: BIDS `channels.tsv` status, `seizure_zone`
  in `electrodes.tsv`, an explicit list, or optional automatic detection. Marked
  at load, removable at preprocessing; marked channels stay out of the reference.
- Events marked `status = bad` in `events.tsv` are dropped at load.
- Results export: TSV (BIDS derivatives style, `n/a` for missing), JSON
  provenance, and a MAT bundle.
- Electrode values on the brain (`ieeglab_topoplot`), also offered in EEGLAB's
  Plot menu, as the intracranial counterpart to scalp topographies.
- `ieeglab_stats_subject` runs N1, CRP, the matrix and export as switchable
  stages; CRP was previously unreachable from the plugin.
- Every step can run without dialogs (`opt` argument), including loading, which
  finds BIDS sidecar files automatically.
- `ieeglab_check_install`, and a startup check that every menu item resolves.
- 50+ headless tests, including 7 that validate CARLA against the reference.

### Fixed issues

#1 dependency check · #2 mesh re-selection · #3 multi-mesh rendering ·
#4 EEGLAB window not updated after preprocessing · #10 filter ringing ·
#12 unclear dialog labels

## 1.0.0 — 2025-12

Initial release.
