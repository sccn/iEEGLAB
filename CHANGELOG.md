# Changelog

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
