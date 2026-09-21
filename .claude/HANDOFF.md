# Handoff, 2026-09-21

Sessions of 19-21 Sep 2026 were run from Claude (Cowork) on the phone, with the repo and data
reached through the desktop app. Committed and merged into main on 21 Sep 2026 (see
worklog/2026-09-21-merge-meeting-doc.md); the meeting document is outside the repo.

## Why this work happened

Preparing a Monday 21 Sep meeting with Dora Hermes (CCEP methods), Arno Delorme (EEGLAB,
contract scope) and Scott Makeig. Scott pushed on four points by email: ERP vs whole-signal
variance, re-referencing and ICA, traveling waves, and multi-model AMICA / RELICA for state
non-stationarity. The analyses answer those and the open questions in Cedric's earlier email to
the group (N1 convention, native-rate blanking/N1, automatic bad-channel detection, validation).

## Final results (HAPwave ds004696 sub-02, 2048 Hz, Python prototype)

Full write-up with corrections: `worklog/2026-09-20-native-rate-numbers.md`. Sections 8-9 override
sections 1, 3, 4, 5 and 6 where they differ. In short:

- N1 vs erdetect 2.6.2: where both detect, latency within 5 ms in 92% (71% within one sample),
  amplitude r = 0.99; detection agreement kappa 0.36 (erdetect 340 pairs, old z-rule 664). 37% of
  responsive pairs have their largest early deflection positive.
- Bad channels: channels.tsv marks 39/237 bad, bimodal (25 quiet, 14 loud). Variance z>5 catches
  12/14 loud and 0/25 quiet; neighbour correlation shows no difference (p = 0.87).
- CRP: responsive pairs 0.70 explained variance vs 0.45 matched surrogate; paired delta +0.17
  (split-half control +0.17); all-pairs delta +0.04. Provisional until the MATLAB cross-check passes.
- ICA: re-referencing changes the decomposition no more than a different random seed (quantiles
  match). Rank drops by one. Removed from the email at Cedric's request (he already made the point).
- Phase: with +-2 s epochs pre-stimulus ITPC sits at chance (0.26-0.29 for n=10). Post-stimulus
  ITPC 0.95, but after subtracting the ERP it falls to baseline: no phase reset in this subject.
- HFB 70-170 Hz: median 1.74x over a distant baseline vs 1.12 surrogate; 24/72 pairs above null.

Figures for the email: `worklog/2026-09-20-figs/email_figB_badchannels.png`,
`email_figC_phase.png`. Email draft: `worklog/2026-09-21-email-draft.md`.

## Code and data changes (uncommitted)

- `functions/ieeglab_detect_n1.m`: new options `polarity` ('abs' | 'negative' | 'positive') and
  `min_baseline_sd` (default 50 uV); baseline outside the epoch is clamped with warning
  `ieeglab_detect_n1:baselineClamped` instead of erroring. Header documents erdetect-equivalent
  settings. NOT yet run in MATLAB.
- `tutorial/dataset_seeg/sub-02_ses-ieeg01_task-ccep_run-01_channels.tsv`: new, from ds004696,
  16 contacts, sampling_frequency set to 128; ROP3, ROP5, ROP6 status 'bad'.
- `tests/test_crp_vs_python.m` + `tests/reference/crp_test_input.csv`, `crp_test_time.csv`,
  `crp_python_reference.json`: run_CRP vs the Python port on fixed synthetic data.
- `TODO.md`, `CHANGELOG.md` (Unreleased), `.claude/LESSONS.md`, `.claude/CLAUDE.md` updated.
- `worklog/2026-09-19-crp-ica-numbers.md` marked superseded.

## Scripts (to reproduce)

- `worklog/2026-09-20-scripts/local_vm/`: extraction from MEF3 (`extract3.py`, `extract_surr.py`,
  `extract_long.py`), `analysis.py` (CRP port, preprocessing, N1 helper), `ica_run.py`, and
  `run_analyses.py` (the interactive steps collected as functions, with the results they gave).
  Needs `pymef numpy scipy pandas python-picard mne erdetect`. They write to `$HOME` and to
  `Documents/iEEGLAB_data/results_sub02/`.
- `worklog/2026-09-20-scripts/cloud/`: the 128 Hz tutorial analyses and figure scripts.

## Data locations (Windows)

- HAPwave full dataset: `C:\Users\ccann\Downloads\v1.0.0` (20 GB, unzipped NEMAR ZIP).
- `C:\Users\ccann\ds004696-download`: deno openneuro clone, metadata only (annex not fetched).
- `C:\Users\ccann\ds004080`: root metadata only. Subject IDs are `sub-ccepAgeUMCU01`..., so
  `openneuro-py download --dataset=ds004080 --include=sub-ccepAgeUMCU01 --target-dir=ds004080`.
- ds004977 (CARLA) and ds004457 (Huang 2023): not downloaded. Check `participants.tsv` for IDs.
- Results: `C:\Users\ccann\Documents\iEEGLAB_data\results_sub02\` (CSV/NPZ per analysis).
- Python on Windows: `%LOCALAPPDATA%\Programs\Python\Python312` (not on PATH by default).
  PowerShell 5.1 needs `[Net.ServicePointManager]::SecurityProtocol = 'Tls12'` for S3.

## Next steps, in order

1. Done 21 Sep: full suite 80/80 (`test_crp_vs_python` passes; four tests updated for the
   tutorial channels.tsv, which `ieeglab_load` now picks up, marking ROP3/5/6 bad). Still to do:
   the unit tests for the new N1 options listed in TODO.md.
2. Run the tutorial end to end once in the GUI (dialogs were never clicked through after the audit).
3. N1 validation in MATLAB against erdetect, with erdetect-equivalent settings, on ds004080
   (download one subject first) and HAPwave sub-02 (needs the MEF3 loader in the menu path).
4. After the meeting: record Dora's answers (N1 polarity default; how often channels.tsv lacks
   status), Arno's scope decisions, and Scott's toolbox answer (Muller wave-matlab or WaveSpace).
5. Time-frequency module design (wavelet spectrogram, HFB), then the phase protocol above.

## Open questions sent to the group

N1 polarity convention; how often BIDS iEEG arrives without status annotations; dataset choice
for validation (ds004696, ds004080, ds004977); MAPA (Tang, Spalding, Cogan, arXiv:2609.13507)
as a downstream consumer or as an idea for cross-subject contact comparability in group analysis.
