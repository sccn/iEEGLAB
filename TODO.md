# iEEGLAB — remaining work

Ordered by what blocks the most. Estimates are with AI-assisted development
(see the contract document for how they were derived).

## Needs a decision first

- [ ] **EEGLAB topoplot routing (with Arno).** `ieeglab_is_ieeg` and
      `ieeglab_topoplot` exist; routing EEGLAB's own `pop_topoplot` to them needs
      a small hook in EEGLAB core, e.g. at the top of `pop_topoplot`:
      ```matlab
      if exist('ieeglab_is_ieeg','file') && ieeglab_is_ieeg(EEG)
          for lat = arg2, ieeglab_topoplot(EEG, lat); end   % one figure per latency
          return
      end
      ```
      That is Arno's call as EEGLAB maintainer. Until then the plugin adds its
      own entry to EEGLAB's Plot menu.
- [ ] **STUDY / group level.** Anatomical grouping (Destrieux labels are already
      in the tutorial electrodes.tsv) versus per-contact MNI. Decide before
      building; consider delegating to MIA.
- [ ] **EEGPrep scope with Dora** — contribute to `ieegprep`/`erdetect`, wrap
      them, or build separately.

## Engineering backlog

- [ ] Triage the ~134 static-audit findings that were never verified
      (concentrated in `ieeglab_gui_load2`, `ieeglab_load_mefd`, `get_elec_coor`).
- [ ] CI: `.github/workflows/tests.yml` is written (Ubuntu, `matlab-actions`,
      EEGLAB with submodules, runs `tests/`). Confirm the first run passes and
      add the badge to the README. Tests are compute-only because graphics hang
      in `-batch` on Windows.
- [ ] Wire the MEF3 loader (`ieeglab_load_mefd`) into the menu; declare matmef.
- [ ] CARLA as its own CCEP-only dialog tab (issue #11). The label and method
      already follow the data mode; the tab is layout work.
- [ ] Vendor `@gifti` and the four vistasoft functions actually used (issue #8).
- [ ] Shrink git history (~190 MB of NIfTI/scratch files were committed before
      `.gitignore` rules existed). Destructive: needs a force-push and a heads-up
      to anyone with a clone.
- [ ] Update the wiki: it still documents aCAR, the 1/f baseline and the old menu.
- [x] Add the real `channels.tsv` from OpenNeuro ds004696 (sub-02) to the
      tutorial (done 2026-09-21, subset to the 16 contacts; ROP3, ROP5, ROP6 are
      status 'bad'). Run the tutorial once to confirm `ieeglab_load` picks it up.

## Methods

- [ ] Basis Profile Curves (`bpc_identify.m`, Miller, Müller & Hermes) —
      clusters stimulation sites by response shape; complements CRP.
- [ ] Bipolar and Laplacian sEEG re-referencing. Not as simple as it looks for
      CCEP: a bipolar derivation such as RA2-RA3 contains the stimulated contact
      RA2 of site RA1-RA2, and a Laplacian averages neighbours that may be
      stimulated. Both must map derived channels back to their contacts so the
      stimulated-contact exclusion still applies. Straightforward for non-CCEP sEEG.
- [ ] Anatomical labelling and MNI coordinates per contact, from
      `mnl_ieegBasics` (prerequisite for any group analysis).
- [ ] Validate N1 detection against `erdetect` on shared data, as was done for
      CARLA. Partly done 2026-09-20 (Python prototype, HAPwave sub-02): latency
      within 5 ms in 92% of jointly detected pairs, amplitude r = 0.99, detection
      kappa 0.36 with the old z-rule. erdetect's rule is now an option of
      `ieeglab_detect_n1` (`polarity`, `min_baseline_sd`). Remaining: MATLAB run of
      `ieeglab_detect_n1` with erdetect-equivalent settings vs erdetect, on
      ds004080 (subject IDs are `sub-ccepAgeUMCU01`...) and on HAPwave sub-02.
- [ ] Unit tests for the new `ieeglab_detect_n1` options: polarity 'negative'
      ignores a positive peak; `min_baseline_sd` floors z on a quiet contact;
      a baseline outside the epoch is clamped with warning `baselineClamped`.
- [ ] Decide the N1 polarity default with Dora (37% of responsive pairs on
      HAPwave sub-02 have their largest early deflection positive).
- [ ] Time-frequency per contact: wavelet spectrogram and 70-170 Hz broadband,
      as in Huang et al. 2023. HAPwave sub-02: HFB rises 1.7x median over a
      distant baseline, a third of pairs above surrogate.
- [ ] Phase: use epochs of at least -2..+2 s (wavelet cone of influence), equal
      trial counts vs surrogate, and an evoked-subtracted (induced) ITPC to test
      phase reset. HAPwave sub-02: no phase reset survives ERP subtraction.
- [ ] Test blanking and N1 on native-rate data (the tutorial extract is 128 Hz,
      too coarse for either to be meaningful — both functions warn about this).
      Done in the Python prototype on HAPwave sub-02 (2048 Hz, see
      worklog/2026-09-20-native-rate-numbers.md); the MATLAB functions themselves
      still need running on it (MEF3: wire `ieeglab_load_mefd`).

## Surfaces and localisation

- [ ] `test_freesurfer.m` is macOS-only (hardcoded paths, POSIX shell) and its CT
      section reads a file written later in the script. Make it cross-platform
      and optional, or replace with FastSurfer (FreeSurfer-compatible surfaces in
      ~1 h, no FreeSurfer licence). Do not grow it into a localisation tool.
