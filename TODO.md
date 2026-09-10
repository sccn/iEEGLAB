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
- [ ] CI on GitHub Actions (`matlab-actions/setup-matlab`, clone EEGLAB, run
      `tests/`). Graphics hang in `-batch` on Windows; run CI on Linux with
      `xvfb-run` or keep tests compute-only as they are now.
- [ ] Wire the MEF3 loader (`ieeglab_load_mefd`) into the menu; declare matmef.
- [ ] CARLA as its own CCEP-only dialog tab (issue #11). The label and method
      already follow the data mode; the tab is layout work.
- [ ] Vendor `@gifti` and the four vistasoft functions actually used (issue #8).
- [ ] Shrink git history (~190 MB of NIfTI/scratch files were committed before
      `.gitignore` rules existed). Destructive: needs a force-push and a heads-up
      to anyone with a clone.
- [ ] Update the wiki: it still documents aCAR, the 1/f baseline and the old menu.
- [ ] Add the real `channels.tsv` from OpenNeuro ds004696 (sub-02) to the
      tutorial so bad-channel handling can be shown on real clinician labels.
      The tutorial currently has seizure-zone labels but no channels.tsv.

## Methods

- [ ] Basis Profile Curves (`bpc_identify.m`, Miller, Müller & Hermes) —
      clusters stimulation sites by response shape; complements CRP.
- [ ] Bipolar and Laplacian sEEG re-referencing, for non-CCEP sEEG.
- [ ] Anatomical labelling and MNI coordinates per contact, from
      `mnl_ieegBasics` (prerequisite for any group analysis).
- [ ] Validate N1 detection against `erdetect` on shared data, as was done for
      CARLA.
- [ ] Test blanking and N1 on native-rate data (the tutorial extract is 128 Hz,
      too coarse for either to be meaningful — both functions warn about this).

## Surfaces and localisation

- [ ] `test_freesurfer.m` is macOS-only (hardcoded paths, POSIX shell) and its CT
      section reads a file written later in the script. Make it cross-platform
      and optional, or replace with FastSurfer (FreeSurfer-compatible surfaces in
      ~1 h, no FreeSurfer licence). Do not grow it into a localisation tool.
