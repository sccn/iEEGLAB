# Changelog

## Unreleased

### New

- **Bipolar re-referencing** (*iEEG re-referencing > Bipolar*,
  `pop_ieeglab_reref(EEG, 'method', 'bipolar')`, `ieeglab_montage`,
  `ieeglab_apply_montage`). Each contact minus the next one on its lead
  (`RA1-RA2`), with the same pairs and weights as FieldTrip's bipolar montage
  per lead on the tutorial sEEG. A'1 is on lead A', not A (`ieeglab_lead_number`).
  Bad contacts are not paired, a missing contact leaves a gap, and with
  coordinates contacts farther apart than twice the lead's median spacing are
  not paired (grid row ends, `'spacing_factor'`; FieldTrip pairs them). Each new
  channel sits at the midpoint of its two contacts and lists its contacts
  (`chanlocs.contacts`), so stimulation sites resolve whatever its name. Works
  on continuous and epoched data.
- **Custom montages** (`'method', 'montage'`): a `.tsv`/`.csv` file or a table
  with columns `positive`, `negative` (one contact, several for a local average,
  or `n/a`; commas between labels that contain a space) and `name`, or a
  FieldTrip montage (`labelold`, `labelnew`, `tra`).
- **Frequency analysis** (*iEEGLAB > Frequency analysis*,
  `pop_ieeglab_spectral`, `ieeglab_spectral`, `ieeglab_plot_spectral`): the
  power spectrum of each contact (Welch, uV^2/Hz; continuous data skip
  discontinuities), and time-frequency power with Hann-tapered windows (fixed
  length or a number of cycles), normalized to a baseline (dB, relative change,
  ratio), with the power of a band (70-150 Hz by default) over time and on the
  brain. Only windows entirely inside the baseline are used for it. On CCEP data
  each contact's own stimulation epochs are left out of its average, windows
  that contain the pulse are left out (`'artifact_window'`, -5 to 15 ms), and
  the spectrum uses only segments away from the pulses (`'stim_exclude'`, -50
  to 1000 ms; on epochs the default window shrinks to the part before the
  pulse). Windows and segments that touch a NaN sample (artifact blanked with
  `'nan'`) are left out for that contact and epoch only. Band power of a
  spectrum is its integral over the band. The spectrum equals MATLAB's `pwelch`
  and time-frequency power equals FieldTrip's `mtmconvol` with a Hann taper, up
  to FieldTrip's constant scaling. From the menu it only draws: the dataset is
  not changed and no new dataset is proposed. The data are transformed once per
  block of epochs and reused for every frequency.
- **Documentation in EEGLAB's Help menu** (*Help > iEEGLAB documentation*, also
  at the bottom of the iEEGLAB menu): the README, the menu-by-menu walkthrough,
  the tutorial script, and the help of each menu's function.

### Changed

- Stimulation sites resolve through the contacts of each channel
  (`ieeglab_label_contacts`): after a bipolar montage, every channel that uses
  a stimulated contact is treated as stimulated (left out of the reference,
  not tested for its own site), and the dataset is still recognized as CCEP.
  Referential data behave as before. A contact listed under two channels
  now marks both. An event that names contacts all removed as bad (on a lead
  of the dataset) still counts as a stimulation, so the dataset stays CCEP.
- The contacts left out of the common average next to a stimulated pair
  (`'neighbors'`) are on the same lead with the prime kept: LA'3 is not a
  neighbour of LA2.
- *iEEG re-referencing* offers what the data allow: on continuous data the
  common average, ICA re-referencing, bipolar and custom montages; after
  segmentation also CARLA and the lowest-covariance subset (CCEP). On continuous
  stimulation data only the montages are offered: the common average and ICA
  wait for epochs, where each pulse's stimulated pair is known and left out.
  On continuous data the common average (ignoring NaN samples, as on epochs)
  and ICA work block by block in the data's own precision (no double copy of
  the recording), and ICA fits its
  unmixing on evenly spaced samples (`max_fit_samples`). Once a montage is
  applied, no further re-referencing is offered.
- **Menu**: the items that need epochs (average heatmap, CCEP analysis, CRP and
  matrix plots, electrode values on the brain, export) are greyed out until the
  data are segmented, and grouped after re-referencing and frequency analysis.
  The heatmap is now *Plot average of epochs, every channel (heatmap)*.
- **Brain surfaces are chosen once** (*iEEGLAB > Choose brain surfaces*, now
  before *Visualize electrodes*; `pop_ieeglab_surfaces`, with `'surf_files'`,
  `'brain', 'template'` or `'brain', 'auto'`). *Visualize electrodes*, *Plot
  electrode values on brain* and the band-power brain map never ask: they use
  the brain chosen, else the one that matches the space of the coordinates.
  A template brain chosen for a dataset is now also used by *Plot electrode
  values on brain* when surfaces were remembered for the participant.
- **Baseline after a montage**: a bipolar or custom montage applied to
  baseline-corrected epochs now applies the baseline correction again, as the
  other re-referencing methods do. The default baseline is a median, and the
  median of a difference is not the difference of the medians: bipolar
  channels of the tutorial kept baselines up to 51 uV away from 0. Montage
  then preprocessing now equals preprocessing then montage for the linear steps.
- **Line noise**: the preprocessing window proposes the power-line frequency and
  its first two harmonics (60, 120, 180 Hz; below Nyquist), still removed by
  fitted sines. The fundamental alone left a 180 Hz peak in the tutorial sEEG.
  The sines of the stimulated contacts are now fitted without the 2 s after
  each of their pulses (`stim_mask_ms`, -10 to 2000 ms): their artifact
  recovery carries tens of uV at the line harmonics for 1 to 2 s, and the fit
  put that back into the rest of the recording (tutorial sEEG, RA4 at 180 Hz:
  75 times the neighbouring frequencies before removal, 338 after; 15 now).
- **Line noise, compared** (help of `ieeglab_line_noise`, with references): on
  the epochs that analyses use (each contact's own stimulation epochs left
  out), a notch changes the tutorial trial average 10-50 ms after the pulse by
  0.76 uV (sEEG) and 0.94 uV (ECoG), the sine fit by 0.32 and 0.44 uV. The
  32 uV given for 1.2.0 included the stimulated contacts' own epochs, where
  the artifact remains. ZapLine and ZapLine-plus, whose spatial components
  include the pulses, changed it by 9 to 400 uV and left most of the 60 Hz on
  the sEEG.

## 1.2.0 (2026-10-25)

### 2026-10-21 to 23

- **Line noise is removed without a filter by default** (`ieeglab_line_noise`,
  `line_method` 'fit'). In windows of 4 s, a sine at the line frequency is fitted
  to the samples away from the pulses (-10 to 300 ms left out) and subtracted; the
  exact mains frequency (within 0.1 Hz of 60 or 50) is found once from each
  continuous stretch of the recording. A notch removes each response's
  own line-frequency content and spreads it up to 1.6 s around the pulse (a 2 Hz
  notch is a 3.3 s filter at 2048 Hz): on the tutorial sEEG it took 32 uV RMS out
  of the trial average around the pulse, the fit 3.6 uV, the same as 600 ms
  before it (line noise 1.7 uV per contact). A wider, shorter notch took twice as
  much. The notch stays available (*Method* in the preprocessing window,
  `line_method` 'notch'). The tutorial results are unchanged.
- *Plot average of all channels (heatmap)* now comes before *iEEG re-referencing*
  in the iEEGLAB menu, in the order of the walkthrough.
- **The brain drawn follows the space of the coordinates**
  (`ieeglab_brain_surfaces`). From `coordsystem.json`: participant coordinates
  (ACPC, ScanRAS, ...) are drawn on the participant's pial surfaces (next to the
  data or in `derivatives/freesurfer/sub-<label>`, GIfTI or FreeSurfer);
  fsaverage coordinates on fsaverage's own pial surfaces, MNI coordinates on them
  as an approximation (fsaverage from FreeSurfer or MNE: `$SUBJECTS_DIR`,
  `$FREESURFER_HOME/subjects`, `~/mne_data`). Without a matching brain the
  window says why and what to do: a template is no longer offered as if it fit
  participant coordinates. The contacts are checked against the surfaces chosen.
  ds004696 sub-02: own pial, 0 of 211 contacts outside; ds004080 (no MRI shared):
  fsaverage, 0 of 88, 0.0 mm. *Plot electrode values on brain* uses the same
  choice when no surfaces were chosen. FreeSurfer surface files (`lh.pial`,
  `rh.pial`) can now be read, besides GIfTI.
- **Brain figures open from a view that suits the electrodes:** sEEG from the
  front, slightly above, showing both hemispheres (a lead running from lateral
  to medial, as in the tutorial, was seen end-on from the side and looked like a
  single contact); ECoG face-on, from outside the brain, along the normal of the
  plane that best fits the contacts (they always opened on the right side, so
  left-hemisphere grids were seen through the brain). Drag to rotate. The figure
  title names the surface files in full.
- The iEEGLAB menu sits right of *Help* on the EEGLAB menu bar (it was before it).
- The average heatmap uses the parula colour map (was blue-white-red), with
  limits symmetric around 0 uV.
- **Import progress.** *Import raw recording* shows a progress window, and the
  command window reports the time elapsed and an estimate of the time left:
  every 16 channels of a MEF3 recording, every 2 minutes of recording (at most 10
  steps) for BrainVision, EDF and BDF files, which are now read in blocks of time
  (the same data as one read).
- The import window says that ECG, EMG, EOG and other channel types can be read
  too (for example for heart-brain analyses with BrainBeats): they stay in the
  dataset, marked non-iEEG, and are left out of the CCEP analysis.
- **Results repeat exactly.** The permutation tests (CRP and N1) and CARLA's
  bootstrap draw from their own random stream with a fixed seed (options `seed`
  and `car_seed`, default 0), so the same data give the same results every run;
  MATLAB's global random state is not touched.
- The preprocessing window proposes the power-line frequency alone (60 or 50 Hz,
  from the recording's `ieeg.json`); harmonics can be typed in.
- The CRP explained-variance histogram outlines the non-significant pairs, so they
  stay visible behind the significant ones.
- **CRP window 15-1000 ms by default** (was 15-500 ms), as in the published CRP
  analyses: with 500 ms, about a quarter of the tutorial's pairs ended exactly at the
  window end. A window longer than the epoch is cut at the epoch end with a
  warning, and the analysis says how many significant pairs last until the end
  of the window.
- **Plot CRP of one site and contact** (`ieeglab_plot_crp`): a new menu item
  that shows one pair as in Miller et al. (2023): the trials and their average
  with the response window and its duration tau_R, the similarity of the trials
  for each candidate duration, the response shape, and its weight and explained
  variance in every trial. The dialog lists each site's contacts, responding
  ones first.
- **CARLA as its authors run it on real CCEP data** (applyCARLARealCCEPs.m):
  sensitive cutoff on, and the 2 contacts on each side of the stimulated pair
  (same lead) left out of the reference. Both were off by default. The
  lowest-variance subset ('varsubset') now defaults to the HAPwave settings:
  20% of the channels, ranked on 15-500 ms.
- **CRP statistics closer to the published analyses.** The correction now runs
  across all tested pairs of the dataset (it ran within each site), as in the
  HAPwave figures. New options: Benjamini-Yekutieli FDR ('by', the published
  correction) and `crp_test` 't-test at tau_R' (the published test, for
  replication). The default test stays the permutation test, which accounts for
  tau_R being chosen from the same data. The results table adds the median
  coefficient of determination (`cod`) and the mean alpha prime
  (`alpha_prime`), the values the HAPwave figures report.
- **CRP on ECoG is off by default** in the CCEP analysis window (N1 is the ECoG
  measure); it can still be ticked. The connectivity matrix passes the new CRP
  options on, and its caption says when the t-test was used.
- **Earlier CCEP history lines now correct across all pairs.** A replayed
  `ieeglab_stats_subject` line from before this version uses the new default
  `correct_scope` 'subject'; add `'correct_scope', 'site'` for the earlier result.
- **N1 amplitude warning at 3000 uV**, where the lab's N1 code takes a peak for a
  saturated recording (it was 1000 uV, reached by real responses next to the
  stimulated pair). 'by' (Benjamini-Yekutieli) is also offered for N1.
- **No more "remove the single-trial ERP baseline" warning** after epoching:
  EEGLAB prints it after every epoching, and iEEGLAB removes the baseline itself.
  Other messages of the epoching step are still shown.
- **Re-referencing defaults belong to each method.** Settings remembered from an
  earlier run apply only when it used the same method, so varsubset after CARLA
  gets 15-500 ms and no neighbours; the history line always names the neighbours.
- **CRP pair figure refuses stale results:** if trials or contacts changed after
  the CCEP analysis, it asks to run the analysis again instead of mixing the two.
  Contact names are matched in any case.
- **Brain surfaces:** FreeSurfer surfaces are moved to scanner coordinates with
  the cras of their volume information (as `mris_convert --to-scanner`), to match
  BIDS electrode coordinates; the coordinate system is read from the
  coordsystem.json paired with the electrodes file that was loaded; other BIDS
  template spaces (Talairach, NIHPD, OASIS30, UNCInfant) get the template brain as
  an approximation; a *Template brain* choice is remembered for both brain
  figures; surfaces remembered for a participant label are offered again only
  within the same dataset.
- The notch frequency of the preprocessing window is taken from the recording's
  own `ieeg.json` when the folder holds several runs (it fell back to 60 Hz).
- **N1 polarity is negative by default** also when `ieeglab_detect_n1` is called
  directly (it was 'abs' there and 'negative' from the CCEP analysis).
- **Epochs -1000 to 1000 ms by default** in the preprocessing dialog (was -500 to
  900 ms), so that the CRP window fits.
- **Contacts whose average is not flat before the pulse are reported.** After
  preprocessing, a contact whose average over epochs deflects before the
  stimulation (3 times the median contact and 3 times what its noise predicts;
  its own stimulation epochs left out) and that the clinicians did not mark bad
  is listed in the command window and in `EEG.ieeglab.baseline_not_flat`. It is
  kept: the channels.tsv labels stay the reference. CRP easily takes such a
  slow, repeated signal for a long response. On the tutorial sEEG this is RA14
  (51 uV against a median of 8 uV); on the full ds004696 sub-02 recording, none.
- **Tutorial sEEG without RB13.** It was the only contact of its lead on the
  brain and looked like an error; the extract now has 17 contacts.

### 2026-10-09

- **Import raw recording** (`pop_ieeglab_import`, `ieeglab_import`): a new menu
  item that reads a recording of a BIDS dataset (MEF3, BrainVision, EDF, BDF,
  .set) with little memory: single precision, a few channels at a time, the
  brain channels only by default (channels.tsv type SEEG, ECOG or DBS), and
  optionally only the stretch around chosen stimulation sites. EEGLAB's BIDS
  import reads every channel as double and saves a copy: 16 minutes and 21 GB
  for OpenNeuro ds004696 sub-02, after which preprocessing ran out of memory on
  a 61 GB machine. `ieeglab_load` shifts the event times when only a stretch
  was read.
- **Fix: MEF3 reader wrongly reported missing.** The check used `exist`, which
  can answer "not a MEX file" for the first function asked about after the
  plugin folder is added to the path. It now looks at the file MATLAB resolves
  (`ieeglab_is_mex`).

### 2026-10-08

- **No high-pass filter by default.** The preprocessing dialog now proposes a
  power-line notch and baseline subtraction only: a zero-phase high-pass can ring
  around the stimulation artifact. The notch defaults to the recording's
  `PowerLineFrequency` (BIDS `ieeg.json`), else 60 Hz.
- **Brain surfaces are asked for once.** *Visualize electrodes* explains in the
  command window and in its dialog what to select (the participant's pial
  surfaces), offers the ones found in `derivatives/freesurfer/sub-<label>` or
  next to the dataset, and remembers the choice for the dataset and the
  participant. Every brain figure then uses them; *Choose brain surfaces*
  changes them. *Plot electrode values on brain* falls back to a template glass
  brain instead of drawing the contacts in empty space.
- **Scalp EEG channels** recorded with the iEEG (channels.tsv type `EEG` next to
  SEEG/ECOG channels) are detected at load and removed unless kept
  (`scalp_eeg`, `'remove'` by default; the dialog asks).
- **Plot average of all channels (heatmap)**: a new menu item
  (`ieeglab_plot_heatmap`), channels by time, usable at any step once the data
  are epoched; each contact's own stimulation epochs are left out.
- **Connectivity matrix.** *Plot connectivity matrix* asks what to show
  (responses, N1 amplitude or latency, CRP duration or explained variance). The
  response view has a key (significant, tested without a response, not tested),
  amplitude keeps its sign on a diverging scale, explained variance is on a
  fixed 0 to 1 scale, and the title names the detection method. The CRP summary
  figure gains an explained-variance panel.
- **Fixes.** The before/after re-referencing figure showed a stimulated contact:
  it now shows the contact with the largest response, own stimulation epochs
  excluded, as mean and SEM. The matrix rebuilt by `ieeglab_ccep_matrix` used
  another N1 polarity (`'abs'`) than *CCEP analysis* (`'negative'`). Figures
  open in their own window without menu bar or toolbar (`ieeglab_figure`).
- **Tutorial folders are flat.** `tutorial/seeg` and `tutorial/ecog` hold the
  `.set` file, its sidecar files and (sEEG) the pial surfaces, with no
  subfolders: open with *File > Load existing dataset*. The T1 MRI is no longer
  shipped. `iEEGLAB_GUI_walkthrough.html` now starts from the recordings as
  downloaded from OpenNeuro.

### 2026-09-30

- **Menu follows the dataset state.** Check installation is the first iEEGLAB
  item and is always available; the other items are greyed out until a dataset
  is loaded (re-referencing until it is epoched). Called without a dataset, the
  functions stop with "no dataset is loaded" instead of a MATLAB indexing error.
- **Clearer dialogs.** Options are greyed out unless the step they belong to
  runs, and each step carries a short description. The event selector lists only
  the columns that can filter, with readable names and the number of events per
  value. Re-referencing has one section per method and a description of the
  selected one. CCEP analysis is grouped into measures, connectivity matrix,
  statistics and output; the matrix's measure (CRP or N1) is now a choice and is
  always run. The seizure-onset-zone option shows how many contacts
  electrodes.tsv labels, and is greyed out when none. Channels bad in
  channels.tsv without a reason read "status bad in channels.tsv, no reason
  given".
- **Seizure-onset zone from a `soz` column.** electrodes.tsv files that mark the
  zone with a yes/no `soz` column (ds004080) are read like the `seizure_zone`
  column (ds004696). On ds004696 sub-02 (native rate), `exclude_soz` removes the 6
  SOZ contacts and `exclude_irritative` the irritative-zone contacts and their
  stimulation trials.
- **Downsampling off by default** in the preprocessing dialog (CCEP measures
  need the native rate).
- **Brain figures.** Sulci are drawn darker than gyri
  (`ieeglab_surface_sulc`, from the mesh alone), so the folding stays visible on
  a transparent surface with MATLAB's current graphics. The figures rotate with
  the mouse and have no toolbars.
- **Tests use the native-rate tutorial datasets** (`tutorial/seeg`,
  `tutorial/ecog`); the 128 Hz `tutorial/dataset_seeg` and `tutorial/dataset_ecog`
  are removed. The sub-02 T1 MRI moved to `tutorial/seeg/sub-02/ses-mri01/anat`.

- **MEF3 check on Apple Silicon.** `ieeglab_check_install` reported the MEF3
  plugin as installed when only its help (.m) files resolved, as on Apple Silicon
  Macs, for which the plugin ships no MEX files. It now requires the compiled
  reader and prints how to build it; `ieeglab_load_mefd` stops with the same hint.
- **Contacts next to the stimulated pair can be left out of the reference.**
  `pop_ieeglab_reref(..., 'neighbors', n)` (a row of the iEEG re-referencing
  dialog; `car_neighbors` in `ieeglab_car`) keeps the contacts of the same lead
  within n of either stimulated contact out of every site's reference, as the
  CARLA publication scripts do with n = 2. With it, CARLA on ds004977 sub-1, site
  RMO8-RMO9, selects the 154 of 206 channels of Huang et al. (2024, Figure 5).
  `EEG.ieeglab.car.groups` records, per site, the reference and excluded channels
  and the CARLA statistics.
- **Tutorial folders renamed** `tutorial/seeg` and `tutorial/ecog`.

- **Electrodes and surfaces checked for a shared space.** `ieeglab_check_coords`,
  called by the electrode plots, warns when more than 10% of contacts lie outside
  the surfaces drawn (`ieeglab_check_coords:spaceMismatch`), which usually means
  the coordinates and surfaces are in different spaces. On ds004696 sub-02 all
  211 contacts are inside the published pial surfaces (median 3.4 mm away).

- **ICA re-referencing** (`ieeglab_icaref`; method `'ica'` in iEEG
  re-referencing), after Michelmann et al. (2018): runica on the good contacts
  with as many components as the data's effective rank (estimated, or given),
  a chi-square test of how uniformly each component's mixing weights spread
  over the contacts, removal of the broad components (p > 0.2), back-projection
  of the local ones. Weights are normalized per component before the test, so
  it does not depend on the data's units (on raw microvolt weights no component
  ever tests as broad); `weights = 'raw'` keeps the literal form. For CCEP data
  the unmixing is estimated on the pre-stimulus samples. runica runs on one
  thread (about 100 times faster at these matrix sizes on a loaded machine) and
  with a fixed seed, so results are reproducible.

- **N1 follows the ECoG convention.** `ieeglab_stats_subject` runs N1 by default
  on ECoG CCEP data only (channel types), counting negative peaks with a 50 uV
  floor on the baseline SD, as erdetect does (`n1_polarity`,
  `n1_min_baseline_sd`). On sEEG the connectivity matrix is built from CRP by
  default, and requesting N1 warns (`ieeglab_stats_subject:n1NotECoG`).

- **iEEG re-referencing is its own step (issue #11).** New menu item iEEGLAB >
  iEEG re-referencing (`pop_ieeglab_reref`): method (CARLA per stimulation site
  for CCEP data; common average; the lowest-covariance subset of Ojeda Valencia
  et al., 2023) and its options in one dialog, where bipolar, Laplacian and
  ICA-based referencing will be added. It runs on preprocessed epochs; when a
  baseline was removed it is removed again afterwards, which gives
  exactly the result of re-referencing inside preprocessing (checked in
  `tests/test_ieeglab_reref.m`, CARLA and CAR). The Preprocess dialog no longer
  re-references; `ieeglab_preprocess` keeps `apply_car` for scripts.

- **Works after EEGLAB's BIDS import (File > Import data > BIDS).** Sidecars are
  found next to the raw BIDS file (from `EEG.BIDS.sourcefile` when EEG-BIDS
  records it, else by mapping its default `derivatives/<pipeline>/sub-...`
  output folder back to the raw dataset). The event column EEG-BIDS moved into
  `EEG.event.type` is followed through `EEG.BIDS.eInfo`, so the stimulation site
  is used instead of `trial_type` (which labelled every pulse the same).
- **Channel labels restored after EEGLAB's BIDS import.** EEG-BIDS assigns
  electrodes.tsv rows to channels by position; on ds004696 sub-02 this gives most
  contacts the wrong label. While every channel is present, labels are restored
  from channels.tsv (data order) with warning `ieeglab_load:bidsImportLabels`,
  and coordinates are matched by name. The fix itself is on branch
  `fix-chanlocs-by-name` of EEG-BIDS (2c03423), to be proposed upstream.
- **No vistasoft needed (issue #8).** `ieeglab_read_gifti` reads .surf.gii files
  (ASCII, base64, gzip; either array order), identical to the @gifti class on all
  tutorial surfaces. `ieeglab_check_install` now lists the MEF3, bva-io and
  EEG-BIDS plugins instead of vistasoft and matmef.
- **Issue #9:** the variance in the legacy fixed-fraction CAR no longer uses
  `var(..., 'omitnan')`, which failed on some MATLAB installations.
- **Native-rate tutorial datasets** `tutorial/seeg` (ds004696
  sub-02, 18 contacts, 3 sites, 2048 Hz) and `tutorial/ecog`
  (ds004080 sub-ccepAgeUMCU02, 20 contacts, 3 sites, 2048 Hz), each a small BIDS
  dataset, built and checked against the source by `tutorial/make_tutorial_datasets.m`.
- **Test: `tests/test_ieeglab_bids_import.m`** (GIfTI reader, sidecar lookup,
  loading after `pop_importbids`).
- **Divisive baseline correction removed.** It is not meaningful for
  time-domain iEEG (the baseline mean is near zero after high-pass filtering);
  only subtraction is offered.
- **Tutorial rewritten** for the native-rate datasets (`ieeglab_tutorial.m` and
  the README), sEEG with CRP and ECoG with N1.
- Events whose type is empty or `n/a` are dropped whichever event field is used.

### 2026-09-21

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

## 1.1.0 (2026-09)

### Scientific fixes: re-run any CCEP analysis made with 1.0

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
  crossing t = 0, too few samples, and overlap with the blanked stimulation
  window. Only subtractive baselines are offered.
- **Removing a contact no longer drops unrelated trials.** Several steps matched
  stimulation sites with a substring test, so removing `RA1` also removed every
  `RA10-RA9` trial.
- **Electrode coordinates are never assigned by row order silently.** When no
  label matched, the loader relabelled every channel by TSV row and reported a
  100% match.

### Other fixes

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

## 1.0.0 (2025-12)

Initial release.
