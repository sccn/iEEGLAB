# iEEGLAB

EEGLAB plugin for intracranial EEG — stereo-EEG and ECoG — with a focus on
**cortico-cortical evoked potentials (CCEP)** from single-pulse electrical
stimulation. Continuous recordings and ordinary stimulus- or response-locked
designs are supported too.

<p align="center">
  <img src="tutorial/images/6_vis_elecs3_opt.gif" width="80%">
  <br>
  <em>Electrodes on the subject's brain surfaces</em>
</p>

## What it does

| Step | Function | Menu |
|---|---|---|
| Load electrode coordinates, events and clinician channel labels from BIDS sidecars | `ieeglab_load` | Load electrode coordinates and events |
| Mark / remove bad channels (BIDS `channels.tsv` status, seizure-zone labels) | `ieeglab_bad_channels` | via Load and Preprocess |
| Blank the stimulation artifact before filtering (CCEP) | `ieeglab_blank_stim` | Preprocess |
| Filter, epoch, baseline-correct | `ieeglab_preprocess` | Preprocess |
| **CARLA** adjusted common average reference (CCEP; validated against the published code) | `ieeglab_carla`, `ieeglab_car` | Preprocess |
| **N1** detection — amplitude, latency, significance | `ieeglab_detect_n1` | CCEP analysis |
| **CRP** — canonical response shape and duration | `run_CRP` | CCEP analysis |
| **Connectivity matrix** — stimulation sites × contacts | `ieeglab_ccep_matrix` | CCEP analysis / Plot connectivity matrix |
| Electrode values on the brain (the iEEG counterpart of topoplot) | `ieeglab_topoplot` | Plot electrode values on brain |
| Export results — TSV, JSON provenance, MAT | `ieeglab_export` | Export results |

Every step runs from the menus **or** from a script with an options struct, so
whole datasets can be processed without clicking — see step 9 of
[`ieeglab_tutorial.m`](ieeglab_tutorial.m):

```matlab
EEG = pop_loadset('sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set');
EEG = ieeglab_load(EEG, struct());                        % finds the BIDS sidecars
EEG = ieeglab_preprocess(EEG, struct('remove_bad_channels',true, 'apply_highpass',true, ...
        'highpass',0.5, 'apply_epoch',true, 'epoch_window',[-500 1000], ...
        'apply_car',true, 'car_method','carla', 'apply_baseline',true));
EEG = ieeglab_stats_subject(EEG, struct('export_dir','results'));   % N1, CRP, matrix, export
```

## Requirements

- MATLAB (Statistics and Signal Processing toolboxes recommended) and
  [EEGLAB](https://github.com/sccn/eeglab) with the firfilt plugin
- [vistasoft](https://github.com/vistalab/vistasoft) for 3D brain surfaces
- 3D electrode coordinates in a BIDS `*_electrodes.tsv`
- For event-related work, events in the dataset or a BIDS `*_events.tsv`

Run `ieeglab_check_install` (or *iEEGLAB > Check installation*) to see what is
present and what each missing piece would enable.

## Tests

```matlab
results = runtests('tests');   % headless; no figures, no dialogs
```

`tests/test_carla_vs_reference.m` checks CARLA against the published
implementation (vendored unmodified in `tests/reference/`).

## Documentation

Tutorial: [`ieeglab_tutorial.m`](ieeglab_tutorial.m) and the
[wiki](https://github.com/sccn/iEEGLAB/wiki). Changes: [CHANGELOG.md](CHANGELOG.md).
Open work: [TODO.md](TODO.md).

## Citing

Please cite iEEGLAB ([CITATION.cff](CITATION.cff)) and the methods you use —
each method prints its reference when it runs. Much of the CCEP methodology
comes from Dora Hermes and the Multimodal Neuroimaging Lab
([github.com/MultimodalNeuroimagingLab](https://github.com/MultimodalNeuroimagingLab)):

- Huang H., et al. (2024). CARLA: Adjusted common average referencing for
  cortico-cortical evoked potential data. *J Neurosci Methods*, 407, 110153.
- Miller K. J., et al. (2023). Canonical Response Parameterization: Quantifying
  the structure of responses to single-pulse intracranial electrical brain
  stimulation. *PLoS Comput Biol*, 19(5), e1011105.
- Ojeda Valencia G., et al. (2023). Signatures of electrical stimulation driven
  network interactions in the human limbic system. *J Neurosci*, 43(39), 6697–6711.
- van Blooijs D., et al. (2018). Evoked directional network characteristics of
  epileptogenic tissue derived from single pulse electrical stimulation.
  *Hum Brain Mapp*, 39(11), 4611–4622.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
