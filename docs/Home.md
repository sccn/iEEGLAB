# iEEGLAB documentation

iEEGLAB is an EEGLAB plugin for intracranial EEG (stereo-EEG and ECoG), with a focus on cortico-cortical evoked
potentials (CCEP) from single-pulse electrical stimulation. Installation, the menus and functions, validation and
references are in the [README](https://github.com/sccn/iEEGLAB#readme).

The tutorials go step by step through the EEGLAB menus. Each step gives what you should get, the window or figure that
opens, and the same step from the command line.

| Tutorial | Data | What it covers |
|---|---|---|
| [sEEG: cortico-cortical evoked potentials](1-sEEG-CCEP-ds004696.md) | ds004696 sub-02: an extract ships with the plugin, or the full recording from OpenNeuro | Import, coordinates, brain surfaces, preprocessing, CARLA, CRP, connectivity matrix, time-frequency power, spectrum, bipolar montage |
| [ECoG: cortico-cortical evoked potentials](2-ECoG-CCEP-ds004080.md) | ds004080 sub-ccepAgeUMCU02: an extract ships with the plugin, or one run from OpenNeuro | The same on fsaverage, with the common average and the N1 |
| [FieldTrip tutorial data](3-FieldTrip-tutorial-data.md) | FieldTrip's SubjectUCI29 (205 MB download) | FieldTrip epochs, common average of the grids, bipolar depth leads, high-frequency activity tested against baseline, the comparison with FieldTrip |
| [Re-referencing compared](4-Re-referencing-compared.md) | ds004696 sub-02, full recording | Six re-referencing schemes against the published CRP of this participant |
| [CARLA on its authors' data](5-CARLA-ds004977.md) | ds004977 sub-1, run 2 (304 MB download) | CARLA with its authors' settings on the recording of their Figure 5, against the common average of every channel |

The same steps, with a checklist, are in `iEEGLAB_GUI_walkthrough.html` in the plugin folder (**Help > iEEGLAB
documentation** in EEGLAB), and as one script in
[`ieeglab_tutorial.m`](https://github.com/sccn/iEEGLAB/blob/main/ieeglab_tutorial.m).

These pages are copied from the `docs/` folder of the [repository](https://github.com/sccn/iEEGLAB/tree/main/docs);
edit them there.
