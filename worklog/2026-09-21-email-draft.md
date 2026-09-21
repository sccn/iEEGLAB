# Email to Dora, Arno and Scott, draft of 2026-09-21

Attach: `worklog/2026-09-20-figs/email_figB_badchannels.png`, `email_figC_phase.png`.

---

Hi Dora, Arno, and Scott,

Short update before we talk. I re-ran the open questions on native-rate data (HAPwave sub-02, 2048 Hz, 237 contacts, 695 pulses) and made three changes to the plugin.

**erdetect is now in iEEGLAB.** Running erdetect 2.6.2 alongside my detector on the same preprocessed averages, the measurement agrees closely: latency within one sample in 71% of pairs and within 5 ms in 92%, amplitudes correlated at 0.99. What differed was the threshold, so I have adopted yours: `ieeglab_detect_n1` now takes a polarity option (negative-peak as in erdetect, or largest deflection of either sign) and a 50 µV floor on the baseline SD, which with the 3.4 factor gives your effective 170 µV criterion. One thing you may want to fix upstream: erdetect's default baseline of -1000 to -100 ms silently returns zero detections when the epoch is shorter than that. iEEGLAB now clamps the baseline to the available pre-stimulus period and warns instead. On the convention itself, 37% of responsive pairs have their largest early deflection positive, so it does change what gets reported for about a third of the data.

**Bad channels: I would drop that item from the list.** Your `channels.tsv` marks 39 of 237 contacts bad, and they are bimodal: 25 quieter than the median good contact and 14 louder. A robust variance criterion catches 12 of the 14 loud ones and none of the 25 quiet ones, and neighbour correlation shows no difference at all (0.92 against 0.93, p = 0.87). So clinical labels encode something other than noise, and `channels.tsv` stays the first source, which is how the plugin was already written. Dora, the question that remains is how often you meet datasets with no status annotation at all, since that is the only case where a fallback would be needed. Related: the tutorial data shipped with iEEGLAB were missing their `channels.tsv`, so three contacts you mark bad were silently included. Fixed, using the real file.

**Phase and high frequency.** I had this wrong on the first pass and it is worth reporting why. With ±0.5 s epochs, the wavelet used to estimate low-frequency phase reaches across the stimulus, so my pre-stimulus baseline was already measuring the response. With ±2 s epochs the pre-stimulus value sits exactly at the chance level for 10 trials. Phase locking after the pulse then reaches 0.95, but once the trial-averaged response is subtracted it falls to the pre-stimulus level and below the surrogate. So in this subject the phase locking in CCEPs is what an additive evoked response produces, with no measurable phase reset. Broadband 70 to 170 Hz power rises by a median factor of 1.7 over a distant baseline, with a third of pairs clearly above the surrogate. That is smaller than I first estimated, but it is real and iEEGLAB computes nothing of the sort, so I would still put time-frequency ahead of Basis Profile Curves.

**Validation.** I added `tests/test_crp_vs_python.m`, which checks `run_CRP.m` against an independent implementation on fixed synthetic data, in the same spirit as the CARLA test. Proposed datasets: ds004696 for sEEG CCEP, ds004080 for ECoG and the erdetect comparison, ds004977 for CARLA. Does that match what you would pick?

**One thing to weigh, for Arno and Dora.** MAPA, a masked autoencoder pretrained on intracranial data (Tang, Spalding and Cogan, Duke; arXiv:2609.13507, Apache 2.0). What caught my attention is not the decoding results but the spatial encoding: each contact carries its brain region and its relative position in the array, which is what lets the model transfer across subjects. That is the same problem as our group-analysis item, where the options on the table were MIA or pooling by anatomical label and MNI coordinate. It is not a drop-in solution, since it is trained for decoding rather than for comparing evoked responses, needs a large unlabelled corpus and a GPU, and returns embeddings rather than statistics. So I see two honest readings: either a downstream consumer of what iEEGLAB already produces, or a source of ideas for how to make contacts comparable across patients. Does either look worth pursuing to you, or is it a distraction from the core?

All numbers are one subject and a Python prototype, so provisional until the MATLAB cross-check passes.

Thanks a lot,
Cedric
