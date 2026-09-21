# Native-rate numbers for the Monday meeting

HAPwave sub-02 (ds004696), full montage, 2048 Hz, 237 SEEG contacts, 695 single pulses across
68 stimulation sites. Python prototype: port of `run_CRP.m`, Picard extended-infomax ICA,
filter-Hilbert for phase and high frequency. Analysis ran on the local machine; figures in
`worklog/2026-09-20-figs/`, per-pair results in `Documents/iEEGLAB_data/results_sub02/`.

This supersedes the 19 Sep note, which used the 128 Hz tutorial extract. Two of its conclusions
changed.

## 1. CRP explained variance, with a paired surrogate null

Responsive pairs (early response |z| > 6 against the pre-stimulus baseline; 1506 of 3921 pairs
across 20 sites) have a median CRP explained variance of 0.70 per trial on the standard 15 to
1000 ms window. The matched surrogate, the same pair and trial count at random times at least 2 s
from any pulse, gives 0.45. The paired difference is +0.17, positive in 69 percent of responsive
pairs. Non-responsive pairs give 0.448 real against 0.451 surrogate, a paired difference of
-0.01, which is the sanity check: the method finds nothing where there is nothing.

The surrogate level is the number worth carrying into the argument. A rank-1 fit of
stimulation-free data already explains 45 percent of each trial, and 389 of 3921 surrogate pairs
exceed 0.9, because spontaneous synchronised events fit a single canonical shape very well. So
explained variance is uninterpretable without a matched null, and the CRP significance test at
tau_R is biased for the same reason already noted in `.claude/LESSONS.md`.

Responsive pairs: early response amplitude median -172 uV, latency median 46 ms, tau_R median
0.21 s.

Figure: `nat_fig1_crp.png`.

## 2. N1 convention, for Dora

Of 1506 responsive pairs, 554 (37 percent) have their largest early deflection positive rather
than negative. Under the current iEEGLAB convention (largest deflection of either sign) those
are reported as positive peaks; under erdetect's negative-peak convention they would be reported
as the largest negative peak instead, and the reported latency differs by more than 5 ms in 551
pairs. For the remaining two thirds the two conventions agree exactly (median latency difference
0 ms). So the choice is not cosmetic for about a third of the data, and it is worth settling
explicitly rather than by default.

## 3. Automatic bad-channel detection does not reproduce clinical labels

The `channels.tsv` marks 39 of 237 contacts bad. They are almost all the highest-numbered
contacts on each shaft (RA15, RB14-15, RO12-15, RX13-15, RZ14-15, RK11-12, RV18, RP18, ROP3,
ROP5-7), that is, contacts outside the brain, not noisy ones. Their filtered baseline SD is
*lower* than that of good contacts (62 uV against 109 uV), and their within-shaft neighbour
correlation is the same (0.92 against 0.93).

Consequently a clean_rawdata-style criterion recovers little: a robust variance threshold at
z > 5 flags 16 contacts and catches 12 of the 39 (31 percent) with 4 false alarms; at z > 3 it
catches 14 (36 percent) with 11 false alarms; a neighbour-correlation threshold below 0.2 catches
1. The answer to the question in the email is therefore that amplitude and neighbour-correlation
criteria are not a substitute for `channels.tsv` in iEEG, and any automatic detector should be
scored against clinical labels before being offered.

One useful continuity check: ROP3, ROP5 and ROP6, the three contacts flagged by variance in the
128 Hz tutorial extract yesterday, are all marked bad in the full `channels.tsv`.

Figure: `nat_fig3_badchan_phase.png`, left panel.

## 4. ICA and re-referencing, with the control that settles it

Infomax components computed before and after a common average, matched one to one by absolute
correlation of their time courses, agree at a median |r| of 0.902 across 198 channels. Two runs
of the same ICA on the same data with different random seeds agree at 0.905. The reference
therefore changes the decomposition no more than the solver's own variability does, which is the
cleanest form of the answer to Scott.

The rank still drops by exactly one (smallest to largest covariance eigenvalue 4.0e-7 as
recorded, 1.8e-16 after the common average), so PCA must be reduced to N-1. What does not
survive from the 128 Hz note is the claim that referencing inflates the apparent spread of the
maps: with 198 contacts the common-average pedestal is 1/198 rather than 1/16, and the median
number of contacts carrying more than half the peak weight is 2 as recorded against 3 after the
common average. That conclusion was an artefact of the 16-contact extract.

Figure: `nat_fig2_ica.png`.

## 5. Phase and high frequency, the gap

On 60 strongly responsive pairs: intertrial phase clustering peaks at 0.996 (4 to 12 Hz) and
0.998 (13 to 30 Hz) in the 11 to 100 ms window, against 0.38 and 0.40 for the matched surrogate,
with 97 percent of pairs above the surrogate 95th percentile. Broadband high-frequency power (70
to 170 Hz) rises by a factor of 6.3 over baseline against 1.0 for the surrogate, in 100 percent of
pairs.

Two honest readings. The high-frequency response is large, unambiguous, and entirely invisible to
iEEGLAB as it stands, which is a straightforward argument for adding it. The phase clustering,
however, saturates because a large evoked deflection produces near-perfect phase alignment by
itself: theta-alpha power also rises 5.3-fold. So these numbers do not separate phase reset from
an additive evoked response, and that separation needs the evoked component removed first. Worth
stating as an open question rather than a finding.

Figure: `nat_fig3_badchan_phase.png`, right panel.

## Caveats

One subject, one run. The CRP port follows `run_CRP.m` line by line but has not been checked
against MATLAB output, so the MATLAB cross-check comes before any number is quoted outside this
group. Picard was stopped at 100 iterations without full convergence, which is why the
seed-to-seed control matters and why it is reported alongside. ICA used 60,000 samples drawn from
the epochs of 20 sites, not the continuous recording. The early-response screen (|z| > 6) is a
simple amplitude criterion, not erdetect's detector.

## 6. N1 detection against erdetect, on the same data

erdetect 2.6.2 was run on the trial-averaged, preprocessed responses for the same 20 sites
(`erdetect.core.detection.ieeg_detect_er`, baseline epoch set to -0.45 to -0.05 s because the
default -1 to -0.1 s falls outside the extracted epoch). On the 1765 pairs common to both
analyses:

Where both flag a response, the measurements agree closely. Latency differs by less than 5 ms in
92 percent of pairs (median difference -0.5 ms) and amplitudes correlate at r = 0.989 (median
-221 uV for erdetect against -226 uV here). So the peak measurement in iEEGLAB and in erdetect
are the same quantity.

What differs is the detection threshold. erdetect flags 340 pairs, the simple |z| > 6 criterion
used here flags 664, and the two agree on presence or absence in 73 percent of pairs. erdetect's
default is deliberately conservative: 3.4 times the baseline SD with a 50 uV floor, so an
effective 170 uV threshold. That is a reason to adopt erdetect's criterion in iEEGLAB rather than
a generic z-threshold, or at least to expose it.

On the convention question, 90 pairs (21 percent of all erdetect detections) are found only when
`detect_positive=True`. That is the same question as in section 2, answered with their own tool.

Per-pair comparison: `results_sub02/erdetect_vs_mine.csv`.

## 7. Cross-check status of the CRP port

MATLAB cannot be reached from this session, so the port was checked in three other ways, all
passing: the kernel-trick PCA component equals the SVD first left singular vector to 1.5e-16; the
explained-variance formula recomputed independently matches `run_crp`'s output exactly; and a
noiseless rank-1 dataset returns explained variance 1.000000.

The MATLAB comparison itself is now a test in the repo. `tests/test_crp_vs_python.m` loads a
fixed synthetic dataset (`tests/reference/crp_test_input.csv`, `crp_test_time.csv`), runs
`run_CRP.m` on it, and compares tau_R, per-trial explained variance, the canonical shape up to
sign, and the alpha weights against `tests/reference/crp_python_reference.json`. It is
compute-only, so it is safe under `-batch`:

    results = runtests('tests/test_crp_vs_python.m')

Until that passes, every CRP number above is provisional.

## 8. Corrections after an adversarial re-check

Every number above was recomputed independently from the saved per-pair files. All reproduce
arithmetically except one, and several need reframing.

Wrong: the "97 percent above the surrogate 95th percentile" in section 5 is 96.7 percent for 4 to
12 Hz and 95.0 percent for 13 to 30 Hz.

Dropped: the phase result. Pre-stimulus ITPC is already 0.59 (4 to 12 Hz) against 0.26 for the
surrogate, which is close to the analytic small-sample bias for 11 trials, so the pre-stimulus
window is contaminated and the post-stimulus value measures the same leakage. ITPC of 0.997 in a
70 to 170 Hz broadband is what one stereotyped evoked transient produces, not oscillatory phase
alignment, and with 11 trials the measure is at ceiling. The 60 pairs used were also the top
responders. The high-frequency power result (6.3-fold) stands; the phase result does not, and
needs the evoked component removed and the filter leakage fixed before it is shown to anyone.

Reframed, section 1: selecting responsive pairs on the real data and then comparing against the
surrogate conditions on an upward fluctuation. The unconditional number across all 3921 pairs is a
median paired delta of +0.039. A split-half control (responsive pairs chosen on odd trials,
explained variance evaluated on even trials, 8 sites, 416 responsive pairs) gives +0.172, so the
conditional effect is not an artefact of trial selection, but the conditional and unconditional
numbers should be quoted together. Note also that 0.70 minus 0.45 is a difference of medians
(0.25) while +0.17 is the median of the paired differences.

Reframed, section 3: the clinician-bad contacts are bimodal, 25 with amplitude below the median of
good contacts and 14 above it. The variance criterion catches 12 of the 14 high-amplitude ones and
0 of the 25 low-amplitude ones, so the headline sensitivity of 31 percent hides a detector that is
structurally blind to flat contacts (the robust z is a monotone transform of the SD itself). The
neighbour-correlation comparison is a null (0.92 against 0.93, Mann-Whitney p = 0.87) and should
be stated as such.

Reframed, section 4: the medians hide a bimodal distribution, with about half the components
matched above 0.9 and a quarter below 0.8. The defensible statement is that the two distributions
are near-identical at every quantile, so re-referencing perturbs the decomposition no more than
the random seed does.

Reframed, section 6: 73 percent presence/absence agreement with erdetect is inflated by the shared
negative base rate (kappa 0.356). The informative asymmetry is that erdetect confirms 39 percent of
my 664 detections while I confirm 77 percent of its 340. The latency and amplitude agreement is
near-tautological, since both read the same extremum of the same averaged waveform: the median
offset is exactly one sample and 71 percent of pairs agree within one sample. It validates the
measurement, not the detection step. The "21 percent" uses 430 (negative or positive detections)
as denominator; against the 340 negative detections it is 26 percent.

Also redundant: in section 2, the 551 pairs whose latency differs between conventions are 551 of
the 554 pairs whose largest deflection is positive, and 0 of the 952 negative ones. It is one fact,
not two.

## 9. Phase, redone properly, and the corrected high-frequency number

The earlier phase result failed because of the cone of influence: in a 0.5 s pre-stimulus window,
a 4 Hz Morlet wavelet with 5 cycles spans 0.625 s, so the "baseline" ITPC was already measuring
the evoked response. Re-extracted epochs of -2 to +2 s for 48 pairs over 6 stimulation sites,
baseline taken at -1.8 to -1.0 s, trial counts equalised between real and surrogate (10 each),
0.5 Hz high-pass, Morlet wavelets 4 to 120 Hz.

Pre-stimulus ITPC is now 0.26 to 0.29 across all frequencies, which is the analytic small-sample
level for 10 trials, and the surrogate sits at the same value. That is the check that the earlier
contamination is gone.

After the pulse, ITPC reaches 0.95 at 4 to 24 Hz, falling to 0.37 by 120 Hz, against 0.20 to 0.29
for the surrogate. Subtracting the trial-averaged response first, post-stimulus ITPC drops to 0.14
to 0.20, its own pre-stimulus level is 0.13 to 0.17, and the surrogate's is 0.21 to 0.38. In other
words, nothing survives removal of the evoked response: in this subject the phase locking in CCEPs
is what an additive evoked response produces, with no measurable phase reset. That is a clean
negative result and it is reportable.

The high-frequency number also changed once the baseline moved away from the pulse and trial
counts were equalised. Broadband 70 to 170 Hz power over the 11 to 100 ms window, referenced to
-1.8 to -1.0 s, rises by a median factor of 1.74 (IQR 1.15 to 5.17) against 1.12 for the surrogate,
and 33 percent of the 72 pairs exceed the surrogate 95th percentile. The earlier "6.3-fold in every
pair" used a baseline immediately before the pulse on the 60 strongest responders and was
inflated. The effect is real but smaller than first reported.

Files: `results_sub02/itpc_fixed.csv`, `results_sub02/hfb_long.csv`.

## 10. Changes made to the plugin

`ieeglab_detect_n1.m` now carries erdetect's rule as options rather than as a separate tool:
`.polarity` ('abs' default, 'negative' as erdetect, 'positive'), and `.min_baseline_sd` (default
50 uV) as a floor on the baseline SD, so the default threshold of 3.4 gives erdetect's effective
170 uV criterion. A baseline window that does not fit the epoch is now clamped to the available
pre-stimulus period with a warning, instead of the silent zero-detection failure erdetect has with
its default of -1000 to -100 ms. erdetect-equivalent settings are documented in the header.

`tutorial/dataset_seeg/sub-02_ses-ieeg01_task-ccep_run-01_channels.tsv` was added, built from the
real HAPwave file and subset to the 16 contacts in the extract. ROP3, ROP5 and ROP6 now carry
status 'bad', so the tutorial exercises the intended path: `ieeglab_bad_channels` already reads
channels.tsv as its first source, the file was simply missing. Both changes need
`runtests('tests')` on a machine with MATLAB before they are trusted.
