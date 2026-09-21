> **Superseded** by `2026-09-20-native-rate-numbers.md`. This note used the 128 Hz tutorial extract;
> its ICA map-spread conclusion was an artefact of 16 contacts and its CRP numbers were inflated by
> downsampling. Kept for the record only.

# Numbers for the Monday meeting

Prototype analysis, 19 Sep 2026, revised after a bad-channel check. Python port of `run_CRP.m`
plus infomax ICA, run on the tutorial sEEG extract shipped with iEEGLAB (`sub-02`, ses-ieeg01,
task-ccep).

## 0. Is the tutorial extract already re-referenced?

No. Three lines of evidence. The file declares `EEG.ref = 'common'`, the setname is `resampled`,
and `EEG.history` is empty apart from the EEGLAB version stamp, so nothing was recorded as
applied. The covariance is full rank (smallest to largest eigenvalue 1.5e-2 on the good
contacts), whereas a common average would have driven one eigenvalue to zero. And in quiet
segments away from stimulation the contacts share a large common component, mean pairwise
correlation 0.47 and median 0.62, which is what a single-contact clinical reference looks like
and is the opposite of what re-referenced data look like.

A caveat worth stating anyway: if a file ever were re-referenced upstream without recording it,
the visible symptoms would be near-zero pairwise correlation in quiet periods and a rank one
below the channel count. Both are cheap to check and are worth adding to `ieeglab_load` as a
warning, since neither BIDS nor the `.set` guarantees honest provenance.

## 1. Bad channels, and why they changed the numbers

Quiet-period variance flags three contacts: ROP5 and ROP6 at robust z of 103 and 102, ROP3 at
27.5, with standard deviations of 1485, 1473 and 422 uV against roughly 30 to 70 uV for the
other thirteen. They are not marked anywhere in the shipped files, which is exactly the gap
noted in TODO: the tutorial has seizure-zone labels but no `channels.tsv`.

Leaving them in mattered twice. They entered the common average, injecting their noise into
every contact. And in the ICA they were the three strongest components, so the first version of
the map figure was showing bad channels rather than brain. Everything below excludes them from
both the reference and the analysis, and this is a concrete argument for the automatic
bad-channel detection item in the email to the group.

## 2. What the data can and cannot support

Sixteen contacts on two shafts (RA1-10, ROP1-6), thirteen usable, 128 Hz, 1140 s, 169 single
pulses across 15 stimulation sites. Nyquist is 64 Hz, so no broadband high frequency and no
time-frequency worth showing. The stimulation artifact is still large past the blanking window:
median absolute voltage across contacts is 657 uV at 15 ms, 237 uV at 50 ms, 84 uV at 100 ms,
against a 19.5 uV baseline. The early response, and therefore N1 and the standard 15 ms CRP
window, cannot be quoted from this file. Everything below uses a 100 to 1000 ms window and drops
pairs whose 99th percentile voltage exceeds 1000 uV.

## 3. CRP explained variance, with a surrogate null

Real stimulation gives a median explained variance of 0.73 per trial (mean 0.70, n = 153 site by
contact pairs). Surrogate epochs, taken at random times at least 2 s from any pulse and matched
for site labels and trial counts, give a median of 0.56 (mean 0.56, n = 152). Forty-five percent
of real pairs exceed the surrogate 95th percentile of 0.759. Median selected response duration
tau_R is 0.48 s for real stimulation against 0.41 s for surrogates.

The surrogate number is the one that matters. A rank-1 fit of non-stimulation data already
explains 56 percent of each trial, because the first principal component of 8 to 25 short,
smooth, low-amplitude traces is far from zero. Explained variance quoted without that null is
not interpretable, which is the same failure mode as the biased t-test at tau_R noted in
`.claude/LESSONS.md`. Note also that the null rose from 0.45 to 0.56 once the noisy contacts
were removed, since what remains is smoother and therefore easier for a rank-1 model to fit.

For the argument with Scott: on responsive pairs the reproducible, phase-locked shape does
dominate, leaving roughly a quarter of each trial unexplained, but the margin over chance
structure is narrower than the raw number suggests. The same analysis on native-rate data is
what should be quoted, since downsampling to 128 Hz has already removed the high-frequency,
non-phase-locked content that the argument is about.

Figures: `fig1_crp_explained_variance.png`, `fig2_crp_example_fits.png`.

## 4. The reference changes the number

On the earlier all-channel version of this analysis, median explained variance was 0.64 as
recorded, 0.93 with a plain common average over all contacts, and 0.75 with the per-site,
per-trial average that excludes the stimulated pair. The plain common average also pushed 64 of
169 pairs over the artifact screen, because the stimulated contacts sit inside the average and
their artifact is subtracted into every channel. That is an independent, quantitative
confirmation of the CARLA fix reported to the group.

## 5. ICA and re-referencing: Scott is right, with one caveat

Rank. The ratio of smallest to largest covariance eigenvalue is 1.5e-2 as recorded and 9.0e-17
after the common average. The rank drops by exactly one, so PCA has to be reduced to N-1 before
ICA or one component is fitted to numerical noise.

Sources. Infomax components computed before and after re-referencing, matched one to one by
absolute correlation of their time courses, agree at a median |r| of 0.971, with 11 of 12 pairs
above 0.9. The decomposition is effectively invariant, and the argument that a trial-varying
reference operator should break the equivalence is not supported: with the per-site, per-trial
reference the agreement was 0.99 on the full channel set.

Maps. The participation ratio of the component maps is 1.55 as recorded and 3.82 after the
common average, a factor of 2.5, while the number of contacts carrying more than half the peak
weight is unchanged at a median of 1. Re-referencing adds a small pedestal to every contact
rather than changing which contacts carry the component. Comparing maps by correlation cannot
detect this at all, since Pearson correlation across contacts removes the across-contact mean,
which is exactly what the common average subtracts. So whether a map looks non-local depends on
the metric and on the reference, while the sources do not.

What this cannot address: thirteen contacts on two shafts cannot reproduce a non-local map on a
32 by 32 grid. That needs Scott's subject or a full montage.

Figure: `fig3_ica_maps.png`.

## 6. ICA time courses versus CRP canonical shapes

On the full channel set, the canonical shape of the best-explained contact at each stimulation
site was matched by at least one IC activation at a median |r| of 0.972, against 0.850 for the
surrogate version of the same test. The effect is real but the metric is weak: in a 100 to
1000 ms window at 128 Hz every slow waveform correlates with every other one, and 4 to 11
components typically exceeded 0.8 at a given site. Suggestive only. The real comparison needs
BPC clustering of stimulation sites against IC groupings, scored with an adjusted Rand index, on
native-rate data.

## Caveats to state before quoting any of this

The CRP port follows `run_CRP.m` line by line (cross-projection profile, tau_R at the maximum of
the mean profile, kernel-trick PCA, explained variance per trial) but has not been checked
against MATLAB output, because the local MATLAB workspace would not start. That check comes
first if any number leaves the room. AMICA itself was not run; infomax is the stand-in, and
multi-model AMICA is not available in this environment. One subject, one run, two shafts,
thirteen usable contacts.

## What needs native-rate data

N1 amplitude and latency, the standard 15 ms CRP window, any time-frequency or 70 to 170 Hz
broadband estimate, intertrial phase clustering, and phase gradients along a shaft. The container
used here cannot reach OpenNeuro, S3, Zenodo or OSF, so one subject and run needs to be fetched
onto the local machine and pointed at this session.
