function out = validate_carla_huang_ds004977(dsDir, carlaDir, mnlDir)
% validate_carla_huang_ds004977() - iEEGLAB's CARLA against the code of Huang
%       et al. (2024) on their Figure 5: ds004977 sub-1 run-02 (sEEG, 4800 Hz),
%       stimulation site RMO8-RMO9. The paper's optimum is 154 of 206 channels.
%
% Usage:
%   out = validate_carla_huang_ds004977('C:\data\ds004977', 'C:\code\CARLA_JNM', 'C:\code\mnl_ieegBasics')
%
% dsDir     the ds004977 BIDS folder (sub-1 and derivatives/event_annotations)
% carlaDir  the authors' code (github.com/hharveygit/CARLA_JNM)
% mnlDir    the lab's toolbox it uses (github.com/MultimodalNeuroimagingLab/mnl_ieegBasics)
%
% (A) Reference, as the authors' applyCARLARealCCEPs.m: high-pass of the whole
%     recording (ieeg_highpass: Butterworth, filtfilt), epochs -1 to 1 s of the
%     events with status good, the good SEEG channels minus the contacts within
%     2 of the stimulated pair (getNeighborChs), channels annotated bad in a
%     trial left out (rmBadTrialsAnnots), rng('default'), then
%     CARLA(tt, data, srate, true). The MEF files are read with iEEGLAB's reader
%     (ieeglab_load_mefd); the authors' readMef3 needs the MEF library.
% (B) The same input: ieeglab_carla on the epochs of (A), rng('default').
% (C) iEEGLAB end to end: ieeglab_import, ieeglab_load, ieeglab_preprocess (0.5 Hz
%     high-pass, which is FIR in iEEGLAB; epochs -1 to 1 s of this site; bad
%     channels removed), pop_ieeglab_reref with CARLA, the sensitive cutoff and
%     2 neighbours left out, rng('default').
%
% out: n_reference and n_total_reference (A), n_same_input and whether (B) chose
% the same channels in the same order (same_input_identical, with the largest
% difference of the zMinMean profile), n_ieeglab and n_total_ieeglab (C), whether
% (C) chose the reference's channels (same_channels, jaccard), the median r of
% the re-referenced trial means of (A) and (C) over 10 to 300 ms, and the
% trial-mean common average subtracted by (C) against (A) (car_r,
% car_max_abs_diff_uv).
%
% Reference:
%   Huang, H., Ojeda Valencia, G., Gregg, N. M., Osman, G. M., Montoya, M. N.,
%   Worrell, G. A., Miller, K. J., & Hermes, D. (2024). CARLA: Adjusted common
%   average referencing for cortico-cortical evoked potential data. Journal of
%   Neuroscience Methods, 407, 110153.
%
% Cedric Cannard, iEEGLAB, 2026

t0 = tic;
addpath(fullfile(carlaDir, 'functions'), fullfile(mnlDir, 'functions'), fullfile(mnlDir, 'utilities'));
cleanup = onCleanup(@() rmpath(fullfile(carlaDir, 'functions'), fullfile(mnlDir, 'functions'), fullfile(mnlDir, 'utilities')));
sub = 'sub-1'; stem = sprintf('%s_ses-ieeg01_task-ccep_run-%02d', sub, 2);
site = 'RMO8-RMO9';
ieeg = fullfile(dsDir, sub, 'ses-ieeg01', 'ieeg');

% ---------------- (A) reference: the authors' pipeline ----------------
channels = ieeg_readtableRmHyphens(fullfile(ieeg, [stem '_channels.tsv']));
events = ieeg_readtableRmHyphens(fullfile(ieeg, [stem '_events.tsv']), 'electrical_stimulation_site', 1);
annots = readtable(fullfile(dsDir, 'derivatives', 'event_annotations', sub, [stem '_eventsAnnot.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
annots = annots.status_description;
goodEv = strcmpi(events.status, 'good');                 % ccep_PreprocessMef keeps the good events
events = events(goodEv, :); annots = annots(goodEv); %#ok<NASGU> used inside evalc below
chMef = readtable(fullfile(ieeg, [stem '_channels.tsv']), 'FileType', 'text', 'Delimiter', '\t');
[meta, dataAll] = ieeglab_load_mefd(fullfile(ieeg, [stem '_ieeg.mefd']), [], chMef.name);
srate = meta.time_series_metadata.section_2.sampling_frequency;
for c0 = 1:32:size(dataAll, 1)                            % per channel, as ieeg_highpass; in blocks for memory
    rows = c0:min(size(dataAll, 1), c0 + 31);
    dataAll(rows, :) = ieeg_highpass(dataAll(rows, :)', srate, 1)';
end
tt = (0:2*srate-1) / srate - 1;
isSite = strcmp(events.electrical_stimulation_site, site);
ev = find(isSite);
data = nan(height(channels), numel(tt), numel(ev));
for k = 1:numel(ev)
    a = events.sample_start(ev(k)) - srate;
    data(:, :, k) = dataAll(:, a + 1:a + numel(tt));
end
clear dataAll
goodChs = strcmp(channels.type, 'SEEG') & strcmp(channels.status, 'good') & ...
    ~ismember(upper(channels.name), getNeighborChs(split(site, '-'), 2));
[~, dataSite] = evalc('rmBadTrialsAnnots(data(goodChs, :, :), channels.name(goodChs), annots(ev))');
clear data
nanChs = any(isnan(dataSite), [2, 3]);
labs = channels.name(goodChs); labs = labs(~nanChs);
X = dataSite(~nanChs, :, :); clear dataSite
rng('default');
[~, Vref, CARref, Sref] = evalc('CARLA(tt, X, srate, true)');
refChs = labs(Sref.chsUsed);
out.n_reference = numel(Sref.chsUsed);
out.n_total_reference = numel(labs);

% ---------------- (B) the same input, iEEGLAB's CARLA ----------------
rng('default');
[~, ~, Sb] = ieeglab_carla(tt, X, srate, struct('sens', true, 'winResp', [0.01 0.3], ...
    'lineFreq', 60, 'notchFirst', true, 'verbose', false));
out.n_same_input = numel(Sb.chsUsed);
out.same_input_identical = isequal(Sb.chsUsed(:), Sref.chsUsed(:)) && isequal(Sb.order(:), Sref.order(:));
out.same_input_zminmean_max_diff = max(abs(Sb.zMinMean(:) - Sref.zMinMean(:)));
clear X

% ---------------- (C) iEEGLAB end to end ----------------
EEG = ieeglab_import(fullfile(ieeg, [stem '_ieeg.mefd']));
EEG = ieeglab_load(EEG, struct('event_field', 'electrical_stimulation_site'));
EEG = ieeglab_preprocess(EEG, struct('event_filters', struct('type', {{site}}), ...
    'apply_highpass', true, 'highpass', 0.5, 'apply_notch', false, 'apply_lowpass', false, ...
    'apply_ds', false, 'apply_blank', false, 'apply_epoch', true, 'epoch_window', [-1000 1000], ...
    'reject_trials', false, 'remove_bad_channels', true, 'apply_baseline', false, 'plot', false));
EEG0 = EEG;
rng('default');
EEG = pop_ieeglab_reref(EEG, 'method', 'carla', 'sensitive', true, 'neighbors', 2);
G = EEG.ieeglab.car.groups;
pool = setdiff({EEG.chanlocs.labels}, G.excluded, 'stable');
out.n_ieeglab = numel(G.reference);
out.n_total_ieeglab = numel(pool);
out.same_channels = isequal(sort(upper(G.reference(:))), sort(upper(refChs(:))));
out.jaccard = numel(intersect(upper(G.reference), upper(refChs))) / numel(union(upper(G.reference), upper(refChs)));

% re-referenced trial means on the same channels, 10 to 300 ms
[~, ia, ib] = intersect(upper(labs), upper({EEG.chanlocs.labels}), 'stable');
assert(EEG.pnts == numel(tt), 'Epoch lengths differ: %d against %d samples.', EEG.pnts, numel(tt));
tw = tt >= 0.01 & tt <= 0.3;                              % the same samples (EEG.times is rounded to ms)
mA = mean(Vref(ia, tw, :), 3); mC = mean(double(EEG.data(ib, tw, :)), 3);
r = arrayfun(@(k) local_r(mA(k, :)', mC(k, :)'), 1:numel(ia));
out.rereferenced_median_r = median(r);
out.rereferenced_min_r = min(r);

% the common average iEEGLAB subtracted: the mean of its reference channels before re-referencing
[~, iref] = ismember(upper(G.reference), upper({EEG0.chanlocs.labels}));
carC = mean(squeeze(mean(double(EEG0.data(iref, :, :)), 1)), 2);
carA = mean(CARref, 2);
out.car_r = local_r(carC, carA);
out.car_max_abs_diff_uv = max(abs(carC - carA));
out.seconds = toc(t0);

fprintf(['\nCARLA on %s (Huang et al., 2024, Figure 5; published: 154 of 206 channels).\n' ...
    '(A) the authors'' code: %d of %d channels.\n' ...
    '(B) iEEGLAB on the same input: %d; same channels in the same order: %d (largest zMinMean difference %.2g).\n' ...
    '(C) iEEGLAB end to end: %d of %d channels; the same channels as (A): %d (Jaccard %.3f).\n' ...
    '    Re-referenced trial means against (A), 10 to 300 ms: median r %.4f (smallest %.4f). Common average ' ...
    'against (A): r %.4f, largest difference %.1f uV.\nDone in %.0f s.\n'], ...
    site, out.n_reference, out.n_total_reference, out.n_same_input, out.same_input_identical, ...
    out.same_input_zminmean_max_diff, out.n_ieeglab, out.n_total_ieeglab, out.same_channels, out.jaccard, ...
    out.rereferenced_median_r, out.rereferenced_min_r, out.car_r, out.car_max_abs_diff_uv, out.seconds);
end

function r = local_r(a, b)
c = corrcoef(a, b); r = c(1, 2);
end
