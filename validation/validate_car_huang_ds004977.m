function T = validate_car_huang_ds004977(dsDir, carlaDir, mnlDir)
% validate_car_huang_ds004977() - iEEGLAB's common average reference against
%       the code of Huang et al. (2024) on their dataset, ds004977 sub-1
%       run-02 (sEEG, 4800 Hz), for every stimulation site of the run.
%
% Usage:
%   T = validate_car_huang_ds004977('C:\data\ds004977', 'C:\code\CARLA_JNM', 'C:\code\mnl_ieegBasics')
%
% dsDir     the ds004977 BIDS folder (sub-1 and derivatives/event_annotations)
% carlaDir  the authors' code (github.com/hharveygit/CARLA_JNM)
% mnlDir    the lab's toolbox it uses (github.com/MultimodalNeuroimagingLab/mnl_ieegBasics)
%
% Reference, as the authors' applyCARLARealCCEPs.m prepares each site for
% CARLA: high-pass of the whole recording (ieeg_highpass: Butterworth, 0.5 Hz
% pass band, filtfilt), epochs -1 to 1 s of the events with status good, the
% good SEEG channels minus the contacts within 2 of the stimulated pair
% (getNeighborChs), trials annotated bad removed and channels bad in any trial
% left out (rmBadTrialsAnnots); then the common average of what remains, as
% their figS3_compareCARvCARLA.m: mean over channels per trial. (Their CARLA,
% on the site of Figure 5, is compared by validate_carla_huang_ds004977.)
% iEEGLAB end to end: ieeglab_import, ieeglab_load, ieeglab_preprocess (0.5 Hz
% high-pass, which is FIR in iEEGLAB; epochs -1 to 1 s; bad channels
% removed), pop_ieeglab_reref('car') with the same 2 neighbours left out.
% iEEGLAB does not read the per-trial annotations in derivatives/: channels
% annotated bad in a trial stay in its average.
% Per site: the channels in each average, and the trial-mean common average
% (same trials): the one iEEGLAB subtracted (data before minus after
% pop_ieeglab_reref, on every re-referenced contact; their spread checks that
% it is one average) against the reference, over the whole epoch and outside
% the pulse (-10 to 20 ms). Where the channels are the same, the unfiltered
% common average of those channels is put through each high-pass (the
% authors' ieeg_highpass and iEEGLAB's FIR filter): if the FIR result matches
% iEEGLAB's, the filter is the whole difference. Both sides read the MEF files
% with iEEGLAB's reader (ieeglab_load_mefd); the authors' readMef3 needs the
% MEF library.
% The public dataset holds only the samples the paper's figures 1 to 6 need;
% the others are 0 (its README). Sites whose epochs are all 0 are marked
% censored: their channels are compared, their averages are not.
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
ieeg = fullfile(dsDir, sub, 'ses-ieeg01', 'ieeg');

% ---------------- reference: the authors' preparation ----------------
channels = ieeg_readtableRmHyphens(fullfile(ieeg, [stem '_channels.tsv']));
events = ieeg_readtableRmHyphens(fullfile(ieeg, [stem '_events.tsv']), 'electrical_stimulation_site', 1);
annots = readtable(fullfile(dsDir, 'derivatives', 'event_annotations', sub, [stem '_eventsAnnot.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
annots = annots.status_description;
goodEv = strcmpi(events.status, 'good');
events = events(goodEv, :); annots = annots(goodEv);
chMef = readtable(fullfile(ieeg, [stem '_channels.tsv']), 'FileType', 'text', 'Delimiter', '\t');
[meta, dataAll] = ieeglab_load_mefd(fullfile(ieeg, [stem '_ieeg.mefd']), [], chMef.name);
srate = meta.time_series_metadata.section_2.sampling_frequency;
% The public dataset holds only what the paper's figures 1-6 need; other samples
% are 0 ("censored", its README). Epochs whose SEEG channels are 0 throughout are
% not compared (the other channels are not censored).
nT = 2 * srate;
seeg = strcmp(chMef.type, 'SEEG');
censored = arrayfun(@(s) ~any(dataAll(seeg, s - srate + 1:s - srate + nT), 'all'), events.sample_start);
sites = unique(events.electrical_stimulation_site, 'stable');
goodSite = cell(numel(sites), 1); carRaw = cell(numel(sites), 1);
for s = 1:numel(sites)                                    % unfiltered common average of each site's channels
    goodSite{s} = strcmp(channels.type, 'SEEG') & strcmp(channels.status, 'good') & ...
        ~ismember(upper(channels.name), getNeighborChs(split(sites{s}, '-'), 2));
    carRaw{s} = mean(dataAll(goodSite{s}, :), 1);
end
for c0 = 1:32:size(dataAll, 1)                            % per channel, as ieeg_highpass; in blocks for memory
    rows = c0:min(size(dataAll, 1), c0 + 31);
    dataAll(rows, :) = ieeg_highpass(dataAll(rows, :)', srate, 1)';
end
tt = (0:2*srate-1) / srate - 1;
ref = struct('site', sites, 'labels', [], 'trials', [], 'n_events', [], 'car', [], 'annotated_out', []);
for s = 1:numel(sites)
    ev = find(strcmp(events.electrical_stimulation_site, sites{s}));
    data = nan(height(channels), numel(tt), numel(ev));
    for k = 1:numel(ev)
        a = events.sample_start(ev(k)) - srate;
        data(:, :, k) = dataAll(:, a + 1:a + numel(tt));
    end
    goodChs = goodSite{s};
    [~, dataSite] = evalc('rmBadTrialsAnnots(data(goodChs, :, :), channels.name(goodChs), annots(ev))');
    keepTr = ~strcmp(upper(strip(annots(ev))), 'ALL');      % rmBadTrialsAnnots drops these trials
    nanChs = any(isnan(dataSite), [2, 3]);
    labs = channels.name(goodChs);
    ref(s).labels = labs(~nanChs);
    ref(s).annotated_out = nnz(nanChs);                      % channels bad in some trial, left out
    ref(s).trials = find(keepTr);                            % positions among the site's good events
    ref(s).n_events = numel(ev);
    ref(s).car = reshape(mean(dataSite(~nanChs, :, :), 1, 'omitnan'), numel(tt), []);   % time x trials (figS3)
end
clear dataAll data dataSite

% ---------------- iEEGLAB end to end ----------------
EEG = ieeglab_import(fullfile(ieeg, [stem '_ieeg.mefd']));
EEG = ieeglab_load(EEG, struct('event_field', 'electrical_stimulation_site'));
EEG = ieeglab_preprocess(EEG, struct('apply_highpass', true, 'highpass', 0.5, 'apply_epoch', true, ...
    'epoch_window', [-1000 1000], 'remove_bad_channels', true, 'plot', false));
EEG0 = EEG;
EEG = pop_ieeglab_reref(EEG, 'method', 'car', 'neighbors', 2);
G = EEG.ieeglab.car.groups;
epSites = ieeglab_epoch_sites(EEG0);
labels0 = {EEG0.chanlocs.labels};

labelsR = {EEG.chanlocs.labels};
away = tt(:) < -0.01 | tt(:) > 0.02;                        % outside the pulse and its artifact
rows = cell(numel(sites), 1);
for s = 1:numel(sites)
    cs = ieeglab_canonical_site(sites(s), EEG0.chanlocs);
    g = find(strcmp({G.site}, char(cs(1))), 1);
    tr = find(epSites == cs(1));
    isCens = all(censored(strcmp(events.electrical_stimulation_site, sites{s})));
    r = struct('site', string(sites{s}), 'censored', isCens, 'events', ref(s).n_events, ...
        'trials_reference', numel(ref(s).trials), 'trials_ieeglab', numel(tr), ...
        'channels_reference', numel(ref(s).labels), 'annotated_out_reference', ref(s).annotated_out, ...
        'channels_ieeglab', NaN, 'channels_both', NaN, 'same_channels', false, ...
        'r', NaN, 'r_outside_pulse', NaN, 'max_abs_diff_uv', NaN, 'max_abs_diff_outside_pulse_uv', NaN, 'rms_uv', NaN, ...
        'applied_spread_uv', NaN, 'r_same_fir', NaN, 'max_abs_diff_same_fir_uv', NaN, 'reference_linear_rel', NaN);
    if g > 0
        r.channels_ieeglab = numel(G(g).reference);
        r.channels_both = numel(intersect(upper(G(g).reference(:)), upper(ref(s).labels)));
        r.same_channels = r.channels_both == r.channels_ieeglab && r.channels_both == numel(ref(s).labels);
    end
    if ~isCens && ~isempty(g) && numel(tr) == ref(s).n_events   % the same epochs, in event order
        % what ieeglab_car subtracted: data before minus after, on every re-referenced contact
        [~, i0] = ismember(upper(G(g).reference(:)), upper(labels0));
        [~, i1] = ismember(upper(G(g).reference(:)), upper(labelsR));
        trK = tr(ref(s).trials);
        sub = double(EEG0.data(i0, :, trK)) - double(EEG.data(i1, :, trK));   % contacts x time x trials
        carI = reshape(mean(sub, 1), EEG0.pnts, []);
        r.applied_spread_uv = max(abs(sub - reshape(carI, 1, EEG0.pnts, [])), [], 'all');
        mR = mean(ref(s).car, 2); mI = mean(carI, 2);
        r.r = local_r(mR, mI); r.r_outside_pulse = local_r(mR(away), mI(away));
        r.max_abs_diff_uv = max(abs(mR - mI)); r.rms_uv = sqrt(mean(mR.^2));
        r.max_abs_diff_outside_pulse_uv = max(abs(mR(away) - mI(away)));
        if r.same_channels && ref(s).annotated_out == 0
            % the filter alone: the unfiltered common average of the same channels through each high-pass
            ev = find(strcmp(events.electrical_stimulation_site, sites{s}));
            ev = ev(ref(s).trials);
            epo = @(x) mean(cell2mat(arrayfun(@(e) x(events.sample_start(e) - srate + (1:numel(tt)))', ...
                ev', 'UniformOutput', false)), 2);
            mB = epo(ieeg_highpass(carRaw{s}', srate, 1)');            % the authors' filter
            F = eeg_emptyset(); F.data = carRaw{s}; [F.nbchan, F.pnts] = size(F.data); F.trials = 1;
            F.srate = srate; F.xmin = 0; F.xmax = (F.pnts - 1) / srate;
            [~, F] = evalc('pop_eegfiltnew(F, ''locutoff'', 0.5, ''usefftfilt'', 1, ''minphase'', false)');
            mF = epo(double(F.data));                                   % iEEGLAB's filter (ieeglab_preprocess)
            r.reference_linear_rel = max(abs(mB - mR)) / max(abs(mR));
            r.r_same_fir = local_r(mF, mI); r.max_abs_diff_same_fir_uv = max(abs(mF - mI));
        end
    end
    rows{s} = struct2table(r);
end
T = vertcat(rows{:});
ok = ~isnan(T.r);
fprintf(['\n%d sites, %d with data (the others are 0 in the public dataset). Same channels in the average: %d of ' ...
    'the sites with data, %d of all sites (the others differ by the channels the reference drops for per-trial ' ...
    'annotations). Trial-mean common average subtracted by iEEGLAB vs the reference: r %s (outside -10 to 20 ms %s), ' ...
    'largest difference %s uV (outside the pulse %s uV) for an RMS of %s uV. With the same FIR high-pass applied to ' ...
    'the unfiltered common average of the same channels: r %s, largest difference %s uV.\nDone in %.0f s.\n'], ...
    height(T), nnz(~T.censored), nnz(T.same_channels & ~T.censored), nnz(T.same_channels), ...
    local_list(T.r(ok), '%.5f'), local_list(T.r_outside_pulse(ok), '%.5f'), local_list(T.max_abs_diff_uv(ok), '%.1f'), ...
    local_list(T.max_abs_diff_outside_pulse_uv(ok), '%.1f'), local_list(T.rms_uv(ok), '%.1f'), ...
    local_list(T.r_same_fir(ok), '%.6f'), local_list(T.max_abs_diff_same_fir_uv(ok), '%.3f'), toc(t0));
disp(T);
end

function r = local_r(a, b)
c = corrcoef(a, b); r = c(1, 2);
end

function s = local_list(x, fmt)
s = strjoin(compose(fmt, x(:)'), ', ');
if isempty(s), s = 'none'; end
end
