function out = validate_fieldtrip_uci29(dataDir)
% validate_fieldtrip_uci29() - iEEGLAB against FieldTrip on FieldTrip's
%       intracranial tutorial dataset (SubjectUCI29; Stolk et al., 2018):
%       re-referencing (common average of the grids, bipolar depth leads),
%       time-frequency power and high-frequency activity, then iEEGLAB's
%       permutation statistics for band power, which the tutorial does not have.
%
% Usage:
%   out = validate_fieldtrip_uci29('C:\data\SubjectUCI29')
%
% dataDir holds, from download.fieldtriptoolbox.org/tutorial/human_ecog/SubjectUCI29/:
%   SubjectUCI29_data.mat         the epochs (26 trials, 152 contacts, 5000 Hz)
%   SubjectUCI29_elec_acpc_f.mat  contact positions (ACPC, mm)
%   SubjectUCI29_freq.mat         FieldTrip's result of the tutorial (ft_freqanalysis)
%
% The tutorial (fieldtriptoolbox.org/tutorial/intracranial/human_ecog): common
% average of the two grids together (LPG*, LTG*), bipolar montage of each depth
% lead, mtmconvol with Hann windows of 0.2 s at 5:5:200 Hz every 10 ms from
% -0.3 to 0.8 s, baseline -0.3 to -0.1 s ('relchange'), HFA 70-150 Hz over 0 to
% 0.8 s. iEEGLAB runs the same steps through its own functions, at FieldTrip's
% frequencies (FieldTrip moves 5:5:200 onto its padded frequency grid). The
% same values are selected on both sides: FieldTrip's baseline holds the 20
% window centres -300 to -110 ms (its time axis has -100 ms at -0.0999...,
% outside [-0.3 -0.1]), which iEEGLAB's baseline [-400 -10] gives with 0.2 s
% windows; ft_selectdata picks the band and latency ends with nearest(), so the
% HFA band starts at 69.98 Hz.
%
% Cedric Cannard, iEEGLAB, 2026

t0 = tic;
f = @(n) fullfile(dataDir, ['SubjectUCI29_' n '.mat']);
types = {'LPG*', 'ECOG'; 'LTG*', 'ECOG'; '*', 'SEEG'};
EEG = ieeglab_import_fieldtrip(f('data'), 'elec', f('elec_acpc_f'), 'types', types);
ft = load(f('freq')); ft = ft.freq;

% ---------------- re-referencing, as the tutorial ----------------
isGrid = startsWith({EEG.chanlocs.labels}, {'LPG', 'LTG'});
G = pop_select(EEG, 'channel', find(isGrid));
G = pop_ieeglab_reref(G, 'method', 'car');
D = pop_select(EEG, 'channel', find(~isGrid));
D = pop_ieeglab_reref(D, 'method', 'bipolar', 'spacing_factor', 0);
R = G;
R.data = cat(1, G.data, D.data);
R.chanlocs = local_cat_locs(G.chanlocs, D.chanlocs);
R.nbchan = size(R.data, 1);
R = eeg_checkset(R);
labels = {R.chanlocs.labels};
[inFt, loc] = ismember(labels, ft.label);
out.n_channels = [numel(labels), numel(ft.label), nnz(inFt)];
fprintf('Channels: iEEGLAB %d, FieldTrip %d, in both %d. Only in iEEGLAB: %s. Only in FieldTrip: %s\n', ...
    out.n_channels, strjoin(labels(~inFt), ' '), strjoin(setdiff(ft.label, labels)', ' '));

% ---------------- time-frequency at FieldTrip's frequencies ----------------
[~, f1] = min(abs(ft.freq - 70)); [~, f2] = min(abs(ft.freq - 150));    % ft_selectdata: nearest()
band = ft.freq([f1 f2]);
S = ieeglab_spectral(R, 'type', 'timefreq', 'freqs', ft.freq, 'window', 0.2, 'timestep', 10, ...
    'baseline', [-400 -10], 'norm', 'relchange', 'band', band, 'band_window', [0 800], 'keeptrials', true);
[tIn, tLoc] = ismember(round(ft.time * 1000), round(S.times));
A = S.tf(inFt, :, tLoc(tIn));                                % iEEGLAB relchange (ch x freq x time)
P = ft.powspctrm(loc(inFt), :, tIn);
tFt = ft.time(tIn);
inBaseFt = tFt >= -0.3 & tFt <= -0.1;                        % ft_freqbaseline's own test
base = mean(P(:, :, inBaseFt), 3, 'omitnan');
B = P ./ base - 1;                                           % FieldTrip ft_freqbaseline 'relchange'
hMs = (S.window_s(1) * S.srate - 1) / 2 / S.srate * 1000;
out.baseline_centres = [nnz(S.times >= -400 + hMs & S.times <= -10 - hMs), nnz(inBaseFt)];
if out.baseline_centres(1) ~= out.baseline_centres(2)
    warning('validate_fieldtrip_uci29:baseline', 'Baselines differ: iEEGLAB %d window centres, FieldTrip %d.', ...
        out.baseline_centres);
end
ok = isfinite(A) & isfinite(B);
r = local_r(A(ok), B(ok));
out.tf_r = r;
out.tf_max_abs_diff = max(abs(A(ok) - B(ok)));
out.tf_median_abs_diff = median(abs(A(ok) - B(ok)));
fprintf('Time-frequency (relchange), %d values, baseline %d window centres on both sides: r = %.6f, median |diff| %.2g, max |diff| %.2g\n', ...
    nnz(ok), out.baseline_centres(1), r, out.tf_median_abs_diff, out.tf_max_abs_diff);

% HFA over 0-0.8 s, one value per channel, as the tutorial's brain map (ends by nearest())
fB = false(size(ft.freq)); fB(f1:f2) = true;
[~, t1] = min(abs(tFt - 0)); [~, t2] = min(abs(tFt - 0.8));
tB = false(size(tFt)); tB(t1:t2) = true;
hA = mean(mean(A(:, fB, tB), 2, 'omitnan'), 3, 'omitnan');
hB = mean(mean(B(:, fB, tB), 2, 'omitnan'), 3, 'omitnan');
out.hfa_band_hz = band;
out.hfa_r = local_r(hA, hB);
out.hfa_max_abs_diff = max(abs(hA - hB));
[~, ia] = sort(hA, 'descend'); [~, ib] = sort(hB, 'descend');
out.hfa_top10_iEEGLAB = labels(find(inFt)) ; out.hfa_top10_iEEGLAB = out.hfa_top10_iEEGLAB(ia(1:10));
fl = ft.label(loc(inFt)); out.hfa_top10_FieldTrip = fl(ib(1:10))';
out.hfa_same_top10 = isequal(out.hfa_top10_iEEGLAB, out.hfa_top10_FieldTrip);
fprintf('HFA %.2f-%.2f Hz per channel: r = %.6f, max |diff| %.3g (relative change); same 10 strongest in the same order %d. Top 5 iEEGLAB: %s; FieldTrip: %s\n', ...
    band, out.hfa_r, out.hfa_max_abs_diff, out.hfa_same_top10, strjoin(out.hfa_top10_iEEGLAB(1:5), ' '), ...
    strjoin(out.hfa_top10_FieldTrip(1:5), ' '));

% ---------------- iEEGLAB's statistics, not in the tutorial ----------------
[T, Rs] = ieeglab_band_stats(S, 'window', [0 800], 'baseline', [-400 -10]);
out.stats = T;
out.stats_info = Rs;
s = T(T.significant, :);
[~, o] = sort(s.band_db, 'descend');
fprintf('HFA vs baseline, sign-flip permutation, FDR: %d of %d channels (%d up, %d down). Largest: %s\n', ...
    height(s), nnz(~isnan(T.p)), nnz(s.band_db > 0), nnz(s.band_db < 0), ...
    strjoin(compose('%s %.2f dB', string(s.channel(o(1:min(end, 8)))), s.band_db(o(1:min(end, 8)))), ', '));
out.seconds = toc(t0);
fprintf('Done in %.0f s.\n', out.seconds);
end

function r = local_r(a, b)
c = corrcoef(a, b); r = c(1, 2);
end

function L = local_cat_locs(a, b)
% Concatenate two chanlocs structs that may not share every field
fn = union(fieldnames(a), fieldnames(b), 'stable');
for k = 1:numel(fn)
    if ~isfield(a, fn{k}), [a.(fn{k})] = deal([]); end
    if ~isfield(b, fn{k}), [b.(fn{k})] = deal([]); end
end
L = [orderfields(a(:), fn); orderfields(b(:), fn)]';
end
