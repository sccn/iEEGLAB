function tests = test_ieeglab_montage
% Bipolar and custom montages (ieeglab_montage, ieeglab_apply_montage,
% pop_ieeglab_reref 'bipolar' / 'montage'), and stimulation sites on derived
% channels. Headless.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
tc.TestData.root = root;
end

function test_bipolar_pairs_and_data(tc)
% Consecutive contacts of each lead; bad contacts, gaps and non-iEEG channels
% are not paired; data = first minus second; coordinates at the midpoint.
labels = {'RA1','RA2','RA3','RA4','RB1','RB2','RB3','RC1','RC2','RC4','ECG1'};
types  = [repmat({'SEEG'}, 1, 10) {'ECG'}];
xyz = [(0:3)' * 3.5, zeros(4, 2); (0:2)' * 3.5, 10 * ones(3, 1), zeros(3, 1); ...
       (0:2)' * 3.5, 20 * ones(3, 1), zeros(3, 1); 0 0 50];
EEG = local_eeg(labels, types, xyz, 500, 3);
EEG.chanlocs(6).status = 'bad';                              % RB2
EEG.icaweights = eye(11); EEG.icasphere = eye(11);
M = ieeglab_montage(EEG, 'bipolar');
tc.verifyEqual(M.labelnew, {'RA1-RA2','RA2-RA3','RA3-RA4','RC1-RC2'});
tc.verifyTrue(all(ismember({'RB1','RB3','RC4','ECG1','RB2'}, cellstr(M.skipped.item))), ...
    'Unpaired, bad and non-iEEG channels must be reported as skipped.');
[~, B] = evalc('ieeglab_apply_montage(EEG, M)');
tc.verifyEqual(B.nbchan, 4);
tc.verifyEqual({B.chanlocs.labels}, M.labelnew);
X = double(EEG.data); Y = double(B.data);
tc.verifyEqual(Y(2, :, :), X(2, :, :) - X(3, :, :), 'AbsTol', 1e-5);
tc.verifyEqual(Y(4, :, :), X(8, :, :) - X(9, :, :), 'AbsTol', 1e-5);
tc.verifyEqual([B.chanlocs(1).X B.chanlocs(1).Y B.chanlocs(1).Z], [1.75 0 0], 'AbsTol', 1e-12);
tc.verifyEqual(B.chanlocs(1).ref, 'bipolar');
tc.verifyEqual(B.ieeglab.montage.labelold, {'RA1','RA2','RA3','RA4','RC1','RC2'});
tc.verifyEqual(B.ieeglab.montage.tra(1, :), [1 -1 0 0 0 0], 'RA1-RA2 = RA1 - RA2');
tc.verifyEqual(B.ref, 'bipolar');
tc.verifyEmpty(B.icaweights, 'An ICA decomposition of the old channels must not survive.');
tc.verifyEqual(B.ieeglab.montage.method, 'bipolar');
end

function test_grid_row_end_not_paired(tc)
% A 2 x 4 grid numbered row by row: G4 and G5 are consecutive numbers but not
% neighbours. With coordinates they are not paired; spacing_factor 0 pairs them.
[gx, gy] = meshgrid(0:3, 0:1);
xyz = [reshape(gx', [], 1) reshape(gy', [], 1) zeros(8, 1)] * 10;
EEG = local_eeg(arrayfun(@(k) sprintf('G%d', k), 1:8, 'UniformOutput', false), repmat({'ECOG'}, 1, 8), xyz, 200, 1);
M = ieeglab_montage(EEG, 'bipolar');
tc.verifyFalse(ismember('G4-G5', M.labelnew));
tc.verifyEqual(numel(M.labelnew), 6);
M0 = ieeglab_montage(EEG, 'bipolar', 'spacing_factor', 0);
tc.verifyTrue(ismember('G4-G5', M0.labelnew));
end

function test_custom_montage_table_file_and_fieldtrip(tc)
% The same custom montage from a table, a TSV file and a FieldTrip struct:
% a bipolar pair, a local average (Laplacian-like) and a contact as it is. A
% row naming a contact not in the dataset is skipped.
labels = {'RA1','RA2','RA3','RA4'};
EEG = local_eeg(labels, repmat({'SEEG'}, 1, 4), (0:3)' * [3.5 0 0], 400, 2);
T = table(["RA1"; "RA2"; "RA4"; "RA9"], ["RA2"; "RA1 RA3"; "n/a"; "RA1"], ["RA1-RA2"; "RA2-local"; ""; ""], ...
    'VariableNames', {'positive', 'negative', 'name'});
M = ieeglab_montage(EEG, T);
tc.verifyEqual(M.labelnew, {'RA1-RA2', 'RA2-local', 'RA4'});
tc.verifyEqual(height(M.skipped), 1);
[~, B] = evalc('ieeglab_apply_montage(EEG, M)');
X = double(EEG.data); Y = double(B.data);
tc.verifyEqual(Y(2, :, :), X(2, :, :) - (X(1, :, :) + X(3, :, :)) / 2, 'AbsTol', 1e-5);
tc.verifyEqual(Y(3, :, :), X(4, :, :), 'AbsTol', 1e-6);

f = [tempname '.tsv'];
c = onCleanup(@() delete(f));
fid = fopen(f, 'w');
fprintf(fid, 'positive\tnegative\tname\nRA1\tRA2\tRA1-RA2\nRA2\tRA1 RA3\tRA2-local\nRA4\tn/a\t\nRA9\tRA1\t\n');
fclose(fid);
Mf = ieeglab_montage(EEG, f);
tc.verifyEqual(Mf.labelnew, M.labelnew);
tc.verifyEqual(Mf.tra, M.tra, 'AbsTol', 1e-15);

ft = struct('labelold', {M.labelold}, 'labelnew', {M.labelnew}, 'tra', M.tra);
[~, Bft] = evalc('ieeglab_apply_montage(EEG, ieeglab_montage(EEG, ft))');
tc.verifyEqual(double(Bft.data), Y, 'AbsTol', 1e-6);
tc.verifyError(@() ieeglab_montage(EEG, table("RA9", "RA8", 'VariableNames', {'positive', 'negative'})), ...
    'ieeglab_montage:nothingMatched');
end

function test_stimulation_sites_on_bipolar_channels(tc)
% A channel that uses a stimulated contact is stimulated; the dataset stays
% CCEP; sites are named by contacts as before.
labels = {'RA1','RA2','RA3','RA4','RA5'};
EEG = local_eeg(labels, repmat({'SEEG'}, 1, 5), (0:4)' * [3.5 0 0], 300, 6, {'RA2-RA3', 'RA3-RA2', 'RA4-RA5'});
[~, B] = evalc('ieeglab_apply_montage(EEG, ieeglab_montage(EEG, ''bipolar''))');
tc.verifyEqual(ieeglab_detect_mode(B), 'ccep');
[sites, stimIdx] = ieeglab_epoch_sites(B);
tc.verifyEqual(sites(1:3), ["RA2-RA3" "RA2-RA3" "RA4-RA5"]);
tc.verifyEqual(sort(stimIdx{1})', [1 2 3], 'RA1-RA2, RA2-RA3 and RA3-RA4 use RA2 or RA3.');
tc.verifyEqual(sort(stimIdx{3})', [3 4], 'RA3-RA4 and RA4-RA5 use RA4 or RA5.');
% referential labels behave as before
[s0, i0] = ieeglab_canonical_site({'RA3-RA2'}, labels);
tc.verifyEqual(s0, "RA2-RA3"); tc.verifyEqual(i0{1}, [2; 3]);
end

function test_reref_dialogless_on_tutorial(tc)
% pop_ieeglab_reref: bipolar on the epoched tutorial and on continuous data,
% history line, and no second re-referencing on top of a montage.
d = fullfile(tc.TestData.root, 'tutorial', 'seeg');
EEG = pop_loadset('filename', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set', 'filepath', d);
[~, EEG] = evalc('ieeglab_load(EEG, struct())');
[~, Bc, comc] = evalc('pop_ieeglab_reref(EEG, ''method'', ''bipolar'')');
tc.verifyEqual(Bc.ref, 'bipolar');
tc.verifyTrue(contains(comc, '''bipolar'''));
pre = struct('apply_highpass',true,'highpass',0.5,'apply_notch',false,'apply_lowpass',false, ...
    'apply_ds',false,'apply_epoch',true,'epoch_window',[-500 1000],'remove_bad_channels',true, ...
    'apply_baseline',true,'baseline_period',[-500 -50],'apply_car',false,'plot',false,'verbose',false);
[~, E] = evalc('ieeglab_preprocess(EEG, pre)');
[~, B, com] = evalc('pop_ieeglab_reref(E, ''method'', ''bipolar'')');
tc.verifyTrue(all(contains({B.chanlocs.labels}, '-')));
EEG = E;                                                      % replay the history line on E
[~, B2] = evalc(com(7:end-1));
tc.verifyEqual(double(B2.data), double(B.data));
i = find(strcmp({B.chanlocs.labels}, B.ieeglab.montage.labelnew{1}));
a = strcmp({E.chanlocs.labels}, B.ieeglab.montage.labelold{find(B.ieeglab.montage.tra(1, :) > 0)}); %#ok<FNDSB>
b = strcmp({E.chanlocs.labels}, B.ieeglab.montage.labelold{find(B.ieeglab.montage.tra(1, :) < 0)}); %#ok<FNDSB>
D = double(E.data(a, :, :)) - double(E.data(b, :, :));
kb = local_baseline_idx(E.times, pre.baseline_period);
tc.verifyEqual(double(B.data(i, :, :)), D - median(D(:, kb, :), 2), 'AbsTol', 1e-3, ...
    'The difference of the two contacts, with the (median) baseline applied again.');
[~, stimIdx] = ieeglab_epoch_sites(B);
tc.verifyTrue(any(~cellfun(@isempty, stimIdx)), 'Stimulated bipolar channels must still be found.');
tc.verifyError(@() pop_ieeglab_reref(B, 'method', 'car'), 'pop_ieeglab_reref:afterMontage');
end

function test_prime_leads_are_separate(tc)
% sEEG writes the other hemisphere's electrode with a prime: A'1 is not on
% lead A. No pair may cross from one to the other.
labels = {'A1','A2','A3','A''1','A''2','A''3'};
xyz = [30 0 0; 30 3.5 0; 30 7 0; -30 0 0; -30 3.5 0; -30 7 0];
EEG = local_eeg(labels, repmat({'SEEG'}, 1, 6), xyz, 200, 1);
M = ieeglab_montage(EEG, 'bipolar');
tc.verifyEqual(sort(M.labelnew), sort({'A1-A2', 'A2-A3', 'A''1-A''2', 'A''2-A''3'}));
[lead, num] = ieeglab_lead_number({'A''1', 'LA''12', 'ROP 3', 'RA1', 'EKG', [char(65) char(8217) '2']});
tc.verifyEqual(lead, ["A'" "LA'" "ROP" "RA" "" "A'"]);
tc.verifyEqual(num, [1 12 3 1 NaN 2]);
end

function test_custom_names_still_stimulated(tc)
% Channels named freely in a custom montage still resolve stimulation sites
% through the contacts they use (chanlocs.contacts).
labels = {'RA1','RA2','RA3','RA4'};
EEG = local_eeg(labels, repmat({'SEEG'}, 1, 4), (0:3)' * [3.5 0 0], 300, 4, {'RA1-RA2', 'RA3-RA4'});
T = table(["RA3"; "RA4"], ["RA2 RA4"; "RA1"], ["RA3_lap"; "far"], 'VariableNames', {'positive', 'negative', 'name'});
[~, B] = evalc('ieeglab_apply_montage(EEG, ieeglab_montage(EEG, T))');
tc.verifyEqual({B.chanlocs.contacts}, {'RA3, RA2, RA4', 'RA4, RA1'});
tc.verifyEqual(ieeglab_detect_mode(B), 'ccep');
[sites, stimIdx] = ieeglab_epoch_sites(B);
tc.verifyEqual(sites(1:2), ["RA1-RA2" "RA3-RA4"]);
tc.verifyEqual(sort(stimIdx{1})', [1 2], 'RA3_lap uses RA2, far uses RA1.');
tc.verifyEqual(sort(stimIdx{2})', [1 2], 'Both use RA3 or RA4.');
end

function test_negatives_with_spaces_and_repeats(tc)
% Labels with a space ('ROP 1') in the negative column, separated by commas;
% a contact listed twice weighs twice.
labels = {'ROP 1', 'ROP 2', 'ROP 3'};
EEG = local_eeg(labels, repmat({'SEEG'}, 1, 3), (0:2)' * [3.5 0 0], 200, 1);
M = ieeglab_montage(EEG, table(["ROP 2"; "ROP 1"], ["ROP 1, ROP 3"; "ROP 2"], 'VariableNames', {'positive', 'negative'}));
tc.verifyEqual(M.labelnew, {'ROP 2-ROP 1-ROP 3', 'ROP 1-ROP 2'});
tc.verifyEqual(M.tra, [-0.5 1 -0.5; 1 -1 0], 'AbsTol', 1e-15);
M2 = ieeglab_montage(EEG, table("ROP 1", "ROP 2, ROP 2, ROP 3", 'VariableNames', {'positive', 'negative'}));
tc.verifyEqual(M2.tra, [1 -2/3 -1/3], 'AbsTol', 1e-15);
end

function test_bipolar_matches_fieldtrip(tc)
% Same pairs and weights as FieldTrip's bipolar montage per lead
% (ft_prepare_montage, refmethod 'bipolar', groupchans 'yes') on the good
% contacts of both tutorials. FieldTrip pairs every consecutive pair; iEEGLAB
% also skips the pairs farther apart than twice the lead's spacing.
tc.assumeNotEmpty(which('ft_prepare_montage'), 'FieldTrip is not on the path.');
sets = {'seeg', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set'; ...
        'ecog', 'sub-ccepAgeUMCU02_ses-1_task-SPESclin_run-041456_ieeg.set'};
for s = 1:2
    EEG = pop_loadset('filename', sets{s, 2}, 'filepath', fullfile(tc.TestData.root, 'tutorial', sets{s, 1}));
    [~, EEG] = evalc('ieeglab_load(EEG, struct())');
    M = ieeglab_montage(EEG, 'bipolar');
    good = {EEG.chanlocs(~strcmpi({EEG.chanlocs.status}, 'bad')).labels};
    cfg = struct('refmethod', 'bipolar', 'groupchans', 'yes', 'channel', {good(:)});
    [~, F] = evalc('ft_prepare_montage(cfg)');
    extra = setdiff(F.labelnew, M.labelnew);
    tc.verifyEmpty(setdiff(M.labelnew, F.labelnew), sprintf('%s: pairs FieldTrip does not make.', sets{s, 1}));
    tc.verifyTrue(all(ismember(extra, cellstr(M.skipped.item(contains(M.skipped.reason, 'apart'))))), ...
        sprintf('%s: FieldTrip pairs that iEEGLAB skips must be the distant ones.', sets{s, 1}));
    [~, io] = ismember(M.labelold, F.labelold);
    [~, in] = ismember(M.labelnew, F.labelnew);
    tc.verifyEqual(M.tra, F.tra(in, io), sprintf('%s: weights differ from FieldTrip''s.', sets{s, 1}));
end
tc.verifyNotEmpty(extra, 'The ECoG tutorial has two distant pairs (FB56-FB57, FB63-FB64).');
end

function test_montage_commutes_with_preprocessing(tc)
% Bipolar then preprocessing = preprocessing then bipolar, for the linear steps
% (epochs, median baseline, high-pass, notch filter, linear blanking). The
% baseline is a median: re-referencing must apply it again after the montage.
% Not linear, so not compared: pchip blanking (differs inside the blanked
% window only) and the sine fit (searches the mains frequency over all channels).
d = fullfile(tc.TestData.root, 'tutorial', 'seeg');
EEG = pop_loadset('filename', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set', 'filepath', d);
[~, EEG] = evalc('ieeglab_load(EEG, struct())');
pre = struct('apply_blank', true, 'blank_method', 'linear', 'apply_highpass', true, 'highpass', 0.5, ...
    'apply_notch', true, 'notch', 60, 'line_method', 'notch', 'apply_epoch', true, 'epoch_window', [-1000 1000], ...
    'remove_bad_channels', true, 'remove_no_coords', true, 'apply_baseline', true, 'baseline_period', [-500 -50], ...
    'plot', false, 'verbose', false);
[~, A] = evalc('ieeglab_preprocess(EEG, pre)');
[~, A] = evalc('pop_ieeglab_reref(A, ''method'', ''bipolar'')');
[~, B] = evalc('pop_ieeglab_reref(EEG, ''method'', ''bipolar'')');
[~, B] = evalc('ieeglab_preprocess(B, pre)');
tc.verifyEqual({A.chanlocs.labels}, {B.chanlocs.labels});
tol = 8 * double(eps(single(max(abs(A.data(:))))));
tc.verifyEqual(double(A.data), double(B.data), 'AbsTol', tol);
k = local_baseline_idx(A.times, pre.baseline_period);
tc.verifyEqual(median(double(A.data(:, k, :)), 2), zeros(A.nbchan, 1, A.trials), 'AbsTol', 1e-3, ...
    'Each bipolar channel''s baseline median must be 0.');
end

function test_analysis_after_bipolar(tc)
% The whole analysis runs on bipolar channels of both tutorials: CCEP statistics
% leave out every channel that uses a stimulated contact, values reach the
% brain plot, export writes the bipolar labels, and band power runs.
sets = {'seeg', 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set', 60; ...
        'ecog', 'sub-ccepAgeUMCU02_ses-1_task-SPESclin_run-041456_ieeg.set', 50};
for s = 1:2
    EEG = pop_loadset('filename', sets{s, 2}, 'filepath', fullfile(tc.TestData.root, 'tutorial', sets{s, 1}));
    [~, EEG] = evalc('ieeglab_load(EEG, struct())');
    pre = struct('apply_blank', true, 'apply_notch', true, 'notch', sets{s, 3}, 'apply_epoch', true, ...
        'epoch_window', [-1000 1000], 'remove_rare_cond', true, 'min_trials', 5, 'remove_no_coords', true, ...
        'remove_bad_channels', true, 'apply_baseline', true, 'baseline_period', [-500 -50], 'plot', false, 'verbose', false);
    [~, E] = evalc('ieeglab_preprocess(EEG, pre)');
    [~, E] = evalc('pop_ieeglab_reref(E, ''method'', ''bipolar'')');
    % each modality's default: CRP for sEEG, N1 for ECoG
    [~, E] = evalc('ieeglab_stats_subject(E, struct(''run_matrix'', true, ''min_trials'', 5, ''verbose'', false))');
    if s == 1
        tc.verifyTrue(isfield(E.ieeglab, 'stats') && height(E.ieeglab.stats.table) > 0, 'sEEG: no CRP statistics.');
    else
        tc.verifyTrue(isfield(E.ieeglab, 'n1') && height(E.ieeglab.n1.table) > 0, 'ECoG: no N1 statistics.');
    end
    tc.verifyTrue(isfield(E.ieeglab, 'ccep_matrix'), sprintf('%s: no connectivity matrix.', sets{s, 1}));
    [sites, stimIdx] = ieeglab_epoch_sites(E);
    [us, iu] = unique(sites, 'stable');
    labels = {E.chanlocs.labels};
    for k = 1:numel(us)
        stim = labels(stimIdx{iu(k)});
        pair = regexp(char(us(k)), '[A-Za-z'' ]+\d+', 'match');
        uses = labels(cellfun(@(l) any(ismember(regexp(l, '[A-Za-z'' ]+\d+', 'match'), pair)), labels));
        tc.verifyEqual(sort(stim(:)), sort(uses(:)), sprintf('%s, site %s: stimulated channels.', sets{s, 1}, us(k)));
    end
    v = ieeglab_topoplot(E, 'in_degree', 'draw', false);
    tc.verifyNumElements(v, E.nbchan);
    out = tempname; mkdir(out); c = onCleanup(@() rmdir(out, 's'));
    [~, files] = evalc('ieeglab_export(E, out, struct(''formats'', {{''tsv''}}, ''verbose'', false))');
    txt = strjoin(cellfun(@fileread, files(endsWith(files, '.tsv')), 'UniformOutput', false), newline);
    tc.verifyTrue(contains(txt, labels{1}), sprintf('%s: export must carry the bipolar labels.', sets{s, 1}));
    [~, S] = evalc('ieeglab_spectral(E, ''freqs'', 70:20:150, ''verbose'', false)');
    tc.verifyEqual(numel(S.band_value), E.nbchan);
    tc.verifyTrue(any(isfinite(S.band_value)));
    clear c
end
end

% ------------------------------------------------------------------ helpers

function k = local_baseline_idx(times, period)
% Baseline samples as ieeglab_rm_baseline takes them: nearest samples to the ends
[~, i0] = min(abs(double(times) - period(1)));
[~, i1] = min(abs(double(times) - period(2)));
k = i0:i1;
end

function EEG = local_eeg(labels, types, xyz, T, N, epochTypes)
% Synthetic dataset: random data, given channels; epoched (N > 1) with one
% event per epoch at time 0 when epochTypes is given.
if nargin < 6, epochTypes = {}; end
C = numel(labels);
rng(7);
EEG = eeg_emptyset();
EEG.srate = 1000; EEG.nbchan = C; EEG.pnts = T; EEG.trials = N;
EEG.data = randn(C, T, N, 'single') * 20;
EEG.xmin = -0.1; EEG.xmax = EEG.xmin + (T - 1) / EEG.srate;
EEG.times = (EEG.xmin + (0:T-1) / EEG.srate) * 1000;
for c = C:-1:1
    cl(c) = struct('labels', labels{c}, 'type', types{c}, 'X', xyz(c, 1), 'Y', xyz(c, 2), 'Z', xyz(c, 3), ...
        'status', 'good', 'theta', 10, 'radius', 0.5); %#ok<AGROW>
end
EEG.chanlocs = cl;
EEG.event = struct('type', {}, 'latency', {}, 'epoch', {});
if ~isempty(epochTypes)
    i0 = find(EEG.times >= 0, 1);
    for n = 1:N
        t = epochTypes{mod(n - 1, numel(epochTypes)) + 1};
        EEG.event(n) = struct('type', t, 'latency', (n - 1) * T + i0, 'epoch', n);
        EEG.epoch(n).event = n; EEG.epoch(n).eventtype = {t}; EEG.epoch(n).eventlatency = {0};
    end
end
end
