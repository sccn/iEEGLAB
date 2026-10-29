function tests = test_ieeglab_reref
% iEEG re-referencing as its own step (pop_ieeglab_reref, issue #11). Headless.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
d = fullfile(root, 'tutorial', 'seeg');
EEG = pop_loadset('filename','sub-02_ses-ieeg01_task-ccep_run-01_ieeg.set','filepath',d);
[~, EEG] = evalc('ieeglab_load(EEG, struct())');
tc.TestData.EEG = EEG;
tc.TestData.pre = struct('apply_highpass',true,'highpass',0.5,'apply_notch',false,'apply_lowpass',false, ...
    'apply_ds',false,'apply_epoch',true,'epoch_window',[-500 1000],'remove_bad_channels',true, ...
    'apply_baseline',true,'baseline_period',[-500 -50],'plot',false,'verbose',false);
end

function test_same_as_inside_preprocessing(tc)
% Re-referencing after baseline correction must equal re-referencing before it
% (the order ieeglab_preprocess uses), for CARLA and for the plain average.
for m = {'carla', 'car'}
    a = tc.TestData.pre; a.apply_car = true; a.car_method = m{1};
    rng(1); [~, A] = evalc('ieeglab_preprocess(tc.TestData.EEG, a)');
    b = tc.TestData.pre; b.apply_car = false;
    [~, B] = evalc('ieeglab_preprocess(tc.TestData.EEG, b)');
    rng(1); [~, B, com] = evalc('pop_ieeglab_reref(B, ''method'', m{1})');
    % single-precision data: tolerance scaled to its rounding step (CI on Linux
    % rounds differently from Windows), near the artifact values reach 3e4 uV
    tol = 8 * double(eps(single(max(abs(A.data(:))))));
    tc.verifyEqual(double(B.data), double(A.data), 'AbsTol', tol, ...
        sprintf('%s after baseline differs from %s inside preprocessing.', m{1}, m{1}));
    tc.verifyEqual(B.ref, A.ref);
    tc.verifyTrue(contains(com, ['''method'', ''' m{1} '''']), 'The history line must replay the method.');
end
end

function test_needs_epochs(tc)
tc.verifyError(@() pop_ieeglab_reref(tc.TestData.EEG, 'method', 'car'), 'pop_ieeglab_reref:notEpoched');
end

function test_rejects_bad_options(tc)
b = tc.TestData.pre; b.apply_car = false;
[~, B] = evalc('ieeglab_preprocess(tc.TestData.EEG, b)');
tc.verifyError(@() pop_ieeglab_reref(B, 'method', 'laplacian'), 'pop_ieeglab_reref:method');
tc.verifyError(@() pop_ieeglab_reref(B, 'method', 'montage'), 'pop_ieeglab_reref:noMontage');
tc.verifyError(@() pop_ieeglab_reref(B, 'method', 'varsubset', 'fraction', 2), 'pop_ieeglab_reref:fraction');
tc.verifyError(@() pop_ieeglab_reref(B, 'window', [10 300]), 'pop_ieeglab_reref:unknownOption');
end

function test_method_defaults_not_inherited(tc)
% Each method's published defaults: values remembered from a run with another
% method (CARLA: 10-300 ms, 2 neighbours) must not carry over to varsubset, and
% the history line replays the neighbours whatever the method.
b = tc.TestData.pre; b.apply_car = false;
[~, B] = evalc('ieeglab_preprocess(tc.TestData.EEG, b)');
[~, R1, com1] = evalc('pop_ieeglab_reref(B, ''method'', ''carla'')');
tc.verifyEqual(R1.ieeglab.car.neighbors, 2);
tc.verifyTrue(contains(com1, '''neighbors'', 2'));
[~, R2, com2] = evalc('pop_ieeglab_reref(R1, ''method'', ''varsubset'')');
tc.verifyEqual(R2.ieeglab.car.neighbors, 0);
tc.verifyEqual(R2.ieeglab.car.timewin_ms, [15 500]);
tc.verifyTrue(contains(com2, '''neighbors'', 0') && contains(com2, '[15 500]'));
end

function test_neighbors_left_out(tc)
% 'neighbors', n keeps the contacts within n of the stimulated pair (same lead)
% out of every site's reference, as the CARLA publication scripts do.
b = tc.TestData.pre; b.apply_car = false;
[~, B] = evalc('ieeglab_preprocess(tc.TestData.EEG, b)');
[~, R0] = evalc('pop_ieeglab_reref(B, ''method'', ''car'')');
[~, R2, com] = evalc('pop_ieeglab_reref(B, ''method'', ''car'', ''neighbors'', 2)');
tc.verifyTrue(contains(com, '''neighbors'', 2'), 'The history line must replay the option.');
[lead, num] = local_lead_number({B.chanlocs.labels});
for g = R2.ieeglab.car.groups
    tok = regexp(g.site, '([A-Za-z]+)(\d+)', 'tokens');
    near = false(size(lead));
    for t = 1:numel(tok)
        near = near | (lead == upper(tok{t}{1}) & abs(num - str2double(tok{t}{2})) <= 2);
    end
    tc.verifyEmpty(intersect(g.reference, {B.chanlocs(near).labels}), ...
        sprintf('Site %s: a contact within 2 of the stimulated pair is in the reference.', g.site));
    g0 = R0.ieeglab.car.groups(strcmp({R0.ieeglab.car.groups.site}, g.site));
    tc.verifyEqual(sort(g.reference(:)), sort(setdiff(g0.reference, {B.chanlocs(near).labels})'), ...
        'Only the neighbours may leave the plain common average.');
end
tc.verifyError(@() pop_ieeglab_reref(B, 'method', 'car', 'neighbors', -1), 'pop_ieeglab_reref:neighbors');
end

function test_continuous_car_and_ica(tc)
% Continuous data without stimulation: the common average and ICA are offered.
% Seven good contacts share a reference signal (equal weight on each) on top
% of their own sources; LA8 is marked bad.
fs = 500; T = 60 * fs; rng(4);
common = 3 * randn(1, T).^3;                      % super-Gaussian, as ICA needs
S = randn(6, T).^3 .* (1 + (0:5)' / 5);              % own sources, 1 to 2 x
X = [S; 0.5 * S(6, :)] + common;                   % 7 contacts, 7 sources
X(8, :) = 10 * randn(1, T);
EEG = eeg_emptyset();
EEG.srate = fs; EEG.nbchan = 8; EEG.pnts = T; EEG.trials = 1;
EEG.data = single(X); EEG.xmin = 0; EEG.xmax = (T - 1) / fs; EEG.times = (0:T-1) / fs * 1000;
for c = 8:-1:1, cl(c) = struct('labels', sprintf('LA%d', c), 'type', 'SEEG', 'status', 'good'); end
cl(8).status = 'bad';
EEG.chanlocs = cl;

[~, B, com] = evalc('pop_ieeglab_reref(EEG, ''method'', ''car'')');
ref = mean(double(EEG.data(1:7, :)), 1);
tc.verifyEqual(double(B.data), double(EEG.data) - ref, 'AbsTol', 1e-4, 'Every channel minus the mean of the 7 good ones.');
tc.verifyEqual(B.ref, 'CAR');
tc.verifyTrue(contains(com, '''car'''));

[~, I] = evalc('pop_ieeglab_reref(EEG, ''method'', ''ica'')');
r0 = abs(corr(double(EEG.data(1:7, :))', common'));
r1 = abs(corr(double(I.data(1:7, :))', common'));
tc.verifyGreaterThan(min(r0), 0.5);
tc.verifyLessThan(max(r1), 0.2, 'The shared reference signal must be removed.');
tc.verifyEqual(I.data(8, :), EEG.data(8, :), 'The bad contact must be left untouched.');
tc.verifyError(@() pop_ieeglab_reref(EEG, 'method', 'carla'), 'pop_ieeglab_reref:notEpoched');
end

function [lead, num] = local_lead_number(labels)
lead = string(upper(regexprep(labels, '\d+$', '')));
num = cellfun(@(s) str2double(regexp(s, '\d+$', 'match', 'once')), labels);
end
