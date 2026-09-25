function tests = test_ieeglab_reref
% iEEG re-referencing as its own step (pop_ieeglab_reref, issue #11). Headless.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
tc.assumeNotEmpty(which('eeglab'), 'EEGLAB is not on the MATLAB path.');
if isempty(which('pop_loadset')), evalc('eeglab nogui'); end
addpath(root); addpath(fullfile(root, 'functions'));
d = fullfile(root, 'tutorial', 'dataset_seeg');
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
tc.verifyError(@() pop_ieeglab_reref(B, 'method', 'bipolar'), 'pop_ieeglab_reref:method');
tc.verifyError(@() pop_ieeglab_reref(B, 'method', 'varsubset', 'fraction', 2), 'pop_ieeglab_reref:fraction');
tc.verifyError(@() pop_ieeglab_reref(B, 'window', [10 300]), 'pop_ieeglab_reref:unknownOption');
end
