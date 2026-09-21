function [EEG, info] = ieeglab_blank_stim(EEG, opt)
% ieeglab_blank_stim() - Remove the stimulation artifact before filtering.
%
% Usage:
%   EEG        = ieeglab_blank_stim(EEG)
%   [EEG, info] = ieeglab_blank_stim(EEG, opt)
%
% Options:
%   .blank_window  [1x2] ms around each stimulation onset to replace.
%                  Default [-1.5 10]. Widen if the artifact or its recovery
%                  transient extends further.
%   .blank_method  'pchip'  - shape-preserving interpolation across the gap (default)
%                  'linear' - straight line across the gap
%                  'nan'    - leave the gap as NaN
%   .verbose       logical. Default true.
%
% Why this is needed
% ------------------
% The stimulation artifact is a near-instantaneous step of tens of millivolts.
% Any subsequent filter sees that as a broadband impulse. A zero-phase FIR
% notch, which is what this pipeline applies, rings SYMMETRICALLY about it, so
% energy from the artifact is smeared both forward into the early evoked
% response and BACKWARD into the pre-stimulus baseline. Baseline-correcting
% against a contaminated baseline then propagates the error through the whole
% epoch. This is the mechanism behind the ringing reported in issue #10.
%
% The fix is ordering: blank the artifact on the CONTINUOUS data first, then
% filter. The filter then never sees the discontinuity.
%
% Blanking is destructive by design - the samples inside the window carry no
% recoverable neural signal, only the stimulus. Any analysis window should
% start after blank_window(2); the CRP and N1 defaults (15 ms and 10 ms) are
% chosen with that in mind.
%
% Cedric Cannard, iEEGLAB, 2026

if nargin < 2 || isempty(opt)
    if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isstruct(EEG.ieeglab.opt)
        opt = EEG.ieeglab.opt;
    else
        opt = struct();
    end
end
def = struct('blank_window',[-1.5 10], 'blank_method','pchip', 'verbose',true);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opt,fn{i}) || isempty(opt.(fn{i})), opt.(fn{i}) = def.(fn{i}); end
end

info = struct('n_blanked',0, 'samples_per_event',0, 'window_ms',opt.blank_window, 'skipped',0);

if ~isfield(EEG,'event') || isempty(EEG.event)
    if opt.verbose
        fprintf('[blank] No events; nothing to blank.\n');
    end
    return
end
if EEG.trials > 1
    warning('ieeglab_blank_stim:epoched', ...
        ['Data is already epoched. Blanking is intended for CONTINUOUS data, before ' ...
         'filtering and epoching - that ordering is the whole point. Blanking now still ' ...
         'removes the artifact but the filter has already rung on it.']);
end

X = EEG.data;
was2D = ismatrix(X);
if was2D, X = reshape(X, size(X,1), size(X,2), 1); end
[C, T, N] = size(X);

% Window in samples, relative to the event sample
w = round(double(opt.blank_window)/1000 * EEG.srate);
if w(2) <= w(1)
    error('ieeglab_blank_stim:badWindow', ...
        'blank_window must be [start stop] in ms with stop > start; got %s.', mat2str(opt.blank_window));
end
offs = w(1):w(2);
info.samples_per_event = numel(offs);
if numel(offs) < 1
    warning('ieeglab_blank_stim:windowTooNarrow', ...
        'Blank window [%g %g] ms is under one sample at %g Hz; nothing blanked.', ...
        opt.blank_window(1), opt.blank_window(2), EEG.srate);
    return
end

% Blanking only means something when the artifact is actually resolved in time.
% Below roughly 500 Hz the stimulation transient has already been smeared by the
% recording system's anti-alias filter, so there is no sharp edge left to cut
% out and no clean sample to interpolate from.
if EEG.srate < 500
    warning('ieeglab_blank_stim:lowSampleRate', ...
        ['Sampling rate is %g Hz, so the %g ms blank window covers only %d sample(s). ' ...
         'At this rate the stimulation artifact is already smeared across samples by ' ...
         'anti-alias filtering and blanking cannot cleanly remove it. Blanking is ' ...
         'intended for data at its native acquisition rate (typically 1-10 kHz); ' ...
         'apply it BEFORE any downsampling.'], ...
         EEG.srate, diff(opt.blank_window), numel(offs));
end

lat = round([EEG.event.latency]);
% For epoched data an event's epoch index says which trial it belongs to
if N > 1 && isfield(EEG.event,'epoch')
    ep = [EEG.event.epoch];
else
    ep = ones(1, numel(lat));
end

nBlank = 0; nSkip = 0;
for e = 1:numel(lat)
    tr = 1;
    base = lat(e);
    if N > 1
        % Epoched: EEGLAB latencies count across the concatenated epochs, so
        % convert to a sample within this event's own epoch.
        tr = min(max(ep(min(e, numel(ep))), 1), N);
        base = lat(e) - (tr - 1) * T;
    end
    idx = base + offs;
    idx = idx(idx >= 1 & idx <= T);
    if numel(idx) < 2
        nSkip = nSkip + 1;
        continue
    end

    if strcmpi(opt.blank_method, 'nan')
        X(:, idx, tr) = NaN;
        nBlank = nBlank + 1;
        continue
    end

    % Anchor the interpolation on up to 3 clean samples either side of the gap
    lo = idx(1) - 1;
    hi = idx(end) + 1;
    anchorIdx = [max(1,lo-2):lo, hi:min(T,hi+2)];
    anchorIdx = unique(anchorIdx(anchorIdx >= 1 & anchorIdx <= T));
    anchorIdx = setdiff(anchorIdx, idx);

    if numel(anchorIdx) >= 2
        % interp1 works column-wise, so pass [samples x channels] and transpose back
        X(:, idx, tr) = interp1(anchorIdx, X(:, anchorIdx, tr).', ...
                                idx, lower(opt.blank_method), 'extrap').';
    elseif ~isempty(anchorIdx)
        % Only one side available (event at a data edge): hold that sample
        X(:, idx, tr) = repmat(X(:, anchorIdx(1), tr), 1, numel(idx));
    else
        nSkip = nSkip + 1;
        continue
    end
    nBlank = nBlank + 1;
end

if was2D, X = reshape(X, C, T); end
EEG.data = X;
info.n_blanked = nBlank;
info.skipped   = nSkip;

EEG.ieeglab.blank = info;

if opt.verbose
    fprintf('[blank] Replaced %g ms (%d samples) around %d stimulation event(s) using %s interpolation.\n', ...
        diff(opt.blank_window), info.samples_per_event, nBlank, opt.blank_method);
    if nSkip > 0
        fprintf('[blank] %d event(s) too close to a data edge to blank.\n', nSkip);
    end
end
end
