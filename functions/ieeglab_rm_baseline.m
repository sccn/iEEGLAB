function EEG = ieeglab_rm_baseline(EEG)
% ieeglab_rm_baseline  Apply baseline correction to EEG.data using GUI choices.
%
% Uses:
%   EEG.data                : [channels x samples x trials]
%   EEG.srate               : sampling rate (Hz)
%   EEG.times               : time (ms), length must match samples
%   EEG.ieeglab.opt fields  :
%       baseline_method     : 'median' | 'mean' | 'trimmed mean'
%       baseline_period     : [start_ms end_ms]
%       baseline_mode       : 'subtract' | 'divide'
%
% Output:
%   EEG.data (baseline-corrected)

% ---------- sanity ----------
assert(isfield(EEG,'data') && ~isempty(EEG.data), 'EEG.data is empty.');
X = EEG.data;
if ndims(X)==2, X = reshape(X, size(X,1), size(X,2), 1); end
[C,T,N] = size(X);

assert(isfield(EEG,'times') && numel(EEG.times)==T, ...
    'EEG.times must exist and match the time dimension of EEG.data.');

% ---------- options from GUI ----------
opt = struct('baseline_method','median', 'baseline_period',[-500 -50], 'baseline_mode','subtract');
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isstruct(EEG.ieeglab.opt)
    f = EEG.ieeglab.opt;
    if isfield(f,'baseline_method') && ~isempty(f.baseline_method)
        opt.baseline_method = f.baseline_method;
    end
    if isfield(f,'baseline_period') && numel(f.baseline_period)>=2
        opt.baseline_period = double(f.baseline_period(1:2));
    end
    if isfield(f,'baseline_mode') && ~isempty(f.baseline_mode)
        opt.baseline_mode = f.baseline_mode;
    end
end

% Normalize names
method = lower(strtrim(string(opt.baseline_method)));
method = strrep(method,'_',' ');                 % allow 'trimmed_mean'
mode   = lower(strtrim(string(opt.baseline_mode)));
valid_methods = ["median","mean","trimmed mean"];

% The '1/f' (aperiodic) option was removed in 2026. It fitted an aperiodic
% slope to a zero-padded baseline FFT and subtracted it from the full-epoch
% magnitude spectrum. That is spectral whitening, not baseline correction, and
% it was wrong on its own terms too: the predicted spectrum was built over
% 0..fs rather than 0..fs/2, so it was not Hermitian-symmetric and
% real(ifft(...)) silently discarded a large imaginary part, while the baseline
% and epoch periodograms were on different energy scales.
% For aperiodic separation, use a method built for it (specparam/FOOOF) in a
% spectral analysis - not in an evoked-response baseline step.
% This check must precede the generic validity assert so the message is useful.
if method=="1/f" || contains(method,"aperiod")
    error('ieeglab_rm_baseline:removedMethod', ...
        ['The "1/f" baseline method has been removed: it performed spectral ' ...
         'whitening rather than baseline correction, and did so incorrectly. ' ...
         'Use "median" (default), "mean" or "trimmed mean". For aperiodic ' ...
         'separation use a dedicated spectral method such as specparam/FOOOF.']);
end

assert(any(method==valid_methods), ...
    'Invalid baseline_method "%s". Valid options: %s', method, strjoin(cellstr(valid_methods), ', '));
assert(any(mode==["subtract","divide"]), 'Invalid baseline_mode: %s', mode);

% ---------- baseline indices from EEG.times (ms) ----------
% The requested window is snapped to the nearest samples. That snapping is
% silent, so the window has to be validated against the epoch FIRST: otherwise
% a window that lies wholly outside the epoch collapses onto one edge sample,
% and a window that lies after the stimulus quietly baselines the response
% against itself. Both were previously possible.
t_ms = double(EEG.times(:));
bsl_ms = sort(double(opt.baseline_period(:))).';

if bsl_ms(2) < t_ms(1) || bsl_ms(1) > t_ms(end)
    error('ieeglab_rm_baseline:windowOutsideEpoch', ...
        ['Baseline window [%g %g] ms lies entirely outside the epoch, which spans ' ...
         '[%g %g] ms. Choose a window inside the epoch.'], ...
         bsl_ms(1), bsl_ms(2), t_ms(1), t_ms(end));
end
if bsl_ms(1) < t_ms(1) || bsl_ms(2) > t_ms(end)
    warning('ieeglab_rm_baseline:windowClipped', ...
        ['Baseline window [%g %g] ms extends past the epoch [%g %g] ms and will be ' ...
         'clipped. The baseline is therefore computed over a shorter interval than requested.'], ...
         bsl_ms(1), bsl_ms(2), t_ms(1), t_ms(end));
end
if bsl_ms(2) > 0
    warning('ieeglab_rm_baseline:windowIncludesPostStim', ...
        ['Baseline window [%g %g] ms extends past t=0 into the post-stimulus period. ' ...
         'The evoked response would then be subtracted from itself, shrinking or ' ...
         'inverting the very effect you are measuring. Use a window that ends before 0 ' ...
         '(for CCEP, before the stimulation artifact, e.g. [-500 -10] ms).'], ...
         bsl_ms(1), bsl_ms(2));
end

[~,i0] = min(abs(t_ms - bsl_ms(1)));
[~,i1] = min(abs(t_ms - bsl_ms(2)));
i0 = max(1, min(T, i0));
i1 = max(1, min(T, i1));
if i1 < i0, [i0,i1] = deal(i1,i0); end
bsl_idx = i0:i1;

% If the stimulation artifact was blanked, the samples inside that window are
% interpolated, not recorded. Averaging them into the baseline mixes a
% synthetic value into the reference level, which then shifts the whole epoch.
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'blank') && ...
        isfield(EEG.ieeglab.blank,'window_ms') && ~isempty(EEG.ieeglab.blank.window_ms)
    bw = EEG.ieeglab.blank.window_ms;
    if bsl_ms(2) > bw(1) && bsl_ms(1) < bw(2)
        warning('ieeglab_rm_baseline:overlapsBlankWindow', ...
            ['Baseline window [%g %g] ms overlaps the blanked stimulation window ' ...
             '[%g %g] ms. Those samples are interpolated, not recorded, so they would ' ...
             'contribute a synthetic value to the baseline estimate. End the baseline ' ...
             'before %g ms.'], bsl_ms(1), bsl_ms(2), bw(1), bw(2), bw(1));
    end
end

if numel(bsl_idx) < 3
    warning('ieeglab_rm_baseline:tooFewSamples', ...
        ['Baseline window [%g %g] ms covers only %d sample(s) at %g Hz. The baseline ' ...
         'estimate will be dominated by noise; widen the window.'], ...
         bsl_ms(1), bsl_ms(2), numel(bsl_idx), EEG.srate);
end
if isempty(bsl_idx)
    warning('Empty baseline index range; no correction applied.');
    EEG.data = X; if ndims(EEG.data)==2, EEG.data = squeeze(EEG.data); end
    return
end

fprintf('Baseline: method="%s", mode="%s", window=[%g %g] ms (idx %d:%d)\n', ...
    method, mode, bsl_ms(1), bsl_ms(2), i0, i1);

% ---------- apply correction ----------
switch method
    case {"median","mean","trimmed mean"}
        % compute baseline statistic per channel × epoch
        switch method
            case "median"
                bsl_vals = median(X(:,bsl_idx,:), 2);   % [C x 1 x N]
            case "mean"
                bsl_vals = mean(X(:,bsl_idx,:), 2);
            case "trimmed mean"
                trim_pct = 20; % %
                trim_amt = round((trim_pct/100) * numel(bsl_idx));
                bsl_vals = nan(C,1,N);
                for ch = 1:C
                    for ep = 1:N
                        vals = squeeze(X(ch,bsl_idx,ep));
                        vals = sort(vals);
                        if numel(vals) > 2*trim_amt
                            vals = vals(trim_amt+1:end-trim_amt);
                        end
                        bsl_vals(ch,1,ep) = mean(vals,'omitnan');
                    end
                end
        end

        switch mode
            case "subtract"
                % The safe operation. Subtraction assumes the baseline is an
                % ADDITIVE offset, which is what a DC/drift offset in a
                % broadband time-domain signal actually is.
                X = X - repmat(bsl_vals, [1 T 1]);

            case "divide"
                % Divisive baselining assumes a MULTIPLICATIVE model: that
                % everything in the epoch scales by a common factor relative to
                % baseline. Gyurkovics et al. (2021, NeuroImage 237:118192)
                % show that when signal and background are instead independent
                % ADDITIVE contributors, divisive/dB baselining misattributes
                % broadband (1/f) changes to narrowband effects. Their analysis
                % is about time-frequency POWER, where the quantity is
                % positive and a ratio is at least well defined.
                %
                % In the time domain it is worse: the signal is signed, and
                % after a high-pass filter the pre-stimulus mean is near zero
                % by construction, so the ratio explodes and its sign is
                % arbitrary. Refuse unless the baseline is large enough for the
                % ratio to mean anything.
                sigSD  = std(X(:), 'omitnan');
                bslMag = median(abs(bsl_vals(:)), 'omitnan');
                if ~isfinite(bslMag) || bslMag < 0.1*sigSD
                    error('ieeglab_rm_baseline:unsafeDivide', ...
                        ['Divisive baseline correction refused: the typical baseline ' ...
                         'magnitude (%.3g) is small relative to the signal (SD %.3g), ' ...
                         'so dividing by it produces exploding and sign-unstable values.\n' ...
                         'This is expected for high-pass filtered time-domain data, ' ...
                         'where the pre-stimulus mean is near zero by construction.\n' ...
                         'Use baseline_mode = "subtract". Divisive and dB baselining ' ...
                         'belong to time-frequency POWER analysis, and even there see ' ...
                         'Gyurkovics et al. (2021) NeuroImage 237:118192 on how they ' ...
                         'misattribute broadband 1/f changes to narrowband effects.'], ...
                         bslMag, sigSD);
                end
                warning('ieeglab_rm_baseline:divisiveBaseline', ...
                    ['Using divisive baseline correction on time-domain data. This ' ...
                     'assumes a multiplicative model; see Gyurkovics et al. (2021) ' ...
                     'NeuroImage 237:118192. "subtract" is the appropriate choice for ' ...
                     'evoked potentials.']);
                denom = max(eps, abs(bsl_vals));
                X = bsxfun(@rdivide, X, denom);
        end
end

% ---------- write back ----------
EEG.data = X;
end
