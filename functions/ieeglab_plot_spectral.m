function ieeglab_plot_spectral(EEG, S, varargin)
% ieeglab_plot_spectral() - Draw the result of ieeglab_spectral.
%
% Usage:
%   ieeglab_plot_spectral(EEG, S)
%   ieeglab_plot_spectral(EEG, S, 'channel', 'RA3')       % one contact's time-frequency map
%   ieeglab_plot_spectral(EEG, S, 'brain', false)
%
% Spectrum: the spectrum of every contact (log power), their median in black,
% bad contacts dashed, the band shaded.
% Time-frequency: the normalized time-frequency map (one contact, or the
% average of the good contacts) and the band power of every contact over time;
% then, when the contacts have coordinates, the band value of each contact on
% the brain (ieeglab_topoplot).
%
% Options:
%   'channel'  contact label or index for the time-frequency map. Default: the
%              average of the good contacts
%   'clim'     colour limits of the time-frequency map
%   'brain'    draw the band value on the brain. Default true
%
% Cedric Cannard, iEEGLAB, 2026

p = inputParser;
p.addParameter('channel', []);
p.addParameter('clim', []);
p.addParameter('brain', true);
p.parse(varargin{:});
o = p.Results;
good = ~S.bad & S.n_trials > 0;
cond = '';
if ~isempty(S.condition), cond = sprintf(' - %s', S.condition); end
bandTxt = '';
if ~isempty(S.band), bandTxt = sprintf('%g-%g Hz', S.band(1), S.band(2)); end

if strcmp(S.type, 'spectrum')
    ieeglab_figure('Name', 'iEEG power spectra', 'Position', [80 80 820 520]);
    ax = axes; hold(ax, 'on');
    if ~isempty(S.band)
        yl = [min(S.power(:), [], 'omitnan') max(S.power(:), [], 'omitnan')];
        if all(isfinite(yl)) && yl(1) > 0
            patch(ax, S.band([1 2 2 1]), yl([1 1 2 2]), [1 0.9 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
    end
    plot(ax, S.freqs, S.power(good, :)', 'Color', [0.55 0.65 0.8]);
    if any(S.bad), plot(ax, S.freqs, S.power(S.bad, :)', '--', 'Color', [0.85 0.35 0.3]); end
    plot(ax, S.freqs, median(S.power(good, :), 1, 'omitnan'), 'k', 'LineWidth', 2);
    set(ax, 'YScale', 'log'); box(ax, 'on'); grid(ax, 'on');
    xlabel(ax, 'Frequency (Hz)', 'FontWeight', 'bold'); ylabel(ax, 'Power (\muV^2/Hz)', 'FontWeight', 'bold');
    title(ax, sprintf('Power spectrum of %d channels (median in black)%s', nnz(good), cond), 'FontWeight', 'normal');
    xlim(ax, [S.freqs(1) S.freqs(end)]);
    what = 'Band power (\muV^2)';
else
    switch S.norm
        case 'db',        unit = 'Power (dB vs baseline)';
        case 'relchange', unit = 'Power (relative change vs baseline)';
        case 'ratio',     unit = 'Power (ratio to baseline)';
        otherwise,        unit = 'Power (\muV^2)';
    end
    if isempty(o.channel)
        M = reshape(mean(S.tf(good, :, :), 1, 'omitnan'), numel(S.freqs), numel(S.times));
        who = sprintf('average of %d channels', nnz(good));
    else
        c = o.channel;
        if ~isnumeric(c), c = find(strcmpi(S.labels, char(c)), 1); end
        if isempty(c), error('ieeglab_plot_spectral:noChannel', 'No contact %s.', char(string(o.channel))); end
        M = reshape(S.tf(c, :, :), numel(S.freqs), numel(S.times));
        who = S.labels{c};
    end
    ieeglab_figure('Name', 'iEEG time-frequency', 'Position', [80 80 1100 480]);
    ax1 = subplot(1, 2, 1);
    imagesc(ax1, S.times, S.freqs, M, 'AlphaData', isfinite(M));
    axis(ax1, 'xy'); set(ax1, 'Color', [0.9 0.9 0.9]);
    cl = o.clim;
    if isempty(cl)
        m = max(abs(M(:)), [], 'omitnan');
        if ~strcmp(S.norm, 'none') && ~strcmp(S.norm, 'ratio') && isfinite(m) && m > 0, cl = [-m m]; end
    end
    if ~isempty(cl), caxis(ax1, cl); end %#ok<CAXIS>
    local_zero_line(ax1);
    cb = colorbar(ax1); cb.Label.String = unit;
    xlabel(ax1, 'Time (ms)', 'FontWeight', 'bold'); ylabel(ax1, 'Frequency (Hz)', 'FontWeight', 'bold');
    title(ax1, sprintf('%s%s', who, cond), 'FontWeight', 'normal');

    ax2 = subplot(1, 2, 2);
    B = S.band_tc;
    imagesc(ax2, S.times, 1:size(B, 1), B, 'AlphaData', isfinite(B));
    set(ax2, 'Color', [0.9 0.9 0.9]);
    m = max(abs(B(:)), [], 'omitnan');
    if ~strcmp(S.norm, 'none') && ~strcmp(S.norm, 'ratio') && isfinite(m) && m > 0, caxis(ax2, [-m m]); end %#ok<CAXIS>
    local_zero_line(ax2);
    if size(B, 1) <= 60, set(ax2, 'YTick', 1:size(B, 1), 'YTickLabel', S.labels, 'FontSize', 7); end
    cb = colorbar(ax2); cb.Label.String = unit;
    xlabel(ax2, 'Time (ms)', 'FontWeight', 'bold'); ylabel(ax2, 'Channel', 'FontWeight', 'bold');
    title(ax2, sprintf('Band %s, every channel', bandTxt), 'FontWeight', 'normal');
    try, colormap(ax1, parula); colormap(ax2, parula); catch, end
    what = sprintf('%s, %g to %g ms', unit, S.band_window(1), S.band_window(2));
end

% ---------------- band value on the brain ----------------
if o.brain && any(isfinite(S.band_value)) && isfield(EEG.chanlocs, 'X') && EEG.nbchan == numel(S.band_value) ...
        && any(arrayfun(@(c) ~isempty(c.X), EEG.chanlocs))
    v = S.band_value;
    if strcmp(S.type, 'spectrum'), v = log10(v); what = 'log10 band power (\muV^2)'; end
    ieeglab_topoplot(EEG, v, 'title', sprintf('Band %s: %s%s', bandTxt, what, cond));
end
end

function local_zero_line(ax)
hold(ax, 'on');
yl = ylim(ax);
plot(ax, [0 0], yl, 'k--', 'LineWidth', 1);
hold(ax, 'off');
end
