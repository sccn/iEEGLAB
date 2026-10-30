function [EEG, com] = ieeglab_import(file, varargin)
% ieeglab_import() - Read a raw iEEG recording of a BIDS dataset into EEGLAB,
%                    using as little memory as possible.
%
% Usage:
%   EEG = ieeglab_import(file)
%   EEG = ieeglab_import(file, 'sites', {'RA2-RA3','RA3-RA4'})
%   EEG = ieeglab_import(file, 'types', 'all', 'time_range', [600 900])
%   info = ieeglab_import(file, 'info', true)      % what the file holds, no data read
%
% file   the recording: a MEF3 session folder (*.mefd), a BrainVision header
%        (*.vhdr), or an EDF, BDF or EEGLAB .set file.
%
% Options (name/value):
%   'types'       channel types to read, from the 'type' column of the
%                 recording's channels.tsv. Default {'SEEG','ECOG','DBS'}: the
%                 brain channels, without scalp EEG, ECG, EMG or trigger
%                 channels. 'all' reads every channel.
%   'channels'    channel labels to read; overrides 'types'.
%   'sites'       stimulation sites (electrical_stimulation_site of events.tsv):
%                 only the stretch from the first to the last pulse at these
%                 sites is read, with 'pad' seconds around it. A site is a
%                 pair of contacts whatever the polarity: 'RA2-RA3' and
%                 'RA3-RA2' are the same site.
%   'time_range'  [start stop] in seconds; overrides 'sites'.
%   'pad'         seconds read before the first and after the last pulse of
%                 'sites'. Default [3 3].
%   'info'        true returns a description of the file instead of the data:
%                 .labels, .types, .srate, .pnts, .sites, .site_counts,
%                 .channels_tsv, .events_tsv.
%   'verbose'     default true. Reports progress in the command window, with the
%                 time elapsed and an estimate of the time left: every 16 channels
%                 for MEF3, every 2 minutes of recording (at most 10 steps) for
%                 BrainVision, EDF and BDF.
%   'waitbar'     also show a progress window. Default false (true from the
%                 iEEGLAB > Import raw recording menu).
%
% Why not File > Import data > From BIDS folder structure: EEGLAB reads the
% whole recording, every channel, as double precision, then saves a copy. A
% long MEF3 session (OpenNeuro ds004696 sub-02: 256 channels, 82 minutes) takes
% 16 minutes and 21 GB that way, and the later steps need several times that.
% Here the data are read in single precision, a few channels at a time, only
% the channels and the time that are asked for, and nothing is written to disk.
%
% The dataset is returned with its channel labels and types only. Coordinates,
% events and channel status come from the sidecar files in the next step
% (iEEGLAB > Load electrode coordinates and events, ieeglab_load), which finds
% them next to the recording and shifts the event times when only a stretch
% of the recording was read.
%
% Cedric Cannard, iEEGLAB, 2026

p = inputParser;
p.addParameter('types', {'SEEG','ECOG','DBS'});
p.addParameter('channels', {});
p.addParameter('sites', {});
p.addParameter('time_range', []);
p.addParameter('pad', [3 3]);
p.addParameter('info', false);
p.addParameter('verbose', true);
p.addParameter('waitbar', false);
p.parse(varargin{:});
o = p.Results;
com = '';

file = char(file);
file = regexprep(file, '[\\/]+$', '');
[folder, name, ext] = fileparts(file);
if isempty(folder), folder = pwd; file = fullfile(folder, [name ext]); end
ext = lower(ext);
if ~(isfolder(file) || exist(file, 'file') == 2)
    error('ieeglab_import:notFound', 'Recording not found: %s', file);
end

% ---------- what the file holds ----------
switch ext
    case '.mefd'
        md = ieeglab_load_mefd(file);
        fileLabels = {md.time_series_channels.name};
        s2 = md.time_series_channels(1).metadata.section_2;
        srate = double(s2.sampling_frequency);
        pnts  = double(s2.number_of_samples);
    case '.vhdr'
        local_need('pop_loadbv', 'bva-io');
        [~, hdr] = evalc('pop_loadbv(folder, [name ext], [], [], true)');
        fileLabels = {hdr.chanlocs.labels}; srate = hdr.srate; pnts = hdr.pnts;
    case {'.edf', '.bdf'}
        local_need('pop_biosig', 'Biosig');
        [~, hdr] = evalc('pop_biosig(file, ''importevent'', ''off'', ''blockrange'', [0 1])');
        fileLabels = {hdr.chanlocs.labels}; srate = hdr.srate;
        h = sopen(file, 'r', 0, 'OVERFLOWDETECTION:OFF'); pnts = double(h.NRec) * double(h.SPR); sclose(h);
    case '.set'
        hdr = pop_loadset('filename', [name ext], 'filepath', folder, 'loadmode', 'info', 'verbose', 'off');
        fileLabels = {hdr.chanlocs.labels}; srate = hdr.srate; pnts = hdr.pnts;
    otherwise
        error('ieeglab_import:format', ...
            'Unknown recording format "%s". Give a .mefd folder, or a .vhdr, .edf, .bdf or .set file.', ext);
end
fileLabels = strtrim(fileLabels(:)');

% ---------- sidecars: channel types and stimulation sites ----------
stub = struct('filepath', folder, 'filename', [name ext]);
chTsv = ieeglab_bids_sibling(stub, 'channels');
evTsv = ieeglab_bids_sibling(stub, 'events');
types = repmat({''}, size(fileLabels));
order = 1:numel(fileLabels);                 % order in which channels are returned
if ~isempty(chTsv)
    ct = readtable(chTsv, 'FileType','text', 'Delimiter','\t', 'TextType','char');
    if all(ismember({'name','type'}, ct.Properties.VariableNames))
        tsvNames = strtrim(cellstr(ct.name))';
        [tf, loc] = ismember(lower(tsvNames), lower(fileLabels));
        order = [loc(tf), setdiff(1:numel(fileLabels), loc(tf), 'stable')];   % channels.tsv order first
        ty = upper(strtrim(cellstr(ct.type)))';
        types(loc(tf)) = ty(tf);
    end
end
labels = fileLabels(order);
types  = types(order);

sites = {}; siteOnset = []; evOnset = [];
if ~isempty(evTsv)
    et = readtable(evTsv, 'FileType','text', 'Delimiter','\t', 'TextType','char');
    if all(ismember({'onset','electrical_stimulation_site'}, et.Properties.VariableNames))
        evOnset = et.onset;
        if ~isnumeric(evOnset), evOnset = str2double(cellstr(evOnset)); end
        siteOnset = strtrim(cellstr(et.electrical_stimulation_site));
        ok = isfinite(evOnset) & ~ismember(lower(siteOnset), {'', 'n/a'});
        evOnset = evOnset(ok); siteOnset = siteOnset(ok);
        % one name per stimulated pair, whatever the polarity (as everywhere in iEEGLAB)
        siteOnset = cellstr(ieeglab_canonical_site(siteOnset, labels))';
        sites = unique(siteOnset, 'stable');
    end
end

if o.info
    EEG = struct('file', file, 'labels', {labels}, 'types', {types}, 'srate', srate, 'pnts', pnts, ...
        'sites', {sites(:)'}, 'site_counts', cellfun(@(s) nnz(strcmp(siteOnset, s)), sites(:)'), ...
        'channels_tsv', chTsv, 'events_tsv', evTsv);
    return
end

% ---------- channels to read ----------
if ~isempty(o.channels)
    want = cellstr(string(o.channels));
    [tf, keep] = ismember(lower(want), lower(labels));
    if ~all(tf)
        error('ieeglab_import:noChannel', 'Not in the recording: %s', strjoin(want(~tf), ', '));
    end
elseif (ischar(o.types) || isstring(o.types)) && strcmpi(char(o.types), 'all')
    keep = 1:numel(labels);
elseif all(cellfun(@isempty, types))
    keep = 1:numel(labels);
    if o.verbose
        fprintf('[import] No channels.tsv with channel types next to the recording: reading all %d channels.\n', numel(labels));
    end
else
    keep = find(ismember(types, upper(cellstr(string(o.types)))));
    if isempty(keep)
        u = unique(types(~cellfun(@isempty, types)));
        error('ieeglab_import:noType', 'No channel of type %s. Types in channels.tsv: %s', ...
            strjoin(cellstr(string(o.types)), '/'), strjoin(u, ', '));
    end
end

% ---------- time to read: samples s0 (0-based, first read) to s1 (one past the last) ----------
s0 = 0; s1 = pnts;
if ~isempty(o.time_range)
    s0 = max(0, floor(o.time_range(1) * srate)); s1 = min(pnts, ceil(o.time_range(2) * srate));
elseif ~isempty(o.sites)
    if isempty(sites)
        error('ieeglab_import:noSites', 'No events.tsv with an electrical_stimulation_site column next to the recording.');
    end
    want = cellstr(ieeglab_canonical_site(cellstr(string(o.sites)), labels));
    unknown = setdiff(want, sites);
    if ~isempty(unknown)
        error('ieeglab_import:unknownSite', 'Site(s) not in events.tsv: %s. Sites: %s', ...
            strjoin(unknown, ', '), strjoin(sites, ', '));
    end
    t = evOnset(ismember(siteOnset, want));
    s0 = max(0, floor((min(t) - o.pad(1)) * srate)); s1 = min(pnts, ceil((max(t) + o.pad(2)) * srate));
end
if s1 <= s0, error('ieeglab_import:emptyRange', 'The time range holds no sample.'); end
n = s1 - s0;
if o.verbose
    fprintf('[import] %s\n[import] %d of %d channels, %.1f of %.1f minutes at %g Hz: %.1f GB in memory (single precision).\n', ...
        [name ext], numel(keep), numel(labels), n/srate/60, pnts/srate/60, srate, numel(keep) * n * 4 / 2^30);
end

% ---------- read ----------
t0 = tic;
data = zeros(numel(keep), n, 'single');
event = [];
wb = [];
if o.waitbar
    wb = waitbar(0, sprintf('Reading %d channels of %s...', numel(keep), [name ext]), 'Name', 'iEEGLAB: import');
    set(findall(wb, 'Type', 'text'), 'Interpreter', 'none');
end
closeWb = onCleanup(@() local_close(wb));
switch ext
    case '.mefd'
        step = 16;                             % channels per read: bounds the double-precision buffer
        for k = 1:step:numel(keep)
            idx = k:min(k + step - 1, numel(keep));
            [~, X] = ieeglab_load_mefd(file, [], labels(keep(idx)), 'samples', [s0 s1]);
            if size(X, 2) ~= n
                error('ieeglab_import:short', 'Read %d samples, expected %d.', size(X, 2), n);
            end
            data(idx, :) = single(X);
            el = toc(t0); left = el / idx(end) * (numel(keep) - idx(end));
            msg = sprintf('%d of %d channels read, %s elapsed, about %s left', ...
                idx(end), numel(keep), local_mmss(el), local_mmss(left));
            if o.verbose, fprintf('[import] %s\n', msg); end
            if ~isempty(wb) && isgraphics(wb), waitbar(idx(end) / numel(keep), wb, msg); end
        end
        clear X
    case {'.vhdr', '.edf', '.bdf'}
        % Read in blocks of time (about 2 minutes each, at most 10) so that the
        % progress can be reported; the markers of each block are shifted back.
        nb = max(1, min(10, ceil(n / (120 * srate))));
        edges = round(linspace(0, n, nb + 1));
        for b = 1:nb
            a0 = s0 + edges(b); a1 = s0 + edges(b+1); m = a1 - a0;
            if strcmp(ext, '.vhdr')
                % every channel in one contiguous read, then the ones kept: asked for a
                % subset, pop_loadbv reads multiplexed files one channel at a time
                % (a strided pass over the block per channel; 88 passes on ds004080)
                [~, E] = evalc('pop_loadbv(folder, [name ext], [a0+1 a1])');
                E.data = E.data(order(keep), :);
            else
                [~, E] = evalc('pop_biosig(file, ''channels'', order(keep), ''blockrange'', [a0 a1] / srate)');
            end
            if size(E.data, 2) < m
                error('ieeglab_import:short', 'Read %d samples of block %d, expected %d.', size(E.data, 2), b, m);
            end
            data(:, edges(b)+1:edges(b+1)) = single(E.data(:, 1:m));
            if ~isempty(E.event)
                ev = E.event;
                for i = 1:numel(ev), ev(i).latency = ev(i).latency + edges(b); end
                if isempty(event), event = ev; else, event = [event(:)' ev(:)']; end %#ok<AGROW>
            end
            clear E
            el = toc(t0); left = el / edges(b+1) * (n - edges(b+1));
            msg = sprintf('%.0f%% read (%.1f of %.1f minutes), %s elapsed, about %s left', ...
                100 * edges(b+1) / n, edges(b+1) / srate / 60, n / srate / 60, local_mmss(el), local_mmss(left));
            if o.verbose, fprintf('[import] %s\n', msg); end
            if ~isempty(wb) && isgraphics(wb), waitbar(edges(b+1) / n, wb, msg); end
        end
    case '.set'
        E = pop_loadset('filename', [name ext], 'filepath', folder, 'verbose', 'off');
        data(:) = single(E.data(order(keep), s0+1:s1)); clear E
end
nBad = nnz(~isfinite(data));
if nBad > 0
    warning('ieeglab_import:nonFinite', '%d sample(s) are NaN or Inf (gaps in the recording?).', nBad);
end

% ---------- dataset ----------
EEG = eeg_emptyset;
EEG.setname = regexprep(name, '_ieeg$', '');
EEG.data    = data; clear data
EEG.srate   = srate;
EEG.nbchan  = numel(keep);
EEG.pnts    = n;
EEG.trials  = 1;
EEG.xmin    = 0;
EEG.chanlocs = struct('labels', labels(keep), 'type', types(keep));
tok = regexp(name, '^(sub-[A-Za-z0-9]+)', 'tokens', 'once');
if ~isempty(tok), EEG.subject = tok{1}; end
if isempty(evTsv) && ~isempty(event), EEG.event = event; end   % markers of the file, when no events.tsv
EEG.BIDS.sourcefile = file;                % where ieeglab_load looks for the sidecar files
EEG.etc.ieeglab_import = struct('source', file, 'first_sample', s0, 'last_sample', s1 - 1, ...
    'onset_shift_s', s0 / srate, 'channels_in_file', numel(labels), 'pnts_in_file', pnts);
EEG = eeg_checkset(EEG);
if o.verbose
    fprintf('[import] Done in %.0f s. Next: iEEGLAB > Load electrode coordinates and events.\n', toc(t0));
end

args = {};
if ~isempty(o.channels), args = [args {'channels', cellstr(string(o.channels))}];
elseif ~isequal(upper(cellstr(string(o.types))), {'SEEG','ECOG','DBS'}), args = [args {'types', o.types}]; end
if ~isempty(o.time_range), args = [args {'time_range', o.time_range}];
elseif ~isempty(o.sites), args = [args {'sites', cellstr(string(o.sites))}];
    if ~isequal(o.pad, [3 3]), args = [args {'pad', o.pad}]; end
end
com = sprintf('EEG = ieeglab_import(%s);', strjoin(cellfun(@local_lit, [{file} args], 'UniformOutput', false), ', '));
end

function s = local_lit(x)
% MATLAB literal of an argument, for the history line (quotes doubled)
if ischar(x) || (isstring(x) && isscalar(x))
    s = ['''' strrep(char(x), '''', '''''') ''''];
elseif iscell(x) || isstring(x)
    q = cellfun(@(v) ['''' strrep(char(v), '''', '''''') ''''], cellstr(x), 'UniformOutput', false);
    s = ['{' strjoin(q, ',') '}'];
else
    s = mat2str(x);
end
end

function s = local_mmss(t)
% seconds as m:ss
t = max(0, round(t));
s = sprintf('%d:%02d', floor(t / 60), mod(t, 60));
end

function local_close(wb)
if ~isempty(wb) && isgraphics(wb), delete(wb); end
end

function local_need(fcn, plugin)
if isempty(which(fcn))
    error('ieeglab_import:missingReader', ...
        'This format needs the %s plugin (File > Manage EEGLAB extensions).', plugin);
end
end
