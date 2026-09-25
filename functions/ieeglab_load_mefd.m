function [metadata, data] = ieeglab_load_mefd(sessPath, password, channels, rangeType, varargin)
% ieeglab_load_mefd() - Read selected channels and samples of a MEF3 session.
%
% Usage:
%   metadata = ieeglab_load_mefd(sessPath)                      % metadata only
%   [metadata, data] = ieeglab_load_mefd(sessPath, [], {'RA1','RA2'}, 'samples', [s0 s1])
%
% - channels: {} (all), cellstr, string, char, or numeric indices
% - rangeType: 'samples' (default) or 'time'
% - ranges as Nx2 or (rangeStart, rangeEnd), 0-based, end exclusive (as meflib)
% - data are in microvolts; a parfor is used when a pool is already open
%
% Unlike pop_MEF3, which loads the whole session, only the requested channels
% and samples are read.
%
% Requires the EEGLAB MEF3 plugin (read_mef_session_metadata, read_mef_ts_data).

metadata = [];
data = [];

narginchk(1,6);

% ---- Fast input normalization
if nargin < 2 || isempty(password),  password = []; end
if nargin < 3, channels = {}; end
if nargin < 4 || isempty(rangeType), rangeType = 'samples'; end

if isstring(sessPath), sessPath = char(sessPath); end
if ~ischar(sessPath) || ~exist(sessPath,'dir')
    error('Invalid session directory: %s', char(sessPath));
end

if ~isempty(password)
    if isstring(password), password = char(password); end
    if ~ischar(password), error('password must be [], string, or char'); end
end

rangeType = lower(rangeType);
if ~ismember(rangeType, {'time','samples'})
    error('rangeType must be ''time'' or ''samples''');
end

% ---- Read metadata once
try
    metadata = read_mef_session_metadata(sessPath, password);
catch e
    error('%s\nUnable to read MEF3 metadata', e.message);
end
if isempty(metadata) ...
        || ~isfield(metadata,'earliest_start_time') || metadata.earliest_start_time < 0 ...
        || ~isfield(metadata,'latest_end_time')     || metadata.latest_end_time < 0
    error('No valid MEF3 metadata found in ''%s''.', sessPath);
end

% Warn if no channels at all
if metadata.number_of_time_series_channels == 0 && metadata.number_of_video_channels == 0
    warning('No channels found in session directory');
end

% ---- Exit early if only metadata requested
if nargout == 1, return; end

% ---- Sort channels by acquisition_channel_number (vectorized)
if metadata.number_of_time_series_channels > 0
    acqChNum = arrayfun(@(c) c.metadata.section_2.acquisition_channel_number, ...
        metadata.time_series_channels);
    [~, ord] = sort(acqChNum);
    metadata.time_series_channels = metadata.time_series_channels(ord);

    if min(acqChNum) ~= 1
        warning('Acquisition channel numbering does not start at 1.');
    end
    if ~isempty(setdiff(min(acqChNum):max(acqChNum), acqChNum))
        warning('Acquisition channel numbering has gaps.');
    end
end

% ---- Build the list of channel indices to load (accept numeric or names)
allNames = {metadata.time_series_channels.name};
nAll = numel(allNames);

if isempty(channels)
    chIdx = 1:nAll; % fast path: all channels in sorted order
else
    if ischar(channels) || isstring(channels)
        channels = cellstr(channels);
    end
    if isnumeric(channels)
        % Convenience: allow numeric indices directly
        chIdx = channels(:).';
        if any(chIdx < 1 | chIdx > nAll)
            error('Channel indices out of bounds (1..%d).', nAll);
        end
    else
        % cellstr of names, case-insensitive match done once
        if ~iscellstr(channels), error('channels must be cellstr, string(s), char, or numeric indices'); end
        lowAll = lower(allNames);
        lowReq = lower(channels(:).');
        [tf, loc] = ismember(lowReq, lowAll);
        if ~all(tf)
            missing = channels(~tf);
            error('Requested channel(s) not found: %s', strjoin(missing, ', '));
        end
        % Preserve requested order
        chIdx = loc;
        % Check duplicates quickly
        if numel(unique(chIdx)) ~= numel(chIdx)
            error('Duplicate channels in request are not allowed.');
        end
    end
end
nCh = numel(chIdx);

% ---- Precompute full file paths for selected channels (no per-iter strcat)
chanPaths = strcat( ...
    {metadata.time_series_channels(chIdx).path}, filesep, ...
    {metadata.time_series_channels(chIdx).name}, '.', ...
    {metadata.time_series_channels(chIdx).extension} );

% ---- Normalize ranges
if numel(varargin) == 0
    ranges = [-1 -1]; % single catch-all
elseif numel(varargin) == 1
    ranges = varargin{1};
    if isempty(ranges) || ~isnumeric(ranges) || size(ranges,2) ~= 2 || size(ranges,1) < 1
        error('ranges must be an Nx2 numeric matrix (>=1 rows).');
    end
    if any(ranges(:) < 0), error('ranges values must be >= 0 (or -1 for defaults).'); end
    ranges = sort(ranges, 2, 'ascend');
    if any(ranges(:,2) - ranges(:,1) <= 0)
        error('Each range must have positive length.');
    end
elseif numel(varargin) == 2
    rangeStart = varargin{1};
    rangeEnd   = varargin{2};
    if ~isscalar(rangeStart) || ~isscalar(rangeEnd)
        error('rangeStart and rangeEnd must be scalars.');
    end
    if (rangeStart ~= -1 && (rangeStart < 0 || mod(rangeStart,1)~=0)) ...
            || (rangeEnd   ~= -1 && (rangeEnd   < 0 || mod(rangeEnd,1)~=0))
        error('rangeStart/end must be whole numbers (>=0) or -1.');
    end
    ranges = sort([rangeStart, rangeEnd], 2, 'ascend');
else
    error('Too many input arguments.');
end
nR = size(ranges,1);

% int64 once (meflib expects 64-bit)
rangeStart64 = int64(ranges(:,1));
rangeEnd64   = int64(ranges(:,2));

% ---- Fast collection strategy:
% Read each (channel, range) into a cell; capture lengths in a matrix; single pack at end
sigCells = cell(nCh, nR);
lenMat   = zeros(nCh, nR);  % double is fine; avoids parfor broadcast issues

% Use parfor only if a pool is already open (otherwise parfor overhead can hurt)
usePar = ~isempty(gcp('nocreate'));

if usePar
    parfor iCh = 1:nCh
        localLen = zeros(1, nR);
        pth = chanPaths{iCh};
        for iR = 1:nR
            s = read_mef_ts_data(pth, password, rangeType, rangeStart64(iR), rangeEnd64(iR), true)';
            sigCells{iCh, iR} = s; %#ok<PFBNS>
            localLen(iR) = numel(s);
        end
        lenMat(iCh, :) = localLen;
    end
else
    for iCh = 1:nCh
        pth = chanPaths{iCh};
        for iR = 1:nR
            s = read_mef_ts_data(pth, password, rangeType, rangeStart64(iR), rangeEnd64(iR), true)';
            sigCells{iCh, iR} = s;
            lenMat(iCh, iR)   = numel(s);
        end
    end
end

% ---- Determine max length per range, allocate once, and pack (pad with NaN)
maxLenPerRange = max(lenMat, [], 1);
if nR == 1
    data = nan(nCh, maxLenPerRange, 'like', sigCells{find(~cellfun('isempty',sigCells),1)}); % keep double
    for iCh = 1:nCh
        s = sigCells{iCh,1};
        if ~isempty(s)
            data(iCh, 1:numel(s)) = s;
        end
    end
else
    data = nan(nCh, max(maxLenPerRange), nR, 'like', sigCells{find(~cellfun('isempty',sigCells),1)});
    for iR = 1:nR
        L = maxLenPerRange(iR);
        for iCh = 1:nCh
            s = sigCells{iCh,iR};
            if ~isempty(s)
                data(iCh, 1:numel(s), iR) = s;
            end
        end
    end
end
end
