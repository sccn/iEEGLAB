function p = ieeglab_abspath(p)
% ieeglab_abspath() - Absolute path of an existing folder or file.
%
% Usage:
%   p = ieeglab_abspath('ds004696/sub-02/ses-ieeg01/ieeg/sub-02_ses-ieeg01_task-ccep_run-01_ieeg.mefd')
%
% A relative path is taken from MATLAB's current folder. The MEF3 reader needs
% this: on macOS and Linux, meflib makes a relative path absolute with the PWD
% environment variable, which MATLAB does not update on cd, so after a cd it
% looks in the folder MATLAB was started from and crashes MATLAB instead of
% returning an error. dir() resolves the path from MATLAB's own current folder.
%
% Cedric Cannard, iEEGLAB, 2026

if isstring(p), p = char(p); end
if ~ischar(p) || isempty(p)
    error('ieeglab_abspath:input', 'The path must be a non-empty character vector or string.');
end
% meflib does not expand '~' either: do it here (macOS and Linux)
if ~ispc && (strcmp(p, '~') || startsWith(p, '~/'))
    p = [getenv('HOME') p(2:end)];
end
if isfolder(p)
    d = dir(p);
    d = d(strcmp({d.name}, '.'));   % the '.' entry's folder is p itself, resolved
    if ~isempty(d)
        p = d(1).folder;
        return
    end
elseif isfile(p)
    d = dir(p);
    if isscalar(d)
        p = fullfile(d.folder, d.name);
        return
    end
else
    error('ieeglab_abspath:notFound', 'No such folder or file: %s', p);
end
% Fallback (a file system that lists no '.' entry): prefix the current folder
% when the path is relative. A Windows drive path ('C:\x', or 'C:x', relative to
% that drive's own current folder) is left as it is.
if ~(startsWith(p, {'/', '\'}) || ~isempty(regexp(p, '^[A-Za-z]:', 'once')))
    p = fullfile(pwd, p);
end
end
