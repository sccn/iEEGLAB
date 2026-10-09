function tf = ieeglab_is_mex(name)
% ieeglab_is_mex() - True when calling 'name' would run a compiled (MEX) file
%                    for this platform.
%
% Usage:
%   tf = ieeglab_is_mex('read_mef_ts_data')
%
% Plugins that ship MEX files often ship help-only .m files of the same name;
% on a platform without the MEX file (the MEF3 plugin on Apple Silicon), the .m
% file is what resolves. exist() is not used: right after a folder is added to
% the path it can answer 2 (a file) for a function whose MEX file is the one
% that runs. which() names the file that will actually be called.
%
% Cedric Cannard, iEEGLAB, 2026

w = which(name);
tf = ~isempty(w) && endsWith(lower(w), ['.' lower(mexext)]);
end
