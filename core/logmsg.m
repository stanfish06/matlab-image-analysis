function logmsg(varargin)
%LOGMSG Print a timestamped progress/log message to the console.
%   LOGMSG(FMT, ...) prints "[HH:MM:SS] " followed by SPRINTF(FMT, ...).
%   LOGMSG(LEVEL, FMT, ...) prepends a level tag when LEVEL is one of
%   'INFO', 'WARN', 'ERROR' (case-insensitive), e.g. "[12:00:00] [WARN] ...".
%
%   Intended for logging critical steps so pipelines are debuggable. Safe to
%   call inside PARFOR (it writes to the client via FPRINTF); include an
%   iteration index in the message so interleaved worker output stays readable,
%   e.g. LOGMSG('position %d/%d done', i, n).
%
%   Logging only: it never changes computation and never throws (a bad format
%   string falls back to printing the raw text).
%
%   Examples:
%       logmsg('stitching montage %d/%d', m, nM);
%       logmsg('WARN', 'position %d has no cells', p);

    levels = {'INFO','WARN','ERROR'};
    if nargin >= 2 && (ischar(varargin{1}) || isstring(varargin{1})) ...
            && any(strcmpi(char(varargin{1}), levels))
        tag = sprintf(' [%s]', upper(char(varargin{1})));
        fmt = varargin{2};
        args = varargin(3:end);
    elseif nargin >= 1
        tag = '';
        fmt = varargin{1};
        args = varargin(2:end);
    else
        return
    end

    try
        msg = sprintf(char(fmt), args{:});
    catch
        msg = char(fmt); % fall back to the raw format string
    end
    ts = char(datetime('now', 'Format', 'HH:mm:ss'));
    fprintf('[%s]%s %s\n', ts, tag, msg);
end
