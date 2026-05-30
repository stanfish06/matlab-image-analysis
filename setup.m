%-----------------------------------------
% set up path (need update)
%-----------------------------------------

% entries of matlab search path
pathentries = regexp(path, pathsep, 'split');

% find old repo entries and remove them from path before adding this checkout
disp('removing old entries from path');
for i = 1:numel(pathentries)
    if ~isempty(regexp(pathentries{i}, 'matlab-image-analysis','once'))...
            || ~isempty(regexp(pathentries{i}, 'matlab-image-analysis','once'))...
            || ~isempty(regexp(pathentries{i}, 'matlab image analysis','once'))
        rmpath(pathentries{i});
    end
end

%%
% Path of the current script is the repository root. Use mfilename so setup
% also works in noninteractive/headless MATLAB sessions.
repopath = fileparts(mfilename('fullpath'));
if isempty(repopath)
    repopath = pwd;
end
disp(['repo path: ' repopath]);

% Add repo paths, excluding archived/dead code so old routines do not shadow
% active implementations.
repoPathEntries = regexp(genpath(repopath), pathsep, 'split');
repoPathEntries = repoPathEntries(~cellfun('isempty', repoPathEntries));
excludeRoots = {fullfile(repopath, 'archive'), fullfile(repopath, '.git')};
keep = true(size(repoPathEntries));
for i = 1:numel(excludeRoots)
    keep = keep & ~strcmp(repoPathEntries, excludeRoots{i}) ...
        & ~startsWith(repoPathEntries, [excludeRoots{i} filesep]);
end
addpath(strjoin(repoPathEntries(keep), pathsep));
disp('added active repository directories to path');

% storing the path in the startup directory
upath = userpath;
if ~isempty(upath)
    upathEntries = regexp(upath, pathsep, 'split');
    savepath(fullfile(upathEntries{1}, 'pathdef.m'));
end
%% add fiji path
setup_fiji = true;
if setup_fiji
    env_path = strsplit(getenv('PATH'), pathsep);
    match = cellfun(@(x)~isempty(regexp(x, 'Fiji', 'once')), env_path);
    if any(match)
        fiji_path = env_path{find(match, 1, 'first')};
        addpath(fullfile(fiji_path, 'scripts'));
    else
        warning('Fiji path not found in PATH; skipping Fiji scripts');
    end
end
