function make_scatter(baseDir)
%MAKE_SCATTER Generate scatter plots and population summaries.

if nargin < 1 || isempty(baseDir)
    baseDir = pwd;
end
originalDir = pwd;
cleanupObj = onCleanup(@() cd(originalDir));
cd(baseDir);
close all;

%%
scriptPath = pwd;
baseDir = scriptPath;
%%
scatterDir = 'scatter';
if ~exist(scatterDir, 'dir')
    mkdir(scatterDir);
end
%%
if exist(fullfile(baseDir,'positions_colony.mat'), 'file')
    load(fullfile(baseDir,'positions_colony.mat'),'positions')
    load(fullfile(baseDir,'meta_colony.mat'),'meta')
    positions = cellfun(@(x)x, positions);
else
    load(fullfile(baseDir,'positions.mat'),'positions')
    load(fullfile(baseDir,'meta.mat'),'meta')
end
%%
nch = length(meta.channelLabel);
stats = cellStats(positions, meta, positions(1).dataChannels);
stats.markerChannels = 1:nch;
% automaticly get thresholds
% confidence = 0.95;
% conditions = 1:numel(meta.conditions);
% whichthreshold = []; %[1 1 2]
% stats.getThresholds(confidence, conditions, whichthreshold);
%%
% set thresholds per channel - adjust for your stain panel
% example ENDO (DAPI,TFAP2C,EOMES,SOX17): thresholds = [100, 250, 300, 150];
% example LM (DAPI,TFAP2C,NANOG,HAND1): thresholds = [100, 250, 200, 200];
% example PM (DAPI,TFAP2C,NANOG,TBX6): thresholds = [100, 250, 200, 200];
% example ECTO (DAPI,TFAP2C,SOX2,NANOG): thresholds = [100, 250, 200, 200];
thresholds = [100, 250, 300, 150];
for ci = 1:nch
    stats.thresholds(ci) = thresholds(ci);
end
%%
dataDir = fullfile(baseDir, scatterDir);
dat = [];

% auto-generate all pairwise channel combos (excluding DAPI ch1)
% to set manually: options.channelCombos = {[4 3 2],[2 3 4],[4 2 3]}; % example ENDO combos
marker_chs = 2:nch;
channelCombos = {};
for ci = 1:length(marker_chs)
    for cj = ci+1:length(marker_chs)
        remaining = setdiff(marker_chs, [marker_chs(ci) marker_chs(cj)]);
        channelCombos{end+1} = [marker_chs(ci) marker_chs(cj) remaining];
    end
end

% build clim struct from channel labels (default [0 3] for all marker channels)
% to set manually: clim_struct = struct('TFAP2C',[0 3],'EOMES',[0 3],'SOX17',[0 3]); % example ENDO
clim_struct = struct();
for ci = 2:nch
    clim_struct.(meta.channelLabel{ci}) = [0 3];
end

% 1:numel(meta.conditions)
for condi = 1:numel(meta.conditions)
    close all;
    options = struct();

    options.selectDisType = ones(1, nch); % 1 - nuc, 2 - cyto, 3 - NCratio, 4 - CNratio
    options.channelThresholds = stats.thresholds;
    options.showThreshold = false(1, nch);
    options.minimal = false;
    options.edgeDistance = false;
    options.clim = clim_struct;
    options.imageType = 'Disordered';
    options.conditionIdx = condi;
    options.channelCombos = channelCombos;
    options.titleChangeLine = false;

    options.axisLabel = meta.channelLabel;
    options.channelMax = exp([8 repelem(3, nch-1)])-1;
    options.log1p = true(1, nch);
    options.conditionsCombined = false;
    options.showThreshold = true(1, nch);
    options.crossConditionColonySummary = true;
    options.savePercent = true;
    options.fs = 18;
    options.pfs = 25;

    result = scatterMicropattern(stats, meta, dataDir, options);
    dat = [dat; result.percentData];
end
%
writetable(dat, fullfile(dataDir, 'scatterPercent.csv'))

end
