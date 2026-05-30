%%
clear; close all; clc
scriptPath = pwd;
baseDir = scriptPath;
%%
scatterDir = 'scatter';
if ~exist(scatterDir, 'dir')
    mkdir(scatterDir);
end
%%
% load(fullfile(baseDir,'positions_colony.mat'),'positions')
% load(fullfile(baseDir,'meta_colony.mat'),'meta')
% positions = cellfun(@(x)x, positions);
load(fullfile(baseDir,'positions.mat'),'positions')
load(fullfile(baseDir,'meta.mat'),'meta')
%%
stats = cellStats(positions, meta, positions(1).dataChannels);
stats.markerChannels = [1 2 3 4];
% automaticly get thresholds
% confidence = 0.95;
% conditions = 1:numel(meta.conditions);
% whichthreshold = []; %[1 1 2]
% stats.getThresholds(confidence, conditions, whichthreshold);
%%
stats.thresholds(1) = 100;
stats.thresholds(2) = 500;
stats.thresholds(3) = 600;
stats.thresholds(4) = 400;
%%
dataDir = fullfile(baseDir, scatterDir);
dat = [];
% 1:numel(meta.conditions)
for condi = 1:numel(meta.conditions)
    close all;
    options = struct();

    options.selectDisType = [1, 1, 1, 1]; % 1 - nuc, 2 - cyto, 3 - NCratio, 4 - CNratio
    options.channelThresholds = stats.thresholds;
    options.showThreshold = [false false false false];
    options.minimal = false;
    options.edgeDistance = false;
    options.clim = struct('TFAP2C',[0 3],'EOMES',[0 3],'SOX17',[0 3]);
    options.imageType = 'Disordered';
    options.conditionIdx = condi;
    options.channelCombos = {[4 3 2],[4 2 3],[2 3 4]};
    options.titleChangeLine = false;

    options.axisLabel = meta.channelLabel;
    options.channelMax = exp([8 3 3 3])-1;
    options.log1p = [true true true true];
    options.conditionsCombined = false;
    options.showThreshold = [true true true true];
    options.crossConditionColonySummary = true;
    options.savePercent = true;

    result = scatterMicropattern(stats, meta, dataDir, options);
    dat = [dat; result.percentData];
end
%
writetable(dat, fullfile(dataDir, 'scatterPercent.csv'))
