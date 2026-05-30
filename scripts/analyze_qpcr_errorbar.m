% ANALYZE_QPCR_ERRORBAR  Analyze qPCR data and plot means with error bars [template].
%   Reads a qPCR export (.txt), groups wells by condition labels, and plots
%   per-condition means with error bars. Edit the input filename(s) and
%   conditionLabels near the top.
clear all; close all;

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% data location, modify if not the same as the location of this script
dataDir = scriptPath; 
cd(dataDir);

filenames = {'20240501 FGFKO shRNA FGFs.txt'};

% filenamesNorm = {'170115 A New TimeSeries ATP5 Oct4 Nanog Sox2_dataX.txt',...
%     '170115 A New TimeSeries ATP5O Oct4 Nanog Sox2 B_data.txt'};

% Conditions: (samples)
conditionLabels = {
'ESI WT GAPDH','','','ESI WT F4','','','ESI WT F8','','','ESI WT F17','','',...
'ESI F4KO GAPDH','','','ESI F4KO F4','','','ESI F4KO F8','','','ESI F4KO F17''','',...
'ESI F17KO GAPDH','','','ESI F17KO F4','','','ESI F17KO F8','','','ESI F17KO F17','','',...
'ESI shScrmbl GAPDH','','','ESI shScrmbl F4','','','ESI shScrmbl F8','','','ESI shScrmbl F17','','',...
'ESI shF4 GAPDH','','','ESI shF4 F4','','','ESI shF4 F8','','','ESI shF4 F17','','',...
'ESI shF4 F4newPrimer','ESI shScrmbl F4newPrimer','ESI F17KO F4newPrimer','ESI F4KO F4newPrimer','ESI WT F4newPrimer'...
'ESI shF4 F4newPrimer','ESI shScrmbl F4newPrimer','ESI F17KO F4newPrimer','ESI F4KO F4newPrimer','ESI WT F4newPrimer'...
'ESI shF4 F4newPrimer','ESI shScrmbl F4newPrimer','ESI F17KO F4newPrimer','ESI F4KO F4newPrimer','ESI WT F4newPrimer'...
};

normName = 'GAPDH';

%% read data

tolerance = 1;
data = combineQPCR(dataDir, filenames, tolerance);

Ntargets = data.Ntargets;
Nsamples = data.Nsamples;
targets = data.targets;
CTmean = data.CTmean;
CTstd = data.CTstd;

% normalization
refIdx = find(strcmp(targets,normName));
if isempty(refIdx)
    disp('using normalization from other file');
	refData =  combineQPCR(dataDir, filenamesNorm);
    refIdx = find(strcmp(refData.targets,normName));
    CTref = refData.CTmean(:,refIdx);
    refIdx = 0; % to make sure this idx isn't excluded from plot later
else
    disp('using normalization from same file');
    CTref = CTmean(:,refIdx);
end

barefname = [targets{:}];

%% normalize

% index into data.samples for ctrl sample (e.g. untreated)
ctrlSampleIdx = 1; 

CTnorm = zeros([Nsamples Ntargets]);
% normalize to normalizing / housekeeping gene
for targeti = 1:Ntargets
    CTnorm(:,targeti) = -(CTmean(:,targeti) - CTref);
end
% normalize to ctrl sample
for targeti = 1:Ntargets
    CTnorm(:,targeti) = CTnorm(:,targeti) - CTnorm(ctrlSampleIdx,targeti);
end

[~,p] = sort(str2double(data.samples));
CTnorm = CTnorm(p,:);

%% make individual bar graphs for each target

labels = categorical(conditionLabels);

targetIndices = 2:data.Ntargets; % skipping 1 here because that is ATP5O

for targetIdx = targetIndices 

    labels = categorical(data.samples);
	labels = reordercats(labels,data.samples);
    
    bar(labels, CTnorm(:, targetIdx))
    ylabel('normalized CT');
    title(data.targets(targetIdx));
    fs = 16;
    cleanSubplot(fs)
%     set(gca, 'LineWidth', 2);
%     set(gca,'FontSize', fs)
    
    saveas(gcf, ['normalizedCT_' data.targets{targetIdx} '.png']);
end
%% make individual bar graphs for each target with error bar

labels = categorical(conditionLabels);

targetIndices = 2:data.Ntargets; % skipping 1 here because that is ATP5O

for targetIdx = targetIndices 

    labels = categorical(data.samples);
	labels = reordercats(labels,data.samples);
    
    bar(labels, CTnorm(:, targetIdx))
    ylabel('normalized CT');
    title(data.targets(targetIdx));
    fs = 16;
    cleanSubplot(fs)
%     set(gca, 'LineWidth', 2);
%     set(gca,'FontSize', fs)
    
    means = CTnorm(:, targetIdx);
    std_devs = CTstd(:, targetIdx);
    hold on;
    %errorbar(labels, means, std_devs, '.m', 'LineWidth', 2); 
    errorbar(labels, means, std_devs, '.', 'Color', [77/255, 190/255, 238/255], 'LineWidth', 4);
    saveas(gcf, ['normalizedCT_' data.targets{targetIdx} '_errorbar.png']);
    hold off;
end
    
%% make a combined bar graph

figure
%instead of making a bar plot with a categorical x axis, use a numerical
%axis and then add the labels in the property axes.XTickLabel
% b = bar(labels, CTnorm(:,targetIndices));
b = bar(1:length(labels), CTnorm(:,targetIndices));
ylabel('normalized CT');
title('FGF2 4 17 control ESI comparison');
fs = 16;
cleanSubplot(fs)
hold on;
% errorbar(labels,  CTnorm(:,targetIndices),  CTstd(:,targetIndices), 'k.', 'LineWidth', 1.5); % Adding black error bars
for bi = 1:length(b)
    xp = b(bi).XEndPoints; yp = b(bi).YEndPoints; errp = CTstd(:,targetIndices(bi));
    errorbar(xp(:),yp(:),errp(:),'k.','LineWidth',1.5);
end
hold off;
legend(b,data.targets(targetIndices),'Location','southwest')
set(gca,'XTick',1:length(labels),'XTickLabel',labels)
xtickangle(30)

saveas(gcf, fullfile(dataDir,'normalizedCT_combined_errorbar.png'));
saveas(gcf, fullfile(dataDir,['normalizedCT_' data.targets{targetIndices} '.png']));

