clear; close all; clc

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
baseDir = scriptPath;

%% setup

dirs = {''};
% dirs = {...
%     '230110_6i9rd_exp20_RD1_SMAD23_pAKT_SOX17',...
%     '230111_6i9rd_exp20_RD2_GATA3_OTX2_LEF1',...
%     };
dataDirs = fullfile(baseDir,dirs);
nrounds = length(dataDirs);
r1 = 1; %first round; should be 1 in general

bare = 'stitched_p%.4d_w%.4d_t%.4d';
channelLabel = cell(1,nrounds);
nucChannels = NaN(1,nrounds);

%load metadata from each round
metas = cell(nrounds,1);
for ri = 1:nrounds
    meta = load(fullfile(dataDirs{ri},'meta.mat'),'meta');
    metas{ri} = meta.meta;
    channelLabel{ri} = metas{ri}.channelLabel;
    nucChannels(ri) = metas{ri}.nucChannel;
end
npos = metas{r1}.nPositions;
nucChannel = nucChannels(r1);

% define imageType as micropatterned 'MP' or 'disordered'
imageType = 'MP';
% radii of micropatterns in micron for each condition (ignored if not 'MP')
radii = 700*ones(1,metas{r1}.nWells)/2; 

%-----------------------------------------------------------------
% workaround to fix positions per condition to format for multiple
% different counts
for ri = 1:nrounds
    if numel(metas{r1}.posPerCondition)==1
        posPerCondition = metas{r1}.posPerCondition;
        metas{r1}.posPerCondition = repmat(posPerCondition, [1 metas{r1}.nWells]);
    end
end
%-----------------------------------------------------------------

%%%% TODO - add code to specify a separate filelist for each condition,
%%%% then update meta.posPerCondition, meta.nPositions, and and
%%%% meta.fileNames accordingly %%%%

segDir = dataDirs{r1};

disp('loading nuclear and cytoplasmic masks')
tic
masks = cell(1,npos); cellData = cell(1,npos); bgmasks = cell(1,npos);
for pi = 1:npos
    [barefname, id] = parseFilename(metas{r1}.fileNames{pi},segDir);
    if isempty(id)
        bare = [barefname,'_w%.4d_t%.4d'];
        prefix = sprintf(bare,nucChannel,0);
    else
        bare = [barefname,'_p%.4d_w%.4d_t%.4d'];
        prefix = sprintf(bare,id,nucChannel,0);
    end
    segname = [prefix,'_masks.mat'];
    mask = load(fullfile(segDir,segname));
    masks{pi} = mask.masks; 
    cellData{pi} = mask.cellData; 
    bgmasks{pi} = mask.bgmask;
end
toc

%% iterate over colonies, folders, and channels, and read out intensities

if strcmp(imageType,'MP')
    positions(metas{r1}.nPositions) = Colony();
elseif strcmp(imageType,'disordered')
    positions(metas{r1}.nPositions) = Position();
end

conditionStartPos = metas{r1}.conditionStartPos;

tic
%parpool(4)
parfor pi = 1:metas{r1}.nPositions

    condi = find(pi >= conditionStartPos,1,'last'); % condition index
          
    fprintf('colony %d of %d\n',pi,npos)

    if strcmp(imageType,'MP')
        positions(pi) = Colony(metas{r1}, pi);
        positions(pi).setRadius(radii(condi), metas{r1}.xres);
        positions(pi).well = condi;
    elseif strcmp(imageType,'Disordered')
        positions(pi) = Position(metas{r1}, pi);
    end

    % shared properties between micropattern and disordered
    positions(pi).dataChannels = 1:positions(pi).nChannels;
    positions(pi).cellData = cellData{pi};
    positions(pi).ncells = length(masks{pi});
    
    if ~isempty(metas{r1}.fileNames)
        [barefname, id] = parseFilename(metas{r1}.fileNames{pi},segDir);
        positions(pi).filenameFormat = metas{r1}.fileNames{pi};
        positions(pi).setID(id);
    else
        positions(pi).filenameFormat = metas{r1}.filenameFormat;
        positions(pi).setID(pi);
    end

    for ri = 1:nrounds
        fprintf('round %d of %d\n',ri,nrounds)
        nchannels = metas{ri}.nChannels;
        ncells = positions(pi).ncells;
        nucLevel = NaN(ncells, nchannels); 
        cytLevel = NaN(ncells, nchannels);
        NCratio = NaN(ncells, nchannels); 
        BG = NaN(1,nchannels);
        for ci = 1:nchannels
            img = positions(pi).loadImage(dataDirs{ri},ci-1,1);
            
            if ri > 1
                %apply shift and warping to the image stack
                img = xyalignImageStack(img,shiftyx);
                %apply z shift + scaling
                img = zShiftScale(img,shiftz,scalez,nz1);
            end
            
            img = correctChromaticAberration(img,ci-1,metas{r1});
            [nL,cL,ncR,bg] = readIntensityValues(img,masks{pi},bgmasks{pi});
            nucLevel(:,ci) = nL; 
            cytLevel(:,ci) = cL; 
            NCratio(:,ci) = ncR;
            BG(ci) = bg;
        end
        if ri==1
            nucLevels = nucLevel;
            cytLevels = cytLevel;
            NCratios = NCratio;
            bgs = BG;
        else
            nucLevels = cat(2, nucLevels, nucLevel); 
            cytLevels = cat(2, cytLevels, cytLevel);
            NCratios = cat(2, NCratios, NCratio);
            bgs = cat(2, bgs, BG);
        end
    end

    positions(pi).cellData.nucLevel = nucLevels;
    positions(pi).cellData.cytLevel = cytLevels;
    positions(pi).cellData.NCratio = NCratios;
    positions(pi).cellData.background = bgs; 

    if strcmp(imageType,'MP')
        % margin around nominal radius to exclude cells from, in micron
        % 96h colonies extend beyond 350 um
        margin = 50; 
        positions(pi).setCenter(margin);
        % positions(pi).makeRadialAvgSeg(); - generates error
        % should we change plotRadialProfiles to replace this and store
        % the result in Colony.radialProfile or remove radialProfile as
        % a propery and recalculate every time?
    end
end
toc

meta = metas{r1};
meta.channelLabel = cat(2,channelLabel{:});

save(fullfile(baseDir,'positions.mat'),'positions')
save(fullfile(baseDir,'meta_combined.mat'),'meta')

%make and save cellStats
stats = cellStats(positions, meta, positions(1).dataChannels);
save(fullfile(baseDir,'stats.mat'),'stats')

%% load data if processing has already been done

load(fullfile(baseDir,'positions.mat'),'positions')
load(fullfile(baseDir,'meta_combined.mat'),'meta')

%% collect colony data and export data to a csv file

%radius of the colony in microns
radiusMicron = 350;

%colony indices for which to save data
pidxs = 1:npos;
np = length(pidxs);

%all channel names
channelLabels = meta.channelLabel;
channelLabels = renameDuplicateChannels(channelLabels);
%choose channels for which to save data
chans = 1:length(channelLabels);
nchan = length(chans);
cls = channelLabels(chans);

%find colony centers
cms = NaN(npos,2);
for pi = 1:npos
    cm = setCenter(positions, pi, radiusMicron, meta.xres);
    cms(pi,:) = cm;
end

%filename to which to save the csv
I = regexp(dirs{1},'_RD[1234567890]+_');
savename = [dirs{1}(1:I),'datamatrix.csv'];

excluded = {'XY','nucLevel','cytLevel','NCratio','background'};

fields = fieldnames(positions(1).cellData);
for fi = 1:length(excluded)
    fields = fields(~strcmp(fields,excluded{fi}));
end
nf = length(fields);

varnames = [fields',{'Colony','X','Y','CenterX','CenterY','RadialDist'}];

%collect data on nucleus geometry, cell position, etc.
geometricdata = cell(np,1);
for pi = 1:np
    pidx = pidxs(pi);
    ncell = size(positions(pidx).cellData.XY,1);
    A = NaN(ncell,nf+6);
    for fi = 1:nf
        A(:,fi) = positions(pidx).cellData.(fields{fi});
    end
    xys = positions(pidx).cellData.XY;
    d = meta.xres*pdist2(xys,cms(pidx,:));
    A(:,nf+1) = pidx;
    A(:,nf+2:nf+3) = xys;
    A(:,nf+4) = cms(pidx,1); A(:,nf+5) = cms(pidx,2);
    A(:,nf+6) = d;
    
    geometricdata{pi} = A;
end
geometricdata = cell2mat(geometricdata);

%nuclear, cytoplasmic, and N:C ratio readouts in each channel
fields = {'nucLevel','cytLevel','NCratio'};
prefs = {'','cyto','nc'};
intnames = [];
for fi = 1:length(fields)
    intnames = [intnames, strcat(prefs{fi},cls)]; %#ok<AGROW>
end

ncc = nchan*length(prefs);
intdata = cell(np,1);
for pi = 1:np
    pidx = pidxs(pi);
    ncell = size(positions(pidx).cellData.XY,1);
    A = NaN(ncell,ncc);
    for fi = 1:length(fields)
        A(:,(fi-1)*nchan + 1:fi*nchan) = positions(pidx).cellData.(fields{fi})(:,chans);
    end
    intdata{pi} = A;
end
intdata = cell2mat(intdata);

%make the collected data matrices into a table
A = array2table([intdata,geometricdata],'VariableNames', [intnames,varnames]);

%write the table to a csv file
writetable(A,fullfile(baseDir,savename))



%% local functions

function cm = setCenter(positions, pidx, radiusMicron, xres, margin)
% setCenter()
% setCenter(margin) % margin in microns

if ~exist('margin','var')
    margin = 20;
end

radiusPixel = radiusMicron/xres;

ntime = positions(pidx).nTime;

cm = zeros(ntime,2);
for ti = 1:ntime
    
    % extractData setting center based on mean cell centroid
    areas = positions(pidx).cellData(ti).nucArea;
    xys = positions(pidx).cellData(ti).XY;
    CM = mean(areas.*xys)/mean(areas);
    cm(ti,:) = CM;
    
    % exclude cells/junk outside colony
    d = sqrt((xys(:,1) - cm(ti,1)).^2 ...
            + (xys(:,2) - cm(ti,2)).^2);
        
    outside = d > radiusPixel + margin/xres;
    xys = xys(~outside,:);
    areas = areas(~outside);

    % recenter after removing junk outside
    CM = mean(areas.*xys)/mean(areas);
    cm(ti,:) = CM;
end

end




