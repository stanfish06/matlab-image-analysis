function quantify(baseDir)
%QUANTIFY Perform single-cell quantification for one plate directory.

if nargin < 1 || isempty(baseDir)
    baseDir = pwd;
end
originalDir = pwd;
cleanupObj = onCleanup(@() cd(originalDir));
cd(baseDir);
close all;

scriptPath = pwd;
baseDir = scriptPath;

%% setup

dirs = {''};

dataDirs = fullfile(baseDir,dirs);
nrounds = length(dataDirs);
r1 = 1; %first round; should be 1 in general

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
imageType = 'disordered';
% radii of micropatterns in micron for each condition (ignored if not 'MP')
radii = 700*ones(1,metas{r1}.nWells)/2; 

segDir = dataDirs{r1}; %segment images in the first directory

logmsg('quantify: %d positions, %d round(s), imageType=%s, nucChannel=%d', npos, nrounds, imageType, nucChannel);
disp('loading nuclear and cytoplasmic masks')
tic
masks = cell(1,npos); cellData = cell(1,npos); bgmasks = cell(1,npos);
for pi = 1:npos
    if isempty(metas{r1}.fileNames)
        [~,bare,~] = fileparts(metas{r1}.filenameFormat);
        id = pi - 1;
        prefix = sprintf(bare,id,nucChannel,0);
    else
        [barefname, id] = parseFilename(metas{r1}.fileNames{pi},segDir);
        if isempty(id)
            bare = [barefname,'_w%.4d_t%.4d'];
            prefix = sprintf(bare,nucChannel,0);
        else
            bare = [barefname,'_p%.4d_w%.4d_t%.4d'];
            prefix = sprintf(bare,id,nucChannel,0);
        end
    end
    segname = [prefix,'_masks.mat'];
    mask = load(fullfile(segDir,segname));
    masks{pi} = mask.masks; 
    % cellData{pi} = mask.cellData; 
    bgmasks{pi} = mask.bgmask;
end
toc

%% iterate over colonies, folders, and channels, and read out intensities

if strcmp(imageType,'MP')
    positions(metas{r1}.nPositions) = Colony();
elseif strcmp(imageType,'disordered')
    positions(metas{r1}.nPositions) = Position();
end

tic
%parpool(4)
parfor pi = 1:metas{r1}.nPositions

    condi = find(pi >= metas{r1}.conditionStartPos,1,'last'); % condition index
          
    logmsg('quantify: position %d/%d start', pi, npos);

    if strcmp(imageType,'MP')
        positions(pi) = Colony(metas{r1}, pi);
        positions(pi).setRadius(radii(condi), metas{r1}.xres);
        positions(pi).well = condi;
    elseif strcmp(imageType,'disordered')
        positions(pi) = Position(metas{r1}, pi);
    end

    % shared properties between micropattern and disordered
    positions(pi).dataChannels = 1:positions(pi).nChannels;
    positions(pi).ncells = length(masks{pi});

    % Compute XY positions, Z, and nuclear area from mask indices
    % Need image size to convert linear indices to subscripts
    % Load one image to get dimensions
    img_for_size = positions(pi).loadImage(dataDirs{r1}, nucChannel, 1);
    img_size = size(img_for_size);
    if numel(img_size) == 2
        img_size = [img_size 1];  % Add z dimension if 2D
    end

    ncells = length(masks{pi});
    XY = zeros(ncells, 2);
    Z = zeros(ncells, 1);
    nucArea = zeros(ncells, 1);

    for ci_mask = 1:ncells
        nucmask_idx = masks{pi}(ci_mask).nucmask;
        if ~isempty(nucmask_idx)
            % Convert linear indices to subscripts [row, col, z]
            [row, col, z] = ind2sub(img_size, nucmask_idx);
            % XY is [col, row] to match image coordinates (x=col, y=row)
            XY(ci_mask, 1) = mean(col);  % X coordinate (column)
            XY(ci_mask, 2) = mean(row);  % Y coordinate (row)
            Z(ci_mask) = mean(z);        % Z coordinate
            nucArea(ci_mask) = length(nucmask_idx);  % Volume in voxels
        end
    end

    % Initialize cellData struct with computed values
    positions(pi).cellData = struct();
    positions(pi).cellData.XY = XY;
    positions(pi).cellData.Z = Z;
    positions(pi).cellData.nucArea = nucArea;

    if ~isempty(metas{r1}.fileNames)
        [barefname, id] = parseFilename(metas{r1}.fileNames{pi},segDir);
        positions(pi).filenameFormat = metas{r1}.fileNames{pi};
        positions(pi).setID(id+1);
    else
        positions(pi).filenameFormat = metas{r1}.filenameFormat;
        positions(pi).setID(pi);
    end

    ncells = positions(pi).ncells;
    channelCounts = cellfun(@(m) m.nChannels, metas);
    totalChannels = sum(channelCounts);
    nucLevels = NaN(ncells, totalChannels);
    cytLevels = NaN(ncells, totalChannels);
    NCratios = NaN(ncells, totalChannels);
    CNratios = NaN(ncells, totalChannels);
    bgs = NaN(1, totalChannels);
    channelOffset = 0;

    for ri = 1:nrounds
        fprintf('round %d of %d\n',ri,nrounds)
        nchannels = metas{ri}.nChannels;
        nucLevel = NaN(ncells, nchannels); 
        cytLevel = NaN(ncells, nchannels);
        NCratio = NaN(ncells, nchannels); 
        CNratio = NaN(ncells, nchannels); 
        BG = NaN(1,nchannels);
        for ci = 1:nchannels
            img = positions(pi).loadImage(dataDirs{ri},ci-1,1);
            
            if ri > 1 %%%TODO - add code to find and load the alignment info for this%%%
                %apply shift and warping to the image stack
                img = xyalignImageStack(img,shiftyx);
                %apply z shift + scaling
                img = zShiftScale(img,shiftz,scalez,nz1);
            end
            
            [nL,cL,ncR,bg] = readIntensityValues(img,masks{pi},bgmasks{pi});
            nucLevel(:,ci) = nL; 
            cytLevel(:,ci) = cL; 
            NCratio(:,ci) = ncR;
            CNratio(:,ci) = 1/ncR;
            BG(ci) = bg;
        end
        channelIdx = channelOffset + (1:nchannels);
        nucLevels(:,channelIdx) = nucLevel;
        cytLevels(:,channelIdx) = cytLevel;
        NCratios(:,channelIdx) = NCratio;
        CNratios(:,channelIdx) = CNratio;
        bgs(channelIdx) = BG;
        channelOffset = channelOffset + nchannels;
    end

    positions(pi).cellData.nucLevel = nucLevels;
    positions(pi).cellData.cytLevel = cytLevels;
    positions(pi).cellData.NCratio = NCratios;
    positions(pi).cellData.CNratio = CNratios;
    positions(pi).cellData.background = bgs;
    logmsg('quantify: position %d/%d done (%d cells)', pi, npos, ncells);

    if strcmp(imageType,'MP')
        % margin around nominal radius to exclude cells from, in micron
        % 96h colonies extend beyond 350 um
        margin = 50; 
        positions(pi).setCenter(margin);
    end
end
toc

meta = metas{r1};
meta.channelLabel = cat(2,channelLabel{:});

save(fullfile(baseDir,'positions.mat'),'positions')
save(fullfile(baseDir,'meta_combined.mat'),'meta')

end
