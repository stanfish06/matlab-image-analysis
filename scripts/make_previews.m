% MAKE_PREVIEWS  Build per-position preview images from MIPs [template].
%   Loads meta.mat, locates MIP files (parse_MIP_names), labels positions by
%   condition, and generates preview images for quick inspection. Run from the
%   data directory after MIPs have been created.
clear; close all; clc

%% setup
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
mipDir = fullfile(dataDir,'MIP');

load(fullfile(dataDir,'meta.mat'),'meta')
conditions = meta.conditions;
channelLabel = meta.channelLabel;
nucChannel = meta.nucChannel;
ti = 1;

[barenames, suffbare] = parse_MIP_names(mipDir,nucChannel,ti);
condlabels = listConditionLabels(meta);

%if you want to label positions based on filenames:
% [~,poslabels] = findRepeatStrings(barenames);

%if you want to label positions based on condition
poslabels = condlabels;

npos = length(barenames);
ncond = length(conditions);
ppc = npos/ncond;
nc = length(channelLabel);

%choose positions for which to make combined previews
condids = 1:meta.nWells; %choose conditions for which to make combined previews
posids = meta.conditionStartPos(condids); %this will use the first position from each condition
% posids = []; %alternatively make a generic list of positions

%set tolerance for each channel
tols = repmat({[0.01,0.99]},1,nc); %use the same tolerance for each channel
%tols = {[0.01 0.99],[0.2 0.99],[0.1 0.9],[0.05 0.995]}; %or manually set tolerance for each channel separately

%for the color overlay of different channels, choose the combination in rgb order
%colorcombo = [2 3 4] will make channel 2 red, channel 3 green, channel 4
%blue; colorcombo = [4 3 2] will make 4 red and 2 blue instead
colorcombo = [2 3 4];

savedir = fullfile(dataDir,'previews','common_contrast');
if ~exist(savedir,'dir'), mkdir(savedir); end

%% load images, downsample, adjust contrast
imgs = cell(npos,nc);
lims = NaN(nc,2);
%optionally specify the position to use for contrast adjustment for each channel
%blank entries use the brightest position for that channel
% e.g., pids = {[],3,5,[]};
pidxs = cell(1,nc); %by default use the brightest image in each channel

for ci = 1:nc
    disp(meta.channelLabel{ci})
    for pidx = 1:npos
        fprintf('.')
        bare = [barenames{pidx},suffbare];
        fname = fullfile(mipDir,sprintf(bare,ci-1));
        im = imread(fname);
        m = size(im,1); n = size(im,2); ss = round(sqrt(m*n/512^2)); %subsampling factor
        
        A = im2double(im);
        A = imfilter(A,ones(ss)/ss^2,'symmetric');
        A = A(1:ss:end,1:ss:end);
        mnew = size(A,1); nnew = size(A,2);
        imgs{pidx,ci} = A;
    end
    fprintf('\n')
    
    means = cellfun(@(x) mean(x,'all'),imgs(:,ci));
    I = ids{ci};
    if isempty(I)
        [~,I] = max(means,[],'all','linear');% [row,col] = ind2sub(size(means),I);
    end
    lim = stitchedlim(imgs{I,ci},tols{ci});
    lims(ci,:) = lim;
    
    for pidx = 1:npos
        imgs{pidx,ci} = imadjust(imgs{pidx,ci},lims(ci,:));
    end
end

%% make and save previews
cfs = 0.01;
fontsize = 24;

C = cell(length(posids),nc+1);

for pidx = 1:npos
%     pidx = posids(pii);
    pii = find(posids == pidx);
    As = cell(1,nc+1);
    for ci = 1:nc
        A = imgs{pidx,ci};
        %channel label
        A = insertText(A,[0.015*nnew, mnew*(1-0.5*cfs)],channelLabel{ci},...
            'BoxColor','black','TextColor','white','BoxOpacity',0.3,...
            'FontSize',fontsize,'Font','Arial Bold','AnchorPoint','LeftBottom');

        %condition label
        if ci == 1
            A = insertText(A,[0.015*nnew, mnew*0.5*cfs],poslabels{pidx},...
                'BoxColor','black','TextColor','white','BoxOpacity',0.3,...
                'FontSize',fontsize,'Font','Arial Bold','AnchorPoint','LeftTop');
        end
        A = [ones(2,nnew+4,3); [ones(mnew,2,3), A, ones(mnew,2,3)]; ones(2,nnew+4,3)];
        
        As{ci} = A;
        if ~isempty(pii)
            C{pii,ci} = A;
        end
    end

    if ~isempty(colorcombo)
        ims = imgs(pidx,colorcombo);
        A = cat(3,ims{:});

        textcols = {'red','green','blue'};
        anchorpoints = {'LeftBottom','CenterBottom','RightBottom'};
        xcoords = [10,nnew/2,nnew-10];


        for ci = 1:3
            A = insertText(A,[xcoords(ci), mnew*(1-0.5*cfs)],channelLabel{colorcombo(ci)},...
                'BoxColor','black','TextColor',textcols{ci},'BoxOpacity',0.3,...
                'FontSize',fontsize,'Font','Arial Bold','AnchorPoint',anchorpoints{ci});
        end
    end
    A = [ones(2,nnew+4,3); [ones(mnew,2,3), A, ones(mnew,2,3)]; ones(2,nnew+4,3)];
    
    As{nc+1} = A;
    if ~isempty(pii)
        C{pii,nc+1} = A;
    end    

    A = cell2mat(As);
    figure('Position',figurePosition(2500,500))
    imshow(A,'InitialMagnification','fit','Border','tight')
    cleanSubplot
    pause(0.5)
    savename = fullfile(savedir,[barenames{pidx},'_previews.jpg']);
    imwrite(A,savename,'Quality',99)
    close all
    
end



%save an overview of all the colonies together
f = figure('WindowState','maximized');
C = cell2mat(C');
imshow(C,'InitialMagnification','fit','Border','tight')
cleanSubplot
pause(0.5)

savename = fullfile(savedir,'combined_positions.jpg');
imwrite(C,savename,'Quality',99)

close(f)
















