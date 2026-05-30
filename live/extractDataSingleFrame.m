function [dbinfo, cellData] = extractDataSingleFrame_stan(img,seg,opts)
% for now, assuming that there is only one channel
dataChannels = opts.dataChannels;
nucmask = seg;
warning('Using inverse dilated nuclear mask as background')
bgmask = imdilate(nucmask,strel('disk',5)) == 0;

if opts.nucShrinkage > 0
    nucmask = bwmorph(nucmask,'shrink',opts.nucShrinkage);
end
nucCC = bwconncomp(nucmask);
nucstats = regionprops(nucCC, 'Area', 'Centroid','PixelIdxList',...
    'Orientation', 'MajorAxisLength','MinorAxisLength',...
    'Circularity');

nCells = numel(nucstats);
cellData.XY = cat(1,nucstats.Centroid);

% nuclear geometry
cellData.nucArea = cat(1,nucstats.Area);
cellData.nucOrientation = cat(1,nucstats.Orientation);
cellData.nucMajorAxis = cat(1,nucstats.MajorAxisLength);
cellData.nucMinorAxis = cat(1,nucstats.MinorAxisLength);
cellData.nucCircularity = cat(1,nucstats.Circularity);

cellData.PixelIdxList = {nucstats.PixelIdxList};
cellData.nucLevel = zeros([nCells numel(opts.dataChannels)]);
cellData.background = zeros([1 numel(opts.dataChannels)]);

if nCells > 0
    for cii = 1:numel(opts.dataChannels)
        imc = img{dataChannels(cii)+1};
        if size(nucmask) ~= size(imc,[1,2])
            error(['nucmask size ' num2str(size(nucmask)) ' does not match image size ' num2str(size(imc))]);
        end
        % current low-tech background subtraction:
        %-------------------------------------------------
        % mean value of segmented empty space in the image
        % or otherwise just min of image
        if ~isempty(bgmask)
            %dont use black buffer regions around the image for
            %background estimation
            bgmask = bgmask & (imc > 0);
            if sum(bgmask,'all') > numel(bgmask)/100
                cellData.background(cii) = mean(imc(bgmask));
            else
                %figure out a better way to handle this
                cellData.background(cii) = [];
            end
        else
            cellData.background(cii) = min(imc(imc>0),[],'all');
        end
        
        for cellidx = 1:nCells
            
            nucPixIdx = nucstats(cellidx).PixelIdxList;
            cellData.nucLevel(cellidx, cii) = mean(imc(nucPixIdx));
            
        end
    end
else
    warning(['------------ NO CELLS' '------------']);
end
dbinfo = 0;
end