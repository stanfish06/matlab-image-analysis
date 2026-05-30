function [barenames, suffbare] = parse_MIP_names(mipDir,nucChannel,ti)
%PARSE_MIP_NAMES Find MIP files and parse their base names.
%   [BARENAMES, SUFFBARE] = PARSE_MIP_NAMES(MIPDIR, NUCCHANNEL, TI) scans
%   MIPDIR for maximum-intensity-projection files of the nuclear channel
%   NUCCHANNEL (0-based) at time TI (default 1), auto-detecting the image
%   extension (.jpg/.tif/.png) and whether names carry a time index.
%   BARENAMES is the per-file prefix (text before _w####); SUFFBARE is the
%   sprintf suffix pattern for the channel/time fields.

if ~exist('ti','var')
    ti = 1;
end

%determine mip extension/filetype
exts = {'.jpg','.tif','.png'};
flag = true;
ii = 1;
while flag
    listing = dir(fullfile(mipDir,sprintf(['*_MIP*_w%.4d*',exts{ii}],nucChannel)));
    if ~isempty(listing)
        ext = exts{ii};
        flag = false;
    elseif ii > length(exts)
        error('no MIP files matching the expected pattern were found')
    end
    ii = ii + 1;
end

%MIP files may or may not end with time index (e.g., end with _t0000.tif),
%depending on preprocessing options
% suffbare = ['_MIP*_w%.4d',ext];
suffbare = ['_w%.4d',ext];
pattern = sprintf(['*_MIP*',suffbare],nucChannel);
listing = dir(fullfile(mipDir,pattern));

if isempty(listing)
%     suffbare = ['_MIP*_w%.4d_t%.4d',ext];
    suffbare = ['_w%.4d',sprintf('_t%.4d',ti-1),ext];
    pattern = sprintf(['*',suffbare],nucChannel);
    listing = dir(fullfile(mipDir,pattern));
end


names = {listing.name};
barenames = cell(size(names));
for ii = 1:length(barenames)
%     I = strfind(names{ii},'_MIP_');
    I = strfind(names{ii},sprintf('_w%.4d',nucChannel));
    barenames{ii} = names{ii}(1:I-1);
end


end