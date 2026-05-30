function newims = flatFieldCorrection(ims,G,D)
%FLATFIELDCORRECTION Apply flat-field (gain/dark) correction to image tiles.
%   NEWIMS = FLATFIELDCORRECTION(IMS, G, D) corrects each image in the cell
%   array IMS as (IM - D).*G + mean(D), where G is the per-pixel gain map and
%   D is the dark/background image. Returns a cell array NEWIMS of uint16
%   corrected images.
Dmean = mean(D,'all');

newims = cell(size(ims));
for ii = 1:size(ims,1)
    for jj = 1:size(ims,2)
        newims{ii,jj} = uint16((double(ims{ii,jj}) - D).*G + Dmean);
    end
end

end
