function [upperleft,links] = registerImageGrid_nominal(imgs, pixelOverlap)
    % register a grid of overlapping images
    % 
    % [upperleft,links] = registerImageGrid(imgs, pixelOverlap)
    %
    %
    % upperleft:    cell array of positions of upperleft corners with
    %               upperleft corner of upperleft image being (1,1)
    % links:        a matrix containing correlation links between images
    %               empty positions will correspond to cells themselves
    %
    % imgs:         cell array of images 
    %               with rows and cols corresponding to grid
    % pixelOverlap: width of overlapping strip in pixels
    %               if left empty, upperleft for a square grid of images
    %               will be returned
    
    
    
    %generate upperleft variable without actual registration - assume the
    %images are in their nominal location based on pixelOverlap values

% find the first actual image on the grid
emptycells = cellfun(@isempty,imgs,'UniformOutput',false);
emptycells = cat(1,emptycells{:});
imsize = size(imgs{find(~emptycells,1,'first')},[1,2]);
M = imsize(1); N = imsize(2);

if numel(pixelOverlap) == 1
    pixelOverlapY = pixelOverlap;
    pixelOverlapX = pixelOverlap;
elseif numel(pixelOverlap) == 2
    pixelOverlapY = pixelOverlap(1); pixelOverlapX = pixelOverlap(2);
end



% make upperleft for square grid of images (no stitching)
grid = size(imgs,[1,2]);
m = grid(1); n = grid(2);

upperleft = cell(grid);
for ii = 1:m
    for jj = 1:n
        upperleft{ii,jj} = [(ii - 1)*(M - pixelOverlapY), (jj - 1)*(N - pixelOverlapX)];
    end
end

end
