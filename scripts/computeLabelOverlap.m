function overlap = computeLabelOverlap(mask1, mask2)
% mask1 and mask2 should be labeled images, not cell array of linear indices
% make sure that for both masks, labels start from 1
% this is a 2D matrix that stores the number of pixels that overlap between the two masks
overlap = zeros([max(mask1(:)), max(mask2(:))]);
sz = size(mask1);
% simply loop through all pixels and count the overlaps
for i = 1:sz(1)
    for j = 1:sz(2)
        if mask1(i,j) > 0 && mask2(i,j) > 0
            overlap(mask1(i,j), mask2(i,j)) = overlap(mask1(i,j), mask2(i,j)) + 1;
        end
    end
end
