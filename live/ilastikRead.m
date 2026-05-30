function seg = ilastikRead(fname,binarize)
%ILASTIKREAD Read an Ilastik HDF5 segmentation/probability export.
%   SEG = ILASTIKREAD(FNAME, BINARIZE) reads the '/exported_data' dataset from
%   the Ilastik output FNAME and returns it as a MATLAB image. Handles 'Simple
%   Segmentation', 'Probabilities' (thresholded at 0.5 for two classes, else
%   argmax) and 'Object Predictions' exports. BINARIZE (default true) controls
%   thresholding. The result is permuted to MATLAB (row,col) order.

if ~exist('binarize','var')
    binarize = true;
end

seg = squeeze(h5read(fname,'/exported_data'));
if contains(fname, 'Simple Segmentation')
    if binarize
        seg = seg == 1;
    end
elseif contains(fname, 'Probabilities')
    %take first ilastik class and threshold if you have probabilities
    nclass = size(seg,1);
    if nclass == 2
        if binarize
            seg = squeeze(seg(1,:,:,:) > 0.5);
        else
            seg = squeeze(seg(1,:,:,:));
        end
    else
        %option for more than two classes that does not discretize classes
        %but is compatible with existing pipeline?
        fprintf('discretizing %d classes\n',nclass)
        [~,seg] = max(seg,[],1);
        seg = squeeze(seg);
    end
elseif contains(fname, 'Object Predictions')
    %no extra step needed for object classification
else
    error('Unknown output type')
end

seg = permute(seg,[2,1,3]);

end
