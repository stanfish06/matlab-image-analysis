function newseg = newNuclearCleanup(seg, fgChannels, divChannels, opts)
%NEWNUCLEARCLEANUP Clean a multi-class (Ilastik) nuclear segmentation.
%   NEWSEG = NEWNUCLEARCLEANUP(SEG, FGCHANNELS, DIVCHANNELS, OPTS) builds a
%   cleaned binary nuclear mask from the labeled segmentation SEG. FGCHANNELS
%   are foreground (nucleus) class labels and DIVCHANNELS the dividing-cell
%   labels; large dividing objects are merged into the nuclear mask, fused
%   nuclei are split (separate_fused), class overlaps removed, and small/large
%   objects filtered via OPTS.decompopts / OPTS.cleanupOptions. Calls
%   nuclearCleanup for the final cleanup pass.

maxArea = opts.decompopts.maxArea;
tau1 = opts.decompopts.tau1;
tau2 = opts.decompopts.tau2;
logmsg('newNuclearCleanup: cleaning (maxArea=%g, tau1=%g, tau2=%g)', maxArea, tau1, tau2);

segChannel = fgChannels(~ismember(fgChannels, divChannels));
nucseg = seg == segChannel;
divmask = ismember(seg, divChannels);

overlap = zeros(size(seg));
for ii = 1:length(fgChannels)
    overlap = overlap + bwmorph(seg == fgChannels(ii),'thicken',2);
end
overlap = overlap > 1;

%transfer large objects from divmask to nucseg
CC = bwconncomp(divmask);
sizes = cellfun(@numel, CC.PixelIdxList);
idxs = CC.PixelIdxList(sizes > maxArea);
idxs = cell2mat(idxs(:));
logmsg('newNuclearCleanup: %d dividing objects, %d large -> moved to nuclei', CC.NumObjects, sum(sizes > maxArea));
divmask(idxs) = 0;
nucseg(idxs) = 1;

newseg = separate_fused(nucseg, tau1, tau2, opts.decompopts);
newseg = (newseg | divmask) & ~overlap;
newseg = nuclearCleanup(newseg, opts.cleanupOptions);

%remove big objects
if isfield(opts.cleanupOptions,'maxArea')
    props = regionprops(newseg,'Area','PixelIdxList');
    idxs = cell2mat({props([props.Area]>opts.cleanupOptions.maxArea).PixelIdxList}');
    newseg(idxs) = 0;
end

end