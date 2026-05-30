function tracks = graphtracks(G, ntime)
%GRAPHTRACKS Convert a tracking digraph into per-track index matrices.
%   TRACKS = GRAPHTRACKS(G, NTIME) splits the tracking digraph G into weakly-
%   connected components (lineages) and, for each, returns an NTIME-by-(number
%   of leaves) matrix whose columns are root-to-leaf paths (node indices placed
%   by frame). Turns graph-based tracking results into time-indexed tracks.

Din = indegree(G);
Dout = outdegree(G);
Frames = G.Nodes.frame;
[bins, binsizes] = conncomp(G,'Type','weak');
ntracks = length(binsizes);

tracks = cell(ntracks,1);

for jj = 1:ntracks
    if mod(jj,round(ntracks/20)) == 0
        fprintf('.')
    end
    list = find(bins == jj);
    start = list(Din(list) == 0);
    ends = list(Dout(list) == 0);
    numcols = length(ends);

    I = NaN(ntime,numcols);

    for ii = 1:numcols
        idxs = shortestpath(G,start,ends(ii));
        I(Frames(idxs),ii) = idxs;
    end
    
    tracks{jj} = I;
end
fprintf('\n')

end