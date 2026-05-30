function knnCellDensity(positions,k)
%KNNCELLDENSITY Add a per-cell k-nearest-neighbor local density to cellData.
%   KNNCELLDENSITY(POSITIONS, K) computes, for every cell at every time point,
%   a local density estimate K/(pi*d^2) where d is the distance to the K-th
%   nearest neighbor, and stores it in POSITIONS(p).cellData(t).density.
%   POSITIONS is a Position array (handle objects, modified in place).

ntime = positions(1).nTime;
npositions = length(positions);

for pidx = 1:npositions
    for ti = 1:ntime
        %get locations of cell centroids in the position and time selected
        XY = positions(pidx).cellData(ti).XY;
        %find squared distance between all pairs of observations and sort 
        %each column in ascending order
        D = sort(squareform(pdist(XY,'squaredeuclidean')));
        %take the kth smallest distance to other points
        d = D(k+1,:);
        %calculate the approximate cell density around each cell
        density = k./(pi*d);
        positions(pidx).cellData(ti).density = density;
    end
end


end