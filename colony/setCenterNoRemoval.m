function [cm, outside] = setCenterNoRemoval(positions, pidx, radiusMicron, xres, margin)
% setCenter()
% setCenter(margin) % margin in microns

if ~exist('margin','var')
    margin = 20;
end

radiusPixel = radiusMicron/xres;

ntime = positions(pidx).nTime;

cm = zeros(ntime,2);
for ti = 1:ntime
    
    % extractData setting center based on mean cell centroid
    areas = positions(pidx).cellData(ti).nucArea;
    xys = positions(pidx).cellData(ti).XY;
    CM = mean(areas.*xys)/mean(areas);
    cm(ti,:) = CM;
    
    % exclude cells/junk outside colony
    d = sqrt((xys(:,1) - cm(ti,1)).^2 ...
            + (xys(:,2) - cm(ti,2)).^2);
        
    outside = d > radiusPixel + margin/xres;
    xys = xys(~outside,:);
    areas = areas(~outside);

    % recenter after removing junk outside
    CM = mean(areas.*xys)/mean(areas);
    cm(ti,:) = CM;
end

end
