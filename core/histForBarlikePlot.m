function [x,y] = histForBarlikePlot(bins,n)
%HISTFORBARLIKEPLOT Build x/y vertices to draw a histogram as a step outline.
%   [X,Y] = HISTFORBARLIKEPLOT(BINS, N) converts bin edges BINS and counts N
%   into coordinates that PLOT renders as a bar-like (staircase) outline,
%   instead of using BAR. Columns of N produce columns of Y.

    x = sort([bins bins(2:end-1)]);
    y = cat(1,n(1:end-1,:),n(1:end-1,:));
    
    y(1:2:end,:) = n(1:end-1,:);
    y(2:2:end,:) = n(1:end-1,:);
end