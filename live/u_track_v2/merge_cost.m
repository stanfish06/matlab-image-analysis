function m = merge_cost(xy_dist, t_dist, iendA, jk_A, jprev_A)
%MERGE_COST Linking cost for a candidate cell-merge event.
%   M = MERGE_COST(XY_DIST, T_DIST, IENDA, JK_A, JPREV_A) returns the cost of
%   merging two tracks given spatial distance XY_DIST, temporal gap T_DIST and
%   areas, penalizing area mismatch via the ratio JK_A/(IENDA + JPREV_A). Used
%   by the 'updated tracking' linkers.

pij = jk_A/(iendA + jprev_A);
if pij > 1
    m = (xy_dist + t_dist^2)*pij;
elseif pij <= 1 && pij > 0
    m = (xy_dist + t_dist^2)/pij^2;
elseif pij <= 0
    disp('Something is wrong')
end

end