function condlabels = listConditionLabels(meta)
%LISTCONDITIONLABELS Condition label for each position.
%   CONDLABELS = LISTCONDITIONLABELS(META) returns a 1-by-nPositions cell
%   array giving the condition name (from META.conditions) for every
%   position, using META.conditionStartPos and META.posPerCondition.

npos = meta.nPositions;
conditionStartPos = meta.conditionStartPos;

if npos ~= sum(meta.posPerCondition)
    error('posPerCondition does not match nPositions')
end

condlabels = cell(1,npos);
for ii = 1:npos
    condi = find(ii >= conditionStartPos,1,'last'); % condition index
    condlabels{ii} = meta.conditions{condi};
end


end