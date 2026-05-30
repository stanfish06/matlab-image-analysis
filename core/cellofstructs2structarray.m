function structarray = cellofstructs2structarray(cellofstructs)
%CELLOFSTRUCTS2STRUCTARRAY Convert a cell array of structs to a struct array.
%   STRUCTARRAY = CELLOFSTRUCTS2STRUCTARRAY(CELLOFSTRUCTS) copies every field
%   of each struct in the cell array CELLOFSTRUCTS into element i of the
%   returned struct array. All structs are assumed to share the fields of the
%   first one.

    tmp = cellofstructs;
    fields = fieldnames(tmp{1});

    clear allmasks
    structarray(length(tmp)) = struct();

    for zi = 1:length(tmp) 
        for fi = 1:length(fields)
            structarray(zi).(fields{fi}) = tmp{zi}.(fields{fi});
        end
    end
end
    