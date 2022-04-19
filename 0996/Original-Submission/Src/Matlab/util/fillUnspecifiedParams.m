function params = fillUnspecifiedParams(params, defparams)
%FILLUNSPECIFIEDPARAMS   Fill unspecified parameters with default values
%   Usage:
%      params = FILLUNSPECIFIEDPARAMS(params, defparams);
    deffields = fieldnames(defparams);
    
    for ii = 1:length(deffields)
        fieldname = deffields{ii};
        if ~isfield(params, fieldname)
            params.(fieldname) = defparams.(fieldname);
        end
    end
end