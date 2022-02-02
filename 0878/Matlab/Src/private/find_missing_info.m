%FIND_MISSING_INFO  Information on missing values
%
%  [KO, RO] = FIND_MISSING_INFO(MISS) returns the number of observed values
%  before time t in KM(t) and the number of observed values at time t in RO(t). 
%
%  [KO, RO, KM, RM] = FIND_MISSING_INFO(MISS) also finds in KM(t) and RM(t) the
%  number of missing values before resp. at time t.
%
%  MISS should be an r by n logical array, true in locations of missing values.

function [ko, ro, km, rm] = find_missing_info(miss)
  ro = sum(~miss, 1);
  ko = [0 cumsum(ro)];
  if nargout > 2
    rm = sum(miss, 1);
    km = [0 cumsum(rm)];
  end
end
