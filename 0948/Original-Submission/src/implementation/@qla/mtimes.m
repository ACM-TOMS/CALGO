function s = mtimes(s, t)
if isa(s, 'numeric') && isa(t, 'qla') % t is qla object
    s = t;
elseif isa(s, 'qla') && isa(t, 'qla') % both are qla objects
    s.offset = min(s.offset, t.offset);
    tmp = max([s.type; t.type]);
    indx = s.type ==2 & t.type==2;
    tmp(indx) = 3;
    s.type = tmp;
end
end