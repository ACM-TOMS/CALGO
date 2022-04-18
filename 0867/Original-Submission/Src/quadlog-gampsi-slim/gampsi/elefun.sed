s/$/ /
s/ +INF / Infinity /g
s/ +Inf / Infinity /g
s/ +inf / Infinity /g
s/ +INFINITY / Infinity /g
s/ +Infinity / Infinity /g
s/ +nan / NaN /g
s/ -nan / NaN /g
s/ INF / Infinity /g
s/ Inf / Infinity /g
s/ inf / Infinity /g
s/ INFINITY / Infinity /g
s/ NAN / NaN /g
s/ nan / NaN /g
s/ NaNQ / NaN /g
s/ NANQ / NaN /g
s/  */ /g
s/ [?][.]00*[Ee][+]0* / NaN /g
s/ [+][.][+]00*[Ee][+]01 / Infinity /g
s/  *$//

