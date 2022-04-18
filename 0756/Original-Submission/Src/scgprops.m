function [trace,tol,v1,v2] = scgprops(fig)
%SCGPROPS (not intended for calling directly by the user)
%	Read current values from SCGUI Properties window.
%	
%	Written by Toby Driscoll.  Last updated 5/23/95.

menus = get(fig,'user');
props = get(menus(4,1),'user');
trace = props(3);
tol = 10^(-props(4));
v1 = props(5);
v2 = props(6);
