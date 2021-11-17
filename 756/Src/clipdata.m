function wc = clipdata(w,lim)
%CLIPDATA (not intended for calling directly by the user)
%	When adding curves one by one to a plot, it may be desirable to
%	override the defualt behavior of redrawing the entire plot each
%	time. This is done with the line's EraseMode property.
%	Unfortunately, in Matlab 4.0, with no-redraw modes, lines are not
%	always clipped to the axes box.  This routine clips manually.  The
%	input vector W should be a vector of closely spaced complex points
%	tracing out a smooth curve.  A cluster of points outside the box is
%	replaced with NaN's, except for the first and last points of the
%	cluster.
%
%       See also HPPLOT, DPLOT, DEPLOT, STPLOT, RPLOT.
%
%       Written by Toby Driscoll.  Last updated 5/24/95.

x = real(w);
y = imag(w);
outside = (x < lim(1)) | (x > lim(2)) | (y < lim(3)) | (y > lim(4));
dout = diff(outside);
kill = outside & [1;dout~=1] & [dout~=-1;1];
wc = w;
junk = NaN;
wc(kill) = junk(ones(size(wc(kill))));


