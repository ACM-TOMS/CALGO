%LOADVPA  Load multiple-precision array
%
function a=loadvpa(filename,dig,n,m)
fid=fopen(filename,'r');
C=textscan(fid,'%s');
C=C{1,1};
fclose(fid);
%
% Define a function that returns 1 when the input
% is not a number represented as a string
%
f=@(x) isnan(str2double(x));
%
% str2double returns NaN if the input string cannot
% cannot be converted to a number
% isnan returns 1 if the input is NaN
%
% Apply the function to every cell of C
%
ind=cellfun(f,C);
%
% Remove cells for which "ind" is 1
%
C(ind)=[];
for k=1:numel(C)
  a(k)=vpa(C{k},dig);
end
a=reshape(a,m,n)';
