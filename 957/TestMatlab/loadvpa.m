function a=loadvpa(filename,nheader_lines,dig,n,m)
%LOADVPA Loading a multi-precision source file.
%LOADVPA(FILENAME,NHEADER_LINES,DIG,N,M) loads a source 
%   text file named FILENAME, which contains NHEADER_LINES 
%   comment lines followed by an NxM array of DIG-digit
%   numbers, into the Matlab working window, ignoring the 
%   comment lines. The file name has to be placed between 
%   single quotes when the function is called.      

fid=fopen(filename,'r');
%Skip first nheader_lines of header
for i =  1:nheader_lines
  fgetl(fid);
end
C = textscan(fid,'%s');
C = C{:};
fclose(fid);
for k=1:numel(C)
  a(k)=vpa(C{k},dig);
end
a=reshape(a,m,n)';
