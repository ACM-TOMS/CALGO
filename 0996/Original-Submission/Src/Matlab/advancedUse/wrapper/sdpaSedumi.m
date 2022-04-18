function [ x,y,info ] = sdpaSedumi( A,b,c,K )
%SDPASEDUMI   Solve the cone problem in sedumi format by SDPA
%
%    Usage:
%       [x,y,info] = sdpaSedumi(A,b,c,K);

    [x,y,info] = sedumiwrap(A,b,c,K);
end

