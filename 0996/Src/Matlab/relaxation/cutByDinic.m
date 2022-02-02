%CUTBYDINIC Compute a min s-t cut by Dinic's max flow algorithm
%Usage:
%   [val, S] = cutByDinic(Capacity, s, t);
%
%Input: 
%  Capacity: Adjacency matrix (in sparse format) of a directed graph
%  s: the index of the source node
%  t: the index of the sink node
%
%Output:
%  val: min s-t cut cost
%  S: indicator vector of the min s-t cut.