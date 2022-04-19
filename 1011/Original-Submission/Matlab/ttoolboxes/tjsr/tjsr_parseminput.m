function [ M, type ] = tjsr_parseminput( type, args )
% [ M, type ] = tjsr_parseminputtjsr( type, args )
% This function belongs to tjsr!
%   Creates the type struct
%   Parses most of the options for tjsr()
%   Some options are parsed by tjsr() itself.
%
%   Input:
%       args        command line arguments
%   Output:
%       type        the type-struct for tjsr
%
%   Notes:
%   See the source code of the file for available options for tjsr()
%
%   See also: tjsr
%
% Written by: tommsch, 2018

%#ok<*ALIGN>

%Pre-processing
    parg = args(2:end); %parg=parse_varargin
    M = args{1};
    %type=struct;
    type.arguments = parg;                %the passed arguments
    
    if( size(M,2)==1 ); 
        M = M'; end; 
    type.M_original = M;  %matrices must be in a row-cell
    
    type.info.errorcode = NaN; %this is a not-valid error code
    type.info.errorinformation = [];      %used to store information in cases of errors
    type.counter.starttime = clock;       %starting time
    type.counter.nummatrix = numel(M);    %number of matrices
    %type.info.infotext = '';              %initialize infotext
    %type.info.errortext = '';             %initialize errortext
    type.info.dim = size(M{1}, 1 );         %the dimension
    type.counter.numstepsmall = 0;        %total number of steps conducted
    type.counter.numstepbig = 0;          %total number of steps conducted including the steps done in linprog
    type.counter.iteration = 0;           %how many times the worker looped. if we make no natural selection, then this number coincides with the true! treedepth
    

    
%preworker options
    [type.JSR,parg] =                    parsem( 'JSR', parg, [], 'expecte',@(x) numel(x)==2 );                  %an inital true! estimate for the JSR. Is used (amongst other things) to search for smp-candidates.
    [type.cyclictree.ordering,parg] =    parsem({'ordering','oo','smp'}, parg, {}, 'expect',@iscell );                %the orderings of the smp-candidates, nearlycandidates. Eg: < 'ordering',{[1 2],[1]} > 
    [type.cyclictree.smpflag,parg] =     parsem( 'smpflag', parg, [], 'expecte',@isvector );                     %(row-vector) Defines wether the candidates are an smp or not. 0=candidate, 1=nearlycandidate, 2=extravertex. If smpflag is given, ordering must be given.
    [type.cyclictree.v0,parg] =          parsem( 'v0', parg, {}, 'expect',@iscell );                            %the corresponding eigenvectors to the candidates. If v0 is given, ordering must be given , v0s should be given.
    [type.cyclictree.v0s,parg] =         parsem( 'v0s', parg, {}, 'expect',@iscell );                           %the corresponding left-eigenvectors to the candidates. If v0s is given, ordering and v0 must be given.
    [type.cyclictree.balancingvector,parg] = parsem( 'balancingvector', parg, [], 'expecte',@isvector );         %balancing-vector. Should have as many entries as ordering has (if there are only unique leading eigenvectors. Otherwise the balancing-vector must have more entries).        eg: < 'balancingvector',[ 1 0.8] >
    [type.cyclictree.extravertex,parg] =  parsem( 'extravertex', parg, {}, 'expect',@iscell );                  %extravertices can be given in 2 ways. Either as a cell array of vectors as argument of 'extravertex', or in the cell array of 'v0' and the corresponding smpflag is set to 2
    [type.cyclictree.multiplicity,parg] = parsem( 'multiplicity', parg, {}, 'expect',@iscell );               %the multiplicity of the correspdoning entries in v0,v0s,...   
    [type.opt.autoextravertex,parg] =    parsem( 'autoextravertex', parg, .01, 'expect',{'clcl',[0 1]} );       %adds a vertex for all directions whose corresponding singular value is less than maximal-singular-value.
    [type.opt.balancingdepth,parg] =     parsem( 'balancingdepth', parg, min(5,floor(10/log(type.counter.nummatrix))), 'expect',{'clcl',[0 10]} );                 
    
    [type.opt.delta,parg] =              parsem( 'delta', parg, 1, 'expect',{'opcl',[0 1]} );                   %Matrices are multiplied by Delta after the construction of the cyclic tree. A value smaller than one leads to faster convergence, but the algorithm cannot return the exact value of the JSR anymore.%level-depth used in balancing of trees.    eg: < 'balancingdepth',4 >
    [type.opt.maxnumcandidate,parg] =    parsem( 'maxnumcandidate', parg, type.counter.nummatrix*10, 'expect',{'clop',[0 inf]} );     %if there are more candidates then this, maxsmpdepth will be reduced
    [type.opt.nearlycandidate,parg] =    parsem( 'nearlycandidate', parg, .9999, 'expect',{'opcl',[0 1]} );     %all candidates with spectral radius> JSR*lambda are used as nearlycandidates
    [type.opt.nobalancing,parg] =        parsem( 'nobalancing', parg, [], 'expecte',{-1 0 1} );                   %no balancing vector is computed. Equivalently: balancing vector = ( 1 1 1 ... 1 ). if 'nobalancing'==-1, a balancing is done if a balancing vector was found.
    if(isempty(type.opt.nobalancing))
        if( type.opt.delta<1); 
            type.opt.nobalancing=-1; 
        else; type.opt.nobalancing=0; end; end;
    [type.opt.noclassify,parg] =         parsem( 'noclassify', parg );                                          %candidates for cyclic trees are not classified. Instead a trivial choice is made.
    [type.opt.invariantsubspace,parg] =  parsem( 'invariantsubspace', parg, 'auto', 'expect',{'auto', 'none', 'perm', 'joint'} );           %whether search for invariantsubspaces or not. Possible arguments are: 'auto','none','perm','joint'
    [type.opt.nomultipleeigenvector,parg] = parsem( 'nomultipleeigenvector', parg );          %Only takes one leading eigenvector per candidate.
    [type.opt.complexeigenvector,parg] = parsem( 'complexeigenvector', parg, 3, 'expect',{'clcl',[0 4]} );             %val=2: Remove complex eigenvectors if there is at least one real eigenvector
    [type.opt.minJSR,parg] =             parsem( 'minJSR', parg, [], 'expecte',@isnumeric );                         %Minimal value for JSR of candidate to be found. This option should not be used.
    %XX Replace options minJSR with option bound
    [type.opt.maxsmpdepth,parg] =        parsem( 'maxsmpdepth', parg, 50, 'expect',{'opop',[0 inf]} );                   %maximal length of matrixproduct tested for smp candidate.      eg: < 'maxsmpdepth',4 >
    [type.opt.findsmp_N,parg] =          parsem( 'findsmp_N', parg, ceil(sqrt(type.info.dim*type.counter.nummatrix))*10, 'expect',{'clcl',[0 inf]} ); %Number of products kept in each step of findsmp algorithm
    
    
    %worker options
    [type.info.algorithm,parg] =         parsem( 'algorithm', parg, [], 'expecte',{0, 'P', 'cone', 1 'R', 'mink', 2, 'C', 'complex'} );                     %the to be chosen algorithm. 0='P'=cone, 1='R'=mink, 2='C'=complex
    if( ~isempty(type.info.algorithm)); 
        switch type.info.algorithm;  
            case {0,1,2}; 
            case {'P','cone'};    type.info.algorithm = 0; 
            case {'R','mink'};    type.info.algorithm = 1; 
            case {'C','complex'}; type.info.algorithm = 2;  
            otherwise; error('Wrong value for ''algorithm''.'); end; end;
    [type.opt.fastnorm,parg] =           parsem( 'fastnorm', parg, 1, 'expect',{0 1 2} );                       %0, 1 or 2. Estimates the norm prior to computing it. If val>=1 and a point is surely inside, its norm is set to TJSR_INSIDE. If val>= 2 and a point is surely outside, its norm is set to TJSR_OUTSIDE==inf. Thus for val=2 thus this gives bad estimates during computation. Default=1.
    [type.opt.naturalselection,parg] =   parsem( 'naturalselection', parg, 64*feature('numcores'), 'expect',{'opcl',[0 inf]} ); %Minimum number of vertices whose norm are computed in each level.
    [type.opt.naturalselectiontype,parg] = parsem( 'naturalselectiontype', parg, inf, 'expect',{inf -inf 1 -1 2 -2 3 -3 100 -100} );        %The type how vertices to test are selected.
                                                                                                    % +/- inf=Auto (3x by estimate, 1x by parent-norm), 
                                                                                                    % +/- 1=norm-estimate, 
                                                                                                    % +/- 2=parent-norm, 
                                                                                                    % +/- 3=rho (is not working good), 
                                                                                                    % +/- 100=-rho (for debug and test reasons). See tjsr_selectvertex.m for more info.
                                                                                              %If the value is positive, then all children of a vertex are selected, if at least on of them is selected.
                                                                                              %If the value is negative, the algorithm can't produce intermediate bounds for the jsr, but the algorithm may be faster.
    [type.opt.testoldvertex,parg] =      parsem( 'testoldvertex', parg, 1, 'expect',{0 1} );                  %Test old vertices after each step. 0: never, 1=Auto, 2=Always.
    [type.opt.simplepolytope,parg] =     parsem( 'simplepolytope', parg, 1, 'expect',{'clcl',[0 1]} );                 %Vertices with distance less than val to other vertices are not used for the polytope in the norm computation. val must be greater than epslinprog in order that this has an effect
    if( type.opt.simplepolytope==1 ); 
        type.opt.simplepolytope = 1e-8; end;               %This makes the polytope smaller, thus the norm larger and the algorithm still reports correct estimates.
    [type.opt.epspolytope,parg] =        parsem( 'epspolytope', parg, max(2*10^(-15+log10(type.info.dim)+1), 2e-9), 'expect',{'clcl',[-1 1e-3]} );                %if polytope norm is larger than 1-eps, then the vertex is considered to be outside the polytope. epspolytope should be larger than epslinprog
    [type.opt.computebound,parg] =       parsem( 'computebound', parg, 0, 'expect',{0 1} );                   %flag if intermediate bounds shall be computed in each level
    [type.opt.numcore,parg] =            parsem({'numcore','numthread'}, parg, feature('numcores'), 'expect',{'clop',[1 inf]} ); %Number of threads used in computation of polytope norm. If numcore==1, then the computation is done in the main thread
    
    %termination options
    [type.opt.maxmemory,parg] =          parsem( 'maxmemory', parg, tavailable_memory-2^20, 'expect',{'clcl',[1 inf]} ); %computation stops if 'info'-field ('type'-field) is larger than 'maxmemory'.
    [type.opt.maxstepnumber,parg] =      parsem( 'maxstepnumber', parg, inf, 'expect',{'clcl',[1 inf]} );                %computation stops if more than 'maxstepnumber' vertices are tested (and the level is completely finished).      eg: < 'maxnumberofsteps',8 >
    [type.opt.maxtreetime,parg] =        parsem( 'maxtreetime', parg, inf, 'expect',{'clcl',[1 inf]} );                  %computation stops after building tree takes more than 'maxtreetime' seconds.                  eg: < 'maxtreedepth',6 >
    [type.opt.maxiteration,parg] =       parsem( 'maxiteration', parg, inf, 'expect',{'clcl',[1 inf]} );                 %computation stops after iterating the algorithm val times.
    [type.opt.maxtime,parg] =            parsem( 'maxtime', parg, inf, 'expect',{'clcl',[1 inf]} );                      %computation stops after 'maxtime' seconds (and the level is completely finished).      eg: < 'maxtime',240 >
    [type.opt.maxvertexnumber,parg] =    parsem( 'maxvertexnumber', parg, inf, 'expect',{'clcl',[1 inf]} );              %computation stops if polytope has more than 'maxvertexnumber' (and the level is completely finished).                  eg: < 'maxvertices',1000 >
    [type.opt.maxtreedepth,parg] =       parsem( 'maxtreedepth', parg, inf, 'expect',{'clcl',[1 inf]} );                 %computation stops after tree has more than 'maxtreedepth' levels.                  eg: < 'maxtreedepth',6 >
    [type.opt.testeigenplane,parg] =     parsem( 'testeigenplane', parg, -inf, 'expect',{'clcl',[-inf inf]} );              %tests whether a product is a possible smp-candidate of not in each step. val is the epsilon used. Negative values are probably prefered. %formerly: noteststoppingcriterion XX
    if( type.opt.testeigenplane==1); 
        type.opt.testeigenplane=-1e-10; end;
    [type.opt.testspectralradius,parg] = parsem( 'testspectralradius', parg, 1, 'expect',{'clcl',[0 1]} );             %tests whether the products which occur during the computation have higher spectral radius then the candidate. Tests are only done for the very first tree. This should be sufficient, since the orderings in the trees are nearly all the same.    
    if( type.opt.testspectralradius==1); 
        type.opt.testspectralradius=-1e-10; end;
    [type.opt.validateupperbound,parg] = parsem( {'validateupperbound','ub'}, parg, -inf, 'expect',{'clcl',[-inf inf]} );   %Algorithm terminates if upper estimate of JSR is less than validateupperbound. Also changes epspolytope. if val<0, then this value is ignored 
    [type.opt.validatelowerbound,parg] = parsem( {'validatelowerbound','lb'}, parg, +inf, 'expect',{'clcl',[-inf inf]} );   %Algorithm terminates if lower estimate of JSR is greater than validatelowerbound.    
    [val,parg] = parsem( {'validatebound','bound','b'}, parg, [], 'expecte',{'clcl',[-inf inf]} );                           %Algorithm terminates if lower estimate of JSR is greater than validatelowerbound.    
    if( ~isempty(val)); 
        type.opt.validateupperbound=val; 
        type.opt.validatelowerbound=val; end;
    
    %Output Options
    [type.opt.plot,parg] =               parsem({'plot','p'}, parg, 'none' );                  %what to plot. See tjsr_plotoutput for help.
    [type.opt.verbose,parg] =            parsem({'verbose','v'}, parg, 1, 'expect',@isnumeric );                    %defines the verbose-level: -1=totally no output, 0=no print output (except severe warnings), 1=little output, 2=standard amount of output, 3=much output
    [type.opt.diary,parg] =              parsem( 'diary', parg );                             %indicates if diary is written. Logic for that is done by tjsr, not by tjsr_worker!    
    [type.opt.profile,parg] =            parsem( 'profile', parg );                           %indicates if we are profiling. Logic for that is done by tjsr, not by tjsr_worker!
    [type.opt.save,parg] =               parsem({'save','s'}, parg, tif(type.opt.diary,2,0), 'expect',{0 1 2 3} ); %how much output (diary, plots, variables) shall be saved: 0= do not save output, 1=save at end, 2=save in between, 3=save in between with different filenames
        
    
    %Constants and Debug-options, shall not be changed by user
    [type.opt.balancing,parg] =          parsem( 'balancing', parg );                         %Flag which indicates that we are balancing.
    [type.opt.memory,parg] =             parsem( 'memory', parg );                            %Sets parameters in a way such that memory is saved.
    [type.opt.alwaysout,parg] =          parsem( 'alwaysout', parg );                         %Assumes a new point is always out, by setting all computed norms to infinity, also the points from the root and nearvertices.    
    [type.opt.epsequal,parg] =           parsem( 'epsequal', parg, 1e-12, 'expect',{'clcl',[-1 1e-3]} );                   %epsilon which is used to compare to floating numbers for equality
    [type.opt.epslinprog,parg] =         parsem( 'epslinprog', parg, 1e-9, 'expect',{'clcl',[-1 1e-3]} );                  %epsilon used for the linear programming part. If Gurobi is used, this value is fixed to 1e-9. If matlab_linprog is used, this value must be greaterequal 1e-10.
    [type.opt.waitafterbalancing,parg] = parsem( 'waitafterbalancing', parg );                %Waits for key-pressing after balancing    
    [type.opt.rholeqval,parg] =          parsem( 'rholeqval', parg, 0, 'expect',{0 1} );                        %Discards all matrix products whos spectral radius is greater than 1+10*eps, and assumes their resp. vertices are inside of the polytope
    [type.opt.showquantity,parg] =       parsem( 'showquantity', parg, type.opt.verbose*25, 'expect',@isnumeric ); %up to which size shall sets of vectors,matrices,etc. be displayed
    [type.opt.solver,parg] =             parsem( 'solver', parg, 'auto', 'expecte',{'auto','gurobi','matlab','a','g','m'} ); %which solver to use
    
    %must be the very last
    [~,parg] =                           parsem( 'help', parg ); %remove 'help'
    
    
    
    %Post-preprocessing
    type.info.infotext = vprintf( 'Arguments: %v\n', type.arguments, 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
    if( ~isempty(parg) ); 
        type.info.infotext = vprintf('Unkown options in argument: %v\nPress any key to continue.\n', parg, 'cpr','err', 'imp',[0 type.opt.verbose], 'str',type.info.infotext );
        type.unkownarguments = parg; 
        pause; end;
    
    type.info.infotext = vprintf( 'Using solver: %s\n', type.opt.solver, 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
    

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 