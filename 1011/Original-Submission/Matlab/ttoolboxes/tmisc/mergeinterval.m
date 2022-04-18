function [varargout] = mergeinterval(varargin)
% [ X ] = mergeinterval( left, right, [options] )
% [ X ] = mergeinterval( v, [options] )
% [ X ] = mergeinterval( C, [options] )
% [ X ] = mergeinterval( M, [options] )
% Merges intervals
%
% Input: 
%   left            1xN vector, left bounds of the intervals
%   right           1xN vector, right bounds of the intervals
%       OR
%   v               1x2N vector, odd entries are the left bounds of the intervals, even entries are the right bounds of the intervals
%       OR
%   C               1xN cell array of 1x2 vector, each cell is one interval
%       OR
%   M               Nx2 matrix, each row is one interval
%
% Option:
%   'switch'        By default, intervals where left>right are dismissed
%                       If this flag is set, those intervals are switched
%   'output',str    default: depends on input
%                       'l'         left/right format for output
%                       'v'         left/right format for output
%                       'C'         left/right format for output
%                       'M'         left/right format for output
%
% Output:           
%   X               default: same format as input.
%                   
%
% Info:
%   The set union{Ii) can be written as a canonical partition by
%       intervals Jk; i.e., union{Ii) = union(Jk), where Jk are M intervals
%       (with M<=N, so the partition is minimum cardinal), and {Jk} are
%       disjoint to each other (their intersections are empty). This function
%       returns Jk = [lower(k),upper(k)], k=1,2,...M, in the ascending sorted
%       order.
%   Algorithm complexity: O(N*log(N))
%
% E.g.: [lower, upper] = mergeinterval([0 1 2 3 4],[1.5 1.6 3.5 3 5])
%           %lower =   0    2  4
%           %upper = 1.6  3.5  5
%       v = mergeinterval([0 1.5 1 1.6 2 3.5 3 3 4 5])
%       C = mergeinterval({[0 1.5];[1 1.6];[2 3.5];[3 3];[4 5]}); vdisp(C)
%       M = mergeinterval([0 1.5; 1 1.6; 2 3.5; 3 3; 4 5])
%       C = mergeinterval([0 1.5; 1 1.6; 2 3.5; 3 3; 4 5],'output','C')
%
%
% Written by: Bruno Luong <brunoluong@yahoo.com>, Original: 25-May-2009
% Changed by: tommsch, 2018

%#ok<*ALIGN> 

    %Parse Input:
    %%%%%%%%%%%%%%%
    [switchflag,varargin]=parsem({'switch','s'},varargin);
    [outputfl,varargin]=parsem({'output','o'},varargin,[]);
    
    if(size(varargin,2)==2);
        if(isempty(outputfl)); 
            outputfl='l'; end; %do nothing %case left,right
        left=varargin{1};
        right=varargin{2};
        if(~isequal(numel(left),numel(right)));
            error( 'mergeinterval:format', 'Wrong input format.'); end
    elseif(size(varargin,2)==0);
        if(isempty(outputfl)); outputfl='l'; end; %do nothing %case left,right
        left=[];
        right=[];
    else
        val=varargin{1};
        if(iscell(val)); %case C %val={[0 3],[4 2],[10 10]};         
            if(isempty(outputfl)); 
                outputfl='C'; end; 
            val=[val{:}];    
        elseif(isrow(val) || isempty(val));  %val=[0 3 4 2 10 10].'; 
            if(isempty(outputfl)); 
                outputfl='v'; end;
            if(mod(numel(val),2)); 
                error( 'mergeinterval:format', 'Wrong input format.'); end
        elseif(ismatrix(val)); %case M %val=[0 3; 4 2; 10 10]; 
            if(isempty(outputfl)); 
                outputfl='M'; end; 
            if(~isequal(size(val,2),2)); 
                error('Wrong input format.'); end;
            val=val.'; 
        else;                            
            error( 'mergeinterval:format', 'Wrong input format.'); end
        val=reshape(val,2,[]);
        left=val(1,:);
        right=val(2,:);
    end

    if(switchflag)
        %switch places where the intervals number are switched
        idx = find(right<left);
        val=left(idx);
        left(idx)=right(idx);
        right(idx)=val;
    end
    % Detect when right < left (empty Ii), and later remove it
    notempty = find(right>=left);
    % sort the rest by left bound
    [left, iorder] = sort(left(notempty));
    right = right(notempty(iorder));

    % Allocate, as we don't know yet the size, we assume the largest case
    if(issym(left))
        lower = sym(zeros(size(left)));
        upper = sym(zeros(size(right)));
    else
        lower = zeros(size(left));
        upper = zeros(size(right));
    end
    % Nothing to do
    if(~isempty(lower))
        % Initialize
        l = left(1);
        u = right(1);
        k = 0;
        % Loop on brakets
        for i=1:length(left)
            if left(i) > u % new Jk detected
                % Stack the old one
                k = k+1;
                lower(k) = l;
                upper(k) = u;
                % Reset l and u
                l = left(i);
                u = right(i);
            else
                u = max(u, right(i));
            end
        end % FOR loop
        % Stack the last one
        k = k+1;
        lower(k) = l;
        upper(k) = u;
        % Remove the tails
        lower(k+1:end) = [];
        upper(k+1:end) = [];
    end
    
    %lower=[0 3 10]; upper=[1 5 10];
    %if(nargout<1); error('Function needs at least one output.'); end;
    switch outputfl
        case 'l'
            if(nargout<2);
                warning( 'mergeinterval:nargout', 'Function needs two outputs.' ); end;
            varargout{1}=lower;
            varargout{2}=upper;
        case 'C'
            varargout{1}=num2cell([lower;upper].',2).';
        case 'v'
            varargout{1}=reshape([lower; upper],1,[]);
        case 'M'
            varargout{1}=[lower; upper].';
        otherwise
            error( 'mergeinterval:format', 'Allowed ''output''-strings: ''l'', ''C'', ''v'', ''M''.' ); end
end 
