% TEXTEDGES  Allows edges to be selected by entering them in the Edges text box.
%            The endpoints of the interval do not need to be enterred.
% 
% Called by: 
%   1) mpt.tex
% Last modified: October 17, 2007

function textEdges(h,ha)

     S(1) = -1;
     sTemp = str2num(get(ha.edgesEditText,'String'));            % e.g., enter -0.5 0 0.5  for a 1x3 matrix
                                                                 %       -1 < S(i) < 1
     n = 1;
     while n <= length(sTemp)
       S(n+1) = sTemp(n);
       n = n + 1;
     end

     S(end+1) = 1;
     if (ha.isPdeSolution & get(ha.chebyButton,'Value')), S = sin( S*asin(ha.alpha) )/ha.alpha;  end
     ha.edgeLocations = S;
     
      T = ['Edges marked at: '];
      if get(ha.chebyButton,'Value') & ha.isPdeSolution
        A = num2str(asin(ha.alpha*S)/asin(ha.alpha), '\n%1.8f');
      else
        A = num2str(S, '\n%1.8f');       % Convert to string array.
      end
      TA = [T, A];

      set(ha.outputTextBox,'string',TA);
      guidata(h,ha);
     
     
