% MANUALEDGES  Allows edge locations to be specified by mouse clicks.  From the MPT menu
%              select Options -> Manual Edges.  A crosshair appears and edges are selected
%              by a left mouse click.  The last edgde is selected with a right mouse click
%              which ends the selection process.  The selected edges appear in the text output
%              box of the GUI.  The endpoints of the interval do not need to be marked.
% 
% Called by: 
%   1) mpt.tex
% Last modified: October 17, 2007


function manualEdges(h,ha) 


  xm = ha.x;
  Np = str2double(get(ha.N,'String'));

  S(1) = -1.0;
  button = 1;
  n = 1;
  while button == 1
          [xi,yi,button] = ginput(1);
    
          edgeIndex=1;
          for i=1:Np-1                          % find closest grid pt to mouse click x
            if xi >= xm(i) & xi <= xm(i+1)
              edgeIndex = i;   
            end      
          end
         
         S(n+1) = xm(edgeIndex);                % add x as the location    
         n = n+1;
    end
    S(end+1) = 1.0;
    
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
   
