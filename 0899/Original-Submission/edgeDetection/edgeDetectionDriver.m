% EDGEDETECTIONDRIVER driver for edge detection functions
%
% Functions called:
%    1) edgeDetectChebyshev.m, 2) edgeDetectFourier.m
% Last modified: July 10, 2008

function edgeDetectionDriver(h,ha)

         polyCh = get(ha.chebyButton,'Value'); 
         JCRIT = str2double(get(ha.edgeDetectCriticalJ,'String'));
             Q = str2double(get(ha.Q    ,'String'));
            NE = str2double(get(ha.nonlinearEnhancementParam   ,'String'));
         concentrationChoice = get(ha.linearConcentrationFactorButton,'Value');

        switch polyCh
            case 1
              [S,uE,uN] = edgeDetectChebyshev(ha.ak,JCRIT,Q,NE,concentrationChoice);
            otherwise
              [S,uE,uN] = edgeDetectFourier(ha.ak,JCRIT,Q,NE,concentrationChoice);
        end
         

   % NOTE: all Chebyshev pseudospectral pde examples have been computed on a mapped CGL grid

           T = ['Edges found at: '];
           
           if get(ha.chebyButton,'Value') & ha.isPdeSolution
              A = num2str(asin(ha.alpha*S)/asin(ha.alpha), '\n%1.8f');
           else
             A = num2str(S, '\n%1.8f');       % Convert to string array.
           end
        
           TA = [T, A];
           set(ha.outputTextBox,'string',TA);
         
         axes(ha.axes2)
         ind = find(uN>0); 

         if (ha.isPdeSolution & get(ha.chebyButton,'Value'))
            plot(ha.xa,uE,'g',ha.xa,uN,'b',asin(ha.alpha*ha.x(ind))/asin(ha.alpha),uN(ind),'k.')
         else 
           plot(ha.x,uE,'g',ha.x,uN,'b',ha.x(ind),uN(ind),'k.')
         end
         if strcmp(get(ha.displayLegends, 'Checked'), 'on'), legend('raw edge','enhanced edge'); end
   
         ha.edgeLocations = S;
         ha.edgeSeries    = uE;
         ha.edgeSeriesNL  = uN;     
         guidata(h,ha);

  
