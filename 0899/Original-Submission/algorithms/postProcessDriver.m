% POSTPROCESSDRIVER  driver for 1d postprocessing functions
%
% Called by:
%   1) mpt.m
% Functions called:
%   1) filterChebyshev.m, 2) filterFourier.m, 3) grp.m,  4) inverseReprojection.m, 5) frp.m,
%   6) chebyshevPade.m, 7) fourierPade.m, 8) digitalTotalVariationFilter.m, 
%   9) digitalTotalVariationFilterPeriodic.m
% Last modified: October 17, 2007

function postProcessDriver(h,ha)

      tic
      polyCh = get(ha.chebyButton,'Value');         % 1 Chebyshev, otherwise 0 and Fourier

   % ----------------------------------------------------------------------------
   % ---------------------- spectral filtering ----------------------------------
   % ----------------------------------------------------------------------------


       if get(ha.filterButton,'Value') == 1

            if get(ha.expFilterButton,'Value') == 1
               fCh = 1;
            elseif get(ha.vandevenFilterButton,'Value') == 1
               fCh = 3;
            else
               fCh = 2;
            end

            filterOrder = str2double(get(ha.spectralFilterOrder,'String'));

            switch polyCh
                case 1
                    if ~ha.isPdeSolution
                       up = filterChebyshev(ha.ak,ha.xa,fCh,filterOrder);
                    else
                       up = filterChebyshev(ha.ak,ha.x,fCh,filterOrder);             % note: gInverse(xm) = xGL, all CPS PDE example on mapped grid
                    end
                otherwise
                     up = filterFourier(ha.ak,ha.xa,fCh,filterOrder);
            end


   % ----------------------------------------------------------------------------
   % ---------------------- Gegenbauer Reconstruction ---------------------------
   % ----------------------------------------------------------------------------

     elseif get(ha.grpButton,'Value') == 1

           LK = str2num(get(ha.gegenbauerL,'String'));
           MK = str2num(get(ha.gegenbauerM,'String'));

           if get(ha.projectButton,'Value') == 1 & ~ha.isPdeSolution
             [up,gh] = grp(ha.edgeLocations,LK,MK,ha.exampleHandle,ha.xa,polyCh);
           elseif ~ha.isPdeSolution
              [up,gh] = grp(ha.edgeLocations,LK,MK,ha.ak,ha.xa,polyCh);
           else
             [up,gh] = grp(ha.edgeLocations,LK,MK,ha.ak,ha.x,polyCh);
           end
       
   % ----------------------------------------------------------------------------
   % ------------- Inverse Gegenbauer Reconstruction ---------------------------
   % ----------------------------------------------------------------------------

      elseif get(ha.igrpButton,'Value') == 1

          MK = str2num(get(ha.gegenbauerM,'String'));
          if get(ha.projectButton,'Value') == 1 & ~ha.isPdeSolution
             [up,condW] = inverseReprojection(ha.edgeLocations,MK,ha.xa,ha.exampleHandle);
          else
             [up,condW] = inverseReprojection(ha.edgeLocations,MK,ha.xa,ha.ak);    % the underlying function is not known
          end

   % ----------------------------------------------------------------------------
   % ---------------------- Freud Reconstruction ---------------------------
   % ----------------------------------------------------------------------------

      elseif get(ha.freudButton,'Value') == 1
      
         
           
           if get(ha.projectButton,'Value') == 1 & ~ha.isPdeSolution
             [up,gh,Lambda,m] = frp(ha.edgeLocations,ha.exampleHandle,ha.xa,polyCh,length(ha.ak));
           elseif ~ha.isPdeSolution
             [up,gh,Lambda,m] = frp(ha.edgeLocations,ha.ak,ha.xa,polyCh,length(ha.ak));
           else
              [up,gh,Lambda,m] = frp(ha.edgeLocations,ha.ak,ha.x,polyCh,length(ha.ak));
           end
           
           set(ha.gegenbauerL,'String',num2str(Lambda));
           set(ha.gegenbauerM,'String',num2str(m));

   % ----------------------------------------------------------------------------
   % ---------------- Generalized Pade approximation ----------------------------
   % ----------------------------------------------------------------------------

        elseif get(ha.padeButton,'Value') == 1

            M = str2double(get(ha.padeM,'String'));
           Nc = str2double(get(ha.PadeNc,'String'));

           switch polyCh
             case 1
               if ~ha.isPdeSolution
                  up = chebyshevPade(ha.ak,M,Nc,ha.xa)';
               else
                  up = chebyshevPade(ha.ak,M,Nc,ha.x)';
               end
             otherwise
                  up = fourierPade(ha.ak,M,Nc,ha.xa);
           end


   % ----------------------------------------------------------------------------
   % ---------------- Digitial Total Variation (DTV) Filtering ------------------
   % ----------------------------------------------------------------------------

         elseif get(ha.dtvFilterButton,'Value') == 1
           iterations = str2double(get(ha.dtvIterations,'String'));          % actually, time steps
           lambda = str2double(get(ha.dtvLambda,'String'));

           switch polyCh
              case 1
                  up = digitalTotalVariationFilter(ha.fa,lambda,iterations);
              otherwise
                  up = digitalTotalVariationFilterPeriodic(ha.fa,lambda,iterations);
           end
           
         end
         
         str = ['Elapsed time:   ', num2str(toc, '\n%4.5f'), ' s'];
         set(ha.outputTextBox,'string',str);

  % ---- plot postprocessed approximation in top image ---------------------------

         ha.fPost = up;
         guidata(h,ha);
         axes(ha.axes1)

         plot(ha.xa,ha.fa,'g',ha.xa,ha.fExact,'r',ha.xa,ha.fPost,'b')
         if strcmp(get(ha.displayLegends, 'Checked'), 'on'), legend('numerical','exact','postprocessed'); end

 % plot the error of the postprocessed approximation in bottom image

      axes(ha.axes2)
      semilogy(ha.xa(2:end-1),abs( ha.fExact(2:end-1) - ha.fa(2:end-1) ),'g',ha.xa(2:end-1),abs( ha.fExact(2:end-1) - ha.fPost(2:end-1) ),'r')
      if strcmp(get(ha.displayLegends, 'Checked'), 'on'), legend('error numerical','error postprocessed'); end
