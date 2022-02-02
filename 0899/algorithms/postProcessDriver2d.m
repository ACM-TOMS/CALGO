% POSTPROCESSDRIVER2d  driver for 2d postprocessing functions
%
% Called by:
%   1) mpt.m
% Functions called:
%   1) filterChebyshev2d.m, 2) filterFourier2d.m, 3) digitalTotalVariationFilter_2d.m,
%   4) digitalTotalVariationFilter_2d_8.m, 5) digitalTotalVariationFilterPeriodic_2d.m,
%   6) digitalTotalVariationFilterPeriodic_2d_8.m
% Last modified: October 17, 2007

function postProcessDriver2d(h,ha)

     tic
     polyCh = get(ha.chebyButton,'Value');

   % ----------------------------------------------------------------------------
   % ---------------- Spectral Filtering ----------------------------------------
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
          Nex = str2double(get(ha.N,'String'));

              switch polyCh
                 case 1

                    sigma = spectralFilter(length(ha.x),fCh,filterOrder,1);
                    up = filterChebyshev2d(ha.ak,ha.xa,sigma);

                 otherwise

                    sigma = spectralFilter(length(ha.x),fCh,filterOrder,0);
                    up = filterFourier2d(ha.ak,length(ha.xa),sigma);

                 end

       end % filtering


   % ----------------------------------------------------------------------------
   % ---------------- Digitial Total Variation (DTV) Filtering ------------------
   % ----------------------------------------------------------------------------

         if get(ha.dtvFilterButton,'Value') == 1
           iterations = str2double(get(ha.dtvIterations,'String'));
           lambda = str2double(get(ha.dtvLambda,'String'));
           
            switch polyCh
                 case 1
                    if get(ha.dtv4_radioButton,'Value') == 1 
                       up = digitalTotalVariationFilter_2d(ha.fa,lambda,iterations);
                    else   
                       up = digitalTotalVariationFilter_2d_8(ha.fa,lambda,iterations);
                    end
                 otherwise
                    if get(ha.dtv4_radioButton,'Value') == 1 
                       up = digitalTotalVariationFilterPeriodic_2d(ha.fa,lambda,iterations);
                    else
                       up = digitalTotalVariationFilterPeriodic_2d_8(ha.fa,lambda,iterations);
                    end
                 end
           
       
         end
         
         str = ['Elapsed time:   ', num2str(toc, '\n%4.5f'), ' s'];
         set(ha.outputTextBox,'string',str);


   % --------------------------------------------------------------------------------


   if strcmp(get(ha.countouMenuItem, 'Checked'), 'on')      % plot post-processed approximation in upper image
      contour(ha.axes1,ha.xa,ha.xa,up);                     % plot post-processed error in botton image
    else
      surf(ha.axes1,ha.xa,ha.xa,up);
    end

    if strcmp(get(ha.countouMenuItem, 'Checked'), 'on')
      contour(ha.axes2,ha.xa,ha.xa,abs(up-ha.fExact));
    else
      surf(ha.axes2,ha.xa,ha.xa,abs(up-ha.fExact));
    end

    ha.fPost = up;
    guidata(h,ha);
