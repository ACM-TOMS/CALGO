function[currentTheta] =  NewtonRaphson(initialTheta,meanvjs,upperTheta,maxit,ThetaFunction,ThetaDevFunction,newton_params)

%%function d = rtsafe(w,b)
%Modified from 'Numerical recipes: the art of scientific computing'
%
% Partially taken from: http://www.mathworks.com/matlabcentral/fileexchange/35977-logicle-histogram/content/logicleHist.m
%
% Created version 01/16/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

%the function
gFunction = ThetaFunction;%@ThetaFunction;
derivFunction = ThetaDevFunction;%@ThetaDevFunction;

X_ACCURACY = 0.001;
MAX_IT = maxit; %100

lowerLimit = initialTheta;
upperLimit = upperTheta;

%if ((gFunction(lowerLimit)>0) && (gFunction(upperLimit)>0)) || ...
%        ((gFunction(lowerLimit)<0) && (gFunction(upperLimit)<0))
%    error('Root must be bracketed');
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (gFunction(lowerLimit,meanvjs,newton_params) == 0)
    currentTheta = initialTheta;
    return;
else if(gFunction(upperLimit,meanvjs,newton_params) == 0)
	currentTheta = upperTheta;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if gFunction(lowerLimit,meanvjs,newton_params)<0
   xLow = lowerLimit;
    xHigh = upperLimit;
else
    xLow = upperLimit;
    xHigh = lowerLimit;
end

%root = (lowerLimit + upperLimit)/2;
root = initialTheta; %<-fijamos un valor inicial para el x.

dxOld = abs(upperLimit - lowerLimit);
dx = dxOld;
g = gFunction(root,meanvjs,newton_params);
dg = derivFunction(root,newton_params);

for i = 0:MAX_IT
    if (((root-xHigh)*dg-g)*((root-xLow)*dg-g) > 0) || ...
            (abs(2*g) > abs(dxOld*dg))
        %Bisect method
        dxOld = dx;
        dx = (xHigh-xLow)/2;
        root = xLow + dx;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (xLow == root) % codigo Josu
            currentTheta = root;
            if (currentTheta < initialTheta)
                currentTheta = initialTheta
            end
            if (currentTheta > upperTheta)
                currentTheta = upperTheta
            end

            return;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %Newton method
        dxOld = dx;
        dx = g/dg;
        temp = root; % codigo Josu
        root = root - dx;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (temp == root) % codigo Josu
            currentTheta = root;
            if (currentTheta < initialTheta)
                currentTheta = initialTheta
            end
            if (currentTheta > upperTheta)
                currentTheta = upperTheta
            end
            return;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    if (abs(dx) < X_ACCURACY)
        currentTheta = root; %%d = root;
        if (currentTheta < initialTheta)
            currentTheta = initialTheta
        end
        if (currentTheta > upperTheta)
            currentTheta = upperTheta
        end
        return; 
    end
    
    g = gFunction(root,meanvjs,newton_params);
    dg = derivFunction(root,newton_params);
    
    if (g < 0)
        xLow = root;
    else
        xHigh = root;
    end
end

currentTheta = (upperTheta+initialTheta)/2;
return%%error('Maximum number of interations exceeded in rtsafe');

end