function [onto6,smooth,callrt,strctr,oscchk,double,reposs]=decid4(problem,linear,rat1,rat2,dfexmx,dfimmx,derivm,dfold,tolval,oldrt1)
%
%   Private function for twpbvpc
%
%  
%
%       Authors:
%
%       Jeff R. Cash 
%            (Department of Mathematics, Imperial College,  London, England.)
%       Davy  Hollevoet 
%            (Vakgroep Toegepaste Wiskunde en Informatica, Universiteit Gent, Belgium.)
%       Francesca Mazzia  
%            (Dipartimento di Matematica, Universita' di Bari, Italy)
%       Abdelhameed Nagy Abdo
%            (Dipartimento di Matematica, Universit\`a di Bari, Italy)
%            (Dept. of Mathematics, Faculty of Sciences, Benha  University,Egypt)
%            
%
    rtst=50;
    derval=50;
    
    onto6=1;
    callrt=0;
    smooth=0;
    oscchk=0;
    strctr=0;
    reposs=0;
    double=0;
    
    stest=1;
    if linear
        stest=dfexmx<0.1*dfold;
    end
    
    if rat2<rtst
        if stest
            smooth=1;
            if problem.debug, disp('smooth'); end
        else
            oscchk=1;
        end
        return
    end
    
    thttol=32*tolval;
    
    if rat1<rtst && dfexmx<thttol
        if stest
            smooth=1;
        else
            oscchk=1;
        end
        return
    end
    
    if rat1<rtst && dfexmx>=thttol
        callrt=1;
        return;
    end
    
    if derivm > derval && dfexmx < thttol
         if stest            % new part from if ....end
            smooth = 1;
         else
            oscchk = 1;
         end

         return
    end                       % 
    if derivm>derval && dfexmx>thttol
        if dfimmx<1
            callrt=1;
        else
            strctr=1;
            if linear
                onto6=0;
                if 2*rat1>=oldrt1
                    double=1;
                end
            end
        end
        return
    end

    
    if linear
        reposs=1;
    end
    
end