function [double, onto6, trst6,smooth,inmsh] = osc(problem,nmsh,incmp,defcor,dfexmx,ratdc,trst6,smooth,double,onto6,inmsh)
%
%   Private function for twpbvpc
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
jsndif = 0;
rmax = 0;
ninter = nmsh - 1;
if problem.debug,  disp('osc'), end
allsum = 0;
smlsum = 0;
bigsum = 0;
ibig = 0;
ism = 0;


for im = 1 : ninter
    abdef = abs(defcor(incmp,im));
    allsum = allsum + abdef;
    if (abdef < 0.5 * dfexmx)
        ism = ism + 1;
        smlsum = smlsum + abdef;
    else
        ibig = ibig + 1;
        bigsum = bigsum + abdef;
    end

    if problem.debug
        disp(sprintf('im %g ratdc %g abdef %g val %g', im, ratdc(im), abdef, 1e-2*dfexmx))
    end
        if (ratdc(im) < 0 && abdef >= 1e-2*dfexmx)
            jsndif = jsndif + 1;

            if (jsndif > 4)
                onto6 = 0;
                double = 1;
                return
            end
            if (abs(ratdc(im))>= rmax)
                rmax = abs(ratdc(im));
                inmsh = im;
            end
        end
end

        if problem.debug  
            disp(sprintf('rmax %g jsndif %g',rmax,jsndif)) 
        end
           
            avsm = 0;
            if (ism > 0), avsm = smlsum /ism;end
             avbg = 0;
             if (ibig > 0), avbg = bigsum /ibig;end
             ave = allsum/ninter;

              if problem.debug
                  disp(sprintf('ave %g  avsm %g avbg %g',ave,avsm,avbg)) 
              end
                         
                if (avsm > avbg*0.1  ||  ave > 0.5*avbg)
                    onto6 = 1;
               elseif (jsndif==0)
                      smooth = 1;
                      onto6 = 1;
               else
                       double = 0;
                       onto6  = 0;
                       trst6  = 0;
                end
                return
end
                
 
