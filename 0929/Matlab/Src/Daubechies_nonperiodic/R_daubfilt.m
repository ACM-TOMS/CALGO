function [lp]=R_daubfilt(k,D)
% Filter coefficient for right boundary scaling function
%Taken from the paper 'Wavelet on the interval and fast wavelet
%transform'' by Cohen daubechies and Vial''


M=D/2;  % Here M is number of vanishing moment

if(M==2)
    if(k==1)
        lp(5)=.4431490496;
        lp(4)=.7675566693;
        lp(3)=.3749553316;
        lp(2)=.1901514184;
        lp(1)=-.1942334074;
    elseif(k==0)
        lp(3)=.2303890438;
        lp(2)=.4348969980;
        lp(1)=.8705087534;
    end
elseif(M==3)
    if(k==0)
        lp(1)=0.9096849943E+00;
        lp(2)=0.3823606559E+00;
        lp(3)= 0.1509872153E+00;
        lp(4)=0.5896101069E-01;
    elseif(k==1)
        lp(1)= -0.2904078511E+00;
       lp(2)= 0.4189992290E+00;
       lp(3)= 0.4969643721E+00;
       lp(4)=0.4907578307E+00;
       lp(5)=0.4643627674E+00;
       lp(6)= 0.1914505442E+00;
   elseif(k==2)
    lp(1)=0.8183541840E-01;
   lp(2)=-0.1587582156E+00;
   lp(3)-0.9124735623E-01;
   lp(4)= 0.6042558204E-03;
   lp(5)= 0.7702933610E-01;
   lp(6)= 0.5200601778E+00;
    lp(7)=0.7642591993E+00;
    lp(8)=0.3150938230E+00;
end
    elseif(M==4)
        if(k==3)
            lp(11)=.32210279e-01;
            lp(10)=-.12598952e-01;
            lp(9)=-.99108040e-01;
            lp(8)=.29771110;
            lp(7)=.80394959;
            lp(6)=.49779209;
            lp(5)=-.30235885e-01;
            lp(4)=-.67659162e-01;
            lp(3)=-.17709184e-01;
            lp(2)=.19132441e-01;
            lp(1)=-.67756036e-02;
        elseif(k==2)
            lp(9)=.32148741e-01;
            lp(8)=-.12574882e-01;
            lp(7)=-.10276635;
            lp(6)=.29864734;
            lp(5)=.81641197;
            lp(4)=.46061686;
            lp(3)=.29213680e-01;
            lp(2)=-.13907160;
            lp(1)=.12900783e-01;
        elseif(k==1)
            lp(7)=.41268408e-01;
            lp(6)=-.16142013e-01;
            lp(5)=-.15813389;
            lp(4)=.39377582;
            lp(3)=.75400048;
            lp(2)=.44880018;
            lp(1)=-.21916264;
        elseif(k==0)
            lp(5)=.64379349e-01;
            lp(4)=-.25191808e-01;
            lp(3)=.59477713e-01;
            lp(2)=.39191428;
            lp(1)=.91547054;
        end
    elseif(M==5)
        if(k==0)
        lp(1)=0.6629994791E+00;
        lp(2)= 0.5590853114E+00;
        lp(3)=-0.3638895712E+00;
        lp(4)=-0.3190292420E+00;
       lp(5)= -0.8577516141E-01;
       lp(6)=    0.7938922949E-01;
   elseif(k==1)
   lp(1)=-0.1615976602E+00;
   lp(2)= 0.7429922432E+00;
    lp(3)=0.5771521705E+00;
    lp(4)=0.2773360273E+00;
    lp(5)=0.1530007177E-01;
   lp(6)=-0.1063942878E+00;
   lp(7)=-0.1216769629E-01;
   lp(8)= 0.1126647052E-01;
elseif(k==2)
   lp(1)=0.2576681838E+00;
   lp(2)=-0.3995181391E-01;
   lp(3)=0.1290658921E-01;
    lp(4)=0.4163708782E+00 ;
   lp(5)= 0.6469374504E+00 ;
   lp(6)= 0.5608421202E+00 ;
   lp(7)= 0.4995700940E-02;
   lp(8)=-0.1569482750E+00;
   lp(9)=-0.2009485730E-01;
   lp(10)= 0.1860648984E-01;
elseif(k==3)
   lp(1)=-0.8955158953E-01;
   lp(2)=0.3895616532E-01;   
   lp(3)= 0.4751157158E-02;
   lp(4)=-0.6596188737E-01;
   lp(5)= 0.1520881337E-01;
   lp(6)= 0.2466644601E+00;
   lp(7)= 0.7199998123E+00;
   lp(8)= 0.6131240289E+00;
   lp(9)= 0.1537319456E-01;
   lp(10)=-0.1721947543E+00;
   lp(11)=-0.2083858871E-01;
   lp(12)= 0.1929513523E-01;
elseif(k==4)
    lp(1)= 0.1546583496E-01;
     lp(2)=-0.6559247482E-02;
    lp(3)=0.2728093541E-02 ;
    lp(4)=0.1830645747E-01;
    lp(5)=0.1708296448E-01;
    lp(6)=0.2156026058E-01;
    lp(7)=-0.3692524465E-01;  
    lp(8)=0.2045334663E+00;
    lp(9)=0.7233635779E+00;      
    lp(10)=0.6327590443E+00;
    lp(11)=0.1649999072E-01;
    lp(12)-0.1751674471E+00;
    lp(13)-0.2109311508E-01;
    lp(14)=0.1953080958E-01;
end
    end
            