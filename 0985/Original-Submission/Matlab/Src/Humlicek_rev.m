
function [ww,dVdx,dVdy] = Humlicek_rev (z)
% Humlicek (JQSRT, 1982) complex probability function: rational approximations
rSqrPi=0.56418958354775628;
ww = NaN(size(z));
t=-1i.*z;
s=abs(real(z))+imag(z);
index_1=(s>=15);
if~isempty(index_1);
ww(index_1)= t(index_1) .* rSqrPi ./ (0.5 + t(index_1).*t(index_1)); %
end;

% original code used: index_2=find(s>=5.5 & s<=15);
%  index_2=find(s>=5.5 & s<=15 & imag(z)>1e-6); %Reformed code uses
 index_2=(s>=5.5 & s<15 ); % Original code
if~isempty(index_2);
u = t(index_2).*t(index_2);
ww(index_2) = (t(index_2) .* (1.410474 + u.*rSqrPi))./ (.75 + (u .*(3.+u))) ;
end;
    
index_3= (s<5.5 & imag(z)>=0.195*abs(real(z))-0.176);
if~isempty(index_3);
ww(index_3)=  (16.4955 + t(index_3) .* (20.20933 + t(index_3) .*...
    (11.96482 + t(index_3) .* (3.778987 + 0.5642236.*t(index_3))))) ./...
    (16.4955 + t(index_3) .* (38.82363 + t(index_3) .* (39.27121 + t(index_3) .*...
    (21.69274 + t(index_3) .* (6.699398 + t(index_3))))));
end;

%Humlcek's Original
index_4=(s<5.5 & imag(z)<(0.195*abs(real(z))-0.176) );
%Zaghloul's correction
% index_4=find(s<5.5 & imag(z)<(0.195*abs(real(z))-0.176) | (s>=5.5 & s<=15 & imag(z)<=1e-6) ); 
if~isempty(index_4);
u=t(index_4).*t(index_4);
% nom=(t(index_4).*(36183.31-u.*(3321.99-u.*(1540.787-u.*(219.031-u.*(35.7668-u.*(1.320522-u.*.56419)))))));
% den=(32066.6-u.*(24322.8-u.*(9022.23-u.*(2186.18-u.*(364.219-u.*(61.5704-u.*(1.84144-u)))))));
ww(index_4)  = exp(u) - (t(index_4).*(36183.31-u.*(3321.99-u.*(1540.787-u.*(219.031-u.*(35.7668-u.*(1.320522-u.*.56419)))))))./...
    (32066.6-u.*(24322.8-u.*(9022.23-u.*(2186.18-u.*(364.219-u.*(61.5704-u.*(1.84144-u)))))));
end;
%-------Calculations of partial derivatives
 dVdx=-2*real(z.*ww);  % Partial derivative of real(w) w.r.t. x
 dVdy=2*imag(z.*ww)-1.128379167095513e+000; % Partial derivative of real(w) w.r.t. y
%dLdx=-dVdy; dLdy=dVdx; % Partial derivatives of imag(w) w.r.t. x&y


