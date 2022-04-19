function [wz,dVdx,dVdy] = wz(z)
one_sqrt_pi=0.5641895835477563;
i_sqrt_pi=1i*one_sqrt_pi;
wz=NaN(size(z));
z_abs_sqr=(real(z).^2+imag(z).^2);
%====================Region I ----> Laplace Cont. Fractions, 1st convergent
index_R1=(z_abs_sqr>=38000.0);
if~isempty(index_R1);
    t=z(index_R1);
    wz(index_R1)=(i_sqrt_pi./t);
end
%====================Region II ----> Laplace Cont. Fractions, 2 convergents
% Also identical to the approximation for regions I in Humlicek’s w4 algorithm
index_R2=(38000.0>z_abs_sqr & z_abs_sqr>=256.0);
if~isempty(index_R2);
    t=z(index_R2);
    wz(index_R2)=i_sqrt_pi.*t./(t.*t-0.5);
end
%====================Region III----> Laplace Cont. Fractions, 3 convergents
index_R3=(256.0>z_abs_sqr & z_abs_sqr>=62.0);
if~isempty(index_R3);
    t=z(index_R3);
    wz(index_R3)=(i_sqrt_pi./t).*(1+0.5./(t.*t-1.5));
end
%====================Region IV----> Laplace Cont. Fractions, 4 convergents
% Also identical to the approximation for regions II in Humlicek’s w4 algorithm
index_R4=(62.0>z_abs_sqr & z_abs_sqr>=30.0 & imag(z).^2>=1e-13 );
if~isempty(index_R4);
    t=z(index_R4);
    tt=t.*t;
    wz(index_R4)=(i_sqrt_pi.*t).*(tt - 2.5)./(tt.*(tt - 3.0)+0.75);
end
%=====================Region V----> Humlicek's w4 (Region IV)
index_R5=(62.0>z_abs_sqr & ~index_R4 & (z_abs_sqr>2.5  &  imag(z).^2<0.072));
if~isempty(index_R5);
    t=z(index_R5);
    u=-t.*t;
    wz(index_R5)=exp(u)+1i*t.*(36183.31-u.*(3321.99-u.*(1540.787-u.*...
        (219.031-u.*(35.7668-u.*(1.320522-u.*one_sqrt_pi))))))./...
        (32066.6-u.*(24322.84-u.*(9022.228-u.*(2186.181-u.*...
        (364.2191-u.*(61.57037-u.*(1.841439-u)))))));
end
%=====================Region VI---->Hui's p-6 Approximation
index_R6=(30.0>z_abs_sqr & ~index_R5);
if~isempty(index_R6);
    t3=-1i*z(index_R6 );
    wz(index_R6)= ((((((one_sqrt_pi*t3+5.9126262).*t3+30.180142).*...
        t3+93.15558).*t3+181.92853).*t3+214.38239).*t3+122.60793)./...
        (((((((t3+10.479857).*t3+53.992907).*t3+170.35400).*t3+348.70392).*...
        t3+457.33448).*t3+352.73063).*t3+122.60793);
end
%===================Calculations of partial derivatives
dVdx=-2*real(z.*wz);  % Partial derivative of real(w) w.r.t. x
dVdy=2*imag(z.*wz)- 1.12837916709551; % Partial derivative of real(w) w.r.t. y
%dLdx=-dVdy; dLdy=dVdx; % Partial derivatives of imag(w) w.r.t. x&y
