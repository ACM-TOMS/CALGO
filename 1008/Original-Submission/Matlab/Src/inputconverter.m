function [out] = inputconverter(real,in,h)

% included in order to make it easier for the user to create a multicomplex 
% number with the step sizes in the correct coefficients. This function 
% inputs the real part, the step size value, and an array with the imaginary 
% component number where you want the step sizes to be in. For example, 
% inputconverter(20,[1,2,3],10^{-10}) creates the multicomplex number 
% 20+10^{-10}i_1+10^{-10}i_2+10^{-10}i_3.

im=max(in);
outp=zeros(1,2^im);
outp(1)=real;

for w=1:length(in)
    
outp((2^(in(w)-1))+1)=h;

end

out=outp;

end

