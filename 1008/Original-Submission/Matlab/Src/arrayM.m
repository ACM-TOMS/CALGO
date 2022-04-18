 function [arrayform] = arrayM(inputMatrix)
 
 % This function just extracts the array form of the multicomplex number
 % from its matrix form. It is used by several methods os the class but
 % should be defined outside of it as it doesnt take multicomplex inputs. 
 
 arrayform=transpose(inputMatrix(:,1)); 
 
 end
 
 
 

