function bb(file_name)
% generate the c code with command mcc -m bb

fid = fopen ( file_name,'r');
X=fscanf(fid,'%f%f');
fclose(fid);


if length(X) ~= 5
   Z = [1e20 1e20 1e20];
else
   Z = [X(5) 0 0];
   for i = 1:5
     Z(2) = Z(2) + (X(i)-1).^2;
     Z(3) = Z(2) + (X(i)+1).^2;
   end
   Z(2) = Z(2) - 25;
   Z(3) = 25 - Z(3);
end

disp(Z);


