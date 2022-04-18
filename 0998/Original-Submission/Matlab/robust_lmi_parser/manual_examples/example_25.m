clear;

Ai = rolmipvar(2,2,'A','full',[2 2],[1 1]);
Pi = rolmipvar(2,2,'P','sym',2,1);

fprintf(1,'explict: \n');
out = texify(Ai'*Pi + Pi*Ai,'explicit',Ai,Pi,{'\alpha','\beta'})

fprintf(1,'explict: \n');
out = texify(Ai'*Pi + Pi*Ai,'implicit',Ai,Pi)

fprintf(1,'polynomial: \n');
out = texify(Ai'*Pi + Pi*Ai,'polynomial',Ai,Pi,{'\alpha','\beta'})