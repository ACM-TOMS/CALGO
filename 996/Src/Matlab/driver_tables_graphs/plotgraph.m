%
% Input POPtable1.txt and POPtable2.txt
% Output figure2.eps -- Fig 2 on page 21 of TOMS's paper
%        figure3.eps -- Fig 3 on page 21 of TOMS's paper
%

% Choose table = 1 for Fig 2 on page 21 of TOMS's paper
% Choose table = 2 for Fig 3 on page 21 of TOMS's paper
table = 2;
%

if (table==1)
   fid = fopen('POPtable1.txt','r');
   ss = '$I_{\mathrm{bin}}=\{1,\ldots,n\}$,\quad   ';
   figure(2);
else
   fid = fopen('POPtable2.txt','r');
   ss = '$I_{\mathrm{box}}=\{1,\ldots,n\}$,\quad    ';
   figure(3);
end

cnt = 0; 
while (true)
     tline = fgetl(fid);
     if ~ischar(tline), break, end
     idx = find(tline=='&');
     if ~isempty(idx)
        cnt = cnt+1;
        L1a(cnt) = str2num(tline(idx(3)+1:idx(4)-1)); %BBC
        T1a(cnt) = str2num(tline(idx(4)+1:idx(5)-1)); %BBC
        L2a(cnt) = str2num(tline(idx(8)+1:idx(9)-1)); %NAL
        T2a(cnt) = str2num(tline(idx(9)+1:idx(10)-1)); %NAL  
        L1b(cnt) = str2num(tline(idx(12)+1:idx(13)-1)); %BBC
        T1b(cnt) = str2num(tline(idx(13)+1:idx(14)-1)); %BBC
        L2b(cnt) = str2num(tline(idx(17)+1:idx(18)-1)); %NAL
        T2b(cnt) = str2num(tline(idx(18)+1:idx(19)-1)); %NAL
     end
     disp(tline)
end
fclose(fid);

fontsize = 20;
linewidth = 1;
idxa = find(abs(L2a)>0);
option = 1;
if (option==1)
   LNAL = L2a(idxa); 
   LBBC = L1a(idxa); 
   TNAL = T2a(idxa); 
   TBBC = T1a(idxa);    
end
subplot(221)
xx = TNAL;
[xx,idx] = sort(xx);
LNAL = LNAL(idx); LBBC = LBBC(idx);
TNAL = TNAL(idx); TBBC = TBBC(idx);
Lratio = (LBBC-LNAL)./abs(LBBC);
Tratio = TNAL./TBBC;
xmin = exp(log(10)*(log10(min(xx))));
xmax = exp(log(10)*(log10(max(xx))));
hh = semilogx(xx,Lratio,'*'); hold on;
semilogx(xx,Lratio,'LineWidth', linewidth); hold off
grid
if (table==1)
   axis([xmin xmax -0.01 0.1]) 
else
   axis([xmin xmax -0.1 1])   
end
xlabel('SDPNAL+ time')
legend(hh,'LBv gap: (BBCPOP-SDPNAL+)/|BBCPOP|','Location','NorthWest');
title(strcat(ss, ' $\mathcal{C}=\emptyset$'), 'FontSize', fontsize, 'Interpreter', 'latex');

subplot(222)
hh = loglog(xx,Tratio,'o'); hold on
loglog(xx,Tratio,'LineWidth', linewidth); hold off
axis([xmin xmax 0.1 100])
xlabel('SDPNAL+ time')
grid
title(strcat(ss, ' $\mathcal{C}=\emptyset$'), 'FontSize', fontsize, 'Interpreter', 'latex');
legend(hh,'Time ratio: SDPNAL+/BBCPOP','Location','NorthWest');

subplot(223)
idxb = find(abs(L2b)>0);
if (option==1)
   LNAL = L2b(idxb);
   LBBC = L1b(idxb);
   TNAL = T2b(idxb);
   TBBC = T1b(idxb);     
end
xx = TNAL; 
[xx,idx] = sort(xx);
LNAL = LNAL(idx); LBBC = LBBC(idx);
TNAL = TNAL(idx); TBBC = TBBC(idx);
Lratio = (LBBC-LNAL)./abs(LBBC);
Tratio = TNAL./TBBC;
xmin = exp(log(10)*(log10(min(xx))));
xmax = exp(log(10)*(log10(max(xx))));
hh = semilogx(xx,Lratio,'*'); hold on;
semilogx(xx,Lratio,'LineWidth', linewidth);  hold off;
xlabel('SDPNAL+ time')
grid
axis([xmin xmax -0.005 0.005])
legend(hh,'LBv gap: (BBCPOP-SDPNAL+)/|BBCPOP|','Location','NorthWest');
title(strcat(ss, ' $\mathcal{C}$: random'), 'FontSize', fontsize, 'Interpreter', 'latex');

subplot(224)
hh = loglog(xx,Tratio,'o'); hold on
loglog(xx,Tratio,'LineWidth', linewidth); hold off
axis([xmin xmax 0.1 100])
xlabel('SDPNAL+ time')
grid
title(strcat(ss, ' $\mathcal{C}$: random'), 'FontSize', fontsize, 'Interpreter', 'latex');
legend(hh,'Time ratio: SDPNAL+/BBCPOP','Location','NorthWest');

if (table==1)
   print -depsc figure2.eps
else
   print -depsc figure3.eps    
end
