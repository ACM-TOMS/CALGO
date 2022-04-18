

tic
for w=1:7
     for i=1:50
        (multicomplex(inputconverter(5,w,10^-10))+multicomplex(inputconverter(5,w,10^-10)));
     end
     tocadd(w)=toc;
     tic
end

tic
for w=1:7
    for i=1:50
        (multicomplex(inputconverter(0.5,w,10^-10)))*(multicomplex(inputconverter(0.5,w,10^-10)));
    end
    tocmult(w)=toc;
    tic
end

for w=1:7
    for i=1:50
        (multicomplex(inputconverter(0.5,w,10^-10)))/(multicomplex(inputconverter(0.5,w,10^-10)));
    end
    tocdiv(w)=toc;
    tic
end


x=[1:7];


semilogy(x,tocadd,'k--')
hold on
semilogy(x,tocmult,'k-^')
hold on
semilogy(x,tocdiv,'k-v')



% % title('Third Derivative error for different h, $$g=\sqrt{e^{2x}+x+e^{x}}$$','FontSize',18,'Interpreter','Latex')
xlabel('Multicomplex Number Order','FontSize',15,'Interpreter','Latex')
ylabel('Computation Time (s)','FontSize',15,'Interpreter','Latex')
legend('plus/minus','multiplication','division','Location','northeast')
% % axis([1 6 10^-3 10^4])
% % xticks([1 2 3 4 5 6])
% % xticklabels({'1','2','3','4','5','6'})

grid on
grid minor
set(gcf,'PaperPosition',[0 0 14 10.5])
print -dpng AddMultDivPerformance






