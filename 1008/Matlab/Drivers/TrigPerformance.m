
for w=1:7
    MC1 = multicomplex(inputconverter(5,w,10^-10))
    tic
    for i=1:50
        sin(MC1);
        
        
    end
    tocsin(w)=toc;
end


for w=1:7
    MC1 = multicomplex(inputconverter(5,w,10^-10))
    tic
    for i=1:50
        cos(MC1);
    end
    toccos(w)=toc;
    
end

for w=1:7
    MC1 = multicomplex(inputconverter(5,w,10^-10))
    tic
    for i=1:50
        tan(MC1);
    end
    toctan(w)=toc;
end


x=[1:7];


semilogy(x,tocsin,'k--')
hold on
semilogy(x,toccos,'k-^')
hold on
semilogy(x,toctan,'k-v')



% % title('Third Derivative error for different h, $$g=\sqrt{e^{2x}+x+e^{x}}$$','FontSize',18,'Interpreter','Latex')
xlabel('Multicomplex Number Order','FontSize',15,'Interpreter','Latex')
ylabel('Computation Time (s)','FontSize',15,'Interpreter','Latex')
legend('sine','cosine','tangent','Location','northwest')
% % axis([1 6 10^-3 10^4])
% % xticks([1 2 3 4 5 6])
% % xticklabels({'1','2','3','4','5','6'})
grid on
grid minor
set(gcf,'PaperPosition',[0 0 14 10.5])
print -dpng TrigPerformance











