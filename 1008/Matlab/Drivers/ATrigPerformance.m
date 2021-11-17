
for w=1:7
    MC1 = multicomplex(inputconverter(5,w,1))
    tic
    for i=1:50
        asin(MC1);
        
        
    end
    tocasin(w)=toc;
end


for w=1:7
    MC1 = multicomplex(inputconverter(5,w,1))
    tic
    for i=1:50
        acos(MC1);
    end
    tocacos(w)=toc;
    
end

for w=1:7
    MC1 = multicomplex(inputconverter(5,w,1))
    MC2 = multicomplex(inputconverter(3,w,1))
    tic
    for i=1:50
        atan2(MC1,MC2);
    end
    tocatan2(w)=toc;
end

for w=1:7
    MC1 = multicomplex(inputconverter(5,w,10^-10))
    MC2 = multicomplex(inputconverter(3,w,10^-10))
    tic
    for i=1:50
        atan2(MC1,MC2);
    end
    tocatan2Sm(w)=toc;
end


x=[1:7];


semilogy(x,tocasin,'k--')
hold on
semilogy(x,tocacos,'k-^')
hold on
semilogy(x,tocatan2,'k-v')
hold on
semilogy(x,tocatan2Sm,'k-')



% % title('Third Derivative error for different h, $$g=\sqrt{e^{2x}+x+e^{x}}$$','FontSize',18,'Interpreter','Latex')
xlabel('Multicomplex Number Order','FontSize',15,'Interpreter','Latex')
ylabel('Computation Time (s)','FontSize',15,'Interpreter','Latex')
legend('sine','cosine','tangent','tangent small angle','Location','northwest')
% % axis([1 6 10^-3 10^4])
% % xticks([1 2 3 4 5 6])
% % xticklabels({'1','2','3','4','5','6'})
grid on
grid minor
set(gcf,'PaperPosition',[0 0 14 10.5])
print -dpng ATrigPerformance











