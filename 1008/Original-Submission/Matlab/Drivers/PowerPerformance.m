

for w=1:7
    MC1 = multicomplex(inputconverter(5,w,10^-10))
    tic
    for i=1:50
        MC1^0.5;
    end
    tocsqrt(w)=toc;
end


for w=1:7
    MC1 = multicomplex(inputconverter(5,w,10^-10))
    tic
    for i=1:50
        MC1^MC1;
    end
    tocsqrt2(w)=toc;
    
end


for w=1:7
    MC1 = multicomplex(inputconverter(5,w,10^-10))
    tic
    for i=1:50
        MC1^2;
    end
    tocsqrt3(w)=toc;
    
end


for w=1:7
    MC1 = multicomplex(inputconverter(5,w,10^-10))
    tic
    for i=1:50
        MC1^20;
    end
    tocsqrt4(w)=toc;
    
end


x=[1:7];


semilogy(x,tocsqrt,'k--')
hold on
semilogy(x,tocsqrt2,'k-^')
hold on
semilogy(x,tocsqrt3,'k-v')
hold on
semilogy(x,tocsqrt4,'k-*')



% % title('Third Derivative error for different h, $$g=\sqrt{e^{2x}+x+e^{x}}$$','FontSize',18,'Interpreter','Latex')
xlabel('Multicomplex Number Order','FontSize',15,'Interpreter','Latex')
ylabel('Computation Time (s)','FontSize',15,'Interpreter','Latex')
legend('z_n^{0.5}','z_n^{z_m}','z_n^{2}','z_n^{20}','Location','northwest')
% % axis([1 6 10^-3 10^4])
% % xticks([1 2 3 4 5 6])
% % xticklabels({'1','2','3','4','5','6'})
grid on
grid minor
set(gcf,'PaperPosition',[0 0 14 10.5])
print -dpng PowerPerformance












