function outputDI(DI)
% Print the dual space DI to the screen.
for n=1:length(DI)
    count=0;
    for m=1:size(DI{n},1)
        if m==1
            if DI{n}(m,1)==1 || abs(DI{n}(m,1)-1)<1e-10
                ftmp=['d',num2str(DI{n}(m,2:end),'%d')];
            elseif DI{n}(m,1)==-1 || abs(DI{n}(m,1)+1)<1e-10
                ftmp=['-d',num2str(DI{n}(m,2:end),'%d')];
            elseif isreal(DI{n}(m,1))
                ftmp=[num2str(DI{n}(m,1),6),'*d',num2str(DI{n}(m,2:end))];
            else
                ftmp=['(',num2str(DI{n}(m,1),6),')*d',num2str(DI{n}(m,2:end))];
            end
            ftmp(ftmp==' ')=[];
        else
            if DI{n}(m,1)==1 || abs(DI{n}(m,1)-1)<1e-10
                ftmp1=['+d',num2str(DI{n}(m,2:end),'%d')];
            elseif DI{n}(m,1)==-1 || abs(DI{n}(m,1)+1)<1e-10
                ftmp1=['-d',num2str(DI{n}(m,2:end),'%d')];
            elseif isreal(DI{n}(m,1)) && DI{n}(m,1)<0
                ftmp1=[num2str(DI{n}(m,1),6),'*d',num2str(DI{n}(m,2:end),'%d')];
            elseif isreal(DI{n}(m,1))
                ftmp1=['+',num2str(DI{n}(m,1),6),'*d',num2str(DI{n}(m,2:end),'%d')];
            else
                ftmp1=['+(',num2str(DI{n}(m,1),6),')*d',num2str(DI{n}(m,2:end),'%d')];
            end
            ftmp1(ftmp1==' ')=[];
            if m<10
                ftmp(end+1:end+length(ftmp1))=ftmp1;
            else
                count=count+length(ftmp1);
            end
        end
        if count==0
            ftmp(length(ftmp)+1)=' ';
        end
        if length(ftmp)>100
            disp(['         ',ftmp])
            ftmp='';
        end
    end
    if count>0
        ftmp1=['... (remaining ', num2str(count),...
            ' terms are not displayed)'];
        ftmp(end+1:end+length(ftmp1))=ftmp1;
    end
    disp(['         ',ftmp])
    disp(' ')
    if n>10
        disp(['         ... remaining ', num2str(length(DI)-n+1),' dual basis are not displayed!'])
        return
    end
end
