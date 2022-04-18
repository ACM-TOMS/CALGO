function write_matrix_file(fid,Term,param)
%Internal use function that writes, in the file identified by fid, the 
%expression contained in the matrix Term of polynomial structures. The
%variable 'param' is either an inequality (if a LMI is to be defined) or
%the name of the auxiliar variable

for contdata=1:length(Term{1,1}.data)
    fprintf(fid,'clear Term; \n');
    for conti=1:size(Term,1)
        for contj = conti:size(Term,2)
            fprintf(fid,'Term{%d,%d} = ',conti,contj);
            expr = Term{conti,contj}.opcode{contdata};
            toprint = [];
            contexpr = 1;
            while (contexpr <= length(expr))
                varnow = [];
                while (~strcmp(expr(contexpr),'#'))
                    varnow = strcat(varnow,expr(contexpr));
                    contexpr = contexpr + 1;
                end
                vartype = expr(contexpr+1);
                contexpr = contexpr + 2;
                if ((vartype ~= 'P') && (vartype ~= 'S') && (vartype ~= 'K'))
                    varnow = strcat(varnow,'{');
                    while ((contexpr <= length(expr)) && ~strcmp(expr(contexpr),'+') && ~strcmp(expr(contexpr),'-') && ~strcmp(expr(contexpr),'*') && ~strcmp(expr(contexpr),char(39)) && ~strcmp(expr(contexpr),')') && ~strcmp(expr(contexpr),'('))
                        varnow = strcat(varnow,expr(contexpr));
                        contexpr = contexpr + 1;
                    end
                    varnow = strcat(varnow,'}');
                else
                    while ((contexpr <= length(expr)) && ~strcmp(expr(contexpr),'+') && ~strcmp(expr(contexpr),'-') && ~strcmp(expr(contexpr),'*') && ~strcmp(expr(contexpr),char(39)) && ~strcmp(expr(contexpr),')') && ~strcmp(expr(contexpr),'('))
                        contexpr = contexpr + 1;
                    end
                end
                while ((contexpr <= length(expr)) && (strcmp(expr(contexpr),'+') || strcmp(expr(contexpr),'-') || strcmp(expr(contexpr),'*') || strcmp(expr(contexpr),char(39)) || strcmp(expr(contexpr),')') || strcmp(expr(contexpr),'(')))
                    varnow = strcat(varnow,expr(contexpr));
                    contexpr = contexpr + 1;
                end
                toprint = strcat(toprint,varnow);
            end
            toprint = strcat(toprint,'; \n');
            fprintf(fid,toprint);
        end
    end

    if (strcmp(param,'>') || strcmp(param,'<') || strcmp(param,'<=') || strcmp(param,'>='))
        ineq = param;
        %Construct the LMI
        fprintf(fid,'\nMatrix = [');
        for conti = 1:size(Term,1)
            if (conti > 1)
                fprintf(fid,'\n');
            end
            for contj = 1:size(Term,2)
                if (conti > contj) %Transpose
                    fprintf(fid,' Term{%d,%d}''',contj,conti);
                else
                    fprintf(fid,' Term{%d,%d}',conti,contj);
                end
            end
            fprintf(fid,';');
        end
        fprintf(fid,'];\n\n');

        %Set the LMI
        expr = strcat('LMIs = [LMIs;  Matrix ',ineq);
        expr = strcat(expr,' 0];\n\n');
        fprintf(fid,expr);
    else %An auxiliar variable is to be defined
        fprintf(fid,'\n%s{%s} = [',param,num2str(contdata));
        for conti = 1:size(Term,1)
            if (conti > 1)
                fprintf(fid,'\n');
            end
            for contj = 1:size(Term,2)
                if (conti > contj) %Transpose
                    fprintf(fid,' Term{%d,%d}''',contj,conti);
                else
                    fprintf(fid,' Term{%d,%d}',conti,contj);
                end
            end
            fprintf(fid,';');
        end
        fprintf(fid,'];\n\n');
    end
end

return
