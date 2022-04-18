function display(poly)
%DISPLAY (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 10
% 
% Update: Cristiano M. Agulhari
% 2016, Feb, 13

cont = 1;
existsdpvar = 0;
while ((cont <= length(poly.data)) && (existsdpvar == 0))
    if (isa(poly.data(cont).value,'sdpvar'))
        existsdpvar = 1;
    else
        cont = cont + 1;
    end
end

if (existsdpvar)%if (isa(poly.data(1).value,'sdpvar'))
    %It is a sdpvar
    poly.data(cont).value
    disp(['Label: ',poly.label]);
    disp(['Vertices: [', num2str(poly.vertices), ']']);
    deg = zeros(1,length(poly.data(1).exponent));
    for cont=1:length(poly.data(1).exponent)
        deg(cont) =  sum(poly.data(1).exponent{cont});
    end
    disp(['Degrees:  [', num2str(deg), ']']);
    disp(' ');
else
    disp(['Label: ',poly.label]);
    disp(['Vertices: [', num2str(poly.vertices), ']']);
    deg = [];
    for cont=1:length(poly.data(1).exponent)
        deg = [deg sum(poly.data(1).exponent{cont})];
    end
    disp(['Degrees:  [', num2str(deg), ']']);
    disp(' ');
    
    for contsimplex = 1:length(poly.data(1).exponent)
        simplexnames{contsimplex} = char('a'+contsimplex-1);
    end
    %simplexnames = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
    
    strmain = [];
    strunder = [];
    rowsunder = size(poly.data(1).value,1) - 1;
    for cont=1:length(poly.data)
        for contsimplex=1:length(poly.data(cont).exponent)
            expon = poly.data(cont).exponent{contsimplex};
            for contexp=1:length(expon)
                if (expon(contexp) == 1)
                    strmain = [strmain simplexnames{contsimplex} num2str(contexp) '*'];
                elseif (expon(contexp) > 1)
                    strmain = [strmain simplexnames{contsimplex} num2str(contexp) '^' num2str(expon(contexp)) '*'];
                end
            end
        end
        tammain = length(strmain);
        strmain = [strmain '[' num2str(poly.data(cont).value(1,:)) '] '];
            
        strunder = [strunder char(32*ones(rowsunder,tammain - size(strunder,2)))]; %Fill with spaces the string under the main string
        auxunder = [];
        for contaux = 1:rowsunder
            thisrow = ['[' num2str(poly.data(cont).value(contaux+1,:)) '] '];
            if (contaux > 1)
                %Compare the length to the prior rows
                difer = size(thisrow,2) - size(auxunder,2);
                if (difer > 0) %The last row is larger
                    auxunder = [auxunder char(32*ones(contaux-1,difer))];
                elseif (difer < 0) %The last row is shorter
                    thisrow = [thisrow char(32*ones(1,-difer))];
                end     
            end
            auxunder(contaux,:) = thisrow;
        end
        strunder = [strunder auxunder];
        
        difer = size(strunder,2) - size(strmain,2);
        if (difer > 0) %The string under is larger
            strmain = [strmain char(32*ones(1,difer))];
        elseif (difer < 0) %The main string is larger
            strunder = [strunder char(32*ones(rowsunder,-difer))];
        end
        
        if (cont < length(poly.data))
            strmain = [strmain ' + '];
            strunder = [strunder char(32*ones(rowsunder,3))];
        end
    end
    disp([strmain; strunder]);
end
                