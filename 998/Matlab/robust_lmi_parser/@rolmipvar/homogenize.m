function h = homogenize(varargin)
%homogenization
%
% Author: Alexandre Felipe
% 2014, Dec, 8
%
% Update: Cristiano Agulhari
% 2016, Feb, 2
%
%  Given a list of elements, returns the
%  same list represented with the same degrees.

vertices = 0;
homogeneous = 1;
if(nargin == 1)
    in = varargin{1};
else
    in = varargin;
end
%  Check the number of vertices

for i = 1:nargin
    if(isa(in{i}, 'rolmipvar'))
        vertices = check_vertices(in{i},vertices);
    end
end
if vertices == 0
    for i = 1:nargin
        if(isa(in{i}, 'rolmipvar'))
            h{i} = in{i};
        else
            h{i} = rolmipvar(in{i}, '<>');
        end
    end
    return
end



% Create exponents for non-rolmipvar, and polynomials of ones for
% the homogenization
for cont = 1:length(vertices)
    expzero{cont} = zeros(1,vertices(cont));
    
    if (vertices(cont) > 0)
        %one = sum_{i=1}^{vertices} alpha_i = 1
        clear M
        degaux = [1 zeros(1,vertices(cont)-1)];
        for cont1 = 1:vertices(cont)
            for cont2 = 1:length(vertices)
                if (cont2 == cont)
                    M{cont1}{cont2} = degaux;
                    degaux = circshift(degaux,[0 1]);
                else
                    M{cont1}{cont2} = [0];
                end
            end
            M{cont1}{cont2+1} = 1;
        end
        
        deg = zeros(1,length(vertices));
        deg(cont) = 1;
        one{cont} = rolmipvar(M,'1',vertices,deg);
    end
    
end
%  Find the maximum degree and
% Homogenize the simplexes
for i = 1:nargin
    if(~isa(varargin{i}, 'rolmipvar')) %It is a constant
        h{i} = rolmipvar(in{i}, '<>', 0, 0);
        h{i} = insert_simplex(h{i},zeros(1,length(vertices)-1));%%%%%
        %h{i} = rolmipvar(in{i}, '<>', vertices, expzero);
        deg_coefs(i,:) = zeros(1,length(vertices));
    else
        h{i} = in{i};
        if (length(h{i}.vertices) < length(vertices))
            h{i} = insert_simplex(h{i},vertices(length(h{i}.vertices)+1:end));
        end
        for j = 1:length(vertices) %Iterates the number of simplexes
            deg_coefs(i,j) = sum(h{i}.data(1).exponent{j});
        end
    end
end
maxexp = max(deg_coefs);
% multiply the terms of lower degree by 1 = sum(alpha)
% in order to get all terms with the same degree.
for i = 1:nargin
    for j = 1:length(vertices) %Iterates the number of simplexes
        while(sum(h{i}.data(1).exponent{j}) < maxexp(j))
            h{i} = h{i} * one{j};
        end
    end
    % disp(['homogenize.m: ',  h{i}.label, ' has ', num2str(h{i}.vertices), ' vertices'])
    

end

  
return

function vertices = check_vertices(atual,vertices)
%Function to verify and update the current number of vertices

if (sum(vertices) == 0)
    vertices = atual.vertices;
else
    if (length(vertices) > length(atual.vertices))
        ind = find(vertices(1:length(atual.vertices)) == 0);
        vertices(ind) = atual.vertices(ind);
        
        ind = find(vertices(1:length(atual.vertices)).*(atual.vertices) ~= 0);
        if (sum(abs(vertices(ind) - atual.vertices(ind))) ~= 0)
            %Error: different number of vertices!!!!!
            error('The number of vertices of the defined simplexes must always be the same!');
        end
    else
        ind = find(vertices == 0);
        vertices(ind) = atual.vertices(ind);
        if (length(vertices) < length(atual.vertices))
            vertices = [vertices atual.vertices(length(vertices)+1:end)];
        end
        
        ind = find(vertices.*atual.vertices ~= 0);
        if (sum(abs(vertices(ind) - atual.vertices(ind))) ~= 0)
            %Error: different number of vertices!!!!!
            error('The number of vertices of the defined simplexes must always be the same!');
        end
    end
end