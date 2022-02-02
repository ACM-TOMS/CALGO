function [h,fv] = pipe_surface(x,y,z,varargin)


h = [];

radius = 0.1;
n = 7;
	
if ~all([length(x)==length(y) length(x)==length(z) length(y)==length(z)])
	error('inputs must all be the same length')
end

render = true;
closed_left = false;
closed_right = false;


while ~isempty(varargin)
	curr_opt = varargin{1};
	varargin = varargin(2:end);
	switch curr_opt
		case 'r'
			radius = varargin{1};
			varargin = varargin(2:end);
		case 'n'
			n = varargin{1};
			varargin = varargin(2:end);
		case 'closed'
			curr_val = varargin{1};
			varargin = varargin(2:end);
			switch curr_val
				case 'left'
					closed_left = true;
				case 'right'
					closed_right = true;
				case 'both'
					closed_left = true;
					closed_right = true;
				case 'none'
					closed_left = false;
					closed_right = false;
				otherwise
					error('bad option %s to option ''closed'' for pipe_surface', curr_val);
			end
		case 'render'
			render = varargin{1};
			varargin = varargin(2:end);
		otherwise
			error('bad option ''%s'' to pipe_surface', curr_opt);
	end
end



fv.faces = [];
fv.vertices = [];

%ring the left end of the edge
left_ring_offset = 2;
ring = left_ring([x(1) y(1) z(1)],[x(left_ring_offset) y(left_ring_offset) z(left_ring_offset)],radius,n);
while any(isnan(ring))
	left_ring_offset = left_ring_offset+1;
	ring = left_ring([x(1) y(1) z(1)],[x(left_ring_offset) y(left_ring_offset) z(left_ring_offset)],radius,n);
end
fv.vertices = [fv.vertices;ring'];



% ring up the middle of the edge
num_rings = 1;
for ii = 1:length(x)-2
	A = [x(ii) y(ii) z(ii)];
	B = [x(ii+1) y(ii+1) z(ii+1)];
	C = [x(ii+2) y(ii+2) z(ii+2)];
	
	ring = center_ring(A,B,C,radius,n);
	if any(isnan(ring))
		continue
	end
	
	fv.vertices = [fv.vertices;ring'];
	
	ell = (num_rings-1)*(n+1)+1:( (num_rings)*(n+1));
	
	a = ell;
	b = circshift(a,[0 1]);
	c = a+n+1;
	
	t1 = [b;...
		a;...
		c];
	
	t2 = [b;...
		c;...
		circshift(c,[0 1])];

	
	fv.faces = [fv.faces;t1';t2'];
	num_rings = num_rings+1;
end

num_rings = num_rings+1;

% put on the ring for the right end.
right_offset = 1;
ring = right_ring([x(end-right_offset) y(end-right_offset) z(end-right_offset)],[x(end) y(end) z(end)],radius,n);
while any(isnan(ring))
	right_offset = right_offset-1;
	ring =  right_ring([x(end-right_offset) y(end-right_offset) z(end-right_offset)],[x(end) y(end) z(end)],radius,n);
end
fv.vertices = [fv.vertices;ring'];



ell = (num_rings-2)*(n+1)+1:( (num_rings-1)*(n+1));
	
	a = ell;
	b = circshift(a,[0 1]);
	c = a+n+1;
	
	t1 = [b;...
		a;...
		c];
	
	t2 = [b;...
		c;...
		circshift(c,[0 1])];
fv.faces = [fv.faces;t1';t2'];


if closed_left
	fv.vertices = [fv.vertices;[x(1) y(1) z(1)]];
	t1 = [2:n+1; 1:n; ones(1,n)*length(fv.vertices(:,1))];
	fv.faces = [fv.faces; t1'];
end


if closed_right
	fv.vertices = [fv.vertices;[x(end) y(end) z(end)]];
	t1 = [(1:n)+(num_rings-1)*(n+1); (2:n+1)+(num_rings-1)*(n+1); ones(1,n)*length(fv.vertices(:,1))];
	fv.faces = [fv.faces; t1'];
end


if render
	h = patch(fv);
	set(h,'FaceAlpha',0.2);
	set(h,'FaceColor',[1 0 0]);
	set(h,'EdgeColor',[0.5 0 0]);
end

fv.left_cap = 1:n+1;
fv.right_cap = (1:n+1)+(num_rings-1)*(n+1);

end



function ring = center_ring(A,B,C,radius,n)

ray = (B-A+C-B)/2;

r = cross(ray,[0;0;1]);
r = r/norm(r);
ux = r(1);
uy = r(2);
uz = r(3);
theta = -acos(ray(3)/norm(ray));
C = cos(theta);
S = sin(theta);
t = 1-cos(theta);

R = [t*ux^2+C	t*ux*uy-S*uz	t*ux*uz+S*uy;...
	 t*ux*uy+S*uz	t*uy^2+C	t*uy*uz-S*ux;...
	 t*ux*uz-S*uy	t*uy*uz+S*ux	t*uz^2+C];

ring = R*base_ring(radius,n);

ring = ring + repmat( reshape(B,[],1),[1 n+1]);

end

function ring = left_ring(A,B,radius,n)

ray = B-A;

r = cross(ray,[0;0;1]);
r = r/norm(r);
ux = r(1);
uy = r(2);
uz = r(3);
theta = -acos(ray(3)/norm(ray));
C = cos(theta);
S = sin(theta);
t = 1-cos(theta);

R = [t*ux^2+C	t*ux*uy-S*uz	t*ux*uz+S*uy;...
	 t*ux*uy+S*uz	t*uy^2+C	t*uy*uz-S*ux;...
	 t*ux*uz-S*uy	t*uy*uz+S*ux	t*uz^2+C];

ring = R*base_ring(radius,n);

ring = ring + repmat( reshape(A,[],1),[1 n+1]);

end

function ring = right_ring(A,B,radius,n)

ray = B-A;

r = cross(ray,[0;0;1]);
r = r/norm(r);
ux = r(1);
uy = r(2);
uz = r(3);
theta = -acos(ray(3)/norm(ray));
C = cos(theta);
S = sin(theta);
t = 1-cos(theta);

R = [t*ux^2+C	t*ux*uy-S*uz	t*ux*uz+S*uy;...
	 t*ux*uy+S*uz	t*uy^2+C	t*uy*uz-S*ux;...
	 t*ux*uz-S*uy	t*uy*uz+S*ux	t*uz^2+C];

ring = R*base_ring(radius,n);

ring = ring + repmat( reshape(B,[],1),[1 n+1]);

end

function ring = base_ring(radius,n)
%make a ring of radius centered at origin in x-y plane, consisting on n
%segments.
t = linspace(0,2*pi,n+1);

ring = radius*[cos(t);sin(t);zeros(1,n+1)];

end
