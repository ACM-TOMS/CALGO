% fv = info2fv(BRinfo)
%
% extract the faces from a bertini_real output.


function [fv] = info2fv(BRinfo)


ind = [1 2 3];

fv.vertices = make_vertices(ind, BRinfo);

if isempty(BRinfo.sampler_data)
	fv.faces = make_faces(BRinfo, ind);
else
	fv.faces = make_surface_faces(BRinfo);
end

degen = any(diff(fv.faces(:,[1:3 1]),[],2)==0,2);
fv.faces(degen,:) = [];

% fv = unifyMeshNormals(fv,'alignTo',1);
end




function faces = make_faces(BRinfo, ind)

%
num_total_faces = 0;
for ii = 1:BRinfo.num_faces
	curr_face = BRinfo.faces(ii);
	num_total_faces = num_total_faces + curr_face.num_left + curr_face.num_right + curr_face.top>=0 + curr_face.bottom>=0;
end
num_total_faces = num_total_faces*2;
faces = zeros(num_total_faces, length(ind));
curr_face_index = 1; %this will be incremented in the loop over the faces




for ii = 1:BRinfo.num_faces
	if BRinfo.faces(ii).midslice_index == -1
		continue
	end
	

	
	pass = 1; left_edge_counter = 1;  right_edge_counter = 1;
	
	
	
	while 1
		
		switch pass
			case 1  %the top edge
				pass = pass+1;
				
				
				if BRinfo.faces(ii).top<0
					continue;
				end
				
				curr_edge = -10;
				if strcmp(BRinfo.faces(ii).system_top,'input_critical_curve')
					curr_edge = BRinfo.crit_curve.edges(BRinfo.faces(ii).top+1,:);
				elseif strcmp(BRinfo.faces(ii).system_top,'input_surf_sphere')
					curr_edge = BRinfo.sphere_curve.edges(BRinfo.faces(ii).top+1,:);
				else
					%do a lookup
					for zz = 1:length(BRinfo.singular_curves)
						if strcmp(BRinfo.singular_names{zz},BRinfo.faces(ii).system_top)
							curr_edge = BRinfo.singular_curves(zz).edges(BRinfo.faces(ii).top+1,:);
						end
					end
					
				end
				
				if curr_edge<0
					continue;
				end
				
				curr_edge = curr_edge([3 2 1]);
				
			case 2  %the bottom edge
				pass = pass+1;
				
				if BRinfo.faces(ii).bottom<0
					continue;
				end
				
				curr_edge = -10;
				if strcmp(BRinfo.faces(ii).system_bottom,'input_critical_curve')
					curr_edge = BRinfo.crit_curve.edges(BRinfo.faces(ii).bottom+1,:);
				elseif strcmp(BRinfo.faces(ii).system_bottom,'input_surf_sphere')
					curr_edge = BRinfo.sphere_curve.edges(BRinfo.faces(ii).bottom+1,:);
				else
					%do a lookup
					for zz = 1:length(BRinfo.singular_curves)
						if strcmp(BRinfo.singular_names{zz},BRinfo.faces(ii).system_bottom)
							curr_edge = BRinfo.singular_curves(zz).edges(BRinfo.faces(ii).bottom+1,:);
						end
					end
					
				end
				
				
				if curr_edge<0
					continue;
				end
				
			case 3
				if left_edge_counter <= BRinfo.faces(ii).num_left
					if BRinfo.faces(ii).left(left_edge_counter)<0 %an error check
						continue;
					end
					
					slice_ind = BRinfo.faces(ii).midslice_index+1; %offset by 1.
					edge_ind = BRinfo.faces(ii).left(left_edge_counter)+1; %offset by 1.
					
					curr_edge = BRinfo.critpoint_slices{slice_ind}.edges(edge_ind,:);
					left_edge_counter = left_edge_counter +1; %increment
					
					
				else
					pass = pass+1;
					continue;
				end
			case 4
				if right_edge_counter <= BRinfo.faces(ii).num_right
					
					if BRinfo.faces(ii).right(right_edge_counter)<0
						continue;
					end
					
					slice_ind = BRinfo.faces(ii).midslice_index+2;
					edge_ind = BRinfo.faces(ii).right(right_edge_counter)+1;
					curr_edge = BRinfo.critpoint_slices{slice_ind}.edges(edge_ind,:);
					right_edge_counter = right_edge_counter +1;
					
					curr_edge = curr_edge([3 2 1]);
					
				else
					pass = pass+1;
					continue;
				end
			case 5
				break;
				
		end
		
		
		faces(curr_face_index,:) = [curr_edge(1) curr_edge(2) BRinfo.faces(ii).midpoint+1];
		faces(curr_face_index+1,:) = [curr_edge(2) curr_edge(3) BRinfo.faces(ii).midpoint+1];
		curr_face_index = curr_face_index+2;
		
		

	end

	
end


end















function [faces] = make_surface_faces(BRinfo)


faces = [];

if ~isempty(BRinfo.sampler_data)
	for ii = 1:length(BRinfo.sampler_data)
		
		
		f = BRinfo.sampler_data{ii}+1;
		f(any(f<=0,2),:) = []; % omit problematic faces.
		faces = [faces;f];
% 		h = patch(fv);
% 		
% 		set(h,'FaceColor',colors(ii,:),'FaceAlpha',0.7,'EdgeColor',0.985*colors(ii,:),'EdgeAlpha',0.5);%,'EdgeColor',0.985*colors(ii,:),'EdgeAlpha',0.5
% 		plot_params.handles.surface_samples(ii) = h;
	end
end


end




















function vertices = make_vertices(ind, BRinfo)



vertices = zeros(BRinfo.num_vertices,length(ind));


for ii=1:BRinfo.num_vertices
	vertices(ii,:) = real(transpose(BRinfo.vertices(ii).point(ind)));
end




end %re: function vertices












