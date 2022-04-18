% a utility for opening fv faces and vertices structures along a periodic
% boundary.  
%
% [fv] = trim_periodic(fv_in, split, nearness_threshold, indices_to_split)
%
% input: fv_in -- a structure with faces and vertices fields.
%        split -- the boundary along which to split.  e.g., to split at
%                    \pm \pi for a trig function, you would input pi
%        nearness_threshold -- don't attempt to split when points are far
%                               away from the split point.  this is for
%                               speed.
%        indices_to_split -- which coordinates in the vertices field to pay
%                             attention to.  e.g. for just x coord, use
%                              [1].
%
% this function assumes the split point is symmetric about 0.
% i have not paid attention to the direction of the produced normal vectors
% -- this may ruin alignment of an aligned surface.

% daniel brake, 2014.  north carolina state university, mathematics.

function fv = trim_periodic(fv_in, split, nearness_threshold, indices_to_split)


fv = fv_in;


for ii = 1:length(indices_to_split)
	
	ind = indices_to_split(ii);
	
	for jj = 1:length(fv.faces)
		curr_face = fv.faces(jj,:);

		curr_ind_vertices = fv.vertices(curr_face,ind);
	
			
		absolute_nearness = abs(abs(curr_ind_vertices)-split);
		absolute_indicator = absolute_nearness>nearness_threshold;
		
		if all(absolute_indicator)
			continue;
		end
		
		% else there is at least one vertex close to a split.
		
		if all(sign(curr_ind_vertices)==1) %all positive, and near the positive split boundary.
			continue%
		elseif all(sign(curr_ind_vertices)==-1) %all negative, and near the negative split boundary
			continue%
		end
		
		if sum(sign(curr_ind_vertices)==1)==2
			%two positive, one negative
			flip = 1;
			[min_coord,base_index] = min(curr_ind_vertices);
		elseif sum(sign(curr_ind_vertices)==-1)==2
			%two negative, one positive.
			flip = -1;
			[max_coord,base_index] = max(curr_ind_vertices);
		else
			'weird'
			jj
			pause
		end
			
		if base_index==1
			base_index = curr_face(1);
			side_index_1 = curr_face(2);
			side_index_2 = curr_face(3);
		elseif base_index==2
			base_index = curr_face(2);
			side_index_1 = curr_face(1);
			side_index_2 = curr_face(3);
		else
			base_index = curr_face(3);
			side_index_1 = curr_face(1);
			side_index_2 = curr_face(2);
		end
		
		
		
		
		B = fv.vertices(base_index,:);
		A = fv.vertices(side_index_1,:);
		C = fv.vertices(side_index_2,:);



		triangle_size = norm(cross( (A-B) ,(C-B) ))/2;

		if flip == -1

			new_coord = B(ind)-2*split;

			T = B;
			T(ind) = new_coord;


			split_triangle_size = norm(cross( (A-T) ,(C-T) ))/2;

			if triangle_size>split_triangle_size % need to flip to negative side
				%norm(A-T) < norm(A-B) 
				side_1 = T-A;
				side_2 = T-C;

				scale_1 =  (-split-A(ind)) / (new_coord-A(ind));
				scale_2 =  (-split-C(ind)) / (new_coord-C(ind));


				E = A + scale_1* side_1;
				E(ind) = -split;
				F = C + scale_2* side_2;
				F(ind) = -split;

				G = E; 
				G(ind) = split;
				H = F; 
				H(ind) = split;


				new_index = size(fv.vertices,1)+1;

				fv.vertices(new_index+0,:) = E;
				fv.vertices(new_index+1,:) = F;

				fv.vertices(new_index+2,:) = G; 
				fv.vertices(new_index+3,:) = H;


				fv.faces(jj,:) = [base_index new_index+2 new_index+3]; %? 

				fv.faces(end+1,:) = [new_index+0 new_index+1 side_index_1];  %AEF
				fv.faces(end+1,:) = [new_index+1 side_index_1 side_index_2 ];  % ACF



			end
		else % flip is positive, so try flipping the most negative coordinate to the other side.
			
			new_coord = B(ind)+2*split;

			T = B;
			T(ind) = new_coord;%[B(1:ind-1) new_coord B(ind+1:end)];

			split_triangle_size = norm(cross( (A-T) ,(C-T) ))/2;

			if triangle_size>split_triangle_size
				%norm(A-T) < norm(A-B) %
				side_1 = T-A;
				side_2 = T-C;

				scale_1 =  (split-A(ind)) / (new_coord-A(ind));
				scale_2 =  (split-C(ind)) / (new_coord-C(ind));


				E = A + scale_1* side_1;
				E(ind) = split;
				F = C + scale_2* side_2;
				F(ind) = split;

				G = E; 
				G(ind) = -split;
				H = F; 
				H(ind) = -split;



				new_index = size(fv.vertices,1)+1;

				fv.vertices(new_index+0,:) = E;
				fv.vertices(new_index+1,:) = F;

				fv.vertices(new_index+2,:) = G; 
				fv.vertices(new_index+3,:) = H;


				fv.faces(jj,:) = [base_index new_index+2 new_index+3]; %? 

				fv.faces(end+1,:) = [new_index+0 new_index+1 side_index_1];  %AEF
				fv.faces(end+1,:) = [new_index+1 side_index_1 side_index_2 ];  % ACF

			end
		end



	end %re jj
end



end