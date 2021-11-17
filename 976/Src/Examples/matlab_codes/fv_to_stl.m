
function fv_to_stl(src,hndl,fv)

if nargin==0
	error('must pass an fv to this function');
end

if nargin==1
	fv = src;
end


[options, user_cancelled] = getopt();
			

if user_cancelled
	return;
end



%delete the bad degenerate faces.  
degen = any(diff(fv.faces(:,[1:3 1]),[],2)==0,2);
fv.faces(degen,:) = [];

num_degen = sum(degen);
if num_degen>0
	display(sprintf('removed %d degenerate faces',num_degen));
end


if options.remove_duplicates
			% now we will look for duplicate points
			num_dup = 0;
			for ii = 1:length(fv.vertices)
				for jj = ii+1:length(fv.vertices)
					if norm(fv.vertices(ii,:)-fv.vertices(jj,:))<0.00001
						num_dup = num_dup+1;
						[row_ind,col_ind] = find(fv.faces==jj);
						if ~isempty(row_ind)
							fv.faces(row_ind,col_ind) = jj;
						end
					end
				end
			end


			%remove unreferenced points.

			F = fv.faces;

			keep = [];
			for ii = 1:length(fv.vertices)
				if ~any(any(fv.faces==ii))
					[a,b] = find(fv.faces>ii);
					for jj = 1:length(a)
						F(a(jj),b(jj)) = F(a(jj),b(jj))-1;
					end

				else
					keep = [keep ii];
				end
			end

			fv.vertices = fv.vertices(keep,:);
			fv.faces = F;


			%delete the bad degenerate faces.  
			degen = any(diff(fv.faces(:,[1:3 1]),[],2)==0,2);
			fv.faces(degen,:) = [];

			fv.faces(1,:) = fv.faces(1,[3 2 1]);
end %re: remove duplicate points


save(sprintf('%s.mat',options.filename),'fv');
stlwrite(sprintf('%s.stl',options.filename),fv);

% Fix non-uniform face orientations
%fv = unifyMeshNormals(fv,'alignTo',1);

if options.align
	fv_unified = unifyMeshNormals(fv,'alignTo','in');
	save(sprintf('%s.mat',options.filename),'fv_unified','-append');
	stlwrite(sprintf('%s_unified.stl',options.filename),fv_unified);
end


end



%http://www.mathworks.com/matlabcentral/fileexchange/25862-inputsdlg--enhanced-input-dialog-box--v2-0-5-


function [answer, cancelled] = getopt()

Title = 'Bertini_real STL options';

%%%% SETTING DIALOG options
% options.WindowStyle = 'modal';
options.Resize = 'on';
options.Interpreter = 'none';
options.CancelButton = 'on';
options.ApplyButton = 'off';
options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration


% prompt = {};
% formats = {};
% default_answers = struct([]);






prompt(1,:) = {'Align faces', 'align','this can take a loooong time'};
formats(1,1).type = 'check';
default_answers.align = false;


prompt(2,:) = {'Remove duplicate points', 'remove_duplicates','this also can take a loooong time'};
formats(2,1).type = 'check';
default_answers.remove_duplicates = false;


prompt(3,:) = {'duplicate dist tolerance', 'duplicate_tolerance',[]};
formats(2,2).type = 'edit';
formats(2,2).format = 'float';
formats(2,2).size = [100 0];
formats(2,2).limits = [0 1e10];
formats(2,2).required = 'on';

default_answers.duplicate_tolerance = 1e-4;



prompt(4,:) = {'Filename', 'filename',[]};
formats(3,1).type = 'edit';
formats(3,1).format = 'text';
formats(3,1).size = [200 0];
formats(3,1).limits = [0 1];
default_answers.filename = 'br_surf';



[answer,cancelled] = inputsdlg(prompt,Title,formats,default_answers,options);


end
