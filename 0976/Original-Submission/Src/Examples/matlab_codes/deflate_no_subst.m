% [] = deflate_no_subst(filename, defl_iteration, minorsize)
%
% computes a deflated system from an input system. 
%
% filename - the name of the file to parse and deflate
% defl_iteration -- a somewhat arbitrary number, used to indicate which
% step in the deflation sequence we are in.  ideally, increment it by 1
% each time you deflate, so that you don't repeat function numbers.
% minorsize -- the size of the minors of which you take the determinant,
% when adding more functions to the system to deflate
%
% dani brake 
% 2017

function deflate_no_subst(filename,defl_iteration, minorsize)

b = bertini_input(filename); % from the Bertini_tropical package

system_variables = b.variable_group;
num_vars = length(system_variables);

num_funcs = size(b.functions,1);

vars = sym(zeros(1,num_vars));
for ii = 1:num_vars
	eval(sprintf('syms %s',system_variables{ii}));  %make the variable be a symbol in matlab memory
	vars(ii) = sym(system_variables{ii}); 
	%vars is used in the jacobian call, so that only those derivatives are computed
end



if ~isempty(b.subfunction)
	subfunc_names = b.subfunction(:,1);
else
	subfunc_names = {};
end
num_subfuncs = length(subfunc_names);


if num_subfuncs > 0
	[var_deps,subfunc_deps] = dependency_graph(b); % determine which subfunctions depend on what, computed by regular expressions
end


absolute_var_deps = zeros(num_subfuncs, num_vars);
for ii= 1:num_subfuncs
	for jj = 1:num_vars
		absolute_var_deps(ii,jj) = depends_on_var(ii,jj, var_deps, subfunc_deps);
	end
end



orig_subfunc_args = cell(num_subfuncs,2);

fprintf('\tsetting up subfunctions in memory\n')

for ii = 1:num_subfuncs
	
	curr_subfunc_args = '';
	curr_subfunc_args_for_regexp = '';
	for jj = 1:num_vars
		if absolute_var_deps(ii,jj)
			curr_subfunc_args = sprintf('%s%s, ',curr_subfunc_args,system_variables{jj});
			curr_subfunc_args_for_regexp = sprintf('%s%s,\\s*',curr_subfunc_args_for_regexp,system_variables{jj});
		end
	end
	curr_subfunc_args = curr_subfunc_args(1:end-2);
	curr_subfunc_args_for_regexp = [curr_subfunc_args_for_regexp(1:end-4) '\s*'];
	
	
	orig_subfunc_args{ii,1} = curr_subfunc_args; %store it.
	orig_subfunc_args{ii,2} = curr_subfunc_args_for_regexp; %store it.
	
	
	currstr = sprintf('syms %s(%s)',b.subfunction{ii,1},curr_subfunc_args);
	eval(currstr);
end



%this block creates the functions, preserving subfunctions as subfunctions,
%without substitution.

fprintf('\tmaking functions in memory\n');
function_names = b.functions(:,1);
f = sym(zeros(length(function_names),1));
for ii = 1:length(function_names)
	curr_func = b.functions{ii,2};
	for jj = 1:num_subfuncs
		oldpattern = sprintf('(\\W|^)%s(\\W|$)',subfunc_names{jj});
		newname = sprintf('%s(%s)',subfunc_names{jj},orig_subfunc_args{jj,1});
		newpattern = sprintf('$1%s$2',newname);
		curr_func = regexprep(curr_func,oldpattern,newpattern); 
	end
	f(ii) = eval(curr_func);
end

J = jacobian(f,vars);







% initialize to empty
new_subfunc_names = {};
is_zero_derivative = zeros(length(subfunc_names),num_vars);

fprintf('\tcomputing partial derivatives of subfunctions\n');
for ii = 1:length(subfunc_names)
	for jj = 1:num_vars
		
		if absolute_var_deps(ii,jj)~=0
			
			curr_subfunc = b.symbol_value(subfunc_names{ii});
			
			for kk = 1:num_subfuncs
				oldpattern = sprintf('(\\W|^)%s(\\W|$)',subfunc_names{kk});
				newname = sprintf('%s(%s)',subfunc_names{kk},orig_subfunc_args{kk,1});
				newpattern = sprintf('$1%s$2',newname);
				curr_subfunc = regexprep(curr_subfunc,oldpattern,newpattern); 
			end
			
			
			current_derivative = diff(eval(curr_subfunc),system_variables{jj});
			if current_derivative==0
				is_zero_derivative(ii,jj) = 1;
			end
			symname = sprintf('DIFF_%s_%s',subfunc_names{ii},system_variables{jj});
			if ~b.is_symbol_declared(symname)
				b.declare_and_define(symname,current_derivative,'subfunction');
				new_subfunc_names{end+1} = symname;
			end
		else
			is_zero_derivative(ii,jj) = 1;
		end
		
		
	end
end

R = nchoosek(1:num_funcs,minorsize);
C = nchoosek(1:num_vars,minorsize);
r = size(R,1);
c = size(C,1);
count = 0;


defl_fns = sym([]);
for j = 1:r
	for k = 1:c
		A = det(J(R(j,:),C(k,:)));
		if A ~= 0
			count = count + 1;
			defl_fns(end+1) = A;
		end
	end
end



fprintf('\tdoing partial derivative substitutions for detjac\n')

for ii = 1:length(defl_fns)
	curr_deflfn = char(defl_fns(ii));
	for jj = 1:length(subfunc_names)
		base = sprintf('diff\\(%s\\(%s\\)',subfunc_names{jj},orig_subfunc_args{jj,2});
		for kk = 1:num_vars

			oldname = sprintf('%s,\\s*%s\\)',base,system_variables{kk});
			oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);
			newname = sprintf('DIFF_%s_%s',subfunc_names{jj},system_variables{kk});
			newpattern = sprintf('$1%s$2',newname);
			curr_deflfn = regexprep(curr_deflfn,oldpattern,newpattern);
		end
	end
	defl_fns(ii) = curr_deflfn;
end



fprintf('\tdoing subfuncion substitutions for defl_fns\n')

%now we substitute away the subfunction f(...) statements, by the names of
%the original subfunctions.


for ii = 1:length(defl_fns)
	curr_deflfn = char(defl_fns(ii));
	for jj = 1:length(subfunc_names)
		oldname = sprintf('%s\\(%s\\)',subfunc_names{jj},orig_subfunc_args{jj,2});
		oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);
		newname = sprintf('%s',subfunc_names{jj});
		newpattern = sprintf('$1%s$2',newname);
		curr_deflfn = regexprep(curr_deflfn,oldpattern,newpattern);
	end
	defl_fns(ii) = curr_deflfn;
end





% finally, we need to substitute away any remaining partial derivative
% statements in the subfunctions.
orig_subfunc_names = subfunc_names;

fprintf('\tdoing substitutions for new subfunctions\n');

for ii = length(new_subfunc_names):-1:1
	ind = b.symbol_index(new_subfunc_names{ii},'subfunction');
	curr_subfunc = char(b.subfunction{ind,2});
	for jj = length(orig_subfunc_names):-1:1
		for kk = 1:num_vars
			oldname = sprintf('diff\\(%s\\(%s\\),\\s*%s)',orig_subfunc_names{jj},orig_subfunc_args{jj,2},system_variables{kk});
			oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);
			
			newname = sprintf('DIFF_%s_%s',orig_subfunc_names{jj},system_variables{kk});
			newpattern = sprintf('$1%s$2',newname);
			curr_subfunc = regexprep(curr_subfunc,oldpattern,newpattern);
		end
		
		oldname = sprintf('%s\\(%s\\)',orig_subfunc_names{jj},orig_subfunc_args{jj,2});
		oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);
		
		newname = orig_subfunc_names{jj};
		newpattern = sprintf('$1%s$2',newname);
		curr_subfunc = regexprep(curr_subfunc,oldpattern,newpattern);
	end
	b.subfunction{ind,2} = curr_subfunc;
end


for ii = 1:length(defl_fns)
	b.declare_and_define(sprintf('defl_%i_%i',defl_iteration,ii),defl_fns(ii),'function');
end

output_filename = sprintf('%s_deflated_%i',filename,defl_iteration);
write_bertini_input_file(b.variable_group, b.functions,'filename',output_filename,'options',b.config,'constants',b.constant,'subfunctions',b.subfunction);


end