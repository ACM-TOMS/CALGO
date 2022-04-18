function nl_opt = prepare_nl_opt_call(obj, x0, lb, ub, options)
  if exist('OCTAVE_VERSION') % running Octave
     % extract MATLAB style options for fmincon
    max_iter = optimget(options, 'MaxIter');
    tol = optimget(options, 'TolFun');
    nl_opt = @() sqp(x0,obj,[],[],lb,ub,max_iter,tol);
  else % running matlab
    nl_opt = @() fmincon(obj,x0,[],[],[],[],lb,ub,[],options);
  end
end
