create( 'series',
  label = "EF.arctanh.power.01",
  booklabel = "11.6.3",
  general = [ (1/(2*k+1))*z^(2*k+1) ],
  constraints = { abs(z) < 1 },
  function = arctanh,
  lhs = arctanh(z), 
  category = "power series"
):

create( 'contfrac',
  label = "EF.arctanh.sfrac.01",
  booklabel = "11.6.8",
  begin = [[z/(1-z^2), 1]],
  even = [(m*(m-1)/((2*m-3)*(2*m-1)))*z^2/(1-z^2),1],
  odd = [((m-1)*(m-2)/((2*m-3)*(2*m-1)))*z^2/(1-z^2),1],
  constraints = { abs(functions:-argument(1-z^2)) < Pi },
  function = arctanh,
  lhs = arctanh(z),
  category = "S-fraction"
):

create( 'contfrac',
  label = "EF.arctanh.sfrac.02",
  booklabel = "11.6.9",
  begin = [[z, 1]],
  general = [ [(-(m-1)^2*z^2/(4*(m-1)^2 - 1)), 1] ],
  constraints = { abs(functions:-argument(1-z^2)) < Pi },
  function = arctanh,
  lhs = arctanh(z),
  category = "S-fraction"
):

