function H = par2h(p)
  switch p(1)
    case 1,  H = [p(2) 1   ; -1   p(5)];
    case 0,  H = [1    p(4); p(3) 1   ];
    case -1, H = [p(2) p(4); p(3) p(5)];
    case -2, H = [1    0   ; 0    1   ];
  end
end