function p = triple(k1,k2,k3)
%
%  test polynomial suggested by Goedecker
%  square of Fibocacci polynomial
%

    p = poly([0.9*ones(1,k1),ones(1,k2),1.1*ones(1,k3)]);
