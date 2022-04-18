ContPhaseVar v1 = DenseContPhaseVar.Erlang(0.8, 3);
ContPhaseVar v2 = DenseContPhaseVar.Erlang(1.5, 2);

ContPhaseVar v3 = v1.sum(v2);
System.out.println("v3:\n"+v3.toString());

		
		