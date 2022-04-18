MomentsACPHFit fitter = new MomentsACPHFit(2, 6, 25); 

ContPhaseVar v1 = fitter.fit();

System.out.println("v1:\n"+v1.toString());

System.out.println("m1:\t"+v1.moment(1));

System.out.println("m2:\t"+v1.moment(2));

System.out.println("m3:\t"+v1.moment(3));
		
		