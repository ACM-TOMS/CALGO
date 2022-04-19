package jphase;

/**
 * This class provides a number of methods to support other classes in the jPhase module. 
 * These include permutations, Gamma functions, factorial, binomial, among others.   
 * @author German Riaño
 */
public class Utils 
{

private	static double cof[] = {76.18009172947146,-86.50532032941677,24.01409824083091,
		-1.231739572450155, 0.1208650973866179e-2,-0.5395239384953e-5};
private static double a[] = new double[101];

private static double B[][] = new double[33][33];
private static long Pers[][] = new long[33][33];

	/**
     * Computes ln( n!/(n-k)! )
	 * @param n
	 * @param k
	 * @return ln( n!/(n-k) )!
	 */
	public static double lnPermut(int n, int k){
		if (k>n) {
			throw new IllegalArgumentException("wrong parameters for permutation function");
			
		}
		if ((n<=32) && (k<=32)) 
			return (B[n][k]!=0)? B[n][k]:
				(B[n][k] = lnFactorial(n) - lnFactorial(n-k) );
		else return lnFactorial(n) - lnFactorial(n-k); 
	}

	/**
     * Computes n!/(n-k)!
	 * @param n
	 * @param k
	 * @return Computes n!/(n-k)!
	 */
	public static double permut(int n, int k){
		if ((k>n)&&(n<0)) {
			throw new IllegalArgumentException("wrong parameters for permutation function");
		}
		if ((n<=32) && (k<=32)){
			if (Pers[n][k]!=0) return Pers[n][k];
			else{
				long mult =1;
				int m = n;
				for (int i=0;i<k;i++) mult *= (m--);
				return (Pers[n][k] = mult );
			}
		}
		else return (long)Math.floor(0.5 + Math.exp(lnPermut(n,k)));
	}

	/**
     * Computes the log of gamma function.
	 * @param xx value
	 * @return lnGamma(xx)
	 */
	public static double lnGamma(double xx){
		double x,y,tmp,ser;
		int j;
		y=x=xx;
		tmp=x+5.5;
		tmp -= (x+0.5)*Math.log(tmp);
		ser = 1.000000000190015;
		for (j=0;j<=5;j++) ser += cof[j]/++y;
		return (float)(- tmp +Math.log(2.506628274631005*ser/x));
	}

	/**
     * Computes the log of Factorial function
	 * @param n
	 * @return ln (n!)
	 */
	public static double lnFactorial(int n){
		if (n<0) 
			throw new ArithmeticException("Negative factorial");
		if (n<=1) return 0.0;
		if (n<=32)return (a[n]!=0)? a[n]: Math.log(fact(n));
		if (n<=100) return (a[n]!=0)? a[n]:(a[n]=lnGamma(n+1.0));
		else return lnGamma(n+1.0);
	}

	private static double fac[] = new double[33];
	private static int topN = 0;
	/**
     * Factorial function ( n! ).
	 * @param n
	 * @return n!
	 */
	public static double fact(int n){
		fac [0] = 1.0;
		int j;
		if (n<0) 
			throw new ArithmeticException("Negative factorial");
		if (n>32) return Math.exp(lnGamma(n + 1.0));
		while (topN<n){
			j = topN++;
			fac[topN]=fac[j]*topN;
		}
		return fac[n];
	}
	
	/**
	 * Incomplete gamma function.
	 * ;
	 * @param a argument
	 * @param x upper limit
	 * @return integ{0,infty,x} e^(-t) t^a / gamma(a) dt
	 */
	public static double gammaP(double a, double x){
		if (x<0.0 || a <= 0.0) 
			throw new IllegalArgumentException("Incorrect arguments for incomplete gamma function ");
		if (x < a+1.0)  // Use series representation
			return gammaSeries(a,x);
		else 
			return (1.0 - gammaContFrac(a,x));
	}
	
	private static double gammaSeries(double a, double x){
		// incomplete gamma using series
		double ap,sum,del;
		double gln = lnGamma(a);
		final int  ITMAX = 100;
		final double EPS = 3.0e-7;
		if (x < 0)
			throw new IllegalArgumentException("x less than 0 in gammaSeies");
		else if (x==0) return 0.0;
		else
		{
			ap = a;
			del = sum = 1.0 /a;
			for (int n=0; n<=ITMAX; n++){
				++ap;
				del *= x / ap;
				sum += del;
				if (Math.abs(del) < Math.abs(sum)*EPS)
					return sum * Math.exp(-x + a*Math.log(x)-gln );
			}
			// did not converge
			throw 
				new ArithmeticException("ITMAX too small, or a to large in routine gammaSeries");
			
		}
	}
	
	private static double gammaContFrac(double a, double x){
		final int  ITMAX = 500;
		final double EPS = 3.0e-7;
		final double MIND = 1e20 * Double.MIN_VALUE;
		
		int i;
		double an,b,c,d,del,h;
		double gln = lnGamma(a);
		b = x + 1.0 - a;
		c = 1.0 / MIND;
		d = 1.0 / b;
		h = d;
		for (i=1; i<=ITMAX; i++){
			an = -i * (i-a);
			b += 2.0;
			d = an*d + b;
			if (Math.abs(d) < MIND) d = MIND;
			c = b + an / c;
			if (Math.abs(c) < MIND) c = MIND;
			d = 1.0 / d;
			del = d*c;
			h *= del;
			if (Math.abs(del -1) < EPS) break;
		}
		if (i>ITMAX)
			throw new ArithmeticException("ITMAX too small, a too large in gammaContFrac");
		return Math.exp(- x + a*Math.log(x) - gln)* h;
	}
	


	/**
     * Binomial coefficient
	 * @param n
	 * @param k
	 * @return n!/ k! (n-k)!
	 */
	public static double binomial(int n,int k){
		return Math.floor(0.5 + Math.exp(lnFactorial(n)-lnFactorial(k)-lnFactorial(n-k)));
	}

	/**
     * ln Binomial coefficient.
	 * @param n
	 * @param k
	 * @return ln [n!/ k! (n-k)!]
	 */
	public static double lnBinomial(int n,int k){
		return (lnFactorial(n)-lnFactorial(k)-lnFactorial(n-k));
	}
	
	/**
     * Power function obtained by multiplying.
	 * @param x
	 * @param n
	 * @return x^n
	 */
	public static double pow(double x, int n){
		double mult = 1.0;
		if (n<0){
			n = -n;
			x = 1.0 / x;
		}
		for (int i=1;i<=n;i++){
			mult *= x;
		}
		return mult;	
	}
	
	/**
     * Euclidean norm btween given arrays
	 * @param v1
	 * @param v2
	 * @return Euclidean norm
	 */
	public static double distance(double v1[],double v2[]){
		int n = v1.length;
		if (n!=v2.length) return -1.0;
		double maxi = -1.0, del, del2;
		for (int i=0;i<n;i++){
			del2 = ( v1[i]-v2[i] );
			del = (v1[i]>0)? del2/v1[i] : del2;
			maxi = Math.max(maxi,del * del);
		}
		return Math.sqrt(maxi);
	}
	
	/**
     * Creates storage for un upper triangular matrix.
	 * @param n
	 * @return An nxn upper triangular matrix.
	 */
	public static double[][] initUpperTriangular(int n){
		double mati[][] = new double[n][];
		for (int i=0;i<n;i++){
			mati[i] = new double[n-i];
		}
		return mati;
	}
	
	
	
	
/*
	public static int combi(int m,int n){
		int nn=Math.min(n,m-n),
			mm=m-nn+1;
		double res = 1.0;
		for (int i=1;i<=nn;i++){
			res *= mm/i;
			mm++;
		}
		return (int)res;
	}
*/
}

