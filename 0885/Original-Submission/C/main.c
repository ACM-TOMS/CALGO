/* This is the driver program for the lnnorm(x), lnanorm(), lnenorm()
 * functions.

   Written by Jean Marie Linhart
              StataCorp LP
              jlinhart@stata.com

   Last modified January 3, 2008
*/

/* if your compiler does not have erf() or log1p() go modify norminc.h */
#include "norminc.h"

#define SEED 314159265
static int N=1000;

double xlin(int i, double Beg, double Dif);
double xrand(int i, double Beg, double Dif);
extern double ren(int i);
void runtests(int NN, double (*xgen)(int i, double b, double d), char *str, unsigned int noisily, char *file, double Beg, double Dif);
#ifdef CALC_TIME
void dowork(double ix, double (*lnn)(double x), double *xtime, double Beg, double Dif);
#endif
void checktable(void);

int main() 
{
	char str[13], file[11];
	unsigned int noisily = 0; /* set to 1 to see when reldifs get bigger*/
	double ix = 1000000;
	double tmp = 0.0e0;
	double Beg;
	double Dif;
#ifdef CALC_TIME
	double cpu_time_used;
#endif

	srand(SEED);

	Beg = -5.0e0;
	Dif = 10.0e0;
	strcpy(str, "evenly spaced");
	strcpy(file, "even.itex");
	runtests(N, xlin, str, noisily, file, Beg, Dif);
	strcpy(str, "random");	
	strcpy(file, "random.itex");
	printf("\n");
	runtests(N, xrand, str, noisily, file, Beg, Dif);	

	Beg = -37.519e0;
	Dif =  2*37.519e0;
	strcpy(str, "evenly spaced");
	strcpy(file, "even2.itex");
	runtests(N, xlin, str, noisily, file, Beg, Dif);
	strcpy(str, "random");	
	strcpy(file, "random2.itex");
	printf("\n");
	runtests(N, xrand, str, noisily, file, Beg, Dif);	

	checktable();
	Beg = -37.519e0;
	Dif =  2*37.519e0;

#ifdef CALC_TIME /* do timings */
	printf("%7.0f iterations \n", ix);
	srand(SEED);
	dowork(ix, lnnorm, &cpu_time_used, Beg, Dif);
	printf("Time for lnnorm  = %20.10f\n", cpu_time_used);  	

	srand(SEED);
	dowork(ix, lnanorm, &cpu_time_used, Beg, Dif);
	printf("Time for lnanorm = %20.10f\n", cpu_time_used);  	

#ifdef HAS_ERF
	srand(SEED);
	dowork(ix, lnenorm, &cpu_time_used, Beg, Dif);
	printf("Time for lnenorm = %20.10f\n", cpu_time_used);  	
#endif
#endif
 	return 0;
}

double xrand(int i, double Beg, double Dif) /* randomly distributed points */
{
	double ans, r;
	r = (double) rand()/ ((double) RAND_MAX);
	ans = Beg + r*Dif;
	return(ans);
}

double xlin(int i, double Beg, double Dif) /* linearly distributed points */
{
	double ans;
	ans = Beg + (Dif/(double)N)*i;
	return(ans);
}

/* check against tabulated values 
   The values are from Table IV from
W. F. Sheppard. 1939. The Probability Integral, 
British Association for the Advancement of Science, Mathematical Tables,
vol 7, Cambridge, University Press, Cambridge, UK
*/
void checktable(void) 
{
	int i, j, istart, iend;
	double x, y, z;
	double lnx, lax;
	double xf, lx;
	double Beg, Dif;
#ifdef HAS_ERF
	double lex;
#endif
	double aeln, reln;
	double aelne, relne;
	double aelna, relna;
	double max_aeln, max_reln, max_aelne, max_relne, max_aelna, max_relna; 
	double mean_aeln, mean_reln, mean_aelne, mean_relne, mean_aelna;
	double mean_relna; 
	FILE *ipf, *opf1;
	char file[14], ofile[11];
	int NN;

	strcpy(file, "sheppard.txt");
	ipf = fopen(file, "r");
	strcpy(ofile, "shepp.itex");
	opf1 = fopen(ofile, "w");
	

	NN = 51;
	for(j=0; j<2; j++) {
		Beg = 0.0e0+ (double)j*5.0e0;
		Dif =  5.0e0;
		max_aeln = max_reln = 0.0e0;
		mean_aeln = mean_reln = 0.0e0;
		max_aelne = max_relne = max_aelna = max_relna = 0.0e0;
		mean_aelne = mean_relne = mean_aelna = mean_relna = 0.0e0;
		for(i=0; i<NN; i++) {
			fscanf(ipf, "%4lf   %24lf\n", &xf, &lx);
			lnx = -lnnorm(-xf);
			aeln = fabs(lx - lnx);
			reln = RELDIF(lx,lnx);
			max_reln = (max_reln > reln ? max_reln : reln);
			mean_reln = mean_reln + reln;
			max_aeln =  (max_aeln > aeln ? max_aeln : aeln);
			mean_aeln = mean_aeln + aeln;
		
			lax = -lnanorm(-xf);
			aelna = fabs(lx - lax);
			relna = RELDIF(lx,lax);
			max_relna = (max_relna > relna ? max_relna : relna);
			mean_relna = mean_relna + relna;
			max_aelna =  (max_aelna > aelna ? max_aelna : aelna);
			mean_aelna = mean_aelna + aelna;
#ifdef HAS_ERF
			lex = -lnenorm(-xf);
			aelne = fabs(lx - lex);
			relne = RELDIF(lx,lex);
			max_relne = (max_relne > relne ? max_relne : relne);
			mean_relne = mean_relne + relne;
			max_aelne =  (max_aelne > aelne ? max_aelne : aelne);
			mean_aelne = mean_aelne + aelne;
#endif
		}

		mean_reln = mean_reln/((double)NN);		
		mean_aeln = mean_aeln/((double)NN);		

		mean_relna = mean_relna/((double)NN);		
		mean_aelna = mean_aelna/((double)NN);		
#ifdef HAS_ERF
		mean_relne = mean_relne/((double)NN);		
		mean_aelne = mean_aelne/((double)NN);		
#endif
		printf("Comparison with tabulated values from %10.4f to %10.4f \n", Beg, Beg+Dif);
		printf("\n lnnorm(): \n");
		printf("Maximum absolute error: %10.5g \n", max_aeln);
		printf("Maximum relative error: %10.5g \n", max_reln);
		printf("Mean absolute error:    %10.5g \n", mean_aeln);
		printf("Mean relative error:    %10.5g \n\n", mean_reln);

		printf("\n lnanorm(): \n");
		printf("Maximum absolute error: %10.5g \n", max_aelna);
		printf("Maximum relative error: %10.5g \n", max_relna);
		printf("Mean absolute error:    %10.5g \n", mean_aelna);
		printf("Mean relative error:    %10.5g \n\n", mean_relna);

#ifdef HAS_ERF
		printf("\n lnenorm(): \n");
		printf("Maximum absolute error: %10.5g \n", max_aelne);
		printf("Maximum relative error: %10.5g \n", max_relne);
		printf("Mean absolute error:    %10.5g \n", mean_aelne);
		printf("Mean relative error:    %10.5g \n\n", mean_relne);
#endif
#ifdef HAS_ERF
		fprintf(opf1, "Comparison with tabulated values from %10.4f to %10.4f \n\n\\medskip \n", Beg, Beg+Dif);
		fprintf(opf1, "\\begin{tabular}{|r|r|r|r|}\n");
		fprintf(opf1, "\\hline\n");
		fprintf(opf1, "Error & {\\tt lnnorm} & {\\tt lnanorm} & {\\tt lnenorm} \\\\ ");
		fprintf(opf1, "\\hline\n");
		fprintf(opf1, "Maximum absolute & %10.5g & %10.5g & %10.5g  \\\\ \n",
			max_aeln, max_aelna, max_aelne);
		fprintf(opf1, "Maximum relative & %10.5g & %10.5g & %10.5g \\\\ \n",
			max_reln, max_relna, max_relne);
		fprintf(opf1, "Mean absolute & %10.5g & %10.5g & %10.5g  \\\\ \n",
			mean_aeln, mean_aelna, mean_aelne);
		fprintf(opf1, "Mean relative & %10.5g & %10.5g & %10.5g ",
			 mean_reln, mean_relna, mean_relne);
		fprintf(opf1, "\\\\ \\hline\n");
		fprintf(opf1, "\\end{tabular}\n\n");
		fprintf(opf1, "\\medskip\n");
#else
		fprintf(opf1, "Comparison with tabulated values from %10.4f to %10.4f \n\n\\medskip \n", Beg, Beg+Dif);
		fprintf(opf1, "\\begin{tabular}{|r|r|r|}\n");
		fprintf(opf1, "\\hline\n");
		fprintf(opf1, "Error type & {\\tt lnnorm} & {\\tt lnanorm} \\\\ ");
		fprintf(opf1, "\\hline\n");
		fprintf(opf1, "Maximum absolute & %10.5g & %10.5g \\\\ \n", max_aeln, max_aelna);
		fprintf(opf1, "Maximum relative & %10.5g & %10.5g \\\\ \n", max_reln, max_relna);
		fprintf(opf1, "Mean absolute & %10.5g & %10.5g  \\\\ \n", mean_aeln, mean_aelna);
		fprintf(opf1, "Mean relative & %10.5g & %10.5g ", mean_reln, mean_relna);
		fprintf(opf1, "\\\\ \\hline\n");
		fprintf(opf1, "\\end{tabular}\n\n");
		fprintf(opf1, "\\medskip\n");
#endif
	}
	fclose(ipf);
	fclose(opf1);
}	

/* This compares the implementations with each other.
   norm(x) and lnnorm(x) are the "base", anorm(x) and enorm(x) are
   compared to norm(x), and lnanorm(x) and lnenorm(x) are compared to
   lnnorm(x).
*/
void runtests(int NN, double (*xgen)(int i, double b, double d), char *str, unsigned int noisily, char *file, double Beg, double Dif)
{
	double x, y, z;
	double lnx, lax, lex;
	double p, an, qn, en;
	double aene, rene;
	double aena, rena;
	double max_aene, max_rene, max_aena, max_rena; 
	double mean_aene, mean_rene, mean_aena, mean_rena; 
	double aelne, relne;
	double aelna, relna;
	double max_aelne, max_relne, max_aelna, max_relna; 
	double mean_aelne, mean_relne, mean_aelna, mean_relna; 
	int i;
	FILE *opf;

	opf = fopen(file, "w");

	NN = NN+1;
	max_aelne = max_relne = max_aelna = max_relna = 0.0e0;
	mean_aelne = mean_relne = mean_aelna = mean_relna = 0.0e0;
	max_aene = max_rene = max_aena = max_rena = 0.0e0;
	mean_aene = mean_rene = mean_aena = mean_rena = 0.0e0;

	for(i=0; i<NN; i++) {
		x = (*xgen)(i, Beg, Dif);
		lnx = (lnnorm(x));
		z = log(norm(x));
		lax = lnanorm(x);
		an = anorm(x);
		qn = norm(x);

		rena = RELDIF(qn,an);
		aena = fabs(qn - an);

		max_rena = (max_rena > rena ? max_rena : rena);
		mean_rena = mean_rena + rena;
		max_aena =  (max_aena > aena ? max_aena : aena);
		mean_aena = mean_aena + aena;

		relna = RELDIF(lnx,lax);
		aelna = fabs(lnx - lax);

		max_relna = (max_relna > relna ? max_relna : relna);
		mean_relna = mean_relna + relna;
		max_aelna =  (max_aelna > aelna ? max_aelna : aelna);
		mean_aelna = mean_aelna + aelna;

#ifdef HAS_ERF
		en = enorm(x);
		lex = lnenorm(x);
	
		rene = RELDIF(qn,en);
		aene = fabs(qn - en);

		max_rene = (max_rene > rene ? max_rene : rene);
		mean_rene = mean_rene + rene;
		max_aene =  (max_aene > aene ? max_aene : aene);
		mean_aene = mean_aene + aene;

		relne = RELDIF(lnx,lex);
		aelne = fabs(lnx-lex);

		max_relne = (max_relne > relne ? max_relne : relne);
		mean_relne = mean_relne + relne;
		max_aelne =  (max_aelne > aelne ? max_aelne : aelne);
		mean_aelne = mean_aelne + aelne;

if (((relna > 1e-14) || (relne > 1e-14) || (rene > 1e-14) || (rena > 1e-14)) && noisily){	

	printf("x = %10.4f  lnnorm(x) =  %22.16e\n", x,lnx);
        printf("                 lnanorm(x) =  %22.16e\n", lax);
        printf("                 lnenorm(x) =  %22.16e\n\n", lex);
        printf("             log(norm(x)) =  %22.16e\n", z);
        printf("              log(anorm(x)) =  %22.16e\n", log(anorm(x)));
        printf("              log(enorm(x)) =  %22.16e\n\n", log(enorm(x)));
        printf("                  norm(x) = %20.16e\n", qn);
        printf("                   anorm(x) = %20.16e\n", an);
        printf("                   enorm(x) = %20.16e\n\n", en);
} 
#endif		

	}

		mean_relna = mean_relna/((double)NN);		
		mean_aelna = mean_aelna/((double)NN);		
		mean_rena = mean_rena/((double)NN);		
		mean_aena = mean_aena/((double)NN);		

		printf("%d %s points on interval from %10.4f to %10.4f \n", NN, str, Beg, Beg+Dif);
		printf("\n anorm(): \n");
		printf("Maximum absolute difference: %10.5g \n", max_aena);
		printf("Maximum relative difference: %10.5g \n", max_rena);
		printf("Mean absolute difference:    %10.5g \n", mean_aena);
		printf("Mean relative difference:    %10.5g \n\n", mean_rena);

#ifdef HAS_ERF
		mean_rene = mean_rene/((double)NN);		
		mean_aene = mean_aene/((double)NN);		
		mean_relne = mean_relne/((double)NN);		
		mean_aelne = mean_aelne/((double)NN);		
		printf("\n enorm(): \n");
		printf("Maximum absolute difference: %10.5g \n", max_aene);
		printf("Maximum relative difference: %10.5g \n", max_rene);
		printf("Mean absolute difference:    %10.5g \n", mean_aene);
		printf("Mean relative difference:    %10.5g \n\n", mean_rene);
		
		fprintf(opf, "%d %s points on interval from %10.4f to %10.4f.\n\n\\medskip \n", NN, str, Beg, Beg+Dif);
		fprintf(opf, "\\begin{tabular}{|r|r|r|r|r|}\n");
		fprintf(opf, "\\hline\n");
		fprintf(opf, "Difference & {\\tt anorm} & {\\tt enorm} & {\\tt lnanorm} & {\\tt lnenorm} \\\\ ");
		fprintf(opf, "\\hline\n");
		fprintf(opf, "Maximum absolute & %10.5g & %10.5g & %10.5g & %10.5g \\\\ \n",
			max_aena, max_aene, max_aelna, max_aelne);
		fprintf(opf, "Maximum relative & %10.5g & %10.5g & %10.5g & %10.5g \\\\ \n",
			max_rena, max_rene, max_relna, max_relne);
		fprintf(opf, "Mean absolute & %10.5g & %10.5g & %10.5g & %10.5g \\\\ \n",
			mean_aena, mean_aene, mean_aelna, mean_aelne);
		fprintf(opf, "Mean relative & %10.5g & %10.5g & %10.5g & %10.5g ", mean_rena,
			mean_rene, mean_relna, mean_relne);
		fprintf(opf, "\\\\ \\hline\n");
		fprintf(opf, "\\end{tabular}\n\n");
#else
		fprintf(opf, "%d %s points on interval from %10.4f to %10.4f. \n\n\\medskip \n", NN, str, Beg, Beg+Dif);
		fprintf(opf, "\\begin{tabular}{|r|r|r|}\n");
		fprintf(opf, "\\hline\n");
		fprintf(opf, "Difference type & {\\tt anorm} & {\\tt lnanorm} \\\\ ");
		fprintf(opf, "\\hline\n");
		fprintf(opf, "Maximum absolute & %10.5g & %10.5g \\\\ \n", max_aena, max_aelna);
		fprintf(opf, "Maximum relative & %10.5g & %10.5g \\\\ \n", max_rena, max_relna);
		fprintf(opf, "Mean absolute & %10.5g & %10.5g  \\\\ \n", mean_aena, mean_aelna);
		fprintf(opf, "Mean relative & %10.5g & %10.5g ", mean_rena, mean_relna);
		fprintf(opf, "\\\\ \\hline\n");
		fprintf(opf, "\\end{tabular}\n\n");
#endif
		printf("\n lnanorm(): \n");
		printf("Maximum absolute difference: %10.5g \n", max_aelna);
		printf("Maximum relative difference: %10.5g \n", max_relna);
		printf("Mean absolute difference:    %10.5g \n", mean_aelna);
		printf("Mean relative difference:    %10.5g \n\n", mean_relna);

#ifdef HAS_ERF

		printf("\n lnenorm(): \n");
		printf("Maximum absolute difference: %10.5g \n", max_aelne);
		printf("Maximum relative difference: %10.5g \n", max_relne);
		printf("Mean absolute difference:    %10.5g \n", mean_aelne);
		printf("Mean relative difference:    %10.5g \n\n", mean_relne);

#endif 
		fclose(opf);

}

#ifdef CALC_TIME
/* This routine is the one that is timed.  There is overhead
   with the timings from the for loop and call to xrand, this overhead
   is timed and subtracted out.
*/
void dowork(double JX, double (*lnn)(double x), double *xtime, double Beg, double Dif)
{
	double x, lnx;
	double ix;
	clock_t start, end;
	double ltime=0.0e0;
	int i = 3141;

	*xtime = 0.0e0;
	start = clock();
	for(ix = 0.0e0; ix < JX; ix += 1.0e0) {
		x = xrand(i, Beg, Dif);
		lnx = lnn(x);
	}
	end = clock();
	*xtime = ((double) (end-start)) / CLOCKS_PER_SEC;

	ltime = 0.0e0;
	start = clock();
	for(ix = 0.0e0; ix < JX; ix += 1.0e0) {
		x = xrand(i, Beg, Dif);
	}
	end = clock();
	ltime = ((double) (end-start)) / CLOCKS_PER_SEC;
	*xtime = *xtime - ltime;
}
#endif
