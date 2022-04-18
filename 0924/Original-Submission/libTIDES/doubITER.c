/****************************************************************************
 libTIDES. Version 1.3.0.
 This file is part of TIDES.
 
 Contributors:
 
 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
 Grupo de Mecanica Espacial
 University of Zaragoza
 SPAIN
 
 http://gme.unizar.es/software/tides
 Contact: <tides@unizar.es>
 
 *****************************************************************************/

#include "commonITER.h"
#include "doubITER.h"


void varDB_init(double var[NVARS+1][NDER][MAX_ORDER+1],  double v[], double t)
{
	int i, j, k;
	for (i=0; i<=NVARS; i++) 
		for (j=0; j<NDER; j++)
			for (k=0; k<=MAX_ORDER; k++)
				var[i][j][k]=0.;
	double_set (&var[0][0][0], t);
	double_set_str (&var[0][0][1], "1.");
	for(i=1; i<=NVARS; i++) 
		double_set (&var[i][0][0], v[i-1]);
}
void parDB_init(double par[NPARS][NDER][MAX_ORDER+1], double p[])
{
	int i, j, k;
	for (i=0; i<NPARS; i++) 
		for (j=0; j<NDER; j++)
			for (k=0; k<=MAX_ORDER; k++)
				par[i][j][k]=0.;
	for(i=0; i<NPARS; i++) double_set (&par[i][0][0], p[i]);
	for(i=0; i <NPARTIALS; i++)
		if(!is_variable(PARTIAL_LIST[i]))
			par[PARTIAL_LIST[i]-NVARS-1][i+1][0] = 1;		
	
}
void linkDB_init(double lk[NLINKS][NDER][MAX_ORDER+1])
{
	int i, j, k;
	for (i=0; i<NLINKS; i++) 
		for (j=0; j<NDER; j++)
			for (k=0; k<=MAX_ORDER; k++)
				lk[i][j][k]=0.;
}


void derDB_init(double var[NVARS+1][NDER][MAX_ORDER+1],
		double par[NPARS][NDER][MAX_ORDER+1], double v[])
{
	int i;
/*	int ncol = NVARS+NDER+NFUNS;  */
	extern int clearPartials;
	if (clearPartials)
	for(i=0; i <NPARTIALS; i++) {
		if(is_variable(PARTIAL_LIST[i])) 
			var[PARTIAL_LIST[i]][i+1][0] = 1.;
	}
	else {
		int j, nvf;
		nvf = NVARS+NFUNS;
		for( i = 0; i < NVARS; i++) 
			for(j = 1; j < NDER; j++) 
					var[i+1][j][0] = v[i+(j*nvf)];
	}

	clearPartials = 0;
}


void	write_solution_DB(double cvf[][MAX_ORDER+1],
		double var[NVARS+1][NDER][MAX_ORDER+1],
		double link[NLINKS][NDER][MAX_ORDER+1])
{
	int i,j,k, nvf;
	double vfun;
	nvf = NVARS+NFUNS;
	for( i = 0; i < NVARS; i++) 
		for(j = 0; j < NDER; j++) 
			for(k = 0; k <= MAX_ORDER; k++)
				double_set (&cvf[i+(j*nvf)][k], var[i+1][j][k]);
	for( i = 0; i < NFUNS; i++) 
		for(j = 0; j < NDER; j++) 
			for(k = 0; k <= MAX_ORDER; k++) {
				if (FUNCTION_LIST[i] > 0) vfun = link[FUNCTION_LIST[i]][j][k];
				else vfun = var[-FUNCTION_LIST[i]][j][k];
				double_set (&cvf[NVARS+i+(j*nvf)][k], vfun);
			}
}

/************************************************************************/

void double_htilde(double h[NDER][MAX_ORDER+1], long j, long v, long i, double *ht, int ORDER_INDEX)
{
	double cero;
	double_init (&cero);
	if(ORDER_INDEX >= 0) {
		if( j == ORDER_INDEX && i==v) double_set(ht, cero);
		else double_set(ht, h[v][j]);
	}
}

/************************************************************************/

void	double_var_t(double f[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	double val; double_init (&val);
	int i;
	if(ORDER_INDEX > 0) {
		for(i = 0; i < NDER; i++) {
			double_div_i(&val, f[i][ORDER_INDEX-1], ORDER_INDEX);
			double_set(&w[i][ORDER_INDEX], val);
		}
	}

}

void	double_var_t_c(char* cs, double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	double c;
	if(ORDER_INDEX == 1) {
		double_init(&c);
		double_set_str(&c, cs);
		double_set(&w[0][1], c);
	}
}
void	double_var_t_cc(double c, double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	if(ORDER_INDEX == 1) {
		double_set(&w[0][1], c);
	}
}

void	double_add_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1],
		 double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < NDER; i++) 
			double_add(&w[i][ORDER_INDEX], u[i][ORDER_INDEX], v[i][ORDER_INDEX]);
	}
}

void	double_sub_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1],
		 double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < NDER; i++) 
			double_sub(&w[i][ORDER_INDEX], u[i][ORDER_INDEX], v[i][ORDER_INDEX]);
	}
}



void	double_add_t_c(char* cs, double u[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	double c; 
	double_init (&c);
	if(ORDER_INDEX >= 0) {
		double_set_str(&c, cs);
		if(ORDER_INDEX == 0) 
			double_add(&w[0][ORDER_INDEX], c, u[0][ORDER_INDEX]);
		else 
			double_set(&w[0][ORDER_INDEX], u[0][ORDER_INDEX]);
		for(i = 1; i < NDER; i++) 
			double_set(&w[i][ORDER_INDEX], u[i][ORDER_INDEX]);
	}

}
void	double_add_t_cc(double c, double u[NDER][MAX_ORDER+1], 
						double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) 
			double_add(&w[0][ORDER_INDEX], c, u[0][ORDER_INDEX]);
		else 
			double_set(&w[0][ORDER_INDEX], u[0][ORDER_INDEX]);
		for(i = 1; i < NDER; i++) 
			double_set(&w[i][ORDER_INDEX], u[i][ORDER_INDEX]);
	}
}


void	double_sub_t_c(char* cs, double u[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	double c; double_init (&c);
	if(ORDER_INDEX >= 0) {
		double_set_str(&c, cs);
		if(ORDER_INDEX == 0) 
			double_sub(&w[0][ORDER_INDEX], c, u[0][ORDER_INDEX]);
		else
			double_mul_i(&w[0][ORDER_INDEX], u[0][ORDER_INDEX],-1);
		for(i = 1; i < NDER; i++) 
			double_mul_i(&w[i][ORDER_INDEX], u[i][ORDER_INDEX],-1);
	}

}
void	double_sub_t_cc(double c, double u[NDER][MAX_ORDER+1], 
						double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) 
			double_sub(&w[0][ORDER_INDEX], c, u[0][ORDER_INDEX]);
		else
			double_mul_i(&w[0][ORDER_INDEX], u[0][ORDER_INDEX],-1);
		for(i = 1; i < NDER; i++) 
			double_mul_i(&w[i][ORDER_INDEX], u[i][ORDER_INDEX],-1);
	}	
}



void	double_mul_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1],
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, vi, j;
	double sumint, sumext, partial;
	if(ORDER_INDEX >= 0) {
		sumint = sumext = partial = 0.;
		for(i = 0; i < NDER; i++) {
			sumext = 0.;
			for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
				sumint = 0.;
				for(j = 0; j <= ORDER_INDEX; j++)
					sumint = sumint + 
			u[PREV_VI[vi]][ORDER_INDEX-j]* v[PREV_IV[vi]][j];
				sumext = sumext + sumint*PREV_COEF[vi];
			}
			w[i][ORDER_INDEX] = sumext;	
		}
	}
}


void	double_mul_t_c(char* cs, double u[NDER][MAX_ORDER+1],
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	double c; double_init (&c);
	if(ORDER_INDEX >= 0) {
		double_set_str(&c, cs);
		for(i = 0; i < NDER; i++) 
			double_mul(&w[i][ORDER_INDEX], c, u[i][ORDER_INDEX]);
	}

}
void	double_mul_t_cc(double c, double u[NDER][MAX_ORDER+1],
					   double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < NDER; i++) 
			double_mul(&w[i][ORDER_INDEX], c, u[i][ORDER_INDEX]);
	}
	
}



void double_div_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	double sumint, sumext, partial, zero, ht;
	if(g[0][0] == 0 ) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		g[0][0]= 1;
	}
	if(ORDER_INDEX >= 0) {
		double_init(&sumint);
		double_init(&sumext);
		double_init(&partial);
		double_init(&zero);
		for(i = 0; i < NDER; i++) {
			double_set(&sumext, zero);
			for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
				double_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					double_htilde(h, ORDER_INDEX-j, PREV_VI[vi], i, &ht, ORDER_INDEX);
					double_mul(&partial, ht, g[PREV_IV[vi]][j]);
					double_add(&sumint, sumint, partial);
				}
				double_mul(&sumint, sumint, PREV_COEF[vi]);
				double_add(&sumext, sumext, sumint);
			}
			double_sub(&sumext, f[i][ORDER_INDEX], sumext);
			double_div(&sumext, sumext, g[0][0]);
			double_set(&h[i][ORDER_INDEX], sumext);
		}
	}
}


void	double_div_t_vc(double u[NDER][MAX_ORDER+1], double c, 
						double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < NDER; i++) 
			double_mul(&w[i][ORDER_INDEX], 1./c, u[i][ORDER_INDEX]);
	}
}

void	double_div_t_cv(double c, double u[NDER][MAX_ORDER+1],
						double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	double sumint, sumext, partial, zero, f, wt;
	double_init (&zero);
	if(double_equal (u[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		double_set_d (&u[0][0], 1.);
	}
	if(ORDER_INDEX >= 0) {
		double_init(&sumint);
		double_init(&sumext);
		double_init(&partial);
		double_init(&f);
		double_init(&wt);
		for(i = 0; i < NDER; i++) {
			if(ORDER_INDEX == 0 && i == 0 ) double_set_d(&f, 1.);
			else double_set(&f, zero);
			double_set(&sumext, zero);
			for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
				double_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					double_htilde(w, ORDER_INDEX-j, PREV_VI[vi], i, &wt, ORDER_INDEX);
					double_div(&wt, wt, c);
					double_mul(&partial, wt, u[PREV_IV[vi]][j]);
					double_add(&sumint, sumint, partial);
				}
				double_mul_i(&sumint, sumint, PREV_COEF[vi]);
				double_add(&sumext, sumext, sumint);
			}
			double_sub(&sumext, f, sumext);
			double_div(&sumext, sumext, u[0][0]);
			double_mul(&sumext, sumext, c);
			double_set(&w[i][ORDER_INDEX], sumext);
		}
	}
}

void double_inv_t(double u[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	double sumint, sumext, partial, zero, f, wt;
	double_init (&zero);
	if(double_equal (u[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		double_set_d (&u[0][0], 1.);
	}
	if(ORDER_INDEX >= 0) {
		double_init(&sumint);
		double_init(&sumext);
		double_init(&partial);
		double_init(&f);
		double_init(&wt);
		for(i = 0; i < NDER; i++) {
			if(ORDER_INDEX == 0 && i == 0 ) double_set_d(&f, 1.);
			else double_set(&f, zero);
			double_set(&sumext, zero);
			for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
				double_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					double_htilde(w, ORDER_INDEX-j, PREV_VI[vi], i, &wt, ORDER_INDEX);
					double_mul(&partial, wt, u[PREV_IV[vi]][j]);
					double_add(&sumint, sumint, partial);
				}
				double_mul_i(&sumint, sumint, PREV_COEF[vi]);
				double_add(&sumext, sumext, sumint);
			}
			double_sub(&sumext, f, sumext);
			double_div(&sumext, sumext, u[0][0]);
			double_set(&w[i][ORDER_INDEX], sumext);
		}
	}
}

void double_exp_t(double u[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	double sumint,sumext, partial, zero;
	if(ORDER_INDEX >= 0) {
		zero = 0.;
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_exp(&sumext, u[0][0]);
					double_set(&w[0][0], sumext);
				} else {
					double_set(&sumint, zero);					
					for( vi = PREVSTAR_ACCUM[i]; vi < PREVSTAR_ACCUM[i+1]; vi++) {
						double_mul(&partial, w[PREVSTAR_VI[vi]][0], u[PREVSTAR_IV[vi]][0]);
						double_mul_i(&partial, partial, PREVSTAR_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_set(&w[i][0], sumint);
				}
			}
		} else {		
			for(i = 0; i < NDER; i++) {
				double_set(&sumext, zero);
				for(j = 0; j < ORDER_INDEX; j++) {
					double_set(&sumint, zero);
					for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
						double_mul(&partial, w[PREV_VI[vi]][j], u[PREV_IV[vi]][ORDER_INDEX-j]);
						double_mul_i(&partial, partial, PREV_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_mul_i(&sumint, sumint, ORDER_INDEX-j);
					double_add(&sumext, sumext, sumint);
				}
				double_div_i(&sumext, sumext, ORDER_INDEX);
				double_set(&w[i][ORDER_INDEX], sumext);
			}
			
		}
	}
}


void double_pow_t_c(double u[NDER][MAX_ORDER+1], char* cs, 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	double sumint, sumext, partial, partialb, zero, wt, c;
	double_init (&zero);
	if(double_equal (u[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		double_set_d (&u[0][0], 1.);
	}
	
	if(ORDER_INDEX >= 0) {
		double_init(&sumint);
		double_init(&sumext);
		double_init(&partial);
		double_init(&partialb);
		double_init(&wt);
		double_init(&c);
		double_set_str(&c, cs);
	
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_pow(&sumext, u[0][0], c);
					double_set(&w[0][0], sumext);
				} else {
					double_set(&sumint, zero);					
					for( vi = PREVSTAR_ACCUM[i]; vi < PREVSTAR_ACCUM[i+1]; vi++) {
						double_mul(&partial, w[PREVSTAR_VI[vi]][0], u[PREVSTAR_IV[vi]][0]);
						double_mul(&partial, partial, c);
						double_htilde(w, 0, PREVSTAR_IV[vi], i, &wt, 0);
						double_mul(&partialb, wt, u[PREVSTAR_VI[vi]][0]);
						double_sub(&partial, partial, partialb);
						double_mul_i(&partial, partial, PREVSTAR_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_div(&sumint, sumint, u[0][0]);
					double_set(&w[i][0], sumint);
				}
			}
		} else {		
			for(i = 0; i < NDER; i++) {
				double_set(&sumext, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					double_set(&sumint, zero);
					for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
						double_htilde(w, j, PREV_VI[vi], i, &wt, ORDER_INDEX);
						double_mul(&partial, wt, u[PREV_IV[vi]][ORDER_INDEX-j]);
						double_mul_i(&partial, partial, PREV_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_add_i(&partial, c, 1);
					double_mul_i(&partial, partial, j);
					double_mul_i(&partialb, c, ORDER_INDEX);
					double_sub(&partial, partialb, partial);
					double_mul(&sumint, sumint, partial);
					double_add(&sumext, sumext, sumint);
				}
				double_div(&sumext, sumext, u[0][0]);
				double_div_i(&sumext, sumext, ORDER_INDEX);
				double_set(&w[i][ORDER_INDEX], sumext);
			}
			
		}
		
	}

	
	
}
void double_pow_t_cc(double u[NDER][MAX_ORDER+1], double c, 
					double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	double sumint, sumext, partial, partialb, zero, wt;
	double_init (&zero);
	if(double_equal (u[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		double_set_d (&u[0][0], 1.);
	}
	
	if(ORDER_INDEX >= 0) {
		double_init(&sumint);
		double_init(&sumext);
		double_init(&partial);
		double_init(&partialb);
		double_init(&wt);
		
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_pow(&sumext, u[0][0], c);
					double_set(&w[0][0], sumext);
				} else {
					double_set(&sumint, zero);					
					for( vi = PREVSTAR_ACCUM[i]; vi < PREVSTAR_ACCUM[i+1]; vi++) {
						double_mul(&partial, w[PREVSTAR_VI[vi]][0], u[PREVSTAR_IV[vi]][0]);
						double_mul(&partial, partial, c);
						double_htilde(w, 0, PREVSTAR_IV[vi], i, &wt, 0);
						double_mul(&partialb, wt, u[PREVSTAR_VI[vi]][0]);
						double_sub(&partial, partial, partialb);
						double_mul_i(&partial, partial, PREVSTAR_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_div(&sumint, sumint, u[0][0]);
					double_set(&w[i][0], sumint);
				}
			}
		} else {		
			for(i = 0; i < NDER; i++) {
				double_set(&sumext, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					double_set(&sumint, zero);
					for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
						double_htilde(w, j, PREV_VI[vi], i, &wt, ORDER_INDEX);
						double_mul(&partial, wt, u[PREV_IV[vi]][ORDER_INDEX-j]);
						double_mul_i(&partial, partial, PREV_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_add_i(&partial, c, 1);
					double_mul_i(&partial, partial, j);
					double_mul_i(&partialb, c, ORDER_INDEX);
					double_sub(&partial, partialb, partial);
					double_mul(&sumint, sumint, partial);
					double_add(&sumext, sumext, sumint);
				}
				double_div(&sumext, sumext, u[0][0]);
				double_div_i(&sumext, sumext, ORDER_INDEX);
				double_set(&w[i][ORDER_INDEX], sumext);
			}
			
		}
		
	}
	
	
	
}

/************************************************************************/

void double_sct_0(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], long i)
{
	long vi;
	double sumint, partial, zero;
	double_init(&zero);
	double_init(&partial);
	double_init(&sumint);
	double_set(&sumint, zero);					
	for( vi = PREVSTAR_ACCUM[i]; vi < PREVSTAR_ACCUM[i+1]; vi++) {
		double_mul(&partial, g[PREVSTAR_VI[vi]][0], f[PREVSTAR_IV[vi]][0]);
		double_mul_i(&partial, partial, PREVSTAR_COEF[vi]);		
		double_add(&sumint, sumint, partial);
	}
	double_set(&h[i][0], sumint);
}

void double_sct_i(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], long i, int ORDER_INDEX)
{
	long vi,j;
	double sumint, sumext, partial, zero;
	double_init(&zero);
	double_init(&partial);
	double_init(&sumint);
	double_init(&sumext);
	double_set(&sumext, zero);
	for(j = 1; j <= ORDER_INDEX; j++) {
		double_set(&sumint, zero);
		for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
			double_mul(&partial, g[PREV_VI[vi]][ORDER_INDEX-j], f[PREV_IV[vi]][j]);
			double_mul_i(&partial, partial, PREV_COEF[vi]);		
			double_add(&sumint, sumint, partial);
		}
		double_mul_i(&sumint, sumint, j);		
		double_add(&sumext, sumext, sumint);
	}
	double_div_i(&sumext, sumext, ORDER_INDEX);		
	double_set(&h[i][ORDER_INDEX], sumext);
}

void    double_sin_cos_t (double f[NDER][MAX_ORDER+1], 
		double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		int ORDER_INDEX) {
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_sin (&s[0][0], f[0][0]);
					double_cos (&c[0][0], f[0][0]);
				} else {
					double_sct_0 (f,c,s,i);
					double_sct_0 (f,s,c,i);
					double_mul_i (&c[i][0], c[i][0], -1);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_sct_i(f,s,c,i,ORDER_INDEX);
				double_sct_i(f,c,s,i,ORDER_INDEX);
				double_mul_i (&c[i][ORDER_INDEX],
					 c[i][ORDER_INDEX], -1);
			}
		}
	}
}



void double_sin_t(double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_sin(&w[0][0], s[0][0]);
				} else {
					double_sct_0(s,c,w,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_sct_i(s,c,w,i,ORDER_INDEX);
			}
		}
	}
}

void double_cos_t(double c[NDER][MAX_ORDER+1], double s[NDER][MAX_ORDER+1],
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_cos(&w[0][0], c[0][0]);
				} else {
					double_sct_0(c,s,w,i);
					double_mul_i(&w[i][0], w[i][0], -1);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_sct_i(c,s,w,i,ORDER_INDEX);
				double_mul_i(&w[i][ORDER_INDEX], w[i][ORDER_INDEX], -1);
			}
		}
	}
}



void    double_sinh_cosh_t (double f[NDER][MAX_ORDER+1], 
		double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		int ORDER_INDEX) {
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_sinh (&s[0][0], f[0][0]);
					double_cosh (&c[0][0], f[0][0]);
				} else {
					double_sct_0 (f,c,s,i);
					double_sct_0 (f,s,c,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_sct_i(f,s,c,i,ORDER_INDEX);
				double_sct_i(f,c,s,i,ORDER_INDEX);
			}
		}
	}
}




void double_sinh_t(double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_sinh(&w[0][0], s[0][0]);
				} else {
					double_sct_0(s,c,w,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_sct_i(s,c,w,i,ORDER_INDEX);
			}
		}
	}
}
void double_cosh_t(double c[NDER][MAX_ORDER+1], double s[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_cosh(&w[0][0], c[0][0]);
				} else {
					double_sct_0(c,s,w,i);
					double_mul_i(&w[i][0], w[i][0], -1);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_sct_i(c,s,w,i,ORDER_INDEX);
			}
		}
	}
}

/***************************************************************/

void double_fgt_0(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], long i)
{
	long vi;
	double sumint, partial, zero;
	double_init(&zero);
	double_init(&partial);
	double_init(&sumint);
	double_set(&sumint, zero);					
	if(double_equal (g[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		double_set_d (&g[0][0], 1.);
	}
	
	for( vi = PREVSTAR_ACCUM[i]; vi < PREVSTAR_ACCUM[i+1]; vi++) {
		if(PREVSTAR_VI[vi]>0) {
			double_mul(&partial, h[PREVSTAR_IV[vi]][0], g[PREVSTAR_VI[vi]][0]);
			double_mul_i(&partial, partial, PREVSTAR_COEF[vi]);
			double_add(&sumint, sumint, partial);
		}
	}
	double_sub(&sumint, f[i][0], sumint);
	double_div(&sumint, sumint, g[0][0]);
	double_set(&h[i][0], sumint);
}
void double_fgt_i(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], long i, int ORDER_INDEX)
{
	long vi,j;
	double sumint, sumext, partial, zero, ht;
	double_init(&zero);
	double_init(&partial);
	double_init(&sumint);
	double_init(&sumext);
	double_init(&ht);
	double_set(&sumext, zero);
	if(double_equal (g[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		double_set_d (&g[0][0], 1.);
	}
	for(j = 0; j < ORDER_INDEX; j++) {
		double_set(&sumint, zero);
		for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
			double_htilde(h, ORDER_INDEX-j, PREV_IV[vi], i, &ht, ORDER_INDEX);
			double_mul(&partial, ht, g[PREV_VI[vi]][j]);
			double_mul_i(&partial, partial, PREV_COEF[vi]);
			double_add(&sumint, sumint, partial);
		}
		double_mul_i(&sumint, sumint, ORDER_INDEX-j);
		double_add(&sumext, sumext, sumint);
	}
	double_div_i(&sumext, sumext, ORDER_INDEX);
	double_sub(&sumext, f[i][ORDER_INDEX], sumext);
	double_div(&sumext, sumext, g[0][0]);
	double_set(&h[i][ORDER_INDEX], sumext);
}


void	double_log_t(double u[NDER][MAX_ORDER+1],
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_log(&w[0][0], u[0][0]);
				} else {
					double_fgt_0(u,u,w,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_fgt_i(u,u,w,i,ORDER_INDEX);
			}
		}
	}
}
void	double_asin_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_asin(&h[0][0], f[0][0]);
				} else {
					double_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	double_acos_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_acos(&h[0][0], f[0][0]);
				} else {
					double_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	double_atan_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_atan(&h[0][0], f[0][0]);
				} else {
					double_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	double_asinh_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_asinh(&h[0][0], f[0][0]);
				} else {
					double_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	double_acosh_t(double f[NDER][MAX_ORDER+1],double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_acosh(&h[0][0], f[0][0]);
				} else {
					double_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	double_atanh_t(double f[NDER][MAX_ORDER+1],double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					double_atanh(&h[0][0], f[0][0]);
				} else {
					double_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				double_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}

