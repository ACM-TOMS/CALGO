/*
 -------------------------------------------------------------------------
 File taputil1.h of ADOL-C version 1.6 as of January 1,   1995          
 Included in ---->
                  adouble.c
		  avector.c
                  fos_reverse.c
                  fov_reverse.c 
                  hos_forward.c
                  hos_reverse.c
                  hov_reverse.c
                  taputil3.c


 ------------------------------------------------------------------------- 

 taputil1.h contains the prototypes for the functions defined in taputil1.c
 These functions provide general information on write statistics and
 provided an interface for writing operations to the 
 tape.
*/


/* Utilities */

extern void get_write_stats(int*, int*, int*, int*);
extern void clear_stats(int);

/* Write routines for scalars */

extern void write_death(locint,locint);
extern void write_int_assign_a(locint,locint);
extern void write_int_assign_d(locint,double);
extern void write_assign_a(locint,locint);
extern void write_assign_d(locint,double);
extern void write_assign_ind(locint); 
extern void write_assign_dep(locint);
extern void write_a_same_arg(unsigned char,locint,locint);
extern void write_d_same_arg(unsigned char,locint,double);
extern void write_args_i_a(unsigned char,locint,int,locint);
extern void write_args_a_i(unsigned char,locint,locint,int);
extern void write_args_d_a(unsigned char,locint,double,locint);
extern void write_args_a_d(unsigned char,locint,locint,double);
extern void write_two_a_rec(unsigned char,locint,locint,locint);
extern void write_single_op(unsigned char,locint,locint);
extern void write_quad(unsigned char,locint,locint,locint);

/* write routines for vectors */

extern void write_intvec_assign_av(locint,locint,locint);
extern void write_assign_vec_dv(locint, locint, double*);
extern void write_assign_av(locint, locint,locint);
extern void write_assign_indvec(locint,locint,double*);
extern void write_assign_depvec(locint size, locint loc1);
extern void write_dv_same_arg(unsigned char,locint,locint, double*);
extern void write_av_same_arg(unsigned char,locint,locint,locint);
extern void write_samearg_av_d(unsigned char,locint,locint,double);
extern void write_two_av_rec(unsigned char,locint,locint,locint,locint);
extern void write_args_dv_av(unsigned char,locint,locint,double*,locint);
extern void write_args_d_av(unsigned char,locint,locint,double,locint);
extern void write_av_a_rec(unsigned char,locint,locint, locint, locint);
extern void write_args_dv_a(unsigned char,locint,locint,double*,locint);


#ifdef conditional
extern void write_condassign(locint result, locint cop, locint r1,locint r2);
extern void write_condassign2(locint result, locint cop, locint r1);
extern void write_associating_value(unsigned char giv_typ, locint location,
				    locint base, locint offset);
extern void write_associating_value_ld(unsigned char giv_typ, double x,
				       locint base, locint offset);
extern void write_associating_vector(unsigned char giv_typ, locint s_loc,
                              locint begin,  locint offset, locint size);
extern void write_associating_vector_ld(double* x, locint begin,  
                                        locint offset, locint size);
#endif

/* End of Include */




