/*
   --------------------------------------------------------------
   File taputil1.c of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
   Provides the definition of the write- interface functions and
   write statistics functions (get_write_stats, clear_stats).
   The majority of these functions are called by the badouble, 
   badoublev operations defined in adouble.c and avector.c.
   These routines write operations to "tape." 
*/

#ifdef __cplusplus
extern "C" {
#endif

#include "dvlparms.h" /* Developers Parameters */

/* Necessary Local Includes */

#include "usrparms.h"
#include "tayutil.h"
#include "oplate.h"
#include "taputil2.h"

/* Static variables for statistics */

static int ind_ptr;
static int dep_ptr;
static int death_ptr;
static int revalso;

/* Store --- external double pointer from adouble.c */
extern double* store;


void get_write_stats(int *x, int *y, int *z, int *w)
{
  *x = dep_ptr;
  *y = ind_ptr;
  *z = death_ptr;
  *w = revalso;
}

void clear_stats(int reverse_flag)
{
  revalso = reverse_flag;
  ind_ptr=0;
  dep_ptr=0;
  death_ptr=0;      
}



/* ---- Write Routines for sequential functions ----- */



void write_death(locint loc1, locint loc2) 
{
  put_op(death_not);
  put_locint(loc1);
  put_locint(loc2);
  death_ptr += loc2-loc1+1;
  if (revalso) 
    {
     write_scaylor(store[loc2]);
     while(loc2>loc1) write_scaylor(store[--loc2]);
     }
}

void write_assign_a( locint loc1,
		     locint loc2)
{
  put_op(assign_a);
  put_locint(loc2);
  put_locint(loc1);
  ++(death_ptr);
  if (revalso) write_scaylor(store[loc1]);
}
void write_int_assign_a( locint loc1,
		         locint loc2)
{
  put_op(int_adb_a);
  put_locint(loc2);
  put_locint(loc1);
  /*  ++(death_ptr); */
}

void write_assign_ind( locint loc1) 
{

  (ind_ptr)++;

  put_op(assign_ind);
  put_locint(loc1);

  ++(death_ptr);
  if(revalso) write_scaylor(store[loc1]);
}

void write_assign_dep( locint loc1)
{
  (dep_ptr)++;

  put_op(assign_dep);
  put_locint(loc1);
}

void write_assign_d( locint loc1, double value)
{
  put_op(assign_d);
  put_locint(loc1);
  put_val(value);
  ++(death_ptr);
  if (revalso)  write_scaylor(store[loc1]);
}

void write_int_assign_d( locint loc1, double value)
{
  put_op(int_adb_d);
  put_locint(loc1);
  put_val(value);
  /* ++(death_ptr); */
  /* if (revalso)  write_scaylor(store[loc1]); */
}
void write_a_same_arg(unsigned char giv_typ,locint result,locint loc1)
{
  put_op(giv_typ);
  put_locint(loc1);
  put_locint(result);
  ++death_ptr;
  if (revalso)  write_scaylor(store[result]);
}
void write_d_same_arg(unsigned char giv_typ,locint result,double value)
{
  put_op(giv_typ);
  put_locint(result);
  put_val(value);
  ++death_ptr;
  if (revalso)  write_scaylor(store[result]);
}
void write_two_a_rec( unsigned char giv_typ,locint result,
		      locint loc_1,locint loc_2)
{
  put_op(giv_typ);
  put_locint(loc_1);
  put_locint(loc_2);
  put_locint(result);
}
void write_args_d_a(unsigned char giv_typ, locint result,
		    double const_val,locint a_loc)
{
  put_op(giv_typ);
  put_locint(a_loc);
  put_locint(result);
  put_val(const_val);
}

void write_single_op(unsigned char giv_typ,
		     locint result,locint old_loc)
{
  put_op(giv_typ);
  put_locint(old_loc);
  put_locint(result);
}
void write_quad(unsigned char giv_typ,
		locint result,locint old_loc,locint deriv_loc)
{
  put_op(giv_typ);
  put_val(store[result]);
  put_locint(old_loc);
  put_locint(deriv_loc);
  put_locint(result);
}




/* -----  Write Routines for vector functions -------- */




void write_assign_av(locint size,locint loc1,locint loc2)
{
  put_op(assign_av);
  put_locint(loc2);
  put_locint(size);
  put_locint(loc1);
  death_ptr+=size;
  if (revalso) write_scaylors((store+loc1),size);
}

void write_intvec_assign_av(locint size,locint result, locint loc1)
{
  put_op(int_av_av);
  put_locint(loc1);
  put_locint(size);
  put_locint(result);
}


void write_assign_indvec(locint size, locint loc1,double *vals)
{
  (ind_ptr)+=size;
  put_op(assign_indvec);
  put_locint(size);
  put_locint(loc1);
  death_ptr+=size;
  if(revalso) write_scaylors((store+loc1),size);
}



void write_assign_depvec(locint size, locint loc1)
{
  (dep_ptr)+=size;
  put_op(assign_depvec);
  put_locint(size);
  put_locint(loc1);
}

void write_assign_vec_dv(locint size, locint loc1, double* vals)
{
  locint space_left, vals_left=size,loc=loc1;
  death_ptr+=size;
  space_left=get_val_space();
  while (space_left<vals_left){
    put_op(assign_dv);
    put_locint(space_left);
    put_locint(loc1);
    put_vals_p(vals,space_left);
    vals+=space_left;
    vals_left-=space_left;
    loc1+=space_left;
    space_left=get_val_space();
  } /* end_while */
  if (vals_left>0){
    put_op(assign_dv);
    put_locint(vals_left);
    put_locint(loc1);
    put_vals_r(vals,vals_left);
  } /* endif */
  if(revalso) write_scaylors((store+loc),size);
}



void write_av_same_arg(unsigned char giv_typ,locint size,
		       locint result,locint loc1)
{
  put_op(giv_typ);
  put_locint(loc1);
  put_locint(size);
  put_locint(result);
  death_ptr+=size;
  if (revalso) write_scaylors((store+result),size);
}
void write_samearg_av_d(unsigned char giv_typ,locint size,
		       locint result,double y)
{
  put_op(giv_typ);
  put_locint(size);
  put_locint(result);
  put_val(y);
  death_ptr+=size;
  if (revalso) write_scaylors((store+result),size);
}



void write_two_av_rec(unsigned char giv_typ, locint size,
			locint result, locint loc_1, locint loc_2)
{
  put_op(giv_typ);
  put_locint(loc_1);
  put_locint(loc_2);
  put_locint(size);
  put_locint(result);
}

void write_av_a_rec(unsigned char giv_typ, locint size,locint result_s,
		    locint x_start,locint a_loc)
{
  put_op(giv_typ);
  put_locint(x_start);
  put_locint(a_loc);
  put_locint(size);
  put_locint(result_s);
}



void write_args_d_av(unsigned char giv_typ,locint size, locint result,
		       double const_val,locint a_loc)
{
  put_op(giv_typ);
  put_locint(a_loc);
  put_locint(size);
  put_locint(result);
  put_val(const_val);
}

#ifdef conditional
void write_condassign(locint result, locint cop, locint r1,locint r2)
{
  put_op(cond_assign);
  put_locint(cop);
  put_locint(r1);
  put_locint(r2);
  put_locint(result);
  ++(death_ptr);
  if (revalso) write_scaylor(store[result]);
}

void write_condassign2(locint result, locint cop, locint r1)
{
  put_op(cond_assign_s);
  put_locint(cop);
  put_locint(r1);
  put_locint(result);
  ++(death_ptr);
  if (revalso) write_scaylor(store[result]);
}

void write_associating_value(unsigned char giv_typ, locint location,
			     locint base, locint offset)
{
  put_op(giv_typ);
  put_locint(base);
  put_locint(offset);
  put_locint(location);
  ++(death_ptr);
  if (revalso) write_scaylor(store[location]);
}
 
void write_associating_value_ld(unsigned char giv_typ, double x,
				locint base, locint offset)
{
  put_op(giv_typ);
  put_val(x);
  put_locint(base);
  put_locint(offset);
  if (giv_typ==subscript_l)
  {
    ++(death_ptr);
    if (revalso) write_scaylor(store[base+(int)store[offset]]);
  }
}

void write_associating_vector(unsigned char giv_typ, locint s_loc,
                              locint begin,  locint offset, locint size)
{
  put_op(giv_typ);
  put_locint(begin);
  put_locint(offset);
  put_locint(size);
  put_locint(s_loc);
  death_ptr+=size;
  if (revalso) write_scaylors((store+s_loc),size); 
}

void write_associating_vector_ld(double* x, locint begin,  locint offset, 
                                 locint size)
{
 locint space_left, vals_left=size,loc=0;
  space_left=get_val_space();
  while (space_left<vals_left){
    put_op(m_subscript_ld);
    put_locint(begin);
    put_locint(offset);
    put_locint(loc);
    put_locint(space_left);
    put_vals_p(x,space_left);
    x+=space_left;
    vals_left-=space_left;
    loc+=space_left;
    space_left=get_val_space();
  } /* end_while */
  if (vals_left>0){
    put_op(m_subscript_ld);
    put_locint(begin);
    put_locint(offset);
    put_locint(loc);
    put_locint(vals_left);
    put_vals_r(x,vals_left);
  } /* endif */
}



#endif


#ifdef __cplusplus
}
#endif

