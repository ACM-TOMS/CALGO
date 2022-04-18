/*
 -------------------------------------------------------------------------
 File taputil2.h of ADOL-C version 1.6 as of January 1,   1995          
 Included in ---->
                  fos_reverse.c
                  fov_reverse.c 
                  hos_forward.c
                  hos_reverse.c
                  hov_reverse.c
		  taputil1.c
                  taputil3.c


 ------------------------------------------------------------------------- 

 taputil2.h contains the prototypes for the functions (and macros) 
 defined in taputil2.c.  The functions and macros provide
 get and put like operation access to the tape.

*/
extern void set_buf_size(int);
extern void set_buffers(char*,unsigned char*,char*, locint*,char*, double *);
extern void close_tape(int*, int);
extern void init_rev_sweep(int);
extern void init_for_sweep(int);
extern void end_sweep();
extern void put_op(unsigned char);
extern void put_val(double);
extern void put_vals_p(double *,int);
extern void put_vals_r(double *,int);
extern int get_val_space();
extern void put_locint(locint);
extern void get_op_block_f();
extern void get_loc_block_f();
extern void get_val_block_f();
extern void get_op_block_r();
extern void get_loc_block_r();
extern void get_val_block_r();
extern double * get_val_v_f(locint);
extern double * get_val_v_r(locint);
extern void reset_val_r();
#ifdef DEBUG
extern unsigned char get_op_f();
#endif
extern unsigned char *op_codes,*g_op_ptr;
extern locint *loc_tape,*g_loc_ptr;
extern double *real_tape,*g_real_ptr;
extern int op_ptr;
extern int loc_ptr;
extern int real_ptr;

#define get_op_r() *(--g_op_ptr)
#define get_locint_r() *(--g_loc_ptr)
#define get_val_r() *(--g_real_ptr)
#ifndef DEBUG
#define get_op_f() *g_op_ptr++ 
#endif
#define get_locint_f() *g_loc_ptr++
#define get_val_f() *g_real_ptr++

