/* 
 -------------------------------------------------------------------------
 File oplate.h of ADOL-C version 1.6 as of January 1,   1995          
 Included in ---->
                  adouble.c
		  avector.c 
                  fos_reverse.c
                  fov_reverse.c 
                  hos_forward.c
                  hos_reverse.c
                  hov_reverse.c
                  taputil1.c
                  taputil2.c
                  taputil3.c

		  
 Provides numeric values for the various opcodes used by ADOL-C.

 ------------------------------------------------------------------------- 

*/

/*
  Opcodes for scalar operations
*/

#define death_not 1
#define assign_ind 2
#define assign_dep 3
#define assign_a 4
#define assign_d 5
#define eq_plus_d 6
#define eq_plus_a 7
#define eq_min_d 8
#define eq_min_a 9
#define eq_mult_d 10
#define eq_mult_a 11
#define plus_a_a 12
#define plus_d_a 13
#define min_a_a 14
#define min_d_a 15
#define mult_a_a 16
#define mult_d_a 17
#define div_a_a 18
#define div_d_a 19
#define exp_op 20
#define cos_op 21
#define sin_op 22
#define atan_op 23
#define log_op 24
#define pow_op 25

/* New as of 4/9/90 */

#define asin_op 26
#define acos_op 27
#define sqrt_op 28
/* removed 1/95
#define eq_div_a 29
#define eq_div_d 30
#define tan_op 31
*/

/* New as of 7/3/90 */

#define ignore_me 0
#define gen_quad 32 /* A General Quadrature */

/* New as of 6/10/93 */

#define int_adb_a 33 /* Initialize an adouble with another adouble */
#define int_adb_d 34 /* Initialize an adouble with a double value  */

/* New as of 7/13/93 */

/* Opcodes for tape delimiters. */

#define end_of_tape 35
#define start_of_tape 36

#define end_of_op 37
#define end_of_int 38
#define end_of_val 39

/* For vector operations */

#define plus_av_av    40
/* removed 1/95
#define plus_dv_av    41
*/
#define sub_av_av     42
/* removed 1/95
#define sub_dv_av     43
#define sub_av_dv     44
*/
#define dot_av_av     45
/* removed 1/95
#define dot_dv_av     46
*/
#define mult_a_av     47
#define mult_d_av     48
/* removed 1/95
#define mult_a_dv     49
*/
#define int_av_av     50
/* removed 1/95
#define int_av_dv     51  
*/
#define assign_av     52
#define assign_dv     53
#define assign_indvec 54
#define assign_depvec 55
/* removed 1/95
#define eq_min_dv     56
*/
#define eq_min_av     57
/* removed 1/95
#define eq_plus_dv    58
*/
#define eq_plus_av    59
#define div_av_a      60
#define eq_mult_av_d  61
#define eq_mult_av_a  62
/* removed 1/95
#define dot_av_dv     63
*/
#define mult_av_a     64

#define cond_assign   70
#define cond_assign_s 71
 
#define m_subscript   72
#define m_subscript_l 73
#define m_subscript_ld 74

#define subscript     75
#define subscript_l   76
#define subscript_ld  77

/* removed 1/95
#define cross_av_av   80
#define mult_cv3_av4  81
*/
/*New as of 8/14/94 */
#define abs_val       101



