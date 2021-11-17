#ifndef BERTINI_HEADERS_H
#define BERTINI_HEADERS_H


/* \file bertini_headers.hpp */

#include <cstddef>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>

#include <mpi.h> // this *cannot* be inside an extern "C"{} wrapper.

extern "C" {
#include <bertini.h>
#include <cascade.h>
#include <parallel.h>
#include <localdim.h>
}

extern "C" {
	
int change_prec_prog_deriv(void const *ED, int prec);
	
void start_system_eval_data_clear_d(start_system_eval_data_d *SSED);//actually lives in bertini library...  testing if this works.

void patch_eval_data_clear_d(patch_eval_data_d *PED);//another which lives in bertini
void patch_eval_data_clear_mp(patch_eval_data_mp *PED);//another which lives in bertini
void changePatchPrec_mp(int new_prec, patch_eval_data_mp *PED); // in bertini


/**
 from the bertini library.  the prototype is not in any header file.
 */
int checkForReal_d(point_d Pt, double realTol);
/**
 from the bertini library.  the prototype is not in any header file.
 */
int checkForReal_mp(point_mp Pt, double realTol);


/**
 from the bertini library.  the prototype is not in any header file.
 */
void findMultSol(post_process_t *endPoints, int num_sols, int num_vars, preproc_data *PPD, double finalTol);
	
	void bcast_prog_t(prog_t *Prog, int MPType, int my_id, int headnode);



/**
 from the bertini library.  the prototype is not in any header file.
 */
void findFiniteSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double maxNorm);

/**
 from the bertini library.  the prototype is not in any header file.
 */
void findRealSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double realTol);


/**
 from the bertini library.  the prototype is not in any header file.
 */
void findSingSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double maxCondNum, double finalTol, int regenToggle);




/**
 from the bertini library.  the prototype is not in any header file.
 */
void getDehomPoint_comp_d(point_d dehomPoint, int *origErrorIsInf, double *origErrorEst, comp_d *sol, int num_vars, preproc_data *PPD, double accuracyEstimate);


/**
 from the bertini library.  the prototype is not in any header file.
 */
void getDehomPoint_comp_mp(point_mp dehomPoint, int *origErrorIsInf, double *origErrorEst, comp_mp *sol, int num_vars, preproc_data *PPD, double accuracyEstimate);

	
	
}


#endif

