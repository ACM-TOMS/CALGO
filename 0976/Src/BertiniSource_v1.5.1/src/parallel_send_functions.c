// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"
#include "regeneration.h"
#include "eqbyeq.h"
#include "dimbydim.h"
#include "pos_dim.h"
#include "regen_pos_dim.h"

// This file contains functions that are used to send data structures over MPI

#ifdef _HAVE_MPI

void bcast_worker_info(worker_info *WI, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts worker_info                                 *
\***************************************************************/
{ 
  MPI_Datatype mpi_wi;
  
  // create the datatypes mpi_wi
  create_worker_info(&mpi_wi);
  
  if (my_id == headnode)
  { // send WI
    MPI_Bcast(WI, 1, mpi_wi, headnode, MPI_COMM_WORLD);
  }
  else // worker process
  { // recv WI
    MPI_Bcast(WI, 1, mpi_wi, headnode, MPI_COMM_WORLD);
  }
 
  // free mpi_wi
  MPI_Type_free(&mpi_wi);

  return;
}

void bcast_tracker_config_t(tracker_config_t *T, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts tracker_config_t                            *
\***************************************************************/
{
  MPI_Datatype mpi_tc;
  tracker_config_t_relevant T_rel;

  // create the datatypes mpi_tc
  create_tracker_config_t(&mpi_tc);

  if (my_id == headnode)
  { // setup T_rel & then broadcast it
    cp_tracker_config_t_relevant(&T_rel, T, 0);
    MPI_Bcast(&T_rel, 1, mpi_tc, headnode, MPI_COMM_WORLD);
  }
  else // worker process
  { // recv T_rel from headnode and then setup T
    MPI_Bcast(&T_rel, 1, mpi_tc, headnode, MPI_COMM_WORLD);
    cp_tracker_config_t_relevant(T, &T_rel, 1);
  }
  
  // free mpi_tc
  MPI_Type_free(&mpi_tc);

  return;
}

void bcast_basic_eval_data_d(basic_eval_data_d *BED, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts basic_eval_data_d                           *
\***************************************************************/
{
  MPI_Datatype mpi_bed_int, mpi_comp_d;
  basic_eval_data_d_int BED_int;
  int size = 0, *prog_inst = NULL, *prog_gp_sizes = NULL, *orig_deg = NULL, *new_deg = NULL, *P = NULL, *W = NULL, *startDeg = NULL, *ppd_type = NULL, *ppd_size = NULL;
  char *progStr = NULL;
  comp_d *st_coeff = NULL, *sq_coeff = NULL, *patch_coeff;

  // create the datatype mpi_bed_int, mpi_comp_d
  create_basic_eval_data_d_int(&mpi_bed_int);
  create_comp_d(&mpi_comp_d);

  if (my_id == headnode)
  { // setup BED_int & the other structures and then broadcast it
    cp_basic_eval_data_d_int(&BED_int, BED, &prog_inst, &prog_gp_sizes, &progStr, 1, &orig_deg, &new_deg, &P, &W, &startDeg, &st_coeff, &sq_coeff, &patch_coeff, &ppd_type, &ppd_size, MPType, 0);
    MPI_Bcast(&BED_int, 1, mpi_bed_int, headnode, MPI_COMM_WORLD);

    // broadcast the prog structures
    MPI_Bcast(prog_inst, BED_int.squareSystem_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(prog_gp_sizes, BED_int.squareSystem_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(progStr, BED_int.squareSystem_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // broadcast the square system structures
    MPI_Bcast(orig_deg, BED_int.squareSystem_int.size_f, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(new_deg, BED_int.squareSystem_int.size_r, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(P, BED_int.squareSystem_int.size_f, MPI_INT, headnode, MPI_COMM_WORLD);
    size = BED_int.squareSystem_int.size_r * (BED_int.squareSystem_int.size_f - BED_int.squareSystem_int.size_r);
    MPI_Bcast(W, size, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(sq_coeff, BED_int.squareSystem_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // broadcast the start system structures
    MPI_Bcast(startDeg, BED_int.startSystem_int.size_r, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(st_coeff, BED_int.startSystem_int.num_coeff, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // broadcast the preproc data
    size = BED_int.preProcData_int.num_hom_var_gp + BED_int.preProcData_int.num_var_gp;
    MPI_Bcast(ppd_type, size, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(ppd_size, size, MPI_INT, headnode, MPI_COMM_WORLD);

    // broadcast the patch structures
    size = BED_int.patch.patchCoeff_rows * BED_int.patch.patchCoeff_cols;
    MPI_Bcast(patch_coeff, size, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // free the memory
    free(prog_inst);
    free(prog_gp_sizes);
    free(orig_deg);
    free(new_deg);
    free(P);
    free(W);
    free(startDeg);
    free(ppd_type);
    free(ppd_size);
    free(progStr);
    free(st_coeff);
    free(sq_coeff);
    free(patch_coeff);
  }
  else // worker process
  { // recv BED_int
    MPI_Bcast(&BED_int, 1, mpi_bed_int, headnode, MPI_COMM_WORLD);
    // allocate and receive the other structures needed to finish setting up BED

    // recv the prog structures
    prog_inst = (int *)bmalloc(BED_int.squareSystem_int.Prog_int.size * sizeof(int));
    MPI_Bcast(prog_inst, BED_int.squareSystem_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);    
    prog_gp_sizes = (int *)bmalloc(BED_int.squareSystem_int.Prog_int.num_var_gps * sizeof(int));
    MPI_Bcast(prog_gp_sizes, BED_int.squareSystem_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
    progStr = (char *)bmalloc(BED_int.squareSystem_int.Prog_int.totalLength * sizeof(char));
    MPI_Bcast(progStr, BED_int.squareSystem_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // recv the square system structures
    orig_deg = (int *)bmalloc(BED_int.squareSystem_int.size_f * sizeof(int));
    MPI_Bcast(orig_deg, BED_int.squareSystem_int.size_f, MPI_INT, headnode, MPI_COMM_WORLD);
    new_deg = (int *)bmalloc(BED_int.squareSystem_int.size_r * sizeof(int));
    MPI_Bcast(new_deg, BED_int.squareSystem_int.size_r, MPI_INT, headnode, MPI_COMM_WORLD);
    P = (int *)bmalloc(BED_int.squareSystem_int.size_f * sizeof(int));
    MPI_Bcast(P, BED_int.squareSystem_int.size_f, MPI_INT, headnode, MPI_COMM_WORLD);
    size = BED_int.squareSystem_int.size_r * (BED_int.squareSystem_int.size_f - BED_int.squareSystem_int.size_r);
    W = (int *)bmalloc(size * sizeof(int));
    MPI_Bcast(W, size, MPI_INT, headnode, MPI_COMM_WORLD);
    sq_coeff = (comp_d *)bmalloc(BED_int.squareSystem_int.num_comp_d * sizeof(comp_d));
    MPI_Bcast(sq_coeff, BED_int.squareSystem_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // recv the start system structures
    startDeg = (int *)bmalloc(BED_int.startSystem_int.size_r * sizeof(int));
    MPI_Bcast(startDeg, BED_int.startSystem_int.size_r, MPI_INT, headnode, MPI_COMM_WORLD);
    st_coeff = (comp_d *)bmalloc(BED_int.startSystem_int.num_coeff * sizeof(comp_d));
    MPI_Bcast(st_coeff, BED_int.startSystem_int.num_coeff, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // recv the preproc data
    size = BED_int.preProcData_int.num_hom_var_gp + BED_int.preProcData_int.num_var_gp;
    ppd_type = (int *)bmalloc(size * sizeof(int));
    MPI_Bcast(ppd_type, size, MPI_INT, headnode, MPI_COMM_WORLD);
    ppd_size = (int *)bmalloc(size * sizeof(int));
    MPI_Bcast(ppd_size, size, MPI_INT, headnode, MPI_COMM_WORLD);

    // recv the patch structures
    size = BED_int.patch.patchCoeff_rows * BED_int.patch.patchCoeff_cols;
    patch_coeff = (comp_d *)bmalloc(size * sizeof(comp_d));
    MPI_Bcast(patch_coeff, size, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // copy to BED - all data is freed in this function call!!!
    cp_basic_eval_data_d_int(BED, &BED_int, &prog_inst, &prog_gp_sizes, &progStr, 1, &orig_deg, &new_deg, &P, &W, &startDeg, &st_coeff, &sq_coeff, &patch_coeff, &ppd_type, &ppd_size, MPType, 1);
  }

  // free mpi_bed_int, mpi_comp_d
  MPI_Type_free(&mpi_bed_int);
  MPI_Type_free(&mpi_comp_d);
  
  return;
}

void bcast_basic_eval_data_mp(basic_eval_data_mp *BED, int setupProg, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: setupProg: 0 - do not setup prog - like for AMP    *
* RETURN VALUES:                                                *
* NOTES: broadcasts basic_eval_data_mp                          *
\***************************************************************/
{
  MPI_Datatype mpi_bed_int;
  basic_eval_data_mp_int BED_int;
  int size = 0, *prog_inst = NULL, *prog_gp_sizes = NULL, *orig_deg = NULL, *new_deg = NULL, *P = NULL, *W = NULL, *startDeg = NULL, *ppd_type = NULL, *ppd_size = NULL;
  char *progStr = NULL, *sqStr = NULL, *patchStr = NULL, *startStr = NULL;

  // create the datatype mpi_bed_int
  create_basic_eval_data_mp_int(&mpi_bed_int);

  if (my_id == headnode)
  { // setup BED_int & the other structures and then broadcast it
    cp_basic_eval_data_mp_int(&BED_int, BED, &prog_inst, &prog_gp_sizes, &progStr, setupProg, &orig_deg, &new_deg, &P, &W, &sqStr, &patchStr, &startDeg, &startStr, &ppd_type, &ppd_size, 1, 0);

    MPI_Bcast(&BED_int, 1, mpi_bed_int, headnode, MPI_COMM_WORLD);

    if (setupProg)
    { // broadcast the prog structures
      MPI_Bcast(prog_inst, BED_int.squareSystem_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
      MPI_Bcast(prog_gp_sizes, BED_int.squareSystem_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
      MPI_Bcast(progStr, BED_int.squareSystem_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    }

    // broadcast the square system structures
    MPI_Bcast(orig_deg, BED_int.squareSystem_int.size_f, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(new_deg, BED_int.squareSystem_int.size_r, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(P, BED_int.squareSystem_int.size_f, MPI_INT, headnode, MPI_COMM_WORLD);
    size = BED_int.squareSystem_int.size_r * (BED_int.squareSystem_int.size_f - BED_int.squareSystem_int.size_r);
    MPI_Bcast(W, size, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(sqStr, BED_int.squareSystem_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // broadcast the patch structures
    MPI_Bcast(patchStr, BED_int.patch_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // broadcast the start system structures
    MPI_Bcast(startDeg, BED_int.startSystem_int.size_r, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(startStr, BED_int.startSystem_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // broadcast the preproc data
    size = BED_int.preProcData_int.num_hom_var_gp + BED_int.preProcData_int.num_var_gp;
    MPI_Bcast(ppd_type, size, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(ppd_size, size, MPI_INT, headnode, MPI_COMM_WORLD);

    // free the memory
    free(prog_inst);
    free(prog_gp_sizes);
    free(progStr);
    free(orig_deg);
    free(new_deg);
    free(P);
    free(W);
    free(sqStr);
    free(patchStr);
    free(startDeg);
    free(startStr);
    free(ppd_type);
    free(ppd_size);
  }
  else // worker process
  { // recv BED_int
    MPI_Bcast(&BED_int, 1, mpi_bed_int, headnode, MPI_COMM_WORLD);
    // allocate and receive the other structures needed to finish setting up BED

    if (setupProg)
    { // recv the prog structures
      prog_inst = (int *)bmalloc(BED_int.squareSystem_int.Prog_int.size * sizeof(int));
      MPI_Bcast(prog_inst, BED_int.squareSystem_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
      prog_gp_sizes = (int *)bmalloc(BED_int.squareSystem_int.Prog_int.num_var_gps * sizeof(int));
      MPI_Bcast(prog_gp_sizes, BED_int.squareSystem_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
      progStr = (char *)bmalloc(BED_int.squareSystem_int.Prog_int.totalLength * sizeof(char));
      MPI_Bcast(progStr, BED_int.squareSystem_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    }

    // recv the square system structures
    orig_deg = (int *)bmalloc(BED_int.squareSystem_int.size_f * sizeof(int));
    MPI_Bcast(orig_deg, BED_int.squareSystem_int.size_f, MPI_INT, headnode, MPI_COMM_WORLD);
    new_deg = (int *)bmalloc(BED_int.squareSystem_int.size_r * sizeof(int));
    MPI_Bcast(new_deg, BED_int.squareSystem_int.size_r, MPI_INT, headnode, MPI_COMM_WORLD);
    P = (int *)bmalloc(BED_int.squareSystem_int.size_f * sizeof(int));
    MPI_Bcast(P, BED_int.squareSystem_int.size_f, MPI_INT, headnode, MPI_COMM_WORLD);
    size = BED_int.squareSystem_int.size_r * (BED_int.squareSystem_int.size_f - BED_int.squareSystem_int.size_r);
    W = (int *)bmalloc(size * sizeof(int));
    MPI_Bcast(W, size, MPI_INT, headnode, MPI_COMM_WORLD);
    sqStr = (char *)bmalloc(BED_int.squareSystem_int.totalLength * sizeof(char));
    MPI_Bcast(sqStr, BED_int.squareSystem_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // recv the patch structures
    patchStr = (char *)bmalloc(BED_int.patch_int.totalLength * sizeof(char));
    MPI_Bcast(patchStr, BED_int.patch_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // recv the start system structures
    startDeg = (int *)bmalloc(BED_int.startSystem_int.size_r * sizeof(int));
    MPI_Bcast(startDeg, BED_int.startSystem_int.size_r, MPI_INT, headnode, MPI_COMM_WORLD);
    startStr = (char *)bmalloc(BED_int.startSystem_int.totalLength * sizeof(char));
    MPI_Bcast(startStr, BED_int.startSystem_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // recv the preproc data
    size = BED_int.preProcData_int.num_hom_var_gp + BED_int.preProcData_int.num_var_gp;
    ppd_type = (int *)bmalloc(size * sizeof(int));
    MPI_Bcast(ppd_type, size, MPI_INT, headnode, MPI_COMM_WORLD);
    ppd_size = (int *)bmalloc(size * sizeof(int));
    MPI_Bcast(ppd_size, size, MPI_INT, headnode, MPI_COMM_WORLD);

    // copy to BED - all data is freed in this function call!!!
    cp_basic_eval_data_mp_int(BED, &BED_int, &prog_inst, &prog_gp_sizes, &progStr, setupProg, &orig_deg, &new_deg, &P, &W, &sqStr, &patchStr, &startDeg, &startStr, &ppd_type, &ppd_size, 1, 1);
  }

  // free mpi_bed_int
  MPI_Type_free(&mpi_bed_int);

  return;
}

void bcast_basic_eval_data_amp(basic_eval_data_d *BED, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts basic_eval_data_d & basic_eval_data_mp (AMP)*
\***************************************************************/
{
  // setup _d structure - MPType is 2 since using AMP
  bcast_basic_eval_data_d(BED, 2, my_id, headnode);

  // setup _mp structure - Prog does not need setup by bcast
  bcast_basic_eval_data_mp(BED->BED_mp, 0, my_id, headnode);

  if (my_id != headnode)
  { // setup Prog in _mp - just pointer to _d Prog
    BED->BED_mp->squareSystem.Prog = BED->squareSystem.Prog;
  }

  return;
}

void bcast_patch_data_d_amp(basic_eval_data_d *BED, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts patch_d if double or patch_mp if AMP        *
\***************************************************************/
{
  if (my_id == headnode)
  { // determine if using _d or _amp
    if (MPType == 0)
    { // using _d
      comp_d *patch_coeff = NULL;
      patch_eval_data_d_int PED_int;
      MPI_Datatype mpi_comp_d, mpi_patch_d_int;

      // setup mpi_comp_d & mpi_patch_d_int
      create_comp_d(&mpi_comp_d);
      create_patch_eval_data_d_int(&mpi_patch_d_int);
      // setup PED_int
      cp_patch_d_int(&PED_int, &BED->patch, &patch_coeff, 0);

      // broadcast patch structures
      MPI_Bcast(&PED_int, 1, mpi_patch_d_int, headnode, MPI_COMM_WORLD);
      MPI_Bcast(patch_coeff, PED_int.patchCoeff_rows * PED_int.patchCoeff_cols, mpi_comp_d, headnode, MPI_COMM_WORLD);

      // free memory
      MPI_Type_free(&mpi_comp_d);
      MPI_Type_free(&mpi_patch_d_int);
      free(patch_coeff);
    }
    else
    { // using _amp
      char *patchStr = NULL;
      patch_eval_data_mp_int PED_int;
      MPI_Datatype mpi_patch_int;

      // setup mpi_patch_int
      create_patch_eval_data_mp_int(&mpi_patch_int);
      // setup PED_int
      cp_patch_mp_int(&PED_int, &BED->BED_mp->patch, &patchStr, 0, 0);

      // send PED_int
      MPI_Bcast(&PED_int, 1, mpi_patch_int, headnode, MPI_COMM_WORLD);
      // send patchStr
      MPI_Bcast(patchStr, PED_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

      // clear memory
      free(patchStr);
      MPI_Type_free(&mpi_patch_int);
    }
  }
  else // worker process
  { // determine if using _d or _amp
    if (MPType == 0)
    { // using _d
      comp_d *patch_coeff = NULL;
      patch_eval_data_d_int PED_int;
      MPI_Datatype mpi_comp_d, mpi_patch_d_int;

      // setup mpi_comp_d & mpi_patch_d_int
      create_comp_d(&mpi_comp_d);
      create_patch_eval_data_d_int(&mpi_patch_d_int);

      // recv patch structures
      MPI_Bcast(&PED_int, 1, mpi_patch_d_int, headnode, MPI_COMM_WORLD);
      // setup patch_coeff
      patch_coeff = (comp_d *)bmalloc(PED_int.patchCoeff_rows * PED_int.patchCoeff_cols * sizeof(comp_d));
      MPI_Bcast(&patch_coeff, PED_int.patchCoeff_rows * PED_int.patchCoeff_cols, mpi_comp_d, headnode, MPI_COMM_WORLD);


      // setup patch
      cp_patch_d_int(&BED->patch, &PED_int, &patch_coeff, 1);

      // free mpi_comp_d & mpi_patch_d_int
      MPI_Type_free(&mpi_comp_d);
      MPI_Type_free(&mpi_patch_d_int);
    }
    else
    { // using _amp
      int i, j;
      char *patchStr = NULL;
      patch_eval_data_mp_int PED_int;
      MPI_Datatype mpi_patch_int;

      // setup mpi_patch_int
      create_patch_eval_data_mp_int(&mpi_patch_int);
      // recv PED_int
      MPI_Bcast(&PED_int, 1, mpi_patch_int, headnode, MPI_COMM_WORLD);
      
      // setup patchStr
      patchStr = (char *)bmalloc(PED_int.totalLength * sizeof(char));
      // recv patchStr
      MPI_Bcast(patchStr, PED_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD); 

      // setup _mp patch
      cp_patch_mp_int(&BED->BED_mp->patch, &PED_int, &patchStr, 1, 1);

      // setup _d patch
      for (i = 0; i < BED->patch.patchCoeff->rows; i++)
        for (j = 0; j < BED->patch.patchCoeff->cols; j++)
        {
          mp_to_d(&BED->patch.patchCoeff->entry[i][j], &BED->BED_mp->patch.patchCoeff->entry[i][j]);
        }

      // free mpi_patch_int (patchStr is freed in cp_patch_mp_int)
      MPI_Type_free(&mpi_patch_int);
    }
  }

  return;
}

void bcast_eqData_t(eqData_t *EqD, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts eqData_t                                    *
\***************************************************************/
{
  MPI_Datatype mpi_eqd_int;
  eqData_t_int EqD_int;
  int *degrees = NULL, *instCount = NULL;
  char *eqdStr = NULL;
  comp_d *coeff_d = NULL;

  // create the datatype mpi_eqd_int
  create_eqData_t_int(&mpi_eqd_int);

  if (my_id == headnode)
  { // setup EqD_int & the other structures and then broadcast it
    cp_eqData_int(&EqD_int, EqD, MPType, &eqdStr, 0, &coeff_d, &degrees, &instCount, 0);
    // send EqD_int
    MPI_Bcast(&EqD_int, 1, mpi_eqd_int, headnode, MPI_COMM_WORLD);
    // send degrees
    MPI_Bcast(degrees, EqD->num_funcs * EqD->num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
    // send instCount
    MPI_Bcast(instCount, EqD_int.numInts, MPI_INT, headnode, MPI_COMM_WORLD);

    if (MPType == 0)
    { // send coeff_d
      MPI_Datatype mpi_comp_d;
      create_comp_d(&mpi_comp_d);

      MPI_Bcast(coeff_d, EqD_int.num_coeff, mpi_comp_d, headnode, MPI_COMM_WORLD);

      MPI_Type_free(&mpi_comp_d);
    }
    else
    { // send eqdStr
      MPI_Bcast(eqdStr, EqD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    }
    free(degrees);
    free(coeff_d);
    free(eqdStr);
    free(instCount);
  }
  else // worker process
  { // recv EqD_int
    MPI_Bcast(&EqD_int, 1, mpi_eqd_int, headnode, MPI_COMM_WORLD);
    // recv degrees
    degrees = (int *)bmalloc(EqD_int.num_funcs * EqD_int.num_var_gps * sizeof(int));
    MPI_Bcast(degrees, EqD_int.num_funcs * EqD_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv instCount
    instCount = (int *)bmalloc(EqD_int.numInts * sizeof(int));
    MPI_Bcast(instCount, EqD_int.numInts, MPI_INT, headnode, MPI_COMM_WORLD);

    if (MPType == 0)
    { // recv coeff_d
      MPI_Datatype mpi_comp_d;
      create_comp_d(&mpi_comp_d);

      coeff_d = (comp_d *)bmalloc(EqD_int.num_coeff * sizeof(comp_d));
      MPI_Bcast(coeff_d, EqD_int.num_coeff, mpi_comp_d, headnode, MPI_COMM_WORLD);

      MPI_Type_free(&mpi_comp_d);
    }
    else
    { // recv eqdStr
      eqdStr = (char *)bmalloc(EqD_int.totalLength * sizeof(char));
      MPI_Bcast(eqdStr, EqD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    }

    // setup EqD - all data is freed in this function call!!!
    cp_eqData_int(EqD, &EqD_int, MPType, &eqdStr, 1, &coeff_d, &degrees, &instCount, 1);
  }

  MPI_Type_free(&mpi_eqd_int);

  return;
}

void bcast_witnessData(eqData_t *EqD, int stage, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts the witness data from 'stage'               *
\***************************************************************/
{
  MPI_Datatype mpi_wd_int;
  witnessData_t_int WD_int;
  char *wdStr = NULL;
  comp_d *coeff_d = NULL;

  // create the datatype mpi_wd_int
  create_witnessData_t_int(&mpi_wd_int);

  if (my_id == headnode)
  { // setup WD_int & the other structures and then broadcast it
    cp_witnessData_int(&WD_int, EqD, stage, MPType, &wdStr, 0, &coeff_d, 0);
    // send WD_int
    MPI_Bcast(&WD_int, 1, mpi_wd_int, headnode, MPI_COMM_WORLD);

    if (MPType == 0)
    { // send coeff_d
      MPI_Datatype mpi_comp_d;
      create_comp_d(&mpi_comp_d);

      MPI_Bcast(coeff_d, WD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);

      MPI_Type_free(&mpi_comp_d);
    }
    else
    { // send wdStr
      MPI_Bcast(wdStr, WD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    }
    free(coeff_d);
    free(wdStr);
  }
  else // worker process
  { // recv WD_int
    MPI_Bcast(&WD_int, 1, mpi_wd_int, headnode, MPI_COMM_WORLD);

    if (MPType == 0)
    { // recv coeff_d
      MPI_Datatype mpi_comp_d;
      create_comp_d(&mpi_comp_d);

      coeff_d = (comp_d *)bmalloc(WD_int.num_comp_d * sizeof(comp_d));
      MPI_Bcast(coeff_d, WD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);

      MPI_Type_free(&mpi_comp_d);
    }
    else
    { // recv wdStr
      wdStr = (char *)bmalloc(WD_int.totalLength * sizeof(char));
      MPI_Bcast(wdStr, WD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    }

    // setup witnessData - all data is freed in this function call!!!
    cp_witnessData_int(EqD, &WD_int, stage, MPType, &wdStr, 1, &coeff_d, 1);
  }

  MPI_Type_free(&mpi_wd_int);

  return;
}

void bcast_stageData(eqData_t *EqD, int stage, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts the stage data from 'stage'                 *
\***************************************************************/
{
  MPI_Datatype mpi_sd_int;
  stageData_t_int SD_int;
  char *sdStr = NULL;
  comp_d *coeff_d = NULL;

  // create the datatype mpi_sd_int
  create_stageData_t_int(&mpi_sd_int);

  if (my_id == headnode)
  { // setup WD_int & the other structures and then broadcast it
    cp_stageData_int(&SD_int, EqD, stage, MPType, &sdStr, 0, &coeff_d, 0);
    // send SD_int
    MPI_Bcast(&SD_int, 1, mpi_sd_int, headnode, MPI_COMM_WORLD);

    if (MPType == 0)
    { // send coeff_d
      MPI_Datatype mpi_comp_d;
      create_comp_d(&mpi_comp_d);

      MPI_Bcast(coeff_d, SD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);

      MPI_Type_free(&mpi_comp_d);
    }
    else
    { // send sdStr
      MPI_Bcast(sdStr, SD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    }
    free(coeff_d);
    free(sdStr);
  }
  else // worker process
  { // recv SD_int
    MPI_Bcast(&SD_int, 1, mpi_sd_int, headnode, MPI_COMM_WORLD);

    if (MPType == 0)
    { // recv coeff_d
      MPI_Datatype mpi_comp_d;
      create_comp_d(&mpi_comp_d);

      coeff_d = (comp_d *)bmalloc(SD_int.num_comp_d * sizeof(comp_d));
      MPI_Bcast(coeff_d, SD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);

      MPI_Type_free(&mpi_comp_d);
    }
    else
    { // recv sdStr
      sdStr = (char *)bmalloc(SD_int.totalLength * sizeof(char));
      MPI_Bcast(sdStr, SD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    }

    // setup stageData - all data is freed in this function call!!!
    cp_stageData_int(EqD, &SD_int, stage, MPType, &sdStr, 1, &coeff_d, 1);
  }

  MPI_Type_free(&mpi_sd_int);

  return;
}

void bcast_mat_d(mat_d A, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts A                                           *
\***************************************************************/
{
  int num_entries;
  MPI_Datatype mpi_mat_d_int, mpi_comp_d;
  mat_d_int A_int;
  comp_d *entries = NULL;

  // create the datatypes mpi_mat_d_int & mpi_comp_d
  create_mat_d_int(&mpi_mat_d_int);
  create_comp_d(&mpi_comp_d);

  if (my_id == headnode)
  { // setup A_int and entries
    cp_mat_d_int(&A_int, A, &entries, 0);
    num_entries = A->rows * A->cols;

    // send A_int
    MPI_Bcast(&A_int, 1, mpi_mat_d_int, headnode, MPI_COMM_WORLD);
    // send entries
    MPI_Bcast(entries, num_entries, mpi_comp_d, headnode, MPI_COMM_WORLD);
  
    // clear entries
    free(entries);
  }
  else
  { // recv A_int
    MPI_Bcast(&A_int, 1, mpi_mat_d_int, headnode, MPI_COMM_WORLD);
    // setup A and entries
    init_mat_d(A, A_int.rows, A_int.cols);
 
    num_entries = A_int.rows * A_int.cols;
    entries = (comp_d *)bmalloc(num_entries * sizeof(comp_d));
    // recv entries
    MPI_Bcast(entries, num_entries, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // setup A
    cp_mat_d_int(A, &A_int, &entries, 1);
  }

  // clear mpi_mat_d_int & mpi_comp_d
  MPI_Type_free(&mpi_mat_d_int);
  MPI_Type_free(&mpi_comp_d);

  return;
}

void bcast_mat_mp(mat_mp A, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts A                                           *
\***************************************************************/
{
  MPI_Datatype mpi_mat_mp_int;
  mat_mp_int A_int;
  char *Astr = NULL;

  // create the datatypes mpi_mat_mp_int
  create_mat_mp_int(&mpi_mat_mp_int);

  if (my_id == headnode)
  { // setup A_int and Astr
    cp_mat_mp_int(&A_int, A, &Astr, 1, 0);

    // send A_int and Astr
    MPI_Bcast(&A_int, 1, mpi_mat_mp_int, headnode, MPI_COMM_WORLD);
    MPI_Bcast(Astr, A_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
  
    // clear Astr
    free(Astr);
  }
  else
  { // recv A_int and Astr
    MPI_Bcast(&A_int, 1, mpi_mat_mp_int, headnode, MPI_COMM_WORLD);
    Astr = (char *)bmalloc(A_int.totalLength * sizeof(char));
    MPI_Bcast(Astr, A_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // setup A and clear Astr
    cp_mat_mp_int(A, &A_int, &Astr, 1, 1);
  }

  // clear mpi_mat_mp_int
  MPI_Type_free(&mpi_mat_mp_int);

  return;
}

void bcast_mat_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts A                                           *
\***************************************************************/
{
  MPI_Datatype mpi_mat_rat;
  mat_rat_int A_int;
  int rows, cols;
  char *ratStr = NULL;

  // create the datatype mpi_mat_rat
  create_mat_rat_int(&mpi_mat_rat);

  if (my_id == headnode)
  { // setup A_int & ratStr
    rows = A_d->rows;
    cols = A_d->cols;
    cp_mat_rat_int(&A_int, A_rat, &ratStr, rows, cols, 1, 0);

    // send A_int & ratStr
    MPI_Bcast(&A_int, 1, mpi_mat_rat, headnode, MPI_COMM_WORLD);
    MPI_Bcast(ratStr, A_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // clear ratStr
    free(ratStr);
  }
  else
  { // recv A_int & ratStr
    MPI_Bcast(&A_int, 1, mpi_mat_rat, headnode, MPI_COMM_WORLD);
    ratStr = (char *)bmalloc(A_int.totalLength * sizeof(char));
    MPI_Bcast(ratStr, A_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // setup A_rat and clear ratStr
    cp_mat_rat_int(A_rat, &A_int, &ratStr, A_int.rows, A_int.cols, 1, 1);

    // setup A_d & A_mp
    for (rows = 0; rows < A_int.rows; rows++)
      for (cols = 0; cols < A_int.cols; cols++)
      {
        mpf_set_q(A_mp->entry[rows][cols].r, A_rat[rows][cols][0]);
        mpf_set_q(A_mp->entry[rows][cols].i, A_rat[rows][cols][1]);
        A_d->entry[rows][cols].r = mpq_get_d(A_rat[rows][cols][0]);
        A_d->entry[rows][cols].i = mpq_get_d(A_rat[rows][cols][1]);
      }
  }

  // clear mpi_mat_rat
  MPI_Type_free(&mpi_mat_rat);

  return;
}

void bcast_vec_d(vec_d b, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts b                                           *
\***************************************************************/
{
  MPI_Datatype mpi_point_d_int, mpi_comp_d;
  point_d_int b_int;
  comp_d *entries = NULL;

  // create the datatype mpi_point_d_int & mpi_comp_d
  create_point_d_int(&mpi_point_d_int);
  create_comp_d(&mpi_comp_d);

  if (my_id == headnode)
  { // setup b_int and entries
    cp_point_d_int(&b_int, b, &entries, 0, 0, 0);

    // send b_int
    MPI_Bcast(&b_int, 1, mpi_point_d_int, headnode, MPI_COMM_WORLD);
    // send entries
    MPI_Bcast(entries, b_int.size, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // clear entries
    free(entries);
  }
  else
  { // recv b_int
    MPI_Bcast(&b_int, 1, mpi_point_d_int, headnode, MPI_COMM_WORLD);

    entries = (comp_d *)bmalloc(b_int.size * sizeof(comp_d));
    // recv entries
    MPI_Bcast(entries, b_int.size, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // setup b
    cp_point_d_int(b, &b_int, &entries, 1, 1, 1);
  }

  // clear mpi_point_d_int & mpi_comp_d
  MPI_Type_free(&mpi_point_d_int);
  MPI_Type_free(&mpi_comp_d);

  return;
}

void bcast_vec_mp(vec_mp b, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts b                                           *
\***************************************************************/
{
  MPI_Datatype mpi_vec_mp_int;
  point_mp_int b_int;
  char *bstr = NULL;

  // create the datatypes mpi_vec_mp_int
  create_point_mp_int(&mpi_vec_mp_int);

  if (my_id == headnode)
  { // setup b_int and bstr
    cp_point_mp_int(&b_int, b, &bstr, 0, 0, 0);

    // send b_int and bstr
    MPI_Bcast(&b_int, 1, mpi_vec_mp_int, headnode, MPI_COMM_WORLD);
    MPI_Bcast(bstr, b_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // clear bstr
    free(bstr);
  }
  else
  { // recv b_int and bstr
    MPI_Bcast(&b_int, 1, mpi_vec_mp_int, headnode, MPI_COMM_WORLD);
    bstr = (char *)bmalloc(b_int.totalLength * sizeof(char));
    MPI_Bcast(bstr, b_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // setup b and clear bstr
    cp_point_mp_int(b, &b_int, &bstr, 1, 1, 1);
  }

  // clear mpi_vec_mp_int
  MPI_Type_free(&mpi_vec_mp_int);

  return;
}

void bcast_vec_rat(mpq_t ***b, int size, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts A                                           *
\***************************************************************/
{
  MPI_Datatype mpi_point_rat;
  point_rat_int b_int;
  char *ratStr = NULL;

  // create the datatype mpi_point_rat
  create_point_rat_int(&mpi_point_rat);

  if (my_id == headnode)
  { // setup b_int & ratStr
    b_int.size = size;
    cp_vec_rat_char(&ratStr, b, &b_int.totalLength, size, 0, 0);

    // send b_int
    MPI_Bcast(&b_int, 1, mpi_point_rat, headnode, MPI_COMM_WORLD);

    // send ratStr
    MPI_Bcast(ratStr, b_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // clear ratStr
    free(ratStr);
  }
  else
  { // recv b_int
    MPI_Bcast(&b_int, 1, mpi_point_rat, headnode, MPI_COMM_WORLD);

    // setup & recv ratStr
    ratStr = (char *)bmalloc(b_int.totalLength * sizeof(char));
    MPI_Bcast(ratStr, b_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // setup b - clears all structures
    cp_vec_rat_char(b, &ratStr, &b_int.totalLength, b_int.size, 1, 1);
  }

  // clear mpi_point_rat
  MPI_Type_free(&mpi_point_rat);

  return;
}

void bcast_comp_d(comp_d c, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts c                                           *
\***************************************************************/
{
  MPI_Datatype mpi_comp_d;
  create_comp_d(&mpi_comp_d);

  MPI_Bcast(c, 1, mpi_comp_d, headnode, MPI_COMM_WORLD);

  MPI_Type_free(&mpi_comp_d);

  return;
}

void bcast_comp_num_d(comp_d *c, int num, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts c                                           *
\***************************************************************/
{
  MPI_Datatype mpi_comp_d;
  create_comp_d(&mpi_comp_d);

  MPI_Bcast(c, num, mpi_comp_d, headnode, MPI_COMM_WORLD);

  MPI_Type_free(&mpi_comp_d);

  return;
}

void bcast_comp_mp(comp_mp c, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts c                                           *
\***************************************************************/
{
  char *str = NULL;
  comp_mp_int c_int;
  MPI_Datatype mpi_comp_mp_int;
  create_comp_mp_int(&mpi_comp_mp_int);

  if (my_id == headnode)
  { // send data
    cp_comp_mp_int(&c_int, c, &str, 0, 0);    
    // send c_int
    MPI_Bcast(&c_int, 1, mpi_comp_mp_int, headnode, MPI_COMM_WORLD);
    // send str
    MPI_Bcast(str, c_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // clear str
    free(str);
  }
  else
  { // recv data
    MPI_Bcast(&c_int, 1, mpi_comp_mp_int, headnode, MPI_COMM_WORLD);
    // setup & recv str
    str = (char *)bmalloc(c_int.totalLength * sizeof(char));
    MPI_Bcast(str, c_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // setup c
    cp_comp_mp_int(c, &c_int, &str, 1, 1);
  }

  MPI_Type_free(&mpi_comp_mp_int);

  return;
}

void bcast_comp_num_mp(comp_mp *c, int num, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts c                                           *
\***************************************************************/
{
  int i, j, total = 0, currLoc = 0;
  comp_mp_int *c_int = (comp_mp_int *)bmalloc(num * sizeof(comp_mp_int));
  char *str = NULL, *tempStr = NULL;
  MPI_Datatype mpi_comp_mp_int;
  create_comp_mp_int(&mpi_comp_mp_int);

  if (my_id == headnode)
  { // send data
    for (i = 0; i < num; i++)
    { // setup c_int[i]
      cp_comp_mp_int(&c_int[i], &c[i], &tempStr, 0, 0);
      // update
      total += c_int[i].totalLength;
      str = (char *)brealloc(str, total * sizeof(char));
      for (j = 0; j < c_int[i].totalLength; j++)
      {
        str[currLoc] = tempStr[j];
        currLoc++;
      }
      free(tempStr);
    }
    // send c_int
    MPI_Bcast(c_int, num, mpi_comp_mp_int, headnode, MPI_COMM_WORLD);
    // send str
    MPI_Bcast(str, total, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // clear data
    free(str);
    tempStr = NULL;
  }
  else
  { // recv data
    MPI_Bcast(c_int, num, mpi_comp_mp_int, headnode, MPI_COMM_WORLD);
    // setup & recv str
    for (i = 0; i < num; i++)
      total += c_int[i].totalLength;
    str = (char *)bmalloc(total * sizeof(char));
    MPI_Bcast(str, total, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // setup c
    for (i = 0; i < num; i++)
    { // setup c[i]
      tempStr = &str[currLoc];
      cp_comp_mp_int(&c[i], &c_int[i], &tempStr, 0, 1);
      // update currLoc
      currLoc += c_int[i].totalLength;
    }

    // clear data
    free(str);
    tempStr = NULL;
  }

  // free data
  free(c_int);
  MPI_Type_free(&mpi_comp_mp_int);

  return;
}

void bcast_comp_rat(mpq_t c[2], int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts c                                           *
\***************************************************************/
{
  comp_rat_int c_int;
  char *str = NULL;
  MPI_Datatype mpi_comp_rat_int;
  create_comp_rat_int(&mpi_comp_rat_int);

  if (my_id == headnode)
  { // send data
    cp_comp_rat_int(&c_int, c, &str, 0, 0);
    // send c_int
    MPI_Bcast(&c_int, 1, mpi_comp_rat_int, headnode, MPI_COMM_WORLD);
    // send str
    MPI_Bcast(str, c_int.length[0] + c_int.length[1], MPI_CHAR, headnode, MPI_COMM_WORLD);

    // clear str
    free(str);
  }
  else
  { // recv data
    MPI_Bcast(&c_int, 1, mpi_comp_rat_int, headnode, MPI_COMM_WORLD);
    // setup & recv str
    str = (char *)bmalloc((c_int.length[0] + c_int.length[1]) * sizeof(char));
    MPI_Bcast(str, c_int.length[0] + c_int.length[1], MPI_CHAR, headnode, MPI_COMM_WORLD);

    // setup c
    cp_comp_rat_int(c, &c_int, &str, 1, 1);
  }

  MPI_Type_free(&mpi_comp_rat_int);

  return;
}

void bcast_comp_num_rat(mpq_t c[][2], int num, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts c                                           *
\***************************************************************/
{
  int i, j, total = 0, currLoc = 0;
  comp_rat_int *c_int = (comp_rat_int *)bmalloc(num * sizeof(comp_rat_int));
  char *str = NULL, *tempStr = NULL;
  MPI_Datatype mpi_comp_rat_int;
  create_comp_rat_int(&mpi_comp_rat_int);

  if (my_id == headnode)
  { // send data
    for (i = 0; i < num; i++)
    { // setup c_int[i]
      cp_comp_rat_int(&c_int[i], c[i], &tempStr, 0, 0);
      // update
      total += c_int[i].length[0] + c_int[i].length[1];
      str = (char *)brealloc(str, total * sizeof(char));
      for (j = 0; j < c_int[i].length[0] + c_int[i].length[1]; j++)
      {
        str[currLoc] = tempStr[j];
        currLoc++;
      }
      free(tempStr);
    }
    // send c_int
    MPI_Bcast(c_int, num, mpi_comp_rat_int, headnode, MPI_COMM_WORLD);
    // send str
    MPI_Bcast(str, total, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // clear str
    free(str);
  }
  else
  { // recv data
    MPI_Bcast(c_int, num, mpi_comp_rat_int, headnode, MPI_COMM_WORLD);
    // setup & recv str
    for (i = 0; i < num; i++)
      total += c_int[i].length[0] + c_int[i].length[1];
    str = (char *)bmalloc(total * sizeof(char));
    MPI_Bcast(str, total, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // setup c
    for (i = 0; i < num; i++)
    { // setup c[i]
      tempStr = &str[currLoc];
      cp_comp_rat_int(c[i], &c_int[i], &tempStr, 0, 1);
      // update currLoc
      currLoc += c_int[i].length[0] + c_int[i].length[1];
    }
  }

  MPI_Type_free(&mpi_comp_rat_int);

  return;
}

void bcast_prog_t(prog_t *Prog, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts CD                                          *
\***************************************************************/
{
  prog_t_int Prog_int;
  int *inst = NULL, *gp_sizes = NULL;
  char *progStr = NULL;

  MPI_Datatype mpi_prog;

  create_prog_t_int(&mpi_prog);

  if (my_id == headnode)
  { // setup data structs
    cp_prog_t_int(&Prog_int, Prog, &inst, &gp_sizes, &progStr, 0, 0);
    // send Prog_int
    MPI_Bcast(&Prog_int, 1, mpi_prog, headnode, MPI_COMM_WORLD);
    // send inst
    MPI_Bcast(inst, Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
    // send gp_sizes
    MPI_Bcast(gp_sizes, Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
    // send progStr
    MPI_Bcast(progStr, Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // clear memory
    free(inst);
    free(gp_sizes);
    free(progStr);
  }
  else
  { // recv data structs
    MPI_Bcast(&Prog_int, 1, mpi_prog, headnode, MPI_COMM_WORLD);
    // recv inst
    inst = (int *)bmalloc(Prog_int.size * sizeof(int));
    MPI_Bcast(inst, Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv gp_sizes
    gp_sizes = (int *)bmalloc(Prog_int.num_var_gps * sizeof(int));
    MPI_Bcast(gp_sizes, Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv progStr
    progStr = (char *)bmalloc(Prog_int.totalLength * sizeof(char));
    MPI_Bcast(progStr, Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // setup Prog
    cp_prog_t_int(Prog, &Prog_int, &inst, &gp_sizes, &progStr, 1, 1);
  }

  MPI_Type_free(&mpi_prog);

  return;
}

void bcast_codim_t(codim_t *CD, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts CD                                          *
\***************************************************************/
{
  MPI_Datatype mpi_codim_int, mpi_comp_d;
  codim_t_int CD_int;
  int size, *degrees = NULL, *ppd_type = NULL, *ppd_size = NULL, *prog_inst = NULL, *prog_gp_sizes = NULL;
  char *codimStr = NULL, *progStr = NULL;
  comp_d *coeff_d = NULL;

  // create the datatypes mpi_codim_int & mpi_comp_d
  create_codim_t_int(&mpi_codim_int);
  create_comp_d(&mpi_comp_d);

  if (my_id == headnode)
  { // setup CD_int & the other structures and then broadcast it
    cp_codim_int(&CD_int, CD, MPType, &codimStr, &progStr, 0, &coeff_d, &degrees, &ppd_type, &ppd_size, &prog_inst, &prog_gp_sizes, 0);
    // send CD_int
    MPI_Bcast(&CD_int, 1, mpi_codim_int, headnode, MPI_COMM_WORLD);
    // send codimStr
    MPI_Bcast(codimStr, CD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send progStr
    MPI_Bcast(progStr, CD_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send coeff_d
    MPI_Bcast(coeff_d, CD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // send degrees
    size = 3 * CD_int.num_funcs;
    MPI_Bcast(degrees, size, MPI_INT, headnode, MPI_COMM_WORLD);
    // send ppd_type & ppd_size
    size = CD_int.PPD_int.num_hom_var_gp + CD_int.PPD_int.num_var_gp;
    MPI_Bcast(ppd_type, size, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(ppd_size, size, MPI_INT, headnode, MPI_COMM_WORLD);
    // send prog_inst
    MPI_Bcast(prog_inst, CD_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
    // send prog_gp_sizes
    MPI_Bcast(prog_gp_sizes, CD_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);

    // clear memory
    free(degrees);
    free(ppd_type);
    free(ppd_size);
    free(prog_inst);
    free(prog_gp_sizes);
    free(codimStr);
    free(progStr);
    free(coeff_d);
  }
  else
  { // recv CD_int
    MPI_Bcast(&CD_int, 1, mpi_codim_int, headnode, MPI_COMM_WORLD);
    // allocate and recv other structures

    // recv codimStr
    codimStr = (char *)bmalloc(CD_int.totalLength * sizeof(char));
    MPI_Bcast(codimStr, CD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv progStr
    progStr = (char *)bmalloc(CD_int.Prog_int.totalLength * sizeof(char));
    MPI_Bcast(progStr, CD_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv coeff_d
    coeff_d = (comp_d *)bmalloc(CD_int.num_comp_d * sizeof(comp_d));
    MPI_Bcast(coeff_d, CD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // recv degrees
    size = 3 * CD_int.num_funcs;
    degrees = (int *)bmalloc(size * sizeof(int));
    MPI_Bcast(degrees, size, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv ppd_type & ppd_size
    size = CD_int.PPD_int.num_hom_var_gp + CD_int.PPD_int.num_var_gp;
    ppd_type = (int *)bmalloc(size * sizeof(int));
    ppd_size = (int *)bmalloc(size * sizeof(int));
    MPI_Bcast(ppd_type, size, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(ppd_size, size, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv prog_inst
    prog_inst = (int *)bmalloc(CD_int.Prog_int.size * sizeof(int));
    MPI_Bcast(prog_inst, CD_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv prog_gp_sizes
    prog_gp_sizes = (int *)bmalloc(CD_int.Prog_int.num_var_gps * sizeof(int));
    MPI_Bcast(prog_gp_sizes, CD_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);

    // setup CD - all data is freed in this function call!
    cp_codim_int(CD, &CD_int, MPType, &codimStr, &progStr, 1, &coeff_d, &degrees, &ppd_type, &ppd_size, &prog_inst, &prog_gp_sizes, 1);
  }

  // free mpi_codim_int & mpi_comp_d
  MPI_Type_free(&mpi_codim_int);
  MPI_Type_free(&mpi_comp_d);

  return;
}

void bcast_codimData_t(codimData_t *CD, int curr_prec, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts CD                                          *
\***************************************************************/
{
  MPI_Datatype mpi_codimData_int, mpi_comp_d;
  codimData_t_int CD_int;
  int *W_send = NULL;
  char *codimStr = NULL;
  comp_d *coeff_d = NULL;

  // create the datatypes mpi_codimData_int & mpi_comp_d
  create_codimData_t_int(&mpi_codimData_int);
  create_comp_d(&mpi_comp_d);

  if (my_id == headnode)
  { // setup CD_int & the other structures and then broadcast it
    cp_codimData_int(&CD_int, CD, curr_prec, MPType, &codimStr, 0, &coeff_d, &W_send, 0);
    // send CD_int
    MPI_Bcast(&CD_int, 1, mpi_codimData_int, headnode, MPI_COMM_WORLD);
    // send codimStr
    MPI_Bcast(codimStr, CD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send coeff_d
    MPI_Bcast(coeff_d, CD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // send W_send
    MPI_Bcast(W_send, CD_int.A_W_rows * CD_int.A_W_cols, MPI_INT, headnode, MPI_COMM_WORLD);

    // clear memory
    free(codimStr);
    free(coeff_d);
  }
  else
  { // recv CD_int
    MPI_Bcast(&CD_int, 1, mpi_codimData_int, headnode, MPI_COMM_WORLD);
    // allocate and recv other structures

    // recv codimStr
    codimStr = (char *)bmalloc(CD_int.totalLength * sizeof(char));
    MPI_Bcast(codimStr, CD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv coeff_d
    coeff_d = (comp_d *)bmalloc(CD_int.num_comp_d * sizeof(comp_d));
    MPI_Bcast(coeff_d, CD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // recv W_send
    W_send = (int *)bmalloc(CD_int.A_W_rows * CD_int.A_W_cols * sizeof(int));
    MPI_Bcast(W_send, CD_int.A_W_rows * CD_int.A_W_cols, MPI_INT, headnode, MPI_COMM_WORLD);

    // setup CD - all data is freed in this function call!
    cp_codimData_int(CD, &CD_int, curr_prec, MPType, &codimStr, 1, &coeff_d, &W_send, 1);
  }

  // free mpi_codimData_int & mpi_comp_d
  MPI_Type_free(&mpi_codimData_int);
  MPI_Type_free(&mpi_comp_d);

  return;
}

void bcast_cascade_t(cascade_t *CD, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts CD                                          *
\***************************************************************/
{
  MPI_Datatype mpi_cascade_int, mpi_comp_d;
  cascade_t_int CD_int;
  int size, *degrees = NULL, *ppd_type = NULL, *ppd_size = NULL, *prog_inst = NULL, *prog_gp_sizes = NULL;
  char *cascadeStr = NULL, *progStr = NULL;
  comp_d *coeff_d = NULL;

  // create the datatypes mpi_cascade_int & mpi_comp_d
  create_cascade_t_int(&mpi_cascade_int);
  create_comp_d(&mpi_comp_d);

  if (my_id == headnode)
  { // setup CD_int & the other structures and then broadcast it
    cp_cascade_int(&CD_int, CD, MPType, &cascadeStr, &progStr, 0, &coeff_d, &degrees, &ppd_type, &ppd_size, &prog_inst, &prog_gp_sizes, 0);
    // send CD_int
    MPI_Bcast(&CD_int, 1, mpi_cascade_int, headnode, MPI_COMM_WORLD);
    // send cascadeStr
    MPI_Bcast(cascadeStr, CD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send progStr
    MPI_Bcast(progStr, CD_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send coeff_d
    MPI_Bcast(coeff_d, CD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // send degrees
    MPI_Bcast(degrees, CD_int.num_int, MPI_INT, headnode, MPI_COMM_WORLD);
    // send ppd_type & ppd_size
    size = CD_int.PPD_int.num_hom_var_gp + CD_int.PPD_int.num_var_gp;
    MPI_Bcast(ppd_type, size, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(ppd_size, size, MPI_INT, headnode, MPI_COMM_WORLD);
    // send prog_inst
    MPI_Bcast(prog_inst, CD_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
    // send prog_gp_sizes
    MPI_Bcast(prog_gp_sizes, CD_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);

    // clear memory
    free(degrees);
    free(ppd_type);
    free(ppd_size);
    free(prog_inst);
    free(prog_gp_sizes);
    free(cascadeStr);
    free(progStr);
    free(coeff_d);
  }
  else
  { // recv CD_int
    MPI_Bcast(&CD_int, 1, mpi_cascade_int, headnode, MPI_COMM_WORLD);
    // allocate and recv other structures

    // recv cascadeStr
    cascadeStr = (char *)bmalloc(CD_int.totalLength * sizeof(char));
    MPI_Bcast(cascadeStr, CD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv progStr
    progStr = (char *)bmalloc(CD_int.Prog_int.totalLength * sizeof(char));
    MPI_Bcast(progStr, CD_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv coeff_d
    coeff_d = (comp_d *)bmalloc(CD_int.num_comp_d * sizeof(comp_d));
    MPI_Bcast(coeff_d, CD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // recv degrees
    degrees = (int *)bmalloc(CD_int.num_int * sizeof(int));
    MPI_Bcast(degrees, CD_int.num_int, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv ppd_type & ppd_size
    size = CD_int.PPD_int.num_hom_var_gp + CD_int.PPD_int.num_var_gp;
    ppd_type = (int *)bmalloc(size * sizeof(int));
    ppd_size = (int *)bmalloc(size * sizeof(int));
    MPI_Bcast(ppd_type, size, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(ppd_size, size, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv prog_inst
    prog_inst = (int *)bmalloc(CD_int.Prog_int.size * sizeof(int));
    MPI_Bcast(prog_inst, CD_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv prog_gp_sizes
    prog_gp_sizes = (int *)bmalloc(CD_int.Prog_int.num_var_gps * sizeof(int));
    MPI_Bcast(prog_gp_sizes, CD_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);

    // setup CD - all data is freed in this function call!
    cp_cascade_int(CD, &CD_int, MPType, &cascadeStr, &progStr, 1, &coeff_d, &degrees, &ppd_type, &ppd_size, &prog_inst, &prog_gp_sizes, 1);
  }

  // free mpi_cascade_int & mpi_comp_d
  MPI_Type_free(&mpi_cascade_int);
  MPI_Type_free(&mpi_comp_d);
  
  return;
}

void bcast_cascadeCodim_t(cascadeCodim_t *CD, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts CD                                          *
\***************************************************************/
{
  MPI_Datatype mpi_cascadeCodim_int;
  cascadeCodim_t_int CD_int;

  // create the datatypes mpi_cascadeCodim_int 
  create_cascadeCodim_t_int(&mpi_cascadeCodim_int);

  if (my_id == headnode)
  { // setup CD_int
    cp_cascadeCodim_int(&CD_int, CD, 0);

    // send CD_int
    MPI_Bcast(&CD_int, 1, mpi_cascadeCodim_int, headnode, MPI_COMM_WORLD);
  }
  else
  { // recv CD_int
    MPI_Bcast(&CD_int, 1, mpi_cascadeCodim_int, headnode, MPI_COMM_WORLD);

    // setup CD
    cp_cascadeCodim_int(CD, &CD_int, 1);
  }

  // free mpi_cascadeCodim_int
  MPI_Type_free(&mpi_cascadeCodim_int);

  return;
}

void bcast_membership_slice_moving_t(membership_slice_moving_t *slice, int sendProg, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts slice                                       *
\***************************************************************/
{
  MPI_Datatype mpi_slice_int, mpi_comp_d;
  membership_slice_moving_t_int slice_int;

  int *prog_inst = NULL, *prog_gp_sizes = NULL;
  char *sliceStr = NULL, *progStr = NULL;
  comp_d *coeff_d = NULL;

  // create the datatype mpi_slice_int & mpi_comp_d
  create_membership_slice_moving_t_int(&mpi_slice_int);
  create_comp_d(&mpi_comp_d);

  if (my_id == headnode)
  { // setup slice_int
    cp_slice_moving_int(&slice_int, &sliceStr, &progStr, &prog_inst, &prog_gp_sizes, &coeff_d, slice, sendProg, MPType, 0, 0);
    // send slice_int
    MPI_Bcast(&slice_int, 1, mpi_slice_int, headnode, MPI_COMM_WORLD);
    // send sliceStr
    MPI_Bcast(sliceStr, slice_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send progStr
    MPI_Bcast(progStr, slice_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send coeff_d
    MPI_Bcast(coeff_d, slice_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // free memory
    free(sliceStr);
    free(progStr);
    free(prog_inst);
    free(prog_gp_sizes);
    free(coeff_d);
  }
  else
  { // recv slice_int
    MPI_Bcast(&slice_int, 1, mpi_slice_int, headnode, MPI_COMM_WORLD);
    // recv sliceStr
    sliceStr = (char *)bmalloc(slice_int.totalLength * sizeof(char));
    MPI_Bcast(sliceStr, slice_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv progStr
    progStr = (char *)bmalloc(slice_int.Prog_int.totalLength * sizeof(char));
    MPI_Bcast(progStr, slice_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv coeff_d
    coeff_d = (comp_d *)bmalloc(slice_int.num_comp_d * sizeof(comp_d));
    MPI_Bcast(coeff_d, slice_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // setup slice
    cp_slice_moving_int(slice, &sliceStr, &progStr, &prog_inst, &prog_gp_sizes, &coeff_d, &slice_int, sendProg, MPType, 1, 1);
  }

  // free mpi_slice_int & mpi_comp_d
  MPI_Type_free(&mpi_slice_int);
  MPI_Type_free(&mpi_comp_d);

  return;
}

void bcast_witnessCodim_t(witnessCodim_t *witCodim, int MPType, int curr_prec, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts witCodim                                    *
\***************************************************************/
{
  MPI_Datatype mpi_witnessCodim, mpi_comp_d, mpi_endpoint_d, mpi_endpoint_mp, mpi_endpoint_amp;
  witnessCodim_t_int witCodim_int;

  int *witW = NULL, *witTypes = NULL, *witMult = NULL;
  char *witStr = NULL;
  comp_d *witComp = NULL;
  endpoint_data_d_int *wit_d = NULL;
  endpoint_data_mp_int *wit_mp = NULL;
  endpoint_data_amp_int *wit_amp = NULL;

  // create the datatypes
  create_comp_d(&mpi_comp_d);
  create_endpoint_data_d_int(&mpi_endpoint_d);
  create_endpoint_data_mp_int(&mpi_endpoint_mp);
  create_endpoint_data_amp_int(&mpi_endpoint_amp);
  create_witnessCodim_t_int(&mpi_witnessCodim);

  if (my_id == headnode)
  { // setup witCodim_int
    cp_witnessCodim_t_int(&witCodim_int, witCodim, &witComp, &witStr, &witW, &witTypes, &witMult, &wit_d, &wit_mp, &wit_amp, curr_prec, 0, MPType, 0);
    // send witCodim_int
    MPI_Bcast(&witCodim_int, 1, mpi_witnessCodim, headnode, MPI_COMM_WORLD);
    // setup witComp
    MPI_Bcast(witComp, witCodim_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // send witStr
    MPI_Bcast(witStr, witCodim_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send witW
    MPI_Bcast(witW, witCodim_int.A_rows * witCodim_int.A_cols, MPI_INT, headnode, MPI_COMM_WORLD);
    // send witTypes
    MPI_Bcast(witTypes, witCodim_int.num_set, MPI_INT, headnode, MPI_COMM_WORLD);
    // send witMult
    MPI_Bcast(witMult, witCodim_int.num_set, MPI_INT, headnode, MPI_COMM_WORLD);
    // send wit_d,_mp,_amp
    if (MPType == 0)
      MPI_Bcast(wit_d, witCodim_int.num_set, mpi_endpoint_d, headnode, MPI_COMM_WORLD);
    else if (MPType == 1)
      MPI_Bcast(wit_mp, witCodim_int.num_set, mpi_endpoint_mp, headnode, MPI_COMM_WORLD);
    else
      MPI_Bcast(wit_amp, witCodim_int.num_set, mpi_endpoint_amp, headnode, MPI_COMM_WORLD);

    // free memory
    free(witW);  free(witTypes); free(witMult);
    free(witStr);
    free(witComp);
    free(wit_d); free(wit_mp); free(wit_amp);
  }
  else
  { // recv witCodim_int
    MPI_Bcast(&witCodim_int, 1, mpi_witnessCodim, headnode, MPI_COMM_WORLD);
    // recv witComp
    witComp = (comp_d *)bmalloc(witCodim_int.num_comp_d * sizeof(comp_d));
    MPI_Bcast(witComp, witCodim_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // recv witStr
    witStr = (char *)bmalloc(witCodim_int.totalLength * sizeof(char));
    MPI_Bcast(witStr, witCodim_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv witW
    witW = (int *)bmalloc(witCodim_int.A_rows * witCodim_int.A_cols * sizeof(int));
    MPI_Bcast(witW, witCodim_int.A_rows * witCodim_int.A_cols, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv witTypes
    witTypes = (int *)bmalloc(witCodim_int.num_set * sizeof(int));
    MPI_Bcast(witTypes, witCodim_int.num_set, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv witMult
    witMult = (int *)bmalloc(witCodim_int.num_set * sizeof(int));
    MPI_Bcast(witMult, witCodim_int.num_set, MPI_INT, headnode, MPI_COMM_WORLD);

    // recv wit_d,_mp,_amp
    if (MPType == 0)
    {
      wit_d = (endpoint_data_d_int *)bmalloc(witCodim_int.num_set * sizeof(endpoint_data_d_int));
      MPI_Bcast(wit_d, witCodim_int.num_set, mpi_endpoint_d, headnode, MPI_COMM_WORLD);
    }
    else if (MPType == 1)
    {
      wit_mp = (endpoint_data_mp_int *)bmalloc(witCodim_int.num_set * sizeof(endpoint_data_mp_int));
      MPI_Bcast(wit_mp, witCodim_int.num_set, mpi_endpoint_mp, headnode, MPI_COMM_WORLD);
    }
    else
    {
      wit_amp = (endpoint_data_amp_int *)bmalloc(witCodim_int.num_set * sizeof(endpoint_data_amp_int));
      MPI_Bcast(wit_amp, witCodim_int.num_set, mpi_endpoint_amp, headnode, MPI_COMM_WORLD);
    }

    // setup witCodim
    cp_witnessCodim_t_int(witCodim, &witCodim_int, &witComp, &witStr, &witW, &witTypes, &witMult, &wit_d, &wit_mp, &wit_amp, witCodim_int.curr_prec, 1, MPType, 1);
  }

  // free datatypes
  MPI_Type_free(&mpi_witnessCodim);
  MPI_Type_free(&mpi_comp_d);
  MPI_Type_free(&mpi_endpoint_d);
  MPI_Type_free(&mpi_endpoint_mp);
  MPI_Type_free(&mpi_endpoint_amp);

  return;
}

void bcast_witness_t(witness_t *W, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts W                                           *
\***************************************************************/
{
  int i;
  MPI_Datatype mpi_witness_int, mpi_comp_d;
  witness_t_int W_int;

  int *progInst = NULL, *gpSizes = NULL, *PPDtype = NULL, *PPDsize = NULL, *origDegs = NULL, *newDegs = NULL, *P = NULL;
  char *progStr = NULL, *witStr = NULL;
  comp_d *witComp = NULL;

  // create the datatypes
  create_comp_d(&mpi_comp_d);
  create_witness_t_int(&mpi_witness_int);

  if (my_id == headnode)
  { // setup W_int
    cp_witness_t_int(&W_int, W, &progInst, &gpSizes, &progStr, &PPDtype, &PPDsize, &origDegs, &newDegs, &P, &witComp, &witStr, 0, MPType, 0);
    // send W_int
    MPI_Bcast(&W_int, 1, mpi_witness_int, headnode, MPI_COMM_WORLD);
    // send progInst
    MPI_Bcast(progInst, W_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
    // send gpSizes
    MPI_Bcast(gpSizes, W_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
    // send progStr
    MPI_Bcast(progStr, W_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send PPDtype
    MPI_Bcast(PPDtype, W_int.PPD_int.num_hom_var_gp + W_int.PPD_int.num_var_gp, MPI_INT, headnode, MPI_COMM_WORLD);
    // send PPDsize
    MPI_Bcast(PPDsize, W_int.PPD_int.num_hom_var_gp + W_int.PPD_int.num_var_gp, MPI_INT, headnode, MPI_COMM_WORLD);
    // send origDegs
    MPI_Bcast(origDegs, W_int.num_funcs, MPI_INT, headnode, MPI_COMM_WORLD);
    // send newDegs
    MPI_Bcast(newDegs, W_int.num_funcs, MPI_INT, headnode, MPI_COMM_WORLD);
    // send P
    MPI_Bcast(P, W_int.num_funcs, MPI_INT, headnode, MPI_COMM_WORLD);
    // send witComp
    MPI_Bcast(witComp, W_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // send witStr
    MPI_Bcast(witStr, W_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // send the codimensions
    for (i = 0; i < W_int.num_codim; i++)
    {
      bcast_witnessCodim_t(&W->codim[i], MPType, W->curr_precision, my_id, headnode);
    }

    // free data
    free(progInst);
    free(gpSizes);
    free(PPDtype);
    free(PPDsize);
    free(origDegs);
    free(newDegs);
    free(P);
    free(progStr);
    free(witStr);
    free(witComp);
  }
  else
  { // recv W_int
    MPI_Bcast(&W_int, 1, mpi_witness_int, headnode, MPI_COMM_WORLD);
    // recv progInst
    progInst = (int *)bmalloc(W_int.Prog_int.size * sizeof(int));
    MPI_Bcast(progInst, W_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv gpSizes
    gpSizes = (int *)bmalloc(W_int.Prog_int.num_var_gps * sizeof(int));
    MPI_Bcast(gpSizes, W_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv progStr
    progStr = (char *)bmalloc(W_int.Prog_int.totalLength * sizeof(char));
    MPI_Bcast(progStr, W_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv PPDtype
    PPDtype = (int *)bmalloc((W_int.PPD_int.num_hom_var_gp + W_int.PPD_int.num_var_gp) * sizeof(int));
    MPI_Bcast(PPDtype, W_int.PPD_int.num_hom_var_gp + W_int.PPD_int.num_var_gp, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv PPDsize
    PPDsize = (int *)bmalloc((W_int.PPD_int.num_hom_var_gp + W_int.PPD_int.num_var_gp) * sizeof(int));
    MPI_Bcast(PPDsize, W_int.PPD_int.num_hom_var_gp + W_int.PPD_int.num_var_gp, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv origDegs
    origDegs = (int *)bmalloc(W_int.num_funcs * sizeof(int));
    MPI_Bcast(origDegs, W_int.num_funcs, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv newDegs
    newDegs = (int *)bmalloc(W_int.num_funcs * sizeof(int));
    MPI_Bcast(newDegs, W_int.num_funcs, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv P
    P = (int *)bmalloc(W_int.num_funcs * sizeof(int));
    MPI_Bcast(P, W_int.num_funcs, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv witComp
    witComp = (comp_d *)bmalloc(W_int.num_comp_d * sizeof(comp_d));
    MPI_Bcast(witComp, W_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // recv witStr
    witStr = (char *)bmalloc(W_int.totalLength * sizeof(char));
    MPI_Bcast(witStr, W_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);

    // setup W
    cp_witness_t_int(W, &W_int, &progInst, &gpSizes, &progStr, &PPDtype, &PPDsize, &origDegs, &newDegs, &P, &witComp, &witStr, 1, MPType, 1);

    // setup the codimensions
    W->codim = (witnessCodim_t *)bmalloc(W_int.num_codim * sizeof(witnessCodim_t));
    for (i = 0; i < W_int.num_codim; i++)
    {
      bcast_witnessCodim_t(&W->codim[i], MPType, W->curr_precision, my_id, headnode);
    }
  }

  // free datatypes
  MPI_Type_free(&mpi_witness_int);
  MPI_Type_free(&mpi_comp_d);

  return;
}

void bcast_witness_structures(prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, int MPType, witness_t *W, int codim_index, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts the witness structures for this codim       *
\***************************************************************/
{
  int i, j, stdProg = -1, compCount = 0, compLoc = 0, totalCount = 0, currLoc = 0;
  char *str = NULL, *tempStr = NULL;
  comp_d *coeff = NULL, *tempCoeff = NULL;
  endpoint_data_d_int *wit_d = NULL;
  endpoint_data_mp_int *wit_mp = NULL;
  endpoint_data_amp_int *wit_amp = NULL;
  MPI_Datatype mpi_comp_d, mpi_endpoint_d, mpi_endpoint_mp, mpi_endpoint_amp;

  // create the datatypes
  create_comp_d(&mpi_comp_d);
  create_endpoint_data_d_int(&mpi_endpoint_d);
  create_endpoint_data_mp_int(&mpi_endpoint_mp);
  create_endpoint_data_amp_int(&mpi_endpoint_amp);

  if (my_id == headnode)
  { // send deflations needed
    MPI_Bcast(W->codim[codim_index].deflations_needed, W->codim[codim_index].num_set, MPI_INT, headnode, MPI_COMM_WORLD);
    // send fullRankProgInfo
    MPI_Bcast(*fullRankProgInfo, W->codim[codim_index].num_set, MPI_INT, headnode, MPI_COMM_WORLD);
    // send the endpoints
    if (MPType == 0)
    { // send endPts_d
      wit_d = (endpoint_data_d_int *)bmalloc(W->codim[codim_index].num_set * sizeof(endpoint_data_d_int));
      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // setup wit_d[i]
        cp_endpoint_data_d_int(&wit_d[i], &(*endPts_d)[i], &tempCoeff, 0, 0, 0);
        // copy to coeff
        compCount += wit_d[i].num_comp_d;
        coeff = (comp_d *)brealloc(coeff, compCount * sizeof(comp_d));
        for (j = 0; j < wit_d[i].num_comp_d; j++)
        {
          set_d(coeff[compLoc], tempCoeff[j]); 
          compLoc++;
        }
        free(tempCoeff);
      }
      MPI_Bcast(wit_d, W->codim[codim_index].num_set, mpi_endpoint_d, headnode, MPI_COMM_WORLD);
      MPI_Bcast(coeff, compCount, mpi_comp_d, headnode, MPI_COMM_WORLD);

      // clear the memory
      free(wit_d);
      free(coeff);
      tempCoeff = NULL;
    }
    else if (MPType == 1)
    { // send endPts_mp
      wit_mp = (endpoint_data_mp_int *)bmalloc(W->codim[codim_index].num_set * sizeof(endpoint_data_mp_int));
      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // setup wit_mp[i]
        cp_endpoint_data_mp_int(&wit_mp[i], &(*endPts_mp)[i], &tempStr, 0, 0, 0);
        // copy to str
        totalCount += wit_mp[i].totalLength;
        str = (char *)brealloc(str, totalCount * sizeof(char));
        for (j = 0; j < wit_mp[i].totalLength; j++)
        {
          str[currLoc] = tempStr[j];
          currLoc++;
        }
        free(tempStr);
      }
      MPI_Bcast(wit_mp, W->codim[codim_index].num_set, mpi_endpoint_mp, headnode, MPI_COMM_WORLD);
      MPI_Bcast(str, totalCount, MPI_CHAR, headnode, MPI_COMM_WORLD);

      // clear the memory
      free(wit_mp);
      free(str);
      tempStr = NULL;
    }
    else
    { // send endPts_amp
      wit_amp = (endpoint_data_amp_int *)bmalloc(W->codim[codim_index].num_set * sizeof(endpoint_data_amp_int));
      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // setup wit_amp[i]
        cp_endpoint_data_amp_int(&wit_amp[i], &(*endPts_amp)[i], &tempCoeff, &tempStr, 0, 0, 0);
        // copy to coeff
        compCount += wit_amp[i].num_comp_d;
        coeff = (comp_d *)brealloc(coeff, compCount * sizeof(comp_d));
        for (j = 0; j < wit_amp[i].num_comp_d; j++)
        {
          set_d(coeff[compLoc], tempCoeff[j]);
          compLoc++;
        }
        free(tempCoeff);
        // copy to str
        totalCount += wit_amp[i].totalLength;
        str = (char *)brealloc(str, totalCount * sizeof(char));
        for (j = 0; j < wit_amp[i].totalLength; j++)
        {
          str[currLoc] = tempStr[j];
          currLoc++;
        }
        free(tempStr);
      }
      MPI_Bcast(wit_amp, W->codim[codim_index].num_set, mpi_endpoint_amp, headnode, MPI_COMM_WORLD);
      MPI_Bcast(coeff, compCount, mpi_comp_d, headnode, MPI_COMM_WORLD);
      MPI_Bcast(str, totalCount, MPI_CHAR, headnode, MPI_COMM_WORLD);

      // clear the memory
      free(wit_amp);
      free(coeff);
      free(str);
      tempCoeff = NULL;
      tempStr = NULL;
    }
    // send fullRankProgs
    for (i = 0; i < W->codim[codim_index].num_set; i++)
      if ((*fullRankProgInfo)[i] == 0 || (*fullRankProgInfo)[i] == 1)
      { // see if this one is a deflated or not
        if (W->codim[codim_index].deflations_needed[i] == 0 && stdProg == -1)
        { // send this SLP
          bcast_prog_t((*fullRankProgs)[i], MPType, my_id, headnode);
          // standard SLP has been sent    
          stdProg = i;
        }
        else if (W->codim[codim_index].deflations_needed[i] > 0)
        { // send this SLP
          bcast_prog_t((*fullRankProgs)[i], MPType, my_id, headnode);
        }
      }
  }
  else
  { // setup deflations_needed
    W->codim[codim_index].deflations_needed = (int *)bmalloc(W->codim[codim_index].num_set * sizeof(int));
    MPI_Bcast(W->codim[codim_index].deflations_needed, W->codim[codim_index].num_set, MPI_INT, headnode, MPI_COMM_WORLD);
    // setup fullRankProgInfo
    *fullRankProgInfo = (int *)bmalloc(W->codim[codim_index].num_set * sizeof(int));
    MPI_Bcast(*fullRankProgInfo, W->codim[codim_index].num_set, MPI_INT, headnode, MPI_COMM_WORLD);
    // setup the endpoints
    if (MPType == 0)
    { // setup endPts_d
      *endPts_d = (endpoint_data_d *)bmalloc(W->codim[codim_index].num_set * sizeof(endpoint_data_d));
      // recv wit_d
      wit_d = (endpoint_data_d_int *)bmalloc(W->codim[codim_index].num_set * sizeof(endpoint_data_d_int));
      MPI_Bcast(wit_d, W->codim[codim_index].num_set, mpi_endpoint_d, headnode, MPI_COMM_WORLD);
      // recv coeff
      for (i = 0; i < W->codim[codim_index].num_set; i++)
        compCount += wit_d[i].num_comp_d;
      coeff = (comp_d *)bmalloc(compCount * sizeof(comp_d));
      MPI_Bcast(coeff, compCount, mpi_comp_d, headnode, MPI_COMM_WORLD);
      // setup endPts_d
      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // setup endPts_d[i]
        tempCoeff = &coeff[compLoc];
        cp_endpoint_data_d_int(&(*endPts_d)[i], &wit_d[i], &tempCoeff, 0, 1, 1);
        // update
        compLoc += wit_d[i].num_comp_d;
      }

      // clear the memory
      free(wit_d);
      free(coeff);
      tempCoeff = NULL;
    }
    else if (MPType == 1)
    { // setup endPts_mp
      *endPts_mp = (endpoint_data_mp *)bmalloc(W->codim[codim_index].num_set * sizeof(endpoint_data_mp));
      // recv wit_mp
      wit_mp = (endpoint_data_mp_int *)bmalloc(W->codim[codim_index].num_set * sizeof(endpoint_data_mp_int));
      MPI_Bcast(wit_mp, W->codim[codim_index].num_set, mpi_endpoint_mp, headnode, MPI_COMM_WORLD);
      // recv str
      for (i = 0; i < W->codim[codim_index].num_set; i++)
        totalCount += wit_mp[i].totalLength;
      str = (char *)bmalloc(totalCount * sizeof(char));
      MPI_Bcast(str, totalCount, MPI_CHAR, headnode, MPI_COMM_WORLD);
      // setup endPts_d
      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // setup endPts_d[i]
        tempStr = &str[currLoc];
        cp_endpoint_data_mp_int(&(*endPts_mp)[i], &wit_mp[i], &tempStr, 0, 1, 1);
        // update
        currLoc += wit_mp[i].totalLength;
      }

      // clear the memory
      free(wit_mp);
      free(str);
      tempStr = NULL;
    }
    else
    { // setup endPts_amp
      *endPts_amp = (endpoint_data_amp *)bmalloc(W->codim[codim_index].num_set * sizeof(endpoint_data_amp));
      // recv wit_amp
      wit_amp = (endpoint_data_amp_int *)bmalloc(W->codim[codim_index].num_set * sizeof(endpoint_data_amp_int));
      MPI_Bcast(wit_amp, W->codim[codim_index].num_set, mpi_endpoint_amp, headnode, MPI_COMM_WORLD);
      // recv coeff
      for (i = 0; i < W->codim[codim_index].num_set; i++)
        compCount += wit_amp[i].num_comp_d;
      coeff = (comp_d *)bmalloc(compCount * sizeof(comp_d));
      MPI_Bcast(coeff, compCount, mpi_comp_d, headnode, MPI_COMM_WORLD);
      // recv str
      for (i = 0; i < W->codim[codim_index].num_set; i++)
        totalCount += wit_amp[i].totalLength;
      str = (char *)bmalloc(totalCount * sizeof(char));
      MPI_Bcast(str, totalCount, MPI_CHAR, headnode, MPI_COMM_WORLD);
      // setup endPts_amp
      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // setup endPts_amp[i]
        tempCoeff = &coeff[compLoc];
        tempStr = &str[currLoc];
        cp_endpoint_data_amp_int(&(*endPts_amp)[i], &wit_amp[i], &tempCoeff, &tempStr, 0, 1, 1);
        // update
        compLoc += wit_amp[i].num_comp_d;
        currLoc += wit_amp[i].totalLength;
      }

      // clear the memory
      free(wit_amp);
      free(coeff);
      free(str);
      tempCoeff = NULL;
      tempStr = NULL;
    }
    // setup fullRankProgs
    *fullRankProgs = (prog_t **)bmalloc(W->codim[codim_index].num_set * sizeof(prog_t *));
    for (i = 0; i < W->codim[codim_index].num_set; i++)
      if ((*fullRankProgInfo)[i] == 0 || (*fullRankProgInfo)[i] == 1)
      { // see if this one is a deflated or not
        if (W->codim[codim_index].deflations_needed[i] == 0)
        {
          if (stdProg == -1)
          { // recv this SLP
            (*fullRankProgs)[i] = (prog_t *)bmalloc(1 * sizeof(prog_t));
            bcast_prog_t((*fullRankProgs)[i], MPType, my_id, headnode);
            // standard SLP has been recv'd
            stdProg = i;
            (*fullRankProgInfo)[i] = 1;
          }
          else
          { // point to the standard SLP
            (*fullRankProgs)[i] = (*fullRankProgs)[stdProg];
          }
        }
        else if (W->codim[codim_index].deflations_needed[i] > 0)
        { // recv this SLP
          (*fullRankProgs)[i] = (prog_t *)bmalloc(1 * sizeof(prog_t));
          bcast_prog_t((*fullRankProgs)[i], MPType, my_id, headnode);
        }
      }
      else
      { // NULL it out
        (*fullRankProgs)[i] = NULL;
      }
  }

  // free datatypes
  MPI_Type_free(&mpi_comp_d);
  MPI_Type_free(&mpi_endpoint_d);
  MPI_Type_free(&mpi_endpoint_mp);
  MPI_Type_free(&mpi_endpoint_amp);


  return;
}

void bcast_trace_structures(vec_d proj_d, vec_mp proj_mp, mpq_t ***proj_rat, vec_d v_d, vec_mp v_mp, mpq_t ***v_rat, comp_d s_d[2], comp_mp s_mp[2], mpq_t s_rat[2][2], comp_d gamma_d, comp_mp gamma_mp, mpq_t gamma_rat[2], int MPType, int curr_prec, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts the witness structures for this codim       *
\***************************************************************/
{
  if (MPType == 0)
  { // _d
    bcast_vec_d(proj_d, my_id, headnode);
    bcast_vec_d(v_d, my_id, headnode);
    bcast_comp_num_d(s_d, 2, my_id, headnode);
    bcast_comp_d(gamma_d, my_id, headnode);
  }
  else if (MPType == 1)
  { // _mp 
    bcast_vec_mp(proj_mp, my_id, headnode);
    bcast_vec_mp(v_mp, my_id, headnode);
    bcast_comp_num_mp(s_mp, 2, my_id, headnode);
    bcast_comp_mp(gamma_mp, my_id, headnode);
  }
  else
  { // _rat
    int i, vec_size = proj_d->size;
    MPI_Bcast(&vec_size, 1, MPI_INT, headnode, MPI_COMM_WORLD);
    bcast_vec_rat(proj_rat, vec_size, my_id, headnode);

    if (my_id != headnode)
    { // setup _d & _mp
      init_vec_d(proj_d, vec_size);
      init_vec_mp2(proj_mp, vec_size, curr_prec);
      proj_d->size = proj_mp->size = vec_size;

      for (i = 0; i < vec_size; i++)
      {
        mpf_set_q(proj_mp->coord[i].r, (*proj_rat)[i][0]);
        proj_d->coord[i].r = mpq_get_d((*proj_rat)[i][0]);

        mpf_set_q(proj_mp->coord[i].i, (*proj_rat)[i][1]);
        proj_d->coord[i].i = mpq_get_d((*proj_rat)[i][1]);
      }
    }

    vec_size = v_d->size;
    MPI_Bcast(&vec_size, 1, MPI_INT, headnode, MPI_COMM_WORLD);
    bcast_vec_rat(v_rat, vec_size, my_id, headnode);

    if (my_id != headnode)
    { // setup _d & _mp
      init_vec_d(v_d, vec_size);
      init_vec_mp2(v_mp, vec_size, curr_prec);
      v_d->size = v_mp->size = vec_size;

      for (i = 0; i < vec_size; i++)
      {
        mpf_set_q(v_mp->coord[i].r, (*v_rat)[i][0]);
        v_d->coord[i].r = mpq_get_d((*v_rat)[i][0]);

        mpf_set_q(v_mp->coord[i].i, (*v_rat)[i][1]);
        v_d->coord[i].i = mpq_get_d((*v_rat)[i][1]);
      }
    }

    bcast_comp_num_rat(s_rat, 2, my_id, headnode);
    bcast_comp_rat(gamma_rat, my_id, headnode);

    if (my_id != headnode)
    { // setup _d & _mp
      init_mp2(s_mp[0], curr_prec);
      init_mp2(s_mp[1], curr_prec);
      init_mp2(gamma_mp, curr_prec);

      for (i = 0; i < 2; i++)
      {
        mpf_set_q(s_mp[i]->r, s_rat[i][0]);
        s_d[i]->r = mpq_get_d(s_rat[i][0]);

        mpf_set_q(s_mp[i]->i, s_rat[i][1]);
        s_d[i]->i = mpq_get_d(s_rat[i][1]);
      }

      mpf_set_q(gamma_mp->r, gamma_rat[0]);
      gamma_d->r = mpq_get_d(gamma_rat[0]);

      mpf_set_q(gamma_mp->i, gamma_rat[1]);
      gamma_d->i = mpq_get_d(gamma_rat[1]);
    }
  }

  return;
}

void bcast_monodromy_structures(int *continueMonodromy, vec_d v_out_d, vec_mp v_out_mp, mpq_t ***v_out_rat, vec_d v_in_d, vec_mp v_in_mp, mpq_t ***v_in_rat, comp_d gamma_out_d, comp_mp gamma_out_mp, mpq_t gamma_out_rat[2], comp_d gamma_in_d, comp_mp gamma_in_mp, mpq_t gamma_in_rat[2], int MPType, int curr_prec, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts the monodromy structures for this round     *
\***************************************************************/
{
  // put everyone on the same page as to whether to continue monodromy loops or to stop
  MPI_Bcast(continueMonodromy, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  if (*continueMonodromy) 
  { // send the structures
    if (MPType == 0)
    { // _d
      bcast_vec_d(v_out_d, my_id, headnode);
      bcast_vec_d(v_in_d, my_id, headnode);
      bcast_comp_d(gamma_out_d, my_id, headnode);
      bcast_comp_d(gamma_in_d, my_id, headnode);
    }
    else if (MPType == 1)
    { // _mp 
      bcast_vec_mp(v_out_mp, my_id, headnode);
      bcast_vec_mp(v_in_mp, my_id, headnode);
      bcast_comp_mp(gamma_out_mp, my_id, headnode);
      bcast_comp_mp(gamma_in_mp, my_id, headnode);
    }
    else
    { // _rat
      int i, vec_size = v_out_d->size;
      MPI_Bcast(&vec_size, 1, MPI_INT, headnode, MPI_COMM_WORLD);
      bcast_vec_rat(v_out_rat, vec_size, my_id, headnode);

      if (my_id != headnode)
      { // setup _d & _mp
        init_vec_d(v_out_d, vec_size);
        init_vec_mp2(v_out_mp, vec_size, curr_prec);
        v_out_d->size = v_out_mp->size = vec_size;

        for (i = 0; i < vec_size; i++)
        {
          mpf_set_q(v_out_mp->coord[i].r, (*v_out_rat)[i][0]);
          v_out_d->coord[i].r = mpq_get_d((*v_out_rat)[i][0]);

          mpf_set_q(v_out_mp->coord[i].i, (*v_out_rat)[i][1]);
          v_out_d->coord[i].i = mpq_get_d((*v_out_rat)[i][1]);
        }
      }

      vec_size = v_in_d->size;
      MPI_Bcast(&vec_size, 1, MPI_INT, headnode, MPI_COMM_WORLD);
      bcast_vec_rat(v_in_rat, vec_size, my_id, headnode);

      if (my_id != headnode)
      { // setup _d & _mp
        init_vec_d(v_in_d, vec_size);
        init_vec_mp2(v_in_mp, vec_size, curr_prec);
        v_in_d->size = v_in_mp->size = vec_size;
  
        for (i = 0; i < vec_size; i++)
        {
          mpf_set_q(v_in_mp->coord[i].r, (*v_in_rat)[i][0]);
          v_in_d->coord[i].r = mpq_get_d((*v_in_rat)[i][0]);

          mpf_set_q(v_in_mp->coord[i].i, (*v_in_rat)[i][1]);
          v_in_d->coord[i].i = mpq_get_d((*v_in_rat)[i][1]);
        }
      }

      bcast_comp_rat(gamma_out_rat, my_id, headnode);
      bcast_comp_rat(gamma_in_rat, my_id, headnode);

      if (my_id != headnode)
      { // setup _d & _mp
        init_mp2(gamma_out_mp, curr_prec);
        init_mp2(gamma_in_mp, curr_prec);

        mpf_set_q(gamma_out_mp->r, gamma_out_rat[0]);
        gamma_out_d->r = mpq_get_d(gamma_out_rat[0]);
 
        mpf_set_q(gamma_out_mp->i, gamma_out_rat[1]);
        gamma_out_d->i = mpq_get_d(gamma_out_rat[1]);

        mpf_set_q(gamma_in_mp->r, gamma_in_rat[0]);
        gamma_in_d->r = mpq_get_d(gamma_in_rat[0]);
 
        mpf_set_q(gamma_in_mp->i, gamma_in_rat[1]);
        gamma_in_d->i = mpq_get_d(gamma_in_rat[1]);
      }    
    }
  }

  return;
}

void bcast_regen_pos_dim_t(regen_pos_dim_t *RPD, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts RPD                                         *
\***************************************************************/
{
  MPI_Datatype mpi_RPD_int, mpi_comp_d;
  regen_pos_dim_t_int RPD_int;
  int size, *degrees = NULL, *ppd_type = NULL, *ppd_size = NULL, *prog_inst = NULL, *prog_gp_sizes = NULL;
  char *rpdStr = NULL, *progStr = NULL;
  comp_d *coeff_d = NULL;

  // create the datatypes mpi_RPD_int & mpi_comp_d
  create_regen_pos_dim_t_int(&mpi_RPD_int);
  create_comp_d(&mpi_comp_d);

  if (my_id == headnode)
  { // setup RPD_int & the other structures and then broadcast it
    cp_regen_pos_dim_int(&RPD_int, RPD, MPType, &rpdStr, &progStr, 0, &coeff_d, &degrees, &ppd_type, &ppd_size, &prog_inst, &prog_gp_sizes, 0);
    // send RPD_int
    MPI_Bcast(&RPD_int, 1, mpi_RPD_int, headnode, MPI_COMM_WORLD);
    // send rpdStr
    MPI_Bcast(rpdStr, RPD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send progStr
    MPI_Bcast(progStr, RPD_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send coeff_d
    MPI_Bcast(coeff_d, RPD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // send degrees
    MPI_Bcast(degrees, RPD_int.num_int, MPI_INT, headnode, MPI_COMM_WORLD);
    // send ppd_type & ppd_size
    size = RPD_int.PPD_int.num_hom_var_gp + RPD_int.PPD_int.num_var_gp;
    MPI_Bcast(ppd_type, size, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(ppd_size, size, MPI_INT, headnode, MPI_COMM_WORLD);
    // send prog_inst
    MPI_Bcast(prog_inst, RPD_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
    // send prog_gp_sizes
    MPI_Bcast(prog_gp_sizes, RPD_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);

    // clear memory
    free(degrees);
    free(ppd_type);
    free(ppd_size);
    free(prog_inst);
    free(prog_gp_sizes);
    free(rpdStr);
    free(progStr);
    free(coeff_d);
  }
  else
  { // recv RPD_int
    MPI_Bcast(&RPD_int, 1, mpi_RPD_int, headnode, MPI_COMM_WORLD);
    // allocate and recv other structures

    // recv rpdStr
    rpdStr = (char *)bmalloc(RPD_int.totalLength * sizeof(char));
    MPI_Bcast(rpdStr, RPD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv progStr
    progStr = (char *)bmalloc(RPD_int.Prog_int.totalLength * sizeof(char));
    MPI_Bcast(progStr, RPD_int.Prog_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv coeff_d
    coeff_d = (comp_d *)bmalloc(RPD_int.num_comp_d * sizeof(comp_d));
    MPI_Bcast(coeff_d, RPD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // recv degrees
    degrees = (int *)bmalloc(RPD_int.num_int * sizeof(int));
    MPI_Bcast(degrees, RPD_int.num_int, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv ppd_type & ppd_size
    size = RPD_int.PPD_int.num_hom_var_gp + RPD_int.PPD_int.num_var_gp;
    ppd_type = (int *)bmalloc(size * sizeof(int));
    ppd_size = (int *)bmalloc(size * sizeof(int));
    MPI_Bcast(ppd_type, size, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(ppd_size, size, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv prog_inst
    prog_inst = (int *)bmalloc(RPD_int.Prog_int.size * sizeof(int));
    MPI_Bcast(prog_inst, RPD_int.Prog_int.size, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv prog_gp_sizes
    prog_gp_sizes = (int *)bmalloc(RPD_int.Prog_int.num_var_gps * sizeof(int));
    MPI_Bcast(prog_gp_sizes, RPD_int.Prog_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);

    // setup RPD - all data is freed in this function call!
    cp_regen_pos_dim_int(RPD, &RPD_int, MPType, &rpdStr, &progStr, 1, &coeff_d, &degrees, &ppd_type, &ppd_size, &prog_inst, &prog_gp_sizes, 1);
  }

  // free mpi_RPD_int & mpi_comp_d
  MPI_Type_free(&mpi_RPD_int);
  MPI_Type_free(&mpi_comp_d);

  return;
}

void bcast_regenCodim_t(regenCodim_t *RCD, int curr_prec, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts CD                                          *
\***************************************************************/
{
  MPI_Datatype mpi_RCD_int, mpi_comp_d;
  regenCodim_t_int RCD_int;
  char *codimStr = NULL;
  comp_d *coeff_d = NULL;

  // create the datatypes mpi_RCD_int & mpi_comp_d
  create_regenCodim_t_int(&mpi_RCD_int);
  create_comp_d(&mpi_comp_d);

  if (my_id == headnode)
  { // setup RCD_int & the other structures and then broadcast it
    cp_regenCodim_int(&RCD_int, RCD, MPType, &codimStr, &coeff_d, curr_prec, 0, 0);
    // send RCD_int
    MPI_Bcast(&RCD_int, 1, mpi_RCD_int, headnode, MPI_COMM_WORLD);
    // send codimStr
    MPI_Bcast(codimStr, RCD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // send coeff_d
    MPI_Bcast(coeff_d, RCD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // clear memory
    free(codimStr);
    free(coeff_d);
  }
  else
  { // recv RCD_int
    MPI_Bcast(&RCD_int, 1, mpi_RCD_int, headnode, MPI_COMM_WORLD);
    // allocate and recv other structures

    // recv codimStr
    codimStr = (char *)bmalloc(RCD_int.totalLength * sizeof(char));
    MPI_Bcast(codimStr, RCD_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv coeff_d
    coeff_d = (comp_d *)bmalloc(RCD_int.num_comp_d * sizeof(comp_d));
    MPI_Bcast(coeff_d, RCD_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // setup RCD - all data is freed in this function call!
    cp_regenCodim_int(RCD, &RCD_int, MPType, &codimStr, &coeff_d, curr_prec, 1, 1);
  }

  // free mpi_RCD_int & mpi_comp_d
  MPI_Type_free(&mpi_RCD_int);
  MPI_Type_free(&mpi_comp_d);

  return;
}

void send_recv_point_data_d(point_data_d **PD, int *numPts, int targetNum, int isSending)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: numPts: > 0 - number of points, otherwise it       *
*   represents some code to the workers                         *
*       isSending: 1 - send, otherwise recv                     *
* RETURN VALUES:                                                *
* NOTES: sends/recvs points to/from targetNum                   *
\***************************************************************/
{
  int i;
  MPI_Datatype mpi_point_data_d_int, mpi_comp_d;
  point_data_d_int PD_int;
  comp_d *coeff = NULL;

  // create the datatype mpi_point_data_d_int & mpi_comp_d
  create_point_data_d_int(&mpi_point_data_d_int);
  create_comp_d(&mpi_comp_d);

  if (isSending)
  { // send numPts & PD to targetNum
    MPI_Send(numPts, 1, MPI_INT, targetNum, TAG_NUM_PTS, MPI_COMM_WORLD);
    // setup PD_int
    for (i = 0; i < *numPts; i++)
    { // setup PD_int
      cp_point_data_d_int(&PD_int, &(*PD)[i], &coeff, 0, 0);

      // send PD_int
      MPI_Send(&PD_int, 1, mpi_point_data_d_int, targetNum, TAG_POINT_DATA_D, MPI_COMM_WORLD);
      // send coeff
      MPI_Send(&PD_int, PD_int.point_int.size, mpi_comp_d, targetNum, TAG_POINT_DATA_D, MPI_COMM_WORLD);
    
      // clear coeff
      free(coeff);
    }
  }
  else 
  { // recv from targetNum
    int tempInt;
    MPI_Status status;

    // recv numPts
    MPI_Recv(&tempInt, 1, MPI_INT, targetNum, TAG_NUM_PTS, MPI_COMM_WORLD, &status);

    if (tempInt > 0)
    { // make sure that PD is the correct size
      if (*numPts != tempInt)
      { // allocate to the correct size
        *PD = (point_data_d *)brealloc(*PD, tempInt * sizeof(point_data_d));
      }
    }
    else
    { // clear PD
      free(*PD);
      *PD = NULL;
    }

    // setup numPts
    *numPts = tempInt;

    // recv PD_int
    for (i = 0; i < tempInt; i++)
    { // recv PD_int
      MPI_Recv(&PD_int, 1, mpi_point_data_d_int, targetNum, TAG_POINT_DATA_D, MPI_COMM_WORLD, &status);
      // setup & recv coeff
      coeff = (comp_d *)bmalloc(PD_int.point_int.size * sizeof(comp_d));    
      MPI_Recv(coeff, PD_int.point_int.size, mpi_comp_d, targetNum, TAG_POINT_DATA_D, MPI_COMM_WORLD, &status);

      // setup PD[i]
      cp_point_data_d_int(&(*PD)[i], &PD_int, &coeff, 1, 1);
    }
  }

  // free mpi_point_data_d_int & mpi_comp_d
  MPI_Type_free(&mpi_point_data_d_int);
  MPI_Type_free(&mpi_comp_d);

  return;
}

int send_recv_endgame_data_t(endgame_data_t **EG, int *numPts, int MPType, int targetNum, int isSending)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: numPts: > 0 - number of points, otherwise it       *
*   represents some code to the workers                         *
*       isSending: 1 - send, otherwise recv                     *
* RETURN VALUES: id number of who it sent/recv to/from          *
* NOTES: sends/recvs endgame_data to/from targetNum             *
\***************************************************************/
{
  int i, j, retVal, strLoc, coeffLoc, totalSize, totalCoeff;
  MPI_Datatype mpi_eg_t, mpi_comp_d;
  endgame_data_t_int *EG_int = NULL; 
  int *coeff_size = NULL;
  char *egStr = NULL, *egStr_temp = NULL;
  comp_d *coeff = NULL, *coeff_temp = NULL;

  // create the datatype mpi_eg_t & mpi_comp_d
  create_endgame_data_t_int(&mpi_eg_t);
  create_comp_d(&mpi_comp_d);

  if (isSending)
  { // send numPts & EG to targetNum
    MPI_Send(numPts, 1, MPI_INT, targetNum, TAG_NUM_PTS, MPI_COMM_WORLD);

    if (*numPts > 0)
    { // setup to send EG
      EG_int = (endgame_data_t_int *)bmalloc((*numPts) * sizeof(endgame_data_t_int));
      coeff_size = (int *)bmalloc((*numPts) * sizeof(int));

      coeffLoc = strLoc = totalSize = totalCoeff = 0;
      for (i = 0; i < *numPts; i++)
      { // setup each of the points
        cp_endgame_data_t_int(&EG_int[i], &(*EG)[i], &egStr_temp, &coeff_temp, 0, 0, 0);
        // add on to the totalSize
        totalSize += EG_int[i].totalLength;
        // add on to the totalCoeff
        coeff_size[i] = 0;
        if (EG_int[i].prec == 52)
        {
          coeff_size[i] = EG_int[i].PD_d_int.point_int.size;
          if (EG_int[i].last_approx_prec == 52)
            coeff_size[i] += EG_int[i].last_approx_d_int.size;
        }
        else if (EG_int[i].last_approx_prec == 52)
          coeff_size[i] = EG_int[i].last_approx_d_int.size;
        totalCoeff += coeff_size[i];

        // setup egStr
        egStr = (char *)brealloc(egStr, totalSize * sizeof(char));
        for (j = 0; j < EG_int[i].totalLength; j++)
        {
          egStr[strLoc] = egStr_temp[j];
          strLoc++;
        }
        free(egStr_temp);

        // setup coeff
        coeff = (comp_d *)brealloc(coeff, totalCoeff * sizeof(comp_d));
        for (j = 0; j < coeff_size[i]; j++)
        {
          set_d(coeff[coeffLoc], coeff_temp[j]);
          coeffLoc++;
        }
        free(coeff_temp);
      }

      // send EG_int
      MPI_Send(EG_int, *numPts, mpi_eg_t, targetNum, TAG_ENDGAME_DATA_T, MPI_COMM_WORLD);
      // send egStr
      MPI_Send(egStr, totalSize, MPI_CHAR, targetNum, TAG_ENDGAME_STR, MPI_COMM_WORLD);
      // send coeff
      MPI_Send(coeff, totalCoeff, mpi_comp_d, targetNum, TAG_ENDGAME_COEFF, MPI_COMM_WORLD);

      // free memory
      free(coeff_size);
      free(coeff);
      free(egStr);
      free(EG_int);
    }

    retVal = targetNum; // who it sent to
  }
  else
  { // recv from targetNum
    int tempInt;
    MPI_Status status;

    // recv numPts 
    MPI_Recv(&tempInt, 1, MPI_INT, targetNum, TAG_NUM_PTS, MPI_COMM_WORLD, &status);

    // determine who sent this packet
    retVal = status.MPI_SOURCE;

    if (tempInt > 0)
    { // setup to recv EG_int
      EG_int = (endgame_data_t_int *)bmalloc(tempInt * sizeof(endgame_data_t_int));
      // recv EG_int
      MPI_Recv(EG_int, tempInt, mpi_eg_t, retVal, TAG_ENDGAME_DATA_T, MPI_COMM_WORLD, &status);

      // setup to recv egStr & coeff
      coeff_size = (int *)bmalloc(tempInt * sizeof(int));
      totalSize = totalCoeff = 0;
      for (i = 0; i < tempInt; i++)
      {
        totalSize += EG_int[i].totalLength;
        coeff_size[i] = 0;
        if (EG_int[i].prec == 52)
        {
          coeff_size[i] = EG_int[i].PD_d_int.point_int.size;
          if (EG_int[i].last_approx_prec == 52)
            coeff_size[i] += EG_int[i].last_approx_d_int.size;
        }
        else if (EG_int[i].last_approx_prec == 52)
          coeff_size[i] = EG_int[i].last_approx_d_int.size;

        totalCoeff += coeff_size[i];
      }
      egStr = (char *)bmalloc(totalSize * sizeof(char));
      coeff = (comp_d *)bmalloc(totalCoeff * sizeof(comp_d));

      // recv egStr from the same source
      MPI_Recv(egStr, totalSize, MPI_CHAR, retVal, TAG_ENDGAME_STR, MPI_COMM_WORLD, &status);
      // recv coeff from the same source
      MPI_Recv(coeff, totalCoeff, mpi_comp_d, retVal, TAG_ENDGAME_COEFF, MPI_COMM_WORLD, &status);

      // make sure the EG is the correct size
      if (tempInt != *numPts)
      { // clear the currently allocated memory if needed
        for (i = *numPts - 1; i >= 0; i--)
          clear_endgame_data(&(*EG)[i]);

        // allocate to the correct size
        *EG = (endgame_data_t *)brealloc(*EG, tempInt * sizeof(endgame_data_t));

        // initialize EG - when doing the copying, they automatically get setup - so only need the ones that are using AMP with double precision
        for (i = 0; i < tempInt; i++)
          init_endgame_data(&(*EG)[i], 64); 
      }

      // setup EG
      strLoc = coeffLoc = 0;
      for (i = 0; i < tempInt; i++)
      { // setup egStr_temp
        egStr_temp = &egStr[strLoc];
        // setup coeff_temp
        coeff_temp = &coeff[coeffLoc];

        // setup each of the points
        cp_endgame_data_t_int(&(*EG)[i], &EG_int[i], &egStr_temp, &coeff_temp, 0, 0, 1);
        // add on to strLoc 
        strLoc += EG_int[i].totalLength;
        // add on to coeffLoc
        coeffLoc += coeff_size[i];
      }

      // free coeff, egStr, coeff_size, EG_int
      free(coeff);
      free(egStr);
      free(coeff_size);
      free(EG_int);
    }
    else
    { // clear the currently allocated memory if needed
      for (i = *numPts - 1; i >= 0; i--)
        clear_endgame_data(&(*EG)[i]);

      // clear EG
      free(*EG);
      *EG = NULL;
    }  
    // setup numPts
    *numPts = tempInt;
  }

  // free mpi_eg_t & mpi_comp_d
  MPI_Type_free(&mpi_eg_t);
  MPI_Type_free(&mpi_comp_d);

  return retVal;
}

int send_recv_corank_data(int *corank, double *sm, double *lg, int numPts, int targetNum, int isSending)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: isSending: 1 - send, otherwise recv                *
* RETURN VALUES: id number of who it sent/recv to/from          *
* NOTES: sends/recvs approximation to/from targetNum            *
\***************************************************************/ 
{
  int retVal = 0;

  if (isSending)
  { // send information
    if (numPts > 0)
    { // send info
      MPI_Send(corank, numPts, MPI_INT, targetNum, TAG_OTHER_DATA, MPI_COMM_WORLD);
      MPI_Send(sm, numPts, MPI_DOUBLE, targetNum, TAG_OTHER_DATA, MPI_COMM_WORLD);
      MPI_Send(lg, numPts, MPI_DOUBLE, targetNum, TAG_OTHER_DATA, MPI_COMM_WORLD);
    }

    retVal = targetNum;
  }
  else
  { // recv information
    retVal = targetNum;

    if (numPts > 0)
    { // recv other info
      MPI_Status status;

      MPI_Recv(corank, numPts, MPI_INT, retVal, TAG_OTHER_DATA, MPI_COMM_WORLD, &status);

      // determine who sent this packet
      retVal = status.MPI_SOURCE;

      MPI_Recv(sm, numPts, MPI_DOUBLE, retVal, TAG_OTHER_DATA, MPI_COMM_WORLD, &status);
      MPI_Recv(lg, numPts, MPI_DOUBLE, retVal, TAG_OTHER_DATA, MPI_COMM_WORLD, &status);
    }
  }

  return retVal;
}

int send_recv_last_approximations(int *prec, point_d *pt_d, point_mp *pt_mp, int *corank, double *sm, double *lg, int numPts, int MPType, int targetNum, int isSending)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: isSending: 1 - send, otherwise recv                *
* RETURN VALUES: id number of who it sent/recv to/from          *
* NOTES: sends/recvs approximation to/from targetNum            *
\***************************************************************/ 
{
  int i, j, retVal, strLoc, coeffLoc, totalSize, totalCoeff;
  char *str = NULL, *str_temp = NULL;
  comp_d *coeff = NULL, *coeff_temp = NULL;
  point_d_int *pt_d_int = NULL;
  point_mp_int *pt_mp_int = NULL;

  if (isSending)
  { // setup to send information
    if (numPts > 0)
    { // setup to send the points
      if (MPType == 0)
      { // all in double
        MPI_Datatype mpi_point_d_int, mpi_comp_d;
        create_point_d_int(&mpi_point_d_int);
        create_comp_d(&mpi_comp_d);

        pt_d_int = (point_d_int *)bmalloc(numPts * sizeof(point_d_int));

        coeffLoc = totalCoeff = 0;
        for (i = 0; i < numPts; i++)
        { // setup each of the points
          cp_point_d_int(&pt_d_int[i], &pt_d[i], &coeff_temp, 0, 0, 0);
          totalCoeff += pt_d_int[i].size;
          coeff = (comp_d *)brealloc(coeff, totalCoeff * sizeof(comp_d));
          for (j = 0; j < pt_d_int[i].size; j++)
          {
            set_d(coeff[coeffLoc], coeff_temp[j]);
            coeffLoc++;
          }
          free(coeff_temp);
        }

        // send pt_d_int
        MPI_Send(pt_d_int, numPts, mpi_point_d_int, targetNum, TAG_POINT_D_INT, MPI_COMM_WORLD);
        // send coeff
        MPI_Send(coeff, totalCoeff, mpi_comp_d, targetNum, TAG_COEFF, MPI_COMM_WORLD);

        // clear memory
        free(coeff);
        free(pt_d_int);

        MPI_Type_free(&mpi_point_d_int);
        MPI_Type_free(&mpi_comp_d);
      }
      else if (MPType == 1)
      { // all in fixed MP
        MPI_Datatype mpi_point_mp_int;
        create_point_mp_int(&mpi_point_mp_int);

        pt_mp_int = (point_mp_int *)bmalloc(numPts * sizeof(point_mp_int));
        
        strLoc = totalSize = 0;
        for (i = 0; i < numPts; i++)
        { // setup each of the points
          cp_point_mp_int(&pt_mp_int[i], &pt_mp[i], &str_temp, 0, 0, 0);
          totalSize += pt_mp_int[i].totalLength;
          str = (char *)brealloc(str, totalSize * sizeof(char));
          for (j = 0; j < pt_mp_int[i].totalLength; j++)
          {
            str[strLoc] = str_temp[j];
            strLoc++;
          }
          free(str_temp);
        }

        // send pt_mp_int
        MPI_Send(pt_mp_int, numPts, mpi_point_mp_int, targetNum, TAG_POINT_MP_INT, MPI_COMM_WORLD);
        // send str
        MPI_Send(str, totalSize, MPI_CHAR, targetNum, TAG_STR, MPI_COMM_WORLD);

        // clear memory
        free(str);
        free(pt_mp_int);

        MPI_Type_free(&mpi_point_mp_int);
      }
      else
      { // using AMP
        MPI_Datatype mpi_point_d_int, mpi_point_mp_int, mpi_comp_d;
        create_point_d_int(&mpi_point_d_int);
        create_point_mp_int(&mpi_point_mp_int);
        create_comp_d(&mpi_comp_d);

        pt_d_int = (point_d_int *)bmalloc(numPts * sizeof(point_d_int));
        pt_mp_int = (point_mp_int *)bmalloc(numPts * sizeof(point_mp_int));
     
        coeffLoc = strLoc = totalSize = totalCoeff = 0;
        for (i = 0; i < numPts; i++)
        { // setup the ith point
          if (prec[i] < 64)
          { // use pt_d
            cp_point_d_int(&pt_d_int[i], &pt_d[i], &coeff_temp, 0, 0, 0);
            totalCoeff += pt_d_int[i].size;
            coeff = (comp_d *)brealloc(coeff, totalCoeff * sizeof(comp_d));
            for (j = 0; j < pt_d_int[i].size; j++)
            {
              set_d(coeff[coeffLoc], coeff_temp[j]);
              coeffLoc++;
            }
            free(coeff_temp);
          }
          else
          { // use pt_mp
            cp_point_mp_int(&pt_mp_int[i], &pt_mp[i], &str_temp, 0, 0, 0);
            totalSize += pt_mp_int[i].totalLength;
            str = (char *)brealloc(str, totalSize * sizeof(char));
            for (j = 0; j < pt_mp_int[i].totalLength; j++)
            {
              str[strLoc] = str_temp[j];
              strLoc++;
            }
            free(str_temp);
          }
        }
   
        // send prec
        MPI_Send(prec, numPts, MPI_INT, targetNum, TAG_OTHER_DATA, MPI_COMM_WORLD);
        // send pt_d_int
        MPI_Send(pt_d_int, numPts, mpi_point_d_int, targetNum, TAG_POINT_D_INT, MPI_COMM_WORLD);
        // send pt_mp_int
        MPI_Send(pt_mp_int, numPts, mpi_point_mp_int, targetNum, TAG_POINT_MP_INT, MPI_COMM_WORLD);
        // send coeff
        MPI_Send(coeff, totalCoeff, mpi_comp_d, targetNum, TAG_COEFF, MPI_COMM_WORLD);
        // send str
        MPI_Send(str, totalSize, MPI_CHAR, targetNum, TAG_STR, MPI_COMM_WORLD);

        // clear memory
        free(coeff);
        free(str);
        free(pt_d_int);
        free(pt_mp_int);

        MPI_Type_free(&mpi_point_d_int);
        MPI_Type_free(&mpi_point_mp_int);
        MPI_Type_free(&mpi_comp_d);
      }

      // send other info
      MPI_Send(corank, numPts, MPI_INT, targetNum, TAG_OTHER_DATA, MPI_COMM_WORLD);
      MPI_Send(sm, numPts, MPI_DOUBLE, targetNum, TAG_OTHER_DATA, MPI_COMM_WORLD);
      MPI_Send(lg, numPts, MPI_DOUBLE, targetNum, TAG_OTHER_DATA, MPI_COMM_WORLD);
    }

    retVal = targetNum; // who it sent to
  }
  else
  { // recv from targetNum

    // initialize retVal
    retVal = targetNum;

    if (numPts > 0)
    { // recv the points
      int tempInt;
      MPI_Status status;

      if (MPType == 0)
      { // all in double
        MPI_Datatype mpi_point_d_int, mpi_comp_d;
        create_point_d_int(&mpi_point_d_int);
        create_comp_d(&mpi_comp_d);

        // setup prec
        for (i = 0; i < numPts; i++)
          prec[i] = 52;

        // setup and recv pt_d_int
        pt_d_int = (point_d_int *)bmalloc(numPts * sizeof(point_d_int));
        MPI_Recv(pt_d_int, numPts, mpi_point_d_int, targetNum, TAG_POINT_D_INT, MPI_COMM_WORLD, &status);

        // determine who sent this packet
        retVal = status.MPI_SOURCE;

        // setup and recv coeff
        totalCoeff = 0;
        for (i = 0; i < numPts; i++)
          totalCoeff += pt_d_int[i].size;
        coeff = (comp_d *)bmalloc(totalCoeff * sizeof(comp_d));
        MPI_Recv(coeff, totalCoeff, mpi_comp_d, retVal, TAG_COEFF, MPI_COMM_WORLD, &status);

        // setup point_d
        coeffLoc = 0;
        for (i = 0; i < numPts; i++)
        { // setup coeff_temp
          coeff_temp = &coeff[coeffLoc];

          // setup each of the points
          cp_point_d_int(&pt_d[i], &pt_d_int[i], &coeff_temp, 0, 0, 1);

          // add on to coeffLoc
          coeffLoc += pt_d_int[i].size;
        }

        // clear memory
        free(pt_d_int);
        free(coeff);

        MPI_Type_free(&mpi_point_d_int);
        MPI_Type_free(&mpi_comp_d);
      }
      else if (MPType == 1)
      { // all in fixed MP
        MPI_Datatype mpi_point_mp_int;
        create_point_mp_int(&mpi_point_mp_int);

        // setup prec
        tempInt = mpf_get_default_prec();
        for (i = 0; i < numPts; i++)
          prec[i] = tempInt;

        // setup and recv pt_mp_int
        pt_mp_int = (point_mp_int *)bmalloc(numPts * sizeof(point_mp_int));
        MPI_Recv(pt_mp_int, numPts, mpi_point_mp_int, targetNum, TAG_POINT_MP_INT, MPI_COMM_WORLD, &status);

        // determine who sent this packet
        retVal = status.MPI_SOURCE;

        // setup and recv str
        totalSize = 0;
        for (i = 0; i < numPts; i++)
          totalSize += pt_mp_int[i].totalLength;
        str = (char *)bmalloc(totalSize * sizeof(char));
        MPI_Recv(str, totalSize, MPI_CHAR, retVal, TAG_STR, MPI_COMM_WORLD, &status);

        // setup point_mp
        strLoc = 0;
        for (i = 0; i < numPts; i++)
        { // setup str_temp
          str_temp = &str[strLoc];

          // setup each of the points
          cp_point_mp_int(&pt_mp[i], &pt_mp_int[i], &str_temp, 0, 0, 1);

          // add on to strLoc
          strLoc += pt_mp_int[i].totalLength;
        }

        // clear memory
        free(pt_mp_int);
        free(str);

        MPI_Type_free(&mpi_point_mp_int);
      }
      else
      { // using AMP
        MPI_Datatype mpi_point_d_int, mpi_point_mp_int, mpi_comp_d;
        create_point_d_int(&mpi_point_d_int);
        create_point_mp_int(&mpi_point_mp_int);
        create_comp_d(&mpi_comp_d);

        // recv prec
        MPI_Recv(prec, numPts, MPI_INT, targetNum, TAG_OTHER_DATA, MPI_COMM_WORLD, &status);

        // determine who sent this
        retVal = status.MPI_SOURCE;

        // setup and recv pt_d_int & pt_mp_int
        pt_d_int = (point_d_int *)bmalloc(numPts * sizeof(point_d_int));
        pt_mp_int = (point_mp_int *)bmalloc(numPts * sizeof(point_mp_int));
        MPI_Recv(pt_d_int, numPts, mpi_point_d_int, retVal, TAG_POINT_D_INT, MPI_COMM_WORLD, &status);
        MPI_Recv(pt_mp_int, numPts, mpi_point_mp_int, retVal, TAG_POINT_MP_INT, MPI_COMM_WORLD, &status);

        // setup and recv coeff & str
        totalSize = totalCoeff = 0;
        for (i = 0; i < numPts; i++)
          if (prec[i] < 64)
            totalCoeff += pt_d_int[i].size;
          else
            totalSize += pt_mp_int[i].totalLength;
        coeff = (comp_d *)bmalloc(totalCoeff * sizeof(comp_d));
        str = (char *)bmalloc(totalSize * sizeof(char));
        MPI_Recv(coeff, totalCoeff, mpi_comp_d, retVal, TAG_COEFF, MPI_COMM_WORLD, &status);
        MPI_Recv(str, totalSize, MPI_CHAR, retVal, TAG_STR, MPI_COMM_WORLD, &status);

        // setup points
        strLoc = coeffLoc = 0;
        for (i = 0; i < numPts; i++)
        { // setup str_temp & coeff_temp
          str_temp = &str[strLoc];
          coeff_temp = &coeff[coeffLoc];

          if (prec[i] < 64)
          { // setup pt_d
            cp_point_d_int(&pt_d[i], &pt_d_int[i], &coeff_temp, 0, 0, 1);
            // add on to coeffLoc
            coeffLoc += pt_d_int[i].size;
          }
          else
          { // setup pt_mp
            cp_point_mp_int(&pt_mp[i], &pt_mp_int[i], &str_temp, 0, 0, 1);
            // add on to strLoc
            strLoc += pt_mp_int[i].totalLength;
          }
        }
        // clear memory
        free(pt_d_int);
        free(pt_mp_int);
        free(coeff);
        free(str);

        MPI_Type_free(&mpi_point_d_int);
        MPI_Type_free(&mpi_point_mp_int);
        MPI_Type_free(&mpi_comp_d);
      }

      // recv other info
      MPI_Recv(corank, numPts, MPI_INT, retVal, TAG_OTHER_DATA, MPI_COMM_WORLD, &status);
      MPI_Recv(sm, numPts, MPI_DOUBLE,  retVal, TAG_OTHER_DATA, MPI_COMM_WORLD, &status);
      MPI_Recv(lg, numPts, MPI_DOUBLE, retVal, TAG_OTHER_DATA, MPI_COMM_WORLD, &status);
    }
  }

  return retVal;
}

int send_recv_trackBack_samples_t(trackBack_samples_t **TB, int *numPts, int MPType, int targetNum, int isSending)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: numPts: > 0 - number of points, otherwise it       *
*   represents some code to the workers                         *
*       isSending: 1 - send, otherwise recv                     *
* RETURN VALUES: id number of who it sent/recv to/from          *
* NOTES: sends/recvs trackBack_samples to/from targetNum        *
\***************************************************************/
{
  int i, j, retVal, egStrLoc = 0, tbStrLoc = 0, egCoeffLoc = 0, tbCoeffLoc = 0, tbDoubleLoc = 0, tbIntLoc = 0;
  int egStrSize = 0, tbStrSize = 0, tbDoubleCount = 0, egCoeffCount = 0, tbCoeffCount = 0, tbIntCount = 0;
  int *egCount = NULL, *tbInt = NULL, **tbInt_temp = NULL;
  char *egStr = NULL, **egStr_temp = NULL, *tbStr = NULL, **tbStr_temp = NULL;
  double *tbDouble = NULL, **tbDouble_temp = NULL;
  comp_d *egCoeff = NULL, **egCoeff_temp = NULL, *tbCoeff = NULL, **tbCoeff_temp = NULL;
  trackBack_samples_t_int *TB_int = NULL;
  MPI_Datatype mpi_comp_d, mpi_tb_t;

  // create the datatypes
  create_comp_d(&mpi_comp_d);
  create_trackBack_samples_t_int(&mpi_tb_t);

  if (isSending)
  { // send to targetNum
    MPI_Send(numPts, 1, MPI_INT, targetNum, TAG_NUM_PTS, MPI_COMM_WORLD);

    if (*numPts > 0)
    { // setup to send TB
      TB_int = (trackBack_samples_t_int *)bmalloc(*numPts * sizeof(trackBack_samples_t_int));
      egStr_temp = (char **)bmalloc(*numPts * sizeof(char *));
      tbStr_temp = (char **)bmalloc(*numPts * sizeof(char *));
      tbDouble_temp = (double **)bmalloc(*numPts * sizeof(double *));
      egCoeff_temp = (comp_d **)bmalloc(*numPts * sizeof(comp_d *));
      tbCoeff_temp = (comp_d **)bmalloc(*numPts * sizeof(comp_d *));
      tbInt_temp = (int **)bmalloc(*numPts * sizeof(int *));
      egCount = (int *)bmalloc(*numPts * sizeof(int));

      // setup each of the points
      for (i = 0; i < *numPts; i++)
      {
        cp_trackBack_samples_t_int(&TB_int[i], &(*TB)[i], &egStr_temp[i], &egCoeff_temp[i], &tbDouble_temp[i], &tbCoeff_temp[i], &tbStr_temp[i], &tbInt_temp[i], 0, 0, 0);

        // increment the sizes
        egStrSize += TB_int[i].EG_int.totalLength;
        tbStrSize += TB_int[i].totalLength;
        tbDoubleCount += TB_int[i].num_double;

        egCount[i] = 0;
        if (TB_int[i].EG_int.prec == 52)
        {
          egCount[i] = TB_int[i].EG_int.PD_d_int.point_int.size;
          if (TB_int[i].EG_int.last_approx_prec == 52)
            egCount[i] += TB_int[i].EG_int.last_approx_d_int.size;
        }
        else if (TB_int[i].EG_int.last_approx_prec == 52)
          egCount[i] = TB_int[i].EG_int.last_approx_d_int.size;

        egCoeffCount += egCount[i];
        tbCoeffCount += TB_int[i].num_comp_d;
        tbIntCount += TB_int[i].numSamples;
      }

      // setup egStr
      egStr = (char *)bmalloc(egStrSize * sizeof(char));
      for (i = 0; i < *numPts; i++)
        for (j = 0; j < TB_int[i].EG_int.totalLength; j++)
        {
          egStr[egStrLoc] = egStr_temp[i][j];
          egStrLoc++;
        }

      // setup tbStr
      tbStr = (char *)bmalloc(tbStrSize * sizeof(char));
      for (i = 0; i < *numPts; i++)
        for (j = 0; j < TB_int[i].totalLength; j++)
        {
          tbStr[tbStrLoc] = tbStr_temp[i][j];
          tbStrLoc++;
        }

      // setup tbDouble
      tbDouble = (double *)bmalloc(tbDoubleCount * sizeof(double));
      for (i = 0; i < *numPts; i++)
        for (j = 0; j < TB_int[i].num_double; j++)
        {
          tbDouble[tbDoubleLoc] = tbDouble_temp[i][j];
          tbDoubleLoc++;
        }

      // setup egCoeff
      egCoeff = (comp_d *)bmalloc(egCoeffCount * sizeof(comp_d));
      for (i = 0; i < *numPts; i++)
        for (j = 0; j < egCount[i]; j++)
        {
          set_d(egCoeff[egCoeffLoc], egCoeff_temp[i][j]);
          egCoeffLoc++;
        }

      // setup tbCoeff
      tbCoeff = (comp_d *)bmalloc(tbCoeffCount * sizeof(comp_d));
      for (i = 0; i < *numPts; i++)
        for (j = 0; j < TB_int[i].num_comp_d; j++)
        {
          set_d(tbCoeff[tbCoeffLoc], tbCoeff_temp[i][j]);
          tbCoeffLoc++;
        }

      // setup tbInt
      tbInt = (int *)bmalloc(tbIntCount * sizeof(int));
      for (i = 0; i < *numPts; i++)
        for (j = 0; j < TB_int[i].numSamples; j++)
        {
          tbInt[tbIntLoc] = tbInt_temp[i][j];
          tbIntLoc++;
        }

      // send the information
      MPI_Send(TB_int, *numPts, mpi_tb_t, targetNum, TAG_TRACKBACK_T, MPI_COMM_WORLD);
      MPI_Send(egStr, egStrSize, MPI_CHAR, targetNum, TAG_ENDGAME_STR, MPI_COMM_WORLD);
      MPI_Send(tbStr, tbStrSize, MPI_CHAR, targetNum, TAG_STR, MPI_COMM_WORLD);
      MPI_Send(tbDouble, tbDoubleCount, MPI_DOUBLE, targetNum, TAG_DOUBLE, MPI_COMM_WORLD);
      MPI_Send(egCoeff, egCoeffCount, mpi_comp_d, targetNum, TAG_ENDGAME_COEFF, MPI_COMM_WORLD);
      MPI_Send(tbCoeff, tbCoeffCount, mpi_comp_d, targetNum, TAG_COEFF, MPI_COMM_WORLD);
      MPI_Send(tbInt, tbIntCount, MPI_INT, targetNum, TAG_INT, MPI_COMM_WORLD);

      // clear memory
      for (i = *numPts - 1; i >= 0; i--)
      {
        free(egStr_temp[i]);
        free(tbStr_temp[i]);
        free(tbDouble_temp[i]);
        free(egCoeff_temp[i]);
        free(tbCoeff_temp[i]);
        free(tbInt_temp[i]);
      }
      free(egStr_temp);
      free(tbStr_temp);
      free(tbDouble_temp);
      free(egCoeff_temp);
      free(tbCoeff_temp);
      free(tbInt_temp);
      free(egStr);
      free(tbStr);
      free(tbDouble);
      free(egCoeff);
      free(tbCoeff);
      free(tbInt);
      free(TB_int);
    }

    retVal = targetNum; // who it sent to
  }
  else
  { // recv from targetNum
    int tempInt;
    MPI_Status status;

    // recv numPts
    MPI_Recv(&tempInt, 1, MPI_INT, targetNum, TAG_NUM_PTS, MPI_COMM_WORLD, &status);
  
    // determine who sent this packet
    retVal = status.MPI_SOURCE;

    if (tempInt > 0)
    { // setup to recv TB_int
      TB_int = (trackBack_samples_t_int *)bmalloc(tempInt * sizeof(trackBack_samples_t_int));
      // recv TB_int
      MPI_Recv(TB_int, tempInt, mpi_tb_t, retVal, TAG_TRACKBACK_T, MPI_COMM_WORLD, &status);

      // setup to recv other structures
      egCount = (int *)bmalloc(tempInt * sizeof(int));
      for (i = 0; i < tempInt; i++)
      {
        egStrSize += TB_int[i].EG_int.totalLength;
        tbStrSize += TB_int[i].totalLength;
        tbDoubleCount += TB_int[i].num_double;

        egCount[i] = 0;
        if (TB_int[i].EG_int.prec == 52)
        {
          egCount[i] = TB_int[i].EG_int.PD_d_int.point_int.size;
          if (TB_int[i].EG_int.last_approx_prec == 52)
            egCount[i] += TB_int[i].EG_int.last_approx_d_int.size;
        }
        else if (TB_int[i].EG_int.last_approx_prec == 52)
          egCount[i] = TB_int[i].EG_int.last_approx_d_int.size;

        egCoeffCount += egCount[i];
        tbCoeffCount += TB_int[i].num_comp_d;
        tbIntCount += TB_int[i].numSamples;
      }
      // allocate structures
      egStr = (char *)bmalloc(egStrSize * sizeof(char));
      tbStr = (char *)bmalloc(tbStrSize * sizeof(char));
      tbDouble = (double *)bmalloc(tbDoubleCount * sizeof(double));
      egCoeff = (comp_d *)bmalloc(egCoeffCount * sizeof(comp_d));
      tbCoeff = (comp_d *)bmalloc(tbCoeffCount * sizeof(comp_d));
      tbInt = (int *)bmalloc(tbIntCount * sizeof(int));

      // recv structures
      MPI_Recv(egStr, egStrSize, MPI_CHAR, retVal, TAG_ENDGAME_STR, MPI_COMM_WORLD, &status);
      MPI_Recv(tbStr, tbStrSize, MPI_CHAR, retVal, TAG_STR, MPI_COMM_WORLD, &status);
      MPI_Recv(tbDouble, tbDoubleCount, MPI_DOUBLE, retVal, TAG_DOUBLE, MPI_COMM_WORLD, &status);
      MPI_Recv(egCoeff, egCoeffCount, mpi_comp_d, retVal, TAG_ENDGAME_COEFF, MPI_COMM_WORLD, &status);
      MPI_Recv(tbCoeff, tbCoeffCount, mpi_comp_d, retVal, TAG_COEFF, MPI_COMM_WORLD, &status);
      MPI_Recv(tbInt, tbIntCount, MPI_INT, retVal, TAG_INT, MPI_COMM_WORLD, &status);

      // make sure the TB is the correct size
      if (tempInt != *numPts)
      { // clear the currently allocated memory if needed
        for (i = *numPts - 1; i >= 0; i--)
        {
          clear_trackBack_sample(&(*TB)[i]);
        }

        // allocate to the correct size
        *TB = (trackBack_samples_t *)brealloc(*TB, tempInt * sizeof(trackBack_samples_t));

        // initialize EG - when doing the copying, they automatically get setup - so only need the ones that are using AMP with double precision
        for (i = 0; i < tempInt; i++)
          init_trackBack_sample(&(*TB)[i], 64);
      }

      // setup TB
      egStr_temp = (char **)bmalloc(1 * sizeof(char *));
      tbStr_temp = (char **)bmalloc(1 * sizeof(char *));
      tbDouble_temp = (double **)bmalloc(1 * sizeof(double *));
      egCoeff_temp = (comp_d **)bmalloc(1 * sizeof(comp_d *));
      tbCoeff_temp = (comp_d **)bmalloc(1 * sizeof(comp_d *));
      tbInt_temp = (int **)bmalloc(1 * sizeof(int *));
      for (i = 0; i < tempInt; i++)
      { // setup the locations
        *egStr_temp = &egStr[egStrLoc];
        *tbStr_temp = &tbStr[tbStrLoc];
        *tbDouble_temp = &tbDouble[tbDoubleLoc];
        *egCoeff_temp = &egCoeff[egCoeffLoc];
        *tbCoeff_temp = &tbCoeff[tbCoeffLoc];
        *tbInt_temp = &tbInt[tbIntLoc];

        // setup TB[i]
        cp_trackBack_samples_t_int(&(*TB)[i], &TB_int[i], egStr_temp, egCoeff_temp, tbDouble_temp, tbCoeff_temp, tbStr_temp, tbInt_temp, 0, 0, 1);
    
        // update locations
        egStrLoc += TB_int[i].EG_int.totalLength;
        tbStrLoc += TB_int[i].totalLength;
        tbDoubleLoc += TB_int[i].num_double;
        egCoeffLoc += egCount[i];
        tbCoeffLoc += TB_int[i].num_comp_d;
        tbIntLoc += TB_int[i].numSamples;
      }

      // free memory
      free(egStr_temp);
      free(tbStr_temp);
      free(tbDouble_temp);
      free(egCoeff_temp);
      free(tbCoeff_temp);
      free(tbInt_temp);
      free(egStr);
      free(tbStr);
      free(tbDouble);
      free(egCoeff);
      free(tbCoeff);
      free(tbInt);
      free(TB_int);
    }
    else
    { // clear the currently allocated memory if needed
      for (i = *numPts - 1; i >= 0; i--)
        clear_trackBack_sample(&(*TB)[i]);

      // clear TB
      free(*TB);
      *TB = NULL;
    }

    // setup numPts
    *numPts = tempInt;
  }

  // free the memory
  MPI_Type_free(&mpi_tb_t);
  MPI_Type_free(&mpi_comp_d);

  return retVal;
}

#endif


