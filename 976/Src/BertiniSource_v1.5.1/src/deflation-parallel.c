// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "diff.h"
#include "ppParse.h"

#define MAX_DEFLATION_ITS 250

typedef struct
{
  int  op;		// operation
  int  memLoc;        	// memory location
  int  in[2];         	// memory locations for the two input arguments
} inst;	 	   	// holds a single SLP instruction

typedef struct
{
  inst *insts;		// the instructions
  int num;	 	// the number of instructions in the list
} inst_list;		// holds all instructions for a single object

void clear_inst_list(inst_list *list, int num_list); // clear memory from list
void setup_insts_from_prog(inst_list **const_insts, int *const_size, inst_list **subfunc_insts, int *subfunc_size, inst_list **func_insts, int *func_size, inst_list **funcjac_insts, int *funcjac_size, inst_list **subfuncjac_insts, int *subfuncjac_size, prog_t *prog); // setup the instructions
void write_insts_to_file(FILE *FOUT, prog_t *prog, inst_list *const_insts, int const_size, inst_list *param_insts, int param_size, inst_list *subfunc_insts, int subfunc_size, inst_list *func_insts, int func_size, int room_before_temps, int total_room_needed); // write all the instructions
void print_insts(FILE *FOUT, inst_list *list, int num_list); // print the instructions in list

int corank_deflation(double *smallest_nonzero_SV, double *largest_zero_SV, point_d point0_d, point_mp point0_mp, int prec0, point_d point1_d, point_mp point1_mp, int prec1, comp_d time_d, comp_mp time_mp, eval_struct_d *e_d, eval_struct_mp *e_mp, tracker_config_t *T, FILE *OUT, prog_t *Prog);

// function to extend matrix K
void extend_deflation_matrix_K(mpq_t ****K_rat, int *K_rows, int *K_cols);

// functions used to extend the point
void extend_deflation_point(point_d output_pt_d, point_mp output_pt_mp, int *output_prec, point_d input_pt_d, point_mp input_pt_mp, mat_d J_d, mat_mp J_mp, mpq_t ***K, int K_rows, int K_cols, mpq_t ***B, int B_rows, int B_cols, int input_prec);
void extend_deflation_point_d(point_d output_pt, point_d input_pt, mat_d Jv, mat_d K, mat_d B);
void extend_deflation_point_mp(point_mp output_pt, point_mp input_pt, mat_mp Jv, mat_mp K, mat_mp B);

int deflation(int *deflations_needed, prog_t *deflatedProg, point_data_d *output_PD_d, point_data_mp *output_PD_mp, int *output_PD_prec, prog_t *origProg, int orig_corank, double smallest_nonzero_SV, double largest_zero_SV, point_data_d *input_PD_d, point_data_mp *input_PD_mp, int input_PD_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, mpq_t ***Kin, int Kin_rows, int Kin_cols, tracker_config_t *T, FILE *OUT, int deflation_iteration_bound)
/************************************************************************\
* USAGE: main deflation function - calls deflator as needed              *
* ARGUMENTS: the SLP structure, the point, K, and stuff for rank finding *
* RETURN VALUES: 0 - success, otherwise - failure                        *
* NOTES:  see deflation dev plan for details                             *
\************************************************************************/
{ 
  int retVal = 0;

  // initialize the number of deflations needed
  *deflations_needed = 0;

  // make sure that it is rank deficient
  if (orig_corank <= 0)
  { // copy origProg to deflatedProg
    cp_prog_t(deflatedProg, origProg); 

    // setup output_PD
    *output_PD_prec = input_PD_prec;
    if (input_PD_prec < 64)
    { // setup output_PD_d
      point_data_cp_d(output_PD_d, input_PD_d);
    }
    else
    { // setup output_PD_mp
      change_prec_point_data_mp(output_PD_mp, *output_PD_prec);
      point_data_cp_mp(output_PD_mp, input_PD_mp);
    }

    // return since we already have full rank
    return retVal;
  }

  // so, we have to do atleast 1 deflation
  int i, j, its, K_rows, K_cols, B_rows, B_cols, new_corank = orig_corank, currPrec = T->MPType == 1 ? T->Precision : 64, maxPrec = T->MPType == 1 ? T->Precision : T->AMP_max_prec;
  int time_prec, temp_approx_prec, temp_PD_prec, corank_offset = Kin_rows - Kin_cols; // change in corank based on linear slicing
  comp_d time_d;
  comp_mp time_mp;
  point_d temp_approx_d;
  point_mp temp_approx_mp;
  point_data_d temp_PD_d;
  point_data_mp temp_PD_mp;
  mat_d B_d;
  mpq_t ***K_rat = NULL, ***B_rat = NULL;
  prog_t input_prog;  
  eval_struct_d e_d;
  eval_struct_mp e_mp;

  init_mp(time_mp);
  init_point_d(temp_approx_d, 0);
  init_point_mp(temp_approx_mp, 0);
  init_point_data_d(&temp_PD_d, 0);
  init_point_data_mp(&temp_PD_mp, 0);
  init_mat_d(B_d, 0, 0);
  init_eval_struct_d(e_d, 0, 0, 0);
  init_eval_struct_mp(e_mp, 0, 0, 0);

  // copy origProg to input_prog to avoid any conflicts with deflatedProg
  cp_prog_t(&input_prog, origProg);

  // setup time
  time_prec = MAX(last_approx_prec, input_PD_prec);
  if (time_prec < 64)
  { // only setup _d
    set_d(time_d, input_PD_d->time);
  }
  else
  { // setup both _d & _mp
    setprec_mp(time_mp, time_prec);
    if (input_PD_prec < 64)
    {
      set_d(time_d, input_PD_d->time);
      d_to_mp(time_mp, time_d);
    }
    else
    {
      set_mp(time_mp, input_PD_mp->time);
      mp_to_d(time_d, time_mp);
    }
  }

  // copy last_approx to temp_approx
  temp_approx_prec = last_approx_prec;
  if (last_approx_prec < 64)
  { // copy to _d
    point_cp_d(temp_approx_d, last_approx_d);
  }
  else
  { // copy to _mp
    setprec_point_mp(temp_approx_mp, temp_approx_prec);
    point_cp_mp(temp_approx_mp, last_approx_mp);
  }

  // copy input_PD to temp_PD
  temp_PD_prec = input_PD_prec;
  if (input_PD_prec < 64)
  { // copy to _d
    point_data_cp_d(&temp_PD_d, input_PD_d);
  }
  else
  { // copy to _mp
    setprec_point_data_mp(&temp_PD_mp, temp_PD_prec);
    point_data_cp_mp(&temp_PD_mp, input_PD_mp);
  }

  // setup K_rat
  K_rows = Kin_rows;
  K_cols = Kin_cols;
  init_mat_rat(K_rat, K_rows, K_cols);
  mat_cp_rat(K_rat, Kin, Kin_rows, Kin_cols);

  //We proceed in a while loop which consists of a deflation followed by a rank check.
  //We continue until the corank is down to 0.
  its = 0;
  do
  { //DEVELOPMENT REMINDER:  K IS DETERMINED FROM THE LINEARS L, AND (FOR DEFLATION) WE ASSUME THAT L IS MOVED IN PARALLEL (ONLY THE CONSTANT MOVES)!!!!

    // setup B
    B_rows = new_corank;
    B_cols = K_cols;
    init_mat_rat(B_rat, B_rows, B_cols); 

    if (T->MPType == 0)
    {
      make_matrix_random_d(B_d, B_rows, B_cols);

      for (i = 0; i < B_rows; i++)
        for (j = 0; j < B_cols; j++)
        {
          mpq_set_d(B_rat[i][j][0], B_d->entry[i][j].r);
          mpq_set_d(B_rat[i][j][1], B_d->entry[i][j].i);
        }
    }
    else
    {
      make_matrix_random_rat(B_d, e_mp.Jp, B_rat, B_rows, B_cols, currPrec, maxPrec, 0, 0);
    }

    // extend the point temp_PD to output_PD - find the original jacobian
    if (temp_PD_prec < 64)
    {
      eval_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, temp_PD_d.point, time_d, &input_prog, evalProg_d_void);
    }
    else
    { // set the current precision
      initMP(temp_PD_prec);
      input_prog.precision = temp_PD_prec + 32;
      setprec_eval_struct_mp(e_mp, temp_PD_prec);

      eval_mp(e_mp.funcVals, e_mp.parVals, e_mp.parDer, e_mp.Jv, e_mp.Jp, temp_PD_mp.point, time_mp, &input_prog, evalProg_mp_void);
    }

    // use Jv, K & B to extend the point
    extend_deflation_point(output_PD_d->point, output_PD_mp->point, output_PD_prec, temp_PD_d.point, temp_PD_mp.point, e_d.Jv, e_mp.Jv, K_rat, K_rows, K_cols, B_rat, B_rows, B_cols, temp_PD_prec);

    // extend the point temp_approx - find the original jacobian
    if (temp_approx_prec < 64)
    {
      eval_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, temp_approx_d, time_d, &input_prog, evalProg_d_void);
    }
    else
    { // set the current precision
      initMP(temp_approx_prec);
      input_prog.precision = temp_approx_prec;
      setprec_eval_struct_mp(e_mp, temp_approx_prec);

      eval_mp(e_mp.funcVals, e_mp.parVals, e_mp.parDer, e_mp.Jv, e_mp.Jp, temp_approx_mp, time_mp, &input_prog, evalProg_mp_void);
    }
    // use Jv, K & B to extend the point
    extend_deflation_point(temp_approx_d, temp_approx_mp, &temp_approx_prec, temp_approx_d, temp_approx_mp, e_d.Jv, e_mp.Jv, K_rat, K_rows, K_cols, B_rat, B_rows, B_cols, temp_approx_prec);

    // do a single-stage deflation
    deflator(deflatedProg, &input_prog, K_rat, K_rows, K_cols, B_rat, B_rows, B_cols);

    // now that we have extended the points, it is time to find the corank
    new_corank = corank_deflation(&smallest_nonzero_SV, &largest_zero_SV, temp_approx_d, temp_approx_mp, temp_approx_prec, output_PD_d->point, output_PD_mp->point, *output_PD_prec, time_d, time_mp, &e_d, &e_mp, T, OUT, deflatedProg);

    // find the rank of this system
    i = MIN(deflatedProg->numFuncs, deflatedProg->numVars) - new_corank;

    // see what happens when the linear slices are added
    i = MIN(i + corank_offset, deflatedProg->numVars);

    // find the corank
    new_corank = MIN(deflatedProg->numFuncs + corank_offset, deflatedProg->numVars) - i;  

    // clear B
    clear_mat_rat(B_rat, B_rows, B_cols);

    // increment its
    its++;

    if (new_corank > 0)
    { // prepare for next iteration
      if (*output_PD_prec < 64)
      { // copy to _d
        point_data_cp_d(&temp_PD_d, output_PD_d);
      }
      else
      { // copy to _mp
        if (temp_PD_prec != *output_PD_prec)
        { // change precision
          setprec_point_data_mp(&temp_PD_mp, *output_PD_prec);
        }
        point_data_cp_mp(&temp_PD_mp, output_PD_mp);
      }
      temp_PD_prec = *output_PD_prec;

      // clear the old prog
      clearProg(&input_prog, T->MPType, 0);
      // setup to be the deflated (but not yet rank deficient) prog
      cp_prog_t(&input_prog, deflatedProg); 

      // extend K
      extend_deflation_matrix_K(&K_rat, &K_rows, &K_cols);
    }
  } while (new_corank > 0 && new_corank <= orig_corank && its < deflation_iteration_bound && its < MAX_DEFLATION_ITS);

  // set the number of deflations used
  *deflations_needed = its;

  // setup retVal
  if (new_corank > 0)
  { // deflation failed to produce a full rank Jacobian matrix
    retVal = 1;
  }
  else
  { // deflation worked properly
    retVal = 0;
  }

  // clear memory
  clear_mp(time_mp);
  clear_point_d(temp_approx_d);
  clear_point_mp(temp_approx_mp);
  clear_point_data_d(&temp_PD_d);
  clear_point_data_mp(&temp_PD_mp);
  clear_mat_d(B_d);
  clear_eval_struct_d(e_d);
  clear_eval_struct_mp(e_mp);
  clear_mat_rat(K_rat, K_rows, K_cols);

  return retVal;
}

void deflator(prog_t *new_prog, prog_t *old_prog, mpq_t ***K, int K_rows, int K_cols, mpq_t ***B, int B_rows, int B_cols)
/************************************************************************\
* USAGE: carries out one stage of deflation, overwriting ED              *
* ARGUMENTS: the SLP structure, the point, K, and the corank             *
* RETURN VALUES: == 0 means success, other means error                   *
* NOTES:  see deflation dev plan (7/25/07) for details                   *
*  started on 7/25/07 by Dan, in London, ON                              *
\************************************************************************/
{ 
  int i, j, k, m; 
  int *addr_conv_table = (int *)bmalloc(old_prog->memSize * sizeof(int)); //Conversion table for workspace array addresses (index is old addr, entry is new addr).
  int old_const_size = 0, old_subfunc_size = 0, old_func_size = 0, old_funcjac_size = 0, old_subfuncjac_size = 0;
  inst_list *old_const_insts = NULL, *old_subfunc_insts = NULL, *old_func_insts = NULL, *old_funcjac_insts = NULL, *old_subfuncjac_insts = NULL;
  int new_const_size = 0, new_subfunc_size = 0, new_func_size = 0;
  inst_list *new_const_insts = NULL, *new_subfunc_insts = NULL, *new_func_insts = NULL;

  int K_num_index, B_num_index, K_const_index, B_const_index, K_addr, B_addr, J_addr;
  int new_workspace_addr, old_workspace_addr;
  int ctr, Psi_dim, tmp_func_addr, tmp_V_addr, ind, room_before_temps;
  FILE *FOUT = NULL, *arrIN = NULL;  
  char ch;

  // setup Psi_dim
  Psi_dim = K_cols;

/////////////////////////////////////////
//  GRABBING THE INSTS FROM THE OLD SLP
/////////////////////////////////////////

  //The idea is that we are going to use the old insts (in old_*_insts structs) to write the new insts (to new_*_insts structs).
  //After everything is written to the new_*_insts structs, we'll throw them together into the SLP, send it off to diff_deflatable(), and call it a day.
  //We will initialize and use the new_*_insts structs later as we fill them.

  // setup the structures for the old SLP
  setup_insts_from_prog(&old_const_insts, &old_const_size, &old_subfunc_insts, &old_subfunc_size, &old_func_insts, &old_func_size, &old_funcjac_insts, &old_funcjac_size, &old_subfuncjac_insts, &old_subfuncjac_size, old_prog);

  //For future reference, the order of the instructions is constants, parameters, subfunctions, functions, subfunc derivs, func derivs.

////////////////////////////////////////
//  SETTING UP THE NEW PROG AND 
//  POPULATING THE CONVERSION TABLE
////////////////////////////////////////

  //At this point, we are done grabbing the instructions from old_prog, so we may start building new_prog, which will contain the deflated system.
  //First we set up the bits of new_prog (the new SLP structure).  
  //new_prog->prog, size, and memSize will be handled after all instructions are determined.
  //nums (the list of numbers) will be dealt with below when we start adding in new numbers for B and K.
  new_prog->precision = old_prog->precision; // The precision at which evaluation should occur 
  new_prog->num_var_gps = old_prog->num_var_gps + 1;  
  new_prog->var_gp_sizes = (int *)bmalloc(new_prog->num_var_gps * sizeof(int));
  for (i=0; i<old_prog->num_var_gps; i++)
    new_prog->var_gp_sizes[i] = old_prog->var_gp_sizes[i];
  new_prog->var_gp_sizes[i] = Psi_dim;
  new_prog->index_of_first_number_for_proj_trans = old_prog->index_of_first_number_for_proj_trans;  //COULD BE WRONG - IS THIS USED???!!!
  //These next few aren't being used yet, so default to old ones for now...needs to be changed later if used!!!  
  new_prog->numInstAtEndUpdate = old_prog->numInstAtEndUpdate; 
  new_prog->numInstAtEndParams = old_prog->numInstAtEndParams; 
  new_prog->numInstAtEndFnEval = old_prog->numInstAtEndFnEval; 
  new_prog->numInstAtEndPDeriv = old_prog->numInstAtEndPDeriv; 
  new_prog->numInstAtEndJvEval = old_prog->numInstAtEndJvEval; 

  //Now for the numbers and locations of various types (plus the new and old workspace starting points).
  //We also populate the conversion table at this point.
  new_workspace_addr = old_workspace_addr = 0;  //We'll increment these as we go through this block to get to the new and old first workspace addrs.

  /////////////////////////VARS
  new_prog->numVars = old_prog->numVars + K_cols;  //Since we're adding in K_cols (=Psi_rows) new variables.
  new_prog->inpVars = old_prog->inpVars; //should just be 0 since vars go first.
  new_workspace_addr += new_prog->numVars;
  old_workspace_addr += old_prog->numVars;
  for (i=0; i<old_prog->numVars; i++)
    addr_conv_table[old_prog->inpVars + i] = new_prog->inpVars + i;
  /////////////////////////PATHVARS
  new_prog->numPathVars = 0;   //Shouldn't be any path vars (not supported) - we'll check that next.
  new_prog->inpPathVars = new_prog->inpVars + new_prog->numVars;
  if (old_prog->numPathVars > 0) 
  {
    printf("WARNING:  Bertini has detected a path variable while deflating the system.  Path variables are not supported in deflation, so any further computation is untrustworthy!\n"); 
    bexit(ERROR_INPUT_SYNTAX);
  }
  new_workspace_addr += new_prog->numPathVars;
  old_workspace_addr += old_prog->numPathVars;
  //Nothing to add to conversion table since there are no pathvars!
  /////////////////////////PARAMS
  new_prog->numPars = 0;   //Shouldn't be any parameters (not supported) - we'll check that next.
  new_prog->evalPars = new_prog->inpPathVars + new_prog->numPathVars;
  if (old_prog->numPars > 0)
  {
    printf("WARNING:  Bertini has detected a parameter while deflating the system.  Parameters are not supported in deflation, so any further computation is untrustworthy!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }
  new_workspace_addr += new_prog->numPars;
  old_workspace_addr += old_prog->numPars;
  //Nothing to add to conversion table since there are no parameters!
  /////////////////////////FUNCS
  new_prog->numFuncs = old_prog->numFuncs + old_prog->numFuncs + B_rows;   //The original functions + the Jacobian times stuff + the third line (B\Psi-I).
  new_prog->evalFuncs = new_prog->evalPars + new_prog->numPars;
  new_workspace_addr += new_prog->numFuncs;
  old_workspace_addr += old_prog->numFuncs;
  /////////////////////////DPARAMS
  new_prog->evalDPars = new_prog->evalFuncs + new_prog->numFuncs;
  new_workspace_addr += new_prog->numPars;  //numPars since numDPars=numPars (one each) - should be 0 anyway!
  old_workspace_addr += old_prog->numPars;
  //Nothing to add to conversion table since there are no parameters (and hence no parameter derivatives)!
  /////////////////////////JFUNCS(V)
  new_prog->evalJVars = new_prog->evalDPars + new_prog->numPars;
  new_workspace_addr += new_prog->numFuncs * new_prog->numVars;
  old_workspace_addr += old_prog->numFuncs * old_prog->numVars;
  /////////////////////////JFUNCS(P)
  new_prog->evalJPars = new_prog->evalJVars + new_prog->numFuncs*new_prog->numVars;
  new_workspace_addr += new_prog->numFuncs * new_prog->numPars;  //numPars=0, so this and next line are worthless
  old_workspace_addr += old_prog->numFuncs * old_prog->numPars;
  //Nothing to add to conversion table since there are no parameters (and hence no derivatives of funcs w.r.t. parameters)!
  /////////////////////////CONSTANTS
  new_prog->numConsts = old_prog->numConsts + K_rows*K_cols + B_rows*B_cols;  //The original consts + consts for K and B.
  new_prog->constAddr = new_prog->evalJPars + new_prog->numFuncs * new_prog->numPars;  
  new_prog->IAddr = new_prog->constAddr;  //since I is always the first constant
  K_addr = new_prog->constAddr + old_prog->numConsts;  //Keeps track of first addr in K.
  B_addr = K_addr + K_rows*K_cols;  //Keeps track of first addr of B.
  new_workspace_addr += new_prog->numConsts;
  old_workspace_addr += old_prog->numConsts;
  for (i=0; i<old_prog->numConsts; i++)
    addr_conv_table[old_prog->constAddr + i] = new_prog->constAddr + i;
  /////////////////////////NUMS
  new_prog->numNums = old_prog->numNums + 2*K_rows*K_cols + 2*B_rows*B_cols;  //The original consts + two for each entry of K and B (real, imag parts).
  new_prog->numAddr = new_prog->constAddr + new_prog->numConsts;
  new_workspace_addr += new_prog->numNums;
  old_workspace_addr += old_prog->numNums;
  for (i=0; i<old_prog->numNums; i++)
    addr_conv_table[old_prog->numAddr + i] = new_prog->numAddr + i;
  /////////////////////////SUBFUNCS
  new_prog->numSubfuncs = old_prog->numSubfuncs + old_prog->numFuncs + old_prog->numVars * old_prog->numFuncs + old_prog->numVars * old_prog->numSubfuncs;   //old subfuncs + one for each old func + one for each entry of old jacobian + one for each derivative of old subfuncs
  new_prog->evalSubs = new_prog->numAddr + new_prog->numNums;
  new_workspace_addr += new_prog->numSubfuncs;
  old_workspace_addr += old_prog->numSubfuncs;
  ctr = new_prog->evalSubs;  //set up ctr as the start addr for subfuncs - we can't do all subfuncs now - see JSUBS(V) next
  for (i=0; i<old_prog->numSubfuncs; i++)  //Just copy in the old subfuncs first
    addr_conv_table[old_prog->evalSubs + i] = ctr++;
  for (i=0; i<old_prog->numFuncs; i++)  //Next we add one to ctr for each of the old funcs 
    addr_conv_table[old_prog->evalFuncs + i] = ctr++;  //We record the function addrs in the appropriate new subfunc locations (for reasons described below) - we handle them in a special way below.
  /////////////////////////JSUBS(V)
  new_prog->evalJSubsV = new_prog->evalSubs + new_prog->numSubfuncs;
  new_workspace_addr += new_prog->numSubfuncs * new_prog->numVars;
  old_workspace_addr += old_prog->numSubfuncs * old_prog->numVars;
  for (i=0; i<old_prog->numSubfuncs; i++)  //Next we add the new subfuncs to the conv table corresp. to old subfunc derivs
    for (j=0; j<old_prog->numVars; j++)  //We loop over old subfuncs (i) and old vars (j)
      addr_conv_table[old_prog->evalJSubsV + old_prog->numVars*i + j] = ctr++;  
  J_addr = ctr;  //J_addr keeps track of first addr of Jacobian.
  for (i=0; i<old_prog->numFuncs; i++)  //Finally we add the old function derivatives - first loop over old funcs
    for (j=0; j<old_prog->numVars; j++)  //Now loop over old vars
      addr_conv_table[old_prog->evalJVars + old_prog->numVars*i + j] = ctr++;  //Just figure out the correct old addr
  /////////////////////////JSUBS(P)
  new_prog->evalJSubsP = new_prog->evalJSubsV + new_prog->numSubfuncs*new_prog->numVars;
  new_workspace_addr += new_prog->numSubfuncs * new_prog->numPars;  //numPars=0, so this and next line are worthless
  old_workspace_addr += old_prog->numSubfuncs * old_prog->numPars;
  //Nothing to add to conversion table since there are no parameters (and hence no derivatives of subfuncs w.r.t. parameters)!  

  //We remember the beginning of the workspace:
  room_before_temps = new_workspace_addr;

  //Finally, the rest of the old workspace was set aside as temp memLocs, so we just copy that into the front of 
  //the new space designated for that purpose....
  //NOTE:  Because of doing this this way, all old temp memLocs will come first, followed by all new ones needed 
  //for the deflated system, even if some memLocs for the deflated system are used before the old ones.
  //That shouldn't matter, but it is worth noting....
  for (i=old_workspace_addr; i<old_prog->memSize; i++)
    addr_conv_table[i] = new_workspace_addr++;

///////////////////////////////////////////
//  CREATING INSTS FOR THE NEW SLP
///////////////////////////////////////////

  //Now we start to populate the SLP with the (old and) new data.
  //Let's start with numbers and constants.  We still have all the old numbers, plus the new ones for K and B.
  //We will store the entries of K and B as actual constants (not just numbers) since that will make life easier 
  //when writing instructions using the entries of K and B.
  //After setting up the numbers, we have all of the old constant update instructions plus the new ones for K and B.

////////////////////////////////////////
//  ADDING NUMBERS TO NEW_PROG
////////////////////////////////////////

  //First initialize new_prog->nums to the new number of numbers:
  new_prog->numNums = old_prog->numNums + 2*K_rows*K_cols + 2*B_rows*B_cols; //old numNums + those for K + those for B (real + im of each!)
  new_prog->nums = (num_t *)bmalloc(new_prog->numNums * sizeof(num_t)); 
  for (i=0; i<old_prog->numNums; i++)  //Just copy in the old numbers
  {
    mpq_init(new_prog->nums[i].rat);
    mpq_set(new_prog->nums[i].rat, old_prog->nums[i].rat);
    new_prog->nums[i].currPrec = old_prog->nums[i].currPrec;
    mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
    mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);
  } 

  //Now we add all new numbers from K (provided) and for B (random rationals), in rational number format
  //Here is K:
  K_num_index = i;  //Remember where in the nums K starts
  for (j=0; j<K_rows; j++)  //Add in the nums for K - 1,1(real) 1,1(im) 1,2(re), etc.
    for (k=0; k<K_cols; k++)
    { //The real part:
      mpq_init(new_prog->nums[i].rat);
      mpq_set(new_prog->nums[i].rat, K[j][k][0]);
      new_prog->nums[i].currPrec = old_prog->nums[0].currPrec;  //pull the precision off of the first (old) number (!!!).
      mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
      mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);
      i++;
      //The imag part:
      mpq_init(new_prog->nums[i].rat);
      mpq_set(new_prog->nums[i].rat, K[j][k][1]);
      new_prog->nums[i].currPrec = old_prog->nums[0].currPrec;  //pull the precision off of the first (old) number (!!!).
      mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
      mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);
      i++;
    }

  //Here is B:
  B_num_index = i;  //Remember where in the nums B starts

  for (j=0; j<B_rows; j++)  //Add in the nums for B - 1,1(real) 1,1(im) 1,2(re), etc.
    for (k=0; k<B_cols; k++)
    { // real part     
      mpq_init(new_prog->nums[i].rat);
      mpq_set(new_prog->nums[i].rat, B[j][k][0]); 
      new_prog->nums[i].currPrec = old_prog->nums[0].currPrec;  //pull the precision off of the first (old) number (!!!).
      mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
      mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);
      i++;
      // imag part
      mpq_init(new_prog->nums[i].rat);
      mpq_set(new_prog->nums[i].rat, B[j][k][1]);
      new_prog->nums[i].currPrec = old_prog->nums[0].currPrec;  //pull the precision off of the first (old) number (!!!).
      mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
      mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);
      i++;
    }

///////////////////////////////////////////
// ADDING CONST EVAL INSTS TO NEW SLP
///////////////////////////////////////////

  //Now that we have all of our numbers chosen and stored in new_prog, it is time to work on the insts of the SLP.
  //First we must handle the consts, starting by preparing new_const_insts:
  new_const_size = new_prog->numConsts - 2; // -2 b/c I & Pi is #0 & #1 (no update insts)
  new_const_insts = (inst_list *)bmalloc(new_const_size * sizeof(inst_list));  

  // copy over all of the old constant evaluations
  for (i=0; i<old_const_size; i++)  
  {
    new_const_insts[i].num = old_const_insts[i].num;
    new_const_insts[i].insts = (inst *)bmalloc(new_const_insts[i].num * sizeof(inst));
    for (j=0; j<new_const_insts[i].num; j++)  //Now just copy over the num insts.
    {
      new_const_insts[i].insts[j].op = old_const_insts[i].insts[j].op;
      new_const_insts[i].insts[j].memLoc = addr_conv_table[old_const_insts[i].insts[j].memLoc];
      new_const_insts[i].insts[j].in[0] = addr_conv_table[old_const_insts[i].insts[j].in[0]];
      new_const_insts[i].insts[j].in[1] = addr_conv_table[old_const_insts[i].insts[j].in[1]];
    }
  }	  
  //Second (and finally for consts), we have const evals for K & B (there are real and imag rational nums for each const)
  //We keep i as the const #.  K first:
  K_const_index = i;
  m = K_num_index;  //m will walk through the indices for the numbers for K, starting (of course) at the first.
  for (j=0; j<K_rows; j++)
  {
    for (k=0; k<K_cols; k++)
    {
      new_const_insts[i].num = 2; // each new constant creates 2 instructions (1) imag_part * I & (2) real_part + imag_part * I
      new_const_insts[i].insts = (inst *)bmalloc(new_const_insts[i].num * sizeof(inst));

      // imag_part * I
      new_const_insts[i].insts[0].op = '*'; // multiply
      new_const_insts[i].insts[0].memLoc = new_workspace_addr;
      new_const_insts[i].insts[0].in[0] = new_prog->numAddr + m + 1; // for imaginary part of the number
      new_const_insts[i].insts[0].in[1] = new_prog->IAddr;

      // real_part + imag_part * I
      new_const_insts[i].insts[1].op = '+'; // add
      new_const_insts[i].insts[1].memLoc = new_prog->constAddr + i + 2; // store result in appropriate addr (+2 b/c of I & Pi)
      new_const_insts[i].insts[1].in[0] = new_workspace_addr;
      new_const_insts[i].insts[1].in[1] = new_prog->numAddr + m; // for real part of the number

      // increment i & m
      i++;
      m += 2;
    }
  }
  //Here we go with B:	    
  B_const_index = i;
  m = B_num_index;
  for (j=0; j<B_rows; j++) 
  {
    for (k=0; k<B_cols; k++)
    {
      new_const_insts[i].num = 2; // each new constant creates 2 instructions (1) imag_part * I & (2) real_part + imag_part * I
      new_const_insts[i].insts = (inst *)bmalloc(new_const_insts[i].num * sizeof(inst));

      // imag_part * I
      new_const_insts[i].insts[0].op = '*'; // multiply
      new_const_insts[i].insts[0].memLoc = new_workspace_addr;
      new_const_insts[i].insts[0].in[0] = new_prog->numAddr + m + 1; // for imaginary part of the number
      new_const_insts[i].insts[0].in[1] = new_prog->IAddr;

      // real_part + imag_part * I
      new_const_insts[i].insts[1].op = '+'; // add
      new_const_insts[i].insts[1].memLoc = new_prog->constAddr + i + 2; // store result in appropriate addr (+2 b/c of I & Pi)
      new_const_insts[i].insts[1].in[0] = new_workspace_addr;
      new_const_insts[i].insts[1].in[1] = new_prog->numAddr + m; // for real part of the number

      // increment i & m
      i++;
      m += 2;
    }
  }

////////////////////////////////////////////////
//  ADDING PARAMETER INSTS TO NEW SLP
////////////////////////////////////////////////
  //There should be no parameters as that case is not yet supported for deflation.
  //This block is here to indicate that this is where such insts would go, if and when we choose to implement that case.

///////////////////////////////////////////////////////
//  ADDING SUBFUNC AND FUNC INSTS TO NEW SLP
///////////////////////////////////////////////////////
  //WARNING:  AT FIRST GLANCE, THE FOLLOWING BUNCH OF INSTS MIGHT SEEM TO BE IN A WEIRD ORDER!
  //The obvious thing to do is to write out the various insts in the following order:
  //old subfuncs, old subfunc derivatives, old func derivatives (all three of those would now be called subfuncs), then
  //old funcs, two types of new funcs (everything on this line would now be called funcs).
  //The problem with this is that the old func derivatives depend upon the old func eval insts.  The problem with this is 
  //that (using the order just mentioned) the old func derivs would depend upon something not yet defined, making 
  //evaluation impossible (and differentiation, for that matter).
  //The resolution that I am using is to add in even more new subfuncs, one for each old func.  The order becomes:
  //old subfuncs, new subfuncs for old funcs, old subfunc derivs, old func derivs (all four are now subfuncs), then
  //old funcs (just set them equal to corresp. new subfuncs for old funcs), and then other new funcs (all called funcs).
  //Use bertini_def_backup_27nov07.tar.gz in case this fails miserably....

/////
/////  ADDING THE OLD SUBFUNCS AS NEW SUBFUNCS
/////

  //There are four kinds of subfuncs now - the old ones, the old funcs, the old JSubsV, and the old JV
  new_subfunc_size = new_prog->numSubfuncs;
  new_subfunc_insts = (inst_list *)bmalloc(new_subfunc_size * sizeof(inst_list));  
  //The first subfunc insts are for the old subfuncs -> we copy over all of the old subfunc evaluations.
  for (i=0; i<old_prog->numSubfuncs; i++)  //i is the new subfunc # on which we are working.
  {
    new_subfunc_insts[i].num = old_subfunc_insts[i].num;
    new_subfunc_insts[i].insts = (inst *)bmalloc(new_subfunc_insts[i].num * sizeof(inst));
    for (j=0; j<new_subfunc_insts[i].num; j++)  //Now just copy over the num insts.
    {
      new_subfunc_insts[i].insts[j].op = old_subfunc_insts[i].insts[j].op;
      new_subfunc_insts[i].insts[j].memLoc = addr_conv_table[old_subfunc_insts[i].insts[j].memLoc];
      new_subfunc_insts[i].insts[j].in[0] = addr_conv_table[old_subfunc_insts[i].insts[j].in[0]];
      new_subfunc_insts[i].insts[j].in[1] = addr_conv_table[old_subfunc_insts[i].insts[j].in[1]];
    }
  }	  

/////
/////  ADDING THE OLD FUNCS AS NEW SUBFUNCS
/////
 
  //The second subfunc insts are for the old funcs -> we copy over all of the old func evaluations.
  //We keep i as the running index in new_subfunc_insts (from the previous block and use k to walk through the old func eval insts.
  for (k=0; k<old_prog->numFuncs; k++)  //k is the old func # on which we are working.
  {
    new_subfunc_insts[i].num = old_func_insts[k].num;
    new_subfunc_insts[i].insts = (inst *)bmalloc(new_subfunc_insts[i].num * sizeof(inst));
    for (j=0; j<new_subfunc_insts[i].num; j++)  //Now just copy over the num insts.
    {
      new_subfunc_insts[i].insts[j].op = old_func_insts[k].insts[j].op;
      new_subfunc_insts[i].insts[j].memLoc = addr_conv_table[old_func_insts[k].insts[j].memLoc];
      new_subfunc_insts[i].insts[j].in[0] = addr_conv_table[old_func_insts[k].insts[j].in[0]];
      new_subfunc_insts[i].insts[j].in[1] = addr_conv_table[old_func_insts[k].insts[j].in[1]];
    }
    i++;
  }	  

/////
/////  ADDING THE OLD SUBFUNC DERIV INSTS AS NEW SUBFUNCS
/////

  //The next type of new subfunc comes from the old subfunc deriv instructions
  //We keep i as the running index in new_subfunc_insts (from the previous block) and use k to walk through the subfunc derivs.

  for (k=0; k<old_prog->numSubfuncs * old_prog->numVars; k++)  //k is the subfunc # on which we are working.
  {
    new_subfunc_insts[i].num = old_subfuncjac_insts[k].num;
    new_subfunc_insts[i].insts = (inst *)bmalloc(new_subfunc_insts[i].num * sizeof(inst));
    for (j=0; j<new_subfunc_insts[i].num; j++)  //Now just copy over the num insts.
    {
      new_subfunc_insts[i].insts[j].op = old_subfuncjac_insts[k].insts[j].op;
      new_subfunc_insts[i].insts[j].memLoc = addr_conv_table[old_subfuncjac_insts[k].insts[j].memLoc];
      new_subfunc_insts[i].insts[j].in[0] = addr_conv_table[old_subfuncjac_insts[k].insts[j].in[0]];
      new_subfunc_insts[i].insts[j].in[1] = addr_conv_table[old_subfuncjac_insts[k].insts[j].in[1]];
    }
    i++;
  }

/////
/////  ADDING THE OLD JACOBIAN INSTS AS NEW SUBFUNCS
/////  

  //Finally we must record the old jacobian instructions as subfunction instructions  
  //We keep i as the running index in new_subfunc insts (from the previous block) and use k to walk through the jacobian, as in the last section.

  for (k=0; k<old_prog->numFuncs * old_prog->numVars; k++)  //k is the subfunc # on which we are working.
  {
    new_subfunc_insts[i].num = old_funcjac_insts[k].num;
    new_subfunc_insts[i].insts = (inst *)bmalloc(new_subfunc_insts[i].num * sizeof(inst));
    for (j=0; j<new_subfunc_insts[i].num; j++)  //Now just copy over the num insts.
    {
      new_subfunc_insts[i].insts[j].op = old_funcjac_insts[k].insts[j].op;
      new_subfunc_insts[i].insts[j].memLoc = addr_conv_table[old_funcjac_insts[k].insts[j].memLoc];
      new_subfunc_insts[i].insts[j].in[0] = addr_conv_table[old_funcjac_insts[k].insts[j].in[0]];
      new_subfunc_insts[i].insts[j].in[1] = addr_conv_table[old_funcjac_insts[k].insts[j].in[1]];
    }
    i++;
  }

/////
/////  ADDING THE OLD FUNCTION INSTS AS NEW FUNCS
/////
  //Now we can just copy over the function instructions from the original system
  new_func_size = new_prog->numFuncs;
  new_func_insts = (inst_list *)bmalloc(new_func_size * sizeof(inst_list));  

  //The first func insts for the deflated system are just the old funcs.
  //The insts to eval the old funcs have already been used above to write out new subfuncs (which just equal the old funcs)
  //As a result, all that remains for the new funcs corresp. to the old ones is to set the new funcs equal to the corresp. 
  //new subfunc addrs.
  for (i=0; i<old_prog->numFuncs; i++)  //i is the func # on which we are working.
  {
    new_func_insts[i].num = 1;
    new_func_insts[i].insts = (inst *)bmalloc(new_func_insts[i].num * sizeof(inst));

    new_func_insts[i].insts[0].op = '=';  // store
    new_func_insts[i].insts[0].memLoc = new_prog->evalFuncs + i;  //This is just the addr of this new func.
    new_func_insts[i].insts[0].in[0] = addr_conv_table[old_prog->evalFuncs+i];  //This gives us the corresp. subfunc addr!
    new_func_insts[i].insts[0].in[1] = 0;  //Standard in[1] for op "=".
  }

/////
/////  ADDING THE NEW MIDDLE DEFLATION FUNC INSTS AS NEW FUNCS
/////

  //Now we write the instructions for the second line of the deflated system (J*K*Psi)

  //The outer loop below (indexed by k) loops through the rows of V := J*K, multiplying to Psi to get the kth func in this part of the deflation system.
  //The next loop in (indexed by j) loops over the columns of V (or K, if  you like) (= dimension of Psi) to actually do this multiplication (V[k,_]*Psi).
  //The innermost loop (indexed by m) loops over the rows of K (= cols of J) to compute J[k,_]*K[_,j] for each j.

  // NOTE: (JDH - 12/18/07) This could be made more efficient by doing K*Psi and then J*(K*Psi) (n^3 + n^2 vs. 2*n^2)
  
  for (k = 0; k < old_prog->numFuncs; k++) //J has the same number of rows as the old # of functions.
  {
    i = old_prog->numFuncs + k; // find the function number

    new_func_insts[i].num = Psi_dim * (2 * K_rows + 1); //See the instructions written below for details.
    new_func_insts[i].insts = (inst *)bmalloc(new_func_insts[i].num * sizeof(inst));
    ind = 0;  //ind keeps track of the instruction #.
    tmp_func_addr = new_workspace_addr++;

    for (j = 0; j < Psi_dim; j++)
    { // compute V[k,j] = (J*K)[k,j] - and store it in tmp_V_addr
      tmp_V_addr = new_workspace_addr;
      new_workspace_addr++;

      // do the first multiplication
      new_func_insts[i].insts[ind].op = '*'; // multiply
      new_func_insts[i].insts[ind].memLoc = tmp_V_addr;
      new_func_insts[i].insts[ind].in[0] = J_addr + old_prog->numVars * k; // J[k,0]
      new_func_insts[i].insts[ind].in[1] = K_addr + j; // K[0,j]
      // increment ind
      ind++;

      // do the other instructions
      for (m = 1; m < K_rows; m++)
      { // J[k,m] * K[m,j]
        new_func_insts[i].insts[ind].op = '*'; // multiply
        new_func_insts[i].insts[ind].memLoc = new_workspace_addr;
        new_func_insts[i].insts[ind].in[0] = J_addr + old_prog->numVars * k + m;  // J[k,m].
        new_func_insts[i].insts[ind].in[1] = K_addr + K_cols * m + j;  // K[m,j].
        // increment
        new_workspace_addr++;
        ind++;

        // add to sum
        new_func_insts[i].insts[ind].op = '+'; // addition
        new_func_insts[i].insts[ind].memLoc = new_workspace_addr;
        new_func_insts[i].insts[ind].in[0] = tmp_V_addr;
        new_func_insts[i].insts[ind].in[1] = new_workspace_addr - 1;
        // update and increment
        tmp_V_addr = new_workspace_addr;
        new_workspace_addr++;
        ind++;
      }

      // so we have computed V[k,j], we just multiply by Psi[j] and add it to the running sum for the new function.
      if (j == 0)
      { // the first multiplication will initialize the total - V[k,0] * Psi[0]
        new_func_insts[i].insts[ind].op = '*'; // multiply
        new_func_insts[i].insts[ind].memLoc = tmp_func_addr;
        new_func_insts[i].insts[ind].in[0] = tmp_V_addr; // V[k,0]
        new_func_insts[i].insts[ind].in[1] = new_prog->inpVars + old_prog->numVars; // Psi[0]
        // increment ind
        ind++;
      }
      else
      { // multiply and add on to sum
        new_func_insts[i].insts[ind].op = '*'; // multiply
        new_func_insts[i].insts[ind].memLoc = new_workspace_addr;
        new_func_insts[i].insts[ind].in[0] = tmp_V_addr; // V[k,j]
        new_func_insts[i].insts[ind].in[1] = new_prog->inpVars + old_prog->numVars + j;  // Psi[j]
        // increment
        ind++;
        new_workspace_addr++;

        // add to sum
        new_func_insts[i].insts[ind].op = '+'; // addition
        new_func_insts[i].insts[ind].memLoc = new_workspace_addr;
        new_func_insts[i].insts[ind].in[0] = tmp_func_addr;
        new_func_insts[i].insts[ind].in[1] = new_workspace_addr - 1;
        // update and increment
        tmp_func_addr = new_workspace_addr;
        new_workspace_addr++;
        ind++;
      }
    }
    // Finally, we set the actual function location equal to the tmp_func_addr:
    new_func_insts[i].insts[ind].op = '='; // store
    new_func_insts[i].insts[ind].memLoc = new_prog->evalFuncs + i; // ith function
    new_func_insts[i].insts[ind].in[0] = tmp_func_addr;
    new_func_insts[i].insts[ind].in[1] = 0;  //since = is unary.
    // increment ind
    ind++;
  }

/////
/////  ADDING THE FINAL NEW DEFLATION FUNC INSTS AS NEW FUNCS
/////

  //We finish the function eval insts by writing the third line of the deflated system (B*Psi-I)
  //This is much simpler than the last block!
  //i keeps going as the index within new_func_insts of the new function

  for (k = 0; k < B_rows; k++) //We just got through each of the B_rows functions (indexed by k).
  { // find the function number
    i = 2 * old_prog->numFuncs + k; 

    new_func_insts[i].num = 2 * (B_cols + 1);  //See the instructions written below for details.
    new_func_insts[i].insts = (inst *)bmalloc(new_func_insts[i].num * sizeof(inst));
    ind = 0;  //ind keeps track of the instruction #.
    tmp_func_addr = new_workspace_addr++;

    // initialize to -1
    new_func_insts[i].insts[ind].op = 'N'; // negate
    new_func_insts[i].insts[ind].memLoc = tmp_func_addr;
    new_func_insts[i].insts[ind].in[0] = new_prog->numAddr + 1; // '1' is always the second 'number' 
    new_func_insts[i].insts[ind].in[1] = 0; // unary negate
    // increment ind
    ind++;
 
    for (j = 0; j < B_cols; j++)
    { // do sum += B[k,j] * Psi[j]
      new_func_insts[i].insts[ind].op = '*'; // multiply
      new_func_insts[i].insts[ind].memLoc = new_workspace_addr;
      new_func_insts[i].insts[ind].in[0] = B_addr + B_cols * k + j;  // B[k,j].
      new_func_insts[i].insts[ind].in[1] = new_prog->inpVars + old_prog->numVars + j;  // Psi[j]
      // increment
      new_workspace_addr++;
      ind++;

      // add to sum
      new_func_insts[i].insts[ind].op = '+'; // addition
      new_func_insts[i].insts[ind].memLoc = new_workspace_addr;
      new_func_insts[i].insts[ind].in[0] = tmp_func_addr;
      new_func_insts[i].insts[ind].in[1] = new_workspace_addr - 1;
      // update & increment
      tmp_func_addr = new_workspace_addr;
      new_workspace_addr++;
      ind++;
    }
    // store to proper location
    new_func_insts[i].insts[ind].op = '='; // store
    new_func_insts[i].insts[ind].memLoc = new_prog->evalFuncs + i;
    new_func_insts[i].insts[ind].in[0] = tmp_func_addr;
    new_func_insts[i].insts[ind].in[1] = 0; // unary store
  }

///////////////////////////////////////////////////
//  COPYING ALL NEW INSTS TO FILE (FOR DIFF)
///////////////////////////////////////////////////
  //Now we write the instructions to a file so that we may call diff (which reads from a file)

  // delete old files
  remove("finalFile.out");
  remove("jacV.out");
  remove("jacP.out");
  remove("arr.out");

  // open up finalFile.out
  FOUT = fopen("finalFile.out", "w");
  // write the instructions to FOUT
  write_insts_to_file(FOUT, new_prog, new_const_insts, new_const_size, NULL, 0, new_subfunc_insts, new_subfunc_size, new_func_insts, new_func_size, room_before_temps, new_workspace_addr);
  // close FOUT
  fclose(FOUT);

////////////////////////////////////////////////
//  CALLING DIFF WITH NEW SLP
////////////////////////////////////////////////
  //We call diff_deflatable() to finish off the new SLP (i.e., to add in the jacobian insts)
  diff_deflatable(1, 0);

////////////////////////////////////////////////
//  COPYING COMPLETE SLP BACK INTO STRUCT
////////////////////////////////////////////////
  //Diff writes to a file all of the evaluation and differentiation instructions, so we just need to read that into our final (new) SLP.

  arrIN = fopen("arr.out", "r");
  if (arrIN == NULL)
  {
    printf("ERROR: The file containing the straight-line program does not exist!\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // get the memSize - move through arr.out until we find the first 'd' in TotalRoomNeeded
  do
  {
    ch = fgetc(arrIN);
  } while (ch != 'd');
  // read in the size
  fscanf(arrIN, "ed %d", &new_prog->memSize);

  // read in the number of instructions - move through arr.out until we find the 'I' in NUMINST
  do
  {
    ch = fgetc(arrIN);
  } while (ch != 'I');
  fscanf(arrIN, "NST %d", &new_prog->size);

  // grab all of the instructions.
  rewind(arrIN);
  new_prog->prog = (int *)bcalloc(new_prog->size, sizeof(int));
  for (i = 0; i < new_prog->size; i++)
    fscanf(arrIN, "%d ", &new_prog->prog[i]);

  // close arr.out
  fclose(arrIN);

  // clear the memory
  free(addr_conv_table);
  clear_inst_list(old_const_insts, old_const_size);
  free(old_const_insts);
  clear_inst_list(old_subfunc_insts, old_subfunc_size);
  free(old_subfunc_insts);
  clear_inst_list(old_func_insts, old_func_size);
  free(old_func_insts);
  clear_inst_list(old_funcjac_insts, old_funcjac_size);
  free(old_funcjac_insts);
  clear_inst_list(old_subfuncjac_insts, old_subfuncjac_size);
  free(old_subfuncjac_insts);

  clear_inst_list(new_const_insts, new_const_size);
  free(new_const_insts);
  clear_inst_list(new_subfunc_insts, new_subfunc_size);
  free(new_subfunc_insts);
  clear_inst_list(new_func_insts, new_func_size);
  free(new_func_insts);

  return;
}

void get_constant_ints(prog_t *prog, int *tmp_insts, inst_list *insts, int startAddr, int numConsts, int *prog_index)
/************************************************************************\
* USAGE: store the instructions for setting up the constants (update)    *
* ARGUMENTS:                                                             *
* RETURN VALUES:                                                         *
* NOTES: since update instructions do not have to follow an order, we    *
* have this special function to make sure that we read them in correctly *
\************************************************************************/
{
  int i = 0, j = *prog_index, k = 0, num_insts = 0, const_count = 0, LHS = -1, endAddr = startAddr + numConsts;
  int *inst_count = (int *)bmalloc(numConsts * sizeof(int));

  while (const_count < numConsts)
  { // increment the number of instructions we have seen
    num_insts++;
    // read in the next instruction to tmp_insts
    tmp_insts[k++] = prog->prog[j++]; // operation
    LHS = tmp_insts[k++] = prog->prog[j++]; // LHS

    // check to see if LHS is a constant location
    if (startAddr <= LHS && LHS < endAddr)
    { // we have found the next constant - save the number of instructions for this one
      inst_count[const_count] = num_insts;
      // reset num_insts to 0 and increment const_count
      num_insts = 0;
      const_count++;
    }

    // read in the operands
    tmp_insts[k++] = prog->prog[j++];
    if (!isUnary(tmp_insts[k-3]))
    { // binary operation so we read in the second operand
      tmp_insts[k++] = prog->prog[j++];
    }
  }

  // save where the index is
  *prog_index = j;

  // now that we have the update instructions, we can setup insts
  k = 0;
  for (i = 0; i < numConsts; i++)
  { // setup the ith set of update instructions
    insts[i].num = inst_count[i];
    insts[i].insts = (inst *)bmalloc(inst_count[i] * sizeof(inst));

    // setup these instructions
    for (j = 0; j < inst_count[i]; j++)
    {
      insts[i].insts[j].op = tmp_insts[k++];
      insts[i].insts[j].memLoc = tmp_insts[k++];
      insts[i].insts[j].in[0] = tmp_insts[k++];
      if (isUnary(insts[i].insts[j].op)) // unary
        insts[i].insts[j].in[1] = 0;
      else
        insts[i].insts[j].in[1] = tmp_insts[k++];
    }
  }

  // clear memory
  free(inst_count);

  return;
}

void get_some_insts(prog_t *prog, int *tmp_insts, inst_list *insts, int insts_index, int target_addr, int *prog_index)
/****************************************************************************\
* USAGE: grabs all insts for defining one object and stores them in "insts"  *
* ARGUMENTS: the SLP, array for storing the insts, target address, etc.      *
* RETURN VALUES: none                                                        *
* NOTES:                                                                     *
\****************************************************************************/
{
  int i, j, k, Done, num_insts;

  j = *prog_index;  
  Done = 0;  //1 triggers the loop to break
  k = 0;  //index within tmp_insts 
  num_insts = 0;  //total number of insts

  while (!Done)
  {
    num_insts++;
    tmp_insts[k++] = prog->prog[j++];  //Grab the next (jth) instructions operation.
    tmp_insts[k++] = prog->prog[j++];  //Mem location of the LHS.
    tmp_insts[k++] = prog->prog[j++];  //Mem location of the first operand.
    if (!isUnary(tmp_insts[k-3])) //If we have a binary operation, we grab the second operand.
    {
      tmp_insts[k++] = prog->prog[j++];  //Mem location of the second operand.
      if (tmp_insts[k-3] == target_addr) //If we have found the target addr as the LHS, we are done!
	Done = 1;
    }
    else //We must check to see if the target addr is the LHS in this case (unary operation) as well.
    {
      if (tmp_insts[k-2] == target_addr)
	Done = 1;
    }
  }

  //Now we just store tmp_insts in the appropriate entry within insts. 
  insts[insts_index].num = num_insts;
  insts[insts_index].insts = (inst *)bcalloc(num_insts, sizeof(inst));
  k = 0;
  for (i=0; i<num_insts; i++)
  {
    insts[insts_index].insts[i].op = tmp_insts[k++];
    insts[insts_index].insts[i].memLoc = tmp_insts[k++];
    insts[insts_index].insts[i].in[0] = tmp_insts[k++];
    if (isUnary(insts[insts_index].insts[i].op))  // unary op
      insts[insts_index].insts[i].in[1] = 0;
    else
      insts[insts_index].insts[i].in[1] = tmp_insts[k++];
  }

  *prog_index = j;

  return;
}

//////////////////////// ADDING A PATCH EQUATION TO THE BOTTOM OF THE SLP ///////////////////////////////////

void add_vec_patch_SLP(prog_t *new_prog, prog_t *old_prog, mpq_t **patch, int patch_size, mpq_t *patch_rhs)
/************************************************************************\
* USAGE: reads in old_prog, appends patch and stores in new_prog         *
* ARGUMENTS:                                                             *
* RETURN VALUES:                                                         *
* NOTES: This is for 1 variable group - i.e. 1 patch                     *
\************************************************************************/
{ // error checking - make sure that patch_size == numVars
  if (patch_size != old_prog->numVars)
  {
    printf("ERROR: The size of the patch must be the same as the number of variables!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // make sure that there are no path variables
  if (old_prog->numPathVars > 0)
  {
    printf("ERROR: Path variables are not allowed when adding a patch to the straight-line program.!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }
  // make sure that there are no parameters
  if (old_prog->numPars > 0)
  {
    printf("ERROR: Parameters are not allowed when adding a patch to the straight-line program.!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  int i, j, k, new_workspace_addr, old_workspace_addr, patch_addr, room_before_temps, patch_num_index, tmp_func_addr;
  int *addr_conv_table = (int *)bmalloc(old_prog->memSize * sizeof(int));
  char ch;
  FILE *FOUT = NULL;

  int old_const_size = 0, old_subfunc_size = 0, old_func_size = 0, old_funcjac_size = 0, old_subfuncjac_size = 0;
  inst_list *old_const_insts = NULL, *old_subfunc_insts = NULL, *old_func_insts = NULL, *old_funcjac_insts = NULL, *old_subfuncjac_insts = NULL;
  int new_const_size = 0, new_subfunc_size = 0, new_func_size = 0;
  inst_list *new_const_insts = NULL, *new_subfunc_insts = NULL, *new_func_insts = NULL;

  // setup the structures
  setup_insts_from_prog(&old_const_insts, &old_const_size, &old_subfunc_insts, &old_subfunc_size, &old_func_insts, &old_func_size, &old_funcjac_insts, &old_funcjac_size, &old_subfuncjac_insts, &old_subfuncjac_size, old_prog);

  // all we are doing is adding some new 'nums', a new function and its derivatives
  // For future reference, the order of the instructions is constants, parameters, subfunctions, functions, subfunc derivs, func derivs.

  new_prog->precision = old_prog->precision; // The precision at which evaluation should occur
  new_prog->num_var_gps = old_prog->num_var_gps;
  new_prog->var_gp_sizes = (int *)bmalloc(new_prog->num_var_gps * sizeof(int));
  for (i = 0; i < old_prog->num_var_gps; i++)
    new_prog->var_gp_sizes[i] = old_prog->var_gp_sizes[i];
  new_prog->index_of_first_number_for_proj_trans = old_prog->index_of_first_number_for_proj_trans;  //COULD BE WRONG - IS THIS USED???!!!
  // These next few aren't being used yet, so default to old ones for now...needs to be changed later if used!!!
  new_prog->numInstAtEndUpdate = old_prog->numInstAtEndUpdate;
  new_prog->numInstAtEndParams = old_prog->numInstAtEndParams;
  new_prog->numInstAtEndFnEval = old_prog->numInstAtEndFnEval;
  new_prog->numInstAtEndPDeriv = old_prog->numInstAtEndPDeriv;
  new_prog->numInstAtEndJvEval = old_prog->numInstAtEndJvEval;

  // now for the numbers and locations of various types
  new_workspace_addr = old_workspace_addr = 0;  // We'll increment these as we go through this block to get to the new and old first workspace addrs.

  /////////////////////////VARS - no changes
  new_prog->numVars = old_prog->numVars;
  new_prog->inpVars = old_prog->inpVars; //should just be 0 since vars go first.
  new_workspace_addr += new_prog->numVars;
  old_workspace_addr += old_prog->numVars;

  for (i = 0; i < old_prog->numVars; i++)
    addr_conv_table[old_prog->inpVars + i] = new_prog->inpVars + i;

  /////////////////////////PATHVARS - no changes
  new_prog->numPathVars = old_prog->numPathVars;
  new_prog->inpPathVars = new_prog->inpVars + new_prog->numVars;
  new_workspace_addr += new_prog->numPathVars;
  old_workspace_addr += old_prog->numPathVars;

  for (i = 0; i < old_prog->numPathVars; i++) // shouldn't be any path variables
    addr_conv_table[old_prog->inpPathVars + i] = new_prog->inpPathVars + i;

  /////////////////////////PARAMS - no changes
  new_prog->numPars = old_prog->numPars;
  new_prog->evalPars = new_prog->inpPathVars + new_prog->numPathVars;
  new_workspace_addr += new_prog->numPars;
  old_workspace_addr += old_prog->numPars;

  for (i = 0; i < old_prog->numPars; i++) // shouldn't be any parameters
    addr_conv_table[old_prog->evalPars + i] = new_prog->evalPars + i;

  /////////////////////////FUNCS - add 1 functions
  new_prog->numFuncs = old_prog->numFuncs + 1; // original functions plus the patch
  new_prog->evalFuncs = new_prog->evalPars + new_prog->numPars;
  new_workspace_addr += new_prog->numFuncs;
  old_workspace_addr += old_prog->numFuncs;

  for (i = 0; i < old_prog->numFuncs; i++)
    addr_conv_table[old_prog->evalFuncs + i] = new_prog->evalFuncs + i;

  /////////////////////////DPARAMS
  new_prog->evalDPars = new_prog->evalFuncs + new_prog->numFuncs;
  new_workspace_addr += new_prog->numPars;
  old_workspace_addr += old_prog->numPars;

  for (i = 0; i < old_prog->numPars; i++) // shouldn't be any parameters
    addr_conv_table[old_prog->evalDPars + i] = new_prog->evalDPars + i;

  /////////////////////////JFUNCS(V)
  new_prog->evalJVars = new_prog->evalDPars + new_prog->numPars;
  new_workspace_addr += new_prog->numFuncs * new_prog->numVars;
  old_workspace_addr += old_prog->numFuncs * old_prog->numVars;

  for (i = 0; i < old_prog->numVars * old_prog->numFuncs; i++)
    addr_conv_table[old_prog->evalJVars + i] = new_prog->evalJVars + i;

  /////////////////////////JFUNCS(P)
  new_prog->evalJPars = new_prog->evalJVars + new_prog->numFuncs * new_prog->numVars;
  new_workspace_addr += new_prog->numFuncs * new_prog->numPars;
  old_workspace_addr += old_prog->numFuncs * old_prog->numPars;

  for (i = 0; i < old_prog->numPars; i++) // shouldn't be any parameters
    addr_conv_table[old_prog->evalJPars + i] = new_prog->evalJPars + i;

  /////////////////////////CONSTANTS - add the number of constants
  new_prog->numConsts = old_prog->numConsts + patch_size + 1;
  new_prog->constAddr = new_prog->evalJPars + new_prog->numFuncs * new_prog->numPars;
  new_prog->IAddr = new_prog->constAddr;  //since I is always the first constant
  patch_addr = new_prog->constAddr + old_prog->numConsts;  // Keeps track of first addr in the patch
  new_workspace_addr += new_prog->numConsts;
  old_workspace_addr += old_prog->numConsts;

  for (i = 0; i < old_prog->numConsts; i++)
    addr_conv_table[old_prog->constAddr + i] = new_prog->constAddr + i;

  /////////////////////////NUMS - add the number of new 'nums'
  new_prog->numNums = old_prog->numNums + 2 * (patch_size + 1);
  new_prog->numAddr = new_prog->constAddr + new_prog->numConsts;
  new_workspace_addr += new_prog->numNums;
  old_workspace_addr += old_prog->numNums;

  for (i = 0; i < old_prog->numNums; i++)
    addr_conv_table[old_prog->numAddr + i] = new_prog->numAddr + i;

  /////////////////////////SUBFUNCS - no changes
  new_prog->numSubfuncs = old_prog->numSubfuncs;
  new_prog->evalSubs = new_prog->numAddr + new_prog->numNums;
  new_workspace_addr += new_prog->numSubfuncs;
  old_workspace_addr += old_prog->numSubfuncs;

  for (i = 0; i < old_prog->numSubfuncs; i++)
    addr_conv_table[old_prog->evalSubs + i] = new_prog->evalSubs + i;

  /////////////////////////JSUBS(V)
  new_prog->evalJSubsV = new_prog->evalSubs + new_prog->numSubfuncs;
  new_workspace_addr += new_prog->numSubfuncs * new_prog->numVars;
  old_workspace_addr += old_prog->numSubfuncs * old_prog->numVars;

  for (i = 0; i < old_prog->numSubfuncs * old_prog->numVars; i++)
    addr_conv_table[old_prog->evalJSubsV + i] = new_prog->evalJSubsV + i;

  /////////////////////////JSUBS(P)
  new_prog->evalJSubsP = new_prog->evalJSubsV + new_prog->numSubfuncs * new_prog->numVars;
  new_workspace_addr += new_prog->numSubfuncs * new_prog->numPars;
  old_workspace_addr += old_prog->numSubfuncs * old_prog->numPars;

  for (i = 0; i < old_prog->numSubfuncs * old_prog->numPars; i++) // shouldn't be any parameters
    addr_conv_table[old_prog->evalJSubsP + i] = new_prog->evalJSubsP + i;

  // remember the beginning of the workspace:
  room_before_temps = new_workspace_addr;

  // the rest of the old workspace was set aside as temp memLocs
  for (i = old_workspace_addr; i < old_prog->memSize; i++)
    addr_conv_table[i] = new_workspace_addr++;

////////////////////// CONSTRUCT THE NEW SLP

  // setup new_prog->nums
  new_prog->numNums = old_prog->numNums + 2 * (patch_size + 1);
  new_prog->nums = (num_t *)bmalloc(new_prog->numNums * sizeof(num_t));
  for (i = 0; i < old_prog->numNums; i++)
  { // copy in the old numbers
    mpq_init(new_prog->nums[i].rat);
    mpq_set(new_prog->nums[i].rat, old_prog->nums[i].rat);

    new_prog->nums[i].currPrec = new_prog->precision;
    mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
    mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);
  }

  // add in the new numbers for the patch
  patch_num_index = i; // remember where the numbers for the patch start
  for (j = 0; j <= patch_size; j++)
    if (j < patch_size)
    { // real part of patch
      mpq_init(new_prog->nums[i].rat);
      mpq_set(new_prog->nums[i].rat, patch[j][0]);

      new_prog->nums[i].currPrec = new_prog->precision;
      mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
      mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);

      // increment i
      i++;

      // imag part of patch
      mpq_init(new_prog->nums[i].rat);
      mpq_set(new_prog->nums[i].rat, patch[j][1]);

      new_prog->nums[i].currPrec = new_prog->precision;
      mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
      mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);

      // increment i
      i++;
    }
    else
    { // real part of patch_rhs
      mpq_init(new_prog->nums[i].rat);
      mpq_set(new_prog->nums[i].rat, patch_rhs[0]);

      new_prog->nums[i].currPrec = new_prog->precision;
      mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
      mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);

      // increment i
      i++;

      // imag part of patch
      mpq_init(new_prog->nums[i].rat);
      mpq_set(new_prog->nums[i].rat, patch_rhs[1]);

      new_prog->nums[i].currPrec = new_prog->precision;
      mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
      mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);

      // increment i
      i++;
    }

///////////////////////////////////////////
// ADDING CONST EVAL INSTS TO NEW SLP
///////////////////////////////////////////

  // First is the constants - prepare new_const_insts
  new_const_size = new_prog->numConsts - 2; // -2 b/c I & Pi is #0 & #1 (no update insts)
  new_const_insts = (inst_list *)bmalloc(new_const_size * sizeof(inst_list));

  // copy over all of the old constant evaluations
  for (i = 0; i < old_const_size; i++)
  {
    new_const_insts[i].num = old_const_insts[i].num;
    new_const_insts[i].insts = (inst *)bmalloc(new_const_insts[i].num * sizeof(inst));
    for (j = 0; j < new_const_insts[i].num; j++)
    { // copy over the instructions
      new_const_insts[i].insts[j].op = old_const_insts[i].insts[j].op;
      new_const_insts[i].insts[j].memLoc = addr_conv_table[old_const_insts[i].insts[j].memLoc];
      new_const_insts[i].insts[j].in[0] = addr_conv_table[old_const_insts[i].insts[j].in[0]];
      new_const_insts[i].insts[j].in[1] = addr_conv_table[old_const_insts[i].insts[j].in[1]];
    }
  }

  // setup the constant evaluations for the patch coefficients and the right hand side of the patch
  k = patch_num_index; // where the numbers for the patch start
  for (i = old_const_size; i < new_const_size; i++)
  {
    new_const_insts[i].num = 2; // each new constant creates 2 instructions (1) imag_part * I & (2) real_part + imag_part * I
    new_const_insts[i].insts = (inst *)bmalloc(new_const_insts[i].num * sizeof(inst));

    // imag_part * I
    new_const_insts[i].insts[0].op = '*'; // multiply
    new_const_insts[i].insts[0].memLoc = new_workspace_addr; 
    new_const_insts[i].insts[0].in[0] = new_prog->numAddr + k + 1; // for imaginary part of the number
    new_const_insts[i].insts[0].in[1] = new_prog->IAddr;

    // real_part + imag_part * I
    new_const_insts[i].insts[1].op = '+'; // add
    new_const_insts[i].insts[1].memLoc = new_prog->constAddr + i + 2; // store result in appropriate addr (+2 b/c of I & Pi)
    new_const_insts[i].insts[1].in[0] = new_workspace_addr;
    new_const_insts[i].insts[1].in[1] = new_prog->numAddr + k; // for real part of the number

    // increment k
    k += 2;
  }
  new_workspace_addr++;

///////////////////////////////////////////
// COPYING SUBFUNC INSTS TO NEW SLP
///////////////////////////////////////////

  // Next is the subfunctions - prepare new_subfunc_insts
  new_subfunc_size = old_subfunc_size;
  new_subfunc_insts = (inst_list *)bmalloc(new_subfunc_size * sizeof(inst_list));

  // copy over all of the old subfunction evaluations
  for (i = 0; i < old_subfunc_size; i++)
  {
    new_subfunc_insts[i].num = old_subfunc_insts[i].num;
    new_subfunc_insts[i].insts = (inst *)bmalloc(new_subfunc_insts[i].num * sizeof(inst));
    for (j = 0; j < new_subfunc_insts[i].num; j++)
    { // copy over the instructions
      new_subfunc_insts[i].insts[j].op = old_subfunc_insts[i].insts[j].op;
      new_subfunc_insts[i].insts[j].memLoc = addr_conv_table[old_subfunc_insts[i].insts[j].memLoc];
      new_subfunc_insts[i].insts[j].in[0] = addr_conv_table[old_subfunc_insts[i].insts[j].in[0]];
      new_subfunc_insts[i].insts[j].in[1] = addr_conv_table[old_subfunc_insts[i].insts[j].in[1]];
    }
  }

///////////////////////////////////////////
// ADDING FUNC INSTS TO NEW SLP
///////////////////////////////////////////

  // Last is the functions - prepare new_func_insts
  new_func_size = old_func_size + 1;
  new_func_insts = (inst_list *)bmalloc(new_func_size * sizeof(inst_list));

  // copy over all of the old function evaluations
  for (i = 0; i < old_func_size; i++)
  {
    new_func_insts[i].num = old_func_insts[i].num;
    new_func_insts[i].insts = (inst *)bmalloc(new_func_insts[i].num * sizeof(inst));
    for (j = 0; j < new_func_insts[i].num; j++)
    { // copy over the instructions
      new_func_insts[i].insts[j].op = old_func_insts[i].insts[j].op;
      new_func_insts[i].insts[j].memLoc = addr_conv_table[old_func_insts[i].insts[j].memLoc];
      new_func_insts[i].insts[j].in[0] = addr_conv_table[old_func_insts[i].insts[j].in[0]];
      new_func_insts[i].insts[j].in[1] = addr_conv_table[old_func_insts[i].insts[j].in[1]];
    }
  }

  // add instructions for evaluating the patch
  i = old_func_size;
  new_func_insts[i].num = 2 * patch_size + 2;
  new_func_insts[i].insts = (inst *)bmalloc(new_func_insts[i].num * sizeof(inst));

  tmp_func_addr = new_workspace_addr++;
  k = 0; // instruction number
  // initialize to the negative of the right-hand size
  new_func_insts[i].insts[k].op = 'N'; // negate
  new_func_insts[i].insts[k].memLoc = tmp_func_addr;
  new_func_insts[i].insts[k].in[0] = patch_addr + patch_size; // right-hand side
  new_func_insts[i].insts[k].in[1] = 0; // unary negate
  // increment k
  k++;

  for (j = 0; j < patch_size; j++)
  { // do sum += patch_coeff[j] * vars[j]
    new_func_insts[i].insts[k].op = '*'; // multiply
    new_func_insts[i].insts[k].memLoc = new_workspace_addr;
    new_func_insts[i].insts[k].in[0] = patch_addr + j; // jth patch coefficient
    new_func_insts[i].insts[k].in[1] = new_prog->inpVars + j; // jth variable
    // increment
    new_workspace_addr++;
    k++;

    // add to the sum
    new_func_insts[i].insts[k].op = '+'; // addition
    new_func_insts[i].insts[k].memLoc = new_workspace_addr;
    new_func_insts[i].insts[k].in[0] = tmp_func_addr;
    new_func_insts[i].insts[k].in[1] = new_workspace_addr - 1;
    // update & increment
    tmp_func_addr = new_workspace_addr;
    new_workspace_addr++;
    k++;
  }
  // store to proper location
  new_func_insts[i].insts[k].op = '='; // store
  new_func_insts[i].insts[k].memLoc = new_prog->evalFuncs + old_prog->numFuncs;
  new_func_insts[i].insts[k].in[0] = tmp_func_addr;
  new_func_insts[i].insts[k].in[1] = 0;

///////////////////////////////////////////////////
//  COPYING ALL NEW INSTS TO FILE (FOR DIFF)
///////////////////////////////////////////////////

  // delete the old files
  remove("finalFile.out");
  remove("jacV.out");
  remove("jacP.out");
  remove("arr.out");

  // open up finalFile.out
  FOUT = fopen("finalFile.out", "w");
  // write the instructions to FOUT
  write_insts_to_file(FOUT, new_prog, new_const_insts, new_const_size, NULL, 0, new_subfunc_insts, new_subfunc_size, new_func_insts, new_func_size, room_before_temps, new_workspace_addr);
  // close FOUT
  fclose(FOUT);

////////////////////////////////////////////////
//  CALLING DIFF WITH NEW SLP
////////////////////////////////////////////////
  diff_deflatable(1, 0);

////////////////////////////////////////////////
//  COPYING COMPLETE SLP BACK INTO STRUCT
////////////////////////////////////////////////

  // read in the differentiated instructions
  FOUT = fopen("arr.out", "r");
  if (FOUT == NULL)
  {
    printf("ERROR: The file containing the straight-line program does not exist!\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // get the memSize - move through arr.out until we find the first 'd' in TotalRoomNeeded
  do
  {
    ch = fgetc(FOUT);
  } while (ch != 'd');
  // read in the size
  fscanf(FOUT, "ed %d", &new_prog->memSize);

  // read in the number of instructions - move through arr.out until we find the 'I' in NUMINST
  do
  {
    ch = fgetc(FOUT);
  } while (ch != 'I');
  fscanf(FOUT, "NST %d", &new_prog->size);

  // grab all of the instructions.
  rewind(FOUT);
  new_prog->prog = (int *)bcalloc(new_prog->size, sizeof(int));
  for (i = 0; i < new_prog->size; i++)
    fscanf(FOUT, "%d ", &new_prog->prog[i]);

  // close arr.out
  fclose(FOUT);

  // clear the memory
  free(addr_conv_table);
  clear_inst_list(old_const_insts, old_const_size);
  free(old_const_insts);
  clear_inst_list(old_subfunc_insts, old_subfunc_size);
  free(old_subfunc_insts);
  clear_inst_list(old_func_insts, old_func_size);
  free(old_func_insts);
  clear_inst_list(old_funcjac_insts, old_funcjac_size);
  free(old_funcjac_insts);
  clear_inst_list(old_subfuncjac_insts, old_subfuncjac_size);
  free(old_subfuncjac_insts);

  clear_inst_list(new_const_insts, new_const_size);
  free(new_const_insts);
  clear_inst_list(new_subfunc_insts, new_subfunc_size);
  free(new_subfunc_insts);
  clear_inst_list(new_func_insts, new_func_size);
  free(new_func_insts);

  return;
}

void write_insts_to_file(FILE *FOUT, prog_t *prog, inst_list *const_insts, int const_size, inst_list *param_insts, int param_size, inst_list *subfunc_insts, int subfunc_size, inst_list *func_insts, int func_size, int room_before_temps, int total_room_needed)
/************************************************************************\
* USAGE: writes the instructions to the FOUT                             *
* ARGUMENTS:                                                             *
* RETURN VALUES:                                                         *
* NOTES:                                                                 *
\************************************************************************/
{
  int i;

  // preamble
  fprintf(FOUT, "RoomBeforeTemps %d;\nTotalRoomNeeded %d;\n", room_before_temps, total_room_needed);
  fprintf(FOUT, "NVAR %d %d;\n", prog->numVars, prog->inpVars);
  fprintf(FOUT, "NPATHVAR %d %d;\n", prog->numPathVars, prog->inpPathVars);
  fprintf(FOUT, "NPAR %d %d %d;\n", prog->numPars, prog->evalPars, prog->evalDPars);
  fprintf(FOUT, "NFCN %d %d %d %d;\n", prog->numFuncs, prog->evalFuncs, prog->evalJVars, prog->evalJPars);
  fprintf(FOUT, "NCON %d %d;\n", prog->numConsts, prog->constAddr);
  fprintf(FOUT, "NNUM %d %d;\n", prog->numNums, prog->numAddr);
  fprintf(FOUT, "SUBFCN %d %d %d %d;\n", prog->numSubfuncs, prog->evalSubs, prog->evalJSubsV, prog->evalJSubsP);
  fprintf(FOUT, "CMPLX %d 0 1;\n", prog->constAddr); //I is the first constant.
  fprintf(FOUT, "ZERO %d;\n", prog->numAddr);  //0 is the first num.
  fprintf(FOUT, "ONE %d;\n", prog->numAddr+1); //1 is the second num.
  fprintf(FOUT, "VARGPS %d;\n", prog->num_var_gps);
  for (i = 0; i < prog->num_var_gps; i++)
    fprintf(FOUT, " %d", prog->var_gp_sizes[i]);
  fprintf(FOUT, ";\n");
  fprintf(FOUT, "RANDINDEX %d;\n",  prog->index_of_first_number_for_proj_trans);

  // constant instructions
  fprintf(FOUT, "BEGIN UPDATE;\n");
  print_insts(FOUT, const_insts, const_size);

  // parameter instructions
  fprintf(FOUT, "BEGIN PARAM;\n");
  print_insts(FOUT, param_insts, param_size);

  // subfunction instructions
  fprintf(FOUT, "BEGIN FUNCTION;\n");
  print_insts(FOUT, subfunc_insts, subfunc_size);
  print_insts(FOUT, func_insts, func_size);

  // the end
  fprintf(FOUT, "END;\n");

  return;
}

void print_insts(FILE *FOUT, inst_list *list, int num_list)
/************************************************************************\
* USAGE: prints the instructions in list to FOUT                         *
* ARGUMENTS:                                                             *
* RETURN VALUES:                                                         *
* NOTES:                                                                 *
\************************************************************************/
{
  int i, j;

  for (i = 0; i < num_list; i++)
    for (j = 0; j < list[i].num; j++)
      if (isUnary(list[i].insts[j].op))
        fprintf(FOUT, "%c %d %d;\n", list[i].insts[j].op, list[i].insts[j].memLoc, list[i].insts[j].in[0]);
      else
        fprintf(FOUT, "%c %d %d %d;\n", list[i].insts[j].op, list[i].insts[j].memLoc, list[i].insts[j].in[0], list[i].insts[j].in[1]);

  return;
}

void clear_inst_list(inst_list *list, int num_list)
/************************************************************************\
* USAGE: clears the memory from list                                     *
* ARGUMENTS:                                                             *
* RETURN VALUES:                                                         *
* NOTES:                                                                 *
\************************************************************************/
{
  int i;

  for (i = num_list - 1; i >= 0; i--)
    free(list[i].insts);

  return;
}

void setup_insts_from_prog(inst_list **const_insts, int *const_size, inst_list **subfunc_insts, int *subfunc_size, inst_list **func_insts, int *func_size, inst_list **funcjac_insts, int *funcjac_size, inst_list **subfuncjac_insts, int *subfuncjac_size, prog_t *prog)
/************************************************************************\
* USAGE: reads in prog, and sets up the instructions                     *
* ARGUMENTS:                                                             *
* RETURN VALUES:                                                         *
* NOTES: allocates the datastructures as needed                          *
\************************************************************************/
{
  int i, j, k;
  int *tmp_insts = (int *)bmalloc(prog->size * sizeof(int)); // too big (obviously) for one func, but oh well

  // allocate the structures needed to read in orig_prog
  *const_size = prog->numConsts - 2; // -2 b/c I & Pi is #0 & #1 (no update insts)
  *const_insts = (inst_list *)bmalloc((*const_size) *  sizeof(inst_list));

  *subfunc_size = prog->numSubfuncs;
  *subfunc_insts = (inst_list *)bmalloc((*subfunc_size) * sizeof(inst_list));

  *func_size = prog->numFuncs;
  *func_insts = (inst_list *)bmalloc((*func_size) * sizeof(inst_list));

  *funcjac_size = prog->numFuncs * prog->numVars;
  *funcjac_insts = (inst_list *)bmalloc((*funcjac_size) * sizeof(inst_list));

  *subfuncjac_size = prog->numSubfuncs * prog->numVars;
  *subfuncjac_insts = (inst_list *)bmalloc((*subfuncjac_size) * sizeof(inst_list));

  // We go through the SLP in prog and grab insts one at a time until the first operand is the location of the target.
  // We store all insts in tmp_insts until done.
  // Then we know the number of them and can plug them into the appropriate inst_list structure.
  // For future reference, the order of the instructions is constants, parameters, subfunctions, functions, subfunc derivs, func derivs.
  // NOTE:  This code assumes no parameters and no path variable (since we are not implementing deflation for user homotopies).
  // The wkspace array goes vars, pathvars, pars, funcs, param ders, J(func/vars), J(func/pars), consts, nums, subfuncs, J(subfunc/vars), J(subfunc/pars).
  // Finally, I (sqrt(-1)) is the first const, 0 is the first num, and 1 is the second num.
  j = 0;  //counter in the SLP
  // Grab the constant eval insts - since constants do not have to be written in order, we have a more elaborate method to read in these instructions
  get_constant_ints(prog, tmp_insts, *const_insts, prog->constAddr + 2, *const_size, &j); // +2 b/c I & Pi are the 2 constants and there are not update instructions for it
 
  // Grab the subfunc eval insts
  for (i = 0; i < *subfunc_size; i++)
  {
    get_some_insts(prog, tmp_insts, *subfunc_insts, i, prog->evalSubs + i, &j);
  }

  // Grab the func eval insts
  for (i = 0; i < *func_size; i++)
  {
    get_some_insts(prog, tmp_insts, *func_insts, i, prog->evalFuncs + i, &j);
  }

  // Now for the subfunc derivative insts (assuming forward AD was used)
  for (i = 0; i < prog->numSubfuncs; i++)
    for (k = 0; k < prog->numVars; k++)  //NOTE: We use k here because j is special.
    {
      get_some_insts(prog, tmp_insts, *subfuncjac_insts, i*prog->numVars + k, prog->evalJSubsV + i*prog->numVars + k, &j);
    }

  // Finally, we grab the only jacobian eval insts (again, assuming forward AD was used to create them!).
  for (i = 0; i < prog->numFuncs; i++)
    for (k = 0; k < prog->numVars; k++)  //NOTE: We use k here because j is special.
    {
      get_some_insts(prog, tmp_insts, *funcjac_insts, i*prog->numVars + k, prog->evalJVars + i*prog->numVars + k, &j);
    }

  // clear local memory
  free(tmp_insts);

  return;
}

void extend_deflation_point(point_d output_pt_d, point_mp output_pt_mp, int *output_prec, point_d input_pt_d, point_mp input_pt_mp, mat_d J_d, mat_mp J_mp, mpq_t ***K, int K_rows, int K_cols, mpq_t ***B, int B_rows, int B_cols, int input_prec)
/************************************************************************\
* USAGE: finds the 'deflation' variables for the given point             *
* ARGUMENTS:                                                             *
* RETURN VALUES:                                                         *
* NOTES:                                                                 *
\************************************************************************/
{
  int i, j;

  // this does not change the precision
  *output_prec = input_prec;

  if (input_prec < 64)
  { // do everything in double precision
    mat_d K_d, B_d;
    // setup K_d
    init_mat_d(K_d, K_rows, K_cols);
    K_d->rows = K_rows;
    K_d->cols = K_cols;
    for (i = 0; i < K_rows; i++)
      for (j = 0; j < K_cols; j++)
      { // setup K_d(i,j)
        K_d->entry[i][j].r = mpq_get_d(K[i][j][0]);
        K_d->entry[i][j].i = mpq_get_d(K[i][j][1]);
      }

    // setup B_d
    init_mat_d(B_d, B_rows, B_cols);
    B_d->rows = B_rows;
    B_d->cols = B_cols;
    for (i = 0; i < B_rows; i++)
      for (j = 0; j < B_cols; j++)
      { // setup B_d(i,j)
        B_d->entry[i][j].r = mpq_get_d(B[i][j][0]);
        B_d->entry[i][j].i = mpq_get_d(B[i][j][1]);
      }

    // extend the point
    extend_deflation_point_d(output_pt_d, input_pt_d, J_d, K_d, B_d);

    // clear
    clear_mat_d(K_d); clear_mat_d(B_d);
  }
  else
  { // do everything in the given precision
    mat_mp K_mp, B_mp;

    // set the precision
    initMP(input_prec);

    // setup K_mp
    init_mat_mp2(K_mp, K_rows, K_cols, input_prec);
    K_mp->rows = K_rows;
    K_mp->cols = K_cols;
    for (i = 0; i < K_rows; i++)
      for (j = 0; j < K_cols; j++)
      { // setup K_mp(i,j)
        mpf_set_q(K_mp->entry[i][j].r, K[i][j][0]);
        mpf_set_q(K_mp->entry[i][j].i, K[i][j][1]);
      }

    // setup B_mp
    init_mat_mp2(B_mp, B_rows, B_cols, input_prec);
    B_mp->rows = B_rows;
    B_mp->cols = B_cols;
    for (i = 0; i < B_rows; i++)
      for (j = 0; j < B_cols; j++)
      { // setup B_mp(i,j)
        mpf_set_q(B_mp->entry[i][j].r, B[i][j][0]);
        mpf_set_q(B_mp->entry[i][j].i, B[i][j][1]);
      }

    // extend the point
    extend_deflation_point_mp(output_pt_mp, input_pt_mp, J_mp, K_mp, B_mp);

    // clear
    clear_mat_mp(K_mp); clear_mat_mp(B_mp);
  }

  return;
}

void extend_deflation_point_d(point_d output_pt, point_d input_pt, mat_d Jv, mat_d K, mat_d B)
/************************************************************************\
* USAGE: finds the 'deflation' variables for the given point using _d    *
* ARGUMENTS:                                                             *
* RETURN VALUES:                                                         *
* NOTES:                                                                 *
\************************************************************************/
{
  int i, j, retVal, size;
  vec_d b;
  mat_d A, R;

  init_vec_d(b, 0);
  init_mat_d(A, 0, 0); init_mat_d(R, 0, 0);

  // setup A = [[Jv*K];[B];]
  mat_mul_d(A, Jv, K); // setup the top of A

  // add rows to bottom of A
  size = A->rows;
  i = A->rows + B->rows;
  increase_rows_mat_d(A, i);
  A->rows = i;
  // setup the bottom of A
  for (i = 0; i < B->rows; i++)
    for (j = 0; j < A->cols; j++)
    {
      set_d(&A->entry[i + size][j], &B->entry[i][j]);
    }

  // setup b = [[0];[1];]
  increase_size_vec_d(b, A->rows);
  b->size = A->rows;
  for (i = 0; i < b->size; i++)
    if (i < Jv->rows)
    { // top is 0
      set_zero_d(&b->coord[i]);
    }
    else
    { // bottom is 1
      set_one_d(&b->coord[i]);
    }

  // setup R - random matrix
  make_matrix_random_d(R, A->cols, A->rows);

  // find R*A & R*b
  mat_mul_d(A, R, A);
  mul_mat_vec_d(b, R, b);

  // solve for the new variables
  retVal = matrixSolve_d(b, A, b);
  if (retVal)
  { // this should never happen
    printf("ERROR: There was a problem extending the variables for deflation!\n");
    bexit(ERROR_OTHER);
  }

  // setup ouptut_pt
  size = input_pt->size;
  increase_size_vec_d(output_pt, size + b->size);
  output_pt->size = size + b->size;
  for (i = 0; i < output_pt->size; i++)
    if (i < size)
    { // leave alone
      set_d(&output_pt->coord[i], &input_pt->coord[i]);
    }
    else
    { // copy in new coordinate
      set_d(&output_pt->coord[i], &b->coord[i - size]);
    }

  clear_vec_d(b);
  clear_mat_d(A); clear_mat_d(R);

  return;
}

void extend_deflation_point_mp(point_mp output_pt, point_mp input_pt, mat_mp Jv, mat_mp K, mat_mp B)
/************************************************************************\
* USAGE: finds the 'deflation' variables for the given point using _mp   *
* ARGUMENTS:                                                             *
* RETURN VALUES:                                                         *
* NOTES:                                                                 *
\************************************************************************/
{
  int i, j, size, retVal;
  vec_mp b;
  mat_mp A, R;

  init_vec_mp(b, 0);
  init_mat_mp(A, 0, 0); init_mat_mp(R, 0, 0);

  // setup A = [[Jv*K];[B];]
  mat_mul_mp(A, Jv, K); // setup the top of A
  // add rows to the bottom of A
  size = A->rows;
  i = A->rows + B->rows;
  increase_rows_mat_mp(A, i);
  A->rows = i;
  // setup the bottom of A
  for (i = 0; i < B->rows; i++)
    for (j = 0; j < A->cols; j++)
    {
      set_mp(&A->entry[i + size][j], &B->entry[i][j]);
    }

  // setup b = [[0];[1];]
  increase_size_vec_mp(b, A->rows);
  b->size = A->rows;
  for (i = 0; i < b->size; i++)
    if (i < Jv->rows)
    { // top is 0
      set_zero_mp(&b->coord[i]);
    }
    else
    { // bottom is 1
      set_one_mp(&b->coord[i]);
    }

  // setup R - random matrix
  make_matrix_random_mp(R, A->cols, A->rows, mpf_get_default_prec());

  // find R*A & R*b
  mat_mul_mp(A, R, A);
  mul_mat_vec_mp(b, R, b);

  // solve for the new variables
  retVal = matrixSolve_mp(b, A, b);
  if (retVal)
  { // this should never happen
    printf("ERROR: There was a problem extending the variables for deflation!\n");
    bexit(ERROR_OTHER);
  }

  // setup ouptut_pt
  size = input_pt->size;
  increase_size_vec_mp(output_pt, size + b->size);
  output_pt->size = size + b->size;
  for (i = 0; i < output_pt->size; i++)
    if (i < size)
    { // leave alone
      set_mp(&output_pt->coord[i], &input_pt->coord[i]);
    }
    else
    { // copy in new coordinate
      set_mp(&output_pt->coord[i], &b->coord[i - size]);
    }

  clear_vec_mp(b);
  clear_mat_mp(A); clear_mat_mp(R);

  return;
}

void extend_deflation_matrix_K(mpq_t ****K_rat, int *K_rows, int *K_cols)
/************************************************************************\
* USAGE: extends the 'deflation' matrix K                                *
* ARGUMENTS:                                                             *
* RETURN VALUES:                                                         *
* NOTES: new_rows = rows + cols, new_cols = cols + cols                  *
*   of the form [[K,0];[0,I]]                                            *
\************************************************************************/
{
  int i, j, K_orig_rows = *K_rows, K_orig_cols = *K_cols;
  mpq_t ***K_orig = *K_rat;

  // setup the new size
  *K_rows = K_orig_rows + K_orig_cols;
  *K_cols = K_orig_cols + K_orig_cols;

  // initialize the new K_rat
  init_mat_rat(*K_rat, *K_rows, *K_cols);

  // setup K_rat = [[K_orig,0][0,I]]
  for (i = 0; i < *K_rows; i++)
    if (i < K_orig_rows)
    { // row is of the form [K_orig[i], 0]
      for (j = 0; j < *K_cols; j++)
        if (j < K_orig_cols)
        { // copy K_orig[i][j]
          set_rat((*K_rat)[i][j], K_orig[i][j]);
        }
        else
        { // set to 0
          set_zero_rat((*K_rat)[i][j]);
        }
    }
    else
    { // row is of the form [0,I]
      for (j = 0; j < *K_cols; j++)
        if (j >= K_orig_cols && i - K_orig_rows == j - K_orig_cols)
        { // set to 1
          set_one_rat((*K_rat)[i][j]);
        }
        else
        { // set to 0
          set_zero_rat((*K_rat)[i][j]);
        }
    }

  // clear the old memory
  clear_mat_rat(K_orig, K_orig_rows, K_orig_cols);

  return;
}

int corank_deflation(double *smallest_nonzero_SV, double *largest_zero_SV, point_d point0_d, point_mp point0_mp, int prec0, point_d point1_d, point_mp point1_mp, int prec1, comp_d time_d, comp_mp time_mp, eval_struct_d *e_d, eval_struct_mp *e_mp, tracker_config_t *T, FILE *OUT, prog_t *Prog)
/************************************************************************\
* USAGE: uses point0 & point1 to determine the corank                    *
* ARGUMENTS:                                                             *
* RETURN VALUES: corank                                                  *
* NOTES:                                                                 *
\************************************************************************/
{
  int corank = 0;
  double CN;

  if (T->MPType == 0)
  { // use only double precision
    corank = Cauchy_corank_d(&CN, smallest_nonzero_SV, largest_zero_SV, point0_d, point1_d, time_d, T, OUT, e_d, Prog, evalProg_d_void);
  }
  else if (T->MPType == 1)
  { // use only the fixed multi-precision
    corank = Cauchy_corank_mp(&CN, smallest_nonzero_SV, largest_zero_SV, point0_mp, point1_mp, time_mp, T, OUT, e_mp, Prog, evalProg_mp_void);
  }
  else
  { // use adaptive precision
    corank = Cauchy_corank_amp(&CN, smallest_nonzero_SV, largest_zero_SV, point0_d, point0_mp, prec0, point1_d, point1_mp, prec1, time_d, time_mp, T, OUT, e_d, e_mp, Prog, Prog, evalProg_d_void, evalProg_mp_void, change_prec_prog);
  }

  return corank;
}

////////////////////////// RANDOMIZING INSIDE OF SLP ///////////////////////

void randomize_SLP(prog_t *new_prog, prog_t *old_prog, mpq_t ***A, int A_rows, int A_cols)
/************************************************************************\
* USAGE: reads in old_prog, sets up new_prog to be [I A]*F               *
* ARGUMENTS:                                                             *
* RETURN VALUES:                                                         *
* NOTES:                                                                 *
\************************************************************************/
{// error checking - make sure that A_rows + A_cols == numFuncs
  if (A_rows + A_cols != old_prog->numFuncs)
  {
    printf("ERROR: The size of the matrix [I A] must be the same as the number of functions!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // make sure that there are no path variables
  if (old_prog->numPathVars > 0)
  {
    printf("ERROR: Path variables are not allowed when adding a patch to the straight-line program.!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }
  // make sure that there are no parameters
  if (old_prog->numPars > 0)
  {
    printf("ERROR: Parameters are not allowed when adding a patch to the straight-line program.!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  int i, j, k, new_workspace_addr, old_workspace_addr, A_addr, room_before_temps, A_num_index, tmp_func_addr;
  int *addr_conv_table = (int *)bmalloc(old_prog->memSize * sizeof(int));
  char ch;
  FILE *FOUT = NULL;

  int old_const_size = 0, old_subfunc_size = 0, old_func_size = 0, old_funcjac_size = 0, old_subfuncjac_size = 0;
  inst_list *old_const_insts = NULL, *old_subfunc_insts = NULL, *old_func_insts = NULL, *old_funcjac_insts = NULL, *old_subfuncjac_insts = NULL;
  int new_const_size = 0, new_subfunc_size = 0, new_func_size = 0;
  inst_list *new_const_insts = NULL, *new_subfunc_insts = NULL, *new_func_insts = NULL;

  // setup the structures
  setup_insts_from_prog(&old_const_insts, &old_const_size, &old_subfunc_insts, &old_subfunc_size, &old_func_insts, &old_func_size, &old_funcjac_insts, &old_funcjac_size, &old_subfuncjac_insts, &old_subfuncjac_size, old_prog); 

  // all we are doing is adding some new 'nums', moving functions to subfunctions, setting f = [I A]*f and Jv accordingly
  // For future reference, the order of the instructions is constants, parameters, subfunctions, functions, subfunc derivs, func derivs.

  new_prog->precision = old_prog->precision; // The precision at which evaluation should occur
  new_prog->num_var_gps = old_prog->num_var_gps;
  new_prog->var_gp_sizes = (int *)bmalloc(new_prog->num_var_gps * sizeof(int));
  for (i = 0; i < old_prog->num_var_gps; i++)
    new_prog->var_gp_sizes[i] = old_prog->var_gp_sizes[i];
  new_prog->index_of_first_number_for_proj_trans = old_prog->index_of_first_number_for_proj_trans;  //COULD BE WRONG - IS THIS USED???!!!
  // These next few aren't being used yet, so default to old ones for now...needs to be changed later if used!!!
  new_prog->numInstAtEndUpdate = old_prog->numInstAtEndUpdate;
  new_prog->numInstAtEndParams = old_prog->numInstAtEndParams;
  new_prog->numInstAtEndFnEval = old_prog->numInstAtEndFnEval;
  new_prog->numInstAtEndPDeriv = old_prog->numInstAtEndPDeriv;
  new_prog->numInstAtEndJvEval = old_prog->numInstAtEndJvEval;

  // now for the numbers and locations of various types
  new_workspace_addr = old_workspace_addr = 0;  // We'll increment these as we go through this block to get to the new and old first workspace addrs.

  /////////////////////////VARS - no changes
  new_prog->numVars = old_prog->numVars;
  new_prog->inpVars = old_prog->inpVars; //should just be 0 since vars go first.
  new_workspace_addr += new_prog->numVars;
  old_workspace_addr += old_prog->numVars;

  for (i = 0; i < old_prog->numVars; i++)
    addr_conv_table[old_prog->inpVars + i] = new_prog->inpVars + i;

  /////////////////////////PATHVARS - no changes
  new_prog->numPathVars = old_prog->numPathVars;
  new_prog->inpPathVars = new_prog->inpVars + new_prog->numVars;
  new_workspace_addr += new_prog->numPathVars;
  old_workspace_addr += old_prog->numPathVars;

  for (i = 0; i < old_prog->numPathVars; i++) // shouldn't be any path variables
    addr_conv_table[old_prog->inpPathVars + i] = new_prog->inpPathVars + i;

  /////////////////////////PARAMS - no changes
  new_prog->numPars = old_prog->numPars;
  new_prog->evalPars = new_prog->inpPathVars + new_prog->numPathVars;
  new_workspace_addr += new_prog->numPars;
  old_workspace_addr += old_prog->numPars;

  for (i = 0; i < old_prog->numPars; i++) // shouldn't be any parameters
    addr_conv_table[old_prog->evalPars + i] = new_prog->evalPars + i;

  /////////////////////////FUNCS - randomizing down to A_rows
  new_prog->numFuncs = A_rows;
  new_prog->evalFuncs = new_prog->evalPars + new_prog->numPars;
  new_workspace_addr += new_prog->numFuncs;
  old_workspace_addr += old_prog->numFuncs;

  /////////////////////////DPARAMS
  new_prog->evalDPars = new_prog->evalFuncs + new_prog->numFuncs;
  new_workspace_addr += new_prog->numPars;
  old_workspace_addr += old_prog->numPars;

  for (i = 0; i < old_prog->numPars; i++) // shouldn't be any parameters
    addr_conv_table[old_prog->evalDPars + i] = new_prog->evalDPars + i;

  /////////////////////////JFUNCS(V)
  new_prog->evalJVars = new_prog->evalDPars + new_prog->numPars;
  new_workspace_addr += new_prog->numFuncs * new_prog->numVars;
  old_workspace_addr += old_prog->numFuncs * old_prog->numVars;

  /////////////////////////JFUNCS(P)
  new_prog->evalJPars = new_prog->evalJVars + new_prog->numFuncs * new_prog->numVars;
  new_workspace_addr += new_prog->numFuncs * new_prog->numPars;
  old_workspace_addr += old_prog->numFuncs * old_prog->numPars;

  /////////////////////////CONSTANTS - add the number of constants
  new_prog->numConsts = old_prog->numConsts + A_rows * A_cols;
  new_prog->constAddr = new_prog->evalJPars + new_prog->numFuncs * new_prog->numPars;
  new_prog->IAddr = new_prog->constAddr;  //since I is always the first constant
  A_addr = new_prog->constAddr + old_prog->numConsts;  // Keeps track of first addr in A
  new_workspace_addr += new_prog->numConsts;
  old_workspace_addr += old_prog->numConsts;

  for (i = 0; i < old_prog->numConsts; i++)
    addr_conv_table[old_prog->constAddr + i] = new_prog->constAddr + i;

  /////////////////////////NUMS - add the number of new 'nums'
  new_prog->numNums = old_prog->numNums + 2 * A_rows * A_cols;
  new_prog->numAddr = new_prog->constAddr + new_prog->numConsts;
  new_workspace_addr += new_prog->numNums;
  old_workspace_addr += old_prog->numNums;

  for (i = 0; i < old_prog->numNums; i++)
    addr_conv_table[old_prog->numAddr + i] = new_prog->numAddr + i;

  /////////////////////////SUBFUNCS - put each old function as a subfunction
  new_prog->numSubfuncs = old_prog->numSubfuncs + old_prog->numFuncs;
  new_prog->evalSubs = new_prog->numAddr + new_prog->numNums;
  new_workspace_addr += new_prog->numSubfuncs;
  old_workspace_addr += old_prog->numSubfuncs;

  k = new_prog->evalSubs;  // setup k as the start addr for subfuncs
  for (i = 0; i < old_prog->numSubfuncs; i++)  // Just copy in the old subfuncs first
    addr_conv_table[old_prog->evalSubs + i] = k++;
  for (i = 0; i < old_prog->numFuncs; i++)  // Next we add one for each of the old funcs
    addr_conv_table[old_prog->evalFuncs + i] = k++;  // We record the function addrs in the appropriate new subfunc locations (for reasons described below) 

  /////////////////////////JSUBS(V)
  new_prog->evalJSubsV = new_prog->evalSubs + new_prog->numSubfuncs;
  new_workspace_addr += new_prog->numSubfuncs * new_prog->numVars;
  old_workspace_addr += old_prog->numSubfuncs * old_prog->numVars;

  for (i = 0; i < old_prog->numSubfuncs * old_prog->numVars; i++) // setup old subfunctions in conv table to correspond to old subfunc derivs
    addr_conv_table[old_prog->evalJSubsV + i] = k++;

  /////////////////////////JSUBS(P)
  new_prog->evalJSubsP = new_prog->evalJSubsV + new_prog->numSubfuncs * new_prog->numVars;
  new_workspace_addr += new_prog->numSubfuncs * new_prog->numPars;
  old_workspace_addr += old_prog->numSubfuncs * old_prog->numPars;

  for (i = 0; i < old_prog->numSubfuncs * old_prog->numPars; i++) // shouldn't be any parameters
    addr_conv_table[old_prog->evalJSubsP + i] = new_prog->evalJSubsP + i;

  // remember the beginning of the workspace:
  room_before_temps = new_workspace_addr;

  // the rest of the old workspace was set aside as temp memLocs
  for (i = old_workspace_addr; i < old_prog->memSize; i++)
    addr_conv_table[i] = new_workspace_addr++;

////////////////////// CONSTRUCT THE NEW SLP

  // setup new_prog->nums
  new_prog->numNums = old_prog->numNums + 2 * A_rows * A_cols;
  new_prog->nums = (num_t *)bmalloc(new_prog->numNums * sizeof(num_t));
  for (i = 0; i < old_prog->numNums; i++)
  { // copy in the old numbers
    mpq_init(new_prog->nums[i].rat);
    mpq_set(new_prog->nums[i].rat, old_prog->nums[i].rat);

    new_prog->nums[i].currPrec = new_prog->precision;
    mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
    mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);
  }

  // add in the new numbers for A
  A_num_index = i; // remember where the numbers for A start
  for (j = 0; j < A_rows; j++)
    for (k = 0; k < A_cols; k++)
    { // real part of A[j][k]
      mpq_init(new_prog->nums[i].rat);
      mpq_set(new_prog->nums[i].rat, A[j][k][0]);

      new_prog->nums[i].currPrec = new_prog->precision;
      mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
      mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);

      // increment i
      i++;

      // imag part of A[j][k]
      mpq_init(new_prog->nums[i].rat);
      mpq_set(new_prog->nums[i].rat, A[j][k][1]);

      new_prog->nums[i].currPrec = new_prog->precision;
      mpf_init2(new_prog->nums[i].real, new_prog->nums[i].currPrec);
      mpf_set_q(new_prog->nums[i].real, new_prog->nums[i].rat);

      // increment i
      i++;
    }

///////////////////////////////////////////
// ADDING CONST EVAL INSTS TO NEW SLP
///////////////////////////////////////////

  // First is the constants - prepare new_const_insts
  new_const_size = new_prog->numConsts - 2; // -2 b/c I & Pi is #0 & #1 (no update insts)
  new_const_insts = (inst_list *)bmalloc(new_const_size * sizeof(inst_list));

  // copy over all of the old constant evaluations
  for (i = 0; i < old_const_size; i++)
  {
    new_const_insts[i].num = old_const_insts[i].num;
    new_const_insts[i].insts = (inst *)bmalloc(new_const_insts[i].num * sizeof(inst));
    for (j = 0; j < new_const_insts[i].num; j++)
    { // copy over the instructions
      new_const_insts[i].insts[j].op = old_const_insts[i].insts[j].op;
      new_const_insts[i].insts[j].memLoc = addr_conv_table[old_const_insts[i].insts[j].memLoc];
      new_const_insts[i].insts[j].in[0] = addr_conv_table[old_const_insts[i].insts[j].in[0]];
      new_const_insts[i].insts[j].in[1] = addr_conv_table[old_const_insts[i].insts[j].in[1]];
    }
  }

  // setup the constant evaluations for A 
  k = A_num_index; // where the numbers for the patch start
  for (i = old_const_size; i < new_const_size; i++)
  {
    new_const_insts[i].num = 2; // each new constant creates 2 instructions (1) imag_part * I & (2) real_part + imag_part * I
    new_const_insts[i].insts = (inst *)bmalloc(new_const_insts[i].num * sizeof(inst));

    // imag_part * I
    new_const_insts[i].insts[0].op = '*'; // multiply
    new_const_insts[i].insts[0].memLoc = new_workspace_addr;
    new_const_insts[i].insts[0].in[0] = new_prog->numAddr + k + 1; // for imaginary part of the number
    new_const_insts[i].insts[0].in[1] = new_prog->IAddr;

    // real_part + imag_part * I
    new_const_insts[i].insts[1].op = '+'; // add
    new_const_insts[i].insts[1].memLoc = new_prog->constAddr + i + 2; // store result in appropriate addr (+2 b/c of I & Pi)
    new_const_insts[i].insts[1].in[0] = new_workspace_addr;
    new_const_insts[i].insts[1].in[1] = new_prog->numAddr + k; // for real part of the number

    // increment k
    k += 2;
  }
  new_workspace_addr++;

///////////////////////////////////////////
// COPYING SUBFUNC INSTS TO NEW SLP
///////////////////////////////////////////

  // Next is the subfunctions - prepare new_subfunc_insts
  new_subfunc_size = new_prog->numSubfuncs;
  new_subfunc_insts = (inst_list *)bmalloc(new_subfunc_size * sizeof(inst_list));

  // copy over all of the old subfunction evaluations
  for (i = 0; i < old_subfunc_size; i++)
  {
    new_subfunc_insts[i].num = old_subfunc_insts[i].num;
    new_subfunc_insts[i].insts = (inst *)bmalloc(new_subfunc_insts[i].num * sizeof(inst));
    for (j = 0; j < new_subfunc_insts[i].num; j++)
    { // copy over the instructions
      new_subfunc_insts[i].insts[j].op = old_subfunc_insts[i].insts[j].op;
      new_subfunc_insts[i].insts[j].memLoc = addr_conv_table[old_subfunc_insts[i].insts[j].memLoc];
      new_subfunc_insts[i].insts[j].in[0] = addr_conv_table[old_subfunc_insts[i].insts[j].in[0]];
      new_subfunc_insts[i].insts[j].in[1] = addr_conv_table[old_subfunc_insts[i].insts[j].in[1]];
    }
  }

/////
/////  ADDING THE OLD FUNCS AS NEW SUBFUNCS
/////

  // The second subfunc insts are for the old funcs -> we copy over all of the old func evaluations.
  for (k = 0; k < old_prog->numFuncs; k++) 
  {
    new_subfunc_insts[i].num = old_func_insts[k].num;
    new_subfunc_insts[i].insts = (inst *)bmalloc(new_subfunc_insts[i].num * sizeof(inst));
    for (j = 0; j < new_subfunc_insts[i].num; j++)  
    { // copy over the instructions
      new_subfunc_insts[i].insts[j].op = old_func_insts[k].insts[j].op;
      new_subfunc_insts[i].insts[j].memLoc = addr_conv_table[old_func_insts[k].insts[j].memLoc];
      new_subfunc_insts[i].insts[j].in[0] = addr_conv_table[old_func_insts[k].insts[j].in[0]];
      new_subfunc_insts[i].insts[j].in[1] = addr_conv_table[old_func_insts[k].insts[j].in[1]];
    }
    i++;
  }

///////////////////////////////////////////
// ADDING FUNC INSTS TO NEW SLP
///////////////////////////////////////////

  // Last is the functions - prepare new_func_insts
  new_func_size = A_rows;
  new_func_insts = (inst_list *)bmalloc(new_func_size * sizeof(inst_list));

  // add instructions for evaluating [I A]*F
  for (i = 0; i < new_func_size; i++)
  {
    new_func_insts[i].num = 2 + 2 * A_cols;
    new_func_insts[i].insts = (inst *)bmalloc(new_func_insts[i].num * sizeof(inst));

    tmp_func_addr = new_workspace_addr++;
    k = 0; // instruction number
    // initialize to the function i
    new_func_insts[i].insts[k].op = '='; // set
    new_func_insts[i].insts[k].memLoc = tmp_func_addr;
    new_func_insts[i].insts[k].in[0] = addr_conv_table[old_prog->evalFuncs + i];
    new_func_insts[i].insts[k].in[1] = 0; // unary set
    // increment k
    k++;

    for (j = 0; j < A_cols; j++)
    { // do sum += A[i][j] * f[A_rows + j]
      new_func_insts[i].insts[k].op = '*'; // multiply
      new_func_insts[i].insts[k].memLoc = new_workspace_addr;
      new_func_insts[i].insts[k].in[0] = A_addr + i * A_cols + j; // A[i][j]
      new_func_insts[i].insts[k].in[1] = addr_conv_table[old_prog->evalFuncs + A_rows + j]; // A_rows + jth function
      // increment
      new_workspace_addr++; 
      k++;

      // add to the sum
      new_func_insts[i].insts[k].op = '+'; // addition
      new_func_insts[i].insts[k].memLoc = new_workspace_addr;
      new_func_insts[i].insts[k].in[0] = tmp_func_addr;
      new_func_insts[i].insts[k].in[1] = new_workspace_addr - 1;
      // update & increment
      tmp_func_addr = new_workspace_addr;
      new_workspace_addr++;
      k++;
    }
    // store to proper location
    new_func_insts[i].insts[k].op = '='; // store
    new_func_insts[i].insts[k].memLoc = new_prog->evalFuncs + i; // store as ith function
    new_func_insts[i].insts[k].in[0] = tmp_func_addr;
    new_func_insts[i].insts[k].in[1] = 0;
  }

///////////////////////////////////////////////////
//  COPYING ALL NEW INSTS TO FILE (FOR DIFF)
///////////////////////////////////////////////////

  // delete the old files
  remove("finalFile.out");
  remove("jacV.out");
  remove("jacP.out");
  remove("arr.out");

  // open up finalFile.out
  FOUT = fopen("finalFile.out", "w");
  // write the instructions to FOUT
  write_insts_to_file(FOUT, new_prog, new_const_insts, new_const_size, NULL, 0, new_subfunc_insts, new_subfunc_size, new_func_insts, new_func_size, room_before_temps, new_workspace_addr);
  // close FOUT
  fclose(FOUT);

////////////////////////////////////////////////
//  CALLING DIFF WITH NEW SLP
////////////////////////////////////////////////
  diff_deflatable(1, 0);

////////////////////////////////////////////////
//  COPYING COMPLETE SLP BACK INTO STRUCT
////////////////////////////////////////////////

  // read in the differentiated instructions
  FOUT = fopen("arr.out", "r");
  if (FOUT == NULL)
  {
    printf("ERROR: The file containing the straight-line program does not exist!\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // get the memSize - move through arr.out until we find the first 'd' in TotalRoomNeeded
  do
  {
    ch = fgetc(FOUT);
  } while (ch != 'd');
  // read in the size
  fscanf(FOUT, "ed %d", &new_prog->memSize);

  // read in the number of instructions - move through arr.out until we find the 'I' in NUMINST
  do
  {
    ch = fgetc(FOUT);
  } while (ch != 'I');
  fscanf(FOUT, "NST %d", &new_prog->size);

  // grab all of the instructions.
  rewind(FOUT);
  new_prog->prog = (int *)bcalloc(new_prog->size, sizeof(int));
  for (i = 0; i < new_prog->size; i++)
    fscanf(FOUT, "%d ", &new_prog->prog[i]);

  // close arr.out
  fclose(FOUT);

  // clear the memory
  free(addr_conv_table);
  clear_inst_list(old_const_insts, old_const_size);
  free(old_const_insts);
  clear_inst_list(old_subfunc_insts, old_subfunc_size);
  free(old_subfunc_insts);
  clear_inst_list(old_func_insts, old_func_size);
  free(old_func_insts);
  clear_inst_list(old_funcjac_insts, old_funcjac_size);
  free(old_funcjac_insts);
  clear_inst_list(old_subfuncjac_insts, old_subfuncjac_size);
  free(old_subfuncjac_insts);

  clear_inst_list(new_const_insts, new_const_size);
  free(new_const_insts);
  clear_inst_list(new_subfunc_insts, new_subfunc_size);
  free(new_subfunc_insts);
  clear_inst_list(new_func_insts, new_func_size);
  free(new_func_insts);

  return;
}


