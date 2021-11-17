/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonsmooth Optimization by Mesh Adaptive Direct search - version 3.4        */
/*                                                                                     */
/*  Copyright (C) 2001-2010  Mark Abramson        - the Boeing Company, Seattle        */
/*                           Charles Audet        - Ecole Polytechnique, Montreal      */
/*                           Gilles Couture       - Ecole Polytechnique, Montreal      */
/*                           John Dennis          - Rice University, Houston           */
/*                           Sebastien Le Digabel - Ecole Polytechnique, Montreal      */
/*                                                                                     */
/*  funded in part by AFOSR and Exxon Mobil                                            */
/*                                                                                     */
/*  Author: Sebastien Le Digabel                                                       */
/*                                                                                     */
/*  Contact information:                                                               */
/*    Ecole Polytechnique de Montreal - GERAD                                          */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada                  */
/*    e-mail: nomad@gerad.ca                                                           */
/*    phone : 1-514-340-6053 #6928                                                     */
/*    fax   : 1-514-340-5665                                                           */
/*                                                                                     */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad               */
/*-------------------------------------------------------------------------------------*/
/**
  \file   defines.hpp
  \brief  Definitions
  \author Sebastien Le Digabel
  \date   2010-03-23
*/
#ifndef __DEFINES__
#define __DEFINES__

#include <string>
#include <iostream>
#include <sstream>
#include <limits>
#include <limits.h>
#include <cstdlib>

// #define DEBUG is obtained by 'make debug'
// #ifndef DEBUG
// #define DEBUG // temporary replacement of 'make debug'
// #endif

/*----------------------------------------------------------------------*/
/*             One of these flags is active in the makefile             */
/*----------------------------------------------------------------------*/
/*  #define GCC_X                 // g++         -- Linux/Unix/Mac-osX  */
/*  #define GCC_WINDOWS           // minGW (gcc) -- Windows             */
/*  #define _MSC_VER              // visual C++  -- Windows             */
/*----------------------------------------------------------------------*/

#ifdef GCC_WINDOWS
#define WINDOWS
#endif

#ifdef _MSC_VER
#define WINDOWS
#endif

namespace NOMAD {

  /// Current version:
#ifdef USE_MPI
  const std::string VERSION = "3.4.2.MPI";
#else
  const std::string VERSION = "3.4.2";
#endif

  // Directory separator, plus LGPL and user guide files
#ifdef WINDOWS
  const char        DIR_SEP = '\\';            ///< Directory separator
  const std::string HOME    = "%NOMAD_HOME%";  ///< Home directory
#else
  const char        DIR_SEP = '/';             ///< Directory separator
  const std::string HOME    = "$NOMAD_HOME";   ///< Home directory
#endif
  
  /// Licence file
  const std::string LGPL_FILE = HOME+DIR_SEP+"src"+DIR_SEP+"lgpl.txt";

  /// User guide file
  const std::string USER_GUIDE_FILE = HOME+DIR_SEP+"doc"+DIR_SEP+"user_guide.pdf";

  /// Examples directory
  const std::string EXAMPLES_DIR = HOME+DIR_SEP+"examples";

  /// Tools directory
  const std::string TOOLS_DIR = HOME+DIR_SEP+"tools";

  /// Tag for valid cache files
#ifdef GCC_X
  const int CACHE_FILE_ID = 77041301;
#else
#ifdef GCC_WINDOWS
  const int CACHE_FILE_ID = 77041302;
#else
#ifdef _MSC_VER
  const int CACHE_FILE_ID = 77041303;
#else
  const int CACHE_FILE_ID = 77041304;
#endif
#endif
#endif

#ifdef USE_MPI
  // MPI constants
  const int   MAX_REQ_WAIT =  3 ;  ///< Maximum time to wait for a request
  const char   STOP_SIGNAL = 'S';  ///< Stop signal
  const char   EVAL_SIGNAL = 'X';  ///< Evaluation signal
  const char  READY_SIGNAL = 'R';  ///< Ready signal
  const char RESULT_SIGNAL = 'O';  ///< Result signal
  const char   WAIT_SIGNAL = 'W';  ///< Wait signal
#endif

  // Mesh index constants
  const int L_LIMITS    = 50;          ///< Limits for the mesh index values
  const int UNDEFINED_L = L_LIMITS+1;  ///< Undefined value for the mesh index
  
  /// Default epsilon used by NOMAD::Double
  /** Use Parameters::set_EPSILON(), or parameter EPSILON,
      or NOMAD::Double::set_epsilon() to change it*/
  const double DEFAULT_EPSILON = 1e-13;
  
  /// Default infinity string used by NOMAD::Double
  /** Use Parameters::set_INF_STR(), or parameter INF_STR,
      or NOMAD::Double::set_inf_str() to change it*/
  const std::string DEFAULT_INF_STR = "inf";
  
  /// Default undefined value string used by NOMAD::Double
  /** Use Parameters::set_UNDEF_STR(), or parameter UNDEF_STR,
      or NOMAD::Double::set_undef_str() to change it*/
  const std::string DEFAULT_UNDEF_STR = "-";

  // VNS constants
  const int VNS_HALTON_INDEX_0 = 997;  ///< Initial Halton index for the VNS
  const int VNS_HALTON_INCR    =  17;  ///< VNS Halton index increment

  /// log(10) (for display widths.)
  const double LOG10 = 2.30258509299;

  const double INF = std::numeric_limits<double>::max(); ///< Infinity

  const double D_INT_MAX = INT_MAX; ///< The INT_MAX constant as a \c double

  /// Default value for parameter POINT_DISPLAY_LIMIT
  /** Use Parameters::set_POINT_DISPLAY_LIMIT() or parameter POINT_DISPLAY_LIMIT
      or Point::set_display_limit() to change it */
  const int DEFAULT_POINT_DISPLAY_LIMIT = 20;

  // Display precisions.
  const int DISPLAY_PRECISION_STD = 10;  ///< Standard display precision
  const int DISPLAY_PRECISION_BB  = 15;  ///< Display precision for blackboxes

  /// Constant for blackbox files #1.
  const std::string BLACKBOX_INPUT_FILE_PREFIX = "nomad"; 

  /// Constant for blackbox files #2.
  const std::string BLACKBOX_INPUT_FILE_EXT = "input";

  /// Constant for blackbox files #3.
  const std::string BLACKBOX_OUTPUT_FILE_PREFIX = "nomad";

  /// Constant for blackbox files #4.
  const std::string BLACKBOX_OUTPUT_FILE_EXT = "output";

  /// Display degree type.
  enum dd_type
    {
      NO_DISPLAY     , ///< No display
      NORMAL_DISPLAY , ///< Normal display
      FULL_DISPLAY     ///< Full display
    };

  /// Types of the variables
  //  (do not modify this order)
  enum bb_input_type
    {
      CONTINUOUS  ,     ///< Continuous variable (default) (R)
      INTEGER     ,     ///< Integer variable              (I)
      CATEGORICAL ,     ///< Categorical variable          (C)
      BINARY            ///< Binary variable               (B)
    };

  /// Blackbox outputs type
  enum bb_output_type
    {
      OBJ         ,    ///< Objective value
      EB          ,    ///< Extreme barrier constraint
      PB          ,    ///< progressive barrier constraint
      // PEB           ///< PB constraint that becomes EB once satisfied
      PEB_P       ,    ///< PEB constraint, state P (PB)
      PEB_E       ,    ///< PEB constraint, state E (EB)
      FILTER      ,    ///< Filter constraint
      CNT_EVAL    ,    ///< Output set to 0 or 1 to specify to count or not the
                       ///<   blackbox evaluation
      STAT_AVG    ,    ///< Stat (average)
      STAT_SUM      ,  ///< Stat (sum)
      UNDEFINED_BBO    ///< Ignored output
    };
  
  /// Formulation for multi-objective optimization
  enum multi_formulation_type
    {
      NORMALIZED            , ///< Normalized formulation
      PRODUCT               , ///< Product formulation
      DIST_L1               , ///< Distance formulation with norm L1
      DIST_L2               , ///< Distance formulation with norm L2
      DIST_LINF             , ///< Distance formulation with norm Linf
      UNDEFINED_FORMULATION   ///< Undefined formulation
    };

  /// Poll type
  enum poll_type
    {
      PRIMARY   , ///< Primary poll
      SECONDARY   ///< Secondary poll
    };
  
  /// Poll center feasibility type
  enum poll_center_type
    {
      FEASIBLE                   ,  ///< Feasible poll center type
      INFEASIBLE                 ,  ///< Infeasible poll center type
      UNDEFINED_POLL_CENTER_TYPE    ///< Undefined poll center type
    };
  
  /// Search type
  enum search_type
    {
      X0_EVAL          ,  ///< Starting point evaluation
      POLL             ,  ///< Poll
      EXTENDED_POLL    ,  ///< Extended poll
      SEARCH           ,  ///< Generic search
      CACHE_SEARCH     ,  ///< Cache search (does not require evals)
      SPEC_SEARCH      ,  ///< MADS speculative search (dynamic order in GPS)
      LH_SEARCH        ,  ///< Latin Hypercube (LH) search
      LH_SEARCH_P1     ,  ///< Latin Hypercube (LH) search during phase one
      VNS_SEARCH       ,  ///< VNS search
      P1_SEARCH        ,  ///< Phase one search
      ASYNCHRONOUS     ,  ///< Parallel asynchronous final evaluations
      USER_SEARCH      ,  ///< User search
      UNDEFINED_SEARCH    ///< Undefined search
    };

  /// Success type of an iteration
  // (do not modify this order)
  enum success_type
    {
      UNSUCCESSFUL    ,  ///< Failure
      PARTIAL_SUCCESS ,  ///< Partial success (improving)
      FULL_SUCCESS       ///< Full success (dominating)
    };

  /// Stopping criteria
  enum stop_type
    {
      NO_STOP                    ,  ///< No stop
      ERROR                      ,  ///< Error
      UNKNOWN_STOP_REASON        ,  ///< Unknown
      CTRL_C                     ,  ///< Ctrl-C
      USER_STOPPED               ,  ///< User-stopped in Evaluator::update_iteration()
      TRIAL_PT_CONSTRUCTION      ,  ///< Problem with trial point construction
      X0_FAIL                    ,  ///< Problem with starting point evaluation
      P1_FAIL                    ,  ///< Problem with phase one
      DELTA_M_MIN_REACHED        ,  ///< Min mesh size
      DELTA_P_MIN_REACHED        ,  ///< Min poll size
      L_MAX_REACHED              ,  ///< Max mesh index
      L_LIMITS_REACHED           ,  ///< Mesh index limits
      MAX_TIME_REACHED           ,  ///< Max time
      MAX_BB_EVAL_REACHED        ,  ///< Max number of blackbox evaluations
      MAX_SGTE_EVAL_REACHED      ,  ///< Max number of surrogate evaluations
      MAX_EVAL_REACHED           ,  ///< Max number of evaluations
      MAX_SIM_BB_EVAL_REACHED    ,  ///< Max number of sim bb evaluations
      MAX_ITER_REACHED           ,  ///< Max number of iterations
      MIN_EPS_REACHED            ,  ///< Min epsilon
      FEAS_REACHED               ,  ///< Feasibility
      F_TARGET_REACHED           ,  ///< F_TARGET
      STAT_SUM_TARGET_REACHED    ,  ///< STAT_SUM_TARGET
      L_CURVE_TARGET_REACHED     ,  ///< L_CURVE_TARGET
      MULTI_MAX_BB_REACHED       ,  ///< Max number of blackbox evaluations (multi obj.)
      MULTI_NB_MADS_RUNS_REACHED ,  ///< Max number of MADS runs (multi obj.)
      MULTI_STAGNATION           ,  ///< Stagnation criterion (multi obj.)
      MULTI_NO_PARETO_PTS        ,  ///< No Pareto points (multi obj.)
      MAX_CACHE_MEMORY_REACHED      ///< Max cache memory
    };

  /// Type of norm used for the computation of h
  enum hnorm_type
    {
      L1   ,  ///< norm L1
      L2   ,  ///< norm L2
      LINF    ///< norm Linf
    };

  /// Types of directions
  // (do not modify this order)
  enum direction_type
    {
      UNDEFINED_DIRECTION    ,  ///< Undefined direction
      NO_DIRECTION           ,  ///< No direction
      ORTHO_1                ,  ///< OrthoMADS 1
      ORTHO_2                ,  ///< OrthoMADS 2
      ORTHO_2N               ,  ///< OrthoMADS 2n
      LT_1                   ,  ///< LT-MADS 1
      LT_2                   ,  ///< LT-MADS 2
      LT_2N                  ,  ///< LT-MADS 2n
      LT_NP1                 ,  ///< LT-MADS n+1
      GPS_BINARY             ,  ///< GPS for binary variables
      GPS_2N_STATIC          ,  ///< GPS 2n static (classic coordinate search)
      GPS_2N_RAND            ,  ///< GPS 2n random
      GPS_NP1_STATIC_UNIFORM ,  ///< GPS n+1 static uniform
      GPS_NP1_STATIC         ,  ///< GPS n+1 static
      GPS_NP1_RAND_UNIFORM   ,  ///< GPS n+1 random uniform
      GPS_NP1_RAND              ///< GPS n+1
    };
  
  /// Type for Eval_Point::check() failures
  enum check_failed_type
    {
      CHECK_OK     ,  ///< Correct check
      LB_FAIL      ,  ///< LB failure
      UB_FAIL      ,  ///< UB failure
      FIX_VAR_FAIL ,  ///< Fixed variables failure
      BIN_FAIL     ,  ///< Binary failure
      CAT_FAIL     ,  ///< Categorical failure
      INT_FAIL        ///< Integer failure
    };
  
  /// Type for cache indexes in Cache:
  enum cache_index_type
    {
      CACHE_1         ,  ///< Cache index #1
      CACHE_2         ,  ///< Cache index #2
      CACHE_3         ,  ///< Cache index #3
      UNDEFINED_CACHE    ///< Undefined cache index
    };
  
  /// Type for DISPLAY_STATS parameter
  // (do not modify this order):
  enum display_stats_type
    {
      DS_OBJ        ,    ///< Objective (f) value
                         //   (keep in first position)
      DS_SIM_BBE    ,    ///< Number of simulated bb evaluations
      DS_BBE        ,    ///< Number of bb evaluations
      DS_SGTE       ,    ///< Number of surrogate evaluations
      DS_BBO        ,    ///< All blackbox outputs
      DS_EVAL       ,    ///< Number of evaluations
      DS_TIME       ,    ///< Wall-clock time
      DS_MESH_INDEX ,    ///< Mesh index
      DS_MESH_SIZE  ,    ///< Mesh size parameter Delta^m_k
      DS_DELTA_M    ,    ///< Same as \c DS_MESH_SIZE
      DS_POLL_SIZE  ,    ///< Poll size parameter Delta^p_k
      DS_DELTA_P    ,    ///< Same as \c DS_POLL_SIZE
      DS_SOL        ,    ///< Solution vector
      DS_VAR        ,    ///< One variable
      DS_STAT_SUM   ,    ///< Stat sum
      DS_STAT_AVG   ,    ///< Stat avg
      DS_UNDEFINED       ///< Undefined value
                         //   (keep in last position)
    };
  
  /// Type for evaluation
  enum eval_type
    {
      TRUTH , ///< Truth
      SGTE    ///< Surrogate
    };
  
  /// Type for an evaluation status
  enum eval_status_type
    {
      EVAL_FAIL        ,  ///< Evaluation failure
      EVAL_OK          ,  ///< Correct evaluation
      EVAL_IN_PROGRESS ,  ///< Evaluation in progress
      UNDEFINED_STATUS    ///< Undefined evaluation status
    };

}

#endif
