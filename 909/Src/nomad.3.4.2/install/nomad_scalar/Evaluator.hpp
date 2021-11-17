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
  \file   Evaluator.hpp
  \brief  Evaluation of blackbox functions (headers)
  \author Sebastien Le Digabel
  \date   2010-04-14
  \see    Evaluator.cpp
*/
#ifndef __EVALUATOR__
#define __EVALUATOR__

#include <csignal>
#include "Priority_Eval_Point.hpp"
#include "Stats.hpp"

namespace NOMAD {

  // forward declarations:
  class Barrier;
  class Pareto_Front;
  class Evaluator_Control;

  /// Evaluation of blackbox functions.
  /**
     This class allows the evaluation of one trial point.
  */
  class Evaluator : private NOMAD::Uncopyable {

  protected:

    /// Parameters.
    const NOMAD::Parameters & _p;

    /// Multi-objective flag.
    /**
       Identifies if a NOMAD::Evaluator object refers to
       a NOMAD::Multi_Obj_Evaluator for more than one
       objective function.
    */
    bool _is_multi_obj;

  private:

    /// Blackbox executable names.
    /**
       - Not the same as Parameters::_bb_exe.
       - NOMAD::Evaluator::_bb_exe is constructed from
         NOMAD::Parameters::_bb_exe.
    */
    std::vector<std::string> _bb_exe;
    
    std::vector<std::string> _sgte_exe;  ///< Surrogate executable names.

    /// Number of outputs for each blackbox executable.
    std::vector<int> _bb_nbo;

    /// Check the constraints to decide if the evaluations have to be stopped.
    /**
       \param x     Point at which the constraints are checked -- \b IN.
       \param h_max Maximal feasibility value \c h_max         -- \b IN.
       \return A boolean equal to \c true if
               the evaluations have to be stopped.
    */
    bool interrupt_evaluations ( const NOMAD::Eval_Point & x     ,
				 const NOMAD::Double     & h_max   ) const;

    /// Process a blackbox executable name.
    /**
       \param bb_exe Executable name -- \b IN/OUT.
    */
    void process_bb_exe_name ( std::string & bb_exe ) const;

  public:

    /// Constructor.
    /**
       \param p Parameters -- \b IN.
    */
    Evaluator ( const NOMAD::Parameters & p );

    /// Destructor.
    virtual ~Evaluator ( void ) {}

    /// Access to the multi-objective flag.
    /**
       \return The multi-objective flag.
    */
    bool is_multi_obj ( void ) const { return _is_multi_obj; }

    /// User pre-processing of a list of points.
    /**
       - This virtual method is called before the evaluation of a list
         of points.
       - It allows the user to pre-process the points to be evaluated.
       \param pts The list of points -- \b IN/OUT.
    */
    virtual void list_of_points_preprocessing
    ( std::set<Priority_Eval_Point> & pts ) const {}

    /// User updates after a success.
    /**
       This virtual method is called every time a new (full) success is made.
       \param stats Stats                 -- \b IN.
       \param x     Last successful point -- \b IN.
    */
    virtual void update_success ( const NOMAD::Stats      & stats ,
				  const NOMAD::Eval_Point & x       ) {}

    /// User updates after an iteration.
    /**
       This virtual method is called every time a MADS iteration is terminated.
       \param success      Success of the iteration              -- \b IN.
       \param stats        Stats                                 -- \b IN.
       \param ev_control   The NOMAD::Evaluator_Control object   -- \b IN.
       \param true_barrier Barrier for true evaluations          -- \b IN.
       \param sgte_barrier Barrier for surrogate evaluations     -- \b IN.
       \param pareto_front Pareto front                          -- \b IN.
       \param stop         Allows the user to stop the algorithm -- \b OUT.
    */
    virtual void update_iteration ( NOMAD::success_type              success      ,
				    const NOMAD::Stats             & stats        ,
				    const NOMAD::Evaluator_Control & ev_control   ,
				    const NOMAD::Barrier           & true_barrier ,
				    const NOMAD::Barrier           & sgte_barrier ,
				    const NOMAD::Pareto_Front      & pareto_front ,
				    bool                           & stop           ) {}

    /// Evaluate the blackbox functions at a given trial point (#1).
    /**
       - Const version.
       - May be user-defined.
       - Surrogate or true evaluation depending on the value of \c x.is_surrogate().
       \param x          The trial point                    -- \b IN/OUT.
       \param h_max      Maximal feasibility value \c h_max -- \b IN.
       \param count_eval Flag indicating if the evaluation has to be counted
                         or not -- \b OUT.
       \return A boolean equal to \c false if the evaluation failed.
    */
    virtual bool eval_x ( NOMAD::Eval_Point   & x          ,
			  const NOMAD::Double & h_max      ,
			  bool                & count_eval   ) const;

    /// Evaluate the blackbox functions at a given trial point (#2).
    /**
       - Non-const version.
       - Calls the const version by default.
       - May be user-defined.
       - Surrogate or true evaluation depending on the value of \c x.is_surrogate().
       \param x          The trial point                    -- \b IN/OUT.
       \param h_max      Maximal feasibility value \c h_max -- \b IN.
       \param count_eval Flag indicating if the evaluation has to be counted
                         or not -- \b OUT.
       \return A boolean equal to \c false if the evaluation failed.
    */
    virtual bool eval_x ( NOMAD::Eval_Point   & x          ,
			  const NOMAD::Double & h_max      ,
			  bool                & count_eval   ) {
      return static_cast<const Evaluator *>(this)->eval_x ( x , h_max , count_eval );
    }

    /// Compute the objective value \c f(x) from the blackbox outputs of a point.
    /**
       - Define a NOMAD::Multi_Obj_Evaluator to treat more than one objective.
       \param x The point -- \b IN/OUT.
    */
    virtual void compute_f ( NOMAD::Eval_Point & x ) const;
    
    /// Compute the feasibility value \c h(x) from the blackbox outputs of a point.
    /**
       \param x The point -- \b IN/OUT.
    */
    void compute_h ( NOMAD::Eval_Point & x ) const;
  };
}

#endif
