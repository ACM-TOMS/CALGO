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
  \file   Priority_Eval_Point.hpp
  \brief  Evaluation point with a priority (headers)
  \author Sebastien Le Digabel
  \date   2010-04-22
  \see    Priority_Eval_Point.cpp
*/
#ifndef __PRIORITY_EVAL_POINT_
#define __PRIORITY_EVAL_POINT_

#include "Set_Element.hpp"
#include "Eval_Point.hpp"

namespace NOMAD {

  /// Evaluation point with a priority.
  class Priority_Eval_Point : public NOMAD::Set_Element<NOMAD::Eval_Point> {

  private:

    NOMAD::Double _h_min;              ///< \c h_min value for comparison operator.
    NOMAD::Double _f_sgte;             ///< Objective surrogate value.
    NOMAD::Double _h_sgte;             ///< Feasibility surrogate value.
    NOMAD::Double _angle_success_dir;  ///< Angle with last successful direction.
    NOMAD::Double _angle_simplex_grad; ///< Angle with simplex gradient.

    /// Affectation operator.
    /**
       \param x The right-hand side object -- \b IN.
    */
    Priority_Eval_Point & operator = ( const Priority_Eval_Point & x );

    /// Compare the \c h and \c f values of two points.
    /**
       The two points to compare are \c x1 and \c x2.
       \param hx1 \c h(x1) -- \b IN.
       \param fx1 \c f(x1) -- \b IN.
       \param hx2 \c h(x2) -- \b IN.
       \param fx2 \c f(x2) -- \b IN.
       \return \c (h(x1),f(x1)) \c < \c (h(x2),f(x2))
               with the following format:
	       -  1: \c x1 best than \c x2.
	       - -1: \c x2 best than \c x1.
	       -  0: undetermined.
    */
    int compare_hf_values ( const NOMAD::Double & hx1 ,
			    const NOMAD::Double & fx1 ,
			    const NOMAD::Double & hx2 ,
			    const NOMAD::Double & fx2   ) const;
  public:
  
    /// Constructor.
    /**
       \param x A pointer to the evaluation point -- \b IN.
       \param h_min \c h_min value                -- \b IN.
    */
    Priority_Eval_Point ( const NOMAD::Eval_Point * x     ,
			  const NOMAD::Double     & h_min   )
      : NOMAD::Set_Element<NOMAD::Eval_Point> ( x     ) ,
	_h_min                                ( h_min )   {}

    /// Copy constructor.
    /**
       \param pep The copied object -- \b IN.
    */
    explicit Priority_Eval_Point ( const Priority_Eval_Point & pep )
      : NOMAD::Set_Element<NOMAD::Eval_Point> ( pep                     ) ,
	_h_min                                ( pep._h_min              ) ,
	_f_sgte                               ( pep._f_sgte             ) ,
	_h_sgte                               ( pep._h_sgte             ) ,
	_angle_success_dir                    ( pep._angle_success_dir  ) ,
	_angle_simplex_grad                   ( pep._angle_simplex_grad )   {}
    
    /// Destructor.
    virtual ~Priority_Eval_Point ( void ) {}
  
    /// Access to specific elements of comparison.
    /**
       - This method is defined virtual in NOMAD::Set_Element so that
         \c operator \c < \c (Set_Element x) can invoke
	 it on \c x (which is in fact a \c Priority_Eval_Point).
       - This avoids an expensive downcast in \c operator \c < .
       \param f_sgte              Objective surrogate value            -- \b OUT.
       \param h_sgte              Feasibility surrogate value          -- \b OUT.
       \param angle_last_succ_dir Angle with last successful direction -- \b OUT.
       \param angle_simplex_grad  Angle with simplex gradient          -- \b OUT.
    */
    virtual void get_priority_criteria ( NOMAD::Double & f_sgte              ,
					 NOMAD::Double & h_sgte              ,
					 NOMAD::Double & angle_last_succ_dir ,
					 NOMAD::Double & angle_simplex_grad    ) const;
    
    /// Comparison operator.
    /**
       \param x The right-hand side object -- \b IN.
       \return A boolean equal to \c true if \c *this \c < \c x.
    */
    virtual bool operator < ( const NOMAD::Set_Element<NOMAD::Eval_Point> & x ) const;

    /// Access to the evaluation point.
    /**
       \return A pointer to the evaluation point.
    */
    const NOMAD::Eval_Point * get_point ( void ) const { return get_element(); }

    /// Set the angle with last successful direction.
    /**
       \param a The angle with last successful direction -- \b IN.
    */
    void set_angle_success_dir ( const NOMAD::Double & a ) { _angle_success_dir = a; }

    /// Set the objective surrogate value.
    /**
       \param f The objective surrogate value -- \b IN.
    */
    void set_f_sgte ( const NOMAD::Double & f ) { _f_sgte = f; }

    /// Set the feasibility surrogate value.
    /**
       \param h The feasibility surrogate value -- \b IN.
    */
    void set_h_sgte ( const NOMAD::Double & h ) { _h_sgte = h; }
  };
}

#endif
