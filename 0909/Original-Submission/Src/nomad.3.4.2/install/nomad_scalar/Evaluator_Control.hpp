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
  \file   Evaluator_Control.hpp
  \brief  Control of the blackbox evaluations (headers)
  \author Sebastien Le Digabel
  \date   2010-04-15
  \see    Evaluator_Control.cpp
*/
#ifndef __EVALUATOR_CONTROL__
#define __EVALUATOR_CONTROL__

#include "Barrier.hpp"
#include "Pareto_Front.hpp"
#include "Slave.hpp"

namespace NOMAD {

  /// Control of the blackbox evaluations.
  /**
     This class allows the evaluation of a list of trial point.
  */
  class Evaluator_Control : private NOMAD::Uncopyable {

  private:

    static bool       _force_quit;  ///< Forces NOMAD to terminate if Ctrl-C is pressed.
    const NOMAD::Parameters &  _p;  ///< Parameters.
    NOMAD::Evaluator        * _ev;  ///< The evaluator.

    /// Flag equal to \c true if the destructor has to delete \c this->_ev.
    bool _del_ev;

    NOMAD::Cache * _cache;       ///< Cache for true function evaluations.
    NOMAD::Cache * _sgte_cache;  ///< Cache for surrogate evaluations.

    /// List of points to be evaluated.
    std::set<NOMAD::Priority_Eval_Point> _eval_lop;

    /// Flag equal to \c true if the destructor has to erase the cache.
    bool _del_cache;

    /// Flag equal to \c true if the destructor has to erase the cache for surrogates.
    bool _del_sgte_cache;

#ifdef USE_MPI
    NOMAD::Eval_Point ** _eval_in_progress;  ///< List of evaluations in progress.
    int                  _nb_in_progress;    ///< Number of evaluations in progress.

    /// Elop tag.
    /**
       Unique tag associated to \c this->eval_lop().
    */
    int _elop_tag;

    /// Elop tag for each slave.
    /**
       _slaves_elop_tags[k] corresponds to the Elop tag for slave \c k.
    */
    int * _slaves_elop_tags;

    NOMAD::Slave * _slave; // Slave object for master process
#endif

    NOMAD::Stats & _stats;  ///< Stats

    /**
       - Lattest tag of a point that has been written for \c display_stats
         or \c stats_file.
       - It is used to avoid any double display at the end of a run.
    */
    mutable int _last_stats_tag;

    ///< Same as \c this->_last_stats_bbe for the number of blackbox evaluations.
    mutable int _last_stats_bbe;

    /// Process an already evaluated Eval_Point.
    /**
       \param x            The point        -- \b IN.
       \param barrier      The barrier      -- \b IN/OUT.
       \param pareto_front The Pareto front -- \b IN/OUT
                           (may be NULL).
    */
    void process_eval_point ( const NOMAD::Eval_Point & x            ,
			      NOMAD::Barrier          & barrier      ,
			      NOMAD::Pareto_Front     * pareto_front   ) const;
    
    /// Save the solution file or update the history file.
    /**
       \param file_name Name of the file -- \b IN.
       \param x         Lattest solution -- \b IN.
       \param is_sol    Flag equal to \c true if the file is
                        a solution file; otherwise it is a
			history file.
    */
    void write_sol_or_his_file ( const std::string        & file_name ,
				 const NOMAD::Eval_Point  & x         ,
				 bool                       is_sol      ) const;

    /// Display evaluation result.
    /**
       \param x                Lattest evaluation                  -- \b IN.
       \param display_degree   Display degree                      -- \b IN.
       \param search           Search type                         -- \b IN.
       \param one_eval_success Success for one evaluation          -- \b IN.
       \param success          Success for a series of evaluations -- \b IN.
    */
    void display_eval_result ( const NOMAD::Eval_Point  & x                ,
			       NOMAD::dd_type             display_degree   ,
			       NOMAD::search_type         search           ,
			       NOMAD::success_type        one_eval_success ,
			       NOMAD::success_type        success            ) const;

    /// Check if evaluations have to be stopped.
    /**
       Checks the opportunistic strategy stopping criterion.
       \param x                The lattest evaluation                     -- \b IN.
       \param search           Search type                                -- \b IN.
       \param k                Evaluation index                           -- \b IN.
       \param nb_points        Number of points to evaluate               -- \b IN.
       \param stop             \c true if the algorithm has to be stopped -- \b IN.
       \param display_degree   Display degree                             -- \b IN.
       \param one_eval_success Success for one evaluation                 -- \b IN.
       \param success          Success for a series of evaluations        -- \b IN.
       \param init_nb_eval     Initial number of evaluations              -- \b IN.
       \param f0               Initial value of the objective function    -- \b IN
                               (may be undefined).
       \param barrier          Barrier                                    -- \b IN.
       \param nb_success       Number of successes                        -- \b IN/OUT.
       \param one_for_luck     \c true if one additional evaluation is made
                               according to parameter \c OPPORTUNISTIC_LUCKY_EVAL
			                                                  -- \b IN/OUT.
    */
    bool stop_evaluations ( const NOMAD::Eval_Point & x                ,
			    NOMAD::search_type        search           ,
			    int                       k                ,
			    int                       nb_points        ,
			    bool                      stop             ,
			    NOMAD::dd_type            display_degree   ,
			    NOMAD::success_type       one_eval_success ,
			    NOMAD::success_type       success          ,
			    int                       init_nb_eval     ,
			    const NOMAD::Double     & f0               ,
			    const NOMAD::Barrier    & barrier          ,
			    int                     & nb_success       ,
			    bool                    & one_for_luck       ) const;
    
    /// Check the opportunistic strategy stopping criterion.
    /**
       \param display_degree   Display degree                             -- \b IN.
       \param one_eval_success Success for one evaluation                 -- \b IN.
       \param init_nb_eval     Initial number of evaluations              -- \b IN.
       \param f0               Initial value of the objective function    -- \b IN
                               (may be undefined).
       \param barrier          Barrier                                    -- \b IN.
       \param nb_success       Number of successes                        -- \b IN/OUT.
       \param one_for_luck     \c true if one additional evaluation is made
                               according to parameter OPPORTUNISTIC_LUCKY_EVAL
			                                                  -- \b IN/OUT.
       \return A boolean equal to \c true to stop the evaluations;
               Otherwise, continue the evaluations.
    */
    bool check_opportunistic_criterion ( NOMAD::dd_type         display_degree   ,
					 NOMAD::success_type    one_eval_success ,
					 int                    init_nb_eval     ,
					 const NOMAD::Double  & f0               ,
					 const NOMAD::Barrier & barrier          ,
					 int                  & nb_success       ,
					 bool                 & one_for_luck ) const;
    
    /// Count the output stats (STAT_SUM and STAT_AVG).
    /**
       \param x Lattest evaluation -- \b IN.
    */
    void count_output_stats ( const NOMAD::Eval_Point & x );
    
    /// Clear the list of evaluation points.
    void clear_eval_lop ( void );

    /// Search a point in the cache.
    /**
       \param x              The point                           -- \b IN/OUT.
       \param true_barrier   Barrier for true evaluations        -- \b IN/OUT.
       \param sgte_barrier   Barrier for surrogate evaluations   -- \b IN/OUT.
       \param pareto_front   A pointer to the Pareto front       -- \b IN/OUT
                             (may be \c NULL).
       \param count_eval     Count or not a simulated evaluation -- \b OUT.
       \param h_max          Maximal feasibility value           -- \b IN.
       \param display_degree Display degree                      -- \b IN.
       \return A boolean equal to \c true if the point is in cache.
    */
    bool cache_check ( const NOMAD::Eval_Point *& x              ,
		       NOMAD::Barrier           & true_barrier   ,
		       NOMAD::Barrier           & sgte_barrier   ,
		       NOMAD::Pareto_Front      * pareto_front   ,
		       bool                     & count_eval     ,
		       const NOMAD::Double      & h_max          ,
		       NOMAD::dd_type             display_degree   ) const;
    
    /// Evaluate a point.
    /*
      \param x            The point                         -- \b IN/OUT.
      \param true_barrier Barrier for true evaluations      -- \b IN/OUT.
      \param sgte_barrier Barrier for surrogate evaluations -- \b IN/OUT.
      \param pareto_front A pointer to the Pareto front     -- \b IN/OUT
                          (may be \c NULL).
      \param count_eval   Count or not the evaluation       -- \b OUT.
      \param stop         Stop flag                         -- \b IN/OUT.
      \param stop_reason  Stop reason                       -- \b OUT.
      \param h_max        Maximal feasibility value         -- \b IN.
    */
    void eval_point ( const NOMAD::Eval_Point * x            ,
		      NOMAD::Barrier          & true_barrier ,
		      NOMAD::Barrier          & sgte_barrier ,
		      NOMAD::Pareto_Front     * pareto_front ,
		      bool                    & count_eval   ,
		      bool                    & stop         ,
		      NOMAD::stop_type        & stop_reason  ,
		      const NOMAD::Double     & h_max          );
    
    /// Check if a search is opportunistic or not.
    /**
       \param t Search type -- \b IN.
       \return A boolean equal to \c true if the search is opportunistic.
    */
    bool is_opportunistic ( NOMAD::search_type t ) const;

    /// Check stopping criteria.
    /**
       \param search       Search type                 -- \b IN.
       \param count_eval   Count or not the evaluation -- \b IN.
       \param x            Lattest evaluation          -- \b IN/OUT.
       \param stop         Stop flag                   -- \b IN/OUT.
       \param stop_reason  Stop reason                 -- \b OUT.
    */
    void check_stopping_criteria ( NOMAD::search_type        search      ,
				   bool                      count_eval  ,
				   const NOMAD::Eval_Point * x           ,
				   bool                    & stop        ,
				   NOMAD::stop_type        & stop_reason   ) const;
    
#ifdef USE_MPI

    /// Receive an evaluation result from a slave.
    /**
       \param search       Search type                       -- \b IN.
       \param x            A pointer to the evaluation       -- \b OUT.
       \param true_barrier Barrier for true evaluations      -- \b IN/OUT.
       \param sgte_barrier Barrier for surrogate evaluations -- \b IN/OUT.
       \param pareto_front A pointer to the Pareto front     -- \b IN/OUT
                           (may be \c NULL).
       \param slave_rank   Slave rank                        -- \b IN.	   
       \param stop         Stop flag                         -- \b IN/OUT.
       \param stop_reason  Stop reason                       -- \b OUT.
    */
    void receive_eval_result ( NOMAD::search_type    search       ,
			       NOMAD::Eval_Point   * x            ,
			       NOMAD::Barrier      & true_barrier ,
			       NOMAD::Barrier      & sgte_barrier ,
			       NOMAD::Pareto_Front * pareto_front ,
			       int                   slave_rank   ,
			       bool                & stop         ,
			       NOMAD::stop_type    & stop_reason    );
    
    /// Check if the evaluation of a point is already in progress.
    /**
       \param x The point -- \b IN.
       \return A boolean equal to \c true if the evaluation
               of the point is already in progress.
    */
    bool already_in_progress ( const NOMAD::Eval_Point & x ) const;
    
#endif

    /// Evaluate a list of points.
    /**
     - "Evaluate a list of points" is abbreviated with "eval_lop".
     - This method is private and is called by the public eval_lop method.
     - This method is defined twice in the .cpp for scalar and parallel versions.
     - The method clears the list of points \c this->_eval_lop.
     \param search         Search type                             -- \b IN.
     \param true_barrier   Barrier for true evaluations            -- \b IN/OUT.
     \param sgte_barrier   Barrier for surrogate evaluations       -- \b IN/OUT.
     \param pareto_front   A pointer to the Pareto front           -- \b IN/OUT
                           (may be \c NULL).
     \param stop           Stop flag                               -- \b IN/OUT.
     \param stop_reason    Stop reason                             -- \b OUT.
     \param new_feas_inc   Pointer to the new feasible incumbent   -- \b OUT.
     \param new_infeas_inc Pointer to the new infeasible incumbent -- \b OUT.
     \param success        Success for these evaluations           -- \b OUT.
     \param evaluated_pts  List of processed points                -- \b OUT.
    */
    void private_eval_list_of_points
    ( NOMAD::search_type                     search         ,
      NOMAD::Barrier                       & true_barrier   ,
      NOMAD::Barrier                       & sgte_barrier   ,
      NOMAD::Pareto_Front                  * pareto_front   ,
      bool                                 & stop           ,
      NOMAD::stop_type                     & stop_reason    ,
      const NOMAD::Eval_Point             *& new_feas_inc   ,
      const NOMAD::Eval_Point             *& new_infeas_inc ,
      NOMAD::success_type                  & success        ,
      std::list<const NOMAD::Eval_Point *> & evaluated_pts    );

  public:

    /// Constructor.
    /**
       \param p          Parameters                     -- \b IN.
       \param stats      Stats                          -- \b IN/OUT.
       \param ev         Pointer to the Evaluator       -- \b IN     (may be \c NULL).
       \param cache      Pointer to the cache           -- \b IN/OUT (may be \c NULL).
       \param sgte_cache Pointer to the surrogate cache -- \b IN/OUT (may be \c NULL).
    */
    Evaluator_Control ( const NOMAD::Parameters & p          ,
			NOMAD::Stats            & stats      ,
			NOMAD::Evaluator        * ev         ,
			NOMAD::Cache            * cache      ,
			NOMAD::Cache            * sgte_cache   );
    
    /// Destructor.
    virtual ~Evaluator_Control ( void );
    
    /// Force quit (called by pressing Ctrl-C).
    static void force_quit ( void ) { Evaluator_Control::_force_quit = true; }

    /// Display in stats file according to parameter \c STATS_FILE.
    /**
       \param file_name Name of the output file             -- \b IN.
       \param x         Pointer to the lattest evaluation   -- \b IN (may be \c NULL).
       \param multi_obj Pointer to several objective values -- \b IN (may be \c NULL).
    */
    void stats_file ( const std::string       & file_name ,
		      const NOMAD::Eval_Point * x         ,
		      const NOMAD::Point      * multi_obj   ) const;
    
    /// Display stats during NOMAD::Mads::run() for minimal and normal display.
    /**
       \param header    Boolean equal to \c true if a header has to be displayed
                                                            -- \b IN.
       \param out       Display                             -- \b IN.
       \param stats     List of stats to display            -- \b IN.
       \param x         Pointer to the lattest evaluation   -- \b IN (may be \c NULL).
       \param multi_obj Pointer to several objective values -- \b IN (may be \c NULL).
    */
    void display_stats ( bool                           header    ,
			 const NOMAD::Display         & out       ,
			 const std::list<std::string> & stats     ,
			 const NOMAD::Eval_Point      * x         ,
			 const NOMAD::Point           * multi_obj ) const;

    /// Display a real according to parameter \c DISPLAY_STATS.
    /**
       \param out    Display             -- \b IN.
       \param d      The real to display -- \b IN.
       \param format The format
                     -- \b IN -- \b optional (default = empty string).
    */
    void display_stats_real ( const NOMAD::Display & out    ,
			      const NOMAD::Double  & d      ,
			      const std::string    & format = ""  ) const;

    /// Display an integer according to parameter \c DISPLAY_STATS.
    /**
       \param out    Display                -- \b IN.
       \param i      The integer to display -- \b IN.
       \param max_i  Maximal value of \c i used to determine the display width
                     -- \b IN -- \b optional (default = \c -1).
       \param format The format
                     -- \b IN -- \b optional (default = empty string).
    */
    void display_stats_int ( const NOMAD::Display & out         ,
			     int                    i           ,
			     int                    max_i  = -1 ,
			     const std::string    & format = ""   ) const;

    /// Display a point according to parameter \c DISPLAY_STATS.
    /**
       \param out           Display                        -- \b IN.
       \param display_stats List of stats to display       -- \b IN.
       \param it            Iterator for the list \c stats -- \b IN/OUT.
       \param x             Pointer to the point           -- \b IN (may be \c NULL).
    */
    void display_stats_point ( const NOMAD::Display                   & out           ,
			       const std::list<std::string>           & display_stats ,
			       std::list<std::string>::const_iterator & it            ,
			       const NOMAD::Point                     * x ) const;
    /// Save the caches.
    /**
       \param overwrite A boolean equal to \c true if the cache files
                        may be overwritten -- \b IN.
       \return A boolean equal to \c true if the caches could be saved.
    */
    bool save_caches ( bool overwrite );

    /// Update the solution file.
    /**
       \param x The lattest solution -- \b IN.
    */
    void write_solution_file ( const NOMAD::Eval_Point & x ) const;
    
    /// Update a barrier from another barrier.
    /**
       Update barrier \c b1 from points in barrier \c b2 and treat
       these points as evaluations (used in VNS search).
       \param b1             First barrier                 -- \b IN/OUT.
       \param b2             Second barrier                -- \b IN.
       \param pareto_front   A pointer to the Pareto front -- \b IN/OUT
                             (may be \c NULL).
       \param display_degree Display degree                -- \b IN.
       \param search         Search type                   -- \b IN.
       \return Success for these simulated evaluations.
    */
    NOMAD::success_type process_barrier_points ( NOMAD::Barrier       & b1             ,
						 const NOMAD::Barrier & b2             ,
						 NOMAD::Pareto_Front  * pareto_front   ,
						 NOMAD::dd_type         display_degree ,
						 NOMAD::search_type     search  ) const;
    /// Access to the cache (#1).
    /**
       Non-const version.
       \return The cache.
    */
    NOMAD::Cache & get_cache ( void ) { return *_cache; }

    /// Access to the cache (#2).
    /**
       Const version.
       \return The cache.
    */
    const NOMAD::Cache & get_cache ( void ) const { return *_cache; }

    /// Access to the surrogate cache (#1).
    /**
       Non-const version.
       \return The surrogate cache.
    */
    NOMAD::Cache & get_sgte_cache ( void ) { return *_sgte_cache; }

    /// Access to the surrogate cache (#2).
    /**
       Const version.
       \return The surrogate cache.
    */
    const NOMAD::Cache & get_sgte_cache ( void ) const { return *_sgte_cache; }

    /// Access to the evaluator.
    /**
       \return The evaluator.
    */
    NOMAD::Evaluator * get_evaluator ( void ) const { return _ev; }

    /// Access to \c _last_stats_tag.
    /**
       \return \c _last_stats_tag.
    */
    int get_last_stats_tag ( void ) const { return _last_stats_tag; }

    /// Access to _last_stats_bbe.
    /**
       \return _last_stats_bbe.
    */
    int get_last_stats_bbe ( void ) const { return _last_stats_bbe; }
    
    /// Access to the list of evaluation points.
    /**
       \return The list of evaluation points.
    */
    const std::set<NOMAD::Priority_Eval_Point> & get_eval_lop ( void ) const
    {
      return _eval_lop;
    }
    
    /// Set the evaluator.
    /**
       \param e A pointer to the evaluator -- \b IN.
    */
    void set_evaluator ( NOMAD::Evaluator * e ) { _ev = e; }
    
    /// Reset.
    void reset ( void ) { _last_stats_tag = _last_stats_bbe = -1; }

    /// Evaluation of a list of points (public version that calls the private version).
    /**
       \param search         Search type                             -- \b IN.
       \param true_barrier   Barrier for true evaluations            -- \b IN/OUT.
       \param sgte_barrier   Barrier for surrogate evaluations       -- \b IN/OUT.
       \param pareto_front   A pointer to the Pareto front           -- \b IN/OUT
                             (may be \c NULL).
       \param stop           Stop flag                               -- \b IN/OUT.
       \param stop_reason    Stop reason                             -- \b OUT.
       \param new_feas_inc   Pointer to the new feasible incumbent   -- \b OUT.
       \param new_infeas_inc Pointer to the new infeasible incumbent -- \b OUT.
       \param success        Success for this series of evaluations  -- \b OUT.
       \param evaluated_pts  List of processed points                -- \b OUT.
    */
    void eval_list_of_points
    ( NOMAD::search_type                     search         ,
      NOMAD::Barrier                       & true_barrier   ,
      NOMAD::Barrier                       & sgte_barrier   ,
      NOMAD::Pareto_Front                  * pareto_front   ,
      bool                                 & stop           ,
      NOMAD::stop_type                     & stop_reason    ,
      const NOMAD::Eval_Point             *& new_feas_inc   ,
      const NOMAD::Eval_Point             *& new_infeas_inc ,
      NOMAD::success_type                  & success        ,
      std::list<const NOMAD::Eval_Point *> * evaluated_pts = NULL );

    /// Add a point to the list of points to be evaluated.
    /**
       - The point has to be a dynamic object.
       - It can be deleted into the method and be \c NULL after that.
       - The point is also snapped to bounds.
       - Periodic variables are checked.
       \param x              The point                           -- \b IN/OUT.
       \param display_degree Display degree                      -- \b IN.
       \param snap_to_bounds Boolean equal to \c true if the
                             point has to be snapped to bounds   -- \b IN.
       \param f_sgte         Objective value for the surrogate   -- \b IN
                             (may be undefined).
       \param h_sgte         Feasibility value for the surrogate -- \b IN
                             (may be undefined).
    */
    void add_eval_point ( NOMAD::Eval_Point  *& x              ,
			  NOMAD::dd_type        display_degree ,
			  bool                  snap_to_bounds ,
			  const NOMAD::Double & f_sgte         ,
			  const NOMAD::Double & h_sgte           );

    /// Display the list of evaluation points \c _eval_lop.
    /**
       \param t Search type -- \b IN
              -- \b optional (default = NOMAD::UNDEFINED_SEARCH ).
    */
    void display_eval_lop ( NOMAD::search_type t = NOMAD::UNDEFINED_SEARCH ) const;

#ifdef USE_MPI

    /// Access to the number of evaluations in progress.
    /**
       \return the number of evaluations in progress.
    */
    int get_nb_eval_in_progress ( void ) const { return _nb_in_progress; }
    
    /// Wait for evaluations in progress.
    /**
       \param search        Search type                       -- \b IN.
       \param true_barrier  Barrier for true evaluations      -- \b IN/OUT.
       \param sgte_barrier  Barrier for surrogate evaluations -- \b IN/OUT.
       \param pareto_front  A pointer to the Pareto front     -- \b IN/OUT
                            (may be \c NULL).
       \param stop          Stop flag                         -- \b IN/OUT.
       \param stop_reason   Stop reason                       -- \b OUT.
       \param success       Success for these evaluations     -- \b OUT.
       \param evaluated_pts List of processed points          -- \b OUT.
    */
    void wait_for_evaluations ( NOMAD::search_type              search         ,
				NOMAD::Barrier                & true_barrier   ,
				NOMAD::Barrier                & sgte_barrier   ,
				NOMAD::Pareto_Front           * pareto_front   ,
				bool                          & stop           ,
				NOMAD::stop_type              & stop_reason    ,
				NOMAD::success_type           & success        ,
				std::list<const NOMAD::Eval_Point *>
				                              & evaluated_pts    );
#endif
  };
}
#endif
