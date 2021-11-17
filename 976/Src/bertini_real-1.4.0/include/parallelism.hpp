#ifndef BR_PARALLELISM_H
#define BR_PARALLELISM_H



/** \file parallelism.hpp */

#include <boost/timer/timer.hpp>

#include "decompositions/surface.hpp"
#include "decompositions/curve.hpp"

#include "nag/nid.hpp"

/**
\defgroup mpienabled MPI-enabled classes

 */

/**
 \defgroup send MPI Sends
 */

/**
 \defgroup receive MPI Receives
 */

/**
 \defgroup bcast MPI Bcasts (sends and receives)
 */


/**
 \brief BR process base class -- holds current state of program and solver.
 
 In order to make the program_options and solve_options 'globally' accessibly, we place them into the containing process.
 */
class Process
{
protected:
	BertiniRealConfig program_options;///< holds the current state of Bertini_real
	SolverConfiguration solve_options; ///< holds the current state of the solver
	
	int MPType; ///< operating MP type.
public:

	
	
	
	
	virtual int main_loop() = 0;
	
	virtual ~Process()
	{
		
	};
	
	
	
};



/**
 \brief Master process, level 0.
 */
class UbermasterProcess : public Process
{
	
public:
	
	UbermasterProcess(BertiniRealConfig & new_options, SolverConfiguration & new_solve_options){
		this->program_options = new_options;
		this->solve_options = new_solve_options;
	}
	
	
	/**
	 \brief Master Bertini_real procedure.
	 
	 Loads the witness_data, tracker config, and decomposes components of the user's choosing.
	 \return An integer flag indicating the success of the loop.
	 */
	int main_loop();

	void bertini_real(WitnessSet & W, vec_mp *pi, VertexSet & V);
	
	
	void critreal(WitnessSet & W, vec_mp *pi, VertexSet & V);
	
	
	
	
	~UbermasterProcess()
	{
		
	}
};



/**
 \brief worker process, level 1.
 */
class WorkerProcess : public Process
{
public:
	
	WorkerProcess(BertiniRealConfig & new_options, SolverConfiguration & new_solve_options){
		this->program_options = new_options;
		this->solve_options = new_solve_options;
	}
	
	
	
	~WorkerProcess()
	{
		
	}
	
	/**
	 \brief the listen-work loop for workers in an MPI ring.
	 
	 \return SUCCESSFUL flag from process.
	 */
	int main_loop();

};








#endif





