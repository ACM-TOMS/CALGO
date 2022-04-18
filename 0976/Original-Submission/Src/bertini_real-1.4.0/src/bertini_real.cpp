#include "bertini_real.hpp"

int main(int argC, char *args[])
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{

	

	////
	//  INITIALIZATION
	////
	MPI_Init(&argC,&args);
	
	
	//instantiate options
	BertiniRealConfig program_options;
	SolverConfiguration solve_options;
	int MPType;
	
	
	program_options.parse_commandline(argC, args); // everybody gets to parse the command line.
	
	
	
	
	
	
	if (program_options.debugwait()) {
		
		if (solve_options.is_head()) {
			std::cout << "in debug mode, waiting to start so you can attach to this process" << std::endl;
			std::cout << "master PID: " << getpid() << std::endl;
			for (int ii=0; ii<30; ii++) {
				std::cout << "starting program in " << 30-ii << " seconds" << std::endl;
				sleep(1);
				std::cout << "\033[F"; // that's a line up movement.  only works with some terminals
			}
		}
		
		
		if (solve_options.use_parallel()) {
			MPI_Barrier(MPI_COMM_WORLD);
		}

	}
	
	
	
	
	if (solve_options.is_head()) {
		// split the input_file.  this must be called before setting up the solver config.
		parse_input_file(program_options.input_filename(), &MPType);
		
		// set up the solver configuration
		get_tracker_config(solve_options,MPType);
		
		solve_options.T.ratioTol = 0.9999999999999999999999999; // manually assert to be more permissive.  i don't really like this.
	}
	else
	{ // catch the bcast from parallel parsing (which cannot be disabled)
		int arbitrary_int;
		MPI_Bcast(&arbitrary_int,1,MPI_INT,0,MPI_COMM_WORLD); // first is declaration that a parse is going to happen.
		MPI_Bcast(&arbitrary_int,1,MPI_INT,0,MPI_COMM_WORLD); // second is the bcast from the parse.
															  //yes, you do have to do both of these.
	}
	
	
	
	
	if (solve_options.use_parallel()) { // everybody participates in this.
		MPI_Bcast(&solve_options.path_number_modulus,1,MPI_INT,0,MPI_COMM_WORLD); // first is declaration that a parse is going to happen.
		bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	}
	
	
	initMP(solve_options.T.Precision); // set up some globals.
	
	
	solve_options.use_midpoint_checker = 0;
	solve_options.verbose_level(program_options.verbose_level());
	solve_options.use_gamma_trick = program_options.use_gamma_trick();
	

	
	
	if (solve_options.is_head()) {
		UbermasterProcess current_process(program_options, solve_options);
		
		try{
			current_process.main_loop();
		}
		catch (std::runtime_error e)
		{
			std::cout << "uncaught runtime exception from ubermaster process caught in bertini_real;" << std::endl;
			std::cout << e.what() << std::endl;
		}
		catch (std::logic_error e)
		{
			std::cout << "uncaught logic exception from ubermaster process caught in bertini_real;" << std::endl;
			std::cout << e.what() << std::endl;
		}
		catch (std::exception e){
			std::cout << "uncaught other exception from ubermaster process caught in bertini_real;" << std::endl;
			std::cout << e.what() << std::endl;
		}
	}
	else{
		WorkerProcess current_process(program_options, solve_options);
		
		try{
			current_process.main_loop();
		}
		catch (std::runtime_error e)
		{
			std::cout << "uncaught runtime exception from ubermaster process caught in bertini_real;" << std::endl;
			std::cout << e.what() << std::endl;
		}
		catch (std::logic_error e)
		{
			std::cout << "uncaught logic exception from ubermaster process caught in bertini_real;" << std::endl;
			std::cout << e.what() << std::endl;
		}
		catch (std::exception e){
			std::cout << "uncaught other exception from worker process caught in bertini_real;" << std::endl;
			std::cout << e.what() << std::endl;
		}
	}
	
	
	clearMP();

	MPI_Finalize();
	
	
	return 0;
}



