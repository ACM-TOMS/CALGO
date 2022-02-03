#ifndef __Manbis_PARAMETERS_H__
#define __Manbis_PARAMETERS_H__
#include <string>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <climits>//Defines constants with the limits of fundamental integral types for the 
//specific system and compiler implementation used.
#include <getopt.h>//Is a C library function used to parse command-line options.
#include <cstdlib>//This library provides the opportunity to use EXIT_FAILURE.
#include <unistd.h>//Provides access to the POSIX operating system API.

class Parameters
{
	// non-copyable
	Parameters(const Parameters &rhs);
	Parameters & operator=(const Parameters &rhs);
	public:
		
		const char* program_name;//Holds the name of the executable
		int files;//Specifies if the user wants the creation 
		//of txt files SubsByLength.txt and SubsByXLeft.txt
		
		int digits;//'digits' is the number of the digits that will be  
		//displayed at the end of the program. 
		
		long double e,x1, x2, fraction;//'e' is the accuracy
		//'x1' is the left point of the interval.
		//'x2' is the right point of the interval.
		//'fraction' will hold the fraction of the roots that the user wants to calculate.
		
		int aqur;//the accuracy in decimal digits.

		double Z, a;//'Z' is a parameter used in the statistics
		//For a=0.05, Z=1.959964.

		int nor;//if nor=1 then the roots will be written in the file Results.txt, else only
		//their number will be printed.

		long double statdiff;//Is the wanted by the user difference between two consecutive statistical levels.

		int totaltime;//if totaltime=1 then the run time of Manbis will be displayed at the end of the program.

		int currentest;//if currentest=1 then the current estimations of Manbis will be displayed at the screen during the run of the program.
		
		
		
		Parameters();
		
		Parameters(int argc, char* argv[]);

 	private:

		//init sets some default values.
		void Init();

		//set_command_line_values gets the arguments from the command line
		void SetCommandLineValues(int &argc, char **argv);

		//checking_values_of_user checks the arguments given by the user.
		void CheckingValuesOfUser();

		//print_usage displays usage instructions 
		void PrintUsage() ;


};//End of the class Parameters.

#endif //__Manbis_PARAMETERS_H__
