#include <iostream>//Library for input and output data.
#include "Parameters.h"
#include "Manbis.h"
#include "Functions.h"

using namespace std;//Library for input and output data.

int main(int argc, char* argv[]){//Beginning of the main program.
	//Its arguments argc and argv are the argument count and array 
	//as passed to the main() function on program invocation.

	Parameters parameters(argc, argv);
	Manbis man(parameters, fun);

	man.FindRoots();

	return 0;//this returning value shows that the program
	//has been terminated successfully.

}//End of main






