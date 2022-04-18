#include "Parameters.h"

using namespace std;

Parameters::Parameters(int argc, char* argv[])
{
	program_name=argv[0];//holds the name of the executable file of the user.

	Init();//sets the default values

	SetCommandLineValues(argc, argv);//gets the arguments from the commant line

	CheckingValuesOfUser();//checks the arguments set by the user.

}
 
//Init sets some default values.
void Parameters::Init() 
{
	totaltime = 0;//Default value for the run time of Manbis.
	//User can change this value at command line. 

	currentest = 0;//Default value for the current estimations of Manbis.
	//User can change this value at command line.

	digits = 8;//Default value for the number of digits when the roots are printed.
	//User can change this value at command line. 

	e = 0.0001, x1 = LLONG_MAX, x2 = LLONG_MIN, fraction = -1;
	//'e' is the accuracy. Default value is equal to 4 decimal points (=0.0001).
	//'x1' is the left point of the interval.
	//'x2' is the right point of the interval.
	//'fraction' will hold the fraction of the roots that the user wants to calculate.

	aqur = 4;//the accuracy in decimal digits.
	//The default value will be 4;

	Z = 1.959964, a = 0.95;//'Z' is a parameter used in the statistics
	//For a=0.95 then Z=1.959964 respectively for a=0.90 then Z=1.644854.

	files = 0;//The variable named "files" by default contains the value 0.
	//If the user wants the files  SubsByLength.txt and SubsByXLeft.txt
	//he/she has to set (at the command line) the -u equals to 1.

	nor = 1;//if nor=1 then the roots will be written in the file Results.txt, else only
		//their number will be printed. Default value is 1.
	
	statdiff = 0.5;//Set the default value for the statistical difference between two consecutives levels.

}//End of function Init

//SetCommandLineValues gets the arguments from the command line
void Parameters::SetCommandLineValues(int &argc, char **argv)
{
	int opt = 0;//Initialization of the iterator opt with the value zero. 
	//It will be used in the getopt() function.

	static struct option long_options[] = {
		{"digits", required_argument, 0, 'a'},
		{"time", required_argument, 0, 't'},
		{"accuracy", required_argument, 0, 'e'},
		{"xleft", required_argument, 0, 'l'},
		{"xright", required_argument, 0, 'r'},
		{"percentage", required_argument, 0, 'q'},
		{"uncutfiles", required_argument, 0, 'u'},
		{"probability", required_argument, 0, 'p'},
		{"noroots", required_argument, 0, 'n'},
		{"statdiff", required_argument, 0, 's'},
		{"currentest", required_argument, 0, 'c'},
		{"help", required_argument, 0, 'h'},
		{0, 0, 0, 0}
	};//Specifying the expected options.
	//The options are:
	//t(run time of Manbis)
	//a(number of digits of the output numbers),
	//e(the accuracy of the roots),
	//l(the left point of the initial interval),
	//r(the right point of the initial inteval)
	//q(the percentage of the roots that the user wants)
	//u(if u is equal to 1, then the program will produce 2 files with the
	//uncut subintervals. The file SubsByLength.txt contains the 
	//uncut subintervals, sorted by length and the file named 
	//SubsByXLeft.txt contains the uncut subintervals, sorted by 
	//the left endpoint.)
	//p(the user can change the confidence interval. User may define 
	//values 0.90, 0.99, 0.80 or 0.95 to obtain 1.644854, 2.575829, 1.281552 or 1.959964	
	//probability respectively)
	//n(set by user. If zero the roots 
	//will not be printed out in the Results.txt file, but only their number.)
	//s(the difference between two consecutive statistical levels)
	//c(current estimetions of Manbis)
	//h(help on the usage of Manbis)
	//The options l, r and q are mandatory to be initialized by the user.
	//The options a ,t ,e, u, p, n, s and c are optional. 
	//The default value for a, t, e, u, p, n, s and c are 8, 0, 4, 0, 0.95, 1 and 0 respectively.

	int long_index = 0;
	while ((opt = getopt_long(argc, argv,"a:t:e:l:r:q:u:p:n:s:c:h", long_options, &long_index )) != -1) {
		switch (opt) {

			case 'a' :  digits = atof(optarg);
						if ( digits < 0 ) {
  			   			    cerr << endl;
							cerr << "Wrong number of --digits (-a)" << endl;
							cerr << "Non negative value required" << endl;
							cerr << "The default value for digits is 8" << endl;
							cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
  			   			    cerr << endl;
							exit(EXIT_FAILURE);
						}


						break;
			case 't' :  totaltime = atof(optarg);
						if ( totaltime != 0 && totaltime != 1 ) {
							cerr << endl;
							cerr << "Wrong number of --time (-t)" << endl;
							cerr << "Valid values are only 0 or 1" << endl;
							cerr << "The default value for time is 0" << endl;
							cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
  			   			    cerr << endl;
							exit(EXIT_FAILURE);
						}

						break;
			case 'e' :  aqur = atof(optarg);
						e = 1;
						for ( int j = 0; j < aqur; j++ ) {
							e /= 10;
						}
						int *accuracy;
						if ( sizeof(accuracy) == 4 && aqur > 15 ) {
  			   			    cerr << endl;
							cerr << "Illegal --accuracy (-e). Maximum accuracy for 32-bit OS is 15." << endl;
							cerr << "Default value for --accuracy is 4 decimal digits" << endl;
							cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
  			   			    cerr << endl;
							exit(EXIT_FAILURE);
						}else if ( sizeof(accuracy) == 8 && aqur > 31 ) {
  			   			    cerr << endl;
							cerr << "Illegal --accuracy (-e). Maximum accuracy for 64-bit OS is 31." << endl;
							cerr << "Default value for --accuracy is 4 decimal digits" << endl;
							cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
  			   			    cerr << endl;
							exit(EXIT_FAILURE);
						}else if ( e <= 0 || e >= 1 ) {
							cerr << endl;
							cerr << "Illegal --accuracy (-e)" << endl;
							cerr << "A valid value can be a positive number." << endl;
							cerr << "The default value for --accuracy is 4 decimal digits" << endl;
							cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
							cerr << endl;
							exit(EXIT_FAILURE);
						}

						break;
			case 'l' : x1 = atof(optarg);//Saving the value of the left end of the interval
					   break;
			case 'r' : x2 = atof(optarg);//Saving the value of the right end of the interval
					   break;
			case 'q' : fraction = atof(optarg);
					   if ( fraction <= 0 || fraction > 100 ) {
  			   			   cerr << endl;
						   cerr << "Wrong value of --percentage (-q)" << endl;
						   cerr << "Valid values are greater that zero and less or equal that 100" << endl;
						   cerr << "Initialization of --percentage is mandatory" << endl;
						   cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
						   cerr << endl;
						   exit(EXIT_FAILURE);
					   }

					   break;
			case 'u' : files = atof(optarg);//Saving the information about the creation or not of the files 
					   //SubsByLength.txt and SubsByXLeft.txt
					   if ( files != 0 && files != 1 ) {
				   	   	   cerr << endl;
					   	   cerr << "Wrong value of --uncutfiles (-u)" << endl;
					   	   cerr << "Valid values are 0 or 1" << endl;
					   	   cerr << "The default value for --uncutfiles is 0" << endl;
   						   cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
					 	   cerr << endl;
   						   exit(EXIT_FAILURE);
					   }
					   break; 
			case 'p' : a = atof(optarg);//Set the required probability. The acceptable values of a
					   // are 0.05 and 0.1.
					   if ( a != 0.95 and a != 0.90 and a != 0.80 and a != 0.99 ) {
						cerr << endl;
						   cerr << "Wrong value of --probability (-p)" << endl;
						   cerr << "Valid values are 0.95 (default), 0.99, 0.90 or 0.80" << endl;
						   cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
       				   	   cerr << endl;
						   exit(EXIT_FAILURE);
					   }
					   break;   
			case 'n': nor = atof(optarg);//if nor=1 then the roots will be written in the file Results.txt, else only
						//their number will be printed
					   if ( nor != 0 && nor != 1 ) {
				   		   cerr << endl;
					   	   cerr << "Wrong value of --noroots (-n)" << endl;
					   	   cerr << "Valid values are 0 or 1" << endl;
					   	   cerr << "Default value of --noroots is 1 (print)" << endl;
   						   cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
					   	   exit(EXIT_FAILURE);
					   }

					   break;
			case 's' :statdiff = atof(optarg);// Sets the difference between two 
					  //consecutive statistical estimations.
						if ( statdiff < 0.1 || statdiff > 1.5 ) {
				   		   cerr << endl;
					   	   cerr << "Wrong value of --statdiff (-s)" << endl;
					   	   cerr << "Valid values are real numbers in the interval [0.1, 1.5]." << endl;
					   	   cerr << "Default value of --statdiff is 0.5." << endl;
   						   cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
					   	   exit(EXIT_FAILURE);
					   }
					   break;
			case 'c' :  currentest= atof(optarg);
						if ( currentest != 0 && currentest != 1 ) {
							cerr << endl;
							cerr << "Wrong number of --currentest (-c)" << endl;
							cerr << "Valid values are only 0 or 1" << endl;
							cerr << "The default value for currentest is 0" << endl;
							cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
  			   			    cerr << endl;
							exit(EXIT_FAILURE);
						}
						break;
			case 'h' :
			default: PrintUsage(); 
					 exit(EXIT_FAILURE);
		}//End of swich
	}//End of while
}//End of function SetCommandLineValues

//CheckingValuesOfUser checks if the arguments are correct.
void Parameters::CheckingValuesOfUser()
{

	if ( LLONG_MAX == x1 ) {
		cerr << endl;
		cerr << "Define the value of the --xleft (-l) of the initial interval." << endl;
		cerr << "Initialization of --xleft is mandatory." << endl;
		cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
		cerr << endl;
		exit(EXIT_FAILURE);
	}//If if(LLONG_MAX==x1) is true means xleft was not set by the user, 
	//so it prints an error message and terminates the process.

	if ( LLONG_MIN == x2 ) {
		cerr << endl;
		cerr << "Define the value of the --xright (-r) of the initial interval." << endl;
		cerr << "The initialization of --xleft is mandatory." << endl;
		cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
		cerr << endl;
		exit(EXIT_FAILURE);
	}//If if(LLONG_MIN==x2) is true means xright was not set by the user, 
	//so it prints an error message and terminates the process.

	if ( x1 >= x2 ) {
		cerr << endl;
		cerr << "Error. Left end of the interval is greater than the right end." << endl;
		cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
		cerr << endl;
		exit(EXIT_FAILURE);
	}//If the if(x1>=x2) is true the endpoints were wrong 
	//so it prints an error message and terminates the process.

	if ( fraction == -1 ) {
		cerr << endl;
		cerr << "Define the value of the --percentage of the roots (-q)." << endl;
		cerr << "Initialization of --percentage is mandatory." << endl;
		cerr << "For more help on the usage of Manbis type ./nameofexecutable.out -h or ./nameofexecutable.out --help" << endl;
		cerr << endl;
		exit(EXIT_FAILURE);
	}//If the if(fraction==-1) is true means fraction was not set by the user
	//so it prints an error message and terminates the process.

	if ( a == 0.90 ) {
		Z = 1.644854;
	} else if ( a == 0.99 ) {
		Z = 2.575829;
	} else if ( a == 0.80 ) { 
		Z = 1.281552;
	}
	fraction /= 100;


}//End of function CheckingValuesOfUser


//PrintUsage displays how to call the Manbis at command line. The usage is printed every time the user 
//makes an unacceptable call of Manbis.
void Parameters::PrintUsage() 
{
	cout << endl << endl;
	cout << "_______________________________USAGE OF Manbis_______________________________" << endl << endl;
	cout << "Usage: " << program_name << " --digits num --accuracy num --xleft num --xright num" << endl;
	cout << " --percentage num --uncutfiles num --probability num --noroots num" << endl;
	cout << "or with the following sortcuts respectively" << endl;
	cout << "Usage: " << program_name << " -a num -e num -l num -r num -q num -u num -p num -n num" << endl;
	cout << endl;
	cout << "--time(-t):       Run time of Manbis." << endl;
	cout << "                  Optional parameter. Valid values are 0 or 1. Default value = 0." << endl;	
	cout << "--digits(-a):     Number of digits that the output numbers will have." << endl;
	cout << "                  Optional parameter. Non-negative integer value. Default value = 8." << endl;
	cout << "--accuracy(-e):   Accuracy of the calculated roots." << endl;
	cout << "                  Optional parameter. Integer value greater than 0 and less than" << endl;
	cout << "                  16 for 32-bit OS, greater than 0 and less than 32 for 64-bit OS." << endl;
	cout << "                  Default value = 4." << endl;
	cout << "--xleft(-l):      Left endpoint of the initial interval." << endl;
	cout <<	"                  Mandatory parameter. Real value." << endl;
	cout << "--xright(-r):     Right endpoint of the initial interval." << endl;
	cout << "                  Mandatory parameter. Real value greater than the value of --xleft." << endl;
	cout << "--uncutfiles(-u): Flag that controls the creations of the files SubsByLength.txt " << endl;
	cout << "                  and SubsByXLeft.txt. Optional parameter." << endl;
	cout << "                  If we set this parameter equal to 1 then SubsByLength.txt will " << endl;
	cout << "                  contain  the uncut subintervals, sorted by length and" << endl;
	cout << "                  SubsByXLeft.txt will contain the uncut subintervals, " << endl;
	cout << "                  sorted by the left endpoint." << endl;
	cout << "                  Default value = 0 (no creation)." << endl;
	cout << "--percentage(-q): Percentage of the roots that the program will attempt calculate." << endl;
	cout << "                  Mandatory parameter." << endl;
	cout << "                  Real value greater than 0 and less or equal than 100."<<endl;
	cout << "--probability(-p):Probability of the confidence interval." << endl;
	cout << "                  Optional parameter. " << endl;
	cout << "                  Valid values are 0.95, 0.99, 0.90 or 0.80." << endl;
	cout << "                  Default value = 0.95" << endl;
	cout << "--noroots(-n):    Value 1 prints the roots in the file Result.txt." << endl; 
	cout << "                  Value 0 prints only their number." << endl;
	cout << "--statdiff(-s):   Difference between two consecutive statistical estimations" << endl;
	cout << "                  Optional parameter. " << endl;
	cout << "                  Valid values positive non zero real numbers in the interval [0.1, 1.5]." << endl;
	cout << "                  Default value = 0.5" << endl;
	cout << "--currentest(-c): Display the progress of the program" << endl;
	cout << "                  Optional parameter. Valid values are 0 or 1. Default value = 0." << endl;	
	cout << "--help(-h):       This help." << endl << endl;
	cout << "_______________________________USAGE OF Manbis_______________________________" << endl << endl << endl;
}//End of function PrintUsage
