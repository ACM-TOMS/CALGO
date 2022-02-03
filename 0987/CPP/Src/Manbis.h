#ifndef __Manbis_ALGORITHM_H__
#define __Manbis_ALGORITHM_H__
#include <cmath>//Library for mathimatical functions as abs(x) and sqrt(x).
#include "Parameters.h"
#include <iostream>//Library for input and output data.
#include <vector>//Since vectors are not a built-in part of the core language, the declaration is mandatory.
#include <map>//Since maps and multimaps are not a built-in part of the core language, the declaration is mandatory.
#include <fstream>//This library gives the potential to make files dynamically.
#include <cstdlib>//This library provides the opportunity to use EXIT_FAILURE.
#include <unistd.h>//Provides access to the POSIX operating system API.
#include <csignal>//This liblary gives the potential to the software to inform running processes of certain events 

#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <iomanip>


class Manbis
{
	// non-copyable
	Manbis(const Manbis &rhs);
	Manbis & operator=(const Manbis &rhs);

	private:
	 	clock_t t;//Holds the total time that Manbis needs, in order to give results.

	 	clock_t ct;//Holds the time that Manbis needs in order to compute a root at the current level.

	 	long int proots;//Holds the amount of the roots that the manbis computed at the previous level.

		long double (*function_)(long double );//Points to the user's formula from the file Functions.h 

		long double xleft, xright, maccuracy, mfraction, mstatdiff;//Creation of three real variables. xleft holds the left point of the
		// initial interval and xright holds the right point. maccuracy holds the accuracy of the roots
		// mfraction holds the percentage of the roots that the user wants.
		// mstatdiff is the wanted by the user difference between two consecutive statistical levels.


		int mdigits, mtotaltime, mcurrentest, mfiles, mnor;//mdigits holds the number of digits when the roots are printed.
		//if mtotaltime=1 then the run time of Manbis will be displayed at the end of the program.
		//if mcurrentest=1 then the currest estimation during the run of Manbis will be displayed at the screen
		//mfiles specifies if the user wants the creation of txt files SubsByLength.txt and SubsByXLeft.txt.
		//if mnor=1 then the roots will be written in the file Results.txt, else only
		//their number will be printed
		double mZ;//'mZ' is a parameter used in the statistics
		//For =0.05, mZ=1.959964.

		long double plower, pupper, nlower, nupper, avg, avg1;//Creation of six real variables.
		//"avg" is equal to fraction*[pow((nlower+nupper),2)] and "avg1" is the statistical estimation 
		//for the previous level.
		
		
	
		long long unsigned int nextcut;//Is a non-zero integer variable, 
		//which determines the level of the tree where the statistics
		//are currently computed.

		long double medial;//Holds the middle point of the subinterval that is being
		//divided at the current iteration.

		int sign, sign2;//Creation of two integer variables.
		//In each iteration they hold the sings of the subintervals  
		//(x1,medial) and (medial,x2), which are the two half's of the "parent"
		//interval $(x1,x2)$, where medial=(x1+x2)/2.
		//Variable "sign" contains the value 1 or the value -1 if the 
		//sings of the subinterval (x1,medial) are positive or negative respectively. 
		//Similarly the variable named "sing2" contains the value 1 or the value -1 if 
		//the sign of the subinterval (medial,x2) is positive or negative respectively.

		std::vector <long long unsigned int> h1;
		std::vector <std::vector<long long unsigned int> > even;
		std::vector <std::vector<long long unsigned int> > odd;
		//Creation of three vectors. The vector "h1" is one-dimensional
		//and consists of two slots. Both of these slots hold the number 99,
		//which is an arbitrary initialization. This vector is then used to create two two-
		//dimensional vectors called "even" and "odd" These vectors 
		//consist of exactly two columns while the number of the rows can change 
		//dynamically as the program proceeds. Vector "even" holds in every row the 
		//information about a subinterval, with an even number of roots. 
		//In its second column the number of times that the initial interval has been sub-divided,
		//is stored. Let m be that number. Then the length of such an interval is obviously l="initiallength"/2^m.
		//Now, in its first column the number of intervals, 
		//with length equal to l that lie between the left point of the initial interval 
		//and the left point of the current subinterval, is stored. Both of these numbers are non-zero integers
		//and can take values from 1 up to (2^64)-1. The actual value depends on the particular
		//system and library implementation. The maximum value allows storing 
		//a large number of subintervals in a very efficient way.
		//The vector "odd" has the same structure as vector "even", the only difference being  
		//the information that is contained in "odd" which now concerns the subintervals
		//with odd number of roots. The use of two separate vectors helps to give priority 
		//to the subintervals, which certainly have roots (vector "odd").
		//For more information see http://www.cplusplus.com/forum/beginner/12409/
	
		int ar, per;//The integer variables named "ar" and "per" are assistant 
		//variables for proper initialization of both vectors "even" and "odd".
	
		std::multimap <long long unsigned int,long long unsigned int> lengtheven;
		std::multimap <long long unsigned int,long long unsigned int> lengthodd;
		std::multimap <long long unsigned int,long long unsigned int>::iterator itmultimap;
		//Creation of two multimaps, which store two integer numbers.
		//These numbers can take values from 1 up to (2^64)-1 (with a 64-bit processor). 
		//Multimap "lengtheven" contains pairs of numbers. Each pair gives the position
		//of the subinterval that will be sub-divided in the future. The "lengtheven" holds
		//the positions (that is, the respective slots in the two dimensional vector "even")
		//for subintervals that are even. Specifically, the first number of the pair 
		//depicts how many times the initial interval has been divided. From this number
		//the length of the current subinterval can be calculated. Furthermore, the second 
		//number of the pair depicts the position of the important information 
		//concerning the current subinterval in the vector "even".
		//The elements in a multimap are always sorted according to the key following a specific strict 
		//weak ordering criterion indicated by its internal comparison object (of type Compare).
		//Multimap containers are generally slower than unordered_multimap containers to 
		//access individual elements by their key, but they allow the direct iteration on subsets
		//based on their order. Apparently $multimaps$ are typically implemented as binary search
		//trees, so the wast time to store a new pair of number or to find a pair of number is 
		//O(log(n)). Furthermore the extent of these forms is changing dynamically.
		//In view of the above properties in every iteration the subinterval that will be 
		//selected for examination can be found very quickly. To clarify this  point, the 
		//next subinterval for examination is the one with the largest length.
		//This derives from the fact that the first pair contains the smallest key 
		//number, so this subinterval is the interval that was created with the fewer cuts from 
		//the initial interval, i.e. it is the subinterval with the largest length.  
		//For the multimap named "lengthodd" the same properties hold, with the difference that it stores
		//the position for the subintervals in the vector odd. 
		//The multimaps "lengtheven" and "lengthodd" give priority to the examinations
		//of the odd subintervals, consulting the multimap "lengthodd", and when the 
		//"lengthodd" is emptied, then the program will start the sub-division of the
		//subintervals with an even number of roots, consulting the multimap "lengtheven".
		//If both "lengthodd" and "lengtheven" are empty the program terminates.
		//The iterator "itmultimap" is used for input and output information
		//both to "lengtheven" and "lengthodd".
		//For more information see http://www.cplusplus.com/reference/map/multimap/


		std::map<long double, long long unsigned int> myroots;
		std::map<long double, long long unsigned int>::iterator itmyroots;// The map
		//"myroots" stores the roots that the program is computing. As for multimaps, 
		//maps are composed of pairs of numbers. Their extents can change dynamically, the 
		//pairs are in order. As opposed to multimaps, maps don' t allow any of the keys to 
		//exist more than once. So even if the program calculates a specific root more than once, the 
		//map "myroots" will hold the value of this root only once. The roots are stored in an 
		//increasing order. The first number of the pair(key) is the value of the roots and the 
		//second number is the iteration of the program that the roots have been computed.
		//The iterator named "itmyroots" is used for input and output information
		//with the map "myroots".
		//For more information see http://www.cplusplus.com/reference/map/map/?kw=map


		std::multimap <long double, long double> subintervals;
		std::multimap <long double, long double>::iterator itsubintervals;// The multimap
		//"subintervals" stores the subintervals that the program has not examinated yet.
		//The subintervals are firstly stored in the multimap "subintervals" by length and then printed out
		//in the file SubsByLength.txt. Then the program cleans the multimap "subintervals"
		//and then the subintervals are stored again in "subintervals", but now
		//sorted by their left endpoint. Then they are printed out in the file SubsByXLeft.txt.  
	
		long double initialength, initialeft, tinylength;//Creation of three Real variables.
		//Variable "initialength" holds the length of the initial interval.
		//Variable "initialeft" holds the left point of the initial interval. 
		//Variable "tinylength" holds the smallest length 
		//that a subinterval can have. If an even subinterval has length smaller 
		//than "tinylength", then the program terminates the examination of this subinterval 
		//and erases it. If an odd subinterval has length smaller than "tinylength", then
		//the program calculates the functional value of the subinterval's middle point, 
		//and if the points (middle-e) and (middle+e) have opposite signs, then the value of 
		//the variable "middle" is stored in the map "myroots" and the program erases
		//the current subinterval.

		long long unsigned int cut, htree;//Creation of three positive variables that 
		//may assume values from 0 up to (2^64)-1 (for 64-bit processor). The variable 
		//"cut" contains the depth of the subdivision tree where the two halfs of the current
		//subinterval lie.
		//The variable "htree" holds the maximum possible depth the subdivision tree may have.



		long long unsigned int numberofit;//The variable "numberofit" holds the 
		//current number of iterations that the program has done up to the current point.
		//By iteration we mean a single subdivision of a subinterval (with either same or opposite signs) 
		//This variable is also used to store the roots in the map "myroots".

		double power, power1;//Variable "power" holds the power of the number 2,
		//which the current  iteration demands lets say equal to 2^k. Then "power1"
		// holds the value 2^(k-1).

		long long unsigned int position, interpolations; //"position" is an integer variable
		//which is holding the position in vector "odd" of the subinterval that will be examined
		//at every iteration.
		//If the initial interval has been sub-divided m times to obtain the current subinterval
		//then "interpolations" is an integer variable, which holds the number of subintervals
		//with length equal to initialength/(2^m) interpolates from the left point
		//of the initial interval up to the left point of the current subinterval.
	
		std::vector <long long unsigned int> h2;
		std::vector <std::vector <long long unsigned int> > evenandodd;///Creation of 
		//two vectors. The vector named "h2" is one-dimensional
		//and is consisted by two slots. Both of these slots are holding the number 0.
		//The creation of this vector is made with the view 
		//to construct the two-dimensional vectors named "evenandodd".
		//This vector has htree+1 rows. The total sum of each row depicts the number of the 
		//subintervals that are having equal length. In each row, the first number is the number 
		//of the subintervals that contain even number of roots and apparently the second number 
		//holds the number of the subintervals that contain an odd number of roots. 
		//For example if the initial interval has been sub-divided 2 times, then the total sum 
		//of the second row of "evenandodd" (starting from zero), eventually will be equal to the 
		//value of 4.

		long long unsigned int number_even, number_odd;//These variables are holding the number
		//of the even and the odd subintervals of the same length, that the Manbis produces in every iteration

		long long unsigned int helphtree;//This variable is assistive for the calculation
		//of the value of the variable named  tinylength. 

		int statistic_stopping_criterion;//This variable holds if the stopping criterion has fulfilled or not.

		long long unsigned int unit;//This variable is assisteive for the calculation of the power 2^i


	public:
	
		//Manbis();
		Manbis(const Parameters &aa, long double (*function_)(long double ));//Constractor

		
		//This function sets some vital initialazations 
		//for the program.
		void PrimordialStep();
		
		//The following function named CleanRoots, 
		//checks if the new root has been calculated 
		//in a previous iteration.
		void CleanRoots(std::map<long double, long long unsigned int>& myroots, long double e, long double medial);

		//This function calculates the statistic, 
		//and if the respective stopping criteria hold, 
		//first it creates the txt files SubsByLength.txt 
		//and SubsByXLeft.txt if the user
		//has asked them and then the function terminates the
		//main loop of the program
		void StatisticStopingCriterion();
		
		//This function cuts the odd uncut subintervals.
		void CutOdds();
	
		//This function cuts the even uncut subintervals.
		void CutEven();
		
		//This function print out the calculated roots.
		void PrintResults();

		int GetExistOdd();

		long long unsigned int GetNumberOfOdd();

		long long unsigned int GetNumberOfEven();

		//This function estimates how much time the Manbis 
		//will take to calculate all the roots in the worst case senario.
		//and print out some usefull information on the user's screen.
		void PrintChoices();

		//This function is the that Manbis calculates the roots
		//using all the above mentioned functions 
		void FindRoots();

	};//End of class Manbis
int Bolzano(long double x1,long double x2, long double (*function_)(long double ) );

void  UserInterruption(int sig);//This function terminates the software if the user interrupt the software
// by the use of Ctrl + C and saves the roots, which have been calculated until the interruption. 

#endif //__Manbis_ALGORITHM_H__
