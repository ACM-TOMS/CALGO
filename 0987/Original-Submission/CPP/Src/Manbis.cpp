#include "Manbis.h"
using namespace std;
volatile sig_atomic_t flag = 0;


//Constractor
Manbis::Manbis(const Parameters &parameters, long double (*f)(long double )){
	function_ = f;//points to the user's formula from the file Functions.h 

	//user's initializations
	xleft = parameters.x1;//left point of the initial interval.
	xright = parameters.x2;//right point of the initial interval.
	mdigits = parameters.digits;//number of digits when the roots are printed.
	mtotaltime = parameters.totaltime;//if mtotaltime == 1 then the run time of Manbis will be displayed 
	//at the end of the program.
	mcurrentest = parameters.currentest;//if mcurrentest == 1 then the current estimations of Manbis will be displayed 
	//on the screen during the run of the software.
	//lolamcurrenttime = parameters.currenttime;//if holds the run time that Manbis needs to compute a root at the current level. 
	//lolaThis will be displayed at the last column an the screen, if the mcurrentest == 1 
	maccuracy = parameters.e;//holds the accuracy of the roots.
	mfraction = parameters.fraction;//holds the percentage of the roots that the user wants.
	mfiles = parameters.files;//Specifies if the user wants the creation of txt files SubsByLength.txt
	//and SubsByXLeft.txt.
	mZ = parameters.Z;//mZ is a parameter used in the statistics. For a == 0.05, mZ == 1.959964.
	mnor = parameters.nor;//if mnor == 1 then the roots will be written in the file Results.txt, else only
	//their number will be printed
	mstatdiff = parameters.statdiff;//Holds the wanted by the user difference between two consecutive 
	//statistical levels.
	//user's initializations

	//Initializations
	number_even = 0;
	number_odd = 0;
	nextcut = 1;
	unit = 1;
	h1.push_back(99);
	h1.push_back(99);
	even.push_back(h1);
	odd.push_back(h1);
	statistic_stopping_criterion = 0;
	ar = 0;
	per = 0;
	cut = 2;
	numberofit = 1;
	initialeft = xleft;
	initialength = xright - xleft;
	htree = ceil( log(initialength / maccuracy) / log(2) ) + 1;
	if ( htree < 64 ) {
		tinylength = initialength / (unit << htree);
	} else {
		helphtree = htree;
		tinylength = initialength;
		while ( helphtree >= 64 ) {
			tinylength /= (unit << 63);
			helphtree -= 63;
		}
		tinylength /= (unit << helphtree);
	}

	h2.push_back(0);
	h2.push_back(0);
	for( long long unsigned int i = 0; i <= (htree + 1); i++ ) { evenandodd.push_back(h2);}
	//End of initializations

}//End of Constractor



//This function sets some vital initializations 

void Manbis::PrimordialStep(){
	medial = (xleft + xright) / 2;//Calculate the middle point of the initial interval (x1,x2).
	//________________________________________________________________________________________
	//If "medial" is a root of the function, then store its value 
	// in the map "myroots".
	if ( Bolzano(medial - maccuracy, medial + maccuracy, function_) <= 0) {
		myroots.insert ( pair<long double, long long unsigned int>(medial, numberofit) );
		numberofit++;
	}
	//________________________________________________________________________________________

	//________________________________________________________________________________________
	//Find the sign of the initial interval (x1 - x2) using the function from the header file 
	//FunctionSubs.h and update vector evenandodd 
	if ( Bolzano(xleft, xright, function_) == 1 ) {
		evenandodd[0][0] += 1;}
	else{
		evenandodd[0][1] += 1;
	}
	//________________________________________________________________________________________

	//________________________________________________________________________________________
	//Find the sign of the interval (x1, medial) and update the appropriate data structures. 
	//This part of the program handles 2 different cases. Only one of them will be true. 
	//The cases are:
	//1.(x1, medial) has opposite signs at its endpoints. In this case the first odd subinterval
	//is found so variable "per" gets the value 1 and sets "odd", "lengthodd" and "evenandodd"
	//are updated
	//2.(x1, medial) has same signs at its endpoints. In this case the first even subinterval
	//is found so variable "ar" gets the value 1 and sets "even", "lengtheven" and "evenandodd"
	//are updated
	sign = Bolzano(xleft, medial, function_);
	if ( sign == -1 ) {
		per = 1;
		odd[0][0] = 0;
		odd[0][1] = 1;
		lengthodd.insert( pair<long long unsigned int, long long unsigned int>(1, 0) );
		evenandodd[1][1]+=1;
	}

	if ( sign == 1 ) {
		ar = 1;
		even[0][0] = 0;
		even[0][1] = 1;
		lengtheven.insert( pair<long long unsigned int, long long unsigned int>(1, 0) );
		evenandodd[1][0] += 1;
	}
	//________________________________________________________________________________________

	//________________________________________________________________________________________
	//Find the sign of the interval (medial, x2) and inform the appropriate data structures. 
	//This part of the program handles 4 different cases. Only one of them will be true.
	//The four cases are the following:

	//1.(x1, medial) > 0 (same signs) and (medial, x2) < 0 (opposite signs). In this case the first
	//odd subinterval is found so the variable "per" gets 
	//the value 1 and sets "odd", "lengthodd" and "evenandodd" are upadated.

	//2.(x1, medial) < 0 (opposite signs) and (medial,x2) < 0 (opposit signs). In this case
	//2 odd subintervals have emerged so the variable "per" is already 1 
	//Second slot of "odd" and sets "lengthodd" and "evenandodd" must be updated. 

	//3.(x1,medial) < 0 (opposit signs) and (medial,x2) > 0 (same signs).  In this case the first
	//even subinterval is found so the variable so the variable "ar" gets 
	//the value of 1 and sets "even", "lengtheven" and "evenandodd" are updated.

	//4.(x1,medial) > 0 (same signs) and (medial,x2) > 0 (same signs). In this case
	//2 even subintervals have emerged so the variable "ar" is already 1 
	//Second slot of "even" and sets "lengtheven" and "evenandodd" must be updated. 
	sign2 = Bolzano(xright, medial, function_);
	if ( sign != -1 && sign2 == -1 ) {
		per = 1;
		odd[0][0] = 1;
		odd[0][1] = 1;
		lengthodd.insert( pair<long long unsigned int, long long unsigned int>(1, 0) );
		evenandodd[1][1] += 1;
	}

	if(sign == -1 && sign2 == -1){
		h1[0] = 1;
		h1[1] = 1;
		odd.push_back(h1);
		lengthodd.insert( pair<long long unsigned int,long long unsigned int>(1, 1) );
		evenandodd[1][1] += 1;
	}

	if ( sign != 1 && sign2 == 1 ) {
		ar = 1;
		even[0][0] = 1;
		even[0][1] = 1;
		lengtheven.insert( pair<long long unsigned int, long long unsigned int>(1, 0) );
		evenandodd[1][0] += 1;
	}


	if ( sign == 1 && sign2 == 1 ) {
		h1[0] = 1;
		h1[1] = 1;
		even.push_back(h1);
		lengtheven.insert( pair<long long unsigned int,long long unsigned int>(1, 1) );
		evenandodd[1][0] += 1;
	}
	//________________________________________________________________________________________

	//________________________________________________________________________________________
	//Calculate the statistics variables for zero cuts of the initial interval:
	//"plower", "pupper", "nlower", "nupper", "avg".
	//Also calculate the number of roots, which have been found up to this point
	plower = evenandodd[0][1] - mZ * sqrt( evenandodd[0][1] * (1 - evenandodd[0][1]) );
	pupper = evenandodd[0][1] + mZ * sqrt( evenandodd[0][1] * (1 - evenandodd[0][1]) );
	nlower = log(1 - 2 * plower) / log(1 - 2);
	nupper = log(1 - 2 * pupper) / log(1 - 2);
	avg = mfraction * (nlower + nupper) / 2;
	//________________________________________________________________________________________
	//-------------------------------END OF PRIMORDIAL STEP------------------------------------- 
	number_even = lengtheven.size();
	number_odd = lengthodd.size();


}//End of function PrimordialStep.


//CleanRoots checks if the new root has been calculated 
//in a previous iteration.
void Manbis::CleanRoots(map<long double, long long unsigned int>& myroots, long double e, long double medial){
	//function consists of two basic cases, as follows:
	//1. if the map named "myroots" contains only 2 roots.
	//2. if the map named "myroots" contains more than 2 roots. 

	long double root1, root2, root3;//Creation of three assist variables.
	map<long double, long long unsigned int>::iterator itmyroots;
	//"itmyroots" is used for input and output information to the map "myroots".

	//Starting point of first case.
	if ( myroots.size() == 2 ) {
		itmyroots = myroots.begin();
		root1 = itmyroots -> first;
		itmyroots++;
		root2 = itmyroots -> first;
		if ( root2 - root1 <= e ) {
			if ( abs( (*function_)(root2) ) > abs( (*function_)(root1) ) ) {
				myroots.erase(itmyroots);
			}//End of if(abs(function_(root2))>abs(function_(root1))).
			else{
				itmyroots--;
				myroots.erase(itmyroots);
			}//End of else.
		}//End of if(root2-root1<=e).
	}//End of if(myroots.size()==2).
	//Ending point of the first case.

	//Starting point of second case.
	//The second case consists of three cases, as follows:
	//a. if the new root is the first element of the map "myroots".
	//b. if the new root is the last element of the map "myroots".
	//c. if the new root is a middle element of the map "myroots".

	if ( myroots.size() > 2 ) {
		//Starting point of case "a".
		//Only one comparison has to be made 
		//with the second element of the map "myroots".
		itmyroots = myroots.find(medial);
		if(itmyroots == myroots.begin()){
			root1 = itmyroots -> first;
			itmyroots++;
			root2 = itmyroots -> first;
			if ( root2 - root1 <= e ) {
				if ( abs( (*function_)(root2) ) > abs( (*function_)(root1) ) ) {
					myroots.erase(itmyroots);
				}//End of if(abs(function_(root2))>abs(function_(root1))).
				else{
					itmyroots--;
					myroots.erase(itmyroots);
				}//if(abs((*function_)(root2))>abs((*function_)(root1)))
			}//End of if(root2-root1<=e). 
		}//End of if(itmyroots==myroots.begin()).
		//End of case "a" .

		//Starting point of the case "b".
		//Only one comparison has to be made 
		//with the previous element of the map "myroots".
		else if ( itmyroots == --myroots.end() ) {
			root2 = itmyroots -> first;
			itmyroots--;
			root1 = itmyroots -> first;
			if ( root2 - root1 <= e ) {
				if ( abs( (*function_)(root1) ) > abs( (*function_)(root2) ) ) {
					myroots.erase(itmyroots);
				}//End of if(abs(function_(root1))>abs(function_(root2))).
				else{
					itmyroots++;
					myroots.erase(itmyroots);
				}//End of if(abs((*function_)(root1))>abs((*function_)(root2)))
			}//End of if(root2-root1<=e).
		}//End of else if(itmyroots==--myroots.end()).
		//End of the case "b" .

		//Starting point of the case "c" .
		else if ( (itmyroots != --myroots.end())  &&  (itmyroots != myroots.begin()) ) {
			itmyroots--;
			root1 = itmyroots -> first;
			itmyroots++;
			root2 = itmyroots -> first;
			itmyroots++;
			root3 = itmyroots -> first;
			itmyroots--;
			if ( (root2 - root1 <= e) && (root3 - root2 <= e) ) {
				if ( abs( (*function_)(root2)) <= abs( (*function_)(root1) ) && abs( (*function_)(root2) ) <= abs( (*function_)(root3) ) ) {
					itmyroots--;
					myroots.erase(itmyroots);
					itmyroots = myroots.find(root3);
					myroots.erase(itmyroots);   
				}//End of if(abs(function_(root2))<=abs(function_(root1)) && abs(function_(root2))<=abs(function_(root3))).
				else if (abs ( (*function_)(root1)) <= abs((*function_)(root2) ) && abs( (*function_)(root1)) <= abs( (*function_)(root3) ) ) {
					myroots.erase(itmyroots);
					itmyroots = myroots.find(root3);
					myroots.erase(itmyroots);
				}//End of if((root2-root1<=e) && (root3-root2<=e)).
				else{
					myroots.erase(itmyroots);
					itmyroots = myroots.find(root1);
					myroots.erase(itmyroots);   
				}//End of else.
			}//End of if((root2-root1<=e) && (root3-root2<=e)).
			else if ( root2 - root1 <= e ) {
				if ( abs( (*function_)(root2) ) <= abs( (*function_)(root1) ) ) {
					itmyroots--;
					myroots.erase(itmyroots);
				}//End of if(abs(function_(root2))<=abs(function_(root1))).
				else{
					myroots.erase(itmyroots);
				}//End of else.
			}//End of else if(root2-root1<=e).
			else if ( root3 - root2 <= e ) {
				if ( abs( (*function_)(root2) ) <= abs( (*function_)(root3) ) ) {
					itmyroots++;
					myroots.erase(itmyroots);
				}//End of if(abs(function_(root2))<=abs(function_(root3))).
				else{
					myroots.erase(itmyroots);
				}//End of else
			}//End of else if(root3-root2<=e)
		}//End of if that checks if the new root is between myroots.
		//End of the case "c" .

	}//End of if((itmyroots!=--myroots.end())  &&  (itmyroots!=myroots.begin())).
	//End of second case.

}//End of function CleanRoots.

//This function calculates the statistics, 
//and if the stopping criteria are fulfilled, 
//creates the txt files SubsByLength.txt 
//and SubsByXLeft.txt if the user
//has asked for them and terminates the
//main loop of the program
void Manbis::StatisticStopingCriterion(){
	//_________Creation of the txt files____________
	ofstream statistika;
	ofstream coutuncut;
	ofstream coutuncut2;
	//_________Creation of the txt files____________

	long double x1,x2;//These variables hold the left and
	//the right point of the current odd subinterval.

	long long unsigned int numberofroots;//Holds the number of roots at the current iteration.
	//*********************Start of the statistics****************************//
	if ( (cut - 1) >= 64 ) {
		power = 1;
		helphtree = cut - 1;
		while ( helphtree >= 64 ) {
			power *= ( unit << 63 );
			helphtree -= 63;
		}
		power *= ( unit << helphtree );

	}else{
		power = ( unit << (cut - 1) );
	}//power of 2^(cut-1) depending on the word size of the OS. "cut" is holding the number
	//of the subintervals sub-divisions. The statistic criterion, is checked 
	//for the current level.  
	if ( power == (evenandodd[cut-1][0] + evenandodd[cut-1][1]) ) {//If all the subintervals
		//of the previous cut have been obtained, then the statistics formulas can be calculated.
		if ( nextcut == cut - 1 ) {//Make sure that the statistic of 
			//each level will be calculated only once.
			//_______________________________________________________________________
			nextcut += 1;
			power1 = power / 2;
			plower = ( evenandodd[cut-2][1] - mZ * sqrt( (evenandodd[cut-2][1] * (power1 - evenandodd[cut-2][1]) ) / power1 ) ) / power1;
			pupper = ( evenandodd[cut-2][1] + mZ * sqrt( (evenandodd[cut-2][1] * (power1 - evenandodd[cut-2][1]) ) / power1 ) ) / power1;
			nlower = log(1 - 2 * plower) / log(1 - 2 * (1 / power1));
			nupper = log(1 - 2 * pupper) / log(1 - 2 * (1 / power1));
			avg1 = mfraction * (nlower + nupper) / 2;

			plower = (evenandodd[cut-1][1] - mZ * sqrt( (evenandodd[cut-1][1] * (power - evenandodd[cut-1][1]) ) / power) ) / power;
			pupper = (evenandodd[cut-1][1] + mZ * sqrt( (evenandodd[cut-1][1] * (power - evenandodd[cut-1][1]) ) / power) ) / power;
			nlower = log(1 - 2 * plower) / log(1 - 2 * (1 / power));
			nupper = log(1 - 2 * pupper) / log(1 - 2 * (1 / power));
			avg = mfraction * (nlower + nupper) / 2;

			if ( mcurrentest == 1 ) {

				proots=myroots.size()-proots;

				if (nextcut-1==1 ){
					cout << "\033[1;31m   Current      \033[0m"<<" | "<<"\033[1;32m   Requested   \033[0m"<<" | "<<"\033[1;33m   Discovered   \033[0m"<<" | "<<"\033[1;34m   Discovered   \033[0m"<<" | "<<"\033[1;36m   Avg time of computing    \033[0m"<<endl;
					cout << "\033[1;31m   estimation   \033[0m"<<" | "<<"\033[1;32m   #roots      \033[0m"<<" | "<<"\033[1;33m   #roots       \033[0m"<<" | "<<"\033[1;34m   fraction     \033[0m"<<" | "<<"\033[1;36m   100 roots    \033[0m"<<endl<<endl;

				}
				ct = clock() - ct;

				if(proots>0 && !isnan(proots) && !isinf(proots)){
					if ( ceil( avg ) >= myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) !=1 && isnan( ceil ( avg ) ) != 1 ){
						cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" << ceil( ( nlower + nupper) / 2 ) << '\t' << '\t' << "\033[1;32m" << ceil( avg ) << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t' << '\t' << "\033[1;34m" ;
						cout<< ceil(( evenandodd[cut-1][1] / ( (nlower + nupper) / 2) ) * 100) << "%" << '\t' << '\t' << "\033[1;36m" << (((float)ct) / CLOCKS_PER_SEC)/proots *100<<"sec"<< "\033[1;37m";
 						cout.flush();
	 				}else if(ceil( avg ) < myroots.size() && ceil( ( nlower + nupper) / 2 ) < myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) != 1 && isnan( ceil( avg ) ) != 1 ){
 						cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" << myroots.size() << '\t' << '\t' << "\033[1;32m" << myroots.size()  << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t'  << '\t' << "\033[1;34m" ;
 						cout<< 100 << "%" << '\t' << '\t' << "\033[1;36m" << (((float)ct) / CLOCKS_PER_SEC)/proots *100 <<"sec"<< "\033[1;37m";
 						cout.flush();
 					}else if(ceil( avg ) < myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) != 1 && isnan( ceil( avg ) ) != 1 ){
 					 	cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" << ceil( ( nlower + nupper) / 2 ) << '\t' << '\t' << "\033[1;32m" << myroots.size()  << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t'  << '\t' << "\033[1;34m" ;
 				 		cout << ceil(( evenandodd[cut-1][1] / ( (nlower + nupper) / 2) ) * 100) << "%" << '\t' << '\t' << "\033[1;36m" << (((float)ct) / CLOCKS_PER_SEC)/proots *100<<"sec"<< "\033[1;37m";
 						cout.flush();
	 				}	
 				}else {
 					if ( ceil( avg ) >= myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) !=1 && isnan( ceil ( avg ) ) != 1 ){
						cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" << ceil( ( nlower + nupper) / 2 ) << '\t' << '\t' << "\033[1;32m" << ceil( avg ) << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t' << '\t' << "\033[1;34m" ;
						cout<< ceil(( evenandodd[cut-1][1] / ( (nlower + nupper) / 2) ) * 100) << "%" << '\t' << '\t' << "\033[1;37m" << string(30, '-') ;
 						cout.flush();
	 				}else if(ceil( avg ) < myroots.size() && ceil( ( nlower + nupper) / 2 ) < myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) != 1 && isnan( ceil( avg ) ) != 1 ){
 						cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" << myroots.size() << '\t' << '\t' << "\033[1;32m" << myroots.size()  << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t'  << '\t' << "\033[1;34m" ;
 						cout<< 100 << "%" << '\t' << '\t' << "\033[1;37m" << string(30, '-');
 						cout.flush();
 					}else if(ceil( avg ) < myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) != 1 && isnan( ceil( avg ) ) != 1 ){
 					 	cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" << ceil( ( nlower + nupper) / 2 ) << '\t' << '\t' << "\033[1;32m" << myroots.size()  << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t'  << '\t' << "\033[1;34m" ;
 				 		cout << ceil(( evenandodd[cut-1][1] / ( (nlower + nupper) / 2) ) * 100) << "%" << '\t' << '\t' << "\033[1;37m" << string(30, '-');
 						cout.flush();
 					}

 				}//end of if(proots>0 && !isnan(proots) && !isinf(proots))

				ct = clock();
 				proots=myroots.size();

			}
			
			
			//At this point MANBIS is checking if the current statistic 
			//estimation is close to the estimation of the previous level.
			if ( abs(avg - avg1) < mstatdiff * avg ) {

				//By now, the statistics have been calculated.
				//_______________________________________________________________________
				if ( myroots.size() < evenandodd[cut-1][1] ) {
					numberofroots = myroots.size();
				}else{
					numberofroots = evenandodd[cut-1][1];
				}
				//----------Start of the stopping criterion of the statistics formulations------------
				if ( numberofroots > avg && avg > 0 && nlower >= 0 && nupper >= 0 && plower >= 0 && pupper >= 0 ) {//if  

					//condition holds, the program produces files SubsByLength.txt
					//and SubsByXLeft.txt, empties both multimaps "evenlength" and "oddlength"
					//If both are empty the main step in main.cpp terminates
					//
					//_________Start of creation of the txt files____________
					if ( mfiles == 1 ) {//If the user wants the files SubsByLength.txt and SubsByXLeft.txt
						//When the program terminates, the files 
						//SubsByLength.txt and SubsByXLeft.txt will 
						//hold the uncut subintervals, which will 
						//be printed out with the number of 
						//digits that the user has defined. 

						//Opening two files SubsByLength.txt and SubsByXLeft.txt
						//The file named "SubsByLength.txt", at the end of the program
						//will contain the subintervals, which were not examined, sorted by length.
						//The file named "SubsByXLeft.txt", at the end of the program,
						//will contain the subintervals, which were not examined, sorted by their left point.
						//These two files will be created only if the user wants.
						coutuncut.open("SubsByLength.txt");
						coutuncut2.open("SubsByXLeft.txt");
						coutuncut.precision(mdigits);
						coutuncut2.precision(mdigits);
						//_________End of creation of the txt files____________
						//_______________________________________________________________________
						//This part of the program outputs the even subintervals sorted by length in the
						//file SubsByLength.txt
						int itsubs = 1;
						coutuncut << "The uncut even subintervals are " << lengtheven.size() << endl;
						coutuncut << "Sorted by length" << endl;
						for ( (itmultimap = lengtheven.begin()); itmultimap != lengtheven.end(); itmultimap++ ) {
						
							//Calculation of the number 2^(even[itmultimap->second][1])
							//for the OS's word size
							if ( (even[itmultimap -> second][1]) >= 64 ) {
								power = 1;
								helphtree = even[itmultimap -> second][1];
								while ( helphtree >= 64 ) {
									power *= (unit << 63);
									helphtree -= 63;
								}
							power *= (unit << helphtree);
							}else{
								power = ( unit << (even[itmultimap -> second][1]) );
							}

							x1 = initialeft + even[itmultimap -> second][0] * (initialength / power);
							x2 = x1 + initialength / power;
							subintervals.insert( pair<long double,long double>(x1, x2) );
							coutuncut << itsubs << ")=" << "\t" << "(" << x1 << " , " << x2 << ") " << "\t" << "length=" << x2 - x1 << endl;
							itsubs++;
						}
						//_______________________________________________________________________

						//_______________________________________________________________________
						//This part of the program outputs the even subintervals sorted by xleft in the
						//file SubsByXLeft.txt
						itsubs = 1;
						coutuncut2 << "The uncut odd subintervals are " << subintervals.size() << endl;
						coutuncut2 << "Sorted by xleft" << endl;
						for ( (itsubintervals = subintervals.begin()); itsubintervals != subintervals.end(); itsubintervals++ ) {
							x1 = itsubintervals -> first;
							x2 = itsubintervals -> second;
							coutuncut2 << itsubs << ")=" << "\t" << "(" << x1 << " , " << x2 << ") " << "\t" << "length=" << x2 - x1 << endl;
							itsubs++;
						}
						subintervals.clear();
						//_______________________________________________________________________
						coutuncut.close();//Closing  file SubsByLength.txt 
						coutuncut2.close();//Closing file SubsByXLeft.txt                     
					}//End of if(mfiles==1)

					lengthodd.clear();
					lengtheven.clear();  
					statistic_stopping_criterion = 1;  
				}//----------End of the stopping criterion of the statistics formulations------------
			}//End of if(nextcut==cut-1)
		}//End of if(abs(avg-avg1)<statdiff*avg)
	}//End of if(power==(evenandodd[cut-1][0]+evenandodd[cut-1][1]))
	//*********************End of the statistics****************************//

}//End of function StatisticStopingCriterion;

//This function divides the odd uncut subintervals.
void Manbis::CutOdds(){//Divides the odd subinterval that is stored
	//in the top element of "lengthodd"
	long double x1, x2;//These variables hold the left and
	//the right point of the current odd subinterval.

	//___________________________Start the "do while" of the odd subintervals____________________________
	while ( lengthodd.size() != 0 && per == 1 ) {//Starting point of the while for odd subintervals.
		//The variable "per", denotes the existence or not of an odd subinterval. 
		//Required because when "lengthodd" is created, it gets random values
		//so even if its  size is 1, no odd subintervals exist. 
		//If "per" is zero, then there is no odd subinterval

		itmultimap = lengthodd.begin();//"itmultimap" holds the location
		//of the subinterval that will be examined at the current iteration.

		cut = itmultimap -> first + 1;//"cut" contains the depth of the subdivision tree where the two halfs of the current
		//subinterval lie.

		position = itmultimap -> second;//"position" holds the position in vector "odd" 
		//of the current odd subinterval. 

		lengthodd.erase(itmultimap);//Erase the current subinterval from "lengthodd"
		//to avoid examinating it for a second time.


		interpolations = odd[position][0];//"interpolations" is holding how many 
		//subintervals with length equal to  initialength/(2^(cut-1)) lie between the 
		//left point of the initial interval and the left point of the current 
		//odd subinterval.

		StatisticStopingCriterion();
		//If the statistics are not ready to be printed out or the stopping criterion of the 
		//statistics doesn't hold, then the program starts to examinate the largest odd subinterval
		//at the following lines of program. 

		if ( statistic_stopping_criterion == 0 ) {

			//_____________________________________________________________________________________
			//Calculation of the number power=2^(odd[position][1])
			//for the word size of the OS.
			if ( (odd[position][1]) >= 64 ) {
				power = 1;
				helphtree = odd[position][1];
				while ( helphtree >= 64 ) {
					power *= (unit << 63);
					helphtree -= 63;
				}
				power *= (unit << helphtree);
			}else{
				power = (unit << (odd[position][1]) );
			}
			//_____________________________________________________________________________________

			x1 = initialeft + odd[position][0] * (initialength / power);//"x1" holds the left endpoint of
			//the current odd subinterval.

			x2 = x1 + initialength / power;//"x2" holds the right endpoint of the current odd subinterval.

			medial = (x1 + x2) / 2;//"medial" holds the middle point of the current odd
			//subinterval (x1,x2).

			//__________________________________________________________________________________________
			//This part of the program checks if the half of the current odd subinterval is less than
			//the length of "tinylength". If it is smaller
			//then the examination of the current subinterval stops, and the value of the 
			//variable "medial" is stored in the map "myroots".
			//Then function CleanRoots is called, in order to eliminate multiple approximations 
			//of the same root.  
			if ( abs(medial - x1) < tinylength ) {
				myroots.insert( pair<long double, long long unsigned int>(medial, numberofit) );
				numberofit++;
				CleanRoots(myroots, maccuracy, medial);	    	
			}//end of the if(abs(medial-x1)<tinylength)
		  	//__________________________________________________________________________________________

			//__________________________________________________________________________________________
			//Next, the signs of both subintervals (x1,medial) and (medial,x2) are calculated if their length
			//is greater that "tinylength".
			else if ( abs(medial - x1) >= tinylength ) {
				//__________________________________________________________________________________________
				//Check if "medial" is a root. If f(medial+e) and
				//f(medial-e) have opossite signs then store the value of the variable
				//"medial", in the map "myroots".
				//Call the function CleanRoots, to eliminate approximations for the same root.
				if ( Bolzano(medial - maccuracy, medial + maccuracy, function_) < 0 ) {
					myroots.insert( pair<long double, long long unsigned int>(medial, numberofit) );
					numberofit++;
					CleanRoots(myroots,  maccuracy, medial);	    	
				}//end of if(Bolzano(medial-e,medial+e)<0)
				//__________________________________________________________________________________________

				//__________________________________________________________________________________________
				//Calculate the signs of the subintervals (x1,medial)
				//and (medial, x2) of the current odd interval.
				sign = Bolzano(x1, medial, function_);//"sign" holds the sign of subinterval (x1,medial).
				sign2 = Bolzano(x2, medial, function_);//"sign2" holds the sign of the subinterval (medial,x2)

				if ( sign == -1 ) {//Check if the subinterval (x1,medial) has opposite 
					//endpoints and update the corresponding sets ("evenandodd", "odd"
					//and "lengthodd").
					evenandodd[cut][1] += 1;
					odd[position][0] = 2 * interpolations;
					odd[position][1] = cut;
					lengthodd.insert( pair<long long unsigned int, long long unsigned int>(cut, position) );
				}//end of if(sign==-1)

				if ( sign == 1 ) {//Check if the subinterval (x1,medial) has same signed endpoints 
					//and update the corresponding sets ("evenandodd", "even" and "lengtheven").
					evenandodd[cut][0] += 1;
					if ( ar == 1 ) {//If ar is equal to 1, then the vector "even" is not empty. 
						//So the subinterval (x1,medial) has to be stored at the end of the
						//vector "even". 
						h1[0] = 2 * interpolations;
						h1[1] = cut;
						even.push_back(h1);
						lengtheven.insert( pair<long long unsigned int, long long unsigned int>(cut, even.size()-1) );
					}
					else if ( ar == 0 ) {//If ar is equal to 0, the vector "even" is empty. 
						//The subinterval(x1,medial) is the first even subinterval that the program 
						//has discovered. So the subinterval (x1,medial) has to be store at the
						//first row of the vector "even".
						ar = 1;
						even[0][0] = 2 * interpolations;
						even[0][1] = cut;
						lengtheven.insert( pair<long long unsigned int, long long unsigned int>(cut, 0) );
					}
				}//end of if(sign==1)

				if ( sign2 == -1 ) {//Check if the subinterval (medial,x2) has opposite signs 
					//at its endpoints and update the corresponding sets ("evenandodd",
					//"odd" and "lengthodd").
					evenandodd[cut][1] += 1;
					odd[position][0] = 2 * interpolations + 1;
					odd[position][1] = cut;
					lengthodd.insert( pair<long long unsigned int, long long unsigned int>(cut, position) );
				}

				if ( sign2 == 1 ) {//Check if the subinterval (medial,x2) has the same signs at its 
					//endpoints and update the corresponding sets ("evenandodd", "even" and 
					//"lengtheven").
					evenandodd[cut][0] += 1;
					if ( ar == 1 ) {//If ar is equal to 1, the vector "even" is not empty. So the subinterval
						// (medial,x2) has to be stored at the end of the vector named
						//"even". 
						h1[0] = 2 * interpolations + 1;
						h1[1] = cut;
						even.push_back(h1);
						lengtheven.insert( pair<long long unsigned int, long long unsigned int>(cut, even.size() - 1) );
					}
					else if ( ar == 0 ) {//If ar is equal to 0, the vector "even" is empty.
						//The subinterval(medial,x2) is the first even subinterval that 
						//the program produces.
						//So the subinterval (medial,x2) has to be stored at the first row 
						//of the vector "even".
						ar = 1;
						even[0][0] = 2 * interpolations + 1;
						even[0][1] = cut;
						lengtheven.insert( pair<long long unsigned int, long long unsigned int>(cut, 0) );
					}
				}//end of if (sign2==1)

				//End of the section, in which the program calculates the signs of the subintervals 
				//(x1,medial) and (medial, x2) of the current odd interval.  
				//__________________________________________________________________________________________

			}//end of if(abs(medial-x1)>=tinylength)

			//End of the section, in which the program calculates the signs of both subintervals (x1,medial) and 
			//(medial,x2), when the current subinterval (x1,x2) is odd and the two halfs (x1,medial) and 
			//(medial,x2) have length greater than the smallest length ("tinylength"), which they can take.
			//__________________________________________________________________________________________
		}//End of statistic_stopping_criterion.
	}//End of the "do while" of odd subintervals

	number_even = lengtheven.size();
	number_odd = lengthodd.size();

}//End of function CutOdds		

//This function divides the even uncut subintervals.
void Manbis::CutEven(){

	if ( lengtheven.size() != 0 ) {
		itmultimap = lengtheven.begin();//"itmultimap" holds the location
		//of the subinterval that will be examined at the current iteration.

		cut = itmultimap -> first + 1;//"cut" contains the depth of the subdivision tree where the two halfs of the current
		//subinterval lie. 

		position = itmultimap -> second;//"position" holds the position in vector "odd" 
		//of the current even subinterval.

		lengtheven.erase(itmultimap);//Erase the current subinterval from "lengtheven"
		//to avoid examinating it for a second time.

		interpolations = even[position][0];//"interpolations" is holding how many 
		//subintervals with length equal to  initialength/(2^(cut-1)) lie between the 
		//left point of the initial interval and the left point of the current 
		//odd subinterval.

		long double x1, x2;//These variables hold the left and
		//the right point of the current even subinterval.

		StatisticStopingCriterion();
		//If the statistics are not ready to be printed out or the stopping criterion of the 
		//statistics doesn't hold, then the program starts to examinate the largest even subinterval
		//at the following lines. 

		if ( statistic_stopping_criterion == 0 ) {

			//________________________________________________________________________________
			//_____________________________________________________________________________________
			//Calculation of the number power=2^(odd[position][1])
			//for the word size of the OS.
			if ( (even[position][1]) >= 64 ) {
				power = 1;
				helphtree = even[position][1];
				while ( helphtree >= 64 ) {
					power *= ( unit << 63 );
					helphtree -= 63;
				}
				power *= ( unit << helphtree );
			} else {
				power = ( unit << ( even[position][1]) );
			}
			//________________________________________________________________________________


			x1 = initialeft + even[position][0] * (initialength / power);//"x1" holds the left endpoint of
			//the current even subinterval.

			x2 = x1 + initialength / power;//"x2" holds the right endpoint of the current even subinterval.

			medial = (x1 + x2) / 2;//"medial" holds the middle point of the current even
			//subinterval (x1,x2).


			//__________________________________________________________________________________________
			//Check if "medial" is a root. If f(medial+e) and
				//f(medial-e) have opossite signs then store the value of the variable
				//"medial", in the map "myroots".
				//Call the function CleanRoots, to eliminate approximations for the same root.
			if ( Bolzano(medial - maccuracy, medial + maccuracy, function_) < 0 ) {
				myroots.insert ( pair<long double, long long unsigned int>(medial, numberofit) );
				numberofit++;
				CleanRoots(myroots, maccuracy, medial);	    	
			}//End of if(Bolzano(medial-e,medial+e)<0)
  
			//__________________________________________________________________________________________

			//__________________________________________________________________________________________
			//At the followings lines, the program calculates the signs of both subintervals (x1,medial) and 
			//(medial,x2), when the current subinterval (x1,x2) is even and the two halfs (x1,medial) and 
			//(medial,x2) have length greater than the smallest length ("tinylength").
			//If the length of the two halfs is lower than the length of "tinylength", then the program will not
			//store the two halfs, as the current subinterval has been examinated thoroughly.

			if ( abs(medial - x1) >= tinylength ) {
				//Calculate the signs of the subintervals (x1,medial)
				//and (medial, x2) of the current even interval.
				sign = Bolzano(x1, medial, function_);//"sign" holds the sign of subinterval (x1,medial).
				sign2 = Bolzano(x2, medial, function_);//"sign2" holds the sign of the subinterval (medial,x2).

				if ( sign == -1 ) {//Check if the subinterval (x1,medial) has opposite 
					//endpoints and update the corresponding sets ("evenandodd", "odd"
					//and "lengthodd").
					evenandodd[cut][1] += 1;

					if ( per == 1 ) {//per=1 means that the vector "odd" is not empty. 
						//So the subinterval (x1,medial) has to be stored at the end of the vector
						//named "odd". 
						h1[0] = 2 * interpolations;
						h1[1] = cut;
						odd.push_back(h1);
						lengthodd.insert( pair<long long unsigned int, long long unsigned int>(cut,odd.size() - 1) );
					}//End of if(per==1)
					if ( per == 0 ) { //per=0 means that the vector "odd" is empty. The subinterval
						// (x1,medial) is the first odd subinterval that the program has discovered.
						//So the subinterval (x1,medial) has to be stored at the first row of the vector 
						//"odd".
						per = 1;
						odd[0][0] = 2 * interpolations;
						odd[0][1] = cut;
						lengthodd.insert( pair<long long unsigned int, long long unsigned int>(cut, 0) );
					}//End of if(per==0).

				}//End of if(sign==-1).

				if ( sign == 1 ) {//Check if the subinterval (x1,medial) has same signed endpoints 
					//and update the corresponding sets ("evenandodd", "even" and "lengtheven").
					evenandodd[cut][0] += 1;
					even[position][0] = 2 * interpolations;
					even[position][1] = cut;
					lengtheven.insert( pair<long long unsigned int, long long unsigned int>(cut, position) );
				}//End of if(sign==1).


				if ( sign2 == -1 ) {//Check if the subinterval (medial,x2) has opposite signs 
					//at its endpoints and update the corresponding sets ("evenandodd",
					//"odd" and "lengthodd").
					evenandodd[cut][1] += 1;
					h1[0] = 2 * interpolations + 1;
					h1[1] = cut;
					odd.push_back(h1);
					lengthodd.insert( pair<long long unsigned int, long long unsigned int>(cut, odd.size()-1) );
				}//End of if(sign2==-1).

				if ( sign2 == 1 ) {//Check if the subinterval (medial,x2) has the same signs at its 
					//endpoints and update the corresponding sets ("evenandodd", "even" and 
					//"lengtheven").
					evenandodd[cut][0] += 1;
					h1[0] = 2 * interpolations + 1;
					h1[1] = cut;
					even.push_back(h1);
					lengtheven.insert( pair<long long unsigned int, long long unsigned int>(cut, even.size() - 1) );
				}//End of if(sign2==1).


			}//end of if(abs(medial-x1)>=tinylength)
			//End of the section in which the program calculates the signs of both subintervals (x1,medial) and 
			//(medial,x2), when the current subinterval (x1,x2) is even and the two halfs (x1,medial) and 
			//(medial,x2) have length greater than the smallest length ("tinylength"), which they can take.
			//__________________________________________________________________________________________
		}//End of statistic_stopping_criterion
		number_even = lengtheven.size();
		number_odd = lengthodd.size();
	}//End of if(lengtheven.size()!=0)
}//End of function CutEven

//This function prints out the calculated roots.
void Manbis::PrintResults(){
	ofstream coutresult;
	if (flag==1) cout<< "Terminating. Results in file Results.txt"<<endl;
	coutresult.open("Results.txt");//Create the file Results.txt, 
	//This file will contain the calculated roots at the end of the program,
	//if the user selected this option

	coutresult.precision(mdigits);//Set the number of decimal digits of the
	//output numbers. Default value is 8.

	//_______________________________________________________________________________________________
	//Print in file Results.txt the values of some parameters 
	coutresult << "Requested percentage is " << mfraction * 100 << "%" << endl;
	coutresult << "Requested accuracy is " << maccuracy << endl;
	coutresult << "Number of the roots that the program has calculated is " << myroots.size() << endl;
	if ( ceil( avg ) >= myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) !=1 && isnan( ceil ( avg ) ) != 1 ){
		coutresult << "Estimated percentage of the calculated roots is " << ceil(( evenandodd[cut-1][1] / ( (nlower + nupper) / 2) ) * 100) << "%" << endl;
		coutresult << "Estimated total number of roots in the interval (" << xleft << "," << xright << ") is " << ceil( ( nlower + nupper) / 2 )<< "." << endl;
	}else if(ceil( avg ) < myroots.size() && ceil( ( nlower + nupper) / 2 ) < myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) != 1 && isnan( ceil( avg ) ) != 1 ){
		coutresult << "Estimated percentage of the calculated roots is " << 100 << "%" << endl;
		coutresult << "Estimated total number of roots in the interval (" << xleft << "," << xright << ") is " << myroots.size() << "." << endl;	
	}else if(ceil( avg ) < myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) != 1 && isnan( ceil( avg ) ) != 1 ){
		coutresult << "Estimated percentage of the calculated roots is " << ceil(( evenandodd[cut-1][1] / ( (nlower + nupper) / 2) ) * 100) << "%" << endl;
		coutresult << "Estimated total number of roots in the interval (" << xleft << "," << xright << ") is " << ceil( ( nlower + nupper) / 2 )<< "." << endl;
	}
	
	if (nextcut-1==1 ){
		cout << "\033[1;31m   Current      \033[0m"<<" | "<<"\033[1;32m   Requested   \033[0m"<<" | "<<"\033[1;33m   Discovered   \033[0m"<<" | "<<"\033[1;34m   Discovered   \033[0m"<<" | "<<"\033[1;36m   Avg time of computing    \033[0m"<<endl;
		cout << "\033[1;31m   estimation   \033[0m"<<" | "<<"\033[1;32m   #roots      \033[0m"<<" | "<<"\033[1;33m   #roots       \033[0m"<<" | "<<"\033[1;34m   fraction     \033[0m"<<" | "<<"\033[1;36m   100 roots    \033[0m"<<endl<<endl;
	}
				
	if ( mcurrentest == 1 ) {
		ct = clock() - ct;
		proots=myroots.size()-proots;
			if(proots>0 && !isnan(proots) && !isinf(proots)){
				if ( ceil( avg ) >= myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) !=1 && isnan( ceil ( avg ) ) != 1 ){
					cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" <<ceil( ( nlower + nupper) / 2 ) << '\t' << '\t' << "\033[1;32m" << ceil( avg ) << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t' << '\t' << "\033[1;34m" ;
					cout<< ceil(( evenandodd[cut-1][1] / ( (nlower + nupper) / 2) ) * 100) << "%" << '\t' << '\t' << "\033[1;36m" << (((float)ct) / CLOCKS_PER_SEC)/proots *100<<"sec"<< "\033[1;37m";
 					cout.flush();
	 			}else if(ceil( avg ) < myroots.size() && ceil( ( nlower + nupper) / 2 ) < myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) != 1 && isnan( ceil( avg ) ) != 1 ){
 					cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" << myroots.size() << '\t' << '\t' << "\033[1;32m" << myroots.size()  << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t'  << '\t' << "\033[1;34m" ;
 					cout<< 100 << "%" << '\t' << '\t' << "\033[1;36m" << (((float)ct) / CLOCKS_PER_SEC)/proots *100 <<"sec"<< "\033[1;37m";
 					cout.flush();
 				}else if(ceil( avg ) < myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) != 1 && isnan( ceil( avg ) ) != 1 ){
 				 	cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" << ceil( ( nlower + nupper) / 2 ) << '\t' << '\t' << "\033[1;32m" << myroots.size()  << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t'  << '\t' << "\033[1;34m" ;
 			 		cout << ceil(( evenandodd[cut-1][1] / ( (nlower + nupper) / 2) ) * 100) << "%" << '\t' << '\t' << "\033[1;36m" << (((float)ct) / CLOCKS_PER_SEC)/proots *100<<"sec"<< "\033[1;37m";
 					cout.flush();
	 			}	
 			}else {
 				if ( ceil( avg ) >= myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) !=1 && isnan( ceil ( avg ) ) != 1 ){
					cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" << ceil( ( nlower + nupper) / 2 ) << '\t' << '\t' << "\033[1;32m" << ceil( avg ) << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t' << '\t' << "\033[1;34m" ;
					cout<< ceil(( evenandodd[cut-1][1] / ( (nlower + nupper) / 2) ) * 100) << "%" << '\t' << '\t' << "\033[1;37m" << string(30, '-') ;
 					cout.flush();
	 			}else if(ceil( avg ) < myroots.size() && ceil( ( nlower + nupper) / 2 ) < myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) != 1 && isnan( ceil( avg ) ) != 1 ){
 					cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" << myroots.size() << '\t' << '\t' << "\033[1;32m" << myroots.size()  << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t'  << '\t' << "\033[1;34m" ;
 					cout<< 100 << "%" << '\t' << '\t' << "\033[1;37m" << string(30, '-');
 					cout.flush();
 				}else if(ceil( avg ) < myroots.size() && isnan( ceil( ( nlower + nupper) / 2 ) ) != 1 && isnan( ceil( avg ) ) != 1 ){
 				 	cout  <<'\r'<< string(100, '-') << '\r' << "---\033[1;31m" << ceil( ( nlower + nupper) / 2 ) << '\t' << '\t' << "\033[1;32m" << myroots.size()  << '\t' << '\t' << "\033[1;33m" << myroots.size() << '\t'  << '\t' << "\033[1;34m" ;
 			 		cout << ceil(( evenandodd[cut-1][1] / ( (nlower + nupper) / 2) ) * 100) << "%" << '\t' << '\t' << "\033[1;37m" << string(30, '-');
 					cout.flush();
 				}
 			}
 	}		

	if ( mtotaltime == 1 ) {
		coutresult << "Run time of Manbis: " << t << " clock ticks (" << ((float)t) / CLOCKS_PER_SEC << "seconds)." << endl;
		cout << endl << "Run time of Manbis: "<< t << " clock ticks (" << ((float)t) / CLOCKS_PER_SEC << "seconds)." << endl;
	}
	if ( mnor == 1 ) {
		coutresult << "The roots are: " << endl;

		for ( itmyroots = myroots.begin(); itmyroots != myroots.end(); ++itmyroots ) {
			coutresult << itmyroots -> first  << '\n';
		}
	}

	coutresult.close();//Close file Results.txt.

}//End of function PrintResults

//Function GetExistOdd returns the value of per. per=1 means that an
//odd subinterval has been created else per=0.
int Manbis::GetExistOdd(){
	return per;
}

//Returns how many odd 
//uncut subintervals exist. 
long long unsigned int Manbis::GetNumberOfOdd(){
	return number_odd;
}

//Returns how many even uncut
//subintervals exist. 
long long unsigned int Manbis::GetNumberOfEven(){
	return number_even;
}

//Print the options to the screen
void Manbis::PrintChoices(){
	cout << endl << "Accuracy: " << maccuracy << endl;
	cout << "Initial interval: " << "(" << xleft << "," << xright << ")" << endl;
	cout << "Percentage: " << mfraction * 100 << "\%" << endl;
	cout << "Number of digits of output roots:" << mdigits << endl << endl;
}

//This function is the that Manbis calculates the roots
//using all the above mentioned functions 
void Manbis::FindRoots(){
	PrintChoices();//Prints the default values of some parameters.

  	t = clock();//Holds the total time that Manbis needs, in order to give results.
  	ct = clock();//Holds the current time that Manbis needs, in order to copmute a root at the every level.


	int exist_odd=GetExistOdd();//This variable checks if at least one uncut odd subinterval exists	

	long long unsigned int number_of_odd = GetNumberOfOdd();//Holds the amount of the 
	//odd uncut subintervals that exist at the current iteration.
	long long unsigned int number_of_even = GetNumberOfEven();//Holds the amount of the 
	//even uncut subintervals that exist at the current iteration.

	PrimordialStep();//Call the function PrimordialStep, which
	//sets some important initializations of the program.

	signal(SIGINT, UserInterruption);   // Register signals 
//---------------------------------START OF MAIN STEP-------------------------------------
//________________________________________________________________________________________

	exist_odd=GetExistOdd();//Update the variable, that checks 
	//if at least one uncut odd subinterval exists	
	number_of_odd=GetNumberOfOdd();//Update the variable, that holds
	//the amount of odd uncut subintervals that exist at the current iteration.
	number_of_even=GetNumberOfEven();//Update the variable, that holds
	//the amount of even uncut subintervals that exist at the current iteration.
	do{//Loop for even subintervals

		//___________________________Start the "do while" of the odd subintervals____________________________
		// Odd subintervals get higher priority than the even subintervals
		while ( number_of_odd != 0 && exist_odd == 1 ) {//Starting point of the while for odd subintervals.
        //The flag "exist_odd", signifies whether an odd subinterval exists. 
        //The  amount_of_odd holds the size of the multimap lengthodd,
        //which exists in the header file Manbis.h.

        	CutOdds();//Call the function CutOdds. 
			//If the program reaches this point, then 
			//at least one odd uncut subinterval needs to be examined
			//and it will be divided once.

		exist_odd = GetExistOdd();//Update the variable, that checks 
		//if at least one uncut odd subinterval exists	
		number_of_odd = GetNumberOfOdd();//Update the variable, that holds
		//the amount of odd uncut subintervals that exist at the current iteration.
		number_of_even=GetNumberOfEven();//Update the variable, that holds
		//the amount of even uncut subintervals that exist at the current iteration.

		}//End of while(amount_of_odd!=0 && exist_odd==1). All odd subintervals
		//have been divided.

		CutEven();//Call the function CutEven. 
			//If at least one even uncut subinterval exists
			//it will be divided once.
	
		exist_odd = GetExistOdd();//Update the variable, that checks 
		//if at least one uncut odd subinterval exists	
		number_of_odd = GetNumberOfOdd();//Update the variable, that holds
		//the amount of odd uncut subintervals that exist at the current iteration.
		number_of_even = GetNumberOfEven();//Update the variable, that holds
		//the amount of even uncut subintervals that exist at the current iteration.

	}while ( number_of_even != 0 && flag==0);//End of "while" even subintervals

 	t = clock() - t;
 	

 	PrintResults();//Call the function PrintResults,
 	//which will print out the calculated roots in the txt file Results.txt.

	//---------------------------------END OF MAIN STEP----------------------------------------
	if(flag){ // my action when signal set it 1
        PrintResults();   

    } 
}

//Function Bolzano returns -1 if f(x1) and f(x2) have opposite signs
//else returns the value 1. The function takes two Real values as arguments
//and one pointer function, which holds the memory address of the user's 
//formula in the Functions.cpp
int Bolzano(long double x1,long double x2, long double (*function_)(long double ) ){

	long double f1, f2;//Set two Real variables.

	f1 = (*function_)(x1);//"f1" contains the functional value of the point named "x1".
	f2 = (*function_)(x2);//"f2" contains the functional value of the point named "x2".

	if ( f1 * f2 < 0 ) 
		return -1;//if f(x1) and f(x2) have opposites signs the function returns -1.
	else
		return 1;//if f(x1) and f(x2) have same signs, the function returns 1.

}//End Bolzano.

void UserInterruption(int sig){ // can be called asynchronously
  flag = 1; // set flag
}


