// Fast computation of the non asymptotic Cochran's Q statistic for heterogeneity detection
// Fahmy and Belletoile
//

#include <math.h>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <ctime>

#define  Eps12   1.0e-12

typedef uint16_t					Tshort;
typedef std::vector<size_t>			vecTsize;
typedef std::vector<Tshort>			vecTshort;
typedef std::vector<vecTshort>		matTshort;
typedef std::pair<vecTshort,double>	pairCij_W;
typedef std::vector<pairCij_W>		vecPairCij_W;
typedef std::vector<vecPairCij_W>	matPairCij_W;

/*****************************************************************************************************************************************************************
*** \brief   Compute the binomial coefficient "choose k among n"
*** \param   [in] n : double value
*** \param   [in] k : double value
*** \return  returns the natural value of the binomial coefficient
*** \author  Arnaud Belletoile
****************************************************************************************************************************************************************/
size_t Cnk(size_t n, size_t k) 
{
    size_t p, hCnk = 0;

	if (n < k) return 0;
    if (k == n || k == 0)
	{
        return 1; 
	}
    p = size_t(n - k);
    if (k < p) 
	{
        p = size_t(k);
		k = n - p;
	}
	hCnk = k + 1;
    for (size_t i = 2; i<= p; i++)
	{
        hCnk = (hCnk * (k + i)) / i;
	}

	return hCnk;
}

/*****************************************************************************************************************************************************************
*** \brief   Compute the log of the binomial coefficient 
*** \param   [in] n : double value
*** \param   [in] k : double value
*** \return  returns the log of the binomial coefficient
*** \author  Arnaud Belletoile
****************************************************************************************************************************************************************/
double LogCnk(double n , double k )
{
    int i=0, p=0;
	double sum=0;
    if (k == n || k == 0) 
	{
        return 0;
	}
    p = int(n) - int(k);
    if (k < p)
	{
        p = int(k);
        k = n - p;
	}
    sum = log(k + 1);
    for (i = 2; i <= p; i++)
	{
        sum = sum + log(k + i) - log(double(i));
	}
    return sum;
}

/*****************************************************************************************************************************************************************
*** \brief   Check that time limit hasn't been exceeded
*** \param   [in] start : initialized at the begining of the program
*** \param   [in] timeLimit : parameter set by user if timeLimit == 0 then no limit
*** \return  returns true if the time limit was exceeded, false otherwise
*** \author  Arnaud Belletoile
****************************************************************************************************************************************************************/
bool exitOnTime(std::clock_t start, double &lastStop, const double timeLimit) 
{
	double q;
	double t = (std::clock() - start) / (double) CLOCKS_PER_SEC; 

	if (lastStop != 0)
		q = t / lastStop;
	else
		q = 1;

	if ((timeLimit > Eps12) && (q * t > timeLimit))
		return true;

	lastStop = t;

	return false;
}

/*****************************************************************************************************************************************************************
*** \brief   Recursive function to examine all possible distributions of S "balls" in maxP "baskets"
***			 The function is called a first time with tab filled with 0
***			 If a solution is found (all S balls could be distributed in the basket(s) without exceeding the upper limit given by tabMax), tab is updated and 0 is returned
***			 Each successive call needs the last tab as input, if another solution is found, tab is updated and 0 is returned.
***			 If all possible combinations (possibly none) has been explored, a value > is returned
*** \param   [in-out] tab : array of size maxP containing the number of balls in each basket 
*** \param   [in] tabMax : array of size maxP containing the upper limit on the number of balls for each basket
*** \param   [in] P : indice of current basket [0; maxP]
*** \param   [in] maxP : number of baskets
*** \param   [in] S : number of success
*** \return  returns the number of balls that haven't been distributed into one basket. 0 means that a new solution was found, >0 that the algorithm couldn't find any solution
*** \author  Arnaud Belletoile
****************************************************************************************************************************************************************/
 int getNextTuple(vecTshort &tab, vecTshort &tabMax, size_t P, const size_t maxP, Tshort S)
{

	if (tab[P] == S)
	{
		tab[P] = 0;
		return(S);
	}

	if (P < maxP - 1)
	{
		while(getNextTuple(tab, tabMax, P+1, maxP, S-tab[P]) > Tshort(0))
		{
			tab[P]++;
			if (tab[P] > tabMax[P])
			{
				tab[P] = 0;
				return(S);
			}
		}
	}
	else if (tab[P] + S <= tabMax[P])
	{
		tab[P] = S;
		return(0);
	}
	else
		return(S);

	return(0);
}

/*****************************************************************************************************************************************************************
*** \brief   Computes the Cochran's Q statistic as
*** \param   [in] Qden : allready computed denominator in the expression of Q 
*** \param   [in] C : vector of size P of sums over rows of the original table
*** \param   [in] P : number of treatment/columns
*** \return  the Q statistic
*** \author  Arnaud Belletoile
****************************************************************************************************************************************************************/
double getCochranQ(double Qden, vecTshort &C, const size_t P)
{
	double Qnum;
	double sumSimple = 0;
	double sumSquare = 0;

	if (Qden == 0)
		return -1;

	for (size_t i = 0; i < P; i++)
	{
		sumSquare += C[i] * C[i];
		sumSimple += C[i];
	}

	Qnum = (P - 1) * ( (P * sumSquare) - sumSimple * sumSimple);
	
	return(Qnum/Qden);
}


/*****************************************************************************************************************************************************************
*** \brief   Merges A and B in A
*** \param   [in-out] A vector of pairs 
*** \param   [in] B vector of pairs 
*** \return  0 on success, exceptions are thrown otherwise
*** \author  Arnaud Belletoile
****************************************************************************************************************************************************************/
int merge(vecPairCij_W &A, vecPairCij_W &B) 
{
	size_t i = 0, j = 0;
	double temp;
	vecPairCij_W C;

	C.reserve(A.size() + B.size());

	if (A.front().first > B.back().first)
	{
		C = move(B);
		C.insert(C.end(), A.begin(), A.end());
		A = move(C);
		return(1);
	}

	if (B.front().first > A.back().first)
	{
		C = move(A);
		C.insert(C.end(), B.begin(), B.end());
		A = move(C);
		return(1);
	}

	while (i < A.size() && j < B.size()) 
	{
		if (A[i].first < B[j].first) 
		{
			C.push_back(A[i]);
            i++;
		}
		else 
		{
			C.push_back(B[j]);
			if(A[i].first == B[j].first)
			{
				temp = exp(A[i].second);
				temp += exp(B[j].second);
				C.back().second = log(temp);
				i++;
			}
			j++;
		}
	}

	if (i < A.size()) 
	{
       for (size_t p = i; p < A.size(); p++) 
	   {
			C.push_back(A[p]);
	   }
    }
	else
	{
        for (size_t p = j; p < B.size(); p++) 
	   {
			C.push_back(B[p]);
	   }
    }

	A = move(C);
	return(0);
}

/*****************************************************************************************************************************************************************
*** \brief   Sorts A according to the first element of the pair, removes duplicates and merges weights (second element of the pair)
*** \param   [in-out] A vector of pairs 
*** \return  0 on success, exceptions are thrown otherwise
*** \author  Arnaud Belletoile
****************************************************************************************************************************************************************/
int	reduce(vecPairCij_W &A)
{
	vecPairCij_W B;

	B.reserve(A.size());

	double temp;

	if (A.size() > 0)
	{
		sort(A.begin(), A.end());
		B.push_back(A[0]);
		for (size_t i = 1; i < A.size(); i++)
		{
			if (A[i].first == B.back().first)
			{
				temp = exp(A[i].second);
				temp += exp(B.back().second);
				B.back().second = log(temp);
			}
			else
			{
				B.push_back(A[i]);
			}
		}

		A = move(B);
	}

	return(0);

}

/*****************************************************************************************************************************************************************
*** \brief   Transforms the solution vector [5 4 4 2 2 2 1] into the number of baskets with sizes [1 2 3 1]
*** \param   [in] Cij_W : solution vector of length P containing the distributed sucesses (ex : [5 4 4 2 2 2 1])
*** \param   [in] P : number of treatment/columns, length of the solution vector
*** \param   [out] K: vector of length the number of "baskets", elements are the size of each basket (ex : [1 2 3 1])
*** \return  void
*** \author  Arnaud Belletoile
****************************************************************************************************************************************************************/
void makeBins(pairCij_W &Cij_W, size_t P,vecTshort &K)
{
	K.assign(1, 1);

	for (Tshort j = 1; j < P; j++)
	{
		if (Cij_W.first[j] != Cij_W.first[j-1])
			K.push_back(1);
		else
			K.back() += 1;

		if (Cij_W.first[j] == 0)
		{
			K.back() += (Tshort) P - j - 1;
			break;
		}
	}
}

/*****************************************************************************************************************************************************************
*** \brief   Performs computations and returns the exact p-value of a Cochran's Q test
*** \param   [in] table [N x P] N-vector of P-vector of results coded with unsigned Tshort int (16 bits): 1 (success) or 0 (failure)
*** \param   [in] timelimit : time limit in seconds (double) > 0 -- optional, if == 0 no limit
*** \return  exact p-values [0; 1] if success, <0 if the program was aborted due to the time limit
*** \error	 Exceptions from the stl should be tried/catched around this function
*** \author  Arnaud Belletoile
****************************************************************************************************************************************************************/
double getPval(matTshort &table, const double timeLimit = 0)
{
	/************************************************************************************************************************************************************/
	/* Check parameters
	************************************************************************************************************************************************************/
	size_t N = table.size();
	if (N == 0) 
		return -3;

	size_t P = table[0].size();
	if (P == 0) 
		return -3.;

	for (size_t i = 1; i < N; ++i)
		if (table[i].size() != P) // not a matrix!
			return -3.;

	if (timeLimit < 0)
		return -3.;

	/************************************************************************************************************************************************************/
	/* Time measurement init + some declarations
	************************************************************************************************************************************************************/
	std::clock_t start = std::clock();
	double lastStop = 0;

	int idx = 0;
	vecTshort R(N, 0), C(P,0);
	vecTshort K, X, classes;
	std::vector<std::vector<double>> logCnk;
	
	matPairCij_W Cij_W; 
	vecPairCij_W temp_Cij_W(1, make_pair(C, 0));

	double Qval = 0;
	double sumSquare = 0;
	double Qden = 0;

	size_t szC = 0, szK = 0;


	/************************************************************************************************************************************************************/
	/* Compute all possible log(Cnk) values once and for all
	************************************************************************************************************************************************************/
	K.reserve(P); 	
	logCnk.reserve(P + 1);
	std::vector<double> tempLogCnk;
	for (size_t i = 0; i < P + 1; i++ )
	{
		for (size_t j = 0; j < i + 1; j++)
			tempLogCnk.push_back(LogCnk(double(i), double(j)));
		logCnk.push_back(tempLogCnk);
		tempLogCnk.clear();
	}

	/************************************************************************************************************************************************************/
	/* Compute sums per rows, sums per columns, overall sums and some usefull values to be used to compute the Q statistic
	************************************************************************************************************************************************************/
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < P; j++)
		{
			R[i] += table[i][j];
			C[j] += table[i][j];
		}
		sumSquare += R[i] * R[i];
		Qden += R[i];
	}

	std::sort(R.begin(), R.end());
	std::reverse(R.begin(), R.end());

	Qden = P * Qden - sumSquare;
	Qval = getCochranQ(Qden, C, P);
	if (Qval < 0 )
		return -4.;

	/************************************************************************************************************************************************************/
	/* Loop on each observation -- This loop can be parallelized provided some minor modifications to prevent race conditions in the merging operations
	************************************************************************************************************************************************************/
	for (size_t i = 0; i < N; i++)
	{
		if (0 < R[i] && R[i] < P) // only value with at least one success and one failure contribute to Q
		{
			Cij_W.resize(temp_Cij_W.size()); 

			for (size_t k = 0; k < temp_Cij_W.size(); k++) 
			{				
				/************************************************************************************************************************************************************/
				/* First transform the solution vector in a number of "baskets" with a size -- ex: Cij = [5 4 4 2 2 2 1] --> K = [1 2 3 1]
				************************************************************************************************************************************************************/
				makeBins(temp_Cij_W[k], P, K);
				
				X.assign(K.size(), 0);

				Cij_W[k].reserve(Cnk((R[i] + K.size() - 1), (K.size() - 1)));
				
				/************************************************************************************************************************************************************/
				/* Explore every possible arrangement of the current number of successes
				************************************************************************************************************************************************************/
				while (getNextTuple(X, K, 0, K.size(), R[i]) == 0)
				{
					idx = 0;
					Cij_W[k].push_back(temp_Cij_W[k]);

					for (size_t l = 0; l < K.size(); l++)
					{
						if (l > 0)
							idx += K[l - 1];

						Cij_W[k].back().second += logCnk[K[l]][X[l]];
						for (Tshort iiClass = 0; iiClass < X[l]; iiClass++)
							Cij_W[k].back().first[idx + iiClass] += 1;
					}

					/************************************************************************************************************************************************************/
					/* Make sure that solution vectors Cij are properly sorted for comparisons when removing duplicates and for the number of "baskets" during the next iteration
					************************************************************************************************************************************************************/
					sort(Cij_W[k].back().first.begin(), Cij_W[k].back().first.end());
					reverse(Cij_W[k].back().first.begin(), Cij_W[k].back().first.end());

				}

				/************************************************************************************************************************************************************/
				/* Remove duplicates
				************************************************************************************************************************************************************/
				reduce(Cij_W[k]);
			}
			
			/************************************************************************************************************************************************************/
			/* Merge branches 
			************************************************************************************************************************************************************/
			szC = Cij_W.size();
			while (szC > 1)
			{
				szK = size_t(floor(0.5*szC));
				/************************************************************************************************************************************************************/
				/* this loop can be parallelized
				************************************************************************************************************************************************************/
				for (size_t k = 0; k < szK; k++)
					merge(Cij_W[k], Cij_W[szC-k-1]);
				szC -= szK;
				Cij_W.resize(szC);
			}

			temp_Cij_W = move(Cij_W[0]);

			if (exitOnTime(start, lastStop, timeLimit)) 
				return -5.;
		}
	}

	/************************************************************************************************************************************************************/
	/* Compute and the return the p-value
	************************************************************************************************************************************************************/
	double Ntot = 0, Nabove = 0;
	for (size_t k = 0; k < temp_Cij_W.size(); k++)
	{
		double tempW = exp(temp_Cij_W[k].second);
		if (getCochranQ(Qden, temp_Cij_W[k].first, P) - Qval > -Eps12)
			Nabove += tempW;
		Ntot += tempW;
	}
	return(Nabove/Ntot);
}

/*****************************************************************************************************************************************************************
*** \brief   Call the getPval sub-function where all computations are performed and handle exceptions
*** \param   [in] table [N x P] N-vector of P-vector of results coded in unsigned Tshort int (16 bits): 1 (success) or 0 (failure)
*** \param   [out] pVal : exact p-value(double)
*** \param   [in] timelimit : time limit in seconds (double) > 0 -- optional, if == 0 no limit
*** \return  0 on success, -1 on bad allocation(memory issue), -2 on out of range (wrong access), -3 bad parameter, -4 unknown error, -5 if time limit is exceeded, -6 unknown exception
*** \author  Arnaud Belletoile
****************************************************************************************************************************************************************/
int xlstatCoqExact(matTshort &table, double &pVal, const double timeLimit = 0)
{
	try
	{
		pVal = getPval(table, timeLimit);

		if (pVal == -3.)
		{
			std::cerr << "Bad parameters: program terminated" << '\n'; 
			return -3;
		}
		else if (pVal == -4.)
		{
			std::cerr << "Unknown error: program terminated" << '\n';
			return -4;
		}
		else if (pVal == -5.)
		{
			std::cerr << "time limit exceeded: program terminated" << '\n';
			return -5;
		}
	}
	catch (const std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << '\n';
		return -1;
	}
	catch (const std::out_of_range& oor) 
	{
		std::cerr << "out_of_range error: " << oor.what() << '\n';
		return -2;
	}
	catch (const std::bad_exception be) 
	{ 
		std::cerr << "caught unknown bad_exception\n"; 
		return -6;
	}
	catch (...)
	{ 
		std::cerr << "caught unknown exception\n"; 
		return -6;
	}


	return(0);
}




