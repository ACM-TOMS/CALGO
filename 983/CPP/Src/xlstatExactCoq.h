// Fast computation of the non asymptotic Cochran's Q statistic for heterogeneity detection
// Fahmy and Belletoile
//

#ifndef XLSTATEXACTCOQ_H
#define XLSTATEXACTCOQ_H

#include <math.h>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <ctime>

#define	Eps12   1.0e-12

typedef uint16_t					Tshort;
typedef std::vector<size_t>			vecTsize;
typedef std::vector<Tshort>			vecTshort;
typedef std::vector<vecTshort>		matTshort;
typedef std::pair<vecTshort,double>	pairCij_W;
typedef std::vector<pairCij_W>		vecPairCij_W;
typedef std::vector<vecPairCij_W>	matPairCij_W;

size_t	Cnk(size_t n, size_t k);
double	LogCnk(double n , double k);
bool	exitOnTime(std::clock_t start, double &lastStop, const double timeLimit);
int		getNextTuple(vecTshort &tab, vecTshort &tabMax, size_t P, const size_t maxP, Tshort S);
double	getCochranQ(double Qden, vecTshort &C, const size_t P);
int		merge(vecPairCij_W &A, vecPairCij_W &B);
int		reduce(vecPairCij_W &A);
void 	makeBins(pairCij_W &Cij_W, size_t P,vecTshort &K);
double	getPval(matTshort &table, const double timeLimit = 0);
int 	xlstatCoqExact(matTshort &table, double &pVal, const double timeLimit = 0);

#endif

