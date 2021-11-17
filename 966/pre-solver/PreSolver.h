/*****************************************************************************
 * This code is part of Art Gallery Solver (AGSol) Package, which aims the 
 * resolution of the Art Gallery Problem With Point Guards.
 *
 * This software version (1.0.2) has been tested under and is compatible with 
 * CGAL 3.9 and GLPK 4.52.
 *
 * Authors:
 *  Davi Colli Tozoni - davi.tozoni@gmail.com
 *  Marcelo Castilho Couto - coutomarcelo@gmail.com
 *
 * AGSol Concept and Design: 
 *  Davi Colli Tozoni, Marcelo Castilho Couto, Pedro Jussieu de Rezende & Cid 
 * Carvalho de Souza.
 * 
 * Other information can be found at:
 *  http://www.ic.unicamp.br/~cid/Problem-instances/Art-Gallery/index.html
 *
 * --
 * 
 *  This program is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation, either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.
 *      
 ****************************************************************************/


#ifndef PRESOLVER_H
#define PRESOLVER_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <tr1/unordered_set>

#ifdef XPRESS
#include "SolverPLIXpress.h"
#endif

#ifdef GLPK
#include "SolverPLIGlpk.h"
#endif

#include "Lagrangian.h"

/* Rotinas para medicao de  tempo. Elas foram implementadas pelo Ionut Aaron (Brown,GSIA) */
#ifndef __SYSTEM_TIME_H
#define __SYSTEM_TIME_H

/* Ionut Aron, ia@cs.brown.edu (1999-2001)
 * created: August, 1999
 * last updated: December 20, 2001
 */

#include <sys/resource.h>
#include <sys/types.h>
#include <time.h>

extern int getrusage();

// this macro produces a time equal to the one produced by clock(),
// but does not suffer from the wraparound problem of clock()

#define CPUTIME(ruse) (\
  getrusage(RUSAGE_SELF,&ruse),\
  ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + \
  1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec)\
)

/* usage of the CPUTIME timing routine:
   ===================================

   struct rusage ruse;
   double t0 = CPUTIME(ruse);
   ...
   double t1 = CPUTIME(ruse);
   printf("CPU time: %f secs.\n", t1-t0);
*/

#endif

#define SMALL_NUMBER_GAP 0.00000001

using namespace std;

class PreSolver {
    private:
        vector< vector<bool> > _matrix;
        vector< boost::dynamic_bitset<> > _iBitsMatrix;
        vector< vector<int> > _groups;
        vector<bool> _cAlive;
        vector<bool> _rAlive;

        // maps variables from the presolved matrix to the original one
        vector<int> _redMap; 

        struct rusage _ruse;
        ofstream _log;

        int _nCand;
        int _nWitn;
        int _nAliveCols;
        int _nAliveRows;
    
        double _extLB;

        // Other variables:
        vector<int> _heurSolution; // indices from reduced matrix
        vector<int> _optSolution; // indices from reduced matrix

        bool _foundOptimal;

    public:
        /**
         * Constructor. Receives the ILP Matrix, the guard candidates groups
         * (which are organized by the Light AVPs where the candidates belong)
         * and the lower bound previously obtained by the solution. Using all
         * this information, we setup the structures used when solving the ILP.
         */
        PreSolver(const vector< vector<bool> > &matrix, vector< vector<int> > &groups, double extLB = 0.0); 

        /**
         * Returns if the solution found is optimal.
         */
        bool isOptimal() { return _foundOptimal; };

        /**
         * Eliminates redundant columns.
         */
        void bitsColsReduction();
        
        /**
         * Eliminates redundant rows.
         */
        void bitsRowsReduction();
        
        /**
         * Checks a group of guard candidates to find redundant columns. 
         */
        void bitsCheckGroup(int index);

#ifdef XPRESS
        /**
         * Prepares the necessary data in order to call the SCP solver library
         * using XPRESS.
         */
        void solveWithXpress(bool preXpress, double solValue, vector<int> sol);
#endif

#ifdef GLPK
        /**
         * Prepares the necessary data in order to call the SCP solver library
         * using GLPK.
         */
        void solveWithGlpk(double solValue, vector<int> sol);
#endif

        /**
         * Solves an SCP using the mode selected by the user.
         */
        void solve(int mode);

        /** 
         * Uses Lagrangian Heuristic to find a good viable solution for the
         * current SCP.
         */
        vector<int> tryHeuristic();

        /**
         * Constructs the optimal solution found using the initial indices.
         */
        vector<int> getOptimalSolution();
};

#endif
