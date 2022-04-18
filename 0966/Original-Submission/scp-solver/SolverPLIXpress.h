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
 * ArtGalleryHope Concept and Design: 
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


#ifndef SOLVER_PLI_XPRESS_H
#define SOLVER_PLI_XPRESS_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include "xprs.h"

using namespace std;

class SolverPLIXpress {
    private:
        static const int MAX_CPU_TIME = 3600;
        
        XPRSprob _prob;
        double* _xStar; // saves best integer solution found
        double  _zStar; // saves best value found
        int _nCols;
        int _problemType;

        double _extLB;
        bool _extLBstop;
        

        SolverPLIXpress(); // avoid this constructor

        void errorMsg(const string, int, int);

    public:
        static const int SET_COVER = 1;
        static const int USE_XPRESS = 1;
        static const int NOT_USE_XPRESS = 0;

        SolverPLIXpress(int, int, int, int, int*, int*, int*, double*, double extLB = 0.0);
        
        ~SolverPLIXpress(){
            /* Frees memory */ 
            int xpressRet = XPRSdestroyprob(_prob);
            if (xpressRet) errorMsg("Main [XPRSdestroyprob]", __LINE__, xpressRet);
            
            xpressRet = XPRSfree();
            if (xpressRet) errorMsg("Main [XPRSfree]", __LINE__, xpressRet);

            freeXStar();
        }

        /** 
         * Loads SCP matrix in Xpress.
         * @param nCols Number of variables
         * @param nRows Number of constraints
         * @param mRwInd First entrance for each a_{i,j} set (not 0); represents "i"
         * @param mStart First entrance for each col j; the offset of the first a_{i,j} not 0 in mRwInd
         * @param mNel First entrance for each col j; represents the total number of non null coeficients in column "j"
         * @param dMatVal First entrance for each a_{?} not 0; represents the value of a_{i,j}
         */   
        int loadSetCoverProblem(int, int, int*, int*, int*, double*);

        /**
         * Callback triggered when an integer solution is found. 
         */
        static void XPRS_CC callbackBestSolution(XPRSprob, void*);

        /**
         * Solves an SCP instance.
         */
        void solveSetCov();

        /**
         * Solves an IP using Xpress.
         */
        void solve();

        /**
         * Verifies if an optimal solution was found.
         */
        bool isOptimal();

        /**
         * Sets an external UB.
         */
        void setUB(double value);

        /**
         * Sets an initial viable solution for the SCP being solved.
         */
        void setViableSolution(vector<int> solution);

        /**
         * Prints the current problem being solved to an output file.
         */
        void writeProblem();

        /**
         * General Getters and Setters.
         */
        double getZStar(){ return _zStar; }
        double getExtLB(){ return _extLB; }
        void setZStar(double zStar){ _zStar = zStar; }
        int getNumCols() { return _nCols; }
        void freeXStar() { if (_xStar != NULL) delete[] _xStar; }
        void setXStar(double* xStar){ _xStar = xStar; }
        double* getXStar() { return _xStar; }
        int getProblemType() { return _problemType; }
        void setLBstop(bool value) { _extLBstop = value;};

};

#endif
#define SOLVER_PLI_XPRESS_H
