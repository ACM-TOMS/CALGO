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

#ifndef LAGRANGEAN_H
#define LAGRANGEAN_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <tr1/unordered_set>
#include <boost/dynamic_bitset.hpp>
#include <math.h>
#include <string.h>

#define MAX_DUAL_IT 100000
#define PARAM_STEP 0.05

using namespace std;

class Lagrangian {
    private:
        vector< vector<bool> > _matrix; // ILP matrix
        vector< vector<int> > _rowsList; // adjacency list
        vector< vector<int> > _colsList; // adjacency list
        int _nRows;
        int _nCols;
        double * _u;  // Lagrangian Multipliers
        double * _c;  // lagrangean costs -> sum_{c}(1 - sum_{w}(a_wc * u)), sendo a_wc 1 ou 0
        double * _s;  // subgradient
        bool * _xl; // solution of relaxation
        bool * _xh; // solution of heuristic procedure
        double _mi; // size of the dual subgradient step
 
        vector<int> _guards;
        vector<int> _optimal; // optimal solution

        bool * _activeConstraints;
        bool * _activeVariables;
        vector<int> _fixedInOne;
        vector<int> _fixedInZero;
        vector<int> _solvedCons;

        double _v; // value of the last lagrangean relaxation iteration
        double _lb;
        int _ub;
        double _extLB;

        int _sameLB;
        int _sameUB;
        int _iterations;
        int _bestUBIt;
    
        unsigned int _sameLBMaxIt;
        unsigned _sameUBMaxIt;
        double _initialMulti;

        int _paramGreedy;
        double _paramStep;

        /** 
         * Sets the initial Lagrangian Multipliers. 
         */
        void setInitialMultipliers();

        /** 
         * Updates Lagrangian Multipliers.
         */
        void updateMultipliers();

        /** 
         * Updates Lagrangian step size.
         */
        void updateStep();

        /** 
         * Computes Lagrangian Costs.
         */
        void calculateCosts();

        /** 
         * Computes Subgradient in order to update the current Lagrangian
         * Multipliers.
         */
        void calculateSubgradient();

        /**
         * Based on some ranking, chooses next variable to be in the guards
         * set.
         */
        int chooseNextSet(int * grades, int paramGreedy);

        /**
         * Verifies if the given solution is a viable one.
         */
        bool isViableSolution(vector<int> solution);

        /**
         * Cleans viable solution, erasing unnecessary guards.
         */
        bool isViableSolutionNew(vector<int> solution);

        /**
         * Cleans viable solution, erasing unnecessary guards.
         */
        void cleaningViable(vector<int> &solution);

    public:
        /**
         * Constructor. Receives the ILP matrix and an external lower bound
         * computed previously.
         */
        Lagrangian(const vector< vector<bool> > &matrix, double extLB = 0.0);
                
        /**
         * Destructor. Frees structures memory.
         */
        ~Lagrangian();
       
        /**
         * Computes lagrangian relaxation.
         */ 
        double relaxation();
        
        /**
         * Iterativelly computes new lower and upper bounds for the ILP.
         */ 
        double dual();
        
        /**
         * Heuristic that finds a viable solution to the SCP from an existing
         * lagrangian relaxation result.
         */
        int heuristicNew(int paramGreedy);
        
        /**
         * If possible, removes variables and/or constraints from the original
         * SCP.
         */
        void problemReduction();

        /**
         * Gets current LB value.
         */
        double getLB();

        /**
         * Gets current UB value.
         */
        int getUB();

        /**
         * Gets iteration where the best UB was found.
         */
        int getBestUBIt();

        /**
         * Gets number of iterations.
         */
        int getIterations();

        /**
         * Gets optimal solution.
         */
        vector<int> getSolution();

        /**
         * Sets optimal solution.
         */
        void setOptimal(int * solution, int size);
};

#endif
