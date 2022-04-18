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


#ifndef SOLVER_PLI_GLPK_H
#define SOLVER_PLI_GLPK_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <glpk.h>

using namespace std;

class SolverPLIGlpk {
    private:
        vector<double> _bestSolution;
        bool _isOptimal;

    public:
        /**
         * Solves SCP using GLPK API.
         */
        int solveSCP (vector< vector<bool> > matrix, vector<int> initSol, double extLB);

        /**
         * Returns best solution found.
         */
        vector<double> getBestSolution() { return _bestSolution; };

        /**
         * Return true if the solution found by GLPK is optimal.
         */
        bool isOptimal() { return _isOptimal; }
};

#endif
#define SOLVER_PLI_GLPK_H
