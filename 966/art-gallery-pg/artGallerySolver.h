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


#ifndef ART_GALLERY_DAVI_H
#define ART_GALLERY_DAVI_H

#define AGPFC_FREQ 1
#define SMALL_NUMBER_GAP 0.00000001
#define MAX_ITERATIONS 30

#define ALL_VERTICES 1
#define CONVEX_VERTICES 2
#define CHWA_POINTS 3
#define CHWA_POINTS_EXTENDED 4

#include "PolygonExt.h"
#include "PolygonWithHolesExt.h"
#include "AVPLightGrid.h"
#include "polAlgorithms.h"
#include "PreSolver.h"
#include "AuxGallery.h"

#include <CGAL/copy_n.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_set_2.h>

#include <limits>
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>
#include <tr1/unordered_map>

typedef CGAL::Polygon_set_2<Extended_kernel, std::list<Point> > PolygonSet;
typedef CGAL::Creator_uniform_2<RT, Point> Creator;
typedef Polygon::Edge_const_iterator Edge_iterator;
typedef Polygon::Vertex_const_iterator Vertex_iterator;
typedef Polygon::Vertex_circulator Vertex_circulator;
typedef Polygon::Edge_const_circulator Edge_circulator;
typedef CGAL::Vector_2<Extended_kernel> Vector;

using CGAL::ORIGIN;

class ArtGallerySolver {

    private:
        /* Used in time calculations */
        struct rusage ruse;
        
        /* Discretization technique for the Witness Set */
        int _mode;
        char * _sMode;
        char * _nameLog;

        /* Time informations for log usage */
        double _procTime;
        double _initDisTime;
        double _insertDisTime;
        double _selectGuardCandTime;
        double _initSolverTime;
        double _runSolverTime;
        double _scpResolTime;
        double _selectNewDisTime;

        /* Solver Blackbox */
        PreSolver * _solverBox;
        int _solverMode;
        
        /* SCP matrix */
        vector< vector<bool> > _matrix;
        
        /* Vertical and Horizontal iterations of the procedure */
        int _mainLoopIterations;
        int _iteration;
        
        /* Lower and Upper Bound */
        int _lowerBound;
        int _upperBound;
        
        /* Polygon informations */
        string _polName; 
        PolygonWithHolesExt _polygon;
        Polygon _pbox;
        
        /* Informations about coverage */
        PolygonSet _notCovered;
        PolygonSet _notCoveredFirst;
        PolygonSet _polCovered;
        RT _areaMiss;
        RT _areaMissFirst;
        RT _remaining;
        RT _zero;
        
        /* Flags */
        int _useXPRESS;
        bool _loadOk;

        void addIterationGridPoints();
        RT solveIteration();

        /* Current solution */
        std::vector<Point> _solution;

        /* Guards information */
        std::vector<Point> _guards;
        std::vector<Point> _guardsFirst;
        std::vector<PolygonExt> _visGuards;
        
        /* Guard candidates information */
        std::vector<PolygonExt> _visGuardCandidates;
        std::vector<Point> _guardCandidates;
        
        /* Arrangemente information */
        vector<Point> _hiddenPoints;
        vector<Point> _hiddenPointsFirst;
        vector<Point> _newHiddenPoints;
        vector<PolygonExt> _visHiddenPoints;
        vector<PolygonExt> _visNewHiddenPoints;
        std::vector<PolygonExt> _arrangement;
        std::vector<PolygonExt> _arrangementFirst;
        std::vector<PolygonExt> _lightAVPs;
        std::vector<PolygonExt> _lightAVPsFirst;
        std::vector< vector<int> > _groups; 
        AVPLightGrid* _guardsGrid;
        std::vector<Point> _gridPoints;
 
        /* Optimization structures (hashtables, cache, etc) */
        map<Point, bool> _witnessTest;
        tr1::unordered_map<pair<Point, Point>, bool, PointPointHash> _visibilityTest;
        tr1::unordered_map<Point, PolygonExt, PointHash> _visCache;

        /* Tamanhos iniciais e finais do conjunto de testemunhas */
        int _initWitnSize;
        int _finalWitnSize;
        int _IPSolved;
        int _maxHorIt;
        int _initCandSize;
        int _finalCandSize;

        bool _isLastOptimal;

        /**
         * Converts mode from string to integer.
         */
        int stringModeToInt(char * mode);

        /**
         * Selects the initial witnesses based on the discretization technique used.
         */
        vector<Point> selectInitialDisPoints();

        /**
         * Based on the uncovered regions, selects new witnesses.
         */
        vector<Point> selectNewDisPoints();

        /**
         * Selects guard candidates for the next iteration. The guard candidates are
         * all the vertices of the polygon and also the vertices of light AVPs obtained 
         * from the visibility arrangement.
         */
        void selectGuardCandidates();

        /**
         * Adds each new witness and its visibility polygon in the Structure.
         */
        void insertDisPointsToArrange(vector<Point> dis);

        /**
         * Initializes the matrix used in SCP resolution.
         */
        void initSolver(PolygonWithHolesExt pol, vector<Point> guardCandidates);

        /**
         * Solves an AGPFC instance using an iterative procedure.
         */
        void runArt();

        /**
         * Iterates in AGPFC algorithm and prints information.
         */
        void stepArt();


    public:
        /** This function implements the full algorithm on a higher level code.
         * The algorithm follows an iterative process of obtaining lower and
         * upper bounds that lead to an optimal solution. These bounds are
         * drawn up through the resolution of instances of the aforementioned
         * AGPW and AGPFC problems using Integer Linear Programming (ILP)
         * techniques.
         */
        int solveProblem(PolygonWithHolesExt pol);

        /**
         * Loads a polygon from file.
         */
        PolygonWithHolesExt loadPol(char* name);

        /**
         * Prints the polygon in .pol (CGAL) format.
         */
        void printPol(PolygonWithHolesExt pol); 

        /**
         * Prints information about the solving process in file.
         */
        void writeLog(char* idFile, std::ostream* logFile);

        /**
         * Prints information about the partial solving process in file.
         */
        void writePartialLog();

        /**
         * Prints best solution found in file.
         */
        void writeBestSol(char* idFile, std::ostream* solFile);

        /* Getters and Setters */
        int getCardinality();
        PolygonSet getNotCovered() { return _notCovered; }
        PolygonSet getCovered() { return _polCovered; }
        std::vector<Point> getGuards() { return _guards; }
        std::vector<Point> getGuardCandidates() { return _guardCandidates; }
        RT getAreaMiss() { return _areaMiss; }
        int getIteration() { return _iteration; }
        bool isGallerySolved() { return (_areaMiss == _zero); }
        bool isLoadOk() { return _loadOk; }
        void setMode(char * mode);
        void setSolverMode(int mode) { _solverMode = mode; }
        void setLogFile(char * nameLog) { _nameLog = nameLog; }
};

#endif
