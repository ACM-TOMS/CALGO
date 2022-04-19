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


#include "PreSolver.h"

/** 
 * Uses Lagrangian Heuristic to find a good viable solution for the current
 * SCP.
 */
vector<int> PreSolver::tryHeuristic() {
    int nRows = _nAliveRows;
    int nCols = _nAliveCols;        
    int tj = 0, ti = 0;
    vector< vector<bool> > tempMatrix(nRows);    
    Lagrangian * g;

    if (nCols != _nCand || nRows != _nWitn) {
        cout << "[PreSolver] Creating temporary Matrix..." << endl;
        
        for (int i=0; i<_nWitn; i++) {
            if (!_rAlive[i])
                continue;

            for (int j=0; j<_nCand; j++) {
                if (!_cAlive[j])
                    continue;

                tempMatrix[ti].push_back(_matrix[i][j]);
            }

            ti++;
        }
        
        for (int j=0; j<_nCand; j++) {
            if (!_cAlive[j])
                continue;

            _redMap[tj] = j;

            tj++;
        }

    }
    else {
        tempMatrix = _matrix; 
    }
    
    g = new Lagrangian(tempMatrix, _extLB);

    g->dual();
    
    double lb = g->getLB();
    int ub = g->getUB();
    int it = g->getIterations();
    int bestIt = g->getBestUBIt();
    _heurSolution = g->getSolution();  

    std::cout << "[PreSolver] Final:" << std::endl;
    std::cout << "[PreSolver]   LB: " << lb << std::endl;
    std::cout << "[PreSolver]   UB: " << ub << std::endl;
    std::cout << "[PreSolver]   Best It: " << bestIt << std::endl;
    std::cout << "[PreSolver]   Total It: " << it << std::endl;
    
    if ((ub - lb) < 1.0) {
        _optSolution = _heurSolution;
        _foundOptimal = true;
    }
    
    if ((ub - _extLB) < 1.0) {
        _optSolution = _heurSolution;
        _foundOptimal = true;
    }

    delete g;

    return _heurSolution;
}

#ifdef XPRESS
/**
 * Prepares the necessary data in order to call the SCP solver library using
 * XPRESS.
 */
void PreSolver::solveWithXpress(bool preXpress = true, double initSol = -1.0, vector<int> sol = vector<int>()) {
    int nRows = _nAliveRows;  
    int nCols = _nAliveCols;                                                              
    int tj = 0, ti = 0;
    int numElem = 0;
    int useXpressPreSolver;    

    SolverPLIXpress* solver;
    
    int* mStart = (int*) calloc(nCols, sizeof(int));                                           
    int* mNel = (int*) calloc(nCols, sizeof(int));                                                   
    
    int startValMatrix = (3 * nCols) * (nRows / 10);                           
    if (startValMatrix < 50) {             
        startValMatrix = 50;
    }                                                 

    int* mRwInd = (int*) calloc(startValMatrix, sizeof(int));  
    double* dMatVal = (double*) calloc(startValMatrix, sizeof(double));

    for (int j=0; j<_nCand; j++) {
        if (!_cAlive[j])
            continue;
        
        mStart[tj] = numElem;
        mNel[tj] = 0;

        ti=0;
        for (int i=0; i<_nWitn; i++) {
            if (!_rAlive[i])
                continue;

            if (_matrix[i][j]) {
                mRwInd[numElem] = ti;
                dMatVal[numElem++] = 1.0;
                mNel[tj]++;

                if (startValMatrix <= numElem) {
                    startValMatrix = 14 * startValMatrix / 10;
                    if (startValMatrix > (nCols * nRows)) {
                        startValMatrix = nCols * nRows;
                    }

                    dMatVal = (double*) realloc(dMatVal, startValMatrix * sizeof(double));
                    mRwInd = (int*) realloc(mRwInd, startValMatrix * sizeof(int));
                }
            }
            ti++;
        }
            
        _redMap[tj] = j;
        tj++;
    }
    
    dMatVal = (double*) realloc(dMatVal, numElem * sizeof(double));
    mRwInd = (int*) realloc(mRwInd, numElem * sizeof(int));

    if (preXpress)
        useXpressPreSolver = SolverPLIXpress::USE_XPRESS;
    else
        useXpressPreSolver = SolverPLIXpress::NOT_USE_XPRESS;

    solver = new SolverPLIXpress(nCols, nRows, SolverPLIXpress::SET_COVER, useXpressPreSolver, mStart, mNel, mRwInd, dMatVal, _extLB);
    
    if (initSol > 0.0) {
        solver->setViableSolution(sol);
    }
    solver->solve();
    double * x = solver->getXStar();

    _optSolution.clear();
    for (int j=0; j<nCols; j++) {
        if (x[j] + SMALL_NUMBER_GAP >= 1.0) {
            _optSolution.push_back(j);
        }
    }

    _foundOptimal = solver->isOptimal(); 

    /* Frees memory */
    free(mStart);
    free(mNel);
    free(mRwInd);
    free(dMatVal);
    delete solver;

    return;
}
#endif

#ifdef GLPK
/**
 * Prepares the necessary data in order to call the SCP solver library using
 * GLPK.
 */
void PreSolver::solveWithGlpk(double initSol = -1.0, vector<int> sol = vector<int>()) {
    int nRows = _nAliveRows;  
    int nCols = _nAliveCols;                                                              
    vector< vector<bool> > reducedMatrix;
    SolverPLIGlpk* solver = new SolverPLIGlpk();
    int ti=0, tj=0;

    reducedMatrix.resize(nRows);
    for (int i=0; i<nRows; i++) {
        reducedMatrix[i].resize(nCols);
    }

    for (int i=0; i<_nWitn; i++) {
        if (!_rAlive[i])
            continue;

        tj = 0;
        for (int j=0; j<_nCand; j++) {
            if (!_cAlive[j])
                continue;

            reducedMatrix[ti][tj] = _matrix[i][j];
            
            tj++;
        }

        ti++;
    }
        
    tj = 0;
    for (int j=0; j<_nCand; j++) {
        if (!_cAlive[j])
            continue;

        _redMap[tj] = j; 
            
        tj++;
    }

    
    solver->solveSCP(reducedMatrix, sol, _extLB);
    
    vector<double> x = solver->getBestSolution();

    _optSolution.clear();
    for (int j=0; j<nCols; j++) {
        if (x[j] + SMALL_NUMBER_GAP >= 1.0) {
            _optSolution.push_back(j);
        }
    }

    _foundOptimal = solver->isOptimal(); 

    delete solver;

    return;
}
#endif

/**
 * Solves an SCP using the mode selected by the user.
 */
void PreSolver::solve(int mode) {
    
    switch (mode) {
        case 1:
#ifdef GLPK
            bitsColsReduction();
            bitsRowsReduction();
            tryHeuristic();
            if (_optSolution.size() == 0)
                solveWithGlpk(_heurSolution.size(), _heurSolution);
            else
                cout << "[PreSolver] Lagrangian Heuristic found optimal solution for SCP!" << endl;
            break;
#else
            goto badoption;
#endif
        case 2:
#ifdef GLPK
            bitsColsReduction();
            bitsRowsReduction();
            solveWithGlpk(-1);
            break;
#else
            goto badoption;
#endif
        case 3:
#ifdef XPRESS
            bitsColsReduction();
            bitsRowsReduction();
            tryHeuristic();
            if (_optSolution.size() == 0)
                solveWithXpress(true, _heurSolution.size(), _heurSolution);
            else
                cout << "[PreSolver] Lagrangian Heuristic found optimal solution for SCP!" << endl;
            break;
#else
            goto badoption;
#endif
        case 4:
#ifdef XPRESS
            bitsColsReduction();
            bitsRowsReduction();
            solveWithXpress(true, _heurSolution.size(), _heurSolution);
            break;
#else
            goto badoption;
#endif
        default:
badoption:
            cout << "[PreSolver] Bad option for solver mode. The solver mode requested is not available!" << endl;
            exit(1);
            break;
    }

    return;
}

/**
 * Constructor. Receives the ILP Matrix, the guard candidates groups
 * (which are organized by the Light AVPs where the candidates belong) 
 * and the lower bound previously obtained by the solution. Using all
 * this information, we setup the structures used when solving the ILP.
 */
PreSolver::PreSolver(const vector< vector<bool> > &matrix, vector< vector<int> > &groups, double extLB) {
    _nWitn = matrix.size(); 
    _nCand = matrix[0].size();
   
    _matrix = matrix;

    int nRows = matrix.size(); 
    int nCols = matrix[0].size();

    _extLB = extLB;

    for (int j=0; j<nCols; j++) {
        boost::dynamic_bitset<> newBitSet;
        for (int i=0; i<nRows; i++) {
            newBitSet.push_back(matrix[i][j]);
        }
        _iBitsMatrix.push_back(newBitSet);
        string s;
        to_string(newBitSet, s);

        newBitSet.clear();
    }

    _groups = groups;

    for (int j=0; j<_nCand; j++) {
        _cAlive.push_back(true);
        _redMap.push_back(j);
    }
    _nAliveCols = _nCand;
    
    for (int i=0; i<_nWitn; i++) {
        _rAlive.push_back(true);
    }
    _nAliveRows = _nWitn;

    return;
}


/**
 * Eliminates redundant rows.
 */
void PreSolver::bitsRowsReduction() {
    int cont = 0;

    /* Creates sets of bits */
    vector< boost::dynamic_bitset<> > bitsRows(_nWitn);
    for (int i=0; i<_nWitn; i++) {
        for (int j=0; j<_nCand; j++) {
            if (_cAlive[j])
                bitsRows[i].push_back(_matrix[i][j]);
        }
    }

    /* Tests the sets */
    boost::dynamic_bitset<> res;

    for (int i=0; i<_nWitn; i++) {
        if (!_rAlive[i])
            continue;

        for (int k=(i+1); k<_nWitn; k++) {
            if (!_rAlive[k])
                continue;

            res = bitsRows[i] | bitsRows[k];

            if (res == bitsRows[i]) {
                _rAlive[i] = false;
                break;
            }
            else if (res == bitsRows[k]) {
                _rAlive[k] = false;
            }
        }
    } 

    /* Counts the alive sets */
    for (int i=0; i<_nWitn; i++) {
        if (_rAlive[i])
            cont++;
    }
    _nAliveRows = cont;

    cout << "[PreSolver] Started with " << _nWitn << " constraints and finished with " << cont << endl;

}

/**
 * Eliminates redundant columns.
 */
void PreSolver::bitsColsReduction() {
    int cont = 0;

    for (int k=0; k<_groups.size(); k++) {
        bitsCheckGroup(k);
    }

    for (int j=0; j<_nCand; j++) {
        if (_cAlive[j])
            cont++;
    }
    _nAliveCols = cont;

    cout << "[PreSolver] Started with " << _nCand << " candidates and finished with " << cont << endl;
}

/**
 * Checks a group of guard candidates to find redundant columns. 
 */
void PreSolver::bitsCheckGroup(int index) {
    int cand1 = -1;
    int cand2 = -1;
    boost::dynamic_bitset<> res;
            
    for (int k=0; k<_groups.at(index).size(); k++) {
        if (!_cAlive[_groups[index].at(k)])
            continue;

        cand1 = _groups[index].at(k);

        for (int h=(k+1); h<_groups.at(index).size(); h++) {
            if (!_cAlive[_groups[index].at(h)])
                continue;

            cand2 = _groups[index].at(h);

            res = _iBitsMatrix[cand1] | _iBitsMatrix[cand2];

            if (res == _iBitsMatrix[cand1]) {
                _cAlive[cand2] = false;
            }
            else if (res == _iBitsMatrix[cand2]) {
                _cAlive[cand1] = false;
                break;
            }
        }
    }

}

/**
 * Constructs the optimal solution found using the initial indices.
 */
vector<int> PreSolver::getOptimalSolution() {
    vector<int> res;
    
    if (!_optSolution.size())
        _optSolution = _heurSolution;

    if (_nAliveCols != _nCand) {
        for (int k=0; k<_optSolution.size(); k++) {
            res.push_back(_redMap[_optSolution[k]]);
        }
    }
    else {
        res = _optSolution;
    }

    return res;
}

