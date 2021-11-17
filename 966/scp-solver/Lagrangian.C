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

#include "Lagrangian.h"

/** 
 * Sets the initial Lagrangian Multipliers. 
 */
void Lagrangian::setInitialMultipliers() {
    for (int i=0; i<_nRows; i++) {
        if (!_activeConstraints[i])
            continue;
        
        _u[i] = _initialMulti;
    }
}

/** 
 * Updates Lagrangian Multipliers.
 */
void Lagrangian::updateMultipliers() {
    for (int i=0; i<_nRows; i++) {
        if (!_activeConstraints[i])
            continue;
        
        _u[i] += _mi * _s[i];
        if (_u[i] < 0) {
            _u[i] = 0;
        }
    } 
}

/** 
 * Updates Lagrangian step size.
 */
void Lagrangian::updateStep() {
    double sum = 0;

    for (int i=0; i<_nRows; i++) {
        if (!_activeConstraints[i])
            continue;
        
        sum += (_s[i] * _s[i]);
    }

    _mi = _paramStep * ((double)_ub - _lb) / sqrt(sum);
}

/** 
 * Computes Lagrangian Costs.
 */
void Lagrangian::calculateCosts() {
    for (int j=0; j<_nCols; j++) {
        if (!_activeVariables[j])
            continue;
        
        double sum = 0.0;
        for (int k=0; k<_colsList[j].size(); k++) {
            int row = _colsList[j].at(k);
            if (!_activeConstraints[row])
                continue;
        
            sum += _u[row];
        }
        _c[j] = 1 - sum; 
    }
}

/** 
 * Computes Subgradient in order to update the current Lagrangian Multipliers.
 */
void Lagrangian::calculateSubgradient() {
    for (int i=0; i<_nRows; i++) {
        if (!_activeConstraints[i])
            continue;

        double sum = 0.0;
        for (int k=0; k<_rowsList[i].size(); k++) {
            int col = _rowsList[i].at(k);
            if (!_activeVariables[col])
                continue;
            
            sum += _xl[col];
        }
        _s[i] = 1 - sum;
    }
}

/**
 * Destructor. Frees structures memory.
 */
Lagrangian::~Lagrangian() {
    free(_u);
    free(_c);
    free(_s);
    free(_xl);
    free(_xh);
    free(_activeConstraints);
    free(_activeVariables);
}

/**
 * Constructor. Receives the ILP matrix and an external lower bound computed
 * previously.
 */
Lagrangian::Lagrangian(const vector< vector<bool> > &matrix, double extLB) {
    srand (time(NULL));
    
    _extLB = extLB;
    
    _matrix = matrix;
    _nRows = matrix.size(); 
    _nCols = matrix[0].size();

    _mi = 0.0;
    _u = (double *) calloc(_nRows, sizeof(double));
    _c = (double *) calloc(_nCols, sizeof(double));
    _s = (double *) calloc(_nRows, sizeof(double));
    _xl = (bool *) calloc(_nCols, sizeof(bool));
    _xh = (bool *) calloc(_nCols, sizeof(bool));
    
    _lb = 0.0;
    _ub = _nCols;

    _sameLB = 0;
    _sameUB = 0;

    _paramGreedy = 2;
    _paramStep = PARAM_STEP;
    
    cout << "[Lagrangian] Param step: " << _paramStep << endl; 

    _sameLBMaxIt = _matrix.size();
    _sameUBMaxIt = 5 * _matrix.size();
   
    _initialMulti = 0.5;

    _activeConstraints = (bool *) calloc(_nRows, sizeof(bool));
    _activeVariables = (bool *) calloc(_nCols, sizeof(bool));
    for (int i=0; i<_nRows; i++) {
       _activeConstraints[i] = true; 
    }
    for (int j=0; j<_nCols; j++) {
       _activeVariables[j] = true; 
    }

    /* Fill in visibility vector */
    _rowsList.resize(_nRows);
    _colsList.resize(_nCols);
    for (int i=0; i<_nRows; i++) {
        for (int j=0; j<_nCols; j++) {
            if (_matrix[i][j]) {
                _rowsList[i].push_back(j);
                _colsList[j].push_back(i);
            }
        }
    }

}

/**
 * Computes lagrangian relaxation.
 */ 
double Lagrangian::relaxation() {
    double value = 0.0;
    
    _sameLB++;
    
    calculateCosts();

    for (int j=0; j<_nCols; j++) {
        if (!_activeVariables[j])
            continue;
        
        if (_c[j] < 0) {
            _xl[j] = 1;
            value += _c[j];
        }
        else {
            _xl[j] = 0;
        }
    }

    for (int i=0; i<_nRows; i++) {
        if (!_activeConstraints[i])
            continue;
        
        value += _u[i];
    } 

    if (_lb < (value + _fixedInOne.size())) {
        _lb = value + _fixedInOne.size();
        _sameLB = 0;
    }

    return value;
}

/**
 * Iterativelly computes new lower and upper bounds for the ILP.
 */ 
double Lagrangian::dual() {
    int ret;
    _iterations = 0;
    int paramChangeIt = _sameUBMaxIt / 10;    

    cout << "[Lagrangian] param change: " << paramChangeIt << endl;

    setInitialMultipliers();
    heuristicNew(1);

    while (_iterations < MAX_DUAL_IT && (_lb + 1) <= _ub && (_extLB + 1) <= _ub && _sameLB <= _sameLBMaxIt && _sameUB <= _sameUBMaxIt) {
        _v = relaxation();

        if (_sameLB == 0) {
            ret = heuristicNew(_paramGreedy);
            if (ret == -1) {
                _lb = _ub;
                goto _end;
            }
        }
        else {
            _sameUB++;
        }

        problemReduction();

        calculateSubgradient();
        updateStep();
        updateMultipliers();
   
        _iterations++;

        if (!(_iterations % paramChangeIt)) {
            _paramStep *= 0.9;
        }
    }

_end:
    if (_lb > (_ub - 1))
        _lb = _ub;

    return _lb;
}

/**
 * Based on some ranking, chooses next variable to be in the guards set.
 */
int Lagrangian::chooseNextSet(int * grades, int paramGreedy) {
    int j;
    int newX = -1;

    if (paramGreedy == 1) {
        int maxGrade = -1;

        for (j=0; j<_nCols; j++) {
            if (!_activeVariables[j])
                continue;
            
            if (maxGrade < grades[j]) {
                newX = j;
                maxGrade = grades[j];
            }
        } 
    }
    else if (paramGreedy == 2) {
        double minGrade = 99999999.9;
        double minC = 0.0;

        for (j=0; j<_nCols; j++) {
            if (_c[j] < minC) {
                minC = _c[j];
            }
        }

        minC = 0.1 - minC;

        for (j=0; j<_nCols; j++) {
            if (!_activeVariables[j])
                continue;
            
            if (grades[j] == 0)
                continue;

            if (minGrade > ((minC + _c[j]) / grades[j])) {
                newX = j;
                minGrade = ((minC + _c[j]) / grades[j]);
            }
        }        
    }

    return newX;
}

/**
 * Verifies if the given solution is a viable one.
 */
bool Lagrangian::isViableSolutionNew(vector<int> solution) {
    tr1::unordered_set<int> covered;

    for (int k=0; k<_solvedCons.size(); k++) {
        covered.insert(_solvedCons[k]);
    } 

    for (int k=0; k<solution.size(); k++) {
        for (int h=0; h<_colsList[solution[k]].size(); h++) {
            int row = _colsList[solution[k]].at(h);
            covered.insert(row);
        }
    }

    if (covered.size() != _nRows)
        return false;


    return true; 
}

/**
 * Cleans viable solution, erasing unnecessary guards.
 */
void Lagrangian::cleaningViable(vector<int> &solution) {
    int index = 0;
    vector<int> sorted;
    vector<int> bingo;

    for (int k=0; k<solution.size(); k++) {
        bingo.push_back(k);
    }

    while (bingo.size()) {
        /* Generates secret number */
        index = rand() % bingo.size();

        sorted.push_back(solution[bingo[index]]);

        bingo.erase(bingo.begin() + index);
    }    
    
    for (int k=0; k<sorted.size(); k++) {
        vector<int> trial;
        trial = sorted;
        trial.erase(trial.begin() + k);
        if (isViableSolutionNew(trial)) {
            sorted = trial;
            k--;
        }
    }

    solution = sorted;
}

/**
 * Heuristic that finds a viable solution to the SCP from an existing
 * lagrangian relaxation result.
 */
int Lagrangian::heuristicNew(int paramGreedy) {
    tr1::unordered_set<int> rest;
    vector<int> solution;
    int * grades = (int *) calloc(_nCols, sizeof(int));
    bool feasibleProblem = true;
    
    _sameUB++;

    for (int i=0; i<_nRows; i++) {
        if (_activeConstraints[i]) {
            rest.insert(i);
        }
    }

    while (rest.size() > 0) {
        feasibleProblem = false;
        
        for (int j=0; j< _nCols; j++) {
            grades[j] = 0;
        }

        for (tr1::unordered_set<int>::iterator it = rest.begin(); it != rest.end(); ++it) {
            for (int k=0; k<_rowsList[*it].size(); k++) {
                int col = _rowsList[*it].at(k);
                if (_activeVariables[col]) {
                    grades[col]++;
                    feasibleProblem = true;
                }
            }
        }

        if (!feasibleProblem) {
            free(grades);
            return -1;
        }


        int newX = chooseNextSet(grades, paramGreedy);

        _xh[newX] = true;
        solution.push_back(newX);
    
        /* Updates restrictions covered */
        for (int k=0; k<_colsList[newX].size(); k++) {
            int row = _colsList[newX].at(k);
            rest.erase(row);
        }
    }
    
    /* Removes unnecessary guards */
    cleaningViable(solution);

    if (_ub > (solution.size() + _fixedInOne.size())) {
        _ub = solution.size() + _fixedInOne.size();
        _sameUB = 0;
        _bestUBIt = _iterations;

        _guards = solution;
        _guards.insert(_guards.end(), _fixedInOne.begin(), _fixedInOne.end() );
    }
    
    free(grades);

    return _ub;
}

/**
 * If possible, removes variables and/or constraints from the original SCP.
 */
void Lagrangian::problemReduction() {
    int cont_1 = 0;
    int cont_0 = 0;
    
    if (_v < (_ub - 2)) 
        return;

    for (int j=0; j<_nCols; j++) {
        if (!_activeVariables[j])
            continue;
        
        if (_xl[j]) {
            if ((_v - _c[j]) > (_ub - 1)) {
                _activeVariables[j] = false;
                _fixedInOne.push_back(j);
                for (int i=0; i<_nRows; i++) {
                    if (!_activeConstraints[i])
                        continue;
                    
                    if (_matrix[i][j]) {
                        _activeConstraints[i] = false;
                        _solvedCons.push_back(i);
                    }
                }
                cont_1++;
            }
        }
        else {
            if ((_v + _c[j]) > (_ub - 1)) {
                _activeVariables[j] = false;
                _fixedInZero.push_back(j);
                cont_0++;
            }
        }
    }

    return;
}

/**
 * Gets current LB value.
 */
double Lagrangian::getLB() { 
    return _lb;
}

/**
 * Gets current UB value.
 */
int Lagrangian::getUB() { 
    return _ub;
}

/**
 * Gets iteration where the best UB was found.
 */
int Lagrangian::getBestUBIt() { 
    return _bestUBIt;
}

/**
 * Gets number of iterations.
 */
int Lagrangian::getIterations() { 
    return _iterations;
}

/**
 * Gets optimal solution.
 */
vector<int> Lagrangian::getSolution() {
    return _guards;
}

/**
 * Sets optimal solution.
 */
void Lagrangian::setOptimal(int * solution, int size) {
    for (int k=0; k<size; k++) {
        _optimal.push_back(solution[k]);
    }
}

