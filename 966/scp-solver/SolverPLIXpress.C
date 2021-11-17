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


#include "SolverPLIXpress.h"

/**
 * Constructor for Xpress SCP solver.
 */
SolverPLIXpress::SolverPLIXpress(int nCols, int nRows, int problemType, int useXPRESS, int* mStart, int* mNel, int* mRwInd, double* dMatVal, double extLB) {
    // Initializes atributes
    _xStar = NULL;
    _nCols = nCols;
    int xpressRet = 0;

    _extLB = extLB;

    // Initializes XPRESS
    xpressRet=XPRSinit(NULL);
    if (xpressRet) errorMsg("Main [initialize XPRESS]",__LINE__,xpressRet);
	
    // Creates the problem
    xpressRet = XPRScreateprob(&_prob);
    if (xpressRet) errorMsg("Main [initialize problem]", __LINE__, xpressRet);

    // Limits maximum runtime
    xpressRet = XPRSsetintcontrol(_prob, XPRS_MAXTIME, MAX_CPU_TIME);
    if (xpressRet) errorMsg("Main [limiting running time]", __LINE__, xpressRet);

    // Enables all XPRESS tricks
    xpressRet = XPRSsetintcontrol(_prob, XPRS_PRESOLVEOPS, 36863);
    if (xpressRet) errorMsg("Main [setting PRESOLVE]", __LINE__, xpressRet);
	
    // Enables (or not) MIPPRESOLVE
    xpressRet = XPRSsetintcontrol(_prob, XPRS_MIPPRESOLVE, useXPRESS);
    if (xpressRet) errorMsg("Main [desabling MIPPRESOLVE]", __LINE__, xpressRet);

    // Enables (or not) PRESOLVE
    xpressRet = XPRSsetintcontrol(_prob, XPRS_PRESOLVE, useXPRESS);
    if (xpressRet) errorMsg("Main [desabling PRESOLVE]", __LINE__, xpressRet);

    // Defines verbose level
    xpressRet = XPRSsetintcontrol(_prob, XPRS_MIPLOG, 3); 
    if (xpressRet) errorMsg("Main [setting MIPLOG]", __LINE__, xpressRet);

    // Modifies MIPABSTOP to stop when it is really close
    xpressRet = XPRSsetdblcontrol(_prob, XPRS_MIPABSSTOP, 0.9);
    if (xpressRet) errorMsg("Main [setting MIPABSSTOP]", __LINE__, xpressRet);

    // Sets the number of threads used
    xpressRet = XPRSsetintcontrol(_prob, XPRS_THREADS, 1);
    if (xpressRet) errorMsg("Main [setting THREADS]", __LINE__, xpressRet);

    // Enables the use of heuristics
    xpressRet = XPRSsetintcontrol(_prob, XPRS_HEURSTRATEGY, 3);
    if (xpressRet) errorMsg("Main [seting HEURSTRATEGY]", __LINE__, xpressRet);

    // Sets the maximum depth where heuristics are allowed 
    xpressRet = XPRSsetintcontrol(_prob, XPRS_HEURDEPTH, 1000);
    if (xpressRet) errorMsg("Main [setting HEURSTRATEGY]", __LINE__, xpressRet);

    // Sets the frequency of heuristics
    xpressRet = XPRSsetintcontrol(_prob, XPRS_HEURFREQ, 10);
    if (xpressRet) errorMsg("Main [setting HEURSTRATEGY]", __LINE__, xpressRet);

    // Sets the maximum number of heuristic solutions
    xpressRet = XPRSsetintcontrol(_prob, XPRS_HEURMAXSOL, 1000000);
    if (xpressRet) errorMsg("Main [setting HEURMAXSOL]", __LINE__, xpressRet);

    // Sets the maximum number of visited nodes
    xpressRet = XPRSsetintcontrol(_prob, XPRS_HEURNODES, 100000);
    if (xpressRet) errorMsg("Main [setting HEURNODES]", __LINE__, xpressRet);

    // Sets the maximum number of Threads for heuristics
    xpressRet = XPRSsetintcontrol(_prob, XPRS_HEURTHREADS, 0);
    if (xpressRet) errorMsg("Main [setting HEURTHREADS]", __LINE__, xpressRet);

    // Enables the use of new cuts
    xpressRet = XPRSsetintcontrol(_prob, XPRS_CUTSTRATEGY, 3);
    if (xpressRet) errorMsg("Main [setting CUTSTRATEGY]", __LINE__, xpressRet);

    // Sets the maximum depth where cuts can be created
    xpressRet = XPRSsetintcontrol(_prob, XPRS_CUTDEPTH, 100);
    if (xpressRet) errorMsg("Main [setting CUTDEPTH]", __LINE__, xpressRet);

    // Sets the frequency of cuts
    xpressRet = XPRSsetintcontrol(_prob, XPRS_CUTFREQ, 10);
    if (xpressRet) errorMsg("Main [setting CUTFREQ]", __LINE__, xpressRet);

    // Sets the DUALGRADIENT
    xpressRet = XPRSsetintcontrol(_prob, XPRS_DUALGRADIENT, 3);
    if (xpressRet) errorMsg("Main [setting DUALGRADIENT]", __LINE__, xpressRet);

    _problemType = problemType;
	
    // Loads the model
    switch (problemType) {
        case SET_COVER:
            _zStar = XPRS_PLUSINFINITY;

            xpressRet = loadSetCoverProblem(nCols, nRows, mStart, mNel, mRwInd, dMatVal);
            if (xpressRet) errorMsg("Main [loading problem]", __LINE__, xpressRet);

            // Sets callback function
            xpressRet = XPRSsetcbintsol(_prob, &SolverPLIXpress::callbackBestSolution, this);
            if (xpressRet) errorMsg("Main [XPRSsetcbintsol]", __LINE__, xpressRet);

            break;
        default: errorMsg("Main [problem type not defined or suported]", __LINE__, xpressRet);
    }
}

void SolverPLIXpress::errorMsg(const string subRotine, int line, int errCode) {

    cout << "[ERR] subroutine " << subRotine << " has failed on line " << line << endl;

    if (errCode != -1) cout << "  with error code " << errCode << endl;
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
int SolverPLIXpress::loadSetCoverProblem(int nCols, int nRows, int* mStart, int* mNel, int* mRwInd, double* dMatVal) {
    int xpressRet;
    double* coefFncObj = new double[nCols]; // coeficients of the objective function 
    double* rhs = new double[nRows]; // First entrance for each restriction; 
                                     //   represents the rhs of restriction
    char* qRType = new char[nRows]; // First entrance for each constraint; Can be 
                                    //   'L', 'G' or 'E' (meaning '<=', '>=' or '=')
    double* dlb = new double[nCols]; // lower bound for each variable
    double* dub = new double[nCols]; // upper bound for each variable
    int* mGCols = new int[nCols]; // First entrance for each variable; represents
                                  //   the variable intex
    char* qGType = new char[nCols]; // First entrance for each variable; represents
                                    //   the type of the variable

    for(int i=0; i < nCols; i++){
        dlb[i] = 0.0;
        dub[i] = 1.0;
        mGCols[i] = i;
        qGType[i] = 'B'; // all variables are booleans
        coefFncObj[i] = 1.0;
    }

    for(int i=0; i < nRows; i++){
        rhs[i] = 1;
        qRType[i] = 'G';
    }

    // Loads model
    xpressRet = XPRSloadglobal(_prob, "GALLERY", nCols, nRows, qRType, rhs, NULL, coefFncObj, \
            mStart, mNel, mRwInd, dMatVal, dlb, dub, nCols, 0, qGType, mGCols, \
            NULL, NULL, NULL, NULL, NULL);

    // Frees memory
    delete[] coefFncObj;
    delete[] rhs;
    delete[] qRType;
    delete[] dlb;
    delete[] dub;
    delete[] mGCols;
    delete[] qGType;

    return xpressRet;	
}

/**
 * Callback triggered when an integer solution is found. 
 */
void XPRS_CC SolverPLIXpress::callbackBestSolution(XPRSprob prob, void* myObj) {
    SolverPLIXpress * solver = (SolverPLIXpress *) myObj;

    int xpressRet;
    int nCols;
    double objVal;
    double* auxCoefFncObj;

    xpressRet = XPRSgetintattrib(prob, XPRS_COLS, &nCols);
    if (xpressRet) solver->errorMsg("callbackBestSolution [XPRSgetintattrib - XPRS_COLS]", __LINE__, xpressRet);

    xpressRet = XPRSgetdblattrib(prob, XPRS_LPOBJVAL, &objVal);
    if (xpressRet) solver->errorMsg("callbackBestSolution [XPRSgetdblattrib - XPRS_LPOBJVAL]", __LINE__, xpressRet);

    if ((objVal < solver->getZStar() && solver->getProblemType() == SET_COVER)) {

        solver->freeXStar();

        auxCoefFncObj = new double[solver->getNumCols()];

        xpressRet = XPRSgetmipsol(prob, auxCoefFncObj, NULL);
        if (xpressRet) solver->errorMsg("callbackBestSolution [XPRSgetsol]", __LINE__, xpressRet);

        solver->setXStar(auxCoefFncObj);
        solver->setZStar(objVal);

        xpressRet = XPRSsetdblcontrol(prob, XPRS_MIPABSCUTOFF, objVal);
        if (xpressRet) solver->errorMsg("callbackBestSolution [XPRSsetdblcontrol", __LINE__, xpressRet);

        if ((objVal - solver->getExtLB()) < 0.999) {
            cout << "[SolverPLIXpress] XPRESS can already stop searching!!" << endl;
            XPRSinterrupt(prob, XPRS_STOP_USER);
            solver->setLBstop(true);
        }
    }
}

/**
 * Solves an SCP instance.
 */
void SolverPLIXpress::solveSetCov(){
    int xpressRet = XPRSminim(_prob,"g"); 

    _extLBstop = false;

    if (xpressRet) errorMsg("Main [XPRSminim]", __LINE__, xpressRet);

    int nCols;
    double objVal;
    double* auxCoefFncObj;

    xpressRet = XPRSgetintattrib(_prob, XPRS_COLS, &nCols);
    if (xpressRet) errorMsg("solve [XPRSgetintattrib - XPRS_COLS]", __LINE__, xpressRet);

    xpressRet = XPRSgetdblattrib(_prob, XPRS_LPOBJVAL, &objVal);
    if (xpressRet) errorMsg("solve [XPRSgetdblattrib - XPRS_LPOBJVAL]", __LINE__, xpressRet);

    int mipSols;
    xpressRet = XPRSgetintattrib(_prob, XPRS_MIPSOLS, &mipSols);
    if (xpressRet) errorMsg("callbackBestSolution [XPRSgetintattrib - XPRS_MIPSOLS]", __LINE__, xpressRet);

    if (objVal < getZStar()) {
        freeXStar();
        auxCoefFncObj = new double[getNumCols()];

        xpressRet = XPRSgetmipsol(_prob, auxCoefFncObj, NULL);
        if (xpressRet) errorMsg("solve [XPRSgetsol]", __LINE__, xpressRet);

        setXStar(auxCoefFncObj);
        setZStar(objVal);
    }
}

/**
 * Solves an IP using Xpress.
 */
void SolverPLIXpress::solve(){
    switch(_problemType) {
        case SET_COVER:
            solveSetCov();
            break;
        default:
            break;	
    }

    return;
}

/**
 * Verifies if an optimal solution was found.
 */
bool SolverPLIXpress::isOptimal() {
    int statusNumber;
    int xpressRet = XPRSgetintattrib(_prob, XPRS_MIPSTATUS, &statusNumber);

    if (xpressRet) 
        errorMsg("XPRS_MIPSTATUS - XPRSgetintattrib",__LINE__,xpressRet);
	
    if (statusNumber == XPRS_MIP_OPTIMAL || _extLBstop) {
        return true;
    } 

    return false;
}

/**
 * Sets an external UB.
 */
void SolverPLIXpress::setUB(double value){
    int xpressRet = XPRSsetdblcontrol(_prob, XPRS_MIPABSCUTOFF, value);
    if (xpressRet) 
        errorMsg("setUB [XPRSsetdblcontrol", __LINE__, xpressRet);
}

/**
 * Sets an initial viable solution for the SCP being solved.
 */
void SolverPLIXpress::setViableSolution(vector<int> solution){
    double * vet = (double *) malloc(_nCols * sizeof(double));
    int status;

    for (int i=0; i<_nCols; i++) {
        vet[i] = 0.0;
    }

    for (int i=0; i<solution.size(); i++) {
        vet[solution.at(i)] = 1.0;
    }

    int xpressRet = XPRSloadmipsol(_prob, vet, &status);
    if (xpressRet) errorMsg("setUB [XPRSloadmipsol", __LINE__, xpressRet);

    if (!status)
        cout << "[SolverPLIXpress] Heuristic solution was accepted!" << endl;
    else
        cout << "[SolverPLIXpress] Xpress didn't accept heuristic solution!" << endl;

    free(vet);

    return;
}

/**
 * Prints the current problem being solved to an output file.
 */
void SolverPLIXpress::writeProblem() {
    // Save file with original ILP model
    int xpressRet = XPRSwriteprob(_prob, "problem", "l");
    if (xpressRet) errorMsg("Main [XPRSwriteprob]", __LINE__, xpressRet);
}
