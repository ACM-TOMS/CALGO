/*****************************************************************************
 * This code is part of Art Gallery Solver (AGSol) Package, which aims the 
 * resolution of the Art Gallery Problem With Point Guards.
 *
 * This software version (1.0.2) has been tested under and is compatible with 
 * CGAL 3.9 and GLPK 4.52.
 *
 * Authors:
 *  Davi Colli Tozoni - davi.tozoni@gmail.com
 *  Marcelo Castilho Couto - cutomarcelo@gmail.com
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


#include "SolverPLIGlpk.h"

double * x;
double _extLB;
int _nCols;

/**
 * Callback that provides an initial viable solution to the solver and also stops whenever the best solution 
 * equals the external lower bound given.
 */
void callback(glp_tree *tree, void *info) {
    if(glp_ios_reason(tree) == GLP_IHEUR && glp_ios_curr_node(tree) == 1) {
        if (x[0] != -1.0) {
            cout << "[SolverPLIGlpk] Providing GLPK with an initial viable solution" << endl;
            
            int res = glp_ios_heur_sol(tree, x);
        
            if (!res)
                cout << "[SolverPLIGlpk] Solution accepted" << endl;
        }
    }
    else if (glp_ios_reason(tree) == GLP_IBINGO) {
        glp_prob *ilp = glp_ios_get_prob(tree);
        double value = glp_mip_obj_val(ilp);
        if (value <= _extLB) {
            cout << "[SolverPLIGlpk] GLPK found best solution possible" << endl;
            glp_ios_terminate(tree);
        }
    }
}

/**
 * Solves SCP using GLPK API.
 */
int SolverPLIGlpk::solveSCP (vector< vector<bool> > matrix, vector<int> initSol, double extLB) {
    glp_prob *ilp;
    int * ia, * ja;
    double * ar;
    int cont = 0;
    double value;
    int status;
    int nCols = matrix[0].size();
    int nRows = matrix.size();
  
    _nCols = nCols;

    /* Saving heuristic solution to help GLPK */
    if (initSol.size()) {
        x = (double *) calloc (nCols+1, sizeof(double));
        
        for (int i=0; i<(nCols+1); i++) {
            x[i] = 0.0;
        }
         
        for (int i=0; i<initSol.size(); i++) {
            x[initSol[i]+1] = 1.0;
        }
    }
    else {
        x = (double *) calloc (1, sizeof(double));
        x[0] = -1.0;
    }

    _extLB = extLB;

    /* Constructs sparse matrix */
    
    ia = (int *) calloc (nCols * nRows + 1, sizeof(int));
    ja = (int *) calloc (nCols * nRows + 1, sizeof(int));
    ar = (double *) calloc (nCols * nRows + 1, sizeof(double));

    for (int i=0; i<nRows; i++) {
        for (int j=0; j<nCols; j++) {
            if (matrix[i][j]) {
                ia[cont + 1] = i+1;
                ja[cont + 1] = j+1;
                ar[cont + 1] = 1.0;
                cont++;
            }
        }
    }
    
    /* Creates problem */
    ilp = glp_create_prob();
    glp_set_prob_name(ilp, "SCP");
    glp_set_obj_dir(ilp, GLP_MIN);
    glp_add_rows(ilp, nRows);
    glp_add_cols(ilp, nCols);

    /* Set rows */
    for (int i=0; i<nRows; i++) {
        glp_set_row_bnds(ilp, i+1, GLP_LO, 1.0, 0.0);
    }
    
    /* Set columns */
    for (int j=0; j<nCols; j++) {
        glp_set_obj_coef(ilp, j+1, 1.0);
        glp_set_col_kind(ilp, j+1, GLP_BV);
    }

    glp_load_matrix(ilp, cont, ia, ja, ar);
    
    /* Configures execution */ 
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.cb_func = callback;
    glp_simplex(ilp, NULL);
    int err = glp_intopt(ilp, &parm);
    
    /* Get the results */
    value = glp_mip_obj_val(ilp);
    status = glp_mip_status(ilp);
    cout << "[SolverPLIGlpk] Result: " << value << endl;

    if (status == GLP_OPT || status == GLP_FEAS) {
        _isOptimal = true;
    }
    else {
        _isOptimal = false;
    }

    /* Saves best solution found */
    _bestSolution.clear();
    for (int j=0; j<nCols; j++) {
       double elem = glp_mip_col_val(ilp, j+1);
       _bestSolution.push_back(elem);
    }

    glp_delete_prob(ilp);
    free(ia);
    free(ja);
    free(ar);
    free(x);

    return 0;
}
