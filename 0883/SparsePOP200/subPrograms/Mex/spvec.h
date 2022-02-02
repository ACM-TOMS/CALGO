/* -------------------------------------------------------------

This file is a component of SparsePOP
Copyright (C) 2007 SparsePOP Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */

#ifndef _spvec_
#define _spvec_

#include "global.h"
using namespace std;

void minkovsum(/*IN*/class spvec_array & vecs1,class spvec_array & vecs2,/*OUT*/class spvec_array & rvecs);
void minkovsum(/*IN*/class spvec_array & vecs,/*OUT*/class spvec_array & rvecs);
void minkSumMinusDiagPart(/*IN*/class spvec_array & vecs,/*OUT*/class spvec_array & rvecs);
void convert_var_pattern(vector<int> pattern,class spvec_array & vecs);

class spvec_array{
    
    public:
		// monomial \x^{\a} \x\in\Real^n, \a\in\Integer^n_+
        vector<vector<int> > vap;//2 \times vap_size
		// vap[0]...index j of nonzero elements of \a
		// vap[1]...value (\a)_j of \a

        int vap_full_size;// vap_full_size >= vap_size
        int vap_size;	//<--- the amounts of nnz of set of exponents
        
        vector<vector<int> > pnz;// 2 \times pnz_size
		// pnz[0]...the order in all monomials of index j which is the first nonzero elements (\a)_j in \a.
		// pnz[1]...nnz(\a) of \x^{\a}
		
        int pnz_size;	//<--- # of monomials
        int pnz_full_size;// pnz_full_size >= pnz_size
        
        int dim;
        int dim2;
    
        spvec_array();
        spvec_array(int size1,int size2);
       	spvec_array(const spvec_array & sa);
        spvec_array& operator=(const spvec_array& sa);
		~spvec_array();
        
        void alloc(int size1,int size2);
        void del();
        void del_vap();
        void disp();
        void disp(int pos1,int pos2);
        void disp(int* slist);
        void write(string fname);
        
        void clean_pnz();
        void clean_vap();
        void clean();
        
        int get_nnz();
        int deg();
};
class poly_info{
    
    public:
        
        int no;
        
        int typeCone;
        int sizeCone;
        int numMs;
		vector<vector<double> > coef;        

        vector<int> mr;
        vector<int> mc;
        
        class spvec_array sup;
        
        poly_info();
		poly_info(const poly_info &);
        ~poly_info();
        
        void alloc_coef(int typecone,int sizecone,int terms,int nnz);
        
        void del();
        
        void disp();
        
};
class bass_info{
    public:
        
        int  dim;
        int  deg;
        vector<int> pattern;
        
        class spvec_array sup;
        
        bass_info();
        bass_info(const bass_info &);
		~bass_info();
        
        void alloc_pattern(int length);
        
        void del();
        
};


#endif //_spvec_
