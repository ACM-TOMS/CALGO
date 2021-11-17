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

#ifndef _sup_
#define _sup_

#include "global.h"
using namespace std;

class sup;
class sup2;
class supSet;
class supSet2;
class supsetSet;

void genLexFixDeg(int k,int n,int W,class sup  Sup,list<class sup> & supList);
void genLexAll(int totalOfVars,int Deg,list<class sup> & supList);
bool comp_sup_for_sort(class sup sup1, class sup sup2);
class sup SupMinusSup(class sup sup1, class sup sup2);
class sup{
    public:
        vector<int> idx;
        vector<int> val;
        
        // these values usually are  not assigned.
        int bij; // for sorting at qsort_psdp in conversion.cpp
        int no; // for sorting at qsort_psdp in conversion.cpp
        
        sup();
        sup(const sup& support);
        ~sup();
        sup& operator=( const sup& support);
        bool operator==(const sup& support) const;
        bool operator!=(const sup& support) const;
        bool operator <(const sup& support) const;
        
        void push(int i,int v);
        void getIdxsVals(vector<int> & Idxs,vector<int> & Vals);
        void changeIndices(vector<int> & o2n_pattern);
        void changeIndices(list<int> o2n_pattern);
        void initByArrayData(const vector<int> & arrayData);
        void clear();
        void erase_end();
        void disp();
        void disp2(int nsize);
        int isEvenSup();
        int dimvar();
        int nnz();
        void assignSupToArray(vector<int> & digit,int type=1);
        int deg();
};
class supSet{
    public:
        
        
        int dimVar;
        int dimVar2;
        
        list<class sup> supList;
        supSet();
        supSet(const supSet& supset);
        ~supSet();
        
        void setSupSet(int v,list<class sup> & suplist);
        void setSupSet(int v,set<class sup> & supset);
        
        void setDimVar(int vnum);
        void changeIndicesAll(vector<int> & o2n_pattern);
        void changeIndicesAll(list<int> o2n_pattern);
        void addSup(class sup & newSup,int type=1);
        void disp();
        int dimvar();
        int deg();
        list<class sup>::iterator begin();
        list<class sup>::iterator end();
        int size();
        void clear();
        //delete specific elments and return the positon of the previous elemens
        list<class sup>::iterator erase(list<class sup>::iterator & supIte);
        bool  doesExist(class sup & Sup);
        void pushSupList(list<class sup> & SupList);
        void pushSupSet(class supSet & newSet);
        void pushSup(class sup & newSup);
        void unique();
        void sort();
        void getEvenSups(class supSet & eSups,int isUnique);
        void out_full(string fName);
        void out_full(int i, string fName);
        int nnz();
};
class sup2{
    private:
	//
	// this sup belong to the l[i]-th position of the r[i]-th BasisSupport.
	//
        vector<int> r;
        vector<int> l;
        
        double ubd;
        double lbd;
        
        vector<int> idx;
        vector<int> val;
        
        public:
            void pushRL(int R,int L);// add the information of r, l
            void pushIdxVal(int Idx,int Val);// add index and degree
            // return the information of index and degree from data of class sup type
            void setSup(class sup Sup);
            void getSup(class sup & Sup);
            void assignSupToArray(vector<int> & Array,int type=1);
            void disp();
            int deg();
            void RL(vector<int> & R,vector<int> & L);
            void setUp(double Ubound);
            void setLow(double Lbound);
            double up();
            double low();
            int compSup(class sup & sup1);
            
};
class supSet2{
    private:
        int dimVar;
        list<class sup2> supList;
        public:
            supSet2(int nVars){
                dimVar=nVars; //give the dimension of variables.
            }
            
            void addSup(int r,int l,class sup & newSup);
            void pushSup(class sup2 & newSup2);
            void disp();
            int dimvar();
            list<class sup2>::iterator begin();
            list<class sup2>::iterator end();
            int size();
            void clear();
};

class supsetSet{
    public:
        vector< class supSet > supsetArray;
        
        void push(class supSet List){
            supsetArray.push_back(List);
        }
        void changeIndicesAll(int nof,vector<int> & o2n_pattern);
        void changeIndicesAll(int nof,list<int> o2n_pattern);
        void getSupSets(vector<class supSet> & supSets);
};

#endif //_sup_
