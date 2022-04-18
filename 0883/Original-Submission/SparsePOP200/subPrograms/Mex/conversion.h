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

#ifndef _PSDPtoLSDP_
#define _PSDPtoLSDP_

#include "global.h"
using namespace std;

#include "polynomials.h"
#include "spvec.h"
#include "mex.h"
#ifdef mex_h
#define printf mexPrintf
#endif

#define RELAXORDER		1
#define MULTICLIFACT    1 //multiCliquesFactor;
#define EQTOLERANCE     1   //eqTolerance;
#define SPARSESW        1   //sparseSW;
#define COMPSW          1   //complementaritySW;
#define BOUNDSW         1   //boundSW;
#define SCALINGSW       1   //scalingSW;
#define PERTURB         1.0E-5  //perturbation;
#define REDUCEMOMENTSW  1 //reduceMomentMatSW;
#define SEDUMISW    	1
#define SEDUMIONSCREEN  1
#define PRINTONSCREEN	1
#define PRINTLEVEL1		2
#define	PRINTLEVEL2		0
#define	SEDUMIEPSILON	1.0E-9
#define	SYMBOLICMATH	1
#define	MEX				1		


/*** conversion ****************************************/
void conversion_part1(
    /*IN*/  class s3r & sr,
    /*OUT*/ double & objconst,
    int & slen, vector<double> & scaluvalue,
    int & blen, vector<double> & bvect,
    int & mlen, vector<double> & permmatrix);

void conversion_part2(
    /*IN*/  class s3r & sr,
    vector<int>  oriidx,
    class SparseMat extofcsp,
	/*OUT*/ class mysdp & sdpdata);

/*******************************************************/

int write_sdp_polyinfo(string outname, class poly_info & polyinfo);
int write_bassinfo(string outname, class spvec_array & bassinfo);

void gen_maxcliques3(int msize,vector<int> oriidx,class SparseMat extofcsp,class cliques & macls);
void perturb_objective(class poly & objpoly,int dimvar,double eps);
void gen_basisindices(int sparsesw,int multifactor,class polysystem & polysys,class cliques & maxcliques,vector<list<int> > & BasisIndices);

void get_moment_matrix(const class spvec_array & bassinfo,class spvec_array & mm_mat);

void info_a_nnz_a_struct_size_eq( class poly_info & polyinfo, class spvec_array & bassinfo,int & info_size,int & nnz_size,int & struct_size);
void info_a_nnz_a_struct_size_ineq_a_ba1( class poly_info & polyinfo, class spvec_array & bassinfo,int & info_size,int & nnz_size,int & struct_size);
void info_a_nnz_a_struct_size_ba1mmt( class spvec_array & bassinfo,int & info_size,int & nnz_size,int & struct_size);
void info_a_nnz_a_struct_size_ineq_a_ba2( class poly_info & polyinfo, class spvec_array & bassinfo,int & info_size,int & nnz_size,int & struct_size);
void info_a_nnz_a_struct_size_ba2mmt( class spvec_array & bassinfo,int & info_size,int & nnz_size,int & struct_size);
void info_a_nnz_a_struct_size_sdp(int mdim, class poly_info & polyinfo, class spvec_array & bassinfo,int & info_size,int & nnz_size,int & struct_size);
void get_info_a_nnz_a_struct_size(int mdim, int msize, vector<class poly_info> polyinfo, vector<class spvec_array> bassinfo,int & info_size,int & nnz_size,int & blst_size);

void convert_obj(  class poly_infon & polyinfo,class mysdp & psdp);
void convert_eq(  class poly_info & polyinfo,  class spvec_array & bassinfo,class mysdp & psdp);
void convert_ineq_a_ba1(  class poly_info & polyinfo,  class spvec_array & bassinfo,class mysdp & psdp);
void convert_ba1mmt(  class spvec_array & bassinfo,class mysdp & psdp);
void convert_ineq_a_ba2(  class poly_info & polyinfo,  class spvec_array & bassinfo,class mysdp & psdp);
void convert_ba2mmt(  class spvec_array & bassinfo,class mysdp & psdp);
void convert_sdp(class poly_info & polyinfo, class spvec_array & bassinfo, class mysdp & psdp);
void get_psdp(/*IN*/int mdim,int msize, vector<class poly_info> polyinfo,vector<class spvec_array> bassinfo ,/*OUT*/class mysdp & psdp);
void gather_diag_blocks(int gbs,vector<vector<int> > ggg,class mysdp & psdp);

void qsort_psdp(vector<int> & slist,class mysdp & psdp);
void qsort_sups(vector<int> & slist,class spvec_array & supset);
void qsort_normal(/*IN*/vector<int> a, int left, int right,/*OUT*/vector<int> sortedorder);
void swap2(vector<int> a,const int & i,const int & j);

bool comp_sup_a_block(class spvec_array & supset,int i,class sup_a_block & supblock);
bool comp_sup_a_block(class sup_a_block & supblock,/* < */ class spvec_array & supset,int i);
bool comp_sup_a_block(class mysdp & psdp,int i,class sup_a_block & supblock);
bool comp_sup_a_block(class sup_a_block & supblock,class mysdp & psdp,int i);

void variable_numbering(class spvec_array & allsups,vector<int> plist,class mysdp & psdp,vector<int> & degOneTerms);
void count_upper_nnz(vector<int> plist,class mysdp & psdp);

void get_lsdp(class spvec_array & allsupset,class mysdp & psdp,vector<int> & degOneTerms);

void pushsups(/*IN*/ class spvec_array insups,/*OUT*/ class spvec_array & outsups);
void simplification(/*IN*/ class spvec_array & vecs);

void remove_sups(class mysdp & psdp,class spvec_array & removesups);
void remove_sups(class spvec_array & removesups,class spvec_array & sups);

//function that write the sparse format of the SDP into the file
void write_sdpa(/*IN*/class mysdp & psdp,/*OUT*/ string sdpafile);

//return the information of POP
void get_poly_a_bass_info(
   /* IN */  class polysystem & polysys, vector<class supSet> & BaSupVect,vector<class supSet> & mmBaSupVect,
   const int mat_size,
   /* OUT */ vector<class poly_info> & polyinfo,vector<class spvec_array> & bassinfo);

void get_subjectto_polys_and_basups(
   /* IN */  class polysystem polysys, vector<list<int> > BaIndices, vector<class supSet> basups,
   /* OUT */ int stsize,vector<class poly_info> & polyinfo_st, vector<class bass_info> & bassinfo_st);
void get_momentmatrix_basups(class polysystem polysys,vector<list<int> > BaIndices, vector<class supSet> & basups, vector<class bass_info> & bassinfo_st);

void get_allsups(int dim,class poly_info & polyinfo_obj, int stsize, vector<class poly_info> polyinfo_st, vector<class bass_info> bassinfo_st,class spvec_array & allsups);
void get_allsups_in_momentmatrix(int dimvar,int mmsize, vector<class bass_info> bassinfo_mm,class spvec_array & mmsups);

//for sorting functions
void sort_infotable(vector<vector<int> > infotable,vector<int> stand, vector<int> infolist,int left,int right);
bool comp_InfoTable(class Vec3 vec1, class Vec3 vec2);
void sortInfoTable(vector<class Vec3 > &infotable);
void sortInfoTable(vector<class Vec3 > infotable, vector<int> & infolist);
bool comp_infotable(vector<vector<int> > infotable,int i,vector<int> stand);
bool comp_infotable(vector<int> stand, vector<vector<int> > infotable, int j);
bool same_edge(class EdgeSet edge1, class EdgeSet edge2);
bool comp_edge(class EdgeSet edge1, class EdgeSet edge2);

//functions that enumerate all monomials
void count_lexall_num_a_nnz(int dimvar,int deg,int & num,int & nnz);
void genLexFixDeg(int k,int n,int W,vector<vector<int> > sup,int nnz,class spvec_array & rsups);
void genLexAll(int totalOfVars,int Deg,class spvec_array & rsups);

//function that extracts the support of complimentarity constraints.
void get_removesups(/*IN*/class polysystem & polysys,/*OUT*/class spvec_array & removesups);

void initialize_supset(/*IN*/class spvec_array & spvecs,/*OUT*/class supSet & supset);
void initialize_spvecs(/*IN*/class supSet & supset,/*OUT*/class spvec_array & spvecs);
void initialize_spvecs(/*IN*/list<class sup> & suplist,/*OUT*/class spvec_array & spvecs);
void initialize_polyinfo(/*IN*/class polysystem & polysys,int nop,/*OUT*/class poly_info & polyinfo);
void copy_polynomialsdpdata(/*IN*/class mysdp & opsdp,/*OUT*/class mysdp & npspd);

//function that reads the parameter file
bool input_params(class pop_params & params);

class mat_info{
    
    public:
        vector<vector<int> > bij;
		vector<double> coef;
		class spvec_array sup;
        
        void del();
        
        mat_info();
        ~mat_info();
        
};
class mysdp{
    public:
        int mDim;
        int nBlocks;
        vector<int> bLOCKsTruct;
        vector<vector<int> > block_info;	//
        class mat_info ele;
        
        int utsize;
        vector<vector<int> > utnnz;
        int idum1,idum2,idum3;
        
        double dtime ;
        
        mysdp();	//constructor
        ~mysdp();	//destructor
        
        void alloc_structs(int blst_size);
        
        void alloc(int blst_size,int ele_size,int nnz_size);
        void del();
        
        void disp();
        void disp(int b1,int b2);
        void disp_lsdp();
        void disp_sparseformat();
        void write(string fname);
        void write_utnnz(string fname);
        
        void input(string fname);
        
        void set_struct_info(int matstruct,int pos,int nnz_size,int typecone);
        
};

class sup_a_block{
    public:
        int block;
        int deg;
        int nnzsize;
        vector<int> vap0;
        vector<int> vap1;
        
        sup_a_block();
        ~sup_a_block();
        
        void alloc(int mdim);
        void input(class mysdp & psdp,int i);
        void input(class spvec_array & supset,int i);
        void disp();
        
};

class pop_params{
    
    public:
        int relax_Order;	//RelaxOrder
        int multiCliquesFactor; //
        double eqTolerance;	//
        int sparseSW;		//
        int complementaritySW;	//
        int boundSW;		//
        int scalingSW;		//
        double perturbation;	//
        int reduceMomentMatSW;	//
        string detailedInfFile; //name of file to write information of POP
        string sdpaDataFile; //name of file to write SDPA sparse format
        int SeDuMiSW;//
        string SeDuMiOutFile;//
        int SeDuMiOnScreen;//
        int SDPASW;//
       	string SDPAOutFile;
		int SDPAOnScreen;
		string printFileName;	
        int printOnScreen;//
		int printLevel1;
		int printLevel2;
		double SeDuMiEpsilon;
		int symbolicMath;
		int mex; 
		int reduceAMatSW;       
		string Method; 
        pop_params();
        void write_parameters(string fname);
};
class cliques{
    public:
        int numnode;
        int numcliques;
		vector<list<int> > clique;
        cliques();
        ~cliques();
        
        void initialize(int nodes,int csize);
        void write_maxCliques(string fname);
        int maxnnz();
};
//Sum Of Squares and Semidifinite Relaxation
class s3r{
    
    public:
        
        int itemp;
       	string problemName; 
        class polysystem Polysys;
        class polysystem OriPolysys;
        
        s3r();//constructor
      	~s3r(){
			timedata1.clear();
			timedata.clear();
			for(int i=0;i<bindices.size();i++){
				bindices[i].clear();
			}
			degOneTerms.clear();
			problemName.clear();
		};
        class cliques maxcliques;
        
        vector<double> timedata1;
        vector<double> timedata;
        
		vector<double> scalevalue;
		vector<double> bvect;       
		vector<double> permmatrix;
 
        vector<list<int> > bindices;
        
        //parameter
        class pop_params param;
        
        double ctime;
        
        string isFull;
        int linearterms;
       	vector<int> degOneTerms; 
        
        string detailedInfFile;
        string sdpaDataFile;
        
        //objective function
        class poly objPoly;
        
        void set_relaxOrder(int Order=2);
        
        void genBasisSupports(class supsetSet & BasisSupports);
        void reduceSupSets(class supsetSet & BasisSupports,class supSet allNzSups);
        void eraseCompZeroSups(class supSet & czSups,vector<class supSet> & BaSups); //delete the supports from Polynomial SDPs via complementarity constraints
       
        void disp_params();
        
        //*** conversion ***
        
        void write_pop(int i,string fname);
        void write_BasisIndices(string fname);
        void write_BasisSupports(int i, string fname, class supsetSet BasisSupports);
        void redundant_ZeroBounds(class supsetSet BasisSupports, class supSet allSup, class supSet & ZeroSup);
        void redundant_OneBounds(class supsetSet BasisSupports, class supSet allSup, class supSet & OneSup);
        
};

class Vec3{
	public:
		vector<int> vec;
		int no;
		
		Vec3(){};
		~Vec3(){
			vec.clear();
		};
		void clear(){
			vec.clear();
		};
		void resize(int s, int i){
			vec.resize(s,i);
		}

};

class SparseMat{
	// Sparse matrix of MATLAB
	public:
		vector<int> ir;
		vector<int> jc;
		
		SparseMat(){};
		~SparseMat(){
			ir.clear();
			jc.clear();
		};
		void clear(){
			ir.clear();
			jc.clear();
		}
		void resizeIr(int n, int val){
			ir.resize(n,val);
		}
		void resizeJc(int m, int val){
			jc.resize(m,val);
		}
};

class EdgeSet{
	public:
		int vertex1;
		int vertex2;
		EdgeSet(){
			vertex1 = 0;
			vertex2 = 0;	
		};
		~EdgeSet(){};
};

#endif //_gensdp_

