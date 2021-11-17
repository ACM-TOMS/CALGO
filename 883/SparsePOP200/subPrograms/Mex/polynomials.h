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

#ifndef _polynomials_
#define _polynomials_

#include "global.h"
using namespace std;

#include "sup.h"

bool comp(class tmpVec vec1, class tmpVec vec2);
class mono multiMonos(class mono Mono1,class mono Mono2);
class poly multiPolys(class poly Poly1,class poly Poly2);
class poly powPoly(class poly Poly,int k);
class poly assignPolyToVar(class poly Poly, vector<class poly>  polyVec, int nov);

// declaration of class mono
class mono{
    public:
		int nDim;
		vector<int> supIdx;
		vector<int> supVal;
    	vector<double>	Coef;	//coefficient vector
   		
		mono(){//constructor
			nDim = -1;
		};
		/*
		mono(const mono& Mono){//copy constructor
			supIdx.resize(Mono.supIdx.size());
			supVal.resize(Mono.supVal.size());
			Coef.resize(Mono.Coef.size());
			nDim = Mono.nDim;
			copy(Mono.supIdx.begin(),Mono.supIdx.end(),supIdx.begin());
			copy(Mono.supVal.begin(),Mono.supVal.end(),supVal.begin());
			copy(Mono.Coef.begin(),Mono.Coef.end(),Coef.begin());
		};
		*/
		~mono(){//destructor
			supIdx.clear();
			supVal.clear();
			Coef.clear();
		};
		/*
		mono Copy(const mono& Mono){//copy function
			if(this == &Mono){
				return *this;
			}
			supIdx.resize(Mono.supIdx.size());
			supVal.resize(Mono.supVal.size());
			Coef.resize(Mono.Coef.size());
			nDim = Mono.nDim;
			copy(Mono.supIdx.begin(),Mono.supIdx.end(),supIdx.begin());
			copy(Mono.supVal.begin(),Mono.supVal.end(),supVal.begin());
			copy(Mono.Coef.begin(),Mono.Coef.end(),Coef.begin());
			return *this;	
		};
		mono operator=(const mono& Mono){//copy operator
			return Copy(Mono);
		};
		*/
	     
	    //initializiation and allocation
    	void allocSupp(int dimvar);		//for support
	    void allocCoef(int length);		//for coefficient
    	void setSupp(int nov,int value);	//Input support
	    void setCoef(double value,int r=0);	//Input coefficient
    	void resizeCoef(int length,double value);
	    void clear(); 
    	    
	    //display
		void printSupp();	//display support
	    void printCoef();	//display coef
    	void writeMono();	//display mono
        
	    //Return the information of mono
		int	getSupp(int nov);//return the position nov of the support
	    double getCoef(int nof);//return the position nof of the coef
    	int	lengthSupp();	//return the length of the support vector
	    int lengthNzSupp();//return the # of nonzero in the support 
    	int	lengthCoef();//return the length of the coef vector
	    int lengthNzCoef();//return the # of nonzero in the coef
		//return the value of mono for the given numeric vector.
	    double evalMono(vector<double> const Var);
    	void getSuppComb(vector<int> & comb);
	    void copyCoef(vector<double> & coef);
    	void getSup(class sup & Sup);
		void sortMono();
};

//the declaration of class poly.
class poly{
    public:
        
	    int noSys;	//return # of poly in POP (obj. func -> 0, i-th constraint -> i)
    	int typeCone;	//return the kind of obj. func and constraints
	    int sizeCone;	//return the size of coef vectors of poly.
    	int dimVar;	//return the dimension of poly.
	    int degree;	//return the degree of poly.
    	int noTerms;	//return the # of monomials in poly.
	    list<class mono> monoList;// monomials in poly.
        
    	double scaleValue;	//max. abs. val in coefs of poly.
        	
	    poly(){//constructor
    		scaleValue=1.0;
        	noTerms=0;
	       	degree =0;
			noSys = -1;
			typeCone = 1;
			sizeCone = 1;	
		};
		/*	
		poly(const poly& P){//copy constructor
			monoList.resize(P.monoList.size());
			noSys = P.noSys;
			typeCone = P.typeCone;
			sizeCone = P.sizeCone;
			dimVar = P.dimVar;
			degree = P.degree;
			noTerms = P.noTerms;
			scaleValue = P.scaleValue;
			copy(P.monoList.begin(),P.monoList.end(),monoList.begin());	
		};
		*/
	   	~poly(){//destructor
			for(list<mono>::iterator it=monoList.begin();it!=monoList.end();++it){
				(*it).clear();
			}
			monoList.clear();
		};
	    /*	
		poly Copy(const poly& P){//copy function
			if(this == &P){
				return *this;
			}
			monoList.resize(P.monoList.size());
			noSys = P.noSys;
			typeCone = P.typeCone;
			sizeCone = P.sizeCone;
			dimVar = P.dimVar;
			degree = P.degree;
			noTerms = P.noTerms;
			scaleValue = P.scaleValue;
			copy(P.monoList.begin(),P.monoList.end(),monoList.begin());	
			return *this;	
		}

		poly operator=(const poly& Poly){//copy operator
			return Copy(Poly);
		}
		*/
		list<class mono>::iterator monoListBegin(){ return monoList.begin(); }
	    list<class mono>::iterator monoListEnd()  { return monoList.end();   }
        
		//functions for monoList
   		void setNoSys(int nosys);	//Input the # of poly.
	    void setDimVar(int nov);	//Input the dimension of poly.
		//Input the typeCone and the sizeCone of poly.
	    void setTypeSize(int typecone,int sizecone=1);
    	void setDegree();			//Compute the degree of poly.
	    void setDegree(int deg);	//Input the degree of poly.
		//Input the # and dimension of poly.
	    void setNosysDimvar(int nosys,int dimvar);
    	void clear(); 
		//output poly.
   		int dimvar();				//return the dimention of poly.	
	    int noterms();				//return the # of monomials
    	int nosys();				//return NO. of poly.
	    int typecone();				//return the typeCone of poly.	
    	int sizecone();				//return the sizeCone of poly.
	    int Degree();				//return the degree of poly.
    	int lengthCoef();
	    void sumSupports(list<int> & sumSupp);
    	void getSups(class supSet & Sups);
	    void pushSupList(list<class sup> & suplist);
		//Investigate whether poly. is complementarity constraint or not.
		//If yes, store the support and the coef into czSup and Coef, respectively.
    	int isComplementarity(class sup & czSup,double & Coef);
	     
	    //functions for monoList.
    	void addMono(class mono monoNew);//add monoNew into monoList
	    void clearMonolist();			//delete all mono in monoList.
    	void resizeCoef(int length,int nop=0);
       
		//display poly.
		void writePolyData();
        
	    //return the values of poly. for the given numeric vector.
    	void evalPoly(vector<double> const var,vector<double> & xVect,vector<double> & maxAbs);
	    //devides all coefficients by scaleValue
		void scalingPoly(double & ScaleValue, double maxCoef);
	    void getConst(vector<double> & constValue,int isErase=YES);
    	//perturb poly. with eps. 
	    void perturbPoly(int dim,int seed,double eps=1.0E-5);
};	

// the declaration of class bounds.
class bounds{
    private:
        vector<double> upper;
        vector<double> lower;
        public:
			bounds(){};
			/*
			bounds(const bounds & B){
				upper.resize(B.upper.size());
				lower.resize(B.lower.size());
			};
			*/
			~bounds(){
				upper.clear();
				lower.clear();
			};
			/*
			bounds Copy(const bounds& B){
				if(this == &B){
					return *this;
				}
				upper.resize(B.upper.size());
				lower.resize(B.lower.size());
				copy(B.upper.begin(),B.upper.end(),upper.begin());	
				copy(B.lower.begin(),B.lower.end(),lower.begin());
				return *this;	
			};
			bounds operator=(const bounds& B){
				return Copy(B);
			};
			*/
            void allocUp(int no);	//allocate the upper bound.
            void setUp(int novar,double uvalue);	//Input the upper bound of novar-th variable 
            void allocLo(int no);	//allocate the lower bound.
            void allocUpLo(int no);	//allocate the upper and lower bound
            void setLow(int novar,double lvalue);	//Input the lower bound of novar-the variable.
            double lbd(int i);	//return the lower bound of i-th variable
            double ubd(int i);	//return the upper bound of i-th variable
            void printdata();	//display bounds.
			void clear();
};


//declaration of class polysystem
class polysystem{
    
    public:
        int numSys;	//the number of obj. and constraints.	
        int dimVar;	//the number of variable in POP
        vector<class poly> polynomial;
        double objConst;	//constant value of obj. function
        vector<list<int> > posOflbds;// indicates the position of lbd(i) in polysys
        vector<list<int> > posOfubds;// indicates the position of ubd(i) in polysys
        
        vector<int> dumvec;
        int itemp;
        //bounds of POP
        class bounds bounds;
        class bounds boundsNew;
        
        //constructor
        polysystem(){
			numSys = -1;
			dimVar = -1;
			objConst = 0;
		};
		/*
		polysystem(const polysystem& P){
			numSys = P.numSys;
			dimVar = P.dimVar;
			objConst = P.objConst;
			polynomial.resize(P.polynomial.size());
			copy(P.polynomial.begin(),P.polynomial.end(),polynomial.begin());	
			posOflbds.resize(P.posOflbds.size());
			posOfubds.resize(P.posOfubds.size());
			copy(P.posOflbds.begin(),P.posOflbds.end(),posOflbds.begin());	
			copy(P.posOfubds.begin(),P.posOfubds.end(),posOfubds.begin());	
			bounds = P.bounds;
			boundsNew = P.boundsNew;
		};
		*/
        ~polysystem(){
			vector<poly>::iterator itBegin=polynomial.begin();
			vector<poly>::iterator itEnd  =polynomial.end();
			for(vector<poly>::iterator it=itBegin;it!=itEnd;++it){
				(*it).clear();
			}
			polynomial.clear();
			vector<list<int> >::iterator it = posOflbds.begin();
			for(;it!=posOflbds.end();++it){
				(*it).clear();
			}
			posOflbds.clear();
			it = posOfubds.begin();
			for(;it!=posOfubds.end();++it){
				(*it).clear();
			}
			posOfubds.clear();
			dumvec.clear();
			bounds.clear();
			boundsNew.clear();	
		};
       /*
		polysystem Copy(const polysystem& P){
			if(this == &P){
				return *this;
			}
			numSys = P.numSys;
			dimVar = P.dimVar;
			objConst = P.objConst;
			polynomial.resize(P.polynomial.size());
			copy(P.polynomial.begin(),P.polynomial.end(),polynomial.begin());	
			posOflbds.resize(P.posOflbds.size());
			posOfubds.resize(P.posOfubds.size());
			copy(P.posOflbds.begin(),P.posOflbds.end(),posOflbds.begin());	
			copy(P.posOfubds.begin(),P.posOfubds.end(),posOfubds.begin());	
			bounds = P.bounds;
			boundsNew = P.boundsNew;
			return *this;	
		};
		polysystem operator=(const polysystem& P){
			return Copy(P);
		};
		*/	
        //function for obj. function.
		void layawayObjConst();//store the constant value of obj. func.
        double getObjScaleValue();	//return the scalevalue of obj. func.
        void copyObjPoly(class poly & objPoly);	//copy the obj. func.
        void perturbObjPoly(int seed,double eps);
        
        //functions that return the information of POP.
		int numsys();	//the # of obj. and constraints
        int dimvar();	//the dimention of variables in POP
        int numofINE();	//the # of inequality constraints
		//store the union of supports of no-th poly. into sumSupp
        void sumSupports(int no,list<int> & sumSupp);
        int maxDeg();   //return the max degree in POP
        
        void allocSys(int numvar,int numsys);	//initialize POP.
        void writePolynomials();				//display POP.
        void addPoly(class poly Poly);			//add Poly into POP.
        //values of constriants of POP for the given numeric vector
        void evalPolynomials(vector<double> const Var,vector<double>& fValue,vector<double> & maxAbs);
        
        int polyDimvar(int nop);	//return the dimension of nop-th poly.
        int polyNoterms(int nop);	//return the # of monomials of nop-th poly.
        int polyNosys(int nop);		//return No. of poly. in POP
        int polyTypeCone(int nop);	//return the typeCone of nop-th poly.
        int polySizeCone(int nop);	//return the sizeCone of nop-th poly.
        int polyDegree(int nop);    //return the degree of nop-th poly.
		//return the information of Complementarity Constraints.
        int polyIsComplimentarity(int nop,class sup & czSup,double & Coef); 
        int polySupNnz(int nop);	//nop th polynomial form's nonzeros of all supports
        int polyCoefNnz(int nop);	//nop th polynomial form's nonzeros of all coefficinet
        int lengthCoef(int nosys);	//the length of coef. of nosys-th poly in POP
        void polyAddMono(int nop,class mono Mono);//add Mono into nop-th poly. in POP
        void polySetDegree(int nop);//Compute the degree of nop-th poly. in POP.
        void polyWritePolyData(int nop);//display the nop-th poly.
        void setPolyNoSys(int nop,int nosys);	// Input the NO of poly. in POP 
        void setPolyDimVar(int nop,int nov);	//Input the dimension of poly. in POP.
		//Input the typeCone and sizeCone into nop-th poly.
		void setPolyTypeSize(int nop,int typecone,int sizecone=1);
        //return the support of nop-th poly.
        void getPolySups(int nop,class supSet & Sups);
		//Input the degree into nop-th poly.
        void setPolyDegree(int nop, int deg); 
        
        //convert bounds into polys. 
        void boundToIneqPolySys(int typeExchange=1);
		class poly createUpLowPoly(int varNo,double lbd,double ubd);
		class poly createUpPoly(int varNo,double ubd);
		class poly createLowPoly(int varNo,double lbd);
        //Applying scaling technique into POP
        void scalingPOP(vector<double> & pMatrix,vector<double> & bVect,vector<double> & scaleVect);
        void scalingAllPolys(vector<double> & scaleVect, vector<double> vscale);
        void addBoundToPOP_simple(class supSet & allNzSups,int & numAdd);//add the contraints on allNzSups into POP.
        void addBoundToPOP(class supSet ZeroSup, class supSet OneSup,int & numAdd);	
       	void relax1EqTo2Ineqs(double eqTolerance);
        
};

class tmpVec{
	public:
		int idx; 
		int val;
		
		tmpVec(){};
		tmpVec(const tmpVec& t){
			idx = t.idx;
			val = t.val;
		};
		~tmpVec(){};
		tmpVec Copy(const tmpVec& t){
			if(this == &t){
				return *this;
			}
			idx = t.idx;
			val = t.val;
			return *this;
		};
		tmpVec operator=(const tmpVec& t){
			return Copy(t);
		};
};


#endif // _polynomials_



