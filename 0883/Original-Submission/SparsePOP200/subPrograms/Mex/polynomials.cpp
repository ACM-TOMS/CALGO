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

#include "polynomials.h"

void mono::clear(){
	supIdx.clear();
	supVal.clear();
	Coef.clear();
}
void poly::clear(){
	list<mono>::iterator it = monoList.begin();
	for(;it!=monoList.end();++it){
		(*it).clear();
	}
	monoList.clear();
	noTerms = 0;
}
//return the variables that have nonzero exponents.
void mono::getSuppComb(vector<int> & comb){
	comb = supIdx;
}
void mono::resizeCoef(int length,double value){
    Coef.resize(length,value);
}
void mono::writeMono(){
	printf(" Supp:");
	printSupp();
	printf("\n");
	printf(" Coef:");
	printCoef();
	printf("\n\n");
}
//return the length of the support
int mono::lengthSupp(){
	return nDim;
}
//return the length of the coef.
int mono::lengthCoef(){
    return Coef.size();
}
double mono::evalMono(vector<double> const Var){
    double value=1;
	for(int i=0;i<supIdx.size();i++){
		value*= pow(Var[supIdx[i]],supVal[i]);
	}	
    return value;
}
int mono::getSupp(int noterm){
	if(noterm <0 || noterm >= nDim){
		printf("*** Access Error for mono::supVal***\n");
		exit(1);
	}
	vector<int>::iterator vt = supVal.begin();
	for(vector<int>::iterator it=supIdx.begin();it!=supIdx.end();++it){
		if(noterm == (*it)){
			return (*vt);
		}
		++vt;
	}
	return 0;
}
double mono::getCoef(int nof){
    if(Coef.empty()){
        return 0;
    }else if(nof < 0 || nof >= Coef.size()){
        cout << " *** mono::getCoef@polynomials.cpp *** " << endl;
        cout << "nof          = " << nof << endl;
        cout << "size of Coef = " << Coef.size() << endl;
    	exit(1);
	}else{
        return Coef[nof];
    }
}
void mono::printSupp(){
	int k;
	bool flag;
	for(int i=0; i<nDim; i++){
		k = 0;
		flag = true;
		while(k < supIdx.size()){
			if(i == supIdx[k]){
				printf(" %d", supVal[k]);
				flag = false;
				break;
			}else{
				k++;
			}
		}
		if(flag){
			printf(" 0");
		}
	}
}
void mono::printCoef(){
    int length=Coef.size();
    for(int i=0;i<length;i++){
        printf(" %10.8e",Coef[i]);
    }
}
void mono::allocSupp(int dimvar){
	nDim = dimvar;
}
//allocate the length of coef.
void mono::allocCoef(int length){
    Coef.resize(length,0.0);
}
void mono::setSupp(int nov,int value){
	
	if(nov < 0 || nov >= nDim){
        cout << " *** mono::setSupp@polynomials.cpp *** " << endl;
       	cout << "nov   = " << nov << endl;
	    cout << "value = " << value << endl;
		exit(1);
	}else{
		vector<int>::iterator itFirst = supIdx.begin();
		vector<int>::iterator itLast  = supIdx.end();
		vector<int>::iterator vt = supVal.begin();
		vector<int>::iterator iret = find(itFirst,itLast,nov);
		if(value == 0){
			if(iret != itLast){
				int k = distance(itFirst,iret);	
				advance(vt,k);
				supIdx.erase(iret);	
				supVal.erase(vt);	
			}	
		}else if(value > 0){
			if(iret != itLast){
				int k = distance(itFirst,iret);	
				advance(vt,k);
				(*vt) = value; 
			}else{
				supIdx.push_back(nov);
				supVal.push_back(value);
			}	
		}else{
			cout << "*** mono::setSupp@polynomial.cpp" << endl;
			cout << "*** value must be nonnegative." << endl;
			exit(1);
		}
		//this->sortMono();	
	}
}
void mono::setCoef(double value,int r){
    if(r < 0 || r >= Coef.size()){
        cout << " *** mono::setCoef@polynomials.cpp *** " << endl;
        cout << " r     = " << r << endl;
        cout << " value = " << value << endl;
    	exit(1);
	}else{
        Coef[r]=value;
    }
}
/*
poly::~poly(){
	this->clear();
}
*/
void poly::scalingPoly(double & scaleValue, double maxCoef){
    
    //get lenght of coefficient matrix
    int coeflength;
    if(this->typeCone != SDP){
        coeflength = this->sizeCone;
    }
    else{
        coeflength = this->sizeCone*this->sizeCone;
    }
    
    //find max value of all coefficent components
    scaleValue = 1;
    list<mono>::iterator Mono;
    //scaling all coefficient of each monomial by max(abs(coefficients))
    if((maxCoef < 1 && fabs(maxCoef) > EPS) || (1 < maxCoef)){
        scaleValue = maxCoef;
		for(Mono=monoList.begin();Mono!=monoList.end();++Mono){
            for(int j=0;j<coeflength;j++){
                (*Mono).Coef[j] = (*Mono).Coef[j]/scaleValue;
                //printf("Coef = %f\n", (*Mono).Coef[j]);
            }
        }
    }
}
void poly::resizeCoef(int length,int nop){
    
    list<class mono>::iterator Mono=monoList.begin();
    
    for(;Mono!=monoList.end();++Mono){
        (*Mono).resizeCoef(length,(*Mono).getCoef(nop));
    }
}
//clear all monomials of poly.
void poly::clearMonolist(){
    monoList.clear();
    noTerms=0;
}
//addition of monomials into poly.
void poly::addMono(class mono monoNew){
    int i;
    for(i=0;i<monoNew.lengthCoef();i++){
        if(fabs(monoNew.getCoef(i)) > EPS){
            break;
        }
    }
    if(i==monoNew.lengthCoef()){
        return;
    }
    if(noTerms==0){
        monoList.push_back(monoNew);
        noTerms++;
        return;
    }
    else{
        list<class mono>::iterator Mono;
        vector<int>::iterator ite1, ite2, ite3;
        bool flag;
        for(Mono = monoList.begin(); Mono != monoList.end();){
			if((*Mono).supIdx.size() != monoNew.supIdx.size()){
				flag = false;
			}else{
				ite1 = (*Mono).supIdx.begin();
				ite2 = (*Mono).supIdx.end();
				ite3 = monoNew.supIdx.begin();
           		flag = equal(ite1, ite2, ite3);
           		if(flag){
					ite1 = (*Mono).supVal.begin();
					ite2 = (*Mono).supVal.end();
					ite3 = monoNew.supVal.begin();
					flag = equal(ite1, ite2, ite3);
				}
			}
			//focus on coef.
            if(flag){
                for(i=0; i<(*Mono).lengthCoef();i++){
                    (*Mono).Coef[i] += monoNew.Coef[i];
                    if(fabs((*Mono).Coef[i]) > EPS){
                        break;
                    }
                }
                if(i == (*Mono).lengthCoef()){
                    monoList.erase(Mono++);
                    noTerms--;
                }
                return;
            }
			++Mono;
        }
        monoList.push_back(monoNew);
        noTerms++;
        return;
    }
    exit(1);
}
int poly::Degree(){
    setDegree();
    return degree;
}
//return the length of coef.
int poly::lengthCoef(){
    if(typeCone==SDP){
        return sizeCone*sizeCone;
    }
    else{
        return sizeCone;
    }
}

//display poly.
void poly::writePolyData(){
    
    if(noSys==0){
        printf("[objective function[%d]]\n",noSys);
    }else{
        printf("[constraint [%d]]\n",noSys);
    }
    printf("   tpyeCone = %d\n",typeCone);
    printf("   sizeCone = %d\n",sizeCone);
    printf("   dimVar   = %d\n",dimVar);
    printf("   degree   = %d\n",degree);
    printf("   noTerms  = %d\n",noTerms);
    printf("   Monos\n");
	int i = 0;
    list<class mono>::iterator Mono =monoList.begin();
	for(;Mono!=monoList.end();++Mono){
        printf("   mono[%d]:\n",i+1);
        printf("supp:");
        (*Mono).printSupp();
        printf(" : ");
        printf("coef:");
        (*Mono).printCoef();
        printf("\n");
		i++;
    }
}
//evaluate the values of poly. for the specific numerical values.
void poly::evalPoly(vector<double> const Var,vector<double> & fValue,vector<double> & maxAbs){
    
    if(Var.size()!=dimVar){
        cout<<"error:not num o variable == dimVar @evalPoly"<<endl;
        exit(1);
    }
    
    // case of  typeCone = 1, -1
    if(typeCone==EQU || typeCone==INE){
        double mValue;
        double dummy;
        
        fValue.resize(sizeCone,0.0);
        maxAbs.resize(sizeCone,0.0);
        
        list<class mono>::iterator Mono=monoList.begin();
    	for(;Mono!=monoList.end();++Mono){
            mValue=(*Mono).evalMono(Var);
            for(int row=0;row<sizeCone;row++){
                dummy=(*Mono).getCoef(row)*mValue;
                fValue[row]+=dummy;
                dummy=fabs(dummy);
                if(maxAbs[row]<dummy){
                    maxAbs[row]=dummy;
                }
            }
        }
    }
    // case of  typeCone = 2
    else if(typeCone==SOC){
        double mValue,mValue2;
        double norm;
        double maxnorm=0;
        
        fValue.resize(sizeCone,0.0);
        maxAbs.resize(1,0.0);
        
        list<class mono>::iterator Mono=monoList.begin();
        for(;Mono!=monoList.end();++Mono){
            mValue=(*Mono).evalMono(Var);
            norm=0;
            for(int row=0;row<sizeCone;row++){
                mValue2=(*Mono).getCoef(row)*mValue;
                fValue[row]+=mValue2;
                norm+=fabs(mValue2);
            }
            if(maxnorm<norm){
                maxnorm=norm;
            }
        }
        
        maxAbs[0]=maxnorm;
        
    }
    ///case of typeCone = 3
    else if(typeCone==SDP){
        
        double norm;
        double maxnorm=0;
        double  mValue,mValue2;
        int row;
        int length=sizeCone*sizeCone;
        
        fValue.resize(length,0.0);
        maxAbs.resize(1,0.0);
        
        list<class mono>::iterator Mono=monoList.begin();
        for(;Mono!=monoList.end();++Mono){
            mValue=(*Mono).evalMono(Var);
            norm=0;
            for(row=0;row<length;row++){
                mValue2=(*Mono).getCoef(row)*mValue;
                fValue[row]+=mValue2;
                norm+=fabs(mValue2);
            }
            if(maxnorm<norm){
                maxnorm=norm;
            }
        }
        maxAbs[0]=maxnorm;
    }
    else {
        cout<< "error@evalPolynomial: the typeCone is not available, typeCone=" <<typeCone<<endl;
        exit(1);
    }
    
}
//Input the NO. of eq and the # of variables.
void poly::setNosysDimvar(int nosys,int dimvar){
    noSys=nosys;
    dimVar=dimvar;
}
//Input the # of constraints.
void poly::setNoSys(int nosys){ noSys=nosys;}
//Input the dimension of variables
void poly::setDimVar(int nov){ dimVar=nov;}
void poly::setTypeSize(int typecone,int sizecone){ typeCone=typecone; sizeCone=sizecone; }
//Input the degree of poly.
void poly::setDegree(){
    int dummy=0;
    int max=dummy;
    list<class mono>::iterator Mono= monoList.begin();
    for(;Mono!=monoList.end();++Mono){
		dummy = accumulate((*Mono).supVal.begin(),(*Mono).supVal.end(),0);
        if(max<dummy){
            max=dummy;
        }
    }
    degree=max;
}
void poly::setDegree(int deg){
    degree=deg;
}
//Substitute the specfic value into one specific variable
//arg: 1st: poly
//     2nd: NO. of variable.
class poly assignPolyToVar(class poly Poly, vector<class poly> polyVec, int nov){
    // We assume that Poly does not have multiple monomials.
    // In mexconv1 and mexconv2, multiple monomials are already removed from Poly.
    if(0 <= nov && nov<Poly.dimVar){

        //Create poly.
        class poly polyNew;
        class mono subMono, aMono;
        bool flag;

        polyNew = Poly;
        list<mono>::iterator Mono,it;
        int val;
        int	min = Poly.sizeCone;
        if(min > polyVec[0].sizeCone){
            min = polyVec[0].sizeCone;
        }
        //double s1,s2;
        //double s = 0;
		list<mono> expandMono;
		list<int> pow;
		for(Mono=polyNew.monoList.begin();Mono != polyNew.monoList.end();){
			val = (*Mono).getSupp(nov);
			if(val>0){
				pow.push_back(val);
				expandMono.push_back((*Mono));
				polyNew.monoList.erase(Mono++);//remove Mono from polyNew.monoList
				polyNew.noTerms--;
			}else{
                 ++Mono;
            }
		}
		list<int>::iterator pt = pow.begin();
		for(Mono=expandMono.begin();Mono != expandMono.end(); ++Mono){
			//
			// Elements in expanedMono have pow-th positive element.
			// Because setSupp overwrites pow-th element, we do not need
			// to apply sortMono into them.
			//
			subMono = (*Mono);
			subMono.setSupp(nov,0);
			if(subMono.supIdx.empty() == false){
				for(it = polyVec[(*pt)-1].monoList.begin(); it != polyVec[(*pt)-1].monoList.end(); ++it){
					aMono = (*Mono);
					val = (*it).getSupp(nov); 
					aMono.setSupp(nov,val);
					for(int i=0;i<min;i++){
						aMono.setCoef((*it).Coef[i]*subMono.Coef[i],i);
					}
					polyNew.addMono(aMono);
				}
			}else{	
				for(it = polyVec[(*pt)-1].monoList.begin(); it != polyVec[(*pt)-1].monoList.end(); ++it){
					aMono = (*it);
					for(int i=0;i<min;i++){
						aMono.setCoef((*it).Coef[i]*subMono.Coef[i],i);
					}
					polyNew.addMono(aMono);
				}
			}
			++pt;
		}
		//printf("addMono time = %f\n",s);
		return polyNew;
	}else{
		cout<<"error@assignPolytoVar: the NumverOfVariable is not included in the Polynomial"<<endl;
		exit(1);
	}
}
int poly::dimvar() { return dimVar;}
int poly::noterms(){ return monoList.size();}
int poly::nosys(){ return noSys;}
int poly::typecone(){ return typeCone; }
int poly::sizecone(){ return sizeCone; }
/*
//constructor
polysystem::polysystem(){
    this->objConst=0;
}
polysystem::~polysystem(){
    polynomial.clear();
    dumvec.clear();
	for(int i=0; i< posOflbds.size(); i++){	
    	posOflbds[i].clear();// indicates the position of lbd(i) in polysys
    	posOfubds[i].clear();// indicates the position of ubd(i) in polysys
	}
	posOflbds.clear();
	posOfubds.clear();
}
*/

void polysystem::scalingAllPolys(vector<double> & scaleVect, vector<double> vscale){
    
    scaleVect.resize(this->numsys(),1.0);
    
    this->itemp = 999;
    
    double ScaleValue=1.0;
    double maxCoef;
    double ratio;
	vector<int>::iterator it,vt;
    list<mono>::iterator ite;
    for(int i=0; i<this->numSys; i++){
        maxCoef = 0;
        scaleVect[i] = 1.0;
        if(polynomial[i].typeCone == 1 || polynomial[i].typeCone == -1){
    		ite = polynomial[i].monoList.begin();
            for(;ite!=polynomial[i].monoList.end();++ite){
                ratio = 1;
				vt = (*ite).supVal.begin();
				for(it = (*ite).supIdx.begin();it != (*ite).supIdx.end();++it){
					ratio *= pow(vscale[(*it)],(*vt));
					++vt;
				}
                for(int k=0; k < polynomial[i].sizeCone;k++){
                    (*ite).Coef[k] = ratio*(*ite).Coef[k];
                    //cout << "maxCoef = " << maxCoef << endl;
                    //cout << "Coef = " << fabs(sup.Coef[k]) << endl;
                    if(maxCoef < fabs((*ite).Coef[k])){
                        maxCoef = fabs((*ite).Coef[k]);
                    }
                }
            }
        }else{
			maxCoef = 1;
		}
        polynomial[i].scalingPoly(scaleVect[i], maxCoef);
    }
}
void bounds::allocUp(int no){
    upper.resize(no,10*MAX);
}
void bounds::setUp(int novar,double uvalue){
    upper[novar-1]=uvalue;
}
void bounds::allocLo(int no){
    lower.resize(no,-10*MAX);
}
void bounds::setLow(int novar,double lvalue){
    lower[novar-1]=lvalue;
}
double bounds::lbd(int i){
    if(lower.empty()){
        return 0;
    }else if( i < 0 || i >= lower.size() ){
        cout << "out of range of lower. " << endl;
        cout << " i = " << i << endl;
        cout << " lower of size = " << lower.size() << endl;
        exit(1);
    }else{
        return lower[i];
    }
}
double bounds::ubd(int i){
    if(upper.empty()){
        return 0;
    }else if( i < 0 || i>= upper.size() ){
        cout << "out of range of upper. " << endl;
        cout << " i = " << i << endl;
        cout << " upper of size = " << upper.size() << endl;
        exit(1);
    }else{
        return upper[i];
    }
}
void bounds::printdata(){
    
    cout<<"[lower bounds for variables]"<<endl;
    for(int i=0;i<lower.size();i++ ){
        cout<<"   x"<<i+1<<">="<<lower[i]<<endl;
    }
    
    cout<<"[upper bounds for variables]"<<endl;
    for(int i=0;i<upper.size();i++ ){
        cout<<"   x"<<i+1<<"<="<<upper[i]<<endl;
    }
    
}
void bounds::allocUpLo(int no){
    this->allocUp(no);
    this->allocLo(no);
}
void bounds::clear(){
    upper.clear();
    lower.clear();
}

void polysystem::scalingPOP(vector<double> & pMatrix,vector<double> & bVect,vector<double> & scaleVect){
    
    //cout<<" ***>polysystem::scalingPOP ---> "<<endl;
    
    //applying the scaling technique via affine trasformation for variables
    //allocate the new upper and lower bounds.
    //boundsNew.allocUpLo(dimVar);
    
    //pMatrix.resize(this->dimVar*this->dimVar,0.0);
    pMatrix.resize(this->dimVar,1);
    bVect.resize(this->dimVar,0.0);
    for(int i=0;i<this->dimVar;i++){
        if(MIN<bounds.lbd(i)&&bounds.ubd(i)<MAX){
            //pMatrix[i+i*this->dimVar] = bounds.ubd(i)-bounds.lbd(i);
            pMatrix[i] = bounds.ubd(i)-bounds.lbd(i);
            bVect[i] = bounds.lbd(i);
        }else if(bounds.lbd(i) <= MIN && bounds.ubd(i) < MAX){
            //pMatrix[i+i*this->dimVar] = -1;
            pMatrix[i] = -1;
            bVect[i]=bounds.ubd(i);
        }else if(bounds.lbd(i) > MIN && bounds.ubd(i) >= MAX){
            //pMatrix[i+i*this->dimVar] = 1;
            pMatrix[i] = 1;
            if(fabs(bounds.lbd(i)) >= EPS){
                bVect[i]=bounds.lbd(i);
            }
        }else{
            //pMatrix[i+i*this->dimVar]=1.0;
            pMatrix[i]=1.0;
        }
	//cout << " bVect[" << i << "] = " << bVect[i] << endl;
    }
    
    //cout << "***Start of maing VarInConst and VarInConstList" << endl;
    //system("top -b -n 1 | grep MATLAB | head -1 |awk '{printf(\"memory = %s\\n\"), $6}' ");
    vector<list<int> > VarInConstList(dimVar);
    vector<int> MaxDegVar(dimVar);
    MaxDegVar.resize(dimVar,0);
	vector<int>::iterator it,vt;
	for(int i=0; i<polynomial.size();i++){
        for(list<mono>::iterator mit = polynomial[i].monoList.begin(); mit!=polynomial[i].monoList.end(); ++mit){
			vt = (*mit).supVal.begin();
			for(it =(*mit).supIdx.begin();it!=(*mit).supIdx.end(); ++it){
               	VarInConstList[(*it)].push_back(i);
               	MaxDegVar[(*it)] = ((*vt)>MaxDegVar[(*it)]) ? (*vt) : MaxDegVar[(*it)];
				++vt;
           	}
        }
    }
	for(int i=0; i<dimVar; i++){
		VarInConstList[i].sort();
		VarInConstList[i].unique();
	}
	/*
    cout << "***Start of affine scaling" << endl;
	system("top -b -n 1 | grep MATLAB | head -1 |awk '{printf(\"memory = %s\\n\"), $6}' ");
	*/
   	int NumOfActiveBounds = 0; 
    list<int>::iterator iteList;
    //double s = 0;
	//int v = 0;
	//double s1 = (double)clock();
    for(int i=0;i<dimVar;i++){
        //cout << " bounds lbd = " << bounds.lbd(i)<< endl;
        //cout << " bounds ubd = " << bounds.ubd(i)<< endl;
        //case:both upper bound and lower bound exist.
        if(MIN<bounds.lbd(i) && bounds.ubd(i)<MAX){
            if(fabs(bounds.lbd(i)) >= EPS || fabs(bounds.ubd(i)-1) >= EPS){
                //Substitute the values into i-th variables of all constraints.
				if(MaxDegVar[i] > 0){
                    //Substitute the values into i-th variables of all constraints.
                    double t0 = (double)clock();
                    if(VarInConstList[i].empty() != true){
                        //Make poly. for affine trasformation.
                        class poly Poly=createUpLowPoly(i,bounds.lbd(i),bounds.ubd(i));
                        vector<poly> polyVec(MaxDegVar[i]);
                        polyVec[0] = Poly;
                        for(int k=1;k<MaxDegVar[i];k++){
                            polyVec[k] = multiPolys(Poly,polyVec[k-1]);
                        }
                        for(iteList = VarInConstList[i].begin();iteList != VarInConstList[i].end(); ++iteList){
                            //s1 = (double)clock();
                            polynomial[(*iteList)] = assignPolyToVar(polynomial[(*iteList)],polyVec,i);
                            //s2 = (double)clock();
                            //s = s + (s2 - s1)/(double)CLOCKS_PER_SEC;
                            //v++;
                        }
                    }
                }
            }
            //store the new bounds.
            boundsNew.setUp(i+1,1.0);
            boundsNew.setLow(i+1,0.0);
   			NumOfActiveBounds += 2; 
        }
        //POP has only lower bounds
        else if(MIN < bounds.lbd(i)){
            if(fabs(bounds.lbd(i)) >= EPS && MaxDegVar[i] > 0){
                //Substitute the values into i-th variables of all constraints.
                if(VarInConstList[i].empty() != true){
                    //Make poly. for affine trasformation.
                    class poly Poly=createLowPoly(i,bounds.lbd(i));
                    vector<poly> polyVec(MaxDegVar[i]);
                    polyVec[0] = Poly;
                    for(int k=1;k<MaxDegVar[i];k++){
                        polyVec[k] = multiPolys(Poly,polyVec[k-1]);
					}
                    for(iteList = VarInConstList[i].begin();iteList != VarInConstList[i].end(); ++iteList){
                        polynomial[(*iteList)] = assignPolyToVar(polynomial[(*iteList)],polyVec,i);
                    	//polynomial[(*iteList)].writePolyData();
					}
                }
            }
            boundsNew.setUp(i+1,bounds.ubd(i));
            boundsNew.setLow(i+1,0);
   			NumOfActiveBounds++; 
        //cout << "boundsNew.lbd( " << i << ") = " << boundsNew.lbd(i) << endl;
        //cout << "boundsNew.ubd( " << i << ") = " << boundsNew.ubd(i) << endl;
        }
        //POP has only upper bounds
        else if(bounds.ubd(i)<MAX){
            //Substitute the values into i-th variables of all constraints.
            if(MaxDegVar[i] > 0){
				if(VarInConstList[i].empty() != true){
    	            //Make poly. for affine trasformation.
        	        class poly Poly=createUpPoly(i,bounds.ubd(i));
            	    vector<poly> polyVec(MaxDegVar[i]);
                	polyVec[0] = Poly;
	                for(int k=1;k<MaxDegVar[i];k++){
    	                polyVec[k] = multiPolys(Poly,polyVec[k-1]);
        	        }
            	    for(iteList = VarInConstList[i].begin();iteList != VarInConstList[i].end(); ++iteList){
                    	polynomial[(*iteList)] = assignPolyToVar(polynomial[(*iteList)],polyVec,i);
                	}
            	}
			}
            //Substitute the values into i-th variables of all constraints.
            boundsNew.setUp(i+1,-bounds.lbd(i));
            boundsNew.setLow(i+1,0);
   			NumOfActiveBounds++; 
        }
        //cout << "boundsNew.ubd( " << i << ") = " << boundsNew.ubd(i) << endl;
        //cout << "boundsNew.lbd( " << i << ") = " << boundsNew.lbd(i) << endl;
    }
    //double s2 = (double)clock();
	//system("top -b -n 1 | grep MATLAB | head -1 |awk '{printf(\"memory = %s\\n\"), $6}' ");
    //cout << "***End of affine scaling" << endl;
    //printf("# of assign = %d\n",v);
    //printf("Time assign = %f\n",(s2-s1)/(double)CLOCKS_PER_SEC);
    //printf("ave. assign = %f\n",(s2-s1)/((double)CLOCKS_PER_SEC*v));
	for(int i=0;i<dimVar;i++){
		VarInConstList[i].clear();
	}
	MaxDegVar.clear();
	VarInConstList.clear();
   	//writePolynomials(); 
    //cout << "***Start of addition of bounds" << endl;
    //system("top -b -n 1 | grep MATLAB | head -1 |awk '{printf(\"memory = %s\\n\"), $6}' ");
    this->itemp = 888;
	vector<poly> tmpPolys(numSys);
	for(int i=0; i< numSys; i++){
		tmpPolys[i] = polynomial[i];
	}	
	polynomial.resize(numSys+NumOfActiveBounds);
	for(int i=0; i< numSys; i++){
		polynomial[i] = tmpPolys[i];
		tmpPolys[i].clear();
	}	
    tmpPolys.clear();
	this->boundToIneqPolySys();
	//system("top -b -n 1 | grep MATLAB | head -1 |awk '{printf(\"memory = %s\\n\"), $6}' ");
    //cout << "***End of addition of bounds" << endl;
    //Extract the constant value of obj. poly.
    this->layawayObjConst();
    
    //cout << "***Start of scaling bounds" << endl;
    //system("top -b -n 1 | grep MATLAB | head -1 |awk '{printf(\"memory = %s\\n\"), $6}' ");
    // The following is the same as scaling.m/
    double absBound = 0;
    vector<double> vscale(dimVar);
    for(int j=0; j<dimVar ; j++){
        vscale[j] = 1;
        if(MIN < boundsNew.lbd(j) && boundsNew.ubd(j) < MAX){
            if(fabs(boundsNew.lbd(j)) < fabs(boundsNew.ubd(j))){
                absBound = fabs(boundsNew.ubd(j));
            }else{
                absBound = fabs(boundsNew.lbd(j));
            }
            if(absBound > EPS){
                if(1 < absBound){
                    vscale[j] = absBound;
                    boundsNew.setUp(j+1,boundsNew.ubd(j)/vscale[j]);
                    boundsNew.setLow(j+1,boundsNew.lbd(j)/vscale[j]);
                }else if(absBound < 1){
                    vscale[j] = absBound;
                    boundsNew.setUp(j+1,boundsNew.ubd(j)/vscale[j]);
                    boundsNew.setLow(j+1,boundsNew.lbd(j)/vscale[j]);
                }
            }
            //pMatrix[j+j*this->dimVar]=pMatrix[j+j*dimVar]*vscale[j];
            pMatrix[j] *= vscale[j];
        }
    }
    //system("top -b -n 1 | grep MATLAB | head -1 |awk '{printf(\"memory = %s\\n\"), $6}' ");
    //cout << "***End of scaling bounds" << endl;
    
    //Scaling all constraints via. absolute values for each poly.
    this->scalingAllPolys(scaleVect,vscale);
	//writePolynomials();
    //cout<<" <--- polysystem::scalingPOP <*** "<<endl<<endl;
}

int polysystem::lengthCoef(int nosys){
    
    int typecone=this->polyTypeCone(nosys);
    int sizecone=this->polySizeCone(nosys);
    
    if( typecone == SDP){
        return sizecone*sizecone;
    }
    else{
        return sizecone;
    }
    
}
void polysystem::setPolyNoSys(int nop,int nosys){
    polynomial[nop].setNoSys(nosys);
}
void polysystem::setPolyDimVar(int nop,int nov){
    polynomial[nop].setDimVar(nov);
}
void polysystem::setPolyTypeSize(int nop,int typecone,int sizecone){
    polynomial[nop].setTypeSize(typecone,sizecone);
}
void polysystem::polyWritePolyData(int nop){
    polynomial[nop].writePolyData();
}
void polysystem::polySetDegree(int nop){
    polynomial[nop].setDegree();
}
void polysystem::polyAddMono(int nop,class mono Mono){
    polynomial[nop].addMono(Mono);
}

int polysystem::polyDimvar(int nop){
    return polynomial[nop].dimvar();
}
int polysystem::polyDegree(int nop){
    return polynomial[nop].Degree();
}
int polysystem::polyNoterms(int nop){
    return polynomial[nop].noterms();
}
int polysystem::polyNosys(int nop){
    return polynomial[nop].nosys();
}
int polysystem::polyTypeCone(int nop){
    return polynomial[nop].typecone();
}
int polysystem::polySizeCone(int nop){
    return polynomial[nop].sizecone();
}

void polysystem::evalPolynomials(vector<double> const Var,vector<double> & fValue,vector<double> & maxAbs){
    
    for(int i=0;i<this->numsys();i++){
        vector<double> fDummy;
        vector<double> maDummy;
        
        polynomial[i].evalPoly(Var,fDummy,maDummy);
        
        for(int j=0;j<fDummy.size();j++){;
        	fValue.push_back(fDummy[j]);
        }
        for(int j=0;j<maDummy.size();j++){
            maxAbs.push_back(maDummy[j]);
        }
    }
}
void polysystem::addPoly(class poly Poly){
    numSys++;
    polynomial.push_back(Poly);
}

int polysystem::numsys(){
    return this->polynomial.size();
}
void polysystem::allocSys(int numvar,int numsys){
    numSys=numsys;
    dimVar=numvar;
    polynomial.resize(numsys);
	bounds.allocLo(numvar);
	bounds.allocUp(numvar);
	boundsNew.allocLo(numvar);
	boundsNew.allocUp(numvar);
}
void polysystem:: writePolynomials(){
    cout <<"[[ Polynomial System ]]"<<endl;
    
    cout<<endl;
    cout<<"   NumOfPolynomials:"<<numSys<<endl;
    cout<<"   NumOfVariables:"<<dimVar<<endl<<endl;
    
    //display all poly.s
    for(int i=0;i<this->numsys();i++){
        polynomial[i].writePolyData();
    }
    //display upper and lower bounds.
    bounds.printdata();
    //boundsNew.printdata();
}
int polysystem::dimvar(){
    return this->dimVar;
}
void mono::copyCoef(vector<double> & coef){
    coef=this->Coef;
}

void poly::sumSupports(list<int> & varList){
	varList.clear();	
	list<mono>::iterator ite;
	vector<int>::iterator vt;
	for(ite = monoList.begin(); ite!=monoList.end(); ++ite){
		for(vt=(*ite).supIdx.begin();vt!=(*ite).supIdx.end();++vt){
			varList.push_back((*vt));
		}
	}
	varList.sort();
	varList.unique();
	/*
		cout << "==="<<endl;
		for(vit=varList.begin();vit!=varList.end();++vit){
			cout << "(*vit) = " << (*vit) <<endl;
		}	
		cout << endl;
	*/
}
int polysystem::numofINE(){
    int num=0;
    for(int i=0;i<this->numSys;i++){
        if( this->polynomial[i].typecone()==INE){
            num++;
        }
    }
    return num;
}
void polysystem::sumSupports(int no,list<int> & varList){
    this->polynomial[no].sumSupports(varList);
}

void polysystem::boundToIneqPolySys(int typeExchange){
    
    //cout<<" ***>polysystem::boundToIneqPolySys ---> "<<endl;
	//cout << "length of vector polynomial = " << polynomial.size() << endl; 
   	//cout << "the number of numSys        = " << numSys << endl; 
    //add lower bounds into class poly
    for(int i=0;i<dimVar;i++){
        //cout << "boundsNew.lbd(" << i << ")=" << bounds.lbd(i) << endl;
        //cout << "MIN                    =" << MIN << endl;
        if(boundsNew.lbd(i)>MIN){
            polynomial[numSys].setNoSys(numSys);
            polynomial[numSys].setDimVar(dimVar);
            polynomial[numSys].setTypeSize(INE,1);
            int length=polynomial[numSys].lengthCoef();
            
            if( fabs(boundsNew.lbd(i)) > EPS){
    			class mono Mono1,Mono2;
                Mono1.allocSupp(dimVar); 
				Mono1.allocCoef(length);
                Mono2.allocSupp(dimVar); 
				Mono2.allocCoef(length);
                Mono2.setSupp(i,1);
                Mono1.setCoef(-boundsNew.lbd(i),0);
                Mono2.setCoef( 1.0,0);
				polynomial[numSys].monoList.push_back(Mono1);
				polynomial[numSys].monoList.push_back(Mono2);
				polynomial[numSys].noTerms = 2;
                polynomial[numSys].setDegree(1);
            }else{
               	class mono Mono1;
                Mono1.allocSupp(dimVar); 
				Mono1.allocCoef(length);
                Mono1.setSupp(i,1);
                Mono1.setCoef( 1,0);
				polynomial[numSys].monoList.push_back(Mono1);
				polynomial[numSys].noTerms = 1;
                polynomial[numSys].setDegree(1);
            }
            //add poly. into POP
            numSys++;
			posOflbds[i].push_back(numSys);
        }
    }
    //add upper bounds into class poly
    for(int i=0;i<dimVar;i++){
        if(boundsNew.ubd(i)<MAX){
            //cout << "boundsNew.ubd(" << i << ")=" << boundsNew.ubd(i) << endl;
            //cout << "MAX              =" << MAX << endl;
            polynomial[numSys].setNoSys(numSys);
            polynomial[numSys].setDimVar(dimVar);
            polynomial[numSys].setTypeSize(INE,1);
            
            int length=polynomial[numSys].lengthCoef();
            class mono Mono1,Mono2;
            Mono1.allocSupp(dimVar); 
			Mono1.allocCoef(length);
            Mono2.allocSupp(dimVar); 
			Mono2.allocCoef(length);
            Mono2.setSupp(i,1);
            Mono1.setCoef( boundsNew.ubd(i),0);
            Mono2.setCoef(-1.0,0);
			if(fabs(boundsNew.ubd(i)) > EPS){
				polynomial[numSys].monoList.push_back(Mono1);
				polynomial[numSys].monoList.push_back(Mono2);
				polynomial[numSys].noTerms = 2;
			}else{
				//printf("boudsNew.ubd(%4d) = %10.8e\n",i,boundsNew.ubd(i));
				polynomial[numSys].monoList.push_back(Mono2);
				polynomial[numSys].noTerms = 1;
			}
            polynomial[numSys].setDegree(1);
            //add poly. into POP
            numSys++;
            posOfubds[i].push_back(numSys);
        }
    }
	//cout << "numSym = " << numSys << endl;
    //cout<<" <--- polysystem::boundToIneqPolysys <*** "<<endl<<endl;
}
//making poly for affine trasformtion.
class poly polysystem::createLowPoly(int varNo,double lbd){
    
    class poly newPoly;
    newPoly.setDimVar(dimVar);
    newPoly.setTypeSize(1,1); //typeCone -1,or,1
    
    class mono Mono1;
    Mono1.allocSupp(dimVar);
    Mono1.allocCoef(1);
    Mono1.setCoef(lbd,0);
    
    class mono Mono2;
    Mono2.allocSupp(dimVar);
    Mono2.allocCoef(1);
    Mono2.setSupp(varNo,1);
    Mono2.setCoef(1,0);
   	
	if(fabs(lbd) >EPS){
		newPoly.monoList.push_back(Mono1);
		newPoly.monoList.push_back(Mono2);
		newPoly.noTerms = 2;
	}else{
		newPoly.monoList.push_back(Mono2);
		newPoly.noTerms = 1;
	}
	newPoly.setDegree(1);	
    return (newPoly);
}
//making poly for affine trasformtion.
class poly polysystem::createUpPoly(int varNo,double ubd){
    
    class poly newPoly;
    newPoly.setDimVar(dimVar);
    newPoly.setTypeSize(1,1); //typeCone??-1,or,1
    
    class mono Mono1;
    Mono1.allocSupp(dimVar);
    Mono1.allocCoef(1);
    Mono1.setCoef(ubd,0);
    
    class mono Mono2;
    Mono2.allocSupp(dimVar);
    Mono2.allocCoef(1);
    Mono2.setSupp(varNo,1);
    Mono2.setCoef(-1,0);
   
	if(fabs(ubd) > EPS){
		newPoly.monoList.push_back(Mono1);
		newPoly.monoList.push_back(Mono2);
		newPoly.noTerms = 2;
	}else{
		newPoly.monoList.push_back(Mono2);
		newPoly.noTerms = 1;
	}
	newPoly.setDegree(1);	
    return (newPoly);
}

//making poly for affine trasformtion.
class poly polysystem::createUpLowPoly(int varNo,double lbd,double ubd){
   
 
    if(lbd<ubd){
        
        class poly newPoly;
        newPoly.setDimVar(dimVar);
        newPoly.setTypeSize(1,1); //typeCone??-1,or,1
        int length=newPoly.lengthCoef();
        
        class mono Mono1;
        Mono1.allocSupp(dimVar);
        Mono1.allocCoef(length);
        Mono1.setCoef(lbd);
        
        class mono Mono2;
        Mono2.allocSupp(dimVar);
        Mono2.allocCoef(length);
        Mono2.setSupp(varNo,1);
        Mono2.setCoef(ubd-lbd);
        
		if(fabs(ubd-lbd) > EPS && fabs(lbd) > EPS){
			newPoly.monoList.push_back(Mono1);
			newPoly.monoList.push_back(Mono2);
			newPoly.noTerms = 2;
		}else if(fabs(ubd-lbd) > EPS){
			newPoly.monoList.push_back(Mono2);
			newPoly.noTerms = 1;
		}
	newPoly.setDegree(1);	
        return (newPoly);
    }
    else if(lbd==ubd){
        
        class poly newPoly;
        
        newPoly.setDimVar(dimVar);
        newPoly.setTypeSize(1,1);		//typeCone??-1,or,1
        
        int length=newPoly.lengthCoef();
        
        class mono Mono1;
        Mono1.allocSupp(dimVar);
        Mono1.allocCoef(length);
        Mono1.setCoef(lbd);
        
		if(fabs(lbd) > EPS){
        	newPoly.monoList.push_back(Mono1);
			newPoly.noTerms = 1;
			newPoly.setDegree(1);	
        }
		return (newPoly);
    }
    else{
	cout << "varNo = " << varNo << endl;
	cout << "lbd   = " << lbd << endl;	
	cout << "ubd   = " << ubd << endl;	
        cout<<"error@createUPLOWPoly: not available ubd,lbd: lbd>ubd";
        exit(1);
    }
}

int polysystem::maxDeg(){
    
    int fsize=polynomial.size();
    int maxdeg=0;
    int deg=0;
    
    for(int i=0;i<fsize;i++){
        deg=polynomial[i].Degree();
        if(maxdeg<deg){
            maxdeg=deg;
        }
    }
    
    return deg;
    
}
void mono::getSup(class sup & Sup){
	Sup.idx = supIdx;
	Sup.val = supVal;
}

void poly::getSups(class supSet & Sups){
    
    list<mono>::iterator Mono = monoList.begin();
    for(;Mono!=monoList.end(); ++Mono){
        class sup Sup;
        (*Mono).getSup(Sup);
        Sups.addSup(Sup);
    }
    
}
void poly::pushSupList(list<class sup> & suplist){
    
    list<mono>::iterator Mono = monoList.begin();
  	for(;Mono!=monoList.end();++Mono){  
        class sup Sup;
        (*Mono).getSup(Sup);
        suplist.push_back(Sup);
    }
    
}
void polysystem::getPolySups(int nop,class supSet & Sups){
    if(nop < 0 || nop >= numSys){
        cout << "number of poly must be in [0," << numSys << ")." << endl;
        return;
    }else{
        Sups.setDimVar(this->dimVar);
        this->polynomial[nop].getSups(Sups);
        return;
    }
}

int mono::lengthNzSupp(){
	return supIdx.size();
    
}
//get nonozeros of coefficient
//[input ] nothing
//[output](return) nonzeros of coefficient
int mono::lengthNzCoef(){
    
    int nzSize=0;
    int size=this->Coef.size();
    for(int i=0;i<size;i++){
        if(this->Coef[i]!=0){
            nzSize++;
        }
    }
    return nzSize;
}
int polysystem::polyIsComplimentarity(int nop,class sup & czSup,double & Coef)
{
    if(nop<this->numsys()){
        if(this->polynomial[nop].isComplementarity(czSup,Coef)==YES){
            return YES;
        }
    }
    else{
        cout<<"error@polysystem::polyIsComplimentarity :: Nop'th polynomial doesn't exit"<<endl;
        exit(1);
    }
    return NO;
}

int poly::isComplementarity(class sup & czSup,double & Coef){
   	int SW = 1; 
    if(typeCone==EQU && sizeCone == 1){
        if(monoList.size() ==2){
            class mono mono1, mono2;
            list<class mono>::iterator ite = monoList.begin();
            mono1 = (*ite);
            int nzSize1 = mono1.lengthNzSupp();
            ++ite;
            mono2 = (*ite);
            int nzSize2 = mono2.lengthNzSupp();
            if(nzSize1==0 && nzSize2>=2){
                Coef=mono1.getCoef(0);
                //cout << " Coef = " << Coef << endl;
                //mono2.writeMono();
                if(fabs(Coef) < EPS){
                    if(SW == 1){
						czSup.idx = mono2.supIdx;
						czSup.val = mono2.supVal;
					}
					//mono2.getSup(czSup);
                    //czSup.disp();
                    return YES;
                }else{
                    return NO;
                }
                
            }else if(nzSize1>=2 && nzSize2==0){
                Coef=mono2.getCoef(0);
                //mono1.writeMono();
                if(fabs(Coef) < EPS){
                    if(SW == 1){
						czSup.idx = mono1.supIdx;
						czSup.val = mono1.supVal;
					}
					//mono1.getSup(czSup);
                    //czSup.disp();
					return YES;
                }else{
                    return NO;
                }
            }
            
        }else if( this->monoList.size()==1){
            class mono mono1;
            list<class mono>::iterator ite = this->monoList.begin();
            mono1 = (*ite);
            //mono1.writeMono();
            int nzSize1= mono1.lengthNzSupp();
            if(nzSize1>=2){
                //mono1.getSup(czSup);
				czSup.idx = mono1.supIdx;
				czSup.val = mono1.supVal;
                Coef=0;
				//czSup.disp();
                return YES;
            }
        }
    }
    
    return NO;
}

void polysystem::addBoundToPOP(class supSet ZeroSup, class supSet OneSup,int & numAdd){
    int dum = 0;
    int zsize = ZeroSup.size();
    int osize = OneSup.size();
    //cout << "zsize = " << zsize << endl;
    //cout << "osize = " << osize << endl;
    list<sup>::iterator ite;
    for(ite = ZeroSup.supList.begin();ite != ZeroSup.supList.end();++ite){
        if((*ite).deg()>=2){
            
            class mono Mono3;
            Mono3.allocSupp(this->dimVar);
            Mono3.allocCoef(1);
            Mono3.setCoef(1.0);
           	Mono3.supIdx = (*ite).idx;
			Mono3.supVal = (*ite).val; 
            class poly lwPoly;
            lwPoly.setDimVar(this->dimVar);
            lwPoly.setTypeSize(INE);
			lwPoly.monoList.push_back(Mono3);
			lwPoly.noTerms = 1;
            lwPoly.setDegree();
            dum++;
            lwPoly.setNoSys(numSys+dum);
            
            this->polynomial.push_back(lwPoly);
        }
    }
    for(ite = OneSup.supList.begin();ite != OneSup.supList.end();ite++){
        if((*ite).deg()>=2){
            
            class mono Mono1,Mono2;
            Mono1.allocSupp(this->dimVar);	Mono2.allocSupp(this->dimVar);
            Mono1.allocCoef(1);	Mono2.allocCoef(1);
			Mono1.supIdx = (*ite).idx;
			Mono1.supVal = (*ite).val;
            Mono1.setCoef(-1.0);
            Mono2.setCoef( 1.0);
            
            class poly uPoly;
            uPoly.setDimVar(this->dimVar);
            uPoly.setTypeSize(INE);
			uPoly.monoList.push_back(Mono1);
			uPoly.monoList.push_back(Mono2);
			uPoly.noTerms = 2;
            uPoly.setDegree();
            dum++;
            uPoly.setNoSys(numSys+dum);
            this->polynomial.push_back(uPoly);
        }
    }
    numAdd = dum;
}

void polysystem::addBoundToPOP_simple(class supSet & allNzSups,int & numAdd){
    vector<int> NZidx;
    vector<int> NZpow;
    int nDim = this->dimVar;
    int i,j,k, nnz;
    bool neglbd = false;
    bool finiteubd = false;
    bool zerlbd = false;
    bool Infubd;
    int currentNumSys = numSys;
    double yUpperBound, yLowerBound, absBd, yBound;
    list<class sup>::iterator ite;
    /*
    for(i=0; i < nDim;i++){
        //cout << "lbd(" << i << ")=" << this->bounds.lbd(i) << endl;
        //cout << "ubd(" << i << ")=" << this->bounds.ubd(i) << endl;
        cout << "Newlbd(" << i << ")=" << this->boundsNew.lbd(i) << endl;
        cout << "Newubd(" << i << ")=" << this->boundsNew.ubd(i) << endl;
    }
     */
	int s=0;
    int p=0;
    for(ite = allNzSups.supList.begin();ite != allNzSups.supList.end();++ite){
        if((*ite).deg()>=2){
            s = 0;
            neglbd = false;
            finiteubd = false;
            zerlbd = true;
            NZidx.clear();
            NZpow.clear();
            nnz = (*ite).nnz();
            NZidx.resize(nnz,0);
            NZpow.resize(nnz,0);
            NZidx = (*ite).idx;
            NZpow = (*ite).val;
            for(i=0;i<(*ite).nnz();i++){
                if(this->boundsNew.lbd(NZidx[i]) < -EPS){
                    neglbd = true;
                }
                if(this->boundsNew.ubd(NZidx[i]) >=  MAX){
                    finiteubd = true;
                }
                if(fabs(this->boundsNew.lbd(NZidx[i])) >= EPS){
                    zerlbd = false;
                }
            }
            if(!neglbd && !finiteubd){
                if(zerlbd){
                    yLowerBound = 0;
                }else{
                    yLowerBound = 1;
                    for(i = 0; i < (*ite).nnz(); i++){
                        yLowerBound = yLowerBound*pow(this->boundsNew.lbd(NZidx[i]), NZpow[i]);
                    }
                }
                yUpperBound = 1;
                for(i = 0; i < (*ite).nnz(); i++){
                    yUpperBound = yUpperBound*pow(this->boundsNew.ubd(NZidx[i]), NZpow[i]);
                }
                class mono Mono1,Mono2;
                Mono1.allocSupp(this->dimVar);	Mono2.allocSupp(this->dimVar);
                Mono1.allocCoef(1);	Mono2.allocCoef(1);
				Mono1.supIdx = (*ite).idx;
				Mono1.supVal = (*ite).val;
                Mono1.setCoef(1.0,0);
                Mono2.setCoef(-yLowerBound,0);
                
                class poly lPoly;
                lPoly.setDimVar(nDim);
                lPoly.setNoSys(currentNumSys+1);
                lPoly.setTypeSize(INE);
				lPoly.monoList.push_back(Mono1);
				lPoly.noTerms = 1;
                if(fabs(yLowerBound) >= EPS){
					lPoly.monoList.push_back(Mono2);
					lPoly.noTerms = 2;
                }
                lPoly.setDegree();
                this->polynomial.push_back(lPoly);
                currentNumSys++;
                
                Mono1.allocSupp(this->dimVar);	Mono2.allocSupp(this->dimVar);
                Mono1.allocCoef(1);	Mono2.allocCoef(1);
                Mono1.setCoef(-1.0,0);
                Mono2.setCoef(yUpperBound,0);
				Mono1.supIdx = (*ite).idx;                
				Mono1.supVal = (*ite).val;                


                class poly uPoly;
                uPoly.setNoSys(currentNumSys+1);
                uPoly.setDimVar(nDim);
                uPoly.setTypeSize(INE);
				uPoly.monoList.push_back(Mono1);
				uPoly.monoList.push_back(Mono2);
				uPoly.noTerms = 2;
                uPoly.setDegree();
                this->polynomial.push_back(uPoly);
                currentNumSys++;
                s = 2;
                numAdd = numAdd + 2;
            }else if(!neglbd){
                if(zerlbd){
                    yLowerBound = 0;
                }else{
                    yLowerBound = 1;
                    for(i = 0; i < (*ite).nnz(); i++){
                        yLowerBound = yLowerBound*pow(this->boundsNew.lbd(NZidx[i]), NZpow[i]);
                    }
                }
                class mono Mono1,Mono2;
                Mono1.allocSupp(this->dimVar);	Mono2.allocSupp(this->dimVar);
                Mono1.allocCoef(1);	Mono2.allocCoef(1);
                Mono1.setCoef(1.0,0);
				Mono1.supIdx = (*ite).idx;
				Mono1.supVal = (*ite).val;
                Mono2.setCoef(-yLowerBound,0);
                
                class poly lPoly;
                lPoly.setNoSys(currentNumSys+1);
                lPoly.setDimVar(nDim);
                lPoly.setTypeSize(INE);
				lPoly.monoList.push_back(Mono1);
				lPoly.monoList.push_back(Mono2);
				lPoly.noTerms =2;
                lPoly.setDegree();
                this->polynomial.push_back(lPoly);
                currentNumSys++;
                numAdd = numAdd + 1;
                s = 1;
            }else{
                absBd = 0;
                for(i = 0; i < (*ite).nnz(); i++){
                    if(absBd < fabs(this->boundsNew.lbd(NZidx[i]))){
                        absBd = fabs(this->boundsNew.lbd(NZidx[i]));
                    }
                    if(absBd < fabs(this->boundsNew.ubd(NZidx[i]))){
                        absBd = fabs(this->boundsNew.ubd(NZidx[i]));
                    }
                }
                Infubd = 0;
                if(absBd > MAX){
                    Infubd = 1;
                }
                j = (*ite).isEvenSup();
                if(Infubd == 0 && j == YES){
                    yBound = 1;
                    double tempCoef = 0;
                    for(i=0; i<(*ite).nnz(); i++){
                        if(tempCoef < fabs(boundsNew.lbd(NZidx[i]))){
                            tempCoef = fabs(boundsNew.lbd(NZidx[i]));
                        }
                        if(tempCoef < fabs(boundsNew.ubd(NZidx[i]))){
                            tempCoef = fabs(boundsNew.ubd(NZidx[i]));
                        }
                        yBound = yBound*pow(tempCoef,NZpow[i]);
                        tempCoef = 0;
                    }
                    class mono Mono1, Mono2;
                    Mono1.allocSupp(this->dimVar);	Mono2.allocSupp(this->dimVar);
                    Mono1.allocCoef(1);	Mono2.allocCoef(1);
                    Mono1.setCoef(-1.0,0);
					Mono1.supIdx = (*ite).idx;
					Mono1.supVal = (*ite).val;
                    Mono2.setCoef(yBound,0);
                    
                    class poly uPoly;
                    uPoly.setNoSys(currentNumSys+1);
                    uPoly.setDimVar(nDim);
                    uPoly.setTypeSize(INE);
					uPoly.monoList.push_back(Mono1);
					uPoly.monoList.push_back(Mono2);
					uPoly.noTerms = 2;
                    uPoly.setDegree();
                    this->polynomial.push_back(uPoly);
                    currentNumSys++;
                    
                    Mono1.allocSupp(this->dimVar);
                    Mono1.allocCoef(1);
					Mono1.supIdx = (*ite).idx;
					Mono1.supVal = (*ite).val;
                    Mono1.setCoef(1.0,0);
                    
                    class poly lPoly;
                    lPoly.setNoSys(currentNumSys+1);
                    lPoly.setDimVar(nDim);
                    lPoly.setTypeSize(INE);
                    lPoly.monoList.push_back(Mono1);
					lPoly.noTerms = 1;
                    lPoly.setDegree();
                    this->polynomial.push_back(lPoly);
                    currentNumSys++;
                    
                    numAdd = numAdd + 2;
                    s = 2;
                    /*
                    cout << " neglbd = " << neglbd << endl;
                    cout << " finiteubd = " << finiteubd << endl;
                    cout << " zerlbd = " << zerlbd << endl;
                     */
                }else if(Infubd == 0 && j == NO){
                    yBound = 1;
                    double tempCoef = 0;
                    for(i=0; i<(*ite).nnz(); i++){
                        if(tempCoef < fabs(boundsNew.lbd(NZidx[i]))){
                            tempCoef = fabs(boundsNew.lbd(NZidx[i]));
                        }
                        if(tempCoef < fabs(boundsNew.ubd(NZidx[i]))){
                            tempCoef = fabs(boundsNew.ubd(NZidx[i]));
                        }
                        //cout << "tempCoef = " << tempCoef << endl;
                        yBound = yBound*pow(tempCoef,NZpow[i]);
                        tempCoef = 0;
                    }
                    //cout << "yBound = " << yBound << endl;
                    class mono Mono1,Mono2;
                    Mono1.allocSupp(this->dimVar);	Mono2.allocSupp(this->dimVar);
                    Mono1.allocCoef(1);	Mono2.allocCoef(1);
                    Mono1.setCoef(-1.0,0);
					Mono1.supIdx = (*ite).idx;
					Mono1.supVal = (*ite).val;
                    Mono2.setCoef(yBound,0);
                    
                    class poly uPoly;
                    uPoly.setNoSys(currentNumSys+1);
                    uPoly.setDimVar(nDim);
                    uPoly.setTypeSize(INE);
					uPoly.monoList.push_back(Mono2);
					uPoly.monoList.push_back(Mono1);
					uPoly.noTerms = 2;
                    uPoly.setDegree();
                    this->polynomial.push_back(uPoly);
                    currentNumSys++;
                    //uPoly.writePolyData();
                    
                    Mono1.allocSupp(this->dimVar);	Mono2.allocSupp(this->dimVar);
                    Mono1.allocCoef(1);	Mono2.allocCoef(1);
                    Mono1.setCoef(1.0,0);
					Mono1.supIdx = (*ite).idx;
					Mono1.supVal = (*ite).val;
                    Mono2.setCoef(yBound,0);
                    
                    class poly lPoly;
                    lPoly.setNoSys(currentNumSys+1);
                    lPoly.setDimVar(nDim);
                    lPoly.setTypeSize(INE);
					lPoly.monoList.push_back(Mono1);
					lPoly.monoList.push_back(Mono2);
					lPoly.noTerms = 2;
                    lPoly.setDegree();
                    this->polynomial.push_back(lPoly);
                    currentNumSys++;
                    
                    s = 2;
                    numAdd = numAdd + 2;
                }
            }
        }else{
            s = 0;
            /*
            cout << " neglbd = " << neglbd << endl;
            cout << " finiteubd = " << finiteubd << endl;
            cout << " zerlbd = " << zerlbd << endl;
             */
        }
    }
    //cout << " size = " << s << endl;
    //cout << " size0 = " << p << endl;
    //cout << "numAdd = " << numAdd << endl;
}
void polysystem::relax1EqTo2Ineqs(double eqTolerance){
    if(fabs(eqTolerance) > EPS){
        int noeq = 0;
        int nDim = dimVar;
        for(int i=0; i<numSys; i++){
            if(polynomial[i].typeCone == -1){
                noeq++;
            }
        }
        class poly Newpoly;
        class mono Zero;
        int pointer = numSys;
        for(int i=1; i<numSys; i++){
            if(polynomial[i].typeCone == -1){
                
                //cout << " i = " << i << endl;
                pointer++;
                polynomial[i].setTypeSize(1,polynomial[i].sizeCone);
                Newpoly.setNosysDimvar(pointer,nDim);
                Newpoly.setTypeSize(1,polynomial[i].sizeCone);
                Newpoly.setDegree(polynomial[i].degree);
                Newpoly.monoList = polynomial[i].monoList;
                
                int doesHaveZero = polynomial[i].noTerms;
                //polynomial[i].writePolyData();
                list<mono>::iterator ite  = polynomial[i].monoList.begin();
                int sum;
                for(int j=0;j <polynomial[i].noTerms;j++){
					sum = accumulate((*ite).supVal.begin(),(*ite).supVal.end(),0);
                    if(sum == 0){
                        doesHaveZero = j;
                        break;
                    }
                }
                list<mono>::iterator ite2 = Newpoly.monoList.begin();
                if(doesHaveZero == polynomial[i].noTerms){
                    for(int ell = 0;  ell < polynomial[i].noTerms; ell++){
                        for(int k=0; k<polynomial[i].sizeCone; k++){
                            (*ite2).Coef[k] = -(*ite).Coef[k];
                        }
                        ++ite;
                        ++ite2;
                    }
                    Zero.allocCoef(polynomial[i].sizeCone);
                    for(int ell=0; ell<polynomial[i].sizeCone; ell++){
                        Zero.Coef[ell] = eqTolerance;
                    }
                    polynomial[i].addMono(Zero);
                    polynomial[i].noTerms = 1+polynomial[i].noTerms;
                    Newpoly.addMono(Zero);
                    //Newpoly.noTerms = 1+polynomial[i].noTerms;
                }else{
                    for(int ell = 0;  ell < polynomial[i].noTerms; ell++){
                        if(ell == doesHaveZero){
                            for(int k=0; k<polynomial[i].sizeCone; k++){
                                (*ite2).Coef[k] = -(*ite).Coef[k]+eqTolerance;
                                (*ite).Coef[k] = (*ite).Coef[k] + eqTolerance;
                            }
                        }else{
                            for(int k=0; k<polynomial[i].sizeCone; k++){
                                (*ite2).Coef[k] = -(*ite).Coef[k];
                            }
                        }
                        ++ite;
                        ++ite2;
                    }
                }
                polynomial.push_back(Newpoly);
            }
        }
        numSys = numSys+noeq;
    }
    return;
}

void poly::getConst(vector<double> & constValue,int isErase){
    
    constValue.resize(sizeCone,0);
    list<class mono>::iterator mmIte = monoList.begin();
    list<class mono>::iterator ite;
    int SW = 1;
	for(int i=0; i < noTerms; i++){
        bool flag = true;
		
		if((*mmIte).supIdx.empty()){
			flag = true; 
		}else{
			flag = false;
		}
        //cout << " flag = " << flag << endl;
        if(flag){
            //If isErase=YES, delete the constant value from obj.
            //cout << " size of Coef = " << (*mmIte).Coef.size() << endl;
            //cout << " size of vec  = " << constValue.size() << endl;
            for(int k=0; k < (*mmIte).Coef.size(); k++){
                constValue[k] = constValue[k] + (*mmIte).Coef[k];
            }
            if(isErase == YES){
                ite = mmIte;
                ++mmIte;
                this->monoList.erase(ite);
                this->noTerms -- ;
            }else{
                ++mmIte;
            }
        }else{
            ++mmIte;
        }
        //cout << " i = " << i << endl;
    }
}
void polysystem::layawayObjConst(){
    vector<double> constValue;
    this->polynomial[0].getConst(constValue,YES);
    if(constValue.size()>=1){
        this->objConst += constValue[0];
    }
}
double polysystem::getObjScaleValue(){
    return this->polynomial[0].scaleValue;
}
void polysystem::copyObjPoly(class poly & objPoly){
    objPoly = this->polynomial[0];
}
//perturb obj. poly.
void poly::perturbPoly(int dim,int seed,double eps){
    
    srand(seed);
    
    for(int i=0;i<dim;i++){
        class mono Mono;
        Mono.allocCoef(1);
        Mono.allocSupp(dim);
        Mono.setCoef(eps* (2*(rand()/(double)RAND_MAX)-1) );
        Mono.setSupp(i,1);
        //cout << eps*(2*(rand()/(double) RAND_MAX)-1) << endl;
        this->addMono(Mono);
        
    }
    
}
void polysystem::perturbObjPoly(int seed,double eps){
    
    this->polynomial[0].perturbPoly(this->dimVar,seed,eps);
    
}
//production of mono1 and mono2
class mono multiMonos(class mono Mono1,class mono Mono2){
    
    
    int slength;
    int clength;
    
    if(Mono1.lengthCoef()==Mono2.lengthCoef()){
        clength=Mono1.lengthCoef();
    }
    else{
        cout << "error@multiMonos: Length of Coef vectors is not same"<<endl;
        exit(1);
    }
    
    if(Mono1.lengthSupp()==Mono2.lengthSupp()){
        slength=Mono1.lengthSupp();
    }
    else{
        cout << "error@multiMonos: Length of Support vectors is not same"<<endl;
        exit(1);
    }
    
    class mono Mono3;
    
    Mono3.allocSupp(slength);
    Mono3.allocCoef(clength);
    
   	Mono3.nDim = Mono1.nDim; 
    //addition of supports
    bool flag;
	vector<int>::iterator it1 = Mono1.supIdx.begin();
	vector<int>::iterator it1E= Mono1.supIdx.end();
	vector<int>::iterator it2, vt2;
	vector<int>::iterator it2E= Mono2.supIdx.end();
	vector<int>::iterator vt1 = Mono1.supVal.begin();
	if(Mono1.supIdx.empty() == false){
		while(it1 != it1E){
			flag = false;
			vt2 = Mono2.supVal.begin();
			it2 = Mono2.supIdx.begin();
			while(it2 != it2E){
				if((*it1) == (*it2)){
					flag = true;
					break;
				}else{
					++it2;
					++vt2;
				}	
			}
			if(flag){
				Mono3.supIdx.push_back((*it1));	
				Mono3.supVal.push_back((*vt1)+(*vt2));	
			}else{
				Mono3.supIdx.push_back((*it1));	
				Mono3.supVal.push_back((*vt1));	
			}
			++it1;
			++vt1;
		}
		it2 = Mono2.supIdx.begin();
		vt2 = Mono2.supVal.begin();
		while(it2 != it2E){
			flag = false;
			it1 = Mono1.supIdx.begin();
			while(it1 != it1E){
				if((*it2) == (*it1)){
					flag = true;
					break;
				}else{
					++it1;
				}	
			}
			if(!flag){
				Mono3.supIdx.push_back((*it2));	
				Mono3.supVal.push_back((*vt2));	
			}
			++it2;
			++vt2;
		}
	}else{
		Mono3 = Mono2;
	}

	/* sort */
	if(Mono3.supIdx.empty() == false){
		Mono3.sortMono();
	}

    //production of coef.s
    for(int i=0;i<clength;i++){
        Mono3.setCoef(Mono1.getCoef(i)*Mono2.getCoef(i),i);
    }
	/*
	cout << "Mono1 " ;
	Mono1.printSupp();
	cout << endl;
	cout << "Mono2 " ;
	Mono2.printSupp();
	cout << endl;
	cout << "Mono3 " ;
	Mono3.printSupp();
	cout << endl;
    */
	return Mono3;
    
}

void mono::sortMono(){
	int kDim = supIdx.size();
	vector<class tmpVec> vec1(kDim);
	for(int i=0; i<kDim;i++){
		vec1[i].idx = supIdx[i];
		vec1[i].val = supVal[i];
	}
	vector<tmpVec>::iterator itB = vec1.begin();	
	vector<tmpVec>::iterator itE = vec1.end();	
	sort(itB,itE,comp);
	for(int i=0; i<kDim;i++){
		supIdx[i] = vec1[i].idx;
		supVal[i] = vec1[i].val;
	}
	/*	
	vector<int>::iterator vt1 = supVal.begin();
	cout << "sort===" << endl;	
	for(vector<int>::iterator it1=supIdx.begin();it1!=supIdx.end();++it1){
		cout << "(*it1) = " << (*it1) << ", (*vt1) = " << (*vt1) << endl;
		++vt1;
	}
	*/
	vec1.clear();
}

bool comp(class tmpVec vec1, class tmpVec vec2){

	if(vec1.idx <= vec2.idx){
		return true;
	}else{	
		return false;
	}
}


// production of poly.s
class poly multiPolys(class poly Poly1,class poly Poly2){
    
    if(Poly1.dimvar()!=Poly2.dimvar()){
        cout << "error@sumPoly: not same Dimension of Variables of Poly1 and Poly2"<<endl;
        exit(1);
    }
    if(Poly1.lengthCoef()!=Poly2.lengthCoef()){
        cout<< "error@sumPoly: not same Length of Coeffient Vector of Poly1 and Poly2"<<endl;
        exit(1);
    }
    
    class poly polyNew;
    polyNew.setDimVar(Poly1.dimvar());
    polyNew.setTypeSize(Poly1.typecone(),Poly1.sizecone());
    
    list<mono>::iterator Mono1 = Poly1.monoList.begin();
    list<mono>::iterator Mono2;
	/*
  	cout << "*** Poly1 ***" << endl; 
	Poly1.writePolyData();
  	cout << "*** Poly2 ***" << endl; 
	Poly2.writePolyData();
 	*/
    for(;Mono1!=Poly1.monoList.end();++Mono1){
        for(Mono2 = Poly2.monoList.begin();Mono2!=Poly2.monoList.end();++Mono2){
            class mono aMono = multiMonos((*Mono1),(*Mono2));
           	/*
			cout << "Mono1==="<<endl;
			(*Mono1).writeMono();
           	cout << "Mono2==="<<endl;
			(*Mono2).writeMono();
			*/ 
			polyNew.addMono( aMono );
        	/*
			cout << "== aMono == " << endl;
			aMono.writeMono();
        	cout << "======== == " << endl;
			*/
		}
    }
    polyNew.setDegree();
   	//polyNew.writePolyData(); 
    return polyNew;
}
//k-th power of poly.
//k must be nonnegative integer.
class poly powPoly(class poly Poly,int k){
    
    if(k==0){
        class poly polyNew;
        class mono Mono;
        Mono.allocCoef(Poly.lengthCoef());
        for(int i=0;i<Poly.lengthCoef();i++){
            Mono.setCoef(1.0,i);
        }
        Mono.allocSupp(Poly.dimvar());
        polyNew.addMono(Mono);
        polyNew.setDimVar(Poly.dimvar());
        polyNew.setDegree();
        return polyNew;
    }
    else if(k>=1){
        
        if(k==1){
            return Poly;
        }
        else{
            class poly polyNew=Poly;
            for(int i=1;i<k;i++){
                polyNew=multiPolys(polyNew,Poly);
            }
            return polyNew;
        }
    }
    else{
        cout<<"error@powPoly:this functon cannot permit you to caluculate Pow(*,k) <- k>=1"<<endl;
        exit(1);
    }
}
void polysystem::setPolyDegree(int nop, int deg){
    this->polynomial[nop].setDegree(deg);
}
//nop th polynomial form's nonzeros of all supports
//[input1  ](int nop) the number of the polynomial form
//[output1 ](return)  the nonzeros of all supports
int polysystem::polySupNnz(int nop){
    int polynnz = 0;
    list<class mono>::iterator ite = polynomial[nop].monoList.begin();
    for(;ite!=polynomial[nop].monoList.end(); ++ite){    
		polynnz += (*ite).supIdx.size();
    }
    return polynnz;
}
//nop th polynomial form's nonzeros of all coefficinet
//[input1  ](int nop) the number of the polynomial form
//[output1 ](return)  the nonzeros of all coefficinet
int polysystem::polyCoefNnz(int nop){
    //nonzeros of all coefficient
    int coefnnz = 0;
    
    list<class mono>::iterator ite = polynomial[nop].monoList.begin();
   	for(;ite!=polynomial[nop].monoList.end();++ite){ 
	    coefnnz += (*ite).Coef.size();
    }
    
    return coefnnz;
    
}

