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

#include "spvec.h"

spvec_array::spvec_array(){
    pnz_full_size = 0;
    vap_full_size = 0;
    pnz_size = 0;
    vap_size = 0;
	alloc(0,0);
}
spvec_array::spvec_array(const spvec_array& sa){
	pnz_full_size = sa.pnz_full_size;
	vap_full_size = sa.vap_full_size;
	alloc(pnz_full_size,vap_full_size);
	pnz_size      = sa.pnz_size;
	vap_size      = sa.vap_size;
	for(int i=0; i<sa.pnz[0].size();i++){
		pnz[0][i] = sa.pnz[0][i];
		pnz[1][i] = sa.pnz[1][i];
	}
	for(int i=0; i<sa.vap[0].size();i++){
		vap[0][i] = sa.vap[0][i];
		vap[1][i] = sa.vap[1][i];
	}
	dim           = sa.dim;
	dim2          = sa.dim2;
}
spvec_array& spvec_array::operator=(const spvec_array& sa){
	if(this == &sa){
		return *this;
	}
	alloc(sa.pnz_full_size,sa.vap_full_size);
	pnz_size      = sa.pnz_size;
	vap_size      = sa.vap_size;
	for(int i=0; i<sa.pnz[0].size();i++){
		pnz[0][i] = sa.pnz[0][i];
		pnz[1][i] = sa.pnz[1][i];
	}
	for(int i=0; i<sa.vap[0].size();i++){
		vap[0][i] = sa.vap[0][i];
		vap[1][i] = sa.vap[1][i];
	}
	return *this;
}
spvec_array::spvec_array(int size1,int size2){
    pnz_size = 0;
    vap_size = 0;
    alloc(size1,size2);
}
spvec_array::~spvec_array(){
    del();
}
void spvec_array::alloc(int size1,int size2){
    pnz_full_size = size1;
   	vap_full_size = size2;
	if(size1 >= 0){
		pnz.resize(2);
		pnz[0].resize(size1,0);
		pnz[1].resize(size1,0);
	}
	if(size2 >= 0){
		vap.resize(2);
		vap[0].resize(size2,0);
		vap[1].resize(size2,0);
	}
	pnz_size = 0;
	vap_size = 0;
}

void spvec_array::del(){
	if(!this->vap.empty()){
		if(!this->vap[0].empty()){
			this->vap[0].clear();
		}
		if(!this->vap[1].empty()){
			this->vap[1].clear();
		}
		this->vap.clear();
	}
	if(!this->pnz.empty()){
		if(!this->pnz[0].empty()){
			this->pnz[0].clear();
		}
		if(!this->pnz[1].empty()){
			this->pnz[1].clear();
		}
		this->pnz.clear();
	}
	if(!vap.empty()||!pnz.empty()){
		cout << "Not deleted " << endl;
	}
}
void spvec_array::disp(int pos1,int pos2){
    
    int i,j;
    
    if(pos2 > this->pnz_size ){
        cout<<"pos2 is limit over. so change ->"<<this->pnz_size<<endl;
        pos2 = this->pnz_size;
    }
    
    if(pos1 <0 || pos2 < 0){
        cout<<" pos1 and pos2 is negative(<0)"<<endl;
        exit(1);
    }
    
    cout<<" SIZE = "<<pos2 - pos1 <<endl;
    for(i=pos1; i<pos2; i++){
        cout<<"[ "<<i<<" ]";
        if(this->pnz[0][i] >=0){
            cout<<" |var|";
            for(j=0;j<this->pnz[1][i];j++){
                cout<<" "<<this->vap[0][ this->pnz[0][i] + j ]+1;
            }
            cout<<" |pow|";
            for(j=0;j<this->pnz[1][i];j++){
                cout<<" "<<this->vap[1][this->pnz[0][i]+j];
            }
        }
        cout<<endl;
    }
    
}

void spvec_array::disp(){
    
    int i,j;
    cout<<" SIZE = "<<this->pnz_size <<endl;
    for(i=0; i<this->pnz_size; i++){
        cout<<"[ "<<i<<" ]";
        if(this->pnz[0][i] >=0){
            cout<<" |var|";
            for(j=0;j<this->pnz[1][i];j++){
                cout<<" "<<this->vap[0][ this->pnz[0][i] + j ]+1;
            }
            cout<<" |pow|";
            for(j=0;j<this->pnz[1][i];j++){
                cout<<" "<<this->vap[1][this->pnz[0][i]+j];
            }
        }
        cout<<endl;
    }
    
    
}

poly_info::poly_info(){
	no = -9999;
	typeCone= -999;
	sizeCone = 0;;
	numMs = 0;
	coef.resize(0);        
	mr.resize(0);
	mc.resize(0);
	sup.alloc(0,0);
}
poly_info::poly_info(const poly_info &polyinfo){
	no       = polyinfo.no;
	typeCone = polyinfo.typeCone;
	sizeCone = polyinfo.sizeCone;
	numMs    = polyinfo.numMs;
	coef     = polyinfo.coef;
	mr       = polyinfo.mr;
	mc       = polyinfo.mc;
	sup      = polyinfo.sup;
}
void poly_info::alloc_coef(int typecone,int sizecone,int terms,int nnz){
    
    if(typecone != SDP){
        //this->mr.resize(nnz,0);
        //this->mc.resize(sizecone * terms + 1);
		coef.resize(terms);
		for(int i=0;i<terms;i++){
			coef[i].resize(sizecone,0);
		}
	}
    else{ //SDP case
        mr.resize(nnz,0);
        mc.resize(sizecone * terms + 1);
		coef.resize(nnz);
        for(int i = 0;i<nnz;i++){
			coef[i].resize(sizecone*sizecone,0);
        }
		
	}
    
}
void poly_info::del(){
	if(!mr.empty()){
    	mr.clear();
	}
	if(!mc.empty()){
    	mc.clear();
	}
	if(!coef.empty()){
		for(int i = 0; i < coef.size(); i++){
			if(!coef[i].empty()){
				coef[i].clear();
			}
		}
		coef.clear();
	}
	sup.del();
}
poly_info::~poly_info(){
    del();
}
int spvec_array::get_nnz(){
    
    int dum=0;
    
    for(int i=0;i<this->pnz_size;i++){
        dum += this->pnz[1][i];
    }
    
    return dum;
    
}
void spvec_array::del_vap(){
	if(!vap.empty()){
		vap[0].clear();
		vap[1].clear();
		vap.clear();
	}
}
void poly_info::disp(){
    
    cout<<" polynomial."<<this->no<<" ";
    if(!coef.empty()){
		cout << endl;
		for(int i=0;i<this->numMs;i++){
    	    cout<<" coef.";
        	for(int j=0;j<this->sizeCone;j++){
            	cout.width(4); cout<<this->coef[i][j];
	        }
    	    cout<<" sup.i.";
        	for(int j=this->sup.pnz[0][i];j<this->sup.pnz[0][i] + this->sup.pnz[1][i];j++){
            	cout<<" "<<this->sup.vap[0][j]+1;
	        }
    	    cout<<" v.";
        	for(int j=this->sup.pnz[0][i];j<this->sup.pnz[0][i] + this->sup.pnz[1][i];j++){
            	cout<<" "<<this->sup.vap[1][j];
	        }	
    	    cout<<endl;
	    }
    }else{
		cout << "### EMPTY! ###" << endl;
	}
}
void spvec_array::write(string fname){
    
    //file opens.
	std::ofstream fout;
    fout.open(fname.c_str(),ios::out);
    fout.close();
    fout.open(fname.c_str(), ios::app );
    if( fout.fail() ){// check whether file is opened or not.
        cout << "error:file not open for output" << endl;
        cout << fname;
        exit(1);
    }
    
    // mDim
    fout<<this->dim<<endl;
    // pnz_size
    fout<<this->pnz_size<<endl;
    // vap_size
    fout<<this->vap_size<<endl;
    cout<<endl;
    
    // pnz[0],pnz[1],bij[0],bij[1],bij[2],coef
    for(int i=0;i<this->pnz_size;i++){
        fout<<this->pnz[0][i]<<" ";
        fout<<this->pnz[1][i]<<" ";
        fout<<endl;
    }
    fout<<endl;
    
    // vap[0],vap[1]
    for(int i=0;i<this->vap_size;i++){
        fout<<this->vap[0][i]<<" ";
        fout<<this->vap[1][i]<<" ";
        fout<<endl;
    }
    fout<<endl;
    
    fout.close();
    
}
bass_info::bass_info(){
    dim = 0;
	deg = 0;
	sup.alloc(0,0);
}
bass_info::bass_info(const bass_info &bassinfo){
	dim     = bassinfo.dim;
	deg     = bassinfo.deg;
	pattern = bassinfo.pattern;
	sup     = bassinfo.sup;
}
void bass_info::del(){
	//cout << "delete bass_info()" << endl;
    this->dim = 0;
	if(!pattern.empty()){
		this->pattern.clear();
	}
    this->sup.del();
}
bass_info::~bass_info(){
	pattern.clear();
    this->del();
}
void minkovsum(/*IN*/class spvec_array & vecs1,class spvec_array & vecs2,/*OUT*/class spvec_array & rvecs){
    
    rvecs.del();
    rvecs.alloc(vecs1.pnz_size*vecs2.pnz_size,vecs2.pnz_size*vecs1.vap_size + vecs1.pnz_size*vecs2.vap_size);
    
    int size1 = vecs1.pnz_size;
    int size2 = vecs2.pnz_size;
    
    int i=0;
    int j=0;
    int k=0;
    
    int idx=0;
    
    while(i<size1){
        int size3=0;
        while(j<size2){
            int pos1 = vecs1.pnz[0][i];
            int max1 = pos1 + vecs1.pnz[1][i];
            int pos2 = vecs2.pnz[0][j];
            int max2 = pos2+vecs2.pnz[1][j];
            
            // summation of supports
            if(pos1 >= 0 && pos2 >= 0){
                rvecs.pnz[0][k] = idx;
                int nnz3=0;
                while(pos1<max1 && pos2<max2){
                    if(vecs1.vap[0][pos1] < vecs2.vap[0][pos2]){
                        rvecs.vap[0][idx] =	vecs1.vap[0][pos1];
                        rvecs.vap[1][idx] = vecs1.vap[1][pos1];
                        idx++;
                        nnz3++;
                        pos1++;
                    }
                    else if(vecs1.vap[0][pos1] > vecs2.vap[0][pos2]){
                        //cout << "pos2 = " << pos2 << endl;
						//cout << "size of vecs2.vap = " << vecs2.vap_full_size << endl;
						rvecs.vap[0][idx] =	vecs2.vap[0][pos2];
                        rvecs.vap[1][idx] = vecs2.vap[1][pos2];
                        idx++;
                        nnz3++;
                        pos2++;
                    }
                    else {
                        rvecs.vap[0][idx] =	vecs1.vap[0][pos1];
                        rvecs.vap[1][idx] = vecs1.vap[1][pos1] + vecs2.vap[1][pos2];
                        idx++;
                        nnz3++;
                        pos1++;
                        pos2++;
                    }
                }
                while(pos1 < max1){
                    rvecs.vap[0][idx] =	vecs1.vap[0][pos1];
                    rvecs.vap[1][idx] = vecs1.vap[1][pos1];
                    idx++;
                    nnz3++;
                    pos1++;
                }
                while(pos2 < max2){
                    rvecs.vap[0][idx] =	vecs2.vap[0][pos2];
                    rvecs.vap[1][idx] = vecs2.vap[1][pos2];
                    idx++;
                    nnz3++;
                    pos2++;
                }
                rvecs.pnz[1][k] = nnz3;
            }
            else if(pos1 >= 0 && pos2 < 0){
                rvecs.pnz[0][k] = idx;
                while(pos1 < max1){
                    rvecs.vap[0][idx] =	vecs1.vap[0][pos1];
                    rvecs.vap[1][idx] = vecs1.vap[1][pos1];
                    idx++;
                    pos1++;
                }
                rvecs.pnz[1][k] = vecs1.pnz[1][i];
            }
            else if(pos2 >= 0 && pos1 < 0){
                rvecs.pnz[0][k] = idx;
                while(pos2 < max2){
                    rvecs.vap[0][idx] =	vecs2.vap[0][pos2];
                    rvecs.vap[1][idx] = vecs2.vap[1][pos2];
                    idx++;
                    pos2++;
                }
                rvecs.pnz[1][k] = vecs2.pnz[1][j];
            }
            else {// pos1 < 0 && pos2 < 0
                rvecs.pnz[0][k] = -1;
                rvecs.pnz[1][k] =  0;
            }
            j++;
            k++;
        }
        i++;
        j=0;
    }
    rvecs.pnz_size = k;
    rvecs.vap_size = idx;
}
void bass_info::alloc_pattern(int length){
    this->pattern.resize(length);
}
void pushsups(/*IN*/ class spvec_array insups,/*OUT*/ class spvec_array & outsups){
// outsups must be intialized before function pushsups.  
 
   	if(insups.pnz_size > 0){  
	    int maxvap = insups.pnz[0][insups.pnz_size-1] + insups.pnz[1][insups.pnz_size-1];
		//supports
	    int i=0;
    	while(i<insups.pnz_size){
        	outsups.pnz[0][outsups.pnz_size] = outsups.vap_size + insups.pnz[0][i];
	        outsups.pnz[1][outsups.pnz_size] = insups.pnz[1][i];
    	    outsups.pnz_size ++ ;
        	i++;
	    }
    
		// nonzero elements
	    i = 0;
		while(i<insups.vap_size){
        	outsups.vap[0][outsups.vap_size] = insups.vap[0][i];
	        outsups.vap[1][outsups.vap_size] = insups.vap[1][i];
    	    outsups.vap_size ++ ;
        	i++;
    	}
	}
}
void spvec_array::disp(int * slist){
    
    int i,j;
    cout<<" SIZE = "<<this->pnz_size <<endl;
    for(i=0; i<this->pnz_size; i++){
        
        cout<<"[ "<<  i     <<" ]";
        cout<<"[ "<<slist[i]<<" ]";
        if(this->pnz[0][slist[i]] >=0){
            cout<<" |var|";
            for(j=0;j<this->pnz[1][slist[i]];j++){
                cout<<" "<<this->vap[0][ this->pnz[0][slist[i]] + j ]+1;
            }
            cout<<" |pow|";
            for(j=0;j<this->pnz[1][slist[i]];j++){
                cout<<" "<<this->vap[1][this->pnz[0][slist[i]]+j];
            }
        }
        cout<<endl;
    }
    
}
void spvec_array::clean_pnz(){
   	if(!pnz.empty()){
		if(!pnz[0].empty()){
			pnz[0].clear();
		}
		if(!pnz[1].empty()){
			pnz[1].clear();
		}
		pnz.clear();
    }
    
}
void spvec_array::clean_vap(){
    
   	if(!vap.empty()){ 
		if(!vap[0].empty()){
			vap[0].clear();
		}
		if(!vap[1].empty()){
			vap[1].clear();
		}
		vap.clear();
    }
}
void spvec_array::clean(){
    
    this->clean_pnz();
    this->clean_vap();
    
}
void minkovsum(/*IN*/class spvec_array & vecs,/*OUT*/class spvec_array & rvecs)
{
    
    int bsize = vecs.pnz_size;
    rvecs.alloc(bsize*(bsize+1)/2,bsize*vecs.get_nnz());
    rvecs.pnz_size = 0;
    rvecs.vap_size = 0;
    int nnz3=0;
    int pos1,max1,pos2,max2;
    int vidx = 0;
    int sidx = 0;
    for(int i=0;i<bsize;i++){
        for(int j=i;j<bsize;j++){
            pos1 = vecs.pnz[0][i];
            max1 = pos1+vecs.pnz[1][i];
            pos2 = vecs.pnz[0][j];
            max2 = pos2+vecs.pnz[1][j];
            
            // summation of supports
			if(pos1 >= 0 && pos2 >= 0){
                nnz3=0;
                rvecs.pnz[0][sidx] = vidx;
                while(pos1<max1 && pos2<max2){
                    if(vecs.vap[0][pos1] < vecs.vap[0][pos2]){
                        rvecs.vap[0][vidx] = vecs.vap[0][pos1];
                        rvecs.vap[1][vidx] = vecs.vap[1][pos1];
                        vidx++;
                        nnz3++;
                        pos1++;
                    }
                    else if(vecs.vap[0][pos1] > vecs.vap[0][pos2]){
						rvecs.vap[0][vidx] = vecs.vap[0][pos2];
                        rvecs.vap[1][vidx] = vecs.vap[1][pos2];
                        vidx++;
                        nnz3++;
                        pos2++;
                    }
                    else {
                        rvecs.vap[0][vidx] = vecs.vap[0][pos1];
                        rvecs.vap[1][vidx] = vecs.vap[1][pos1] + vecs.vap[1][pos2];
                        vidx++;
                        nnz3++;
                        pos1++;
                        pos2++;
                    }
                }
                while(pos1 < max1){
                    rvecs.vap[0][vidx] = vecs.vap[0][pos1];
                    rvecs.vap[1][vidx] = vecs.vap[1][pos1];
                    vidx++;
                    nnz3++;
                    pos1++;
                }
                while(pos2 < max2){
                    rvecs.vap[0][vidx] = vecs.vap[0][pos2];
                    rvecs.vap[1][vidx] = vecs.vap[1][pos2];
                    vidx++;
                    nnz3++;
                    pos2++;
                }
                rvecs.pnz[1][sidx] = nnz3;
            }
            else if(pos1 >= 0 && pos2 < 0){
                rvecs.pnz[0][sidx] = vidx;
                while(pos1 < max1){
                    rvecs.vap[0][vidx] = vecs.vap[0][pos1];
                    rvecs.vap[1][vidx] = vecs.vap[1][pos1];
                    vidx++;
                    pos1++;
                }
                rvecs.pnz[1][sidx] = vecs.pnz[1][i];
            }
            else if(pos2 >= 0 && pos1 < 0){
                rvecs.pnz[0][sidx] = vidx;
                while(pos2 < max2){
                    rvecs.vap[0][vidx] = vecs.vap[0][pos2];
                    rvecs.vap[1][vidx] = vecs.vap[1][pos2];
                    vidx++;
                    pos2++;
                }
                rvecs.pnz[1][sidx] = vecs.pnz[1][j];
            }
            else {// pos1 < 0 && pos2 < 0
                rvecs.pnz[0][sidx] = -1;
                rvecs.pnz[1][sidx] =  0;
            }
            sidx++;
        }
    }
    rvecs.pnz_size = sidx;
    rvecs.vap_size = vidx;
}

// Let G be a set of supports. This function returns (G+G)\(2G).
void minkSumMinusDiagPart(/*IN*/class spvec_array & vecs,/*OUT*/class spvec_array & rvecs)
{
    int bsize = vecs.pnz_size;
    rvecs.alloc(bsize*(bsize-1)/2,bsize*vecs.get_nnz());
    rvecs.pnz_size = 0;
    rvecs.vap_size = 0;
    int nnz3=0;
    int pos1,max1,pos2,max2;
    int vidx = 0;
    int sidx = 0;
    for(int i=0;i<bsize;i++){
        for(int j=i+1;j<bsize;j++){
            pos1 = vecs.pnz[0][i];
            max1 = pos1+vecs.pnz[1][i];
            pos2 = vecs.pnz[0][j];
            max2 = pos2+vecs.pnz[1][j];
            //summation of supports
			if(pos1 >= 0 && pos2 >= 0){
                nnz3=0;
                rvecs.pnz[0][sidx] = vidx;
                while(pos1<max1 && pos2<max2){
                    if(vecs.vap[0][pos1] < vecs.vap[0][pos2]){
                        rvecs.vap[0][vidx] = vecs.vap[0][pos1];
                        rvecs.vap[1][vidx] = vecs.vap[1][pos1];
                        vidx++;
                        nnz3++;
                        pos1++;
                    }
                    else if(vecs.vap[0][pos1] > vecs.vap[0][pos2]){
                        rvecs.vap[0][vidx] = vecs.vap[0][pos2];
                        rvecs.vap[1][vidx] = vecs.vap[1][pos2];
                        vidx++;
                        nnz3++;
                        pos2++;
                    }
                    else {
                        rvecs.vap[0][vidx] = vecs.vap[0][pos1];
                        rvecs.vap[1][vidx] = vecs.vap[1][pos1] + vecs.vap[1][pos2];
                        vidx++;
                        nnz3++;
                        pos1++;
                        pos2++;
                    }
                }
                while(pos1 < max1){
                    rvecs.vap[0][vidx] = vecs.vap[0][pos1];
                    rvecs.vap[1][vidx] = vecs.vap[1][pos1];
                    vidx++;
                    nnz3++;
                    pos1++;
                }
                while(pos2 < max2){
                    rvecs.vap[0][vidx] = vecs.vap[0][pos2];
                    rvecs.vap[1][vidx] = vecs.vap[1][pos2];
                    vidx++;
                    nnz3++;
                    pos2++;
                }
                rvecs.pnz[1][sidx] = nnz3;
            }
            else if(pos1 >= 0 && pos2 < 0){
                rvecs.pnz[0][sidx] = vidx;
                while(pos1 < max1){
                    rvecs.vap[0][vidx] = vecs.vap[0][pos1];
                    rvecs.vap[1][vidx] = vecs.vap[1][pos1];
                    vidx++;
                    pos1++;
                }
                rvecs.pnz[1][sidx] = vecs.pnz[1][i];
            }
            else if(pos2 >= 0 && pos1 < 0){
                rvecs.pnz[0][sidx] = vidx;
                while(pos2 < max2){
                    rvecs.vap[0][vidx] = vecs.vap[0][pos2];
                    rvecs.vap[1][vidx] = vecs.vap[1][pos2];
                    vidx++;
                    pos2++;
                }
                rvecs.pnz[1][sidx] = vecs.pnz[1][j];
            }
            else {// pos1 < 0 && pos2 < 0
                rvecs.pnz[0][sidx] = -1;
                rvecs.pnz[1][sidx] =  0;
            }
            sidx++;
        }
    }
    rvecs.pnz_size = sidx;
    rvecs.vap_size = vidx;
}
int spvec_array::deg(){
    
    int degree=0;
    int dum;
    
    for(int i=0;i<this->pnz_size;i++){
        
        dum = 0;
        int a  = this->pnz[0][i];
        int ma = a + this->pnz[1][i];
        
        while(a >=0 && a < ma){
            dum += this->vap[1][a];
            a++;
        }
        
        if( degree < dum ){
            degree = dum;
        }
        
    }
    
    return degree;
    
}
void convert_var_pattern(vector<int> pattern,class spvec_array & vecs){
    
    for(int i=0;i<vecs.vap_size;i++){
        vecs.vap[0][i] = pattern[ vecs.vap[0][i] ];
    }
    
}
