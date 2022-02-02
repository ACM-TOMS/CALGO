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

#include "conversion.h"

bool comp_sup_for_unique(class sup sup1, class sup sup2){
    return sup1 < sup2;
}
void s3r::redundant_OneBounds(class supsetSet BasisSupports, class supSet allSup, class supSet & OneSup){
    int nDim = Polysys.dimVar;
    int mDim = Polysys.numSys;
    int kDim = BasisSupports.supsetArray.size()-mDim;
    class spvec_array spvecs,tarray,spSet, cspv,ei;
    class supSet momentSup,redSup;
    int i,j,k,m;
    int ssize, nnz,idx, tsize, tnnz;
	list<int>::iterator it;
    //pick up diagonal parts of all moment matrices
    
    nnz = 0;
    ssize = 0;
    for(i = 0;i < kDim; i++){
        tnnz  =  BasisSupports.supsetArray[i+mDim].nnz();
        nnz  = nnz + tnnz*tnnz;
        tsize =  BasisSupports.supsetArray[i+mDim].size();
        ssize = ssize + tsize*(tsize+1)/2;
    }
    for(i=0;i < nDim;i++){
        // find local mat. having index i.
		for(it=Polysys.posOfubds[i].begin();it!=Polysys.posOfubds[i].end();++it){
			//nnzVec[i]++;
        	idx = (*it);
	        nnz   = nnz   + BasisSupports.supsetArray[idx-1].nnz();
    	    nnz   = nnz   + BasisSupports.supsetArray[idx-1].size();
        	ssize = ssize + BasisSupports.supsetArray[idx-1].size();
		}		
    }
    cspv.alloc(ssize,nnz);
    cspv.dim = nDim;
    nnz = 0;
    ssize = 0;
    //pick up diagonal parts of all moment matrices
    for(i = 0;i < kDim; i++){
        momentSup = BasisSupports.supsetArray[i+mDim];
        initialize_spvecs(momentSup,spvecs);
        minkSumMinusDiagPart(spvecs,tarray);
        pushsups(tarray, cspv);
        nnz = nnz + tarray.vap_size;
        ssize = ssize + tarray.pnz_size;
    }
    
    //pick up diagonal parts of all localizin mat of type ((1-x_i)* u_iu_i^T)
    ei.alloc(1,1);
    ei.pnz_size = 1;
    ei.vap_size = 1;
    ei.vap[1][0] = 1;
    ei.pnz[0][0] = 0;
    ei.pnz[1][0] = 1;
    for(i=0;i < nDim;i++){
        // make e_i vector with sparse format.
        ei.vap[0][0] = i;
		it = Polysys.posOfubds[i].begin();
		for(;it !=Polysys.posOfubds[i].end();++it){
			idx =(*it);
            momentSup = BasisSupports.supsetArray[idx-1];
           	initialize_spvecs(momentSup,spvecs);
	        for(k = 0;k <spvecs.vap_size ; k++){
    			spvecs.vap[1][k] = spvecs.vap[1][k]  << 1;
            }
           	minkovsum(spvecs,ei,tarray);
	        pushsups(tarray,cspv);
    	    nnz = nnz + tarray.vap_size;
            ssize = ssize + tarray.pnz_size;
        }
    }
    cspv.vap_size = nnz;
    cspv.pnz_size = ssize;
    initialize_supset(cspv,redSup);
    redSup.unique();
    list<sup>::iterator b1 = allSup.begin();
    list<sup>::iterator e1 = allSup.supList.end();
    list<sup>::iterator b2 = redSup.begin();
    list<sup>::iterator e2 = redSup.supList.end();
    list<sup>::iterator r1 = OneSup.begin();
    set_difference(b1,e1,b2,e2,inserter(OneSup.supList,r1),comp_sup_for_unique);
    OneSup.sort();
    //cout << "Size of OneSup = " << OneSup.size() << endl;
    //OneSup.disp();
}

void s3r::redundant_ZeroBounds(class supsetSet BasisSupports, class supSet allSup, class supSet & ZeroSup){
    int nDim = Polysys.dimVar;
    int mDim = Polysys.numSys;
    int kDim = BasisSupports.supsetArray.size()-mDim;
    class supSet redSup, momentSup;
    class spvec_array spvecs, cspv,aspv,ei, tarray;
    //vector<int> nnzVec(nDim);
    int i,j,k,m,ssize, nnz,idx;
    list<sup>::iterator ite;
    list<int>::iterator it;
    
    // counting the number of monomials in supsets.
    ssize = 0;
    nnz = 0;
    for(i = 0;i < kDim; i++){
        nnz   = nnz   + BasisSupports.supsetArray[i+mDim].nnz();
        ssize = ssize + BasisSupports.supsetArray[i+mDim].size();
    }
    for(i=0;i < nDim;i++){
        // find local mat. having index i.
		for(it = Polysys.posOflbds[i].begin();it!=Polysys.posOflbds[i].end();++it){
                idx = (*it);
                nnz   = nnz   + BasisSupports.supsetArray[idx-1].nnz();
                nnz   = nnz   + BasisSupports.supsetArray[idx-1].size();
                ssize = ssize + BasisSupports.supsetArray[idx-1].size();
		}	
    }
    cspv.alloc(ssize,nnz);
    cspv.dim = nDim;
    nnz = 0;
    //pick up diagonal parts of all moment matrices
    for(i = 0;i < kDim; i++){
        momentSup = BasisSupports.supsetArray[i+mDim];
        initialize_spvecs(momentSup,spvecs);
        for(j = 0;j <spvecs.vap_size ; j++){
            spvecs.vap[1][j] = spvecs.vap[1][j]  << 1;
        }
        nnz = nnz + spvecs.vap_size;
        pushsups(spvecs, cspv);
    }
    //pick up diagonal parts of all localizin mat of type (x_i u_iu_i^T)
    ei.alloc(1,1);
    ei.pnz_size = 1;
    ei.vap_size = 1;
    ei.vap[1][0] = 1;
    ei.pnz[0][0] = 0;
    ei.pnz[1][0] = 1;
    for(i=0;i < nDim;i++){
        // make e_i vector with sparse format.
        ei.vap[0][0] = i;
		it = Polysys.posOflbds[i].begin();
        for(;it!= Polysys.posOflbds[i].end();++it){
			idx = (*it);
            momentSup = BasisSupports.supsetArray[idx-1];
            initialize_spvecs(momentSup,spvecs);
            for(k = 0;k <spvecs.vap_size ; k++){
                spvecs.vap[1][k] = spvecs.vap[1][k]  << 1;
            }
            minkovsum(spvecs,ei,tarray);
            pushsups(tarray, cspv);
            nnz = nnz + tarray.vap_size;
        }
    }
    cspv.vap_size = nnz;
    cspv.pnz_size = ssize;
    initialize_supset(cspv,redSup);
    redSup.unique();
    list<sup>::iterator b1 = allSup.begin();
    list<sup>::iterator e1 = allSup.supList.end();
    list<sup>::iterator b2 = redSup.begin();
    list<sup>::iterator e2 = redSup.supList.end();
    list<sup>::iterator r1 = ZeroSup.begin();
    set_difference(b1,e1,b2,e2,inserter(ZeroSup.supList,r1),comp_sup_for_unique);
    ZeroSup.sort();
    //cout << "Size of ZeroSup = " << ZeroSup.size() << endl;
    //ZeroSup.disp();
}

mat_info::mat_info(){
    //this->bij = NULL;
}

void mat_info::del(){
    this->bij.clear();
    this->coef.clear();
    this->sup.del();
}
mat_info::~mat_info(){
    this->del();
}

void info_a_nnz_a_struct_size_eq(class poly_info & polyinfo,class spvec_array & bassinfo,int & info_size,int & nnz_size,int & blst_size){
    
    info_size=0;
    nnz_size =0;
    
    int mcsize = polyinfo.numMs;
    
    for(int i=0;i<polyinfo.sizeCone;i++){
        for(int k=0;k<mcsize;k++){
            if(polyinfo.coef[k][i] != 0){
                info_size++;
            }
        }
    }
    
    nnz_size   = ( polyinfo.sup.pnz_size*bassinfo.get_nnz() +bassinfo.pnz_size*polyinfo.sup.get_nnz() ) << 1;
    info_size *= (bassinfo.pnz_size << 1);
    blst_size  = 1;
    
}
void info_a_nnz_a_struct_size_ineq_a_ba1( class poly_info & polyinfo, class spvec_array & bassinfo,int & info_size,int & nnz_size,int & blst_size){
    
    info_size=0;
    nnz_size=0;
    
    int mcsize = polyinfo.numMs;
    
    for(int i=0; i<mcsize; i++){
        for(int j=0; j<polyinfo.sizeCone; j++){
            if(polyinfo.coef[i][j] != 0){
                info_size++;
                nnz_size += polyinfo.sup.pnz[1][i];
            }
        }
    }
    
    blst_size =1;
    
}
void info_a_nnz_a_struct_size_ba1mmt( class spvec_array & bassinfo,int & info_size,int & nnz_size,int & blst_size){
    
    if( bassinfo.pnz[1][0] == 0 ){
        info_size = 0;
        nnz_size  = 0;
        blst_size = 0;
    }
    else{
        
        info_size = 1;
        nnz_size  = bassinfo.pnz[1][0];
        blst_size = 1;
    }
    
}
void info_a_nnz_a_struct_size_ineq_a_ba2( class poly_info & polyinfo, class spvec_array & bassinfo,int & info_size,int & nnz_size,int & blst_size){
    
    info_size=0;
    nnz_size =0;
    
    int ndum=0;
    
    int mna = bassinfo.pnz_size * bassinfo.vap_size;
    int mmsize = bassinfo.pnz_size * ( bassinfo.pnz_size + 1 ) / 2;
    int numms;
    
    for(int i=0;i<polyinfo.sizeCone;i++){
        numms=0;
        ndum=0;
        for(int j=0;j<polyinfo.numMs;j++){
            if(polyinfo.coef[j][i] != 0){
                
                info_size++;
                
                numms++;
                ndum += polyinfo.sup.pnz[1][j];
            }
        }
        nnz_size += numms * mna + mmsize * ndum;
    }
    
    info_size *= bassinfo.pnz_size*(bassinfo.pnz_size+1)/2;
    blst_size =polyinfo.sizeCone;
    
}
void info_a_nnz_a_struct_size_ba2mmt( class spvec_array & bassinfo,int & info_size,int & nnz_size,int & blst_size){
    
    info_size = bassinfo.pnz_size*(bassinfo.pnz_size+1)/2;
    nnz_size  = bassinfo.pnz_size * bassinfo.vap_size;
    blst_size = 1;
    
}
void info_a_nnz_a_struct_size_sdp(int mdim, class poly_info & polyinfo, class spvec_array & bassinfo,int & info_size,int & nnz_size,int & blst_size){
    
    //info_size = number of all supports
    //nnz_size  = nonzeros of all supports
    
    info_size=0;
    nnz_size =0;
    
    int k,s,totalcols;
    
    //get number of all supports
    int nnzMx,nnzDiagMx,nnzTotalMx;
    int nnzAa,nnzDiagAa;
    
    nnzDiagMx = bassinfo.pnz_size;
    nnzMx = nnzDiagMx*nnzDiagMx;
    
    k=0;
    
    nnzDiagAa = 0;
    totalcols = polyinfo.sizeCone * polyinfo.numMs;
    for(int i=0;i<polyinfo.numMs;i++){
        for(int j=0;j<polyinfo.sizeCone;j++){
            
            if( polyinfo.mc[k+1] - polyinfo.mc[k] ){
                if( polyinfo.mr[polyinfo.mc[k]] == j){
                    nnzDiagAa++;
                }
            }
            k++;
        }
    }
    
    nnzAa =  2 * polyinfo.mc[k] - nnzDiagAa;
    
    info_size = (nnzMx*nnzAa + nnzDiagMx*nnzDiagAa)/2;
    
    //get number of vairables of which the moment vector consist
    int maxlength = bassinfo.dim2;
    
    //get nonzeros fo all supports
    int totalnnz;
    k=0;
    for(int i=0;i<polyinfo.numMs;i++){
        
        nnzMx=0;
        nnzDiagMx=0;
        for(int j=0;j<bassinfo.pnz_size;j++){
            for(int t=j;t<bassinfo.pnz_size;t++){
                
                totalnnz = polyinfo.sup.pnz[1][i] + bassinfo.pnz[1][j] + bassinfo.pnz[1][t];
                if(totalnnz > maxlength)
                    totalnnz = maxlength;
                
                if(j==t){
                    nnzDiagMx += totalnnz;
                }
                nnzMx += totalnnz;
            }
        }
        nnzTotalMx = 2*nnzMx - nnzDiagMx;
        
        nnzAa=0;
        nnzDiagAa=0;
        for(int s=0;s<polyinfo.sizeCone;s++){
            if( polyinfo.mc[k+1] - polyinfo.mc[k] ){
                if( polyinfo.mr[polyinfo.mc[k]] == s){
                    nnzDiagAa++;
                }
                nnzAa += polyinfo.mc[k+1]-polyinfo.mc[k];
            }
            k++;
        }
        nnzAa = nnzAa - nnzDiagAa;
        
        nnz_size += nnzDiagAa*nnzMx + nnzAa*nnzTotalMx;
        
    }
    
    blst_size = 1;
    
}
void info_a_nnz_a_struct_size_obj(class poly_info & polyinfo,int & info_size,int & nnz_size,int & blst_size){
    
    info_size = polyinfo.sup.pnz_size;
    nnz_size  = polyinfo.sup.vap_size;
    blst_size = 1;
    
}

void get_info_a_nnz_a_struct_size(int mdim, int msize, vector<class poly_info> polyinfo, vector<class spvec_array> bassinfo,int & info_size,int & nnz_size,int & blst_size){
    
    info_size   = 0;
    nnz_size    = 0;
    blst_size   = 0;
    
    int idum;
    int ndum;
    int bdum;
    
    info_a_nnz_a_struct_size_obj(polyinfo[0],idum,ndum,bdum);
    info_size += idum;
    nnz_size  += ndum;
    blst_size += bdum;
    
    for(int i=1;i<msize;i++){
        idum = 0;
        ndum = 0;
        bdum = 0;
        
        if(polyinfo[i].typeCone == EQU){
            info_a_nnz_a_struct_size_eq(polyinfo[i],bassinfo[i],idum,ndum,bdum);
            //cout<<" size[ "<<i-1<<" ] info = "<<idum<<" nz = "<<ndum<<" @info_size_eq"<<endl;
        }
        else if(polyinfo[i].typeCone == INE && bassinfo[i].pnz_size == 1){
            info_a_nnz_a_struct_size_ineq_a_ba1(polyinfo[i],bassinfo[i],idum,ndum,bdum);
            //cout<<" size[ "<<i-1<<" ] info = "<<idum<<" nz = "<<ndum<<" @ineq_a_ba1"<<endl;
        }
        else if(polyinfo[i].typeCone == INE && bassinfo[i].pnz_size >= 2){
            info_a_nnz_a_struct_size_ineq_a_ba2(polyinfo[i],bassinfo[i],idum,ndum,bdum);
            //cout<<" size[ "<<i-1<<" ] info = "<<idum<<" nz = "<<ndum<<" @ineq_a_ba2"<<endl;
        }
        else if(polyinfo[i].typeCone == SDP){
            info_a_nnz_a_struct_size_sdp(mdim,polyinfo[i],bassinfo[i],idum,ndum,bdum);
            //cout<<" size[ "<<i-1<<" ] info = "<<idum<<" nz = "<<ndum<<" @sdp"<<endl;
        }
        else if(bassinfo[i].pnz_size == 1){
            info_a_nnz_a_struct_size_ba1mmt(bassinfo[i],idum,ndum,bdum);
            //cout<<" size[ "<<i-1<<" ] info = "<<idum<<" nz = "<<ndum<<" @ba1mmt"<<endl;
        }
        else if(bassinfo[i].pnz_size >= 2){
            info_a_nnz_a_struct_size_ba2mmt(bassinfo[i],idum,ndum,bdum);
            //cout<<" size[ "<<i-1<<" ] info = "<<idum<<" nz = "<<ndum<<" @ba2mmt"<<endl;
            
        }
        else {
            cout<<"error@get_info_size :: not available combinatio of typeCone and basupsize"<<endl;
            exit(1);
        }
        
        info_size += idum;
        nnz_size  += ndum;
        blst_size += bdum;
        
    }
    
}
void mysdp::set_struct_info(int matstruct,int pos,int nnz_size,int typecone){
    
    this->bLOCKsTruct[this->nBlocks]   = matstruct;
    this->block_info[0][this->nBlocks] = pos;
    this->block_info[1][this->nBlocks] = nnz_size;
    this->block_info[2][this->nBlocks] = typecone;
    
}
void convert_eq(  class poly_info & polyinfo,  class spvec_array & bassinfo,class mysdp & psdp){
    
    int size1 = polyinfo.sup.pnz_size;
    int size2 = bassinfo.pnz_size;
    int nnz3=0;
    
    int pos1,max1,pos2,max2;
    
    int idx  = psdp.ele.sup.vap_size;
    int idx2 = psdp.ele.sup.pnz_size;
    
    int	sno_f = idx2;
    int vno_f = idx;
    
    psdp.nBlocks++;
    
    int i,j;
    
    for(int s=0;s<polyinfo.sizeCone;s++){
        for(j=0;j<size2;j++){
            for(i=0;i<size1;i++){
                if(polyinfo.coef[i][s] != 0){
                    pos1 = polyinfo.sup.pnz[0][i];
                    max1 = pos1+polyinfo.sup.pnz[1][i];
                    
                    pos2 = bassinfo.pnz[0][j];
                    max2 = pos2+bassinfo.pnz[1][j];
                    
                    //computation of summation--
                    //In case that both supports are not zero.
                    if(pos1 >= 0 && pos2 >= 0){
                        
                        nnz3=0;
                        
                        psdp.ele.sup.pnz[0][idx2] = idx;
                        
                        while(pos1<max1 && pos2<max2){
                            
                            if(polyinfo.sup.vap[0][pos1] < bassinfo.vap[0][pos2]){
                                
                                psdp.ele.sup.vap[0][idx] = polyinfo.sup.vap[0][pos1];
                                psdp.ele.sup.vap[1][idx] = polyinfo.sup.vap[1][pos1];
                                idx++;
                                nnz3++;
                                
                                pos1++;
                                
                            }
                            else if(polyinfo.sup.vap[0][pos1] > bassinfo.vap[0][pos2]){
                                
                                psdp.ele.sup.vap[0][idx] = bassinfo.vap[0][pos2];
                                psdp.ele.sup.vap[1][idx] = bassinfo.vap[1][pos2];
                                idx++;
                                nnz3++;
                                
                                pos2++;
                                
                            }
                            else {
                                
                                psdp.ele.sup.vap[0][idx] = polyinfo.sup.vap[0][pos1];
                                psdp.ele.sup.vap[1][idx] = polyinfo.sup.vap[1][pos1] + bassinfo.vap[1][pos2];
                                idx++;
                                nnz3++;
                                
                                pos1++;
                                pos2++;
                                
                            }
                            
                        }
                        while(pos1 < max1){
                            
                            psdp.ele.sup.vap[0][idx] = polyinfo.sup.vap[0][pos1];
                            psdp.ele.sup.vap[1][idx] = polyinfo.sup.vap[1][pos1];
                            idx++;
                            nnz3++;
                            
                            pos1++;
                            
                        }
                        while(pos2 < max2){
                            
                            psdp.ele.sup.vap[0][idx] = bassinfo.vap[0][pos2];
                            psdp.ele.sup.vap[1][idx] = bassinfo.vap[1][pos2];
                            idx++;
                            nnz3++;
                            
                            pos2++;
                            
                        }
                        
                        psdp.ele.sup.pnz[1][idx2] = nnz3;
                        
                    }
                    else if(pos1 >= 0 && pos2 < 0){
                        psdp.ele.sup.pnz[0][idx2] = idx;
                        while(pos1 < max1){
                            psdp.ele.sup.vap[0][idx] = polyinfo.sup.vap[0][pos1];
                            psdp.ele.sup.vap[1][idx] = polyinfo.sup.vap[1][pos1];
                            idx++;
                            
                            pos1++;
                        }
                        psdp.ele.sup.pnz[1][idx2] = polyinfo.sup.pnz[1][i];
                    }
                    else if(pos2 >= 0 && pos1 < 0){
                        
                        psdp.ele.sup.pnz[0][idx2] = idx;
                        
                        while(pos2 < max2){
                            
                            psdp.ele.sup.vap[0][idx] = bassinfo.vap[0][pos2];
                            psdp.ele.sup.vap[1][idx] = bassinfo.vap[1][pos2];
                            idx++;
                            
                            pos2++;
                            
                        }
                        
                        psdp.ele.sup.pnz[1][idx2] = bassinfo.pnz[1][j];
                    }
                    else {// pos1 < 0 && pos2 < 0
                        psdp.ele.sup.pnz[0][idx2] = -1;
                        psdp.ele.sup.pnz[1][idx2] =  0;
                    }
                    
                    psdp.ele.bij[0][idx2]=psdp.nBlocks;
                    psdp.ele.bij[1][idx2]=j + 1 +s ;
                    psdp.ele.bij[2][idx2]=j + 1 +s ;
                    
                    psdp.ele.coef[idx2] = polyinfo.coef[i][s];
                    
                    idx2++;
                    
                }
                
            }
            
        }
    }
    
    //copy supports
    int size = idx - vno_f;
    for(int i=idx;i<idx+size;i++){
        psdp.ele.sup.vap[0][ i ] = psdp.ele.sup.vap[0][i - size ];
        psdp.ele.sup.vap[1][ i ] = psdp.ele.sup.vap[1][i - size ];
    }
    idx += size;
    //copy sup, bij and val
    size = idx2 - sno_f;
    int move_size = size2 + polyinfo.sizeCone -1;
    for(int i=idx2;i<idx2+size;i++){
        psdp.ele.sup.pnz[0][i] = psdp.ele.sup.pnz[0][i - size];
        psdp.ele.sup.pnz[1][i] = psdp.ele.sup.pnz[1][i - size ];
        
        psdp.ele.bij[0][i] = psdp.ele.bij[0][i-size] ;
        psdp.ele.bij[1][i] = psdp.ele.bij[1][i-size] +  move_size;
        psdp.ele.bij[2][i] = psdp.ele.bij[2][i-size] +  move_size;
        
        psdp.ele.coef[i] = -psdp.ele.coef[i-size];
        
    }
    
    idx2 += size;
    
    psdp.set_struct_info(-2*polyinfo.sizeCone*bassinfo.pnz_size,sno_f,2*size,EQU);
    
    psdp.ele.sup.vap_size = idx;
    psdp.ele.sup.pnz_size = idx2;
    
    //psdp.disp(psdp.nBlocks -1, psdp.nBlocks);
    
    //cout<<"*****************************************"<<endl;
    //psdp.ele.sup.disp(sno_f,idx2);
    //cout<<"*****************************************"<<endl;
    
}
void convert_ineq_a_ba1(  class poly_info & polyinfo,  class spvec_array & bassinfo,class mysdp & psdp){
    
    int fsize = polyinfo.sup.pnz_size;
    int nnz3=0;
    
    int pos1,max1,pos2,max2;
    
    int vidx  = psdp.ele.sup.vap_size;
    int sidx = psdp.ele.sup.pnz_size;
    
    int	sno_f = sidx;
    
    psdp.nBlocks++;
    
    int i;
    
    for(int s=0;s<polyinfo.sizeCone;s++){
        for(i=0;i<fsize;i++){
            if(polyinfo.coef[i][s] != 0){
                pos1 = polyinfo.sup.pnz[0][i];
                max1 = pos1+polyinfo.sup.pnz[1][i];
                
                pos2 = bassinfo.pnz[0][0];
                max2 = pos2+bassinfo.pnz[1][0];
                
                //computation of summation--
                //In case that both supports are not zero.
                if(pos1 >= 0 && pos2 >= 0){
                    nnz3=0;
                    
                    psdp.ele.sup.pnz[0][sidx] = vidx;
                    
                    while(pos1<max1 && pos2<max2){
                        
                        if(polyinfo.sup.vap[0][pos1] < bassinfo.vap[0][pos2]){
                            
                            psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                            psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1];
                            vidx++;
                            nnz3++;
                            
                            pos1++;
                            
                        }
                        else if(polyinfo.sup.vap[0][pos1] > bassinfo.vap[0][pos2]){
                            
                            psdp.ele.sup.vap[0][vidx] = bassinfo.vap[0][pos2];
                            psdp.ele.sup.vap[1][vidx] = bassinfo.vap[1][pos2];
                            vidx++;
                            nnz3++;
                            
                            pos2++;
                            
                        }
                        else {
                            
                            psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                            psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1] + bassinfo.vap[1][pos2];
                            vidx++;
                            nnz3++;
                            
                            pos1++;
                            pos2++;
                            
                        }
                        
                    }
                    while(pos1 < max1){
                        
                        psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                        psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1];
                        vidx++;
                        nnz3++;
                        
                        pos1++;
                        
                    }
                    while(pos2 < max2){
                        
                        psdp.ele.sup.vap[0][vidx] = bassinfo.vap[0][pos2];
                        psdp.ele.sup.vap[1][vidx] = bassinfo.vap[1][pos2];
                        vidx++;
                        nnz3++;
                        
                        pos2++;
                        
                    }
                    
                    psdp.ele.sup.pnz[1][sidx] = nnz3;
                    
                }
                else if(pos1 >= 0 && pos2 < 0){
                    psdp.ele.sup.pnz[0][sidx] = vidx;
                    while(pos1 < max1){
                        psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                        psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1];
                        vidx++;
                        
                        pos1++;
                    }
                    psdp.ele.sup.pnz[1][sidx] = polyinfo.sup.pnz[1][i];
                }
                else if(pos2 >= 0 && pos1 < 0){
                    
                    psdp.ele.sup.pnz[0][sidx] = vidx;
                    
                    while(pos2 < max2){
                        
                        psdp.ele.sup.vap[0][vidx] = bassinfo.vap[0][pos2];
                        psdp.ele.sup.vap[1][vidx] = bassinfo.vap[1][pos2];
                        vidx++;
                        
                        pos2++;
                        
                    }
                    
                    psdp.ele.sup.pnz[1][sidx] = bassinfo.pnz[1][0];
                }
                else {// pos1 < 0 && pos2 < 0
                    psdp.ele.sup.pnz[0][sidx] = -1;
                    psdp.ele.sup.pnz[1][sidx] =  0;
                }
                
                psdp.ele.bij[0][sidx]=psdp.nBlocks;
                psdp.ele.bij[1][sidx]= s + 1 ;
                psdp.ele.bij[2][sidx]= s + 1 ;
                
                psdp.ele.coef[sidx] = polyinfo.coef[i][s];
                
                sidx++;
            }
        }
    }
    
    psdp.set_struct_info(-1*polyinfo.sizeCone,sno_f,sidx-sno_f,INE);
    psdp.ele.sup.pnz_size = sidx;
    psdp.ele.sup.vap_size = vidx;
    
}
void convert_ba1mmt(  class spvec_array & bassinfo,class mysdp & psdp){
    
    if(bassinfo.pnz[1][0] != 0){
        int vidx  = psdp.ele.sup.vap_size;
        int sidx = psdp.ele.sup.pnz_size;
        
        psdp.nBlocks++;
        
        psdp.ele.sup.pnz[0][sidx] = vidx;
        psdp.ele.sup.pnz[1][sidx] = bassinfo.pnz[1][0];
        
        psdp.ele.bij[0][sidx] =psdp.nBlocks;
        psdp.ele.bij[1][sidx] =1;
        psdp.ele.bij[2][sidx] =1;
        
        psdp.ele.coef[sidx] = 1;
        
        for(int i=0;i<bassinfo.pnz[1][0];i++){
            psdp.ele.sup.vap[0][vidx + i] = bassinfo.vap[0][i];
            psdp.ele.sup.vap[1][vidx + i] = 2*bassinfo.vap[1][i];
        }
        
        psdp.set_struct_info(-1,sidx,bassinfo.pnz[1][0],INE);
        psdp.ele.sup.vap_size += bassinfo.pnz[1][0];
        psdp.ele.sup.pnz_size++;
        
    }
    
}
void get_moment_matrix(  class spvec_array & bassinfo,class spvec_array & mm_mat){
    
    int bsize = bassinfo.pnz_size;
    
    mm_mat.alloc(bsize*(bsize+1)/2,bsize*bassinfo.get_nnz());
    mm_mat.pnz_size = 0;
    mm_mat.vap_size = 0;
    
    int nnz3=0;
    int pos1,max1,pos2,max2;
    
    int vidx = 0;
    int sidx = 0;
    
    int i;
    
    for(i=0;i<bsize;i++){
        for(int j=i;j<bsize;j++){
            pos1 = bassinfo.pnz[0][i];
            max1 = pos1+bassinfo.pnz[1][i];
            
            pos2 = bassinfo.pnz[0][j];
            max2 = pos2+bassinfo.pnz[1][j];
            
            //computation of summation--
            //In case that both supports are not zero.
            if(pos1 >= 0 && pos2 >= 0){
                
                nnz3=0;
                
                mm_mat.pnz[0][sidx] = vidx;
                
                while(pos1<max1 && pos2<max2){
                    
                    if(bassinfo.vap[0][pos1] < bassinfo.vap[0][pos2]){
                        
                        mm_mat.vap[0][vidx] = bassinfo.vap[0][pos1];
                        mm_mat.vap[1][vidx] = bassinfo.vap[1][pos1];
                        vidx++;
                        nnz3++;
                        
                        pos1++;
                        
                    }
                    else if(bassinfo.vap[0][pos1] > bassinfo.vap[0][pos2]){
                        
                        mm_mat.vap[0][vidx] = bassinfo.vap[0][pos2];
                        mm_mat.vap[1][vidx] = bassinfo.vap[1][pos2];
                        vidx++;
                        nnz3++;
                        
                        pos2++;
                        
                    }
                    else {
                        
                        mm_mat.vap[0][vidx] = bassinfo.vap[0][pos1];
                        mm_mat.vap[1][vidx] = bassinfo.vap[1][pos1] + bassinfo.vap[1][pos2];
                        vidx++;
                        nnz3++;
                        
                        pos1++;
                        pos2++;
                        
                    }
                    
                }
                while(pos1 < max1){
                    
                    mm_mat.vap[0][vidx] = bassinfo.vap[0][pos1];
                    mm_mat.vap[1][vidx] = bassinfo.vap[1][pos1];
                    vidx++;
                    nnz3++;
                    
                    pos1++;
                    
                }
                while(pos2 < max2){
                    
                    mm_mat.vap[0][vidx] = bassinfo.vap[0][pos2];
                    mm_mat.vap[1][vidx] = bassinfo.vap[1][pos2];
                    vidx++;
                    nnz3++;
                    
                    pos2++;
                    
                }
                
                mm_mat.pnz[1][sidx] = nnz3;
                
            }
            else if(pos1 >= 0 && pos2 < 0){
                mm_mat.pnz[0][sidx] = vidx;
                while(pos1 < max1){
                    mm_mat.vap[0][vidx] = bassinfo.vap[0][pos1];
                    mm_mat.vap[1][vidx] = bassinfo.vap[1][pos1];
                    vidx++;
                    
                    pos1++;
                }
                mm_mat.pnz[1][sidx] = bassinfo.pnz[1][i];
            }
            else if(pos2 >= 0 && pos1 < 0){
                
                mm_mat.pnz[0][sidx] = vidx;
                
                while(pos2 < max2){
                    
                    mm_mat.vap[0][vidx] = bassinfo.vap[0][pos2];
                    mm_mat.vap[1][vidx] = bassinfo.vap[1][pos2];
                    vidx++;
                    
                    pos2++;
                    
                }
                
                mm_mat.pnz[1][sidx] = bassinfo.pnz[1][j];
            }
            else {// pos1 < 0 && pos2 < 0
                mm_mat.pnz[0][sidx] = -1;
                mm_mat.pnz[1][sidx] =  0;
            }
            sidx++;
        }
    }
    
    mm_mat.pnz_size = sidx;
    mm_mat.vap_size = vidx;
    
}
void convert_ineq_a_ba2(  class poly_info & polyinfo,  class spvec_array & bassinfo,class mysdp & psdp){
    
    int bsize = bassinfo.pnz_size;
    
    //generate upper part of moment matrix
    class spvec_array mm_mat;
    get_moment_matrix(bassinfo,mm_mat);
    
    //set result of products of polynomial form and moment matrix
    int fsize = polyinfo.sup.pnz_size;
    int nnz3=0;
    
    int pos1,max1,pos2,max2;
    
    int i,j,k,u;
    
    for(int s=0;s<polyinfo.sizeCone;s++){
        
        int vidx  = psdp.ele.sup.vap_size;
        int sidx = psdp.ele.sup.pnz_size;
        
        int	snof = sidx;
        
        psdp.nBlocks++;
        u=0;
        for(j=0;j<bsize;j++){
            for(k=j;k<bsize;k++){
                
                for(i=0;i<fsize;i++){
                    if(polyinfo.coef[i][s] != 0){
                        
                        pos1 = polyinfo.sup.pnz[0][i];
                        max1 = pos1+polyinfo.sup.pnz[1][i];
                        
                        pos2 = mm_mat.pnz[0][u];
                        max2 = pos2+mm_mat.pnz[1][u];
                        
                        //computation of summation--
                        //In case that both supports are not zero.
                        if(pos1 >= 0 && pos2 >= 0){
                            
                            nnz3=0;
                            
                            psdp.ele.sup.pnz[0][sidx] = vidx;
                            
                            while(pos1<max1 && pos2<max2){
                                
                                if(polyinfo.sup.vap[0][pos1] < mm_mat.vap[0][pos2]){
                                    
                                    psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                                    psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1];
                                    vidx++;
                                    nnz3++;
                                    
                                    pos1++;
                                    
                                }
                                else if(polyinfo.sup.vap[0][pos1] > mm_mat.vap[0][pos2]){
                                    
                                    psdp.ele.sup.vap[0][vidx] = mm_mat.vap[0][pos2];
                                    psdp.ele.sup.vap[1][vidx] = mm_mat.vap[1][pos2];
                                    vidx++;
                                    nnz3++;
                                    
                                    pos2++;
                                    
                                }
                                else {
                                    
                                    psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                                    psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1] + mm_mat.vap[1][pos2];
                                    vidx++;
                                    nnz3++;
                                    
                                    pos1++;
                                    pos2++;
                                    
                                }
                                
                            }
                            while(pos1 < max1){
                                
                                psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                                psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1];
                                vidx++;
                                nnz3++;
                                
                                pos1++;
                                
                            }
                            while(pos2 < max2){
                                
                                psdp.ele.sup.vap[0][vidx] = mm_mat.vap[0][pos2];
                                psdp.ele.sup.vap[1][vidx] = mm_mat.vap[1][pos2];
                                vidx++;
                                nnz3++;
                                
                                pos2++;
                                
                            }
                            
                            psdp.ele.sup.pnz[1][sidx] = nnz3;
                            
                        }
                        else if(pos1 >= 0 && pos2 < 0){
                            psdp.ele.sup.pnz[0][sidx] = vidx;
                            while(pos1 < max1){
                                psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                                psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1];
                                vidx++;
                                
                                pos1++;
                            }
                            psdp.ele.sup.pnz[1][sidx] = polyinfo.sup.pnz[1][i];
                        }
                        else if(pos2 >= 0 && pos1 < 0){
                            
                            psdp.ele.sup.pnz[0][sidx] = vidx;
                            
                            while(pos2 < max2){
                                
                                psdp.ele.sup.vap[0][vidx] = mm_mat.vap[0][pos2];
                                psdp.ele.sup.vap[1][vidx] = mm_mat.vap[1][pos2];
                                vidx++;
                                
                                pos2++;
                                
                            }
                            
                            psdp.ele.sup.pnz[1][sidx] = mm_mat.pnz[1][u];
                        }
                        else {// pos1 < 0 && pos2 < 0
                            psdp.ele.sup.pnz[0][sidx] = -1;
                            psdp.ele.sup.pnz[1][sidx] =  0;
                        }
                        
                        psdp.ele.bij[0][sidx]=psdp.nBlocks;
                        psdp.ele.bij[1][sidx]= j + 1;
                        psdp.ele.bij[2][sidx]= k + 1;
                        
                        psdp.ele.coef[sidx] = polyinfo.coef[i][s];
                        
                        sidx++;
                        
                    }
                    
                }
                u++;
            }
            
        }
        
        psdp.set_struct_info(bsize,snof,sidx - snof,SDP);
        psdp.ele.sup.pnz_size = sidx;
        psdp.ele.sup.vap_size = vidx;
    }
    
}
void convert_ba2mmt(  class spvec_array & bassinfo,class mysdp & psdp){
    
    int bsize = bassinfo.pnz_size;
    
    int vidx  = psdp.ele.sup.vap_size;
    int sidx = psdp.ele.sup.pnz_size;
    
    int	snof = sidx;
    
    psdp.nBlocks++;
    
    int nnz3=0;
    int pos1,max1,pos2,max2;
    
    int i;
    
    for(i=0;i<bsize;i++){
        for(int j=i;j<bsize;j++){
            pos1 = bassinfo.pnz[0][i];
            max1 = pos1+bassinfo.pnz[1][i];
            
            pos2 = bassinfo.pnz[0][j];
            max2 = pos2+bassinfo.pnz[1][j];
            
            //computation of summation--
            //In case that both supports are not zero.
            if(pos1 >= 0 && pos2 >= 0){
                
                nnz3=0;
                
                psdp.ele.sup.pnz[0][sidx] = vidx;
                
                while(pos1<max1 && pos2<max2){
                    
                    if(bassinfo.vap[0][pos1] < bassinfo.vap[0][pos2]){
                        
                        psdp.ele.sup.vap[0][vidx] = bassinfo.vap[0][pos1];
                        psdp.ele.sup.vap[1][vidx] = bassinfo.vap[1][pos1];
                        vidx++;
                        nnz3++;
                        
                        pos1++;
                        
                    }
                    else if(bassinfo.vap[0][pos1] > bassinfo.vap[0][pos2]){
                        
                        psdp.ele.sup.vap[0][vidx] = bassinfo.vap[0][pos2];
                        psdp.ele.sup.vap[1][vidx] = bassinfo.vap[1][pos2];
                        vidx++;
                        nnz3++;
                        
                        pos2++;
                        
                    }
                    else {
                        
                        psdp.ele.sup.vap[0][vidx] = bassinfo.vap[0][pos1];
                        psdp.ele.sup.vap[1][vidx] = bassinfo.vap[1][pos1] + bassinfo.vap[1][pos2];
                        vidx++;
                        nnz3++;
                        
                        pos1++;
                        pos2++;
                        
                    }
                    
                }
                while(pos1 < max1){
                    
                    psdp.ele.sup.vap[0][vidx] = bassinfo.vap[0][pos1];
                    psdp.ele.sup.vap[1][vidx] = bassinfo.vap[1][pos1];
                    vidx++;
                    nnz3++;
                    
                    pos1++;
                    
                }
                while(pos2 < max2){
                    
                    psdp.ele.sup.vap[0][vidx] = bassinfo.vap[0][pos2];
                    psdp.ele.sup.vap[1][vidx] = bassinfo.vap[1][pos2];
                    vidx++;
                    nnz3++;
                    
                    pos2++;
                    
                }
                
                psdp.ele.sup.pnz[1][sidx] = nnz3;
                
            }
            else if(pos1 >= 0 && pos2 < 0){
                psdp.ele.sup.pnz[0][sidx] = vidx;
                while(pos1 < max1){
                    psdp.ele.sup.vap[0][vidx] = bassinfo.vap[0][pos1];
                    psdp.ele.sup.vap[1][vidx] = bassinfo.vap[1][pos1];
                    vidx++;
                    
                    pos1++;
                }
                psdp.ele.sup.pnz[1][sidx] = bassinfo.pnz[1][i];
            }
            else if(pos2 >= 0 && pos1 < 0){
                
                psdp.ele.sup.pnz[0][sidx] = vidx;
                
                while(pos2 < max2){
                    
                    psdp.ele.sup.vap[0][vidx] = bassinfo.vap[0][pos2];
                    psdp.ele.sup.vap[1][vidx] = bassinfo.vap[1][pos2];
                    vidx++;
                    
                    pos2++;
                    
                }
                
                psdp.ele.sup.pnz[1][sidx] = bassinfo.pnz[1][j];
            }
            else {// pos1 < 0 && pos2 < 0
                psdp.ele.sup.pnz[0][sidx] = -1;
                psdp.ele.sup.pnz[1][sidx] =  0;
            }
            
            psdp.ele.bij[0][sidx]=psdp.nBlocks;
            psdp.ele.bij[1][sidx]= i + 1;
            psdp.ele.bij[2][sidx]= j + 1;
            
            psdp.ele.coef[sidx] = 1.0;
            
            sidx++;
        }
    }
    
    psdp.set_struct_info(bsize,snof,sidx - snof,SDP);
    psdp.ele.sup.pnz_size = sidx;
    psdp.ele.sup.vap_size = vidx;
    
    //psdp.disp(psdp.nBlocks-1,psdp.nBlocks);
    
}
int write_sdp_polyinfo(string outname, class poly_info & polyinfo){
    
    if(polyinfo.typeCone != SDP)
        return false;
    
    //file open
    std::fstream fout( outname.c_str(), ios::out );
    if( fout.fail() )
        return false;
    
    //data out
    fout<<polyinfo.typeCone<<endl;
    fout<<polyinfo.sizeCone<<endl;
    fout<<polyinfo.numMs<<endl;
    fout<<polyinfo.mc[polyinfo.numMs*polyinfo.sizeCone]<<endl;
    fout<<polyinfo.numMs*polyinfo.sizeCone+1<<endl;
    for(int i=0;i<polyinfo.mc[polyinfo.numMs*polyinfo.sizeCone];i++){
        fout<<polyinfo.coef[i][0]<<" ";
    }
    fout<<endl;
    for(int i=0;i<polyinfo.mc[polyinfo.numMs*polyinfo.sizeCone];i++){
        fout<<polyinfo.mr[i]<<" ";
    }
    fout<<endl;
    for(int i=0;i<polyinfo.numMs*polyinfo.sizeCone+1;i++){
        fout<<polyinfo.mc[i]<<" ";
    }
    fout<<endl;
    
    if(!write_bassinfo("polysups",polyinfo.sup))
        return false;
    
    //file close
    fout.close();
    
    return true;
}
int write_bassinfo(string outname, class spvec_array & bassinfo){
    
    //file open
    std::fstream fout( outname.c_str(), ios::out );
    if( fout.fail() )
        return false;
    
    //data out
    fout<<bassinfo.dim<<endl;
    fout<<bassinfo.dim2<<endl;
    fout<<bassinfo.pnz_size<<endl;
    fout<<bassinfo.vap_size<<endl;
    if(bassinfo.pnz_size>0){
        for(int i=0;i<bassinfo.pnz_size;i++){
            fout<<bassinfo.pnz[0][i]<<" ";
        }
        fout<<endl;
        for(int i=0;i<bassinfo.pnz_size;i++){
            fout<<bassinfo.pnz[1][i]<<" ";
        }
        fout<<endl;
        if(bassinfo.vap_size>0){
            for(int i=0;i<bassinfo.vap_size;i++){
                fout<<bassinfo.vap[0][i]<<" ";
            }
            fout<<endl;
            for(int i=0;i<bassinfo.vap_size;i++){
                fout<<bassinfo.vap[1][i]<<" ";
            }
            fout<<endl;
        }
    }
    
    //file close
    fout.close();
    
    return true;
    
}
void convert_sdp(  class poly_info & polyinfo,  class spvec_array & bassinfo,class mysdp & psdp){
    
    int rowlength = bassinfo.pnz_size;
    
    //generate upper part of moment matrix
    class spvec_array momentmatrix;
    get_moment_matrix(bassinfo,momentmatrix);
    
    //get product of polynomial form and moment matrix;
    psdp.nBlocks++;
    
    int vidx  = psdp.ele.sup.vap_size;
    int sidx = psdp.ele.sup.pnz_size;
    int	snof = sidx;
    
    int orividx,orisidx;
    
    int nnz3=0;
    int pos1,max1,pos2,max2;
    int u,t,r;
    int rowsize;
    int colsize;
    u=0;
    t=0;
    for(int j=0;j<rowlength;j++){
        rowsize = j*polyinfo.sizeCone;
        for(int k=j;k<rowlength;k++){
            colsize = k*polyinfo.sizeCone;
            
            r=0;
            t=0;
            for(int i=0;i<polyinfo.sup.pnz_size;i++){
                
                //get product of support from polynomial form and support from moment matrix
                pos1 = polyinfo.sup.pnz[0][i];
                max1 = pos1+polyinfo.sup.pnz[1][i];
                
                pos2 = momentmatrix.pnz[0][u];
                max2 = pos2 + momentmatrix.pnz[1][u];
                
                //computation of summation--
                //In case that both supports are not zero.
                if(pos1 >= 0 && pos2 >= 0){
                    
                    nnz3=0;
                    
                    psdp.ele.sup.pnz[0][sidx] = vidx;
                    
                    while(pos1<max1 && pos2<max2){
                        
                        if(polyinfo.sup.vap[0][pos1] < momentmatrix.vap[0][pos2]){
                            
                            psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                            psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1];
                            vidx++;
                            nnz3++;
                            
                            pos1++;
                            
                        }
                        else if(polyinfo.sup.vap[0][pos1] > momentmatrix.vap[0][pos2]){
                            
                            psdp.ele.sup.vap[0][vidx] = momentmatrix.vap[0][pos2];
                            psdp.ele.sup.vap[1][vidx] = momentmatrix.vap[1][pos2];
                            vidx++;
                            nnz3++;
                            
                            pos2++;
                            
                        }
                        else {
                            
                            psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                            psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1] + momentmatrix.vap[1][pos2];
                            vidx++;
                            nnz3++;
                            
                            pos1++;
                            pos2++;
                            
                        }
                        
                    }
                    while(pos1 < max1){
                        
                        psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                        psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1];
                        vidx++;
                        nnz3++;
                        pos1++;
                    }
                    while(pos2 < max2){
                        
                        psdp.ele.sup.vap[0][vidx] = momentmatrix.vap[0][pos2];
                        psdp.ele.sup.vap[1][vidx] = momentmatrix.vap[1][pos2];
                        vidx++;
                        nnz3++;
                        pos2++;
                        
                    }
                    
                    psdp.ele.sup.pnz[1][sidx] = nnz3;
                    
                }
                else if(pos1 >= 0 && pos2 < 0){
                    psdp.ele.sup.pnz[0][sidx] = vidx;
                    while(pos1 < max1){
                        psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
                        psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1];
                        vidx++;
                        
                        pos1++;
                    }
                    psdp.ele.sup.pnz[1][sidx] = polyinfo.sup.pnz[1][i];
                }
                else if(pos2 >= 0 && pos1 < 0){
                    
                    psdp.ele.sup.pnz[0][sidx] = vidx;
                    
                    while(pos2 < max2){
                        
                        psdp.ele.sup.vap[0][vidx] = momentmatrix.vap[0][pos2];
                        psdp.ele.sup.vap[1][vidx] = momentmatrix.vap[1][pos2];
                        vidx++;
                        
                        pos2++;
                        
                    }
                    
                    psdp.ele.sup.pnz[1][sidx] = momentmatrix.pnz[1][u];
                }
                else {// pos1 < 0 && pos2 < 0
                    psdp.ele.sup.pnz[0][sidx] = -1;
                    psdp.ele.sup.pnz[1][sidx] =  0;
                }
                
                //set data for each nonzero value of coefficient matrix of polynomial form
                orisidx = sidx;
                orividx = psdp.ele.sup.pnz[0][sidx];
                max1 = orividx + psdp.ele.sup.pnz[1][sidx];
                if(orividx >= 0){
                    vidx = orividx;
                }
                if(j==k){
                    
                    for(int s=0;s<polyinfo.sizeCone;s++){
                        if( polyinfo.mc[t+1] - polyinfo.mc[t] > 0){
                            
                            while(r<polyinfo.mc[t+1]){
                                
                                //set block_number , row_index and col_index
                                psdp.ele.bij[0][sidx] = psdp.nBlocks;
                                psdp.ele.bij[1][sidx] = polyinfo.mr[r] + rowsize +1;
                                psdp.ele.bij[2][sidx] = s + colsize +1;
                                
                                psdp.ele.coef[sidx] = polyinfo.coef[r][0];
                                
                                //set(copy) value of indecies and nonzeros of supports
                                if(orividx >= 0){
                                    psdp.ele.sup.pnz[0][sidx] = vidx;
                                    psdp.ele.sup.pnz[1][sidx] = psdp.ele.sup.pnz[1][orisidx];
                                    
                                    for(int q = orividx; q < max1; q ++ ){
                                        psdp.ele.sup.vap[0][vidx] = psdp.ele.sup.vap[0][q];
                                        psdp.ele.sup.vap[1][vidx] = psdp.ele.sup.vap[1][q];
                                        vidx++;
                                    }
                                }
                                else{
                                    psdp.ele.sup.pnz[0][sidx] = -1;
                                    psdp.ele.sup.pnz[1][sidx] =  0;
                                }
                                sidx ++ ;
                                
                                r++;
                                
                            }
                            
                        }
                        t++;
                    }
                }
                else{
                    
                    max1 = orividx + psdp.ele.sup.pnz[1][sidx];
                    for(int s=0;s<polyinfo.sizeCone;s++){
                        if( polyinfo.mc[t+1] - polyinfo.mc[t] > 0){
                            while(r<polyinfo.mc[t+1]){
                                
                                //set block_number , row_index and col_index
                                psdp.ele.bij[0][sidx] = psdp.nBlocks;
                                psdp.ele.bij[1][sidx] = polyinfo.mr[r] + rowsize +1;
                                psdp.ele.bij[2][sidx] = s + colsize + 1;
                                
                                psdp.ele.coef[sidx] = polyinfo.coef[r][0];
                                
                                //set(copy) value of indecies and nonzeros of supports
                                if(orividx >= 0){
                                    psdp.ele.sup.pnz[0][sidx] = vidx;
                                    psdp.ele.sup.pnz[1][sidx] = psdp.ele.sup.pnz[1][orisidx];
                                    for(int q = orividx; q < max1; q ++ ){
                                        psdp.ele.sup.vap[0][vidx] = psdp.ele.sup.vap[0][q];
                                        psdp.ele.sup.vap[1][vidx] = psdp.ele.sup.vap[1][q];
                                        vidx++;
                                    }
                                }
                                else{
                                    psdp.ele.sup.pnz[0][sidx] = -1;
                                    psdp.ele.sup.pnz[1][sidx] =  0;
                                }
                                
                                sidx ++ ;
                                
                                if(polyinfo.mr[r] != s){
                                    //set block_number , row_index and col_index
                                    psdp.ele.bij[0][sidx] = psdp.nBlocks;
                                    psdp.ele.bij[1][sidx] = s + rowsize + 1;
                                    psdp.ele.bij[2][sidx] = polyinfo.mr[r] + colsize + 1;
                                    
                                    psdp.ele.coef[sidx] = polyinfo.coef[r][0];
                                    
                                    //set(copy) value of indecies and nonzeros of supports
                                    if( orividx >= 0 ){
                                        psdp.ele.sup.pnz[0][sidx] = vidx;
                                        psdp.ele.sup.pnz[1][sidx] = psdp.ele.sup.pnz[1][orisidx];
                                        for(int q = orividx; q < max1; q ++ ){
                                            psdp.ele.sup.vap[0][vidx] = psdp.ele.sup.vap[0][q];
                                            psdp.ele.sup.vap[1][vidx] = psdp.ele.sup.vap[1][q];
                                            vidx++;
                                        }
                                    }
                                    else{
                                        psdp.ele.sup.pnz[0][sidx] = -1;
                                        psdp.ele.sup.pnz[1][sidx] =  0;
                                    }
                                    
                                    sidx ++ ;
                                    
                                }
                                
                                r++;
                                
                            }
                        }
                        t++;
                    }
                    
                }
                
            }
            u++;
        }
    }
    
    psdp.set_struct_info(bassinfo.pnz_size*polyinfo.sizeCone,snof,sidx - snof,SDP);
    
    psdp.ele.sup.pnz_size = sidx;
    psdp.ele.sup.vap_size = vidx;
    
    
}
void convert_obj(  class poly_info & polyinfo,class mysdp & psdp){
    
    int fsize = polyinfo.sup.pnz_size;
    int pos1,max1;
    
    int vidx = 0;
    int sidx = 0;
    
    int	sno_f = 0;
    
    int i;
    
    psdp.nBlocks=0;
    
    for(i=0;i<fsize;i++){
        pos1 = polyinfo.sup.pnz[0][i];
        max1 = pos1+polyinfo.sup.pnz[1][i];
        
        psdp.ele.sup.pnz[0][sidx] = vidx;
        while(pos1 < max1){
            psdp.ele.sup.vap[0][vidx] = polyinfo.sup.vap[0][pos1];
            psdp.ele.sup.vap[1][vidx] = polyinfo.sup.vap[1][pos1];
            vidx++;
            pos1++;
        }
        psdp.ele.sup.pnz[1][sidx] = polyinfo.sup.pnz[1][i];
        psdp.ele.bij[0][sidx]= 0;
        psdp.ele.bij[1][sidx]= 1 ;
        psdp.ele.bij[2][sidx]= 1 ;
        psdp.ele.coef[sidx] = polyinfo.coef[i][0];
        sidx++;
    }
    
    psdp.set_struct_info(1,sno_f,sidx-sno_f,-999);
    
    /*
    cout << "nBlocks = " << psdp.nBlocks << endl;
    cout << psdp.block_info[0][psdp.nBlocks] << endl;
    cout << psdp.block_info[1][psdp.nBlocks] << endl;
    cout << psdp.block_info[2][psdp.nBlocks] << endl;
     */
    psdp.ele.sup.pnz_size = sidx;
    psdp.ele.sup.vap_size = vidx;
    
}

void get_psdp(/*IN*/int mdim,int msize, vector<class poly_info> polyinfo,vector<class spvec_array> bassinfo ,/*OUT*/class mysdp & psdp){
    
    
    //set number of variables
    psdp.mDim = mdim;
    
    int info_size;
    int nnz_size;
    int blst_size;
    
    //get number and nonzeros of all supports, and blocks.
    get_info_a_nnz_a_struct_size(mdim,msize,polyinfo,bassinfo,info_size,nnz_size,blst_size);
    psdp.idum1 = info_size;
    psdp.idum2 = nnz_size;
    psdp.idum3 = blst_size;
    
    psdp.alloc(blst_size,info_size,nnz_size);
    
    vector<vector<int> > ggg(2);
    int ng = 2;
    for(int i=0;i<ng;i++){
        ggg[i].resize(2,0);
        ggg[i][0] = -1;
        ggg[i][1] =  0;
    }
    
    convert_obj(polyinfo[0],psdp);
    
    for(int i=1;i<msize;i++){
        //int dum = info_size;
        //(negative) equality constraints
        if(polyinfo[i].typeCone == EQU){
            polyinfo[i].no = i;
            if(ggg[0][0] == -1){
                ggg[0][0] = psdp.nBlocks;
            }
            convert_eq(polyinfo[i],bassinfo[i],psdp);
            ggg[0][1] = psdp.nBlocks - ggg[0][0];
            //polyinfo[i].disp();
        }
        //(negative) inequality constraints with the size of basis supports 1
        else if(polyinfo[i].typeCone == INE && bassinfo[i].pnz_size == 1){
            polyinfo[i].no = i;
            if(ggg[1][0] == -1){
                ggg[1][0] = psdp.nBlocks;
            }
            convert_ineq_a_ba1(polyinfo[i],bassinfo[i],psdp);
            ggg[1][1] = psdp.nBlocks - ggg[1][0];
            //cout<<"mat[ "<<i<< " ] ineq_a_ba1"<<endl;
            //polyinfo[i].disp();
        }
        //(positive) inequality constraints with the size of basis supports more than 1
        else if(polyinfo[i].typeCone == INE && bassinfo[i].pnz_size >= 2){
            polyinfo[i].no = i;
            convert_ineq_a_ba2(polyinfo[i],bassinfo[i],psdp);
            //cout<<"mat[ "<<i<< " ] ineq_a_ba2"<<endl;
        }
        //(positive) SDP
        else if(polyinfo[i].typeCone == SDP){
            polyinfo[i].no = i;
            convert_sdp(polyinfo[i],bassinfo[i],psdp);
            //cout<<"mat[ "<<i<< " ] sdp"<<endl;
        }
        //(negative) Moment matrices with size 1
        else if(bassinfo[i].pnz_size == 1){
            //cout << bassinfo[i].pnz[0][0] << endl;
            //cout << bassinfo[i].pnz[1][0] << endl;
            polyinfo[i].no = i;
            if(ggg[1][0] == -1 && bassinfo[i].pnz[1][0] > 0){
                ggg[1][0] = psdp.nBlocks;
            }
            convert_ba1mmt(bassinfo[i],psdp);
            if(bassinfo[i].pnz[1][0] > 0){
                ggg[1][1] = psdp.nBlocks - ggg[1][0];
            }
            //cout<<"mat[ "<<i<< " ] ba1mmt"<<endl;
        }
        //(positive) Moment matrices with size > 1
        else if(bassinfo[i].pnz_size >= 2){
            polyinfo[i].no = i;
            convert_ba2mmt(bassinfo[i],psdp);
            //cout<<"mat[ "<<i<< " ] ba2mmt"<<endl;
        }
        // Error handling
        else {
            cout<<"error@get_psdp :: not available combinatio of typeCone and basupsize"<<endl;
            exit(1);
        }
        
    }
    //psdp.disp();
    /*
    cout << "ggg info = " << endl;
    cout << " ggg[0][0] = " << ggg[0][0] << endl;
    cout << " ggg[0][1] = " << ggg[0][1] << endl;
    cout << " ggg[1][0] = " << ggg[1][0] << endl;
    cout << " ggg[1][1] = " << ggg[1][1] << endl;
     */
    psdp.nBlocks ++;
    //psdp.disp(0,1);
    gather_diag_blocks(ng,ggg,psdp);
    //psdp.disp(0,1);
    //cout<<" <--- get_psdp <*** "<<endl;
    
}
void mysdp::alloc(int blst_size,int ele_size,int nnz_size){
    
    this->bLOCKsTruct.resize(blst_size,0);
    block_info.resize(3);
    block_info[0].resize(blst_size,0);
    block_info[1].resize(blst_size,0);
    block_info[2].resize(blst_size,0);
    ele.bij.resize(3);
    ele.bij[0].resize(ele_size,0);
    ele.bij[1].resize(ele_size,0);
    ele.bij[2].resize(ele_size,0);
    this->ele.sup.alloc(ele_size,nnz_size);
    this->ele.coef.resize(ele_size,0);
}
void mysdp::del(){
    //cout << "delete mysdp" <<endl;
    this->bLOCKsTruct.clear();
    this->block_info.clear();
    this->ele.del();
}
mysdp::mysdp(){
    this->nBlocks = 0;
}
mysdp::~mysdp(){
    this->del();
}
void mysdp::disp(){
    
    cout<<"psdp.data out--------------------------------"<<endl;
    cout<<"block[ "<<0<<" ] --> [ "<<this->nBlocks<<" ] "<<endl;
    cout<<"nBlocks = "<<this->nBlocks<<endl;
    
    for(int i=0;i<this->nBlocks+1;i++){
        cout<<"--- BLOCK [ "<<i<<" ] ----BlockStruct = "<< this->bLOCKsTruct[i]<<"--"<<endl;
        
        for(int j=this->block_info[0][i];j<this->block_info[0][i]+this->block_info[1][i];j++){
            cout<<"[ "<<j<<" ] ";
            cout<<"bij= ";
            cout<<this->ele.bij[0][j]<<" ";
            cout<<this->ele.bij[1][j]<<" ";
            cout<<this->ele.bij[2][j]<<" ";
            
            cout<<"coef= ";
            cout.precision(3);
            cout.width(7);
            cout<<this->ele.coef[j]<<" ";
            cout.width();
            cout<<"sup.idx ";
            for(int k=this->ele.sup.pnz[0][j];k<this->ele.sup.pnz[0][j]+this->ele.sup.pnz[1][j];k++){
                cout<<this->ele.sup.vap[0][k]<<" ";
            }
            cout<<"vap ";
            for(int k=this->ele.sup.pnz[0][j];k<this->ele.sup.pnz[0][j]+this->ele.sup.pnz[1][j];k++){
                cout<<this->ele.sup.vap[1][k]<<" ";
            }
            cout<<endl;
        }
    }
    
    cout<<"-----------------------------"<<endl<<endl;
    
}

void mysdp::disp(int b1,int b2){
    
    
    cout<<"psdp.data out-----------------"<<endl;
    
    cout<<"block[ "<<b1<<" ] --> [ "<<b2<<" ] "<<endl;
    cout<<"nBlocks = "<<this->nBlocks<<endl;
    if(b2 > this->nBlocks ){
        cout<<"error@mysdp::disp :: b2(="<<b2<<") is over(>"<<this->nBlocks<<")"<<endl;
        exit(1);
    }
    
    for(int i=b1;i<b2;i++){
        cout<<"--- BLOCK [ "<<i<<" ] ----------"<<endl;
        
        for(int j=this->block_info[0][i];j<this->block_info[0][i]+this->block_info[1][i];j++){
            cout<<"[ "<<j<<" ] ";
            cout<<"bij= ";
            cout<<this->ele.bij[0][j]<<" ";
            cout<<this->ele.bij[1][j]<<" ";
            cout<<this->ele.bij[2][j]<<" ";
            
            cout<<"coef= ";
            cout.precision(3);
            cout.width(7);
            cout<<this->ele.coef[j]<<" ";
            cout.width();
            cout<<"sup.idx ";
            for(int k=this->ele.sup.pnz[0][j];k<this->ele.sup.pnz[0][j]+this->ele.sup.pnz[1][j];k++){
                cout<<this->ele.sup.vap[0][k]<<" ";
            }
            cout<<"vap ";
            for(int k=this->ele.sup.pnz[0][j];k<this->ele.sup.pnz[0][j]+this->ele.sup.pnz[1][j];k++){
                cout<<this->ele.sup.vap[1][k]<<" ";
            }
            cout<<endl;
        }
    }
    cout<<"-------------------"<<endl<<endl;
    
}

void gather_diag_blocks(int gbs,vector<vector<int> > dibs,class mysdp & psdp){
    //cout<<" ***> gather_diag_blocks ---> "<<endl;
    int matsize;
    //psdp.disp();
    //gather diagonal parts in each block
    for(int i=0;i<gbs;i++){
        if(dibs[i][0] >=0 ){
            //In case of gathering blocks corresponding to free variables
            if(i==0){
                // copy information on position of supports and save other memroy.
                int totalnum=0;
                
                for(int b=dibs[i][0];b<=dibs[i][0]+dibs[i][1];b++){
                    totalnum += psdp.block_info[1][b];
                }
                vector<vector<int> > temppnz(2);
                temppnz[0].clear();
                temppnz[1].clear();
                temppnz[0].resize(totalnum,0);
                temppnz[1].resize(totalnum,0);
                
                for(int s =  psdp.block_info[0][dibs[i][0]+1];s <  psdp.block_info[0][dibs[i][0]+dibs[i][1]]+psdp.block_info[1][dibs[i][0]+dibs[i][1]];s ++){
                    temppnz[0][s] = psdp.ele.sup.pnz[0][s];
                    temppnz[1][s] = psdp.ele.sup.pnz[1][s];
                }
                
                //shift information of sup,bij,coef into other memory
                matsize = 0;
                int numsup = psdp.block_info[1][0] + (psdp.block_info[1][dibs[i][0]+1] >> 1);
                int tmp = 0;
                for(int b=dibs[i][0]+2;b<=dibs[i][0]+dibs[i][1];b++){
                    tmp = -psdp.bLOCKsTruct[b-1];
                    tmp = (tmp >> 1);
                    matsize += tmp;
                    //matsize += abs(psdp.bLOCKsTruct[b-1])/2;
                    for(int j=psdp.block_info[0][b];j<psdp.block_info[0][b]+(psdp.block_info[1][b] >> 1);j++){
                        
                        //sup
                        psdp.ele.sup.pnz[0][numsup] = temppnz[0][j];
                        psdp.ele.sup.pnz[1][numsup] = temppnz[1][j];
                        
                        //bij[1],bij[2]
                        psdp.ele.bij[1][numsup] = psdp.ele.bij[1][j] + matsize;
                        psdp.ele.bij[2][numsup] = psdp.ele.bij[2][j] + matsize;
                        
                        //coef
                        psdp.ele.coef[numsup] = psdp.ele.coef[j];
                        
                        numsup++;
                    }
                }
                tmp = -psdp.bLOCKsTruct[dibs[i][0]+dibs[i][1]];
                tmp = (tmp >> 1);
                matsize += tmp;
                //matsize += abs(psdp.bLOCKsTruct[dibs[i][0]+dibs[i][1]])/2;
                int numsup2 = psdp.block_info[1][0];
                for(int b=dibs[i][0]+1;b<=dibs[i][0]+dibs[i][1];b++){
                    for(int j=psdp.block_info[0][b] + (psdp.block_info[1][b] >> 1);j<psdp.block_info[0][b]+psdp.block_info[1][b];j++){
                        
                        //sup
                        psdp.ele.sup.pnz[0][numsup] = temppnz[0][j];
                        psdp.ele.sup.pnz[1][numsup] = temppnz[1][j];
                        
                        //bij[1],bij[2]
                        psdp.ele.bij[1][numsup] = psdp.ele.bij[1][numsup2] + matsize ;
                        psdp.ele.bij[2][numsup] = psdp.ele.bij[2][numsup2] + matsize ;
                        
                        //coef
                        psdp.ele.coef[numsup] = -1*psdp.ele.coef[numsup2];
                        
                        numsup ++;
                        numsup2++;
                    }
                    
                    psdp.block_info[0][b] = -1;
                    psdp.block_info[1][b] =  0;
                }
                
                psdp.block_info[0][dibs[i][0]+1] = psdp.block_info[1][0];
                psdp.block_info[1][dibs[i][0]+1] = numsup - psdp.block_info[1][0];
                psdp.bLOCKsTruct[dibs[i][0]+1] = -(matsize << 1);
                psdp.block_info[2][dibs[i][0]+1] = EQU;
            }
            //In case of gathering blocks corresponding to nonnegative variables
            else{
                matsize = 0;
                for(int b=dibs[i][0]+2;b<=dibs[i][0]+dibs[i][1];b++){
                    matsize += abs(psdp.bLOCKsTruct[b-1]);
                    for(int j=psdp.block_info[0][b];j<psdp.block_info[0][b]+psdp.block_info[1][b];j++){
                        psdp.ele.bij[1][j] += matsize;
                        psdp.ele.bij[2][j] += matsize;
                    }
                    
                    psdp.block_info[1][dibs[i][0]+1] += psdp.block_info[1][b];
                    
                    psdp.block_info[0][b] = -1;
                    psdp.block_info[1][b] =  0;
                    
                }
                matsize += abs(psdp.bLOCKsTruct[dibs[i][0]+dibs[i][1]]);
                /*
                cout << "size of psdp.bLOCKsTruct = " << psdp.bLOCKsTruct.size() << endl;
                cout << "index = " << dibs[i][0]+1 << endl;
                 */
                psdp.bLOCKsTruct[dibs[i][0]+1] = -1*matsize;
                psdp.block_info[2][dibs[i][0]+1] = INE;
                /*
                cout << "matsize = " << matsize << endl;
                for(int i=0; i<psdp.bLOCKsTruct.size(); i++){
                    cout << psdp.bLOCKsTruct[i] << endl;
                }
                 */
            }
        }
    }
    
    //renumbering no. of blocks
    int a = 1;
    for(int i=1;i<psdp.nBlocks;i++){
        if(psdp.block_info[1][i] > 0){
            psdp.block_info[0][a] = psdp.block_info[0][i];
            psdp.block_info[1][a] = psdp.block_info[1][i];
            psdp.block_info[2][a] = psdp.block_info[2][i];
            psdp.bLOCKsTruct[a]   = psdp.bLOCKsTruct[i];
            
            for(int j=psdp.block_info[0][a];j<psdp.block_info[0][a]+psdp.block_info[1][a];j++){
                psdp.ele.bij[0][j] = a;
            }
            
            a++;
        }
    }
    
    psdp.nBlocks = a;
    //psdp.disp();
    //cout<<" <--- gather_diag_blocks <*** "<<endl<<endl;
}
//void get_lsdp(class spvec_array & allsups,class mysdp & psdp,int & linearterms){
void get_lsdp(class spvec_array & allsups,class mysdp & psdp,vector<int> & linearterms){
    
    //double t1, t2, t3, t4, t5, t6, t7;
    //t1 = (double)clock();
    
    //cout<<" ***> get_lsdp ---> "<<endl;
    
    // allocate memory to save order of vector
    vector<int> slist(psdp.ele.sup.pnz_size);
    for(int i=0;i<psdp.ele.sup.pnz_size;i++){
        slist[i] = i;
    }
    //t2 = (double)clock();
    
    //sort PSDP
    qsort_psdp(slist,psdp);
    //t3 = (double)clock();
    
    // numbering variables
    //in this function, check whether some one order monomials are deleted or not.
    variable_numbering(allsups,slist,psdp,linearterms);
    //t4 = (double)clock();
    
    //computation of nonnzero elements of upper triangle matrices
    count_upper_nnz(slist,psdp);
    psdp.nBlocks --;
    //t5 = (double)clock();
    
    
    /*
    cout << (double)(t2-t1)/(double)CLOCKS_PER_SEC << "sec" << endl;
    cout << (double)(t3-t2)/(double)CLOCKS_PER_SEC << "sec" << endl;
    cout << (double)(t4-t3)/(double)CLOCKS_PER_SEC << "sec" << endl;
    cout << (double)(t5-t4)/(double)CLOCKS_PER_SEC << "sec" << endl;
     */
    //cout<<" <--- get_lsdp <*** "<<endl<<endl;
}
void mysdp::disp_lsdp(){
    
    cout<<"--- LSDP ------------------------------------"<<endl;
    cout<<"block[ "<<0<<" ] --> [ "<<this->nBlocks<<" ] "<<endl;
    cout<<"nBlocks = "<<this->nBlocks<<endl;
    cout<<"bLOCKsTRUCT: "<<endl;
    for(int i=0;i<=this->nBlocks;i++){
        cout<<" "<<this->bLOCKsTruct[i];
    }
    cout<<endl;
    
    for(int i=0;i<=this->nBlocks;i++){
        cout<<"--- BLOCK [ "<<i<<" ] ----BlockStruct = "<< this->bLOCKsTruct[i]<<"-----"<<endl;
        for(int j=this->block_info[0][i];j<this->block_info[0][i]+this->block_info[1][i];j++){
            cout<<" [ "<<j<<" ] ";
            cout<<" var= "<<this->ele.sup.pnz[0][j]<<" ";
            cout<<"bij= ";
            cout<<this->ele.bij[0][j]<<" ";
            cout<<this->ele.bij[1][j]<<" ";
            cout<<this->ele.bij[2][j]<<" ";
            
            cout<<"coef= ";
            cout.precision(3);
            cout.width(7);
            cout<<this->ele.coef[j]<<" ";
            cout.width();
            cout<<endl;
        }
    }
    
    cout<<"----------------------------------"<<endl<<endl;
    
}
void mysdp::disp_sparseformat(){
    
    cout<<" *** sparse format *** "<<endl;
    cout<<" mDim = "<<this->mDim<<endl;
    cout<<" nBlocks = "<<this->nBlocks-1;
    cout<<" bLOCKsTRUCT. "<<endl;
    for(int i=1;i<this->nBlocks;i++){
        cout.width(7); cout<<this->bLOCKsTruct[i];
        if( i% 5 == 0){
            cout<<endl;
        }
    }
    cout<<endl;
    cout<<" cvect."<<endl;
    for(int i=0;i<this->block_info[1][0];i++){
        cout.width(10); cout<<this->ele.sup.pnz[0][i];
        cout.width(5); cout<<this->ele.coef[i]<<endl;
    }
    cout<<endl;
    
    cout<<" sparse data."<<endl;
    int num=1;
    for(int i=1;i<this->nBlocks;i++){
        for(int j=this->block_info[0][i];j<this->block_info[0][i]+this->block_info[1][i];j++){
            cout.width(5); cout<<" "<<num++<<". ";
            cout.width(5); cout<<this->ele.sup.pnz[0][j];
            cout.width(5); cout<<this->ele.bij[0][j];
            cout.width(5); cout<<this->ele.bij[1][j];
            cout.width(5); cout<<this->ele.bij[2][j];
            
            cout.width(10); cout.precision(3); cout<<this->ele.coef[j];
            cout<<endl;
        }
    }
    
    
}

void mysdp::write(string fname){
    // This function does not output sdpa sparse format.
    // Function write_sdpa outputs sdpa sparse format for an obtained SDP.
    //file open
    std::ofstream fout;
    fout.open(fname.c_str(),ios::out);
    fout.close();
    fout.open(fname.c_str(), ios::app );
    if( fout.fail() ){//error check of opening file
        cout << "error:file not open for output" << endl;
        cout << fname;
        exit(1);
    }
    
    // mDim;
    fout<<this->mDim<<endl;
    // nBlocks;
    fout<<this->nBlocks<<endl;
    // pnz_size
    fout<<this->ele.sup.pnz_size<<endl;
    // vap_size
    fout<<this->ele.sup.vap_size<<endl;
    cout<<endl;
    //* bLOCKsTruct,block_info[0],block_info[1]
    for(int i=0;i<this->nBlocks;i++){
        fout<<this->bLOCKsTruct[i]<<" ";
        fout<<this->block_info[0][i]<<" ";
        fout<<this->block_info[1][i]<<" ";
        fout<<endl;
    }
    fout<<endl;
    
    // pnz[0],pnz[1],bij[0],bij[1],bij[2],coef
    for(int i=0;i<this->ele.sup.pnz_size;i++){
        fout<<this->ele.sup.pnz[0][i]<<" ";
        fout<<this->ele.sup.pnz[1][i]<<" ";
        fout<<this->ele.bij[0][i]<<" ";
        fout<<this->ele.bij[1][i]<<" ";
        fout<<this->ele.bij[2][i]<<" ";
        fout<<this->ele.coef[i]<<endl;
    }
    fout<<endl;
    
    // vap[0],vap[1]
    for(int i=0;i<this->ele.sup.vap_size;i++){
        fout<<this->ele.sup.vap[0][i]<<" ";
        fout<<this->ele.sup.vap[1][i]<<" ";
        fout<<endl;
    }
    fout<<endl;
    
    fout.close();
    
}
void mysdp::write_utnnz(string fname){
    
    // file open
    std::ofstream fout;
    fout.open(fname.c_str(),ios::out);
    fout.close();
    fout.open(fname.c_str(), ios::app );
    if( fout.fail() ){//error check
        cout << "error:file not open for output" << endl;
        cout << fname;
        exit(1);
    }
    
    fout<<this->utsize<<endl;
    for(int i=0;i<this->utsize;i++){
        fout<<this->utnnz[0][i]<<" "<<this->utnnz[1][i]<<" "<<this->utnnz[2][i]<<endl;
    }
    
    // mDim;
    fout.close();
    
}

bool comp_sup(class sup sup1,class sup sup2){
    
    int t = 0;
    int s = 0;
    
    if(sup1.nnz() != 0 && sup2.nnz() != 0){
        if(sup1.deg() < sup2.deg()){
            return true;
        }else if(sup1.deg() > sup2.deg()){
            return false;
        }else{
            while(t < sup1.nnz()&& s < sup2.nnz()){
                if(sup2.idx[s] < sup1.idx[t]){
                    return false;
                }else if(sup2.idx[s] > sup1.idx[t]){
                    return true;
                }else {
                    if(sup2.val[s] < sup1.val[t]){
                        return true;
                    }else if(sup2.val[s] > sup1.val[t]){
                        return false;
                    }else{
                        t++;
                        s++;
                    }
                }
            }
        }
        if(sup1.nnz() < sup2.nnz()){
            return true;
        }else if(sup1.nnz() > sup2.nnz()){
            return false;
        }else if(sup1.bij < sup2.bij){
            return true;
        }
        return false;
    }else if(sup1.nnz() != 0){
        return false;
    }else if(sup2.nnz() != 0){
        return true;
    }else if(sup1.bij < sup2.bij){
        return true;
    }
    return false;
}


void qsort_psdp(vector<int> & slist,class mysdp & psdp){
    
    class supSet supsets;
    initialize_supset(psdp.ele.sup,supsets);
    
    list<class sup>::iterator ite;
    int i = 0;
    for(ite=supsets.begin(); ite!=supsets.supList.end(); ++ite){
        (*ite).bij = psdp.ele.bij[0][i];
        (*ite).no = i;
        i++;
    }
    supsets.supList.sort(comp_sup);
    i=0;
    for(ite=supsets.begin(); ite!=supsets.supList.end(); ++ite){
        slist[i] = (*ite).no;
        i++;
    }
    supsets.clear();
}

void variable_numbering(class spvec_array & allsups,vector<int> plist,class mysdp & psdp,vector<int> & linearterms){
    //if this value is i+1, index i is a monomial with degree 1.
    //allsups has been already sorted and is unique.
    //cout<<" ***> variable_numbering ---> "<<endl;
    //psdp.disp();
    int as=0;
    int pp=0;
    /*
    allsups.disp();
    cout << "allsups.pnz_size = " << allsups.pnz_size << endl;
    cout << "allsups.pnz = " << endl;
    for(int i=0; i<allsups.pnz[0].size(); i++){
        cout << allsups.pnz[0][i] << ", " << allsups.pnz[1][i] << endl;
    }
    cout << endl;
    cout << "allsups.vap_size = " << allsups.vap_size << endl;
    cout << "allsups.vap = " << endl;
    for(int i=0; i<allsups.vap[0].size(); i++){
        cout << allsups.vap[0][i] << ", " << allsups.vap[1][i] << endl;
    }
    cout << endl;
     */
    int novar=0;
    
    if(allsups.pnz[1][0]!=0){
        novar=1;
    }
    //cout << "nover = " << novar << endl;
    int deg1terms = 0;
    int deg = 0;
    int apos,ppos;
    int amax,pmax;
    int idx;
    bool issame = true;
    bool doesuse = false;
    
    while(as < allsups.pnz_size && pp < psdp.ele.sup.pnz_size){
        //compare supports
        if( allsups.pnz[1][as] == psdp.ele.sup.pnz[1][plist[pp]]){
            
            apos = allsups.pnz[0][as];
            amax = allsups.pnz[0][as] + allsups.pnz[1][as];
            ppos = psdp.ele.sup.pnz[0][plist[pp]];
            pmax = psdp.ele.sup.pnz[0][plist[pp]] + psdp.ele.sup.pnz[1][plist[pp]];
            
            issame = true;
            
            while(apos < amax && ppos <pmax){
                //compare index
                if(allsups.vap[0][apos] == psdp.ele.sup.vap[0][ppos]){
                    //compare value
                    if(allsups.vap[1][apos] == psdp.ele.sup.vap[1][ppos]){
                        apos++;
                        ppos++;
                    }
                    else{
                        issame = false;
                        break;
                    }
                }
                else{
                    issame = false;
                    break;
                }
            }
        }
        else{
            issame = false;
        }
        
        //if same...
        if(issame){
            // substitute no. of variables
            psdp.ele.sup.pnz[0][plist[pp]] = novar;
            if(novar == 0){
                //cout << " obj.val. = " << psdp.ele.coef[plist[pp]]<< endl;
                psdp.ele.coef[plist[pp]] *= -1;
            }
            doesuse = true;
            pp++;
        }
        //otherwise
        else{
            if(doesuse == true){
                idx = allsups.pnz[0][as];
                if(allsups.pnz[1][as] == 1 && allsups.vap[1][idx] == 1){
                    linearterms[deg1terms] = allsups.vap[0][idx] + 1;
                    deg1terms++;
                }
                //increase no. of variables
                novar++;
                doesuse = false;
            }
            as++;
        }
    }
    
    if(doesuse == true){
        idx = allsups.pnz[0][as];
        if(allsups.pnz[1][as] == 1 && allsups.vap[1][idx] == 1){
            linearterms[deg1terms] = allsups.vap[0][idx] + 1;
            deg1terms++;
        }
        novar++;
    }
    
    psdp.mDim = novar-1;
    //cout<<" <--- variable_numbering <*** "<<endl<<endl;
    
}
void count_upper_nnz(vector<int> plist,class mysdp & psdp){
    
    //cout<<" ***> count_upper_nnz ---> "<<endl;
    
    vector<vector<int> > temp(3);
    for(int i=0;i<3;i++){
        temp[i].resize(psdp.ele.sup.pnz_size,0);
    }
    psdp.utnnz = temp;
    psdp.utsize=0;
    
    int nonobj_size=0;
    for(int i=0;i<psdp.ele.sup.pnz_size;i++){
        if(psdp.ele.bij[0][plist[i]] !=0){
            plist[nonobj_size] = plist[i];
            nonobj_size++;
        }
    }
    
    psdp.utnnz[0][psdp.utsize] = psdp.ele.sup.pnz[0][plist[0]];
    psdp.utnnz[1][psdp.utsize] = psdp.ele.bij[0][plist[0]];
    psdp.utnnz[2][psdp.utsize] = 1;
    
    
    
    for(int i=1;i<nonobj_size;i++){
        if(psdp.ele.sup.pnz[0][plist[i]] == psdp.ele.sup.pnz[0][plist[i-1]]){
            if(psdp.ele.bij[0][plist[i]] == psdp.ele.bij[0][plist[i-1]]){
                psdp.utnnz[2][psdp.utsize] ++ ;
            }
            else {
                psdp.utsize++;
                psdp.utnnz[0][psdp.utsize] = psdp.ele.sup.pnz[0][plist[i]];
                psdp.utnnz[1][psdp.utsize] = psdp.ele.bij[0][plist[i]];
                psdp.utnnz[2][psdp.utsize] = 1;
            }
        }
        else{
            psdp.utsize++;
            psdp.utnnz[0][psdp.utsize] = psdp.ele.sup.pnz[0][plist[i]];
            psdp.utnnz[1][psdp.utsize] = psdp.ele.bij[0][plist[i]];
            psdp.utnnz[2][psdp.utsize] = 1;
        }
    }
    
    psdp.utsize ++;
    
    //cout<<" <--- count_upper_nnz <*** "<<endl<<endl;
    
}
sup_a_block::sup_a_block(){
    
    //this->vap0 = NULL;
    //this->vap1 = NULL;
    
}
sup_a_block::~sup_a_block(){
    this->vap0.clear();
    this->vap1.clear();
}
void sup_a_block::alloc(int mdim){
    this->vap0.resize(mdim,0);
    this->vap1.resize(mdim,0);
}
void sup_a_block::input(class mysdp & psdp,int i){
    if(this->vap0.empty()){
        this->vap0.resize(psdp.mDim,0);
    }
    if(this->vap1.empty()){
        this->vap1.resize(psdp.mDim,0);
    }
    this->block = psdp.ele.bij[0][i];
    
    this->nnzsize = psdp.ele.sup.pnz[1][i];
    
    this->deg=0;
    for(int k=0;k<psdp.ele.sup.pnz[1][i];k++){
        //		cout<<":s:";
        this->vap0[k] = psdp.ele.sup.vap[0][k+psdp.ele.sup.pnz[0][i]];
        //		cout<<" "<<this->vap0[k]+1;
        //		cout<<" :v:";
        this->vap1[k] = psdp.ele.sup.vap[1][k+psdp.ele.sup.pnz[0][i]];
        //		cout<<" "<<this->vap1[k];
        this->deg += this->vap1[k];
    }
    //	cout<<endl;
    
}
void sup_a_block::input(class spvec_array & supset,int i){
    
    if(this->vap0.empty()){
        this->vap0.resize(supset.dim,0);
    }
    if(this->vap1.empty()){
        this->vap1.resize(supset.dim,0);
    }
    
    this->nnzsize = supset.pnz[1][i];
    
    this->deg=0;
    for(int k=0;k<supset.pnz[1][i];k++){
        //		cout<<":s:";
        this->vap0[k] = supset.vap[0][k+supset.pnz[0][i]];
        //		cout<<" "<<this->vap0[k]+1;
        //		cout<<" :v:";
        this->vap1[k] = supset.vap[1][k+supset.pnz[0][i]];
        //		cout<<" "<<this->vap1[k];
        this->deg += this->vap1[k];
    }
    //	cout<<endl;
    
}
bool comp_sup_spvecs(class sup sup1,class sup sup2){
    
    int t = 0;
    int s = 0;
    
    if(sup1.nnz() != 0 && sup2.nnz() != 0){
        if(sup1.deg() < sup2.deg()){
            return true;
        }else if(sup1.deg() > sup2.deg()){
            return false;
        }else{
            while(t < sup1.nnz()&& s < sup2.nnz()){
                if(sup2.idx[s] < sup1.idx[t]){
                    return false;
                }else if(sup2.idx[s] > sup1.idx[t]){
                    return true;
                }else {
                    if(sup2.val[s] < sup1.val[t]){
                        return true;
                    }else if(sup2.val[s] > sup1.val[t]){
                        return false;
                    }else{
                        t++;
                        s++;
                    }
                }
            }
        }
        if(sup1.nnz() < sup2.nnz()){
            return true;
        }
        return false;
    }else if(sup1.nnz() != 0){
        return false;
    }else if(sup2.nnz() != 0){
        return true;
    }
    return false;
}
void qsort_sups(vector<int> & slist,class spvec_array & spvecs){
    
    class supSet supsets;
    initialize_supset(spvecs,supsets);
    list<class sup>::iterator ite;
    int i = 0;
    for(ite=supsets.begin(); ite!=supsets.supList.end(); ++ite){
        (*ite).no = i;
        i++;
    }
    supsets.supList.sort(comp_sup_spvecs);
    i=0;
    for(ite=supsets.begin(); ite!=supsets.supList.end(); ++ite){
        slist[i] = (*ite).no;
        i++;
    }
    
}

void simplification(/*IN*/ class spvec_array & vecs){
    
    if(vecs.pnz_size > 0){
        
        int ssize = vecs.pnz_size;
        vector<int> slist(ssize);
        for(int i=0;i<ssize;i++){
            slist[i] = i;
        }
        //cout << "before sorting" << endl;
        qsort_sups(slist,vecs);
        
        //delete the same elements from the first if exists
        int k=0;
        int s=1;
        int a,b;
        int ma;
        bool flag;
        
        //cout << "after sorting" << endl;
        
        while(s<ssize){
            
            flag = true;
            
            if( vecs.pnz[1][slist[k]] == vecs.pnz[1][slist[s]] ){
                
                a  = vecs.pnz[0][slist[k]];
                ma = a + vecs.pnz[1][slist[k]];
                
                b = vecs.pnz[0][slist[s]];
                
                while(a < ma){
                    if( vecs.vap[0][a] != vecs.vap[0][b]){
                        flag = false;
                        break;
                    }
                    else if(vecs.vap[1][a] != vecs.vap[1][b]){
                        flag = false;
                        break;
                    }
                    a++;
                    b++;
                }
                
            }
            else{
                flag = false;
            }
            if(flag == false){
                //k-th supports not equal to s^th supports => changes s-th supports into k+1-th.
                k++;
                slist[k] = slist[s];
            }
            s++;
        }
        vector<vector<int> > dumarray(2);
        dumarray[0].clear();
        dumarray[0].resize(vecs.pnz_size,0);
        dumarray[1].clear();
        dumarray[1].resize(vecs.pnz_size,0);
        for(int i=0;i<vecs.pnz_size;i++){
            dumarray[0][i] = vecs.pnz[0][i];
            dumarray[1][i] = vecs.pnz[1][i];
        }
        vecs.pnz_size = k+1;
        for(int i=0;i<vecs.pnz_size;i++){
            vecs.pnz[0][i] = dumarray[0][slist[i]];
            vecs.pnz[1][i] = dumarray[1][slist[i]];
        }
        dumarray.clear();
    }
    
    
}
void get_removelist(/*IN*/ int * alist, class spvec_array & allsups,class spvec_array & removesups,/*OUT*/ int * rlist)
{
    
    int rp,ap,ma;
    int numsame;
    
    for(int i=0;i<removesups.pnz_size;i++){
        for(int j=0;j<allsups.pnz_size;i++){
            //check whether supports of allsups included in removesup exist or not.
            if(removesups.pnz[1][i] < allsups.pnz[1][alist[j]]){
                
                rp = removesups.pnz[0][i];
                
                ap = allsups.pnz[0][alist[j]];
                ma = ap + allsups.pnz[1][alist[j]];
                
                numsame = 0;
                
                while(ap < ma){
                    //equal to index
                    if(removesups.vap[0][rp] > allsups.vap[0][ap]){
                        ap ++ ;
                    }
                    else if(removesups.vap[0][rp] < allsups.vap[0][ap]){
                        break;
                    }
                    else{
                        // degree is included
                        if(removesups.vap[1][rp] <= allsups.vap[1][ap]){
                            numsame ++;
                            rp++;
                        }
                        else{
                            break;
                        }
                        ap++;
                    }
                }
                
                // includes!
                if(numsame == removesups.pnz[1][i]){
                    rlist[alist[j]] = false;
                }
            }
        }
        
    }
    
}
void remove_sups(class mysdp & psdp,class spvec_array & removesups){
    
    //cout<<"***> remove_sups ( for psdp) ---> "<<endl;
    
    int dum = 0;
    
    int movesize = 0;
    int loss;
    int nodiag;
    int firstIdx,maxIdx;
    int idum;
    
    psdp.ele.sup.pnz_size = 0;
    
    int maxblock = 0;
    for(int i=0;i<psdp.nBlocks;i++){
        if(maxblock < psdp.block_info[1][i]){
            maxblock = psdp.block_info[1][i];
        }
    }
    
    vector<bool> remainIdx(maxblock);
    remainIdx.clear();
    
    int noblock = 0;
    int nnznoblock = 0;
    
    vector<int> onpattern ;
    int temp,temp2;
    
    while(noblock < psdp.nBlocks){
        
        firstIdx = psdp.block_info[0][noblock];
        maxIdx   = firstIdx + psdp.block_info[1][noblock];
        
        int s,k;
        
        //1. remove unnecesarry supports
        for(s=0; s<psdp.block_info[1][noblock]; s++){
            remainIdx[s] = true;
        }
        
        //1-1.find positions of support in sups, included in removesups
        int r,a,ma;
        int numsame;
        
        for(r=0;r<removesups.pnz_size;r++){
            temp=-999;
            temp2 = -1;
            
            for( s=firstIdx; s<maxIdx; s++){
                
                if(remainIdx[s - firstIdx] == true && removesups.pnz[1][r] <= psdp.ele.sup.pnz[1][s]){
                    
                    k = removesups.pnz[0][r];
                    
                    a  = psdp.ele.sup.pnz[0][s];
                    ma = a + psdp.ele.sup.pnz[1][s];
                    
                    numsame = 0;
                    
                    while(a < ma){
                        //equal to index
                        if(removesups.vap[0][k] > psdp.ele.sup.vap[0][a]){
                            a ++ ;
                        }
                        else if(removesups.vap[0][k] < psdp.ele.sup.vap[0][a]){
                            break;
                        }
                        else{
                            //removesup.degree is less than or equal to sup.degree
                            if(removesups.vap[1][k] <= psdp.ele.sup.vap[1][a]){
                                numsame ++;
                                k++;
                            }
                            //removesup.degree is more than sup.degree
                            else{
                                break;
                            }
                            a++;
                        }
                    }
                    
                    // Yes it is included
                    if(numsame == removesups.pnz[1][r]){
                        
                        //1-2.remove unnecesarry supports
                        remainIdx[s - firstIdx] = false; dum++;
                        
                    }
                }
            }
        }
        
        onpattern.clear();
        onpattern.resize(abs(psdp.bLOCKsTruct[noblock])+1,-999);
        
        temp  = -1;
        temp2 =  1;
        
        for(s=firstIdx; s<maxIdx; s++){
            if(remainIdx[s-firstIdx] == true){
                if(temp != psdp.ele.bij[1][s] ){
                    
                    temp = psdp.ele.bij[1][s];
                    
                    onpattern[temp] = temp2;
                    temp2++;
                    
                }
            }
        }
        
        //2. Set size of matrices
        if(psdp.bLOCKsTruct[noblock] < 0){
            psdp.bLOCKsTruct[noblock] = 1 - temp2;
        }
        else{
            psdp.bLOCKsTruct[noblock] = temp2 -1 ;
        }
        
        //3.arrange no. of row and column
        loss = 0;	//<---defference of no. of row and col.
        nodiag = 1;	//<---no. of diagonal elements
        idum = movesize;
        
        psdp.block_info[0][noblock] -= movesize;
        
        temp = psdp.ele.bij[1][firstIdx];
        
        //2-1.reassign no. of row and col.
        for(s=firstIdx; s<maxIdx; s++){
            if(remainIdx[s-firstIdx] == true){
                
                psdp.ele.sup.pnz[0][s - movesize] = psdp.ele.sup.pnz[0][s];
                psdp.ele.sup.pnz[1][s - movesize] = psdp.ele.sup.pnz[1][s];
                
                psdp.ele.bij[0][s - movesize] = nnznoblock;
                psdp.ele.bij[1][s - movesize] = onpattern[psdp.ele.bij[1][s]];
                psdp.ele.bij[2][s - movesize] = onpattern[psdp.ele.bij[2][s]];
                
                psdp.ele.coef[s - movesize] = psdp.ele.coef[s];
                
            }
            else{
                movesize ++ ;
                psdp.block_info[1][noblock]--;
            }
        }
        
        if(psdp.block_info[1][noblock] !=0){
            psdp.bLOCKsTruct[nnznoblock]   = psdp.bLOCKsTruct[noblock];
            psdp.block_info[0][nnznoblock] = psdp.block_info[0][noblock];
            psdp.block_info[1][nnznoblock] = psdp.block_info[1][noblock];
            psdp.block_info[2][nnznoblock] = psdp.block_info[2][noblock];
            psdp.ele.sup.pnz_size += psdp.block_info[1][noblock];
            
            nnznoblock++;
        }
        noblock++;
    }
    psdp.nBlocks = nnznoblock;
    //cout<<"<--- remove_sups ( for psdp) <*** "<<endl<<endl;
}
void remove_sups(class spvec_array & removesups,class spvec_array & sups){
    
    //cout<<" ***> remove_sups( for spvecarray ) ---> "<<endl;
    
    int sv,mv;
    int rp,rv;
    int sp;
    int numsame=0;
    
    for(rp=0;rp<removesups.pnz_size;rp++){
        for(sp=0;sp<sups.pnz_size;sp++){
            if(sups.pnz[0][sp] >= 0){
                
                rv = removesups.pnz[0][rp];
                
                sv = sups.pnz[0][sp];
                mv = sv + sups.pnz[1][sp];
                
                numsame = 0;
                
                while(sv < mv){
                    //equal to index
                    if(removesups.vap[0][rv] > sups.vap[0][sv]){
                        sv ++ ;
                    }
                    else if(removesups.vap[0][rv] < sups.vap[0][sv]){
                        break;
                    }
                    else{
                        //removesup.degree is less than or equal to sup.degee
                        if(removesups.vap[1][rv] <= sups.vap[1][sv]){
                            numsame ++;
                            rv ++;
                        }
                        //removesup.degree is more than sup.degree
                        else{
                            break;
                        }
                        sv ++;
                    }
                }
                //Yes it is included
                if(numsame == removesups.pnz[1][rp]){
                    sups.pnz[0][sp] = -1;
                }
            }
        }
    }
    
    //arrange supports
    rp = 0;
    for(sp=0;sp<sups.pnz_size;sp++){
        if(sups.pnz[0][sp] >= 0){
            
            sups.pnz[0][rp] = sups.pnz[0][sp];
            sups.pnz[1][rp] = sups.pnz[1][sp];
            
            rp ++ ;
            
        }
    }
    
    sups.pnz_size = rp;
    
    //cout<<" <--- remove_sups( for spvecarray ) <*** "<<endl<<endl;
    
}



void s3r::set_relaxOrder(int Order){
    
    //cout<<" ***>s3r::set_relaxOrder ---> "<<endl;
    
    //find max degree
    int minOrder=(int)ceil((double)(this->Polysys.maxDeg())/2.0);
    
    if(Order<minOrder){
        this->param.relax_Order=minOrder;
    }
    else{
        param.relax_Order=Order;
    }
    
    //cout<<" <--- s3r::set_relaxOrder <*** "<<endl<<endl;
    
}

void gen_basisindices(
int sparsesw,int multifactor,
class polysystem & polysys,
class cliques & maxcliques,
vector<list<int> > & BasisIndices)
{
/*
    printf("numcliques = %2d\n", maxcliques.numcliques);
    printf("numnode    = %2d\n", maxcliques.numnode);
    for(int i=0; i<maxcliques.clique.size(); i++){
    	printf("clique %2d = ", i);
        for(list<int>::iterator lit=maxcliques.clique[i].begin();lit!=maxcliques.clique[i].end();++lit){
    		printf("%2d ", (*lit));
        }
    	printf("\n");
    }
   	printf("\n");
 */
    //vector<double> cpuTime(10);
    //cpuTime.resize(10,0);
    //cout<<" ***>s3r::genBasicIndices ---> "<<endl;
    
    int nDim=polysys.dimvar();
    int rowSize=polysys.numsys();
    
    if(sparsesw==0){
		BasisIndices.resize(rowSize+1);
		for(int i=0; i<rowSize+1;i++){
			for(int j=0; j<nDim; j++){
				BasisIndices[i].push_back(j);
			}
		}
    }else{
        //cpuTime[0] = (double)clock();
        int nClique=maxcliques.numcliques;
        BasisIndices.resize(rowSize+nClique);
       	bool flag;
        int maxSizeClique=maxcliques.maxnnz();
        int maxSize, nnz;
		list<int> candidates, varList,nzIndicator,dummyNzIndicator,clique;
		list<int>::iterator lit1,lit2;
        //cpuTime[1] = (double)clock();
        for(int i=1;i<rowSize;i++){
			varList.clear();
			candidates.clear();
            polysys.sumSupports(i,varList);
            //double t1,t2,t3,t4;
            //t1 = (double)clock();
            for(int j=0;j<nClique;j++){
                nnz=0;
				lit1 = varList.begin();
				lit2 = varList.end();
                flag = includes(maxcliques.clique[j].begin(),maxcliques.clique[j].end(),lit1,lit2);
				if(flag){
                    candidates.push_back(j);
                }
            }
            //t2 = (double)clock();
            //cpuTime[3] = cpuTime[3] + (t2-t1);
            //t1 = (double)clock();
			lit1 = candidates.begin();
			nzIndicator = maxcliques.clique[(*lit1)];
            maxSize=maxSizeClique*multifactor;
            //expanding nzIndicator until its size does not exceed maxSize
            nnz = nzIndicator.size();
           	dummyNzIndicator = nzIndicator; 
            if(nnz < maxSize){
                advance(lit1,1);
				for(;lit1 != candidates.end();++lit1){
					clique = maxcliques.clique[(*lit1)];
					dummyNzIndicator.merge(clique);
					dummyNzIndicator.unique();
					nnz = dummyNzIndicator.size();
                    if(nnz <= maxSize){
                        nzIndicator = dummyNzIndicator;
                        nnz = nzIndicator.size();
                        if(nnz >= maxSize){
                            break;
                        }
                    }
                    else{
                        dummyNzIndicator = nzIndicator;
                    }
                }
            }
			BasisIndices[i] = nzIndicator;
			//t2 = (double)clock();
            //cpuTime[4] = cpuTime[4] + (t2-t1);
        }
        for(int i=0;i<nClique;i++){
			BasisIndices[i+rowSize] = maxcliques.clique[i];
        }
        //s2 = (double)clock();
        //cpuTime[5] = s2 -s1;
        //cpuTime[2] = (double)clock();
    }
    /*
    printf("cpu time Part1 = %f\n",(cpuTime[1]-cpuTime[0])/(double)CLOCKS_PER_SEC);
    printf("cpu time Part2 = %f\n",(cpuTime[2]-cpuTime[1])/(double)CLOCKS_PER_SEC);
    printf("cpu time for 1 = %f\n",(cpuTime[3])/(double)CLOCKS_PER_SEC);
    printf("cpu time for 2 = %f\n",(cpuTime[4])/(double)CLOCKS_PER_SEC);
    printf("cpu time for 3 = %f\n",(cpuTime[5])/(double)CLOCKS_PER_SEC);
     */
}
void s3r::write_pop(int ell, string fname){
    
    std::ofstream fout;
    fout.open(fname.c_str(),ios::app);
    if(fout.fail()){
        cout << "error@write_pop:file does not open for output" << endl;
        cout << fname;
        exit(1);
    }
    if(ell == 1){
        fout << "# Scaled and modified POP to be solved " << endl;
    }else if(ell == 0){
        fout << "# POP to be solved";
    }else{
        fout << "Error MSG" << endl;
        exit(1);
    }
    int nDim = this->Polysys.dimVar;
    
    for(int k=0;k < this->Polysys.numSys;k++){
        if(k==0){
            fout << "objPoly" << endl;
        }else{
            fout << "ineqPolySys" << endl;
            fout << "Polynomial   " << k << endl;
        }
        fout <<	"typeCone = " << this->Polysys.polyTypeCone(k);
        fout << " sizeCone = " << this->Polysys.polySizeCone(k);
        fout << " dimVar = " << nDim;
        fout << " degree = " << this->Polysys.polyDegree(k);
        fout << " noTerms = " << this->Polysys.polyNoterms(k) << endl;
        fout << "supportSet = " << endl;
       
	int j; 
        list<class mono>::iterator Mono = this->Polysys.polynomial[k].monoList.begin();
        int length = this->Polysys.polynomial[k].monoList.size();
		bool flag;
		for(;Mono!=Polysys.polynomial[k].monoList.end();++Mono){
			for(int i=0; i<(*Mono).nDim; i++){
				j = 0;
				flag = true;
				while(j < (*Mono).supIdx.size()){
					if(i == (*Mono).supIdx[j]){
						fout << (*Mono).supVal[j] << " ";
						flag = false;
						break;
					}else{
						j++;
					}
				}
				if(flag){
						fout << "0 ";
				}
			}
			fout << endl;
		}
        Mono = this->Polysys.polynomial[k].monoList.begin();
        fout << "coefficient = " << endl;
        if(Polysys.polynomial[k].sizeCone != 1){
            for(int i =0; i < length; i++){
                fout << i+1 << ":";
                for(int j = 0; j < (*Mono).lengthCoef(); j++){
                    fout << (*Mono).Coef[j] << " ";
                }
                fout << endl;
                ++Mono;
            }
        }else{
            for(int i =0; i < length; i++){
                fout << i+1 << ":";
                for(int j = 0; j < (*Mono).lengthCoef(); j++){
                    fout << (*Mono).Coef[j] << " ";
                }
                if((i+1) % 10 == 0){
                    fout << endl;
                }
                ++Mono;
            }
        }
        fout << endl;
        fout << endl;
    }
    fout << endl;
    fout << "lbd = " << endl;
    for(int i = 0; i < nDim; i++){
        fout << i << ":";
        if(ell == 1){
            fout << this->Polysys.boundsNew.lbd(i);
        }else if(ell== 0){
            fout << this->Polysys.bounds.lbd(i);
        }
        fout << " ";
        if((i+1) % 10 == 0){
            fout << endl;
        }
    }
    fout << endl;
    fout << "ubd = " << endl;
    for(int i = 0; i < nDim; i++){
        fout << i << ":";
        if(ell == 1){
            fout << this->Polysys.boundsNew.ubd(i);
        }else if(ell== 0){
            fout << this->Polysys.bounds.ubd(i);
        }
        fout << " ";
        if((i+1) % 10 == 0){
            fout << endl;
        }
    }
    fout << endl;
    fout << endl;
    fout.close();
}
void s3r::write_BasisSupports(int i,string fname,class supsetSet bSup){
    std::ofstream fout;
    fout.open(fname.c_str(),ios::app);
    if(fout.fail()){
        cout << "error@write_BasisSupports:file does not open for output" << endl;
        cout << fname;
        exit(1);
    }
    if(i == 1){
        fout << "# basisSupports after reduction" << endl;
    }else if(i == 0){
        fout << "";
    }else{
        fout << "Error MSG" << endl;
        exit(1);
    }
    fout << "# basisSupports --- the support set" << endl;
    fout << "  used for each equality and inequality" << endl;
    
    int mDim;
    class supSet sup;
    mDim = bSup.supsetArray.size();
    for(int i=0;i<mDim;i++){
        sup = bSup.supsetArray[i];
        sup.out_full(i,fname);
    }
    
    fout << endl;
    fout.close();
}

void s3r::write_BasisIndices(string fname){
    std::ofstream fout;
    //fout.open(fname.c_str(),ios::out);
    //fout.close();
    fout.open(fname.c_str(),ios::out|ios::app);
    if(fout.fail()){
        cout << "error@write_BasisIndices:file does not open for output" << endl;
        cout << fname;
        exit(1);
    }
    fout << "# basisIndices --- the set of nonzero coodinates" << endl;
    fout << "  used for each equality and inequality" << endl;
    
    int i,j,k,nDim, mDim;
    nDim = this->Polysys.dimVar;
    mDim = (int)(this->bindices.size());
	list<int>::iterator lit;
    for(i=1;i<mDim;i++){
        fout << i << " : ";
		lit = bindices[i].begin();
		for(;lit != bindices[i].end();++lit){
			j = (*lit);
            fout << j+1 << " ";
		}
        /*
        for(j=0;j<nDim;j++){
			if(this->bindices[i*nDim+j] == 1){
                fout << j+1 << " ";
            }
			
        }
		*/
        fout << endl;
    }
    fout << endl;
    fout.close();
}

//genBasisSupports
void s3r::genBasisSupports(class supsetSet & BasisSupports){
    
    //cout<<" ***>s3r::genBasisSupports ---> "<<endl;
    
    int i,j;
    int rowSize=bindices.size();
    int nDim=this->Polysys.dimvar();
    
    int nVars;
    list<class sup> List;
    class supSet Moment;
    int sosDim;
    BasisSupports.push(Moment);
    for(i=1;i<rowSize;i++){
        nVars=0;
		/*
		for(list<int>::iterator vit =o2n_pattern.begin(); vit != o2n_pattern.end(); ++vit){
			printf("%2d ", (*vit));
		}
		printf("\n");
		*/
		//nVars = o2n_pattern.size();
		nVars = bindices[i].size();
        List.clear();
        if(i<this->Polysys.numsys()){
            
            sosDim=this->param.relax_Order-(int)ceil((double)(this->Polysys.polyDegree(i))/2.0);
            if(this->Polysys.polyTypeCone(i)==EQU){
                //cout<<"Equality constraints "<<i<<" sosDim="<<sosDim<<endl;
                
                genLexAll(nVars,2*sosDim,List);
            }
            else{
                //cout<<"Inequality constraints"<<i<<" sosoDim="<<sosDim<<endl;
                genLexAll(nVars,sosDim,List);
            }
        }
        else{
            //cout<<"Moment matrices "<<i<<" sosDim="<<sosDim<<endl;
            genLexAll(nVars,this->param.relax_Order,List);
        }
        Moment.dimVar = this->Polysys.dimVar;
        Moment.setSupSet(nDim,List);
        Moment.changeIndicesAll(bindices[i]);
        BasisSupports.push(Moment);
    }
	
    //cout<<" <--- s3r::genBasisSupports <*** "<<endl<<endl;
}

//Delete unnecessarry supports by exploiting complimentarity constraints.
void s3r::eraseCompZeroSups(class supSet & czSups,vector<class supSet> & BaSups)
{
    //cout<<" ***> s3r::eraseCompZeroSups ---> "<<endl;
    
    int sizeP=BaSups.size();
    int sizeCZ=czSups.size();
    list<class sup>::iterator czIte=czSups.begin();
    vector<int> czIdx,czVal;
    
    for(int i=0;i<sizeCZ;i++){
        
        czIdx.clear();
        czVal.clear();
        (*czIte).getIdxsVals(czIdx,czVal);
        
        for(int j=1;j<sizeP;j++){
            
            int sizeBA=BaSups[j].size();
            list<class sup>::iterator baIte=BaSups[j].begin();
            vector<int> baIdx,baVal;
            
            for(int k=0;k<sizeBA;k++){
                
                baIdx.clear();
                baVal.clear();
                (*baIte).getIdxsVals(baIdx,baVal);
                
                int sizeCID=czIdx.size();
                int sizeBID=baIdx.size();
                
                if(sizeBID >= sizeCID){
                    int h=0;
                    for(int s=0;s<sizeBID && h<sizeCID;s++){
                        if(baIdx[s] == czIdx[h]){
                            if(baVal[s] >= czVal[h]){
                                h++;
                            }
                        }
                    }
                    if(h==sizeCID){
                        //	cout<<"Erased: czSup--> ";(*czIte).disp();
                        //	cout<<"        baSup--> ";(*baIte).disp();
                        baIte=BaSups[j].erase(baIte);
                        sizeBA--;
                        //	cout<<"        from BaSups["<<j<<"]"<<endl;
                    }
                    
                }
                
                baIte++;
                
            }
        }
        czIte++;
    }
    
    //cout<<" <--- s3r::eraseCompZeroSups <*** "<<endl<<endl;
    
}

void Div2(class sup & sup1){
	for(vector<int>::iterator it=sup1.val.begin();it!=sup1.val.end();++it){
		(*it) >>= 1;
	}
}

void Multi2(class sup & sup1){
	for(vector<int>::iterator it=sup1.val.begin();it!=sup1.val.end();++it){
		(*it) <<= 1;
	}
}

bool IsNonNegative(class sup sup1){
	for(vector<int>::iterator it=sup1.val.begin();it!=sup1.val.end();++it){
		if((*it) < 0){
			return false;
		}	
	}
	return true;
}
bool IsZero(class sup sup1){
	for(vector<int>::iterator it=sup1.val.begin();it!=sup1.val.end();++it){
		if((*it) != 0){
			return false;
		}
	}
	return true;
}

void s3r::reduceSupSets(class supsetSet & BasicSupports,class supSet allNzSups){
    
    //cout<<" ***>s3r::reduceSupSets --->"<<endl;
    
    //Find supports with all even coordinates
    class supSet Fe;
    Fe.setDimVar(Polysys.dimvar());
    allNzSups.getEvenSups(Fe,YES);
    class sup zerosup, Sup, Sup1, Sup2, Sup3;
    Fe.pushSup(zerosup);
    
    //Fe=Fe/2
    for_each(Fe.begin(),Fe.supList.end(), Div2);
    
    int size = Polysys.numsys();
    int ABSsize = BasicSupports.supsetArray.size();
    int Csize = ABSsize-size;
    
    vector<vector<int> > checkList(Csize);
    class supSet bSupSet;
    int bsize, tempsize, glsize;
    
    class supSet2 A(Polysys.dimvar());
    list<class sup>::iterator SupIte, SupRIte, SupWIte;
    
    int oldAsize = 0;
    for(int i=0;i<Csize;i++){
        bSupSet=BasicSupports.supsetArray[i+size];
        SupIte=bSupSet.begin();
        checkList[i].resize(bSupSet.size(),1);
        bsize=bSupSet.size();
    	oldAsize += bsize;    
        for(int j=0;j<bsize;j++){
            if(!Fe.doesExist(*SupIte)){
                A.addSup(i,j,(*SupIte));
            }
            ++SupIte;
        }
    }
    list<class sup2>::iterator SupIte2;
    vector<int> Indices, gl, no;
    int newAsize = oldAsize+1;
	vector<int>::iterator it1,it2,vt1,vt2;   
	bool gonext,reduce;
    while(oldAsize < newAsize){
        newAsize = oldAsize;
        for(SupIte2=A.begin(); SupIte2 != A.end(); ++SupIte2){
            gl.clear();
            no.clear();
            (*SupIte2).RL(gl,no);
            if(gl.empty()==true){
                cout<<"error@reduceSupSets: gl.size is ZERO"<<endl;
                exit(1);
            }
            gonext = true;
            reduce = false;
            glsize = gl.size();
            for(int i=0;i<glsize;i++){
                SupRIte=BasicSupports.supsetArray[gl[i]+size].begin();
                tempsize = BasicSupports.supsetArray[gl[i]+size].size();
                for(int j=0;j<tempsize;j++){
                    if(checkList[gl[i]][j]!=0 && j!=no[i]){
						(*SupIte2).getSup(Sup);
                        Multi2(Sup);
						gonext = includes(Sup.idx.begin(),Sup.idx.end(),(*SupRIte).idx.begin(),(*SupRIte).idx.end());
						if(gonext){
							Sup1 = SupMinusSup(Sup,(*SupRIte));
							gonext = IsNonNegative(Sup1);
						}
                        if(gonext == false){
                            ++SupRIte;
                            continue;
                        }
                        SupWIte  = BasicSupports.supsetArray[gl[i]+size].begin();
                       	advance(SupWIte,j); 
                        for(int k=j;k<tempsize;k++){
                            if(checkList[gl[i]][k]!=0 && k!=no[i]){
								Sup2 = Sup1;
								reduce = includes(Sup1.idx.begin(),Sup1.idx.end(),(*SupWIte).idx.begin(),(*SupWIte).idx.end());
								if(reduce){
									Sup3 = SupMinusSup(Sup2,(*SupWIte));
									reduce = IsZero(Sup3);
								}
								if(reduce){
									break;
								}
                            }
                            ++SupWIte;
                        }
                    }
                    if(reduce){
                        break;
                    }
                    ++SupRIte;
                }
                if(reduce){
                    break;
                }
            }
            if(reduce == false){
                for(int i=0;i<glsize;i++){
                	checkList[gl[i]][no[i]]=0;
                }
            }
        }
        oldAsize = 0;
        for(int i=0;i<Csize;i++){
            for(int j=0;j<checkList[i].size();j++){
                if(checkList[i][j]!=0){
                    oldAsize = oldAsize +1;
                }
            }
        }
    }
    class supSet dumSupSet;
    for(int i=0;i<Csize;i++){
        dumSupSet.clear();
        dumSupSet.setDimVar(Polysys.dimvar());
        SupIte=BasicSupports.supsetArray[i+size].begin();
        for(int j=0;j<checkList[i].size();j++){
            if(checkList[i][j]!=0){
                dumSupSet.addSup(*SupIte);
            }
            ++SupIte;
        }
        BasicSupports.supsetArray[i+size]=dumSupSet;
    }
    //cout<<" <--- s3r::reduceSupSets <*** "<<endl<<endl;
}

void get_subjectto_polys_and_basups(
   /* IN */  class polysystem polysys, vector<list<int> > BaIndices, vector<class supSet> basups, int stsize,
   /* OUT */ vector<class poly_info> & polyinfo_st, vector<class bass_info> & bassinfo_st)
{
    
    int ndim = polysys.dimvar();
    int bdim;
   	list<int>::iterator lit;
    for(int i=1;i<stsize+1;i++){
        initialize_polyinfo(polysys,i,polyinfo_st[i-1]);
		bdim = BaIndices[i].size();
        bassinfo_st[i-1].dim = bdim;
        bassinfo_st[i-1].deg = basups[i].deg();
        bassinfo_st[i-1].alloc_pattern(bdim);
        bdim = 0;
		lit = BaIndices[i].begin();
		for(;lit != BaIndices[i].end(); ++lit){
			bassinfo_st[i-1].pattern[bdim] = (*lit);
			bdim++;
		}
        initialize_spvecs(basups[i].supList,bassinfo_st[i-1].sup);
        polyinfo_st[i-1].no = i;
    }
}
void count_lexall_num_a_nnz(/*IN*/int dimvar,int deg,/*OUT*/int & num,int & nnz){
    
    num = 1;
    nnz = 1;
    
    double dum=1;
    
    if(deg != 0){
        
        for(double k=1;k<=(double)deg;k++){
            dum *= (dimvar+k)/k;
        }
        
        int idum = (int)dum;
        if( dum - (double)dum < 1.0){
            idum += 1;
        }
        
        num = idum;
        
        
        dum = 1;
        if(dimvar < deg - 1){
            for(double k=0; k <= (double)dimvar-1; k++){
                dum *= (k+(double)deg)/(k+1.0);
            }
        }
        else{
            for(double k=1; k <=(double)deg-1; k++){
                dum *= (dimvar+k)/k;
            }
        }
        
        idum = (int)dum;
        if( dum - (double)idum < 1.0 ){
            idum += 1;
        }
        
        nnz = dimvar * idum + 1;
        
    }
    else{
        num = 1;
        nnz = 0;
    }
    
    
}
void get_allsups(int dim,class poly_info & polyinfo_obj, int stsize,vector<class poly_info> polyinfo_st, vector<class bass_info> bassinfo_st,class spvec_array & allsups){
    
    //cout<<" ***> get_allsups ---> "<<endl;
    
    allsups.dim = dim;
    if(stsize > 0){
        vector<class Vec3> InfoTable(stsize);
        for(int i=0; i<stsize; i++){
            InfoTable[i].clear();
            InfoTable[i].resize(3,0);
            InfoTable[i].vec[0] = polyinfo_st[i].typeCone;
            InfoTable[i].vec[1] = bassinfo_st[i].dim;
            InfoTable[i].vec[2] = bassinfo_st[i].deg;
            InfoTable[i].no = i;
        }
        
        
        vector<int> stand(3);
        vector<int> infolist(stsize);
        for(int i=0;i<stsize;i++){
            infolist[i] = i;
        }
        sortInfoTable(InfoTable,infolist);
        
        int nzele=0;
        int numele = 0;
        nzele  += polyinfo_obj.sup.pnz_size;
        numele += polyinfo_obj.sup.vap_size;
        
        int k=0;
        while(k < stsize && InfoTable[infolist[k]].vec[0] == EQU){
            nzele += polyinfo_st[infolist[k]].sup.vap_size * bassinfo_st[k].sup.pnz_size
            + bassinfo_st[k].sup.vap_size * polyinfo_st[infolist[k]].sup.pnz_size;
            numele += polyinfo_st[infolist[k]].sup.pnz_size * bassinfo_st[infolist[k]].sup.pnz_size;
            k++;
        }
        int bpsize, bvsize;
        while(k<stsize){
            count_lexall_num_a_nnz(InfoTable[infolist[k]].vec[1],2*InfoTable[infolist[k]].vec[2],bpsize,bvsize);
            
            nzele += polyinfo_st[infolist[k]].sup.vap_size * bpsize
            + bvsize * polyinfo_st[infolist[k]].sup.pnz_size;
            numele += polyinfo_st[infolist[k]].sup.pnz_size * bpsize;
            
            k++;
        }
        
        allsups.alloc(numele*2,nzele*2);
        
        pushsups(polyinfo_obj.sup,allsups);
        
        k=0;
        class spvec_array minsups;
        
        while(k < stsize && InfoTable[infolist[k]].vec[0] == EQU){
            minkovsum(polyinfo_st[infolist[k]].sup,bassinfo_st[infolist[k]].sup,minsups);
            pushsups(minsups,allsups);
            k++;
        }
        
        
        //class spvec_array lexallsups;
        class spvec_array lexallsups;
        
        bool issame = false;
        vector<int> onpattern(dim);
        while(k<stsize){
            
            if(issame == false){
                genLexAll(InfoTable[infolist[k]].vec[1],2*InfoTable[infolist[k]].vec[2],lexallsups);
                for(int i=0;i<InfoTable[infolist[k]].vec[1];i++){
                    onpattern[i] = bassinfo_st[infolist[k]].pattern[i];
                }
            }
            else{
                for(int i=0;i<InfoTable[infolist[k]].vec[1];i++){
                    onpattern[bassinfo_st[infolist[k-1]].pattern[i]] = bassinfo_st[infolist[k]].pattern[i];
                }
            }
            
            for(int i=0;i<lexallsups.vap_size;i++){
                lexallsups.vap[0][i] = onpattern[lexallsups.vap[0][i]];
            }
            /*
            if(k == 2){
                cout << "onpattern = ";
                for(int i=0; i<dim; i++){
                    cout << onpattern[i] << " ";
                }
                cout << endl;
                cout << "k      = " << k << endl;
                cout << "InfoTable.vec[1] = " << InfoTable[infolist[k]].vec[1] << endl;
                cout << "InfoTable.vec[2] = " << InfoTable[infolist[k]].vec[2] << endl;
                cout << "issame = " << issame << endl;
                lexallsups.disp();
            }
             */
            minkovsum(polyinfo_st[infolist[k]].sup,lexallsups,minsups);
            pushsups(minsups,allsups);
            
            k++;
            
            if(k<stsize){
                if(InfoTable[infolist[k]].vec[1] == InfoTable[infolist[k-1]].vec[1]){
                    
                    if(InfoTable[infolist[k]].vec[2] == InfoTable[infolist[k-1]].vec[2]){
                        issame = true;
                    }else{
                        issame = false;
                        count_lexall_num_a_nnz(InfoTable[infolist[k]].vec[1],InfoTable[infolist[k]].vec[2],lexallsups.pnz_size,lexallsups.vap_size);
                    }
                }
                else{
                    issame = false;
                }
            }
            
        }
        InfoTable.clear();
        lexallsups.del();
        //allsups.disp();
        //simplification(allsups);
    }
    else{
        //case:POP constains no constraints
        //get maximal size and nonzeros of objective function's supports
        int nzele  = polyinfo_obj.sup.pnz_size;
        int numele = polyinfo_obj.sup.vap_size;
        
        //allocate memory
        allsups.alloc(nzele*2,numele*2);
        
        //get all supports of objective function
        pushsups(polyinfo_obj.sup,allsups);
    }
    
    //allsups.clean();
    //cout<<" <--- get_allsups <*** "<<endl;
}

bool comp_InfoTable(class Vec3  vec1, class Vec3 vec2){
    if(vec1.vec[0] == EQU){
        if(vec2.vec[0] != EQU){
            return true;
        }
        return false;
    }else if(vec2.vec[0] == EQU){
        return false;
    }else{
        if(vec1.vec[1] < vec2.vec[1]){
            return true;
        }else if(vec1.vec[1] > vec2.vec[1]){
            return false;
        }else{
            if(vec1.vec[2] > vec2.vec[2]){
                return true;
            }
            return false;
        }
    }
    
}
void sortInfoTable(vector<class Vec3> InfoTable,vector<int> & infolist){
    vector<class Vec3> temp(InfoTable.size());
    copy(InfoTable.begin(),InfoTable.end(),temp.begin());
    sort(temp.begin(), temp.end(),comp_InfoTable);
    for(int i=0; i<infolist.size(); i++){
        infolist[i] = temp[i].no;
    }
    temp.clear();
}


void genLexFixDeg(int k,int n,int W,vector<vector<int> > sup,int nnz,class spvec_array & rsups){
    int d;
    for(int i=W;i>0;i--){
        //Sup.push(k,i);
        sup[0][nnz] = k;
        sup[1][nnz] = i;
        nnz ++ ;
        if(W-i>0){
            for(int j=k+1;j<n;j++){
                genLexFixDeg(j,n,W-i,sup,nnz,rsups);
            }
        }
        else{
            rsups.pnz[0][rsups.pnz_size] = rsups.vap_size;
            rsups.pnz[1][rsups.pnz_size] = nnz;
            rsups.pnz_size ++ ;
            d = 0;
            while(d < nnz){
                rsups.vap[0][rsups.vap_size] = sup[0][d];
                rsups.vap[1][rsups.vap_size] = sup[1][d];
                rsups.vap_size ++ ;
                d++;
            }
            if(k==n-1){ 
				break; 
			}
        }
        nnz -- ;
    }
}
void genLexAll(int totalOfVars,int Deg,class spvec_array & rsups){
    
    int nnz = 0;
    
    vector<vector<int> > sup(2);
    sup[0].clear();
    sup[1].clear();
    sup[0].resize(totalOfVars,0);
    sup[1].resize(totalOfVars,0);
    
    int psize,vsize;
    count_lexall_num_a_nnz(totalOfVars,Deg,psize,vsize);
    rsups.alloc(psize,vsize);
    
    rsups.pnz[0][0] = -1;
    rsups.pnz[1][0] =  0;
    
    rsups.pnz_size = 1;
    rsups.vap_size = 0;
    
    if(Deg > 1){
        for(int W=1;W<=Deg;W++){
            for(int k=0;k<totalOfVars;k++){
                genLexFixDeg(k,totalOfVars,W,sup,nnz,rsups);
            }
        }
    }
    //rsups.disp();
}

void get_removesups(class polysystem & polysys,class spvec_array & removesups){
    
    //cout<<" ***> get_removesups ---> "<<endl;
    
    list<class sup> czlist;
    
    int size=polysys.numsys();
    
    for(int i=1;i<size;i++){
        class sup czSup;
        double Coef;
        /*
        cout << "i = " << i << endl;
        cout << " val = " << polysys.polyIsComplimentarity(i,czSup,Coef) << endl;
         */
        if(polysys.polyIsComplimentarity(i,czSup,Coef)==YES){
            //cout << "Coef = " << Coef << endl;
            if(fabs(Coef) < EPS){
                czlist.push_back(czSup);
                //czSup.disp();
            }
        }
    }
    //cout << "End of for loop " << endl;
    
    removesups.dim = polysys.dimVar;
    initialize_spvecs(czlist,removesups);
    simplification(removesups);
    
    //cout<<" <--- get_removesups <*** "<<endl;
}
void get_allsups_in_momentmatrix(int dimvar,int mmsize, vector<class bass_info> bassinfo_mm,class spvec_array & mmsups){
    
    mmsups.dim = dimvar;
    
    int num = 0;
    int nnz = 0;
    for(int i=0;i<mmsize;i++){;
	    num += (bassinfo_mm[i].sup.pnz_size * ( bassinfo_mm[i].sup.pnz_size + 1) ) << 1;
    	nnz += bassinfo_mm[i].sup.pnz_size * bassinfo_mm[i].sup.get_nnz();
    }
    
    mmsups.alloc(num,nnz);
    mmsups.pnz_size = 0;
    mmsups.vap_size = 0;
    
    class spvec_array lexallsups;
    for(int i=0;i<mmsize;i++){
        minkovsum(bassinfo_mm[i].sup,lexallsups);
        pushsups(lexallsups,mmsups);
    }
    lexallsups.del();
    simplification(mmsups);
}
void get_momentmatrix_basups(class polysystem polysys,vector<list<int> > BaIndices, vector<class supSet> & basups, vector<class bass_info> & bassinfo_mm)
{
    int bdim;
    int ndim = polysys.dimvar();
   	list<int>::iterator lit; 
    int mmsize = basups.size() - polysys.numsys();
    for(int i=0;i<mmsize;i++){
		bdim = BaIndices[i].size();
        bassinfo_mm[i].dim = bdim;
        bassinfo_mm[i].deg = basups[i+polysys.numsys()].deg();
        bassinfo_mm[i].alloc_pattern(bdim);
        
        bdim = 0;
        lit = BaIndices[i].begin();
		for(;lit != BaIndices[i].end();++lit){
			bassinfo_mm[i].pattern[bdim] = (*lit);
			bdim++;
		}
        initialize_spvecs(basups[i+polysys.numsys()],bassinfo_mm[i].sup);
    	//bassinfo_mm[i].sup.disp();
	}
}

void get_poly_a_bass_info(
   /* IN */  class polysystem & polysys, vector<class supSet> & BaSupVect,vector<class supSet> & mmBaSupVect,
   const int mat_size,
   /* OUT */ vector<class poly_info> & polyinfo,vector<class spvec_array> & bassinfo){
       
       //cout<<" ***> set_mat_info_data ---> "<<endl;
       
       int no_poly = 0;
       
       initialize_polyinfo(polysys,0,polyinfo[no_poly]);
       //	cout << "0 no_poly = " << no_poly << endl;
       no_poly++;
       for(int i=1;i<polysys.numsys();i++){
           if(polysys.polyTypeCone(i)==EQU){
               initialize_polyinfo(polysys,i,polyinfo[no_poly]);
               initialize_spvecs(BaSupVect[i],bassinfo[no_poly]);
               no_poly++;
           }
       }
       //cout << "1 no_poly = " << no_poly << endl;
       for(int i=1;i<polysys.numsys();i++){
           if(polysys.polyTypeCone(i)==INE && BaSupVect[i].size()==1){
               //cout << "***0***" << endl;
               initialize_polyinfo(polysys,i,polyinfo[no_poly]);
               
               //cout << "***1***" << endl;
               initialize_spvecs(BaSupVect[i],bassinfo[no_poly]);
               
               //cout << "***2***" << endl;
               no_poly++;
           }
       }
       
       //cout << "2 no_poly = " << no_poly << endl;
       int mmsize=mmBaSupVect.size();
       for(int i=0;i<mmsize;i++){
           if(mmBaSupVect[i].size() == 1){
               polyinfo[no_poly].typeCone = -999;
               initialize_spvecs(mmBaSupVect[i],bassinfo[no_poly]);
               no_poly++;
           }
       }
       
       for(int i=0;i<mmsize;i++){
           if(mmBaSupVect[i].size() > 1){
               polyinfo[no_poly].typeCone = -999;
               initialize_spvecs(mmBaSupVect[i],bassinfo[no_poly]);
               no_poly++;
           }
       }
       //cout << "4 no_poly = " << no_poly << endl;
       
       for(int i=1;i<polysys.numsys();i++){
           if(polysys.polyTypeCone(i)==INE && BaSupVect[i].size()>1){
               initialize_polyinfo(polysys,i,polyinfo[no_poly]);
               initialize_spvecs(BaSupVect[i],bassinfo[no_poly]);
               no_poly++;
           }
       }
       //cout << "5 no_poly = " << no_poly << endl;
       
       for(int i=1;i<polysys.numsys();i++){
           if(polysys.polyTypeCone(i)==SDP){
               initialize_polyinfo(polysys,i,polyinfo[no_poly]);
               initialize_spvecs(BaSupVect[i],bassinfo[no_poly]);
               no_poly++;
           }
       }
       //cout<<" <--- set_mat_info_data <*** "<<endl<<endl;
       
}
void initialize_supset(class spvec_array & spvecs,class supSet & supset){
    
    supset.supList.clear();
    supset.dimVar = spvecs.dim;
    
    for(int i=0;i<spvecs.pnz_size;i++){
        class sup Sup;
        int a = spvecs.pnz[0][i];
        //cout << " i = " << i << endl;
        for(int j=0;j<spvecs.pnz[1][i];j++){
            Sup.idx.push_back(spvecs.vap[0][a]);
            Sup.val.push_back(spvecs.vap[1][a]);
            a++;
        }
        //cout << " i = " << i << ", spvecs.pnz_size = " << spvecs.pnz_size << endl;
        supset.supList.push_back(Sup);
    }
}
void initialize_spvecs(/*IN*/class supSet & supset,/*OUT*/class spvec_array & spvecs){
    
    int size1 = supset.size();
    int size2 = supset.nnz();
    spvecs.dim  = supset.dimVar;
    spvecs.alloc(size1,size2);
    spvecs.pnz_size = size1;
    list<int> vars;
    
    list<class sup>::iterator ite = supset.begin();
    spvecs.vap_size=0;
    vector<int> Idx,Val;
    for(int i=0;i<size1;i++){
        if((*ite).nnz() ==0){
            spvecs.pnz[0][i] = -1;
            spvecs.pnz[1][i] =  0;
        }
        else{//if( (*ite).nnz > 0 )
            spvecs.pnz[0][i] =  spvecs.vap_size;
            spvecs.pnz[1][i] = (*ite).nnz();
            
            Idx.clear();
            Val.clear();
            
            (*ite).getIdxsVals(Idx,Val);
            
            for(int j=0; j<spvecs.pnz[1][i]; j++){
                spvecs.vap[0][spvecs.vap_size] = Idx[j];
                spvecs.vap[1][spvecs.vap_size] = Val[j];
                spvecs.vap_size++;
                vars.push_back(Idx[j]);
            }
        }
        
        ++ite;
    }
    vars.unique();
    spvecs.dim2 = vars.size();
    
}
void initialize_spvecs(/*IN*/list<class sup> & suplist,/*OUT*/class spvec_array & spvecs){
    
    int size1 = suplist.size();
    int size2=0;
    
    list<class sup>::iterator ite = suplist.begin();
    for(int i=0;i<size1;i++){
        size2 += (*ite).nnz();
        ++ite;
    }
    
    spvecs.alloc(size1,size2);
    spvecs.pnz_size = size1;
    
    ite = suplist.begin();
    spvecs.vap_size=0;
    vector<int> Idx,Val;
    for(int i=0;i<size1;i++){
        if((*ite).nnz() ==0){
            spvecs.pnz[0][i] = -1;
            spvecs.pnz[1][i] =  0;
        }
        else{//if( (*ite).nnz > 0 )
            spvecs.pnz[0][i] =  spvecs.vap_size;
            spvecs.pnz[1][i] = (*ite).nnz();
            
            Idx.clear();
            Val.clear();
            
            (*ite).getIdxsVals(Idx,Val);
            
            for(int j=0; j<spvecs.pnz[1][i]; j++){
                spvecs.vap[0][spvecs.vap_size] = Idx[j];
                spvecs.vap[1][spvecs.vap_size] = Val[j];
                spvecs.vap_size++;
            }
        }
        
        ++ite;
    }
    //spvecs.disp();
    
}
void initialize_polyinfo(/*IN*/class polysystem & polysys,int nop,/*OUT*/class poly_info & polyinfo){
    
    //set typeCone, sizeCone and number of Monomials.
    polyinfo.typeCone = polysys.polynomial[nop].typeCone;
    polyinfo.sizeCone = polysys.polynomial[nop].sizeCone;
    polyinfo.numMs = polysys.polynomial[nop].monoList.size();
    polyinfo.no = nop;
    
    //set data of all supports
    list<class sup> suplist;
    polysys.polynomial[nop].pushSupList(suplist);
    initialize_spvecs(suplist,polyinfo.sup);
    
    //alloc memory to set all nonzeros of the coefficients
    if(polyinfo.typeCone != SDP){
        polyinfo.alloc_coef(polyinfo.typeCone,polyinfo.sizeCone,polyinfo.numMs,0);
        //cout << "size of polyinfo.coef = " << polyinfo.coef.size() << endl;
    }else{
        polyinfo.alloc_coef(polyinfo.typeCone,polyinfo.sizeCone,polyinfo.numMs,polysys.polyCoefNnz(nop));
    }
    
    
    //set data of all coefficients
    list<class mono>::iterator ite = polysys.polynomial[nop].monoList.begin();
    int nummonos = polysys.polynomial[nop].monoList.size();
    vector<double> Co;
    int lencoef;
    
    if(polyinfo.typeCone != SDP){
        Co.resize(polyinfo.sizeCone,0);
        lencoef = polyinfo.sizeCone;
        //cout << "nummonos = " << nummonos << endl;
        for(int i=0;i<nummonos;i++){
            Co.clear();
            (*ite).copyCoef(Co);
            for(int j=0;j<lencoef;j++){
                polyinfo.coef[i][j] = Co[j];
            }
            ++ite;
        }
    }
    else{
        
        int idx = 0;
        int stidx,edidx;
        int matsize;
        
        for(int i=0;i<nummonos;i++){
            
            Co.clear();
            (*ite).copyCoef(Co);
            
            matsize = i*polyinfo.sizeCone;
            
            for(int j=0;j<polyinfo.sizeCone;j++){
                
                polyinfo.mc[matsize+j] = idx;
                
                stidx = j*polyinfo.sizeCone;
                edidx = stidx + j + 1;
                
                for(int k= stidx; k < edidx ; k ++){
                    if(Co[k] != 0.0){
                        polyinfo.mr  [idx]    = k - stidx;
                        polyinfo.coef[idx][0] = Co[k];
                        idx ++ ;
                    }
                }
            }
            ++ite;
        }
        polyinfo.mc[nummonos*polyinfo.sizeCone] = idx;
        
    }
}

void rescale_sol(int dimvar,vector<double> & pMat,vector<double> & bVec,double * & sol){
    
    //cout<<" ***> rescale_sol ---> "<<endl;
    vector<double> yyy(dimvar);
    for(int i=0; i<dimvar; i++){
        yyy[i] = sol[i];
    }
    
    //x=A(y*)+b
    for(int i=0;i<dimvar;i++){
        sol[i]=0.0;
        for(int j=0;j<dimvar;j++){
            if(pMat[i*dimvar+j] != 0.0 ){
                sol[i] += yyy[j]*pMat[i*dimvar+j];
            }
        }
        if(pMat[i*dimvar+i] != 0.0){
            sol[i] -= bVec[i];///pMat[i*dimvar+i];
        }
    }
    //cout<<" <--- rescale_sol <*** "<<endl<<endl;
    
}
pop_params::pop_params(){
    relax_Order		= RELAXORDER;
    sparseSW		= SPARSESW;
    multiCliquesFactor	= MULTICLIFACT;
    scalingSW		= SCALINGSW;
    boundSW			= BOUNDSW;
    eqTolerance		= EQTOLERANCE;
    perturbation		= PERTURB;
    reduceMomentMatSW       = REDUCEMOMENTSW;
    complementaritySW	= COMPSW;
    //SeDuMiSW		= SEDUMISW;
    //SeDuMiOnScreen		= SEDUMIONSCREEN;
    printOnScreen = PRINTONSCREEN;//
    printLevel1 = PRINTLEVEL1;
    printLevel2 = PRINTLEVEL2;
    SeDuMiEpsilon = SEDUMIEPSILON;
    symbolicMath = SYMBOLICMATH;
    mex = MEX;
    SDPASW = 1;
    SDPAOutFile = "";
    SDPAOnScreen = 1;
    Method = "";
}
void pop_params::write_parameters(string fname){
    std::ofstream fout;
    fout.open(fname.c_str(),ios::out|ios::app);
    if(fout.fail()){
        cout << "error@write_parameters:file does not open for output" << endl;
        cout << fname;
        exit(1);
    }
    fout << "# parameters:" << endl;
    fout << "  relaxOrder         = " << this->relax_Order << endl;
    fout << "  sparseSW           = " << this->sparseSW << endl;
    fout << "  multiCliquesFactor = " << this->multiCliquesFactor << endl;
    fout << "  scalingSW          = " << this->scalingSW << endl;
    fout << "  boundSW            = " << this->boundSW << endl;
    fout << "  eqTolerance        = " << this->eqTolerance << endl;
    fout << "  perturbation       = " << this->perturbation << endl;
    fout << "  reduceMomentMatSW  = " << this->reduceMomentMatSW << endl;
    fout << "  complementaritySW  = " << this->complementaritySW << endl;
    fout << "  SeDuMiSW           = " << this->SeDuMiSW << endl;
    if(this->SeDuMiOnScreen == 1){
        fout << "  SeDuMiOutFile      = " << 1 << endl;
    }else if(this->SeDuMiOnScreen == 0 && this->SeDuMiOutFile.empty() == false){
        fout << "  SeDuMiOutFile      = " << this->SeDuMiOutFile << endl;
    }else if(this->SeDuMiOnScreen == 0 && this->SeDuMiOutFile.empty() == true){
        fout << "  SeDuMiOutFile      = 0" << endl;
    }
    if(this->detailedInfFile.empty() == false){
        fout << "  detailedInfFile    = " << this->detailedInfFile << endl;
    }else{
        fout << "  detailedInfFile    = " << endl;
    }
    if(this->sdpaDataFile.empty() == false){
        fout << "  sdpaDataFile       = " << this->sdpaDataFile << endl;
    }else{
        fout << "  sdpaDataFile       = " << endl;
    }
    if(this->printOnScreen == 1){
        fout << "  printFileName      = " << 1 << endl;
    }else if(this->printOnScreen == 0 && this->printFileName.empty() == false){
        fout << "  printFileName      = " << this->printFileName << endl;
    }else if(this->printOnScreen == 0 && this->printFileName.empty() == true){
        fout << "  printFileName      = 0" << endl;
    }
    fout << "  printLevel         = [" << printLevel1 << ", " << printLevel2 << "]" << endl;
    fout << "  SeDuMiEpsilon      = " << SeDuMiEpsilon << endl;
    fout << "  symbolicMath       = " << symbolicMath << endl;
    fout << "  mex                = " << mex << endl;
    fout << endl;
}
bool input_params(class pop_params & params){
    
    string paramfile("sppop.param");
    ifstream fin(paramfile.c_str());
    
    if(fin.is_open()){
        
        double ddum;
        fin >> ddum; params.relax_Order = (int)ddum;
        fin >> ddum; params.sparseSW = (int)ddum;
        fin >> ddum; params.multiCliquesFactor = (int)ddum;
        fin >> ddum; params.scalingSW = (int)ddum;
        fin >> ddum; params.boundSW = (int)ddum;
        fin >> ddum; params.eqTolerance = ddum;
        fin >> ddum; params.perturbation = ddum;
        fin >> ddum; params.reduceMomentMatSW = (int)ddum;
        fin >> ddum; params.complementaritySW = (int)ddum;
        fin >> ddum; params.SeDuMiSW = (int)ddum;
        fin >> ddum; params.SeDuMiOnScreen = (int)ddum;
        //fin >> ddum; params.SeDuMiOutFile = (string)ddum;
        //fin >> ddum; params.detailedInfFile = (string)ddum;
        //fin >> ddum; params.sdpaDataFile = (string)ddum;
        fin >> ddum; params.printOnScreen = (int)ddum;
        //fin >> ddum; params.printFileName = (string)ddum;
        fin >> ddum; params.printLevel1 = (int)ddum;
        fin >> ddum; params.printLevel2 = (int)ddum;
        fin >> ddum; params.SeDuMiEpsilon = ddum;
        fin >> ddum; params.symbolicMath = (int)ddum;
        fin >> ddum; params.mex = (int)ddum;
        
        fin.close();
        
        return true;
    }
    else{
        return false;
    }
    
}
void s3r::disp_params(){
    
    cout<<"# Parameters: "<<endl;
    cout << "# parameters:" << endl;
    cout << "  relaxOrder         = " << param.relax_Order << endl;
    cout << "  sparseSW           = " << param.sparseSW << endl;
    cout << "  multiCliquesFactor = " << param.multiCliquesFactor << endl;
    cout << "  scalingSW          = " << param.scalingSW << endl;
    cout << "  boundSW            = " << param.boundSW << endl;
    cout << "  eqTolerance        = " << param.eqTolerance << endl;
    cout << "  perturbation       = " << param.perturbation << endl;
    cout << "  reduceMomentMatSW  = " << param.reduceMomentMatSW << endl;
    cout << "  complementaritySW  = " << param.complementaritySW << endl;
    cout << "  SeDuMiSW           = " << param.SeDuMiSW << endl;
    if(param.SeDuMiOnScreen == 1){
        cout << "  SeDuMiOutFile      = " << 1 << endl;
    }else if(param.SeDuMiOnScreen == 0 && param.SeDuMiOutFile.empty() == false){
        cout << "  SeDuMiOutFile      = " << param.SeDuMiOutFile << endl;
    }else if(param.SeDuMiOnScreen == 0 && param.SeDuMiOutFile.empty() == true){
        cout << "  SeDuMiOutFile      = 0" << endl;
    }
    if(param.detailedInfFile.empty() == false){
        cout << "  detailedInfFile    = " << param.detailedInfFile << endl;
    }else{
        cout << "  detailedInfFile    = " << endl;
    }
    if(param.sdpaDataFile.empty() == false){
        cout << "  sdpaDataFile       = " << param.sdpaDataFile << endl;
    }else{
        cout << "  sdpaDataFile       = " << endl;
    }
    if(param.printOnScreen == 1){
        cout << "  printFileName      = " << 1 << endl;
    }else if(param.printOnScreen == 0 && param.printFileName.empty() == false){
        cout << "  printFileName      = " << param.printFileName << endl;
    }else if(param.printOnScreen == 0 && param.printFileName.empty() == true){
        cout << "  printFileName      = 0" << endl;
    }
    cout << "  printLevel = [" << param.printLevel1 << ", " << param.printLevel2 << "]" << endl;
    cout << "  SeDuMiEpsilon = " << param.SeDuMiEpsilon << endl;
    cout << "  symbolicMath     = " << param.symbolicMath << endl;
    cout << "  mex              = " << param.mex << endl;
    cout<<endl;
    
}

void write_sdpa(/*IN*/class mysdp & psdp,/*OUT*/ string sdpafile){
    
    //cout<<"[Write SDPAsparseformat data]"<<endl;
    
    //cout<<" ***> write_sdpa --->"<<endl;
    /*
    std::ofstream fout;
    fout.open(sdpafile.c_str(), ios::out|ios::out );
    if( fout.fail() ){
        cout << "error@write_sdpa:file not open for output" << endl;
        cout << sdpafile;
        exit(1);
    }
     */
    FILE *fp;
    fp = fopen(sdpafile.c_str(),"w+");
    if(fp == NULL){
        printf("file open error\n");
        exit(1);
    }
    
    //psdp.disp();
    //int size = psdp.block_info[0][psdp.nBlocks-1]+psdp.block_info[1][psdp.nBlocks-1]-1;
    //int size = psdp.block_info[1][psdp.nBlocks-1];
    //int size = accumulate(psdp.block_info[1].begin()+1,psdp.block_info[1].end()-1,0);
    int size = 0;
    for(int i=1; i<psdp.nBlocks+1;i++){
        size = size + psdp.block_info[1][i];
    }
    
    
    /*
    cout << "# Output : SDPA sparse format data"  << endl;
    cout << "  File name = "  << sdpafile  << endl;
    cout << "  mDim = " << psdp.mDim << " nBlock = " << psdp.nBlocks << endl;
    cout << "  size of bVect = 1 * "  << psdp.mDim << endl;
    cout << "  size of sparseMatrix = " << size << " * " << 5 << endl;
    fout << "* SDPA sparse format data"  << endl;
    fout << "* File name = "  << sdpafile  << endl;
    fout << "* mDim = " << psdp.mDim << " nBlock = " << psdp.nBlocks << endl;
    fout << "* size of bVect = 1 * "  << psdp.mDim  << endl;
    fout << "* size of sparseMatrix = " << size << " * " << 5 << endl;
     */
    /*
    mexPrintf("* SDPA sparse format data\n");
    mexPrintf("* File name = %s\n", sdpafile.c_str());
    mexPrintf("* mDim = %3d, nBlock = %2d\n",psdp.mDim,psdp.nBlocks);
    mexPrintf("* size of bVect = 1 * %3d\n",psdp.mDim);
    mexPrintf("* size of sparseMatrix = %4d * 5\n",size);
     */
    printf("* SDPA sparse format data\n");
    printf("* File name = %s\n", sdpafile.c_str());
    printf("* mDim = %3d, nBlock = %2d\n",psdp.mDim,psdp.nBlocks);
    printf("* size of bVect = 1 * %3d\n",psdp.mDim);
    printf("* size of sparseMatrix = %4d * 5\n",size);
    
    fprintf(fp,"* SDPA sparse format data\n");
    fprintf(fp,"* File name = %s\n", sdpafile.c_str());
    fprintf(fp,"* mDim = %3d, nBlock = %2d\n",psdp.mDim,psdp.nBlocks);
    fprintf(fp,"* size of bVect = 1 * %3d\n",psdp.mDim);
    fprintf(fp,"* size of sparseMatrix = %4d * 5\n",size);
    /*
    fout<<psdp.mDim<<endl;
    fout<<psdp.nBlocks<<endl;
     */
    fprintf(fp,"%3d\n",psdp.mDim);
    fprintf(fp,"%3d\n",psdp.nBlocks);
    
    for(int i=1;i<psdp.nBlocks+1;i++){
        //fout<<" "<<psdp.bLOCKsTruct[i];
        fprintf(fp,"%2d ",psdp.bLOCKsTruct[i]);
    }
    fprintf(fp,"\n");
    //fout<<endl;
    
    vector<double> obj_coef(psdp.mDim+1);
    obj_coef.clear();
    obj_coef.resize(psdp.mDim+1,0);
    
    for(int j=psdp.block_info[0][0];j<psdp.block_info[0][0]+psdp.block_info[1][0];j++){
        obj_coef[psdp.ele.sup.pnz[0][j]-1] = psdp.ele.coef[j];
    }
    
    for(int i=0;i<psdp.mDim;i++){
        //fout<<obj_coef[i]<<" ";
        fprintf(fp,"%15.10f ",obj_coef[i]);
    }
    fprintf(fp,"\n");
    //fout<<endl;
    
    for(int i=1;i<psdp.nBlocks+1;i++){
        for(int j=psdp.block_info[0][i];j<psdp.block_info[0][i]+psdp.block_info[1][i];j++){
    /*
            fout<<psdp.ele.sup.pnz[0][j]<<" ";
            fout<<psdp.ele.bij[0][j]<<" ";
            fout<<psdp.ele.bij[1][j]<<" ";
            fout<<psdp.ele.bij[2][j]<<" ";
            fout<<psdp.ele.coef[j]<<" ";
            fout<<endl;
     */
            
            fprintf(fp,"%3d ",psdp.ele.sup.pnz[0][j]);
            fprintf(fp,"%3d ",psdp.ele.bij[0][j]);
            fprintf(fp,"%3d ",psdp.ele.bij[1][j]);
            fprintf(fp,"%3d ",psdp.ele.bij[2][j]);
            fprintf(fp,"%15.10f\n",psdp.ele.coef[j]);
        }
    }
    fclose(fp);
    //fout.close();
    
    //cout<<"OK"<<endl;
    
    //cout<<" <--- write_sdpa <***"<<endl<<endl;
}

void perturb_objective(class poly & objpoly,int dimvar,double eps){
    class spvec_array sups;
    //incomplete function!!
}
//class s3r's constructor
s3r :: s3r(){
    this->timedata1.clear();
    this->timedata1.resize(10,0);
    this->timedata.clear();
    this->timedata.resize(21,0);
}

void conversion_part1(
    /*IN*/  class s3r & sr,
    /*OUT*/ double & objconst,
    int & slen, vector<double> & scalevalue,
    int & blen, vector<double> & bvect,
    int & mlen, vector<double> & permmatrix)
{
    if(sr.param.detailedInfFile.empty() == false){
        sr.param.write_parameters(sr.param.detailedInfFile);
    }
    sr.timedata1[0] = (double)clock();
    
    if(sr.param.detailedInfFile.empty() == false){
        sr.write_pop(0,sr.param.detailedInfFile);
    }
    sr.timedata1[1] = (double)clock();
	blen = sr.Polysys.dimvar();
	mlen = sr.Polysys.dimvar();
    
    sr.timedata1[2] = (double)clock();
    //sr.Polysys.writePolynomials();
    if(sr.param.scalingSW == 1){
        sr.Polysys.itemp = 222;
        //double t1 = (double)clock();
        sr.Polysys.scalingPOP(permmatrix,bvect,scalevalue);
	slen = scalevalue.size();
        //double t2 = (double)clock();
        //printf("Scaling Time = %f\n",(t2-t1)/(double)CLOCKS_PER_SEC );
        /*
        for(int i=0; i<sr.Polysys.dimvar(); i++){
            cout <<"lbd(" << i << ")   =" << sr.Polysys.bounds.lbd(i) << endl;
            cout <<"ubd(" << i << ")   =" << sr.Polysys.bounds.ubd(i) << endl;
            cout <<"Newlbd(" << i << ")=" << sr.Polysys.boundsNew.lbd(i) << endl;
            cout <<"Newubd(" << i << ")=" << sr.Polysys.boundsNew.ubd(i) << endl;
        }
         */
    }
    else{
        //set boudsNew.lbd and ubd bound.lbd and and ubd.
        //sr.Polysys.boundsNew.allocUpLo(sr.Polysys.dimvar());
        //cout << sr.Polysys.dimvar() << endl;
        int NumOfActiveBounds = 0;
	for(int i=0;i<sr.Polysys.dimvar();i++){
            sr.Polysys.boundsNew.setUp (i+1,sr.Polysys.bounds.ubd(i));
            sr.Polysys.boundsNew.setLow(i+1,sr.Polysys.bounds.lbd(i));
		if(sr.Polysys.bounds.ubd(i) > MIN && sr.Polysys.bounds.ubd(i) < MAX){
			NumOfActiveBounds++;
		}
		if(sr.Polysys.bounds.lbd(i) > MIN && sr.Polysys.bounds.lbd(i) < MAX){
			NumOfActiveBounds++;
		}
            /*
            cout <<"lbd(" << i << ")   =" << sr.Polysys.bounds.lbd(i) << endl;
            cout <<"ubd(" << i << ")   =" << sr.Polysys.bounds.ubd(i) << endl;
            cout <<"Newlbd(" << i << ")=" << sr.Polysys.boundsNew.lbd(i) << endl;
            cout <<"Newubd(" << i << ")=" << sr.Polysys.boundsNew.ubd(i) << endl;
             */
        }
    	//add upper and lower bounds as constraints to constraints system
	int numSys = sr.Polysys.numSys;	
	vector<poly> tmpPolys(numSys);
	for(int i=0; i< numSys; i++){
		tmpPolys[i] = sr.Polysys.polynomial[i];
	}	
	sr.Polysys.polynomial.resize(numSys+NumOfActiveBounds);
	for(int i=0; i< numSys; i++){
		sr.Polysys.polynomial[i] = tmpPolys[i];
		tmpPolys[i].clear();
	}	
	tmpPolys.clear();
	sr.Polysys.boundToIneqPolySys();
    	//eliminate constant from objective function
        sr.Polysys.layawayObjConst();
        //set permutation matrix and constant vector
        permmatrix.resize(mlen,0);
        for(int j=0;j<sr.Polysys.dimvar();j++){
            permmatrix[j]=1.0;
        }
        bvect.resize(blen,0);
	slen = sr.Polysys.numSys;
        scalevalue.resize(slen,1);
        for(int j=0;j<sr.Polysys.dimvar();j++){
            scalevalue[j]=1.0;
        }
    }
	
    sr.timedata1[3] = (double)clock();
    //cout << "***Start of after scaling" << endl;
    //system("top -b -n 1 | grep MATLAB | head -1 |awk '{printf(\"memory = %s\\n\"), $6}' ");
    
    //get data of objective function' constant
    objconst = sr.Polysys.objConst;
    sr.timedata1[4] = (double)clock();
    
    //perturbate objective function.
    if( sr.param.perturbation > 1.0E-12 ){
        sr.Polysys.perturbObjPoly(3201,sr.param.perturbation);
    }
    sr.timedata1[5] = (double)clock();
    
    if(fabs(sr.param.eqTolerance) > EPS){
        sr.Polysys.relax1EqTo2Ineqs(sr.param.eqTolerance);
    }
    sr.timedata1[6] = (double)clock();
    
    if(sr.param.detailedInfFile.empty() == false){
        if(sr.param.scalingSW == 1 || abs(sr.param.eqTolerance) > EPS || abs(sr.param.perturbation) > EPS){
            sr.write_pop(1,sr.param.detailedInfFile);
        }
    }
    sr.timedata1[7] = (double)clock();
    //system("top -b -n 1 | grep MATLAB | head -1 |awk '{printf(\"memory = %s\\n\"), $6}' ");
    //cout << "***End of conversion_part1" << endl;
  	sr.timedata1[8] = sr.timedata1[7];  
  	sr.timedata1[9] = sr.timedata1[7];  
    //sr.Polysys.writePolynomials();
}

void conversion_part2(
    /*IN*/  class s3r & sr,
    vector<int> oriidx,
    class SparseMat extofcsp,
    /*OUT*/ class mysdp & sdpdata)
{
    int stsize;
    int mmsize;
    int mmsetSize;
    
    vector<class supSet> mmBaSupVect;
    
    class supsetSet BasisSupports;
    class poly_info polyinfo_obj;
    class spvec_array allsups_st;
    class spvec_array removesups;
    class supSet czSups;
    class spvec_array mmsups;
    class spvec_array allsups;
    
    sdpdata.dtime = (double)clock();
    sr.timedata[0] = (double)clock();
    //cout << "0 " << sr.timedata[0] << endl;
    
    //generate max cliques
    if(sr.param.sparseSW == 0){
		sr.maxcliques.initialize(sr.Polysys.dimVar,1);
		for(int i=0;i<sr.Polysys.dimVar ;i++){
			sr.maxcliques.clique[0].push_back(i);
		}
	}else if(sr.param.sparseSW == 1){
		gen_maxcliques3(sr.Polysys.dimVar,oriidx,extofcsp,sr.maxcliques);
    }else{
		printf("sparseSW should be 0 or 1.\n");
	}
	if(sr.param.detailedInfFile.empty() == false){
        sr.maxcliques.write_maxCliques(sr.param.detailedInfFile);
    }
    sr.timedata[1] = (double)clock();
    //cout << "1 " << sr.timedata[1] << endl;
    
    //generate basis indices.
    gen_basisindices(sr.param.sparseSW,sr.param.multiCliquesFactor,sr.Polysys,sr.maxcliques,sr.bindices);
    
	if(sr.param.detailedInfFile.empty() == false){
        sr.write_BasisIndices(sr.param.detailedInfFile);
    }
    sr.timedata[2] = (double)clock();
    //cout << "2 " << sr.timedata[2] << endl;
    
    //generate sets of basis supports
    sr.genBasisSupports(BasisSupports);
    if(sr.param.detailedInfFile.empty() == false){
        sr.write_BasisSupports(0,sr.param.detailedInfFile,BasisSupports);
    }
    
    //sr.dumsups = BasisSupports.supsetArray[1];
    sr.timedata[3] = (double)clock();
    //cout << "3 " << sr.timedata[3] << endl;
    
    //get polyinfo_obj( array data-type to have polynomial form data )
    initialize_polyinfo(sr.Polysys,0,polyinfo_obj);
    sr.timedata[4] = (double)clock();
    //cout << "4 " << sr.timedata[4] << endl;
    
    stsize = sr.Polysys.numsys() -1 ;
    //generate all supports being consisted polynomial sdp without moment matrices
    vector<class poly_info> polyinfo_st;
    vector<class bass_info> bassinfo_st;
    if(stsize >= 1){
        polyinfo_st.resize(stsize);
        bassinfo_st.resize(stsize);
        get_subjectto_polys_and_basups(sr.Polysys,sr.bindices,BasisSupports.supsetArray,stsize,polyinfo_st,bassinfo_st);
    }
    sr.timedata[5] = (double)clock();
    //cout << "5 " << sr.timedata[5] << endl;
    get_allsups(sr.Polysys.dimvar(),polyinfo_obj,stsize,polyinfo_st,bassinfo_st,allsups_st);
    sr.timedata[6] = (double)clock();
    //cout << "6 " << sr.timedata[6] << endl;
    
    class supSet allSups;
    initialize_supset(allsups_st,allSups);
    sr.timedata[7] = (double)clock();
    //cout << "7 " << sr.timedata[7] << endl;
    
    //reduce each basis supports
    if(sr.param.reduceMomentMatSW == 1){
        sr.reduceSupSets(BasisSupports,allSups);
    }
    if(sr.param.detailedInfFile.empty() == false){
        sr.write_BasisSupports(1,sr.param.detailedInfFile,BasisSupports);
    }
    sr.timedata[8] = (double)clock();
    //cout << "8 " << sr.timedata[8] << endl;
    
    //eliminate supports of each basis supports, using special complementarity x(a)=0
    if(sr.param.complementaritySW==YES){
        get_removesups(sr.Polysys,removesups);
        initialize_supset(removesups,czSups);
        if(czSups.size()>0){
            sr.eraseCompZeroSups(czSups,BasisSupports.supsetArray);
        }
    }
    sr.timedata[9] = (double)clock();
    //cout << "9 " << sr.timedata[9] << endl;
    
    mmsize = BasisSupports.supsetArray.size() - sr.Polysys.numsys();
    vector<class bass_info> bassinfo_mm(mmsize);
    get_momentmatrix_basups(sr.Polysys,sr.bindices,BasisSupports.supsetArray,bassinfo_mm);
    sr.timedata[10] = (double)clock();
    //cout << "10 " << sr.timedata[10] << endl;
    
    get_allsups_in_momentmatrix(sr.Polysys.dimvar(),mmsize,bassinfo_mm,mmsups);
    sr.timedata[11] = (double)clock();
    //cout << "11 " << sr.timedata[11] << endl;
    
    //generate all supports being consisted POP
    allsups.alloc(allsups_st.pnz_size + mmsups.pnz_size,
    allsups_st.vap_size + mmsups.vap_size );
    allsups.pnz_size = 0;
    allsups.vap_size = 0;
    allsups.dim = allsups_st.dim;
    pushsups(allsups_st,allsups);
    if(mmsups.vap_size > 0){
        pushsups(mmsups,allsups);
    }
    simplification(allsups);
    sr.timedata[12] = (double)clock();
    //cout << "12 " << sr.timedata[12] << endl;
    
    //eliminate vain supports of allsupports, using special complementary supports x(a) = 0
    if(sr.param.complementaritySW == YES  && removesups.pnz_size > 0){
        remove_sups(removesups,allsups);
    }
    sr.timedata[13] = (double)clock();
    //cout << "13 " << sr.timedata[13] << endl;
    
    initialize_supset(allsups,allSups);
    sr.timedata[14] = (double)clock();
    //cout << "14 " << sr.timedata[14] << endl;
    
    bool flag = true;
    for(int i=0;i<sr.Polysys.dimVar ;i++){
        if(fabs(sr.Polysys.boundsNew.lbd(i)) > EPS || fabs(sr.Polysys.boundsNew.ubd(i)-1) > EPS){
            flag = false;
            break;
        }
    }
    int numofbds=0;
    if(sr.param.boundSW == 1 && flag){
        allSups.unique();
        class supSet OneSup, ZeroSup;
        ZeroSup = allSups;
        OneSup  = allSups;
        sr.Polysys.addBoundToPOP(ZeroSup,OneSup,numofbds);
    }else if(sr.param.boundSW == 2 && flag){
        allSups.unique();
        class supSet OneSup, ZeroSup;
        sr.redundant_OneBounds(BasisSupports,allSups,OneSup);
        sr.redundant_ZeroBounds(BasisSupports,allSups,ZeroSup);
        sr.Polysys.addBoundToPOP(ZeroSup,OneSup,numofbds);
    }else if(sr.param.boundSW == 2 && !flag){
        allSups.unique();
        allSups.sort();
        sr.Polysys.addBoundToPOP_simple(allSups,numofbds);
        //cout << "numofbds = " << numofbds << endl;
    }else if(sr.param.boundSW == 1 && !flag){
        allSups.unique();
        allSups.sort();
        sr.Polysys.addBoundToPOP_simple(allSups,numofbds);
        //cout << "numofbds = " << numofbds << endl;
    }
    sr.timedata[15] = (double)clock();
    //cout << "15 " << sr.timedata[15] << endl;
    
    
    mmsetSize=BasisSupports.supsetArray.size()-(sr.Polysys.numsys()-numofbds);
    for(int i=0;i<mmsetSize;i++){
        mmBaSupVect.push_back(BasisSupports.supsetArray[i+(sr.Polysys.numsys()-numofbds)]);
    }
    //eliminate basis supports of moment matrices from BasisSupports data;
    BasisSupports.supsetArray.resize(sr.Polysys.numsys()-numofbds);
    if(sr.param.boundSW == 1 || sr.param.boundSW == 2){
        for(int i=0;i<numofbds;i++){
            class supSet SE1Set;
            class sup    constSup;
            SE1Set.setDimVar(sr.Polysys.dimvar());
            SE1Set.pushSup(constSup);
            BasisSupports.push(SE1Set);
        }
    }
    sr.timedata[16] = (double)clock();
    //cout << "16 " << sr.timedata[16] << endl;
    
    //polyinfo,bassinfo
    int msize = sr.Polysys.numsys() + mmBaSupVect.size();
    vector<class poly_info> polyinfo(msize);
    vector<class spvec_array> bassinfo(msize);
    get_poly_a_bass_info(sr.Polysys,BasisSupports.supsetArray,mmBaSupVect,msize,polyinfo,bassinfo);
    sr.timedata[17] = (double)clock();
    //cout << "17 " << sr.timedata[17] << endl;
    
    //generate olynomial sdp
    get_psdp(sr.Polysys.dimvar(),msize,polyinfo,bassinfo,sdpdata);
    sr.timedata[18] = (double)clock();
    //cout << "18 " << sr.timedata[18] << endl;
    
    if(sr.param.complementaritySW == YES && removesups.pnz_size > 0){
        remove_sups(sdpdata,removesups);
    }
    sr.timedata[19] = (double)clock();
    //cout << "19 " << sr.timedata[19] << endl;
    
    //linearize polynomial sdp
    get_lsdp(allsups,sdpdata,sr.degOneTerms);
    if(sr.param.sdpaDataFile.empty() == false){
        write_sdpa(/*IN*/sdpdata,/*OUT*/ sr.param.sdpaDataFile);
    }
    sr.timedata[20] = (double)clock();
    //cout << "20 " << sr.timedata[20] << endl;
}

cliques::cliques(){
    numnode = 0;
    numcliques = 0;
}
cliques::~cliques(){
	for(int i=0;i<clique.size();i++){
		clique[i].clear();
	}
	clique.clear();
}
void cliques::initialize(int nodes,int csize){
    this->numnode = nodes;
    this->numcliques = csize;
	clique.resize(csize);	
}
int cliques::maxnnz(){
    
    int maxsize = 0;
	vector<list<int> >::iterator it=clique.begin();
	for(;it!=clique.end();++it){
		if(maxsize < (*it).size()){
			maxsize = (*it).size();
		}
	}
    return maxsize;
    
}
void cliques::write_maxCliques(string fname){
    std::ofstream fout;
    fout.open(fname.c_str(),ios::app);
    if(fout.fail()){
        cout << "errror@write_maxCliques:file does not open for output" << endl;
        cout << fname;
        exit(1);
    }
    
    int i,maxC,minC;
    maxC = clique[0].size();
    minC = clique[0].size();
    for(i = 1; i< this->numcliques;i++){
        if(maxC < clique[i].size()){
            maxC = clique[i].size();
        }
        if(minC > clique[i].size()){
            minC = clique[i].size();
        }
    }
	list<int>::iterator it;
    fout << "#Cliques = " << this->numcliques << ", maxC = " << maxC << ", minC = " << minC << endl;
    for(int j=0; j<this->numcliques;j++){
        fout << "clique " << j + 1 << " : ";
		it = clique[j].begin();
		for(;it!=clique[j].end();++it){
			fout<< (*it)+1 << " ";
		}
        fout << endl;
    }
    fout << endl;
    fout.close();
}


void gen_maxcliques3(int msize,vector<int> oriidx,class SparseMat extofcsp,class cliques & macls){
    
    int clsize = 1;
    bool doesinc;
	int nonzeros;   
    vector<bool> rowdummyIdx(msize);
    
    rowdummyIdx[0] = true;
   	nonzeros = extofcsp.jc[1] - extofcsp.jc[0];
	vector<int>::iterator it1,it2,it3,it4;
    for(int i=1;i<msize;i++){
        //check( inclusion relation between two cliques)
		it1 = extofcsp.ir.begin();
        doesinc = false;
		it1 = it1 + extofcsp.jc[i];
		it2 = it1 + extofcsp.jc[i+1] - extofcsp.jc[i];
        for(int j=0;j<i;j++){
			it3 = extofcsp.ir.begin();
			it3 = it3 + extofcsp.jc[j];
			it4 = it3 + extofcsp.jc[j+1] - extofcsp.jc[j];		
			doesinc = includes(it3,it4,it1,it2);
			if(doesinc == true){
				break;
			}	
        }
        
        //it's max clique
        if(doesinc == false){
            rowdummyIdx[i] = true;
            clsize ++;
            nonzeros += extofcsp.jc[i+1]-extofcsp.jc[i];
        }else{
			rowdummyIdx[i] = false;
		}
    }
    //generate max cliques
    macls.initialize(msize, clsize);
    
    clsize = 0;
    for(int i=0;i<msize;i++){
        if(rowdummyIdx[i] == true){
        	for(int j=extofcsp.jc[i];j<extofcsp.jc[i+1];j++){
				macls.clique[clsize].push_back(oriidx[extofcsp.ir[j]]);
        	}
			macls.clique[clsize].sort();	
            clsize++;
		}
    }
}

void copy_polynomialsdpdata(/*IN*/class mysdp & opsdp,/*OUT*/class mysdp & npsdp){
    
    npsdp.mDim = opsdp.mDim;
    npsdp.nBlocks = opsdp.nBlocks;
    
    npsdp.alloc(opsdp.nBlocks,opsdp.ele.sup.pnz_size,opsdp.ele.sup.vap_size);
    npsdp.ele.sup.pnz_size = opsdp.ele.sup.pnz_size;
    npsdp.ele.sup.vap_size = opsdp.ele.sup.vap_size;
    
    for(int i=0;i<npsdp.nBlocks;i++){
        npsdp.bLOCKsTruct[i]   = opsdp.bLOCKsTruct[i];
        npsdp.block_info[0][i] = opsdp.block_info[0][i];
        npsdp.block_info[1][i] = opsdp.block_info[1][i];
        npsdp.block_info[2][i] = opsdp.block_info[2][i];
    }
    
    for(int i=0;i<opsdp.nBlocks;i++){
        for(int j=opsdp.block_info[0][i];j<opsdp.block_info[0][i]+opsdp.block_info[1][i];j++){
            
            npsdp.ele.bij[0][j] = opsdp.ele.bij[0][j];
            npsdp.ele.bij[1][j] = opsdp.ele.bij[1][j];
            npsdp.ele.bij[2][j] = opsdp.ele.bij[2][j];
            
            npsdp.ele.sup.pnz[0][j] = opsdp.ele.sup.pnz[0][j];
            npsdp.ele.sup.pnz[1][j] = opsdp.ele.sup.pnz[1][j];
            
            npsdp.ele.coef[j]   = opsdp.ele.coef[j];
            
            for(int k=opsdp.ele.sup.pnz[0][j];k<opsdp.ele.sup.pnz[0][j]+opsdp.ele.sup.pnz[1][j];k++){
                npsdp.ele.sup.vap[0][k] = opsdp.ele.sup.vap[0][k];
                npsdp.ele.sup.vap[1][k] = opsdp.ele.sup.vap[1][k];
            }
            
        }
    }
    
    
}
void qsort_normal(/*IN*/vector<int> a, int left, int right,/*OUT*/vector<int> sortedorder)
{
    int i, j, stand;
    if (left < right){
        i = left; j = right; stand = a[sortedorder[left]];
        while(i <= j){
            while(a[sortedorder[i]]<stand) i++;
            while(a[sortedorder[j]]>stand) j--;
            if(i <= j) {
                swap2(sortedorder, i, j);
                i++;
                j--;
            }
        }
        qsort_normal(a, left, j,sortedorder);
        qsort_normal(a, i, right,sortedorder);
    }
}
void swap2(vector<int> a,const int & i,const int & j){
    int dum ;
    dum = a[i];
    a[i] = a[j];
    a[j] = dum;
}
void swap(int* a,const int & i,const int & j){
    int dum ;
    dum = a[i];
    a[i] = a[j];
    a[j] = dum;
}
bool same_edge(class EdgeSet edge1, class EdgeSet edge2){
			
	if(edge1.vertex1 == edge2.vertex1 && edge1.vertex2 == edge2.vertex2){
		return true;
	}
	if(edge1.vertex1 == edge2.vertex2 && edge1.vertex2 == edge2.vertex1){
		return true;
	}
	return false;
}
bool comp_edge(class EdgeSet edge1, class EdgeSet edge2){
			
	if(edge1.vertex2 < edge2.vertex2){
		return true;
	}else if(edge1.vertex2 == edge2.vertex2){
		if(edge1.vertex1 <= edge2.vertex1){
			return true;
		}else{
			return false;
		} 
	}else if(edge1.vertex2 > edge2.vertex2){
		return false;
	}
}
