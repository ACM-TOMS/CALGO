// IT_LP object coincides with the original index and dimension settings,
// but a bit waste at index 0 since level 0 is dealed by an "L0_IML" object.

#ifndef _IT_LP_H_
#define _IT_LP_H_

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <memory.h>

// Terms
//   IT:   Index Tree
//   LP:   Linear Programming multi-linked list
//   VB:   vertical branch in IT
//   HB:   horizontal branch in IT
//   cell: indices in HB when full
//------------------------------------------------------------------------------
struct LPdata
{
   int dim;
   double* xxx;
   double** INV;
   int* JJJ;
   
   LPdata(int &n, int* J, double* x, double** A): dim(n)
   {
      int i;
      JJJ = new int [n];
      xxx = new double [n];
      INV = new double* [n];
      INV[0] = new double [n*n];
      for(i=1; i<n; i++) INV[i] = INV[i-1] + n;
      for(i=0; i<n; i++)
      {
         JJJ[i] = J[i];
         xxx[i] = x[i];
         for(int j=0; j<n; j++) INV[i][j] = A[i][j];
      }
   };
};
//------------------------------------------------------------------------------
struct IndexNode
{
   int      idx;
   LPdata*   info;
   IndexNode*   S; // point south
   IndexNode(int i=0, LPdata* ptr=0): idx(i), info(ptr) { S=0; };
};
//------------------------------------------------------------------------------
struct LPPL // LP pointer link
{
   LPdata* addr;
   LPPL*   next;
   LPPL( LPdata* A=0, LPPL* N=0 ): addr(A), next(N) {};
};

/******************************************************************************/
class IT_LP
{
   public:
      IT_LP( int &, int* ); // Constructor
      ~IT_LP();             // Destructor
      void Add   ( int &n, int* sJ, int &N, int* J, double*, double** );
      void Add   ( int &oneidx, int &N, int* J, double*, double** );
      bool NextLevel();
      void FreeIT();
      void FreeLP();
      void StepBack() { NP[CurLevel--]=0; }
      void PrintIT(ostream &);
      void PrintLP(ostream &);
      //void PrintLPPL(ostream &);

      bool IsEmpty() { return CurLevel<1/*0?*/; } // for both IT and LP
      bool IsFull()  { return CurLevel+1>=MaxLevels ; } // for both IT and LP
      int* Cell() { return cell + 1; } // return the ptr containing HB indices(a cell)
      int& CurSptIdx() { return InSpt[CurLevel]; } // return idx of spt where cur node lies.
      int& MinNumPt()  { return minNP[CurLevel]; } // min # pts required for moving to next level
      int& NumPt()     { return NP[CurLevel]; } // # pts in cur VB of IT
      int& Level()     { return CurLevel; }
      int& CurLPdim()  { return DIM[CurLevel]; } // ??
      IndexNode* FixedIdxNdPtr() { return (CurLevel>=0? IT[CurLevel]: 0); } // cur VB head
      IndexNode* ResetCurLevelTo1();
      void RenewNP1();

   private:
      IndexNode** IT; // array of index links(of IndexNode type) with useless IT[0] 
      IndexNode** last; // last[] points to the last valid node in IT[]
      LPPL**      LP; // array of LP address links(of LPPL type) with useless LP[0] & LP[1] 
      LPPL*       LPlast; // points to the last valid node in LP[CurLevel]
      // WARNING: IT[0]=0 always, IT[1] initially points to dummy but will be replaced
      IndexNode*  curr; // Current IT node ptr, for Add(-).
      IndexNode*  prev; // parent of curr, for Add(-).
      int* DIM;     // dimension of LP in each level.// ??
      int* NP;      // NP[L]= # of nodes in IT[L], including the fixed node.
      int* cell;    // NEW FEATURE
      int* InSpt;   // NEW FEATURE
      int* minNP;   // minimal # of nodes in IT[L], excluding the fixed node.
      int  MaxLevels; // maximal # of total levels for IT[]
      int  CurLevel;  // Current level index.
      //// int  count; // # IT nodes has been compared after Find() is called
      bool Find( int & ); // Find an index in V-branch
};

/******************************************************************************/

IT_LP::IT_LP( int &nSpt, int* type ) // Constructor
{
   int i, j, itmp=0;
   for(i=0; i<nSpt; i++) itmp += type[i];
   MaxLevels = itmp + nSpt + 1; // "+1" is for unused level_0

   IT    = new IndexNode* [MaxLevels];
   last  = new IndexNode* [MaxLevels];
   LP    = new LPPL*      [MaxLevels];
   NP    = new int        [MaxLevels];
   DIM   = new int        [MaxLevels];
   cell  = new int        [MaxLevels]; // NEW FEATURE
   InSpt = new int        [MaxLevels]; // NEW FEATURE
   minNP = new int        [MaxLevels];
   
   memset( IT,    0, MaxLevels * sizeof(IndexNode*) );
   memset( last,  0, MaxLevels * sizeof(IndexNode*) );
   memset( LP,    0, MaxLevels * sizeof(LPPL*) );
   memset( NP,    0, MaxLevels * sizeof(int) );
   memset( DIM,   0, MaxLevels * sizeof(int) ); // ??
   memset( cell,  0, MaxLevels * sizeof(int) ); // NEW FEATURE
   memset( InSpt, 0, MaxLevels * sizeof(int) ); // NEW FEATURE
   memset( minNP, 0, MaxLevels * sizeof(int) );

   int sum=0;
   DIM[sum] = itmp++;
   for(i=0; i<nSpt; i++)
   {
      minNP[sum] = type[i]+1;
      InSpt[sum] = i;
      for(j=1; j<=type[i]; j++)
      {
         DIM[sum+j] = --itmp;
         minNP[sum+j] = type[i]+1-j;
      }
      sum += type[i]+1;
      if(sum < MaxLevels)  DIM[sum] = itmp;
   }

   last[1] = prev = curr = IT[1] = new IndexNode(-1); // points to dummy node
   NP[1] = 1;
   CurLevel = 1;

   for(i=0; i<MaxLevels; i++) LP[i] = new LPPL;
   // for(i=2; i<MaxLevels-1; i++) LP[i] = new LPPL;
   // create dummies for all links which will be used

   //cout << "MaxLevels=" << MaxLevels <<", minNP=[";
   //for(i=0; i<MaxLevels; i++) cout << ' ' << minNP[i];
   //cout << " ]\n";

   //cout << "LP DIM=[";
   //for(i=0; i<MaxLevels; i++) cout <<' '<< DIM[i];
   //cout << " ]\n";
}

//------------------------------------------------------------------------------
IT_LP:: ~IT_LP() // Destructor
{
   FreeIT();
   FreeLP(); 
   delete [] IT;
   delete [] last;
   delete [] LP;
   delete [] NP;
   delete [] DIM;
   delete [] cell;
   delete [] InSpt;
   delete [] minNP;
}

//------------------------------------------------------------------------------
void IT_LP :: FreeIT() // PUBLIC
{
   CurLevel = MaxLevels-1;
   while( CurLevel > 1 )
   {
      //cerr << "Free level " << CurLevel << ":\n";
      prev = IT[CurLevel];
      curr = prev->S;
      while( curr )
      {
         prev->S = curr->S;
         delete curr;
         curr = prev->S;
      }
      //IT[CurLevel]= 0;
      //NP[CurLevel]= 0;
      --CurLevel;
   }
   for(CurLevel=0; CurLevel<MaxLevels; CurLevel++) delete IT[CurLevel];
}

//------------------------------------------------------------------------------
void IT_LP :: FreeLP() // PUBLIC
{
   LPdata *lpp;
   CurLevel = MaxLevels-1;
   while( CurLevel > 1 )
   {
      //cerr << "Free level " << CurLevel << ":\n";
      LPlast = LP[CurLevel]->next;
      while( LPlast )
      {
         LP[CurLevel]->next = LPlast->next;
         lpp = LPlast->addr;
         if( lpp )
         {
            delete [] lpp->JJJ;
            delete [] lpp->xxx;
            delete [] (lpp->INV)[0];
            delete [] lpp->INV;
         }
         // delete lpp;
         delete LPlast;
         LPlast = LP[CurLevel]->next;
      }
      //IT[CurLevel]= 0;
      //NP[CurLevel]= 0;
      --CurLevel;
   }
   for(CurLevel=0; CurLevel<MaxLevels; CurLevel++) delete LP[CurLevel];
}

//------------------------------------------------------------------------------
void IT_LP :: PrintIT(ostream & OUT) // PUBLIC
{
   int L;
   IndexNode* p;
   OUT << "IT:\n";
   for(L=1; L< MaxLevels; L++)
   {
      p= IT[L];
      if( L==CurLevel ) OUT << " ->"; else OUT << "   "; OUT << L << ')';
      OUT << setw(2) << minNP[L]
         << setw(2) << NP[L];
      if( p ) OUT   << ' ' << p->info ;
      while( p )
      {
         OUT <<"- "<< setw(2) << p->idx ;
         if( p == last[L] ) OUT << "*"; else OUT << " ";
         p= p->S;
      }
      OUT << endl;
   }
}

//------------------------------------------------------------------------------
void IT_LP :: PrintLP(ostream & OUT) // PUBLIC
{
   IndexNode* q;
   int i, k;
   OUT << "LP:\n";
   for(k=1; k< MaxLevels; k++)
   {
      if( k==CurLevel ) OUT << " ->"; else OUT << "   "; OUT << k << ')';
      q = IT[k]; 
      if( q ){ OUT << q->idx << ": " << q->info;   q = q->S; }
      while( q )
      {
         OUT <<" -[";
         for(i=0; i<DIM[k] && q->info; i++) {
            OUT <<' '<< q->info->JJJ[i];
            if( q->info->JJJ[i] == q->idx )
               OUT <<'#';
         }
         OUT <<"] "<< q->info;
         q = q->S;
      }
      OUT << endl;
   }
   OUT << endl;
}
//------------------------------------------------------------------------------

IndexNode* IT_LP :: ResetCurLevelTo1()
{
   prev = IT[1];
   curr = prev->S;
   while( curr )
   {
      prev->S = curr->S;
      delete curr;
      curr = prev->S;
   }
   NP[1]=1;
   CurLevel=1;
   return IT[1];
}

//------------------------------------------------------------------------------
void IT_LP::RenewNP1() // assume IT[1] != 0 
{
   for(prev=IT[1]; prev->S; ++NP[1],prev=prev->S);
   last[1] = prev;
   cell[1] = IT[1]->idx;
}
//------------------------------------------------------------------------------
bool IT_LP::NextLevel() // PUBLIC
{
   if( CurLevel+1 < MaxLevels )
   {
      if( NP[ CurLevel ] <= minNP[ CurLevel ] )
      {
           //cout << "NextLevel: Not enough pts to go on.\n";
           // OUT << "NextLevel: Not enough pts to go on.\n";
         return false; // case 1
      }

      /*** Now IT[ CurLevel ] has enough point to go on ***/
      if( IT[ CurLevel+1 ] ) // backtracking (next level non-empty)
      {
         IndexNode* tmp = IT[ CurLevel+1 ];
         IT[ CurLevel+1 ] = tmp->S;
         tmp->S = last[CurLevel]->S;
         last[CurLevel]->S = tmp;
         //
         tmp = IT[ CurLevel ]->S;
         IT[ CurLevel ]->S = tmp->S;
         tmp->S = IT[ CurLevel+1 ];
         IT[ CurLevel+1 ] = tmp;
      }
      else // forward
      {
         IndexNode* tmp = IT[ CurLevel ]->S;
         IT[ CurLevel+1 ] = tmp;
         IT[ CurLevel ]->S = tmp->S;
         tmp->S = 0;
      }
      if( NP[CurLevel] == 2 ) last[CurLevel] = IT[CurLevel];
      --NP[CurLevel++];
      ++NP[CurLevel];
      last[CurLevel] = IT[CurLevel];

      curr = IT[CurLevel];
        cell[CurLevel] = curr->idx; // NEW FEATURE
      // cell[CurLevel] = IT[CurLevel]->idx; // alternative
      
      LPlast = LP[CurLevel]; 
      //cout <<"NextLevel: <"<< curr->idx <<"> moved successfully.\n";
      // PrintIT(cout);
      // OUT <<"NextLevel: <"<< curr->idx <<"> moved successfully.\n";
      // PrintIT(OUT); 
      return true; // case 2
   }
   else
   {
      //cerr << "NextLevel: Levels are full.\n";
      // OUT << "NextLevel: Levels are full.\n";
      return false; // case 3
   }
}

//------------------------------------------------------------------------------
bool IT_LP::Find( int & IDX ) // PRIVATE, compare curr with initialized prev
{
   //// for( curr = prev->S; count< NP[CurLevel]; count++, curr = curr->S )
   for( curr=prev->S; curr !=last[CurLevel]->S; curr=curr->S )
   {
      if ( IDX <= curr->idx )
      {
         if ( IDX == curr->idx )
            return true; // curr->idx = IDX
         else
            return false; // prev->idx < IDX < curr->idx
      }
      prev = prev->S;
   }
   return false; // => all indices are smaller than IDX,
   // i.e. *** < prev->idx < IDX ( prev->S = curr = 0 )
} 

//==================================================================//
// Base index JJ (size nn) is not assumed sorted *****. PUBLIC
void IT_LP :: Add( int &n, int* J, int &nn, int* JJ, double* X, double** A )
{
   //if( CurLevel == 1 ) { cerr << "\a\a\aNot allowed to add index!\n"; return; }

   int i, j, k;
   prev = IT[CurLevel]; // moved out of for loop *****
   //// count = 1;
   bool AddIdx = false;
   LPdata* lpp;
    
   // Create or Overwrite LP:
   for(i=0; i<n; i++) // searching for the 1st idx of J not in level 0
   {
      if( !Find(J[i]) ) // add J[i] after prev -----------------------|
      {
         if( LPlast->next )
         {
            lpp = LPlast->next->addr;
            for(k=0; k<nn; k++)
            {
               lpp->JJJ[k] = JJ[k];
               lpp->xxx[k] = X[k];
               for(j=0; j<nn; j++) lpp->INV[k][j] = A[k][j];
            }
         }
         else
         {
            lpp = new LPdata(nn,JJ,X,A);
            LPlast->next = new LPPL(lpp);
         }
         LPlast = LPlast->next;
         AddIdx = true;
         break;
      }
   }
    //-------------------------------------------------------------------
   // Create or Overwrite IT:
   //for(i=0; i<n; i++) // loop thru all pts in J
   while( AddIdx )
   {
      //cerr << "Adding J[" << i << "]=" << J[i] << endl;
      //if( !Find(J[i]) ) // add J[i] after prev ---------------------------------------------|
      //{
         if( ! last[CurLevel]->S ) // IT[CurLevel] is full
         {
            //lpp = new LPdata(nn,JJ,X,A);
            curr = new IndexNode( J[i], lpp );
            curr->S = prev->S;
            prev = prev->S = curr;
            if( last[CurLevel]==IT[CurLevel] || last[CurLevel]->idx < J[i] ) last[CurLevel] = curr;
                //            ==   by G
         }
         // Having spare slots for LP data ....
         else if( last[CurLevel] == prev ) // after *last[CurLevel]
         {
            IndexNode* ptr = prev->S;
            ptr->idx = J[i];
            ptr->info = lpp;
            /*
            LPdata* tmp = ptr->info;
            for(k=0; k<nn; k++)
            {
               tmp->JJJ[k] = JJ[k];
               tmp->xxx[k] = X[k];
               for(j=0; j<nn; j++) tmp->INV[k][j] = A[k][j];
            }
            */
            last[CurLevel] = ptr;
            prev = ptr;
         }
         else // intermediate position
         {
            IndexNode* ptr = last[CurLevel]->S;
            ptr->idx = J[i];
            ptr->info = lpp;
            /*
            LPdata* tmp = ptr->info;
            for(k=0; k<nn; k++)
            {
               tmp->JJJ[k] = JJ[k];
               tmp->xxx[k] = X[k];
               for(j=0; j<nn; j++) tmp->INV[k][j] = A[k][j];
            }
            */
            prev->S = ptr;
            last[CurLevel]->S = ptr->S;
            ptr->S = curr;
            prev = ptr;
         }
         ++NP[CurLevel];
         //// ++count;
      //}
      AddIdx = false;
      for(++i; i<n; i++) 
      {
         //cerr << i << " " << J[i] << endl;
         if( !Find(J[i]) ) 
         { 
            AddIdx = true;
            break;
         }
      }
   }
}

/**************************************************************************/
// adding 1 index only 
void IT_LP :: Add( int &oneidx, int &nn, int* JJ, double* X, double** A ) // PUBLIC
{
   int i,j;
   prev = IT[CurLevel]; // move out of this loop *****
   //// count = 1; 
   LPdata* lpp;

   if( !Find(oneidx) ) // add oneidx after prev ---------------------------------------------|
   {
      //cerr << "1-Adding point " << oneidx  << endl;
      if( LPlast->next )
      {
         lpp = LPlast->next->addr;
         for(i=0; i<nn; i++)
         {
            lpp->JJJ[i] = JJ[i];
            lpp->xxx[i] = X[i];
            for(j=0; j<nn; j++) lpp->INV[i][j] = A[i][j];
         }
      }
      else
      {
         lpp = new LPdata(nn,JJ,X,A);
         LPlast->next = new LPPL(lpp);
      }
      LPlast = LPlast->next;

      if( ! last[CurLevel]->S ) // IT[CurLevel] is full
      {
         //LPdata* LP = new LPdata(nn,JJ,X,A);
         curr = new IndexNode(oneidx,lpp);
         curr->S = prev->S;
         prev = prev->S = curr;
         if( last[CurLevel]->idx < oneidx ) 
             //need last[CurLevel]==IT{CurLevel] || in if ( ) ??
            last[CurLevel] = curr;
      }
      // Having spare slots for LP data ....
      else if( last[CurLevel] == prev ) // after *last[CurLevel]
      {
         IndexNode* ptr = prev->S;
         ptr->idx = oneidx;
         ptr->info = lpp;
         /*
         LPdata* tmp = ptr->info;
         for(i=0; i<nn; i++)
         {
            tmp->JJJ[i] = JJ[i];
            tmp->xxx[i] = X[i];
            for(j=0; j<nn; j++) tmp->INV[i][j] = A[i][j];
         }
         */
         last[CurLevel] = ptr;
         prev = ptr;
      }
      else // intermediate position before last, though last->S !=0
      {
         IndexNode* ptr = last[CurLevel]->S;
         ptr->idx = oneidx;
         ptr->info = lpp;
         /*
         LPdata* tmp = ptr->info;
         for(i=0; i<nn; i++)
         {
            tmp->JJJ[i] = JJ[i];
            tmp->xxx[i] = X[i];
            for(j=0; j<nn; j++) tmp->INV[i][j] = A[i][j];
         }
         */
         prev->S = ptr;
         last[CurLevel]->S = ptr->S;
         ptr->S = curr;
         prev = ptr;
      }
      ++NP[CurLevel];
   }
}

#endif
