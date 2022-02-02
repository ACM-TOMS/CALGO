#ifndef _L0_IML_H_
#define _L0_IML_H_

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <memory.h>

struct LPdata;
struct IndexNode;
struct LPPL; // a structure for LP ptr link

struct L0IdxNode
{
   int      idx;
   IndexNode*   R; // point right
   L0IdxNode*   D; // point down
   L0IdxNode(int i=0): idx(i) { R=0; D=0; };
};


/*
L0head -> -1-0
           |
         v-h-h-h-0
           |
         v-h-h-0
           |
         v-h-0
*/
/******************************************************************************/

class L0_IML // level 0 index multi-linked list
{
   public:
      L0_IML(); // Constructor
      ~L0_IML(); // Destructor
      void Add( int n, int* sJ, int &N, int* J, double*, double** );
      void Add( int* sJ, int &N, int* J, double*, double** );
      bool Migrate( IndexNode* inp );
      void PrintIT( ostream & );
      void PrintLP( ostream & );
      void PrintLPPL( ostream & );
      void Free();

   private:
      L0IdxNode* L0head; 
      L0IdxNode* L0curr; // Current IT node ptr, for Add(-).
      L0IdxNode* L0prev; // parent of curr, for Add(-).
      IndexNode* curr; // Current IT node ptr, for Add(-).
      IndexNode* prev; // parent of curr, for Add(-).
      LPPL* LP1; // dummy head of the LP address link (used in level 1)

      bool FindInD( int & ); // Find an index in V-branch
      bool FindInR( int & ); // Find an index in H-branch
};

L0_IML :: L0_IML() // Constructor
{
   L0prev= L0curr= L0head= new L0IdxNode(-1);
   LP1= new LPPL;
}

L0_IML :: ~L0_IML() // Destructor
{
   Free();
   delete LP1;
   // L0head is removed by Migrate when empty.
}

//------------------------------------------------------------------------------
void L0_IML :: Free()
{
   // Clean LP1 link
   LPPL* P = LP1->next;
   while( P )
   {
      LP1->next = P->next;
      delete [] P->addr->JJJ;
      delete [] P->addr->xxx;
      delete [] (P->addr->INV)[0];
      delete [] P->addr->INV;
      //delete P->addr;
      delete P;
      P = LP1->next;
   }
   // Multi index links and L0head would be removed by Migrate(-) !
}
//------------------------------------------------------------------------------
void L0_IML :: PrintIT(ostream & OUT) // PUBLIC
{
   OUT << "I0:0)\n";
   L0curr = L0head->D;
   IndexNode* p;
   while( L0curr )
   {
      OUT << "     " << setw(2) << L0curr->idx;
      p = L0curr->R;
      while( p ) { OUT << "- " << setw(2) << p->idx;   p=p->S; }
      OUT << endl;
      L0curr = L0curr->D;
   }
   OUT << endl;
}
//------------------------------------------------------------------------------
void L0_IML :: PrintLP(ostream & OUT) // PUBLIC
{
   OUT << "L0:0)\n";
   L0curr = L0head->D;
   IndexNode* p;
   while( L0curr )
   {
      OUT <<"     "<< setw(2) << L0curr->idx <<' ';
      p = L0curr->R;
      while( p ) { OUT <<"- "<< setw(2) << p->idx <<':'<< p->info; p=p->S; }
      OUT << endl;
      L0curr = L0curr->D;
   }
   OUT << endl;
}
//------------------------------------------------------------------------------
void L0_IML :: PrintLPPL(ostream &OUT)
{
   OUT << "LP addresses passed to L1:\n  ";
   for(LPPL* P=LP1->next; P; P=P->next)
      OUT <<' '<< P->addr;
   OUT << endl <<endl;
}
//------------------------------------------------------------------------------
bool L0_IML :: Migrate( IndexNode* inp ) // PUBLIC, assume *inp has been created
{
   if( !inp ) { cerr << "\aMigrate: inp = 0!\n"; exit(1); }
   if( L0head->D )
   {
      L0prev = L0head->D;
      inp->idx  = L0prev->idx;
      // inp->info = 0;
      inp->S    = L0prev->R;

      L0head->D = L0prev->D;
      delete L0prev;
      //cout <<"Migrate: <"<< inp->idx <<"> moved successfully.\n";
      return true;
   }
   else
   {
      delete L0head;
         //cout << "Migrate: empty.\n";
      return false;
   }
}

//------------------------------------------------------------------------------
bool L0_IML :: FindInR( int & IDX ) // only compare curr in ->R with initial prev settled
{
   for( curr = prev->S; curr; curr = curr->S )
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

//------------------------------------------------------------------------------
bool L0_IML :: FindInD( int & IDX ) // only compare L0curr in VB with initial L0prev settled
{
   for( L0curr = L0prev->D; L0curr; L0curr = L0curr->D )
   {
      if ( IDX <= L0curr->idx )
      {
         if ( IDX == L0curr->idx )
            return true; // L0curr->idx = IDX
         else
            return false; // L0prev->idx < IDX < L0curr->idx
      }
      L0prev = L0prev->D;
   }
   return false; // => all indices are smaller than IDX,
   // i.e. *** < L0prev->idx < IDX ( L0prev->D = L0curr = 0 )
} 

//==================================================================//
// Assume LP (LP data pointer) is prepared and base index J (size n)
// is sorted in ascending order.
void L0_IML :: Add( int n, int* J, int &N, int* I, double* X, double** A )
{
   bool LPnotused = true;
   LPdata*   lpp = new LPdata(N,I,X,A);
   L0prev = L0head;
   int i;
   for(i=0; i<n; i++) // loop thru all pts in J
   {
      //cerr << "Adding J[" << i << "]=" << J[i] <<endl;
      int j=i+2;
      if( FindInD(J[i]) ) // then L0curr != 0
      {
         if( i+1<n )
         {
            if( L0curr->R )
            {
               L0prev = L0curr;
                 curr = L0curr->R;
               if( curr->idx > J[i+1] )
               {
                  curr = new IndexNode( J[i+1], lpp );
                  curr->S = L0prev->R;
                  L0prev->R = curr;
                  LPnotused = false;
               }
               else if( curr->idx < J[i+1] )
                  j=i+1;
               for( ; j<n; j++)
               {
                  prev = curr;
                   if( ! FindInR(J[j]) )
                  {
                     // Add2S( J[j], lpp );
                     curr = new IndexNode( J[j], lpp );
                     curr->S = prev->S;
                     prev->S = curr;
                     LPnotused = false;
                  }
               }
            }
            else
            {
               curr =  L0curr->R =  new IndexNode( J[i+1], lpp ); // Add2E( J[i+1], lpp );
               LPnotused = false;
               for( ; j<n; j++) 
               {
                  prev = curr;
                  // Add2S( J[j], lpp );
                  curr = new IndexNode( J[j], lpp );
                  curr->S = prev->S;
                  prev->S = curr;
               }
            }
         }
      }
      else // add J[i] after L0prev ---------------------------------------------|
      {
         L0curr = new L0IdxNode( J[i] );
         LPnotused = false;
         L0curr->D = L0prev->D;
         L0prev->D = L0curr; 
         if( i+1<n )
         {
            curr =  L0curr->R =  new IndexNode( J[i+1], lpp ); // Add2E( J[i+1], lpp );
            for( ; j<n; j++)
            {
               prev = curr;
               // Add2S( J[j], lpp );
               curr = new IndexNode( J[j], lpp );
               curr->S = prev->S;
               prev->S = curr;
            }
         }
      }
      L0prev = L0curr;
   }

   if( LPnotused )
   {
      delete [] lpp->JJJ;
      delete [] lpp->xxx;
      delete [] (lpp->INV)[0];
      delete [] lpp->INV;
      delete lpp;
   }
   else
      LP1->next = new LPPL(lpp,LP1->next);
}

//==================================================================//
// Assume LP (LP data pointer) is prepared and base index J (size n)
// is sorted in ascending order. #(J) = 2.
void L0_IML :: Add( int* J, int &N, int* I, double* X, double** A )
{
   LPdata*   lpp;
   L0prev = L0head;
   int i;
   for(i=0; i<2; i++) // loop thru all pts in J
   {
      //cerr << "Adding J[" << i << "]=" << J[i] <<endl;
      if( FindInD(J[i]) ) // then L0curr != 0
      {
         if( i==0 )
         {
            if( L0curr->R )
            {
               L0prev = L0curr;
                 curr = L0curr->R;
               if( curr->idx > J[1] )
               {
                  lpp = new LPdata(N,I,X,A);
                  LP1->next = new LPPL(lpp,LP1->next);
                  curr = new IndexNode( J[1], lpp );
                  curr->S = L0prev->R;
                  L0prev->R = curr;
               }
               else if( curr->idx < J[1] )
               {
                  prev = curr;
                   if( ! FindInR(J[1]) )
                  {
                     // Add2S( J[1], lpp );
                     lpp = new LPdata(N,I,X,A);
                     LP1->next = new LPPL(lpp,LP1->next);
                     curr = new IndexNode( J[1], lpp );
                     curr->S = prev->S;
                     prev->S = curr;
                  }
                  else return; // J[1] is in the right branch of J[0]
               }
            }
            else
            {
               lpp = new LPdata(N,I,X,A);
               LP1->next = new LPPL(lpp,LP1->next);
               L0curr->R = new IndexNode( J[1], lpp ); // Add2E( J[1], lpp );
            }
         }
      }
      else // add J[i] after L0prev ---------------------------------------------|
      {
         if( i==0 )
         {
            lpp = new LPdata(N,I,X,A);
            LP1->next = new LPPL(lpp,LP1->next);
            L0curr = new L0IdxNode( J[i] );
            L0curr->D = L0prev->D;
            L0prev->D = L0curr; 
            L0curr->R = new IndexNode( J[1], lpp ); // Add2E( J[1], lpp );
         }
         else
         {
            L0curr = new L0IdxNode( J[i] );
            L0curr->D = L0prev->D;
            L0prev->D = L0curr; 
         }
      }
      L0prev = L0curr;
   }
}

#endif
