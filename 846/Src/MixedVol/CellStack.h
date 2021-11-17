#ifndef _CELL_STACK_H_
#define _CELL_STACK_H_
#include <iostream>

class CellStack
{
   private:
      struct cell
      {
         int* idx;
         cell* next;

         cell(): next(0) {};
         // cell( int* J=0, cell* ptr=0 ): idx(J), next(ptr) {};
         cell( int & n, int* J, cell* ptr=0 ): next(ptr)
         {
            idx = new int [n];
            for(int i=0; i< n; i++) idx[i] = J[i];
         };
      };
      int size;
      int count;
      cell* top;
      cell* cur;
      
   public:
      CellStack(int & n) { size=n; count=0; cur=top=0; };
      ~CellStack() { MakeEmpty(); };
      // void Push( int* &J ) { top = new cell(J, top); };
      void Push( int* J )
      {
         cur = top = new cell(size, J, top);
         ++count;
      };
      void Pop()
      {
         cell *ptr = top;
         //if( cur == top ) cur = cur->next;
         cur = top = top->next;
         delete [] ptr->idx;
         delete ptr;
         --count;
      };
      
      bool Next()
      {
         //if( cur && cur->next )
         if( cur->next )
         { cur = cur->next; return true; }
         else
            return false;
      };

      int* Cur() { return cur->idx; }
      void Top() { cur = top; };
      void MakeEmpty()
      {
         while( !IsEmpty() ) Pop();
         count = 0;
      };
      
      bool IsEmpty() { return top == 0; };
      
      void Print()
      {
         cerr << "CELL::Print:\n";
         cell* p = top;
         while( p )
         {
            cout << "->[";
            for( int i=0; i <size; i++) cout <<' '<< p->idx[i];
            cout << " ] ";
            p= p->next;
         }
         cout << endl;
      };

      friend ostream & operator <<( ostream & OUT, const CellStack &CS )
      {
         OUT << "CELLS ("<< CS.count <<"):\n";
         cell* p;
         int c=0;
         for(p=CS.top; p; p=p->next )
         {
            OUT << "->[";
            for( int i=0; i <CS.size; i++) OUT <<' '<< p->idx[i];
            OUT << " ] ";
            if( !(++c % 6) ) OUT << endl;
         }
         OUT << endl;
         return OUT;
      };
      int & Count() { return count; }
      
      int* operator [](int &i)
      {
         if( -1< i && i< count )
         {
            int c=0;
            cur = top;
            while( c<i ) { cur=cur->next; c++;}
            return cur->idx;
         }
         else
            return 0;
      }
};
#endif
