#ifndef _Sorting_
#define _Sorting_  

void Sort(int n, int* J)
{
   int j, itmp;
   for ( int i=1; i<n; i++ )
   {
      itmp = J[i];
      for ( j = i; j>0 && itmp<J[j-1]; j--)
         J[j] = J[j-1];
      J[j] = itmp;
   }
}

#endif
