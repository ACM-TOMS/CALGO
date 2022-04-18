/**************************************************************
                           Class Net
Defines the Class Net which includes the variables and functions
for creating and defining a netcdf file and for storing values
in the netcdf file.
****************************************************************/

#include <stdio.h>
#include <iostream.h>


class Net {
  private:
         int    M,       //dimensions of the hypervolume
             ncid,       //netcdf file number
     *var1_dimids,       //array containing axes' sizes
          var1_id;       //netcdf file descriptor

/*     long *start;
     long *count;

   double var1_vals;
*/

  public:
    //Net Constructor
    Net(int,int *,float[][3],float **,int,int);
   ~Net();

//    void Add_Knodes(Interpolator *);
//    void All_Knodes(long*,float);
    void All_Knodes(long*,long*,float*); // Stores values in file

 };
      

