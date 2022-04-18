/**************************************************************
                          Class Net
This class contains the variables and functions necessary to
create a netcdf file and store all the values (true and
interpolated) in netcdf format.
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <netcdf.h>
#include <string.h>
#include <strstream.h>
#include "Interp.h"
#include "Net.h"

//Net constructor - defines and creates a netcdf file

Net::Net(int DIM,int *Axes_Dim,float Axes_Ranges[][3],float **Axes_Tickes,int suff1,int suff2) : M(DIM){

  int i,j,k;      // loop indices
  // names of axes
  char *names[] = { "O3", "CO2", "N", "TMP", "PCP" }; 
  ostrstream Name;  // name of netcdf file if testing

  //define the dimensions of the hypervolume to be represented
  // in the netcdf file
  var1_dimids = new int[M]; 

  // determine if test case or not and name the netcdf file
  if((suff1==0)&&(suff2==0)){
    ncid  = nccreate("Interp_view.nc",NC_NOCLOBBER);
   }
  else {
    Name << "Test" << suff1 << suff2 << ".nc" << ends;
    char *Filename = Name.str();
    ncid  = nccreate(Filename,NC_NOCLOBBER);
   }

  //define each axis
  for(i=0;i<M;i++){
    var1_dimids[i] = ncdimdef(ncid, names[i],(long)Axes_Dim[i+1]);
   }

  //finish defining the netcdf variable
  if((var1_id = ncvardef(ncid,"var1", NC_FLOAT, M, var1_dimids))<0)
      printf("Error creating variable: var1_id = %d\n",var1_id);

  ncendef(ncid);  //done defining the netcdf file

 } //end of Net

/*
void
Net::Add_Knodes(Interpolator *interp){
   
    int i,j,k,
          myN,
          res;
    float **myX,
           *myW;


    myN = interp->getN();
    myX = interp->getX();
    myW = interp->getW();

cout << "myN = " << myN << endl;
    for(i=0;i<myN;i++){
      for(k=0;k<M;k++){
        start[k] = (long)myX[k+1][i+1] - 1; 
        cout << start[k] << " "; 
       }
      var1_vals = (double)myW[i+1];
      cout << "W[i] = " << var1_vals << endl;
      res = ncvarput1(ncid,var1_id,start,&var1_vals);
      if(res == -1)   {
         cout << "Could not Add Node " << i << endl ;
       }
     }
   cout << start[0] << " " << start[1] << " " << start[2] << " " << var1_vals << endl;
      res = ncvarput1(ncid,var1_id,start,&var1_vals);
      if(res == -1)   {
         cout << "Could not Add Node " << i << endl ;
       }
    var1_vals = 100;
    ncvarget1(ncid,var1_id,start,&var1_vals);
    cout << "Just did a get - var1_vals = " << var1_vals << endl;
  }

void
Net::All_Knodes(long *PP,float value){
   
   int i;
   long *P;
   double put_value;
   static int counter=0;

   put_value = (double)value;
   P = new long[M];
   for(i=0;i<M;i++){
    P[i] = PP[i+1]-1;
   }
   if(counter < 10){
     cout << P[0] << " " << P[1] << " " << P[2] << " " << value << " " << put_value << endl;
    }
   if(ncvarput1(ncid,var1_id,P,&put_value) == -1)
     {
        cout << "No go for ncvarput1 in All_Knodes" << endl;
     }
   counter++;
 }

*/

/**************************************************************
                             All_Knodes
This function receives a vector of floats along with the
coordinates and length of the vector so it can be stored in
the netcdf file.
  Input Parameters:
      coords        -array of coordinates where int_vector should
                       be stored
      edges         -array of lengths defining int_vector
      int_vector    -array holding values to be stored
  Functions Called:
      ncvarput      -netcdf function that stores a vector of
                     values of size defined by edges at a
                     location defined by coords
****************************************************************/
void
Net::All_Knodes(long *coords, long *edges, float *int_vector){

   if(ncvarput(ncid,var1_id,coords,edges,int_vector) == -1)
     {
        cout << "No go for ncvarput1 in All_Knodes" << endl;
     }
 }


//Destructor for Class Net
Net::~Net(){

  ncclose(ncid);  //close the netcdf file

 }
