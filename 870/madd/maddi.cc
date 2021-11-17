///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    MADDI                                                                  //
//    A line command user interface for the MADD library                     //
//    Ver. 1.0,  09/2005                                                     //
//                                                                           //
//    Author: Leonidas Linardakis                                            //
//    email:  lxlina@wm.edu                                                  //
//                                                                           //
//---------------------------------------------------------------------------//
//    This interface executes commands for reading and decomposing           
//    planar domains. The basic commands are ([] indicate an optional        
//    argument):
// 
//    read <filename>   : reads a file describing a geormtry in a .poly 
//                        format
//    dec [<number_of_subdomains>] : decomposes the geometry into 
//                        number_of_subdoains. if no argument is given, then 
//                        an area criterion is used.
//    read_weights <filename> : reads a file describing a weighted background 
//                        grid. The file is in poly format with one          
//                        attribute, a float number indicating the density   
//                        weight of the point.                                 
//    read_areas <filename> : reads a file describing an area background     
//                        grid. The file is in poly format with one          
//                        attribute, a float number indicating the required   
//                        area of the subdomain containing this point.         
//    setWeightFunction <library_name> <function> : sets a weight function.  
//                        The function is loaded dynamically from            
//                        <library_name>. An example library is given in      
//                        weightLibrary.cc                                    
//    setAreaFunction <library_name> <function> : sets an area function.     
//                        The function is loaded dynamically from            
//                        <library_name>. An example library is given in      
//                        weightLibrary.cc     
//    write <filename>  : writes to a file the current geometry including the
//                        separators.       
//    exit              : exit
//                         
//    set <param>=<value> : set a parameter. The values for <param> are:     
//            f : is the min angle created by the separators (in rads)       
//            g : the grid step for a structured background grid. This grid  
//                will be used to evaluate the weight or area function. If   
//                it is not set, then a value will be calculated             
//                automatically.                                             
//            r : the level of uniform refinement. A value of 2 will double  
//                the number of refining segments, 3 will triple it, and so 
//                on.
//            a : the gradation factor for density weights. The density 
//                weight of the subdomains is multiplied by a, while the 
//                area is multiplied by r. The total sum gives the density 
//                weight of the subdomain. A value a=0 will result no 
//                gradation.       
//            p : the interpolation factor for density weights. The weights 
//                of the new points are computed by linear interpolation 
//                of the weights of the existing end points, multipleid by 
//                p. A value p=0 eliminates the interpolation.              
//            c : is the max area to required area ratio for a subdomain to 
//                be decomposed. For example, if c=1, all the subdomains 
//                will be decomposed until none has area greater than the
//                required area, as defined by the background grid or the
//                area function.
//            d : defines the number of smoothing iterations before a 
//                separator is inserted. 
//            v : defines the number of vertices examined during the 
//                smoothing procedure.
//            i : (0.6 - 0.9) defines the max acceptable imbalance during 
//                the smoothing procedure.   
//
//---------------------------------------------------------------------------
//    An example of a simple set of commands:
//
//read pipe.poly
//dec 100
//write pipe100.poly
//exit  
//                                                                           
///////////////////////////////////////////////////////////////////////////////


using namespace std;


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include "L_Timer.h"
#include "madd.h"
#include <sys/time.h>


//---------------------------------------------------------------------------
// prototypes
void printParamHelp();
void readFile(char *name, char *param2, char *param3, int noPar);
void readWeightsFile(char *param1, char *param2, char *param3, int noPar);
void readAreasFile(char *param1, char *param2, char *param3, int noPar);
void writeFile(char *name, char *param2, char *param3, int noPar);
void writeSubdomainFile(char *param1, char *param2, char *param3, int noPar);
void writeSubpartitionFile(char *param1, char *param2, char *param3, int noPar);
void setParam(char *param1, char *param2, char *param3, int noPar);
void dec(char *param1, char *param2, char *param3, int noPar);
void weight(char *param1, char *param2, char *param3, int noPar);
void makeVerbal(char *param1, char *param2, char *param3, int noPar);
void getStats(char *param1, char *param2, char *param3, int noPar);
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------
// the MADD object
madd myPSLG;

// auxiliary global variables
L_Timer myTime;
long double decTime;
int noOfSubdomains;


//---------------------------------------------------------------------------
// the main menu
//---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
   
   char* inFileName, *outFileName, *outFileName2;
   int NoOfPieces;
   double decoupleArea, decoupleLength;
   char* triangle_param;
   int in_memory;
   DFloat breakCoeff;

   bool keepOn = true;
   char line[190];
   char command[40];
   char param1[80];
   char param2[80];
   char param3[80];

   while (keepOn){
	 cin.getline(line, 188);
	 command[0] = param1[0] = param2[0] = param3[0] = 0;
	 int noPar = sscanf(line, "%s%s%s%s", command, param1, param2, param3);

	 if (noPar <= 0)
	   continue;
	 
	 if (command[0] == '#'){
	   // it is a comment
	   continue;
	 }
	 
	 noPar--;
	 if (strcmp(command, "read") == 0){
	   readFile(param1, param2, param3, noPar);
	   continue;
	 }
	 if (strcmp(command, "read_weights") == 0){
	   readWeightsFile(param1, param2, param3, noPar);
	   continue;
	 }
	 if (strcmp(command, "read_areas") == 0){
	   readAreasFile(param1, param2, param3, noPar);
	   continue;
	 }
	 if (strcmp(command, "write") == 0){
	   writeFile(param1, param2, param3, noPar);
	   continue;
	 }
	 if (strcmp(command, "writeSubdomains") == 0){
	   writeSubdomainFile(param1, param2, param3, noPar);
	   continue;
	 }
	 if (strcmp(command, "writeSubpartitions") == 0){
	   writeSubpartitionFile(param1, param2, param3, noPar);
	   continue;
	 }
	 if (strcmp(command, "dec") == 0){
	   dec(param1, param2, param3, noPar);
	   continue;
	 }
	 if (strcmp(command, "weight") == 0){
	   weight(param1, param2, param3, noPar);
	   continue;
	 }
	 if (strcmp(command, "verbal") == 0){
	   makeVerbal(param1, param2, param3, noPar);
	   continue;
	 }
	 if (strcmp(command, "setWeightFunction") == 0){
	   int ret = myPSLG.setWeightFunction(param1, param2);
	   if (ret == 0)
		 cout << "Weight function is set." << endl;
	   continue;
	 }
	 if (strcmp(command, "setAreaFunction") == 0){
	   int ret = myPSLG.setAreaFunction(param1, param2);
	   if (ret == 0)
		 cout << "Area function is set." << endl;
	   continue;
	 }
	 if (strcmp(command, "set") == 0){
	   setParam(param1, param2, param3, noPar);
	   continue;
	 }
	 if (strcmp(command, "statistics") == 0){
	   getStats(param1, param2, param3, noPar);
	   continue;
	 }
	 if (strcmp(command, "exit") == 0){
	    keepOn = false;
	   continue;
	 }

	 cout << "?" << endl;  cout.flush();
   }// while

   cout << "\n Bye.\n"; cout.flush();

   return 0;
}


//---------------------------------------------------------------------------
// get statistics
//---------------------------------------------------------------------------
void getStats(char *param1, char *param2, char *param3, int noPar)
{
  double minAngle, meanAngle, phi;
  double minWeightedArea, maxWeightedArea, meanWeightedArea; 
  double minCommunication, maxCommunication, meanCommunication;
  long int finalSegments;
  double angleDistribution[20];

  myPSLG.getSepStatistics(angleDistribution, &minAngle, &meanAngle, &phi, &finalSegments,
						  &minWeightedArea, &maxWeightedArea, &meanWeightedArea, 
						  &minCommunication, &maxCommunication,  &meanCommunication);
 

  //  for (int i=0; i < 19; i++){
  //	cout << i << ". " << angleDistribution[i] << endl;
  // }
  // write the angle distribution
  FILE* anglefile = fopen("angles.out", "w");
  for (int i=0; i < 19; i++)
	fprintf(anglefile, "%d %lf\n", i * 10 + 5, angleDistribution[i]);
  fclose(anglefile);

  if (noPar > 0){
	// write to a file
	cout << "\n Stats to file " << param1; cout.flush();
	FILE* statfile = fopen(param1, "a");
	fprintf(statfile, "%d %Lf %f %f %f %ld %f %f %f %f %f %f\n", 
			noOfSubdomains, decTime, phi, minAngle, meanAngle, finalSegments,
			minWeightedArea, maxWeightedArea, meanWeightedArea,
			minCommunication, maxCommunication, meanCommunication);
	fclose(statfile);
	return;
  }

  // display to screen
  cout << endl << " Subdomains: " << noOfSubdomains
	   << "   Decomposing time: " << decTime
	   << "   Phi: " << phi
	   << "   Min angle: " << minAngle
	   << "   Mean angle: " << meanAngle << endl; 
  cout.flush();
  return;
}


//---------------------------------------------------------------------------
// decompose the geometry
//---------------------------------------------------------------------------
void dec(char *param1, char *param2, char *param3, int noPar)
{

  int noOfPieces, retValue;

  if (noPar <=0){
	noOfPieces = 0;
	cout << " decomposing using the min area criterion..." << endl;
  }
  else{
	// get the number of subdomains
	sscanf(param1, "%d", &noOfPieces);
  }

  myTime.start();
  retValue = myPSLG.decompose(noOfPieces);
  decTime = myTime.stop();

   if (retValue < 0){
	 cout << endl << "Error: Decomposition failed: " << retValue << endl;
   }
   else{
	 cout << endl << retValue << " subdomains were created in " << decTime << " sec." << endl;
   }
   noOfSubdomains = retValue;

}


//---------------------------------------------------------------------------
// set a parameter
//---------------------------------------------------------------------------
void setParam(char *param1, char *param2, char *param3, int noPar)
{
  if (noPar < 1){
	cout << "?" << endl;
	return;
  }
  double dvalue;
  int ivalue;
  int accepted = 0;
  sscanf(param1+2, "%lf", &dvalue);
  sscanf(param1+2, "%d", &ivalue);
  switch (param1[0]){
  case 'r': accepted = myPSLG.setUniformRefineLevel(dvalue); break;
  case 'a': accepted = myPSLG.setAdaptiveRefineLevel(dvalue); break;
  case 'f': accepted = myPSLG.setPhi(dvalue); break;
  case 'g': accepted = myPSLG.setGridStep(dvalue); break;
  case 'i': accepted = myPSLG.setMaxImbalanceLevel(dvalue); break;
  case 'p': accepted = myPSLG.setInterpolationWeightLevel(dvalue); break;
  case 'd': accepted = myPSLG.setMaxSmoothLevel(ivalue); break;
  case 's': accepted = myPSLG.setSmoothVertices(ivalue); break;
  case 'n': accepted = myPSLG.setMetiseMaxNode(ivalue); break;
  case 'e': accepted = myPSLG.setMetiseMaxEdge(ivalue); break;
  case 'c': accepted = myPSLG.setDecAreaRatio(dvalue); break;
  default: printParamHelp();
  }//switch

  if (accepted != 0){
	cout << endl << "Argument not accepted: " << dvalue << endl;
	printParamHelp(); 
  }
}


//---------------------------------------------------------------------------
// set the verbal level
//---------------------------------------------------------------------------
void makeVerbal(char *param1, char *param2, char *param3, int noPar)
{
  if (noPar < 1){
	myPSLG.setVerbal(1);
	return;
  }
  int verb_level;
  sscanf(param1, "%d", &verb_level);
  myPSLG.setVerbal(verb_level);
  return;
}


//---------------------------------------------------------------------------
// read the geometry
//---------------------------------------------------------------------------
void readFile(char *param1, char *param2, char *param3, int noPar)
{
  if (noPar < 1){
	cout << "?" << endl;
	return;
  }

   if (myPSLG.readPolyFile(param1) !=0)
	 cout << endl << " Error: Failed to read " << param1 << "." << endl;
   else
	 cout << endl << param1 << " was read." << endl;   
   cout.flush();
}


//---------------------------------------------------------------------------
// read a weighted background grid
//---------------------------------------------------------------------------
void readWeightsFile(char *param1, char *param2, char *param3, int noPar)
{
  if (noPar < 1){
	cout << "?" << endl;
	return;
  }

   if (myPSLG.readBackgroundWeights(param1) !=0)
	 cout << endl << " Error: Failed to read " << param1 << "." << endl;
   else
	 cout << endl << param1 << " was read." << endl;   
}


//---------------------------------------------------------------------------
// read an area background grid
//---------------------------------------------------------------------------
void readAreasFile(char *param1, char *param2, char *param3, int noPar)
{
  if (noPar < 1){
	cout << "?" << endl;
	return;
  }

   if (myPSLG.readBackgroundAreas(param1) !=0)
	 cout << endl << " Error: Failed to read " << param1 << "." << endl;
   else
	 cout << endl << param1 << " was read." << endl;   
}


//---------------------------------------------------------------------------
// write the geometry to a file
//---------------------------------------------------------------------------
void writeFile(char *param1, char *param2, char *param3, int noPar)
{
  if (noPar < 1){
	cout << "?" << endl;
	return;
  }

  if (noPar == 1){
	if (myPSLG.writePolyFile(param1, -1, 'n', 'y', 'n', -1) == 0){
	  cout << endl << param1 << " is written." << endl;   
	}
	return;
  }

  if (strcmp(param2, "ext") == 0){
  	myPSLG.writePolyFile(param1, -2, 'n', 'y', 'n', -1);
	return;
  }
  if (strcmp(param2, "all") == 0){
  	myPSLG.writePolyFileAll(param1);
	return;
  }

  int subdom;
  sscanf(param2, "%d", &subdom);
  if (noPar == 2){
	cout << endl << "Writing subdomain : " << subdom << endl; cout.flush(); 
  	myPSLG.writePolyFile(param1, subdom, 'n', 'n', 'n', -1);
	return;
  }
  
 cout << "?" << endl;
   
}

//---------------------------------------------------------------------------
// write the decomposition to a file
//---------------------------------------------------------------------------
void writeSubdomainFile(char *param1, char *param2, char *param3, int noPar)
{
  if (noPar < 1){
	cout << "?" << endl;
	return;
  }

  if (myPSLG.writeSubdomainFile(param1) == 0)
	  cout << endl << param1 << " is written." << endl;   
  
  return;
}


//---------------------------------------------------------------------------
// write the decomposition to a file
//---------------------------------------------------------------------------
void writeSubpartitionFile(char *param1, char *param2, char *param3, int noPar)
{
  if (noPar < 1){
	cout << "?" << endl;
	return;
  }

  if (myPSLG.writeSubpartitions(param1) == 0)
	  cout << endl << param1 << " is written." << endl;   
  
  return;
}


//---------------------------------------------------------------------------
// assign weights to all points
//---------------------------------------------------------------------------
void weight(char *param1, char *param2, char *param3, int noPar)
{
  if (noPar < 1){
	cout << "?" << endl;
	return;
  }

  int value;
  sscanf(param1, "%d", &value);
  if (noPar == 1){
	myPSLG.weightAllPoints(value);
	return;
  }

  int subdom;
  sscanf(param2, "%d", &subdom);
  myPSLG.weightSubdomainPoints(value, subdom);

   
}


//---------------------------------------------------------------------------
// print a help for the avaliable parameters
//---------------------------------------------------------------------------
void printParamHelp()
{

  cout << endl 
	   << endl << "variables:  r=<refine_level>,  double, from 0 to 10."
	   << endl << "            a=<adaptive_level>,  double, from 0 to 10."
	   << endl << "            f=<min_angle>,   double, from 0.2 to 0.4."
	   << endl << "            c=<decompose_area_Ratio>,   double, > 0"
	   << endl << "            i=<max_imbalance>,   double, from 0.6 to 0.9."
	   << endl << "            p=<interpolation_weight>,   double, from 0 to 1."
	   << endl << "            d=<smooth_depth>,   int, from 0 to 10."
	   << endl << "            v=<smooth_vertices>,   int, from 0 to 10."
	   << endl << "            n=<metis_max_node_weight>,   int, from 10 to 10000."
	   << endl << "            e=<metis_max_edge_weight>,   int, from 10 to 10000."; 
  cout.flush();

}

