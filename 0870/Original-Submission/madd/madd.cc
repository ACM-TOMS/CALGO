
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//          MM   MM      A     DD    DD                                      //
//          M M M M    A   A   D DD  D DD                                    //
//          M  M  M   A A A A  D  DD D  DD                                   //
//          M     M   A     A  D DD  D DD                                    //
//          M     M   A     A  DDD   DDD                                     //
//                                                                           //
//                                                                           //
//    A two-dimensional geometric domain decomposer                          //
//    Ver. 1.0,  09/2005                                                     //
//                                                                           //
//    Author: Leonidas Linardakis                                            //
//    email:  lxlina@wm.edu                                                  //
//                                                                           //
//                                                                           //
//---------------------------------------------------------------------------//
//    This program requires the following libraries:                         //
//                                                                           //
//    Triangle, by J.Shewchuk: http://www.cs.cmu.edu/~quake/triangle.html    //
//        is used to produce Delaunay triangulations                         //
//                                                                           //
//    Metis, by G.Karypis et al.: http://www-users.cs.umn.edu/~karypis/metis///
//        is used for graph partitioning                                     //
//                                                                           //
//    These libraries are avaliable from the corresponding sites,            //
//    under the conditions stated there and in their source code             //
//---------------------------------------------------------------------------//
//                                                                           //
//    The format of the input and output files describing a geometry is the  //
//    same as the poly files used by Triangle, see                           //
//    http://www.cs.cmu.edu/~quake/triangle.poly.html                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

using namespace std;

//---------------------------------------------------------------------------
//  headers
#include <iomanip>
#include <stdio.h>
#include <dlfcn.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#ifdef __sun
#include <floatingpoint.h>
#define isnormal(f) (f != HUGE_VAL)
//#include <sunmath.h>
#endif
#include "madd.h"
#define REAL double // this is for triangle.h
#include "triangle.h"
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Pi constants
#define TWO_PI_MINUS (2 * M_PI - 0.001) 
#define PI_MINUS (M_PI - 0.001)
#define TWO_PI (2 * M_PI)
#define PI_OVER_TWO (M_PI / 2)
#define PI_OVER_TWO_MINUS (M_PI / 2 - 0.001)
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// basic math shortcuts
#define sqr(x) ((x) * (x))
#define min(x, y) ((x) > (y)? (y): (x))
#define max(x, y) ((x) > (y)? (x): (y))
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// basic array access
#define x(p) (*p)                              // the x coordinate of a point
#define y(p) (*(p+1))                          // the y coordinate of a point
#define point(i) (thePointsList+((i) << 1))         // the DFloat* pointer of the i point
#define point_x(i) thePointsList[((i) << 1)]        // the x coordinate of point i
#define point_y(i) thePointsList[((i) << 1)+1]      // the x coordinate of point i
#define backgroundNode(i) (backgroundNodeList+((i) << 1))   // the DFloat* pointer of the i background node
#define background_x(i) (backgroundNodeList[((i) << 1)])    // the x coordinate of the i background node
#define background_y(i) (backgroundNodeList[((i) << 1)+1])  // the y coordinate of the i background node 
#define weightOfNode(i) backgroundAreaList[i]            // the weight parameter (area) of the i background node
#define hole(i) (theHolesList+((i) << 1))           // the DFloat* pointer of the i hole
#define hole_x(i) theHolesList[((i) << 1)]          // the x coordinate of the i hole
#define hole_y(i) theHolesList[((i) << 1)+1]        // the y coordinate of the i hole
#define segm(i) (theSegmsList+((i) << 1))           // the LInt* pointer of segment i
#define segm_point1(i) theSegmsList[((i) << 1)]     // the first point of segment i
#define segm_point2(i) theSegmsList[((i) << 1)+1]   // the second point of segment i
#define segm_colour(i) segmColourList[(i)]     // the color (after partitioning) of work segement i
#define tria_point1(i) theTriaList[(i)*3]      // the first point (LInt) of triangle i
#define tria_point2(i) theTriaList[(i)*3+1]    // the second point (LInt) of triangle i
#define tria_point3(i) theTriaList[(i)*3+2]    // the third point (LInt) of triangle i
#define tria(i) (theTriaList+(i)*3)            // Lint* pointer to triangle i 
#define neighborTria(i, j) triaNeighborList[3*(i)+(j)]  // the adjacent triangle to i triangle, throught the j (0-2) face
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// some short routines
// compute the weight value of subdomain s
#define  subdomainValue(s) (DOMAIN_AREA_WEIGHT_COEF * subdomainAreaCoef * subdomainArea[s]\
                         +  POINT_WEIGHT_COEF * subdomainWeightCoef * subdomainWeight[s])
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// global parameters
DFloat PHI  = (M_PI/3 * 0.9);         // the min angle 
DFloat AREA_WEIGHT_COEF = 2;          // the uniform decomposition factor 
DFloat DOMAIN_AREA_WEIGHT_COEF;       // the uniform decomposition factor, >= 1
DFloat POINT_WEIGHT_COEF = 0;         // the weight factor
DFloat INTERPOLATION_WEIGHT_COEF = 0; // the interpolation parameter for the weights
DFloat GRID_STEP = 0;                 // the uniform grid step
DFloat DEC_AREA_RATIO = 1;            // the min area ratio for decomposing a subdomain
int MAX_NO_OF_SUBDOMAINS = 1400;      // the default max number of subdomains
int VERBAL = 0;                       // the Verbal level
bool  BACKGROUND_WEIGHTS = false;     // if we use background weights
bool  BACKGROUND_AREAS = false;       // if we use background areas
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------*
//                                                                           *
//     START OF THE MAIN DECOMPOSITION ROUTINES                              *
//                                                                           *
//---------------------------------------------------------------------------*

//---------------------------------------------------------------------------
// the main decomposition routine
// input the number of subdomains to be created
// if subDomains <= 0 use the area criterion
// returns the number of the created subdomains
//---------------------------------------------------------------------------
int
madd::decompose(int subDomains)
{
  if (subDomains > 0){
	// decompose into the given number of subdomains
	totalNoOfSubdomains = subDomains;
  }
  else{
	// we have an area defined decomposition
	if (!BACKGROUND_AREAS && areaFunction == NULL){
	  // no area criteria are given; return error
	  return -1;
	}
	totalNoOfSubdomains = MAX_NO_OF_SUBDOMAINS;
  	if (AREA_WEIGHT_COEF < 3)
	  AREA_WEIGHT_COEF = 3;
  }//if

  //--------------------------------------
  // compute the refinement 
  DOMAIN_AREA_WEIGHT_COEF = max(AREA_WEIGHT_COEF, 1);
  int refPoints = (int) (sqrt((double)totalNoOfSubdomains) * 8);
  computeDomainLengths();
  maxPointWeight = computeMaxPointWeight();
  // get the number of total required points
  LInt refinePoints = (noOfPoints < refPoints) ? refPoints: noOfPoints; 
  if (refinePoints < 50)
	refinePoints = 50;
  if (refinePoints > 1000)
	refinePoints = 1000;
  // the the required refining length of the boundary segments
  refineLength = totalBoundaryLength / refinePoints;
#ifdef __VERBAL2
  cout << endl << "refineLength : " << refineLength 
	   << "   totalBoundaryLength : " << totalBoundaryLength 
	   << "   refinePoints : " << refinePoints  
	   << endl; cout.flush();
	   
  cout << endl << " Refining boundary... "; cout.flush();
#endif
  // break all segemnts according to the refineLength
  // set the basic flag for the new points to 2 
  breakAllSegm(refineLength, 2);
  // refinement computed 
  //--------------------------------------

  // this is the number of the boundary points and segments
  noOfExtPoints = noOfPoints;
  noOfExtSegms = noOfSegms;
#ifdef __VERBAL2 
  cout << endl << " allocateSubdomainStr... "; cout.flush();
#endif
  // allocate the subdomains
  allocateSubdomainStr();

#ifdef __VERBAL2 
  cout << endl << " createInitialSubdomain... "; cout.flush();
#endif
  createInitialSubdomain();
#ifdef __VERBAL2  
  cout << endl << " initSubdomain... "; cout.flush();
#endif
  initSubdomain(0);

#ifdef __WRITE_GRID
  // just write the uniform grid for debuging
  DFloat theMin, theMax, theTotal;
  gridFunction = areaFunction;
  if (GRID_STEP <= 0){
	// the grid step will be defined automatically to  200 at the longest coordinate
	DFloat longest = max(domainMaxX-domainMinX, domainMaxY-domainMinY);
	GRID_STEP = (float)longest / 200;
  }// gridstep defined
  evaluateGrid(0, &theMin, &theMax, &theTotal);
  writePolyFile((char *)"grid.poly", -1, 'n', 'y', 'y', -1, -1);
  return 0;
#endif

  //--------------------------------------
  // start the decomposition 
  //--------------------------------------
  if (subDomains <=0){
	//--------------------------------------
	// we will decompose on an area ratio criterion
	//--------------------------------------
	// the grid function is the area function
	gridFunction = areaFunction;
	// check if the grid step is set
	if (GRID_STEP <= 0){
	  // grid step is not set
	  // the grid step will be defined automatically to  200 at the longest coordinate
	  DFloat longest = max(domainMaxX-domainMinX, domainMaxY-domainMinY);
	  GRID_STEP = (float)longest / 200;
	}
	// gridstep defined
	// partition the initial domain
	partition(0);

	//--------------------------------------
	// keep decomposingh while there are 
	// subdomains with larger area ratio than 
	// required or the maximum number of 
	// subdomains has been reached
	//--------------------------------------
	bool keepDecompose = true;
	while(keepDecompose){
	  keepDecompose = false;
	  //--------------------------------------
	  // get the relatively bigger subdomain
	  int maxSubdomain  = 0;
	  if (subdomainMinArea[0] <= 0)
		subdomainMinArea[0] = 1;
	  DFloat maxRatio = subdomainArea[0] / subdomainMinArea[0];
	  //	  cout << "0. Area: " <<  subdomainArea[0] << "   MinArea: " << subdomainMinArea[0] << endl;
	  for (int i=1; i < noOfSubdomains; i++){
		if (subdomainMinArea[i] <= 0)
		  subdomainMinArea[i] = 1;
		//	cout << i << ". Area: " <<  subdomainArea[i] << "   MinArea: " << subdomainMinArea[i] << endl;
		if (subdomainArea[i] / subdomainMinArea[i] > maxRatio){
		  maxRatio = subdomainArea[i] / subdomainMinArea[i];
		  maxSubdomain = i;
		}
	  }//for all subdomains
	  //  bigger subdomain is found in maxSubdomain
	  //--------------------------------------

	  // if the area ratio is greater than rquired and 
	  // number of subdomains area less than max
	  // partition the subdomains and keep decomposing
	  if (maxRatio > DEC_AREA_RATIO && noOfSubdomains < totalNoOfSubdomains){
		if (VERBAL > 0){
		  cout << endl << endl << " ------------------------------------- " ;
		  cout << endl << " Creating the " << noOfSubdomains << " subdomain.";
		  cout << endl << " decomposing subdomain: " << maxSubdomain << " with area " << subdomainArea[maxSubdomain]
			   << ".  The min area for this subdomain is " <<  subdomainMinArea[maxSubdomain]; cout.flush();
		  // partition  the bigger subdomain
		  cout << endl << " ------------------------------------- "  << endl; cout.flush();
		}

		keepDecompose = true;
		partition(maxSubdomain);
		// subdomain is partitoned
		//--------------------------------------		
	  }// if
	  
	}// while(keepDecompose)
	//--------------------------------------
	// decomposition is finished
	// return the number of subdomains 
	//--------------------------------------
	return noOfSubdomains;
  }// if subDomains <=0
  // end of are-wise decompostion
  //--------------------------------------

  //--------------------------------------
  // decompose into totalNoOfSubdomains subdomains
  //--------------------------------------
  if (BACKGROUND_AREAS){
	//--------------------------------------
	// change the background areas to weights 
	DFloat domainArea = areaOfSubdomain(0);
	unsigned int startWeight = 0;
	unsigned int endWeight = backgroundListSize;
	//	cout << endl << DOMAIN_AREA_WEIGHT_COEF << "   " << POINT_WEIGHT_COEF << endl;
	for (unsigned int i=0; i < endWeight; i++){
	     weightOfNode(i) = (domainArea - weightOfNode(i));
	  // weightOfNode(i) = 1 /  weightOfNode(i);
	  //	  cout << endl << i+1 << "    " <<  weightOfNode(i);
	}
	//--------------------------------------
  }

  while(noOfSubdomains < totalNoOfSubdomains){


	//#if 0	
	//--------------------------------------
	// find the max area, weight, length of subdomains
	// deactivated 
	DFloat maxArea =  subdomainArea[0];
	DFloat maxLength =  subdomainLength[0];
	DFloat  maxWeight =  subdomainWeight[0];
	for (int i=1; i < noOfSubdomains; i++){
	  if (maxArea  < subdomainArea[i]){
		maxArea = subdomainArea[i];
	  }
	  if (maxLength < subdomainLength[i]){
		maxLength =  subdomainLength[i];
	  }
	  if (maxWeight  < subdomainWeight[i]){
		maxWeight =  subdomainWeight[i];
	  }
	}
	// find the normalizing coefficients. the max value will be 1000.
	if (maxWeight <= 0)
	  maxWeight = 250;
	DFloat subdomainAreaCoef = 1000 / maxArea;
	DFloat subdomainLengthCoef = 1000 / maxLength;
	DFloat subdomainWeightCoef =  (DFloat)1000 / (DFloat)maxWeight;
	//subdomainAreaCoef = 1;
	//subdomainLengthCoef = 1;
	//subdomainWeightCoef =  1;
	// max area, weight, length of subdomains found
	//--------------------------------------
	//#endif


#ifdef __VERBAL2
	cout << endl << "Initial Domain Parameters: " <<  subdomainArea[0] << ",  " <<  subdomainLength[0] 
		 << ",   " << subdomainWeight[0];    
	cout << endl << "Coefficiens: " <<  subdomainAreaCoef << ",  " <<  subdomainLengthCoef 
		 << ",   " << subdomainWeightCoef << endl;    
 
	cout << endl << "maxWeight: " << maxWeight << ",  subdomainWeightCoef: " << subdomainWeightCoef << endl; cout.flush();
#endif

	//--------------------------------------
	// find the bigger subdomain
	DFloat maxValue = subdomainValue(0);
	int maxSubdomain = 0;
	for (int i=0; i < noOfSubdomains; i++){
#ifdef __VERBAL
  	  cout << endl << " Subdomain: " << i << "  has area " << subdomainArea[i] 
		   << ",  length " << subdomainLength[i] << "  and weight " << subdomainWeight[i] 
		   <<  " giving a value of " << subdomainValue(i); cout.flush();
#endif
	  
	  if (maxValue <  subdomainValue(i)){
		maxValue =  subdomainValue(i);
		maxSubdomain = i;
	  }
	}
	// maxSubdomain is the bigger subdomain
	//--------------------------------------
	//	cout << endl << DOMAIN_AREA_WEIGHT_COEF << "   " << POINT_WEIGHT_COEF << endl;
	if (VERBAL > 0){
	  cout << endl << endl << " ------------------------------------- " ;
	  cout << endl << " Creating the " << noOfSubdomains << " subdomain.";
	  cout << endl << " Max Subdomain Value: " << maxValue << ",  area: " << subdomainArea[maxSubdomain] << 
		"   weight: " << subdomainWeight[maxSubdomain]; cout.flush();
	  cout << endl << " Max Subdomain: " << maxSubdomain; cout.flush();
	  // partition  the bigger subdomain
	  cout << endl << " ------------------------------------- "  << endl; cout.flush();
	}
	// partiton the  bigger subdomain
    partition(maxSubdomain);
  }//while(noOfSubdomains < totalNoOfSubdomains)

  //--------------------------------------
  // decomposition is done
  // return the number of subdomains
  //--------------------------------------
  return noOfSubdomains;
}//decompose()

 

//---------------------------------------------------------------------------
// the main partiton routine
// it partitions theSubDomain into two subdomains
// returns 0 if successful
//---------------------------------------------------------------------------
int
madd::partition(int theSubDomain)
{
  //#define __DEBUG_MODE
#ifdef __DEBUG_MODE
  writePolyFile((char *)"work.poly", theSubDomain, 'n', 'n', 'n', -1);
#endif
#undef __DEBUG_MODE

  // get theSubDomain into the work subdomains
  initSubdomain(theSubDomain);
  // create the work holes
  createWorkHoles(theSubDomain);

  if (VERBAL > 1){
	cout << endl << "Triangulating..."; cout.flush();
  }
  // get the Delauany triangulation
  triangulateDomain(NULL);

  if (VERBAL > 1){
	cout << endl << "Creating MADD graph and separators..."; cout.flush();
  }
  // invoke the MADD procedure
  if (MADDpartition()) {
	ErrorExit((char *)"Unable to create the MADD Graph.");
  }

  if (VERBAL > 1){
	cout << endl << "Idenifying the subdomains..."; cout.flush();
  }

#ifdef __VERBAL2
  cout <<  endl << "Refine separator..."; cout.flush();
#endif

  // refine the inserted separators
  refineThisSepLength();


  //--------------------------------------
  // create the two subdomains
  //--------------------------------------
#ifdef __VERBAL
  cout <<  endl << "Create Subdomains..."; cout.flush();
#endif
  // get the number of segments for each subdomain
  LInt noOfSep = madSepEnd - madSepStart;
  LInt noOfSubSegms1 = edgesPerColour[0] + noOfSep;
  LInt noOfSubSegms2 = edgesPerColour[1] + noOfSep;
#ifdef __VERBAL2
  cout << endl << " edgesPerColour " << edgesPerColour[0] << ",   " << edgesPerColour[1]; cout.flush();
#endif
  // allocate the subdomains segment list
  LInt* subSegments1 = (LInt*)malloc((noOfSubSegms1)*sizeof(LInt) + 1);
  LInt* subSegments2 = (LInt*)malloc((noOfSubSegms2)*sizeof(LInt) + 1);
  if (subSegments1 == NULL || subSegments2 == NULL)
	InsufficientMem((char *)"Allocating subdomains segment lis");
  // get the indeces of the two subdomains
  int subdomain1 = workSubdomain;
  int subdomain2 = noOfSubdomains++;
  // initialize the weights
  subdomainLength[subdomain1] = subdomainLength[subdomain2] = 0;
  subdomainWeight[subdomain1] = subdomainWeight[subdomain2] = 0;

  //--------------------------------------
  // identify which segment belongs to which subdomain
  //--------------------------------------
  LInt sg1 = 0;
  LInt sg2 = 0;
  LInt segm, p1, p2;
  DFloat *point1, *point2;
  DFloat midweight = 0;
  // go throu all the work segments
  for (LInt i = 0; i < workNoOfSegms; i++){
	// get the points and the weights of this segment
	segm =  workSegm[i];
	p1 = segm_point1(segm);
	p2 = segm_point2(segm);
	point1 = point(p1);
	point2 = point(p2);
	midweight = (pointWeight[p1] + pointWeight[p2]) / 2.0;
	//	 cout << endl << "Point " << p1 << " weight " <<  (int)pointWeight[p1]
	//		  << "Point " << p2 << " weight " <<  (int)pointWeight[p2];
	if ( segm_colour(i) == 0){
	  // the segment belongs to the first subdomain
	  subSegments1[sg1++] = segm;
	  subdomainLength[subdomain1] += distance(point(p1), point(p2)); 
	  subdomainWeight[subdomain1] += midweight; 

	}
	else{
	  if ( segm_colour(i) == 1){
	  // the segment belongs to the second subdomain
		subSegments2[sg2++] = segm;
		subdomainLength[subdomain2] += distance(point(p1), point(p2)); 
		subdomainWeight[subdomain2] += midweight; 
	  }//if
	}//if

  }// for all the work segments
  //  the work segments have been assigned 
  //  to the two subdomains
  //--------------------------------------

#ifdef __VERBAL2
  cout <<  endl << "       insert separators in both Subdomains..."; cout.flush();
#endif

  //--------------------------------------
  // insert separators to both subdomains
  for (LInt i = madSepStart; i < madSepEnd; i++){
	//	 cout <<  endl << "           " << sg1 << ",  " << sg2 << "  := " << i; cout.flush();
	subSegments1[sg1++] = i;
	subSegments2[sg2++] = i;
	p1 = segm_point1(i);
	p2 = segm_point2(i);
	point1 = point(p1);
	point2 = point(p2);
	DFloat dist =  distance(point1, point2);
	midweight = (pointWeight[p1] + pointWeight[p2]) / 2.0;
	subdomainLength[subdomain1] += dist;
	subdomainLength[subdomain2] += dist;
	subdomainWeight[subdomain1] += midweight; 
	subdomainWeight[subdomain2] += midweight; 
  }
  // separators have been inserted
  //--------------------------------------
  // all segments have been assgned to the 
  // two subdomains 
  //--------------------------------------

#ifdef __VERBAL2
  cout <<  endl << " the number of segments:." << noOfSubSegms1 << " =? " 
	   << sg1 << ",   " << noOfSubSegms2 << " =? " << sg2  ; cout.flush();
#endif

  //--------------------------------------
  // update the subdomain strucure
  //--------------------------------------
  subdomainArea[subdomain1] = partArea[0];
  subdomainArea[subdomain2] = partArea[1];
  subdomainSegments[subdomain1] = subSegments1;
  subdomainSegments[subdomain2] = subSegments2;
  noOfSubSegms[subdomain1] = sg1;
  noOfSubSegms[subdomain2] = sg2;
  if (sg1 != noOfSubSegms1 || sg2 != noOfSubSegms2){
	//	 ErrorExit((char *)"Error calculating the subdomain segments.");
	cout << endl << "Error calculating the subdomain segments."; cout.flush();
  }
  if (sg1 > noOfSubSegms1 || sg2 > noOfSubSegms2){
	ErrorExit((char *)"Error calculating the subdomain segments.");
  }
	 
  // delete previous subdomain
  free(workSegm);
#ifdef __VERBAL
  cout << endl << "clear partition..."; cout.flush();
#endif
  // delete triangulation
  clearPartition();

  //--------------------------------------
  //  compute the minareas and weights
  //  for the two subdomanis 
  //  evaluate the background grids
  //--------------------------------------
  int parentSubdomain = subdomain1;
  subdomainMinArea[subdomain1] = subdomainArea[subdomain1] + 10;
  subdomainMinArea[subdomain2] = subdomainArea[subdomain2] + 10;

  // get the non strucutred background grid
  if (backgroundListSize > 0) 
	getBackgroundWeights(parentSubdomain, subdomain2);

  //--------------------------------------
  // get the structured background grid
  if (GRID_STEP > 0 && gridFunction != NULL){
	DFloat theMin, theMax, theMean, theTotal;
	// create and evaluate the grid for the first subdomains
	evaluateGrid(subdomain1, &theMin, &theMax, &theMean, &theTotal);
	subdomainWeight[subdomain1] += theTotal;  
	subdomainMinArea[subdomain1] = min(theMin, subdomainMinArea[subdomain1]); 
#ifdef __VERBAL2
	cout << endl << "  -> The min of the grid is " << theMin 
		 <<  "  -> The max of the grid is " << theMax;
#endif

	// create and evaluate the grid for the first subdomains
	evaluateGrid(subdomain2, &theMin, &theMax, &theMean, &theTotal);
	subdomainWeight[subdomain2] += theTotal;  
	subdomainMinArea[subdomain2] = min(theMin, subdomainMinArea[subdomain2]);  
#ifdef __VERBAL2
	cout << endl << "  -> The min of the grid is " << theMin 
		 <<  "  -> The max of the grid is " << theMax;
#endif
  }// if grid
  // structured ackground grid evaluated
  //--------------------------------------
  //  minareas and weights have been 
  //  computed for the two subdomanis 
  //--------------------------------------

  //--------------------------------------
  //  assign the holes to the two subdomains
  //--------------------------------------
  {
#ifdef __VERBAL
	cout << endl << "define subdomain holes...."; cout.flush();
#endif
	// store parent holes
	uchar* subHole = subdomainHoles[parentSubdomain];	
	for (int i=0; i < noOfHoles; i++){
	  parentSubdomainHoles[i] = subHole[i];
	}

	// check the hole of the first domain
   	findSubdomainHoles(subdomain1);

	// remove from parent the ones we found
	subHole = subdomainHoles[subdomain1];	
	for (int i=0; i < noOfHoles; i++){
	  parentSubdomainHoles[i] -= subHole[i];
	}

	// check the hole of the second domain
   	findSubdomainHoles(subdomain2);
  }
  //--------------------------------------
  //  the holes have been assigned to the two subdomains
  //--------------------------------------

#ifdef __VERBAL
  cout << endl << "--------------------------------------------------";
  cout << endl << " Subdomain Areas of " << subdomain1 << ": " <<  subdomainArea[subdomain1]
	   << ",   and of " << subdomain2 << ":  " << subdomainArea[subdomain2]; cout.flush();
  cout << endl << " Subdomain Lengths of " << subdomain1 << ": " <<  subdomainLength[subdomain1]
	   << "  and of " << subdomain2 << ":  " << subdomainLength[subdomain2]; cout.flush();
  cout << endl << " Subdomain Weights of " << subdomain1 << ": " <<  subdomainWeight[subdomain1]
	   << "  and of " << subdomain2 << ":  " << subdomainWeight[subdomain2]; cout.flush();
  cout << endl << "--------------------------------------------------" << endl;
#endif
#undef __VERBAL

  //--------------------------------------
  // the main partiton routine is done
  //--------------------------------------
  return 0;
}// partition()



//---------------------------------------------------------------------------*
//                                                                           *
//     END OF THE MAIN DECOMPOSITION ROUTINES                                *
//                                                                           *
//---------------------------------------------------------------------------*


//---------------------------------------------------------------------------*
//     THE DELAUNAY TRIANGULATION                                            *
//---------------------------------------------------------------------------*
// the Delaunay triangulation routine
// triangulates the work subdomain using Shewchuck's Triangle
// gets a string of flags passed to Triangle
//---------------------------------------------------------------------------
int   
madd::triangulateDomain(char *inFlags)
{

#ifdef  __VERBAL2
  cout << endl << "      Triangulating. Creating the input structure..."; cout.flush();
#endif
  struct triangulateio in, out;
  //--------------------------------------
  // initialize input structure
  // get input points and segments. 
  LInt *fromPoint;
  in.numberofpoints  = createPointsSegms(&in.pointlist, &in.segmentlist, &fromPoint);
  in.numberofsegments = workNoOfSegms;
  // get the holes
  in.numberofholes = workNoOfHoles;
  in.holelist = workHolesList;

#ifdef  __VERBAL2
  cout << endl << " Input Points: " << in.numberofpoints
	   << " Input Segments: " << in.numberofsegments
	   << " Input Holes: " << in.numberofholes << endl; cout.flush();
#endif
  //--------------------------------------
  // zero the output structure
  out.numberofpoints =  in.numberofpoints;
  out.pointlist = in.pointlist;
  out.numberofsegments = in.numberofsegments;
  out.segmentlist  = in.segmentlist; 
  out.holelist =in.holelist;

  out.regionlist = in.regionlist = NULL;
  out.edgelist = in.edgelist = NULL;             
  out.segmentmarkerlist = in.segmentmarkerlist = NULL;
  out.edgemarkerlist = in.edgemarkerlist = NULL;       
  out.normlist = in.normlist = NULL;             
  out.trianglelist = in.trianglelist = NULL;                                 
  out.triangleattributelist = in.triangleattributelist = NULL;                        
  out.trianglearealist = in.trianglearealist = NULL;                             
  out.neighborlist = in.neighborlist = NULL;
  out.pointattributelist = in.pointattributelist = NULL;
  out.pointmarkerlist = in.pointmarkerlist = NULL;
  out.normlist =  in.normlist = NULL;
                              
  out.numberofpoints =  in.numberofpoints;
  out.numberofholes = in.numberofholes;

  out.numberofpointattributes = in.numberofpointattributes = 0;
  out.numberofregions = in.numberofregions = 0;
  out.numberofedges =  in.numberofedges = 0;        
  out.numberoftriangles =in.numberoftriangles =  0;                                  
  out.numberofcorners = in.numberofcorners = 0;                                 
  out.numberoftriangleattributes =in.numberoftriangleattributes =  0;   
  out.numberofedges = in.numberofedges = 0;                    /* Out only */
                   
  out.trianglelist = NULL;            
  out.pointlist = NULL;

  //--------------------------------------
  // initialize the Triangle flags
  char flags[20];
  for (int j=0; j < 20; j++)
	flags[j] = ' ';

  flags[0] = 'p';  // triangulate poly file
  flags[1] = 'P';  // no output poly file
  flags[2] = 'N';  // no output  node 
  flags[3] = 'Q';  // quiet
  //     flags[3] = 'V';  // verbose
  flags[4] = 'B';  // no boundary markers
  flags[5] = 'z';  // start from zero
  flags[6] = 'n';  // create neighborlist
  //  flags[7] = 'V';  // more verbose
  //  flags[8] = 'V';  // more verbose

  int startFlags = 9;
  int j;
  for (j = startFlags; j < 20; j++)
	flags[j] = 0;

  if (inFlags != NULL) {
	for (j=0; j < 10 &&  inFlags[j] != 0; j++) {
	  flags[startFlags++] = inFlags[j];
	  if (inFlags[j] == 'q' || inFlags[j] == 'a')
		flags[2] = ' ';
	}//for
  }//if
  // Triangle flags are initialized
  //--------------------------------------

#ifdef  __VERBAL2
  cout << endl << "      Triangle started..."; cout.flush();
#endif

  //--------------------------------------
  // triangulate
  triangulate(flags, &in, &out, (struct triangulateio *) NULL);
  //--------------------------------------

#ifdef  __VERBAL2
  cout << "finished."; cout.flush();
#endif

  //--------------------------------------
  // get the results and 
  // check if we got the triangulation
  triaNeighborList  = out.neighborlist;
  noOfTria = out.numberoftriangles;
  theTriaList  = out.trianglelist;
  if (noOfTria <= 0)
	ErrorExit((char *)"No Triangles in the list!"); 
  if (out.trianglelist == NULL)
	InsufficientMem((char *)" running triangle");
  if (flags[2] == ' ' && out.pointlist == NULL)
	InsufficientMem((char *)" running triangle");
  //--------------------------------------


#ifdef  __VERBAL2
  cout << endl << "     getting the triangles..."; cout.flush();
  cout << endl << "     assign the actual point values to the triangles..."; cout.flush();
#endif
  //--------------------------------------
  // assign the actual point values to the triangles
  for (int i = 0; i < noOfTria; i++){
	/*
	  cout << endl << "tria: " << i << " -> " <<  tria_point1(i) << ",  "
	  <<  tria_point2(i) << ",  "<<  tria_point3(i); cout.flush();
	  cout << endl << "      original points: " << " -> " <<  fromPoint[tria_point1(i)] << ",  "
	  <<  fromPoint[tria_point2(i)] << ",  "<<  fromPoint[tria_point3(i)]; cout.flush();
	*/
	tria_point1(i) = fromPoint[tria_point1(i)];
	tria_point2(i) = fromPoint[tria_point2(i)];
	tria_point3(i) = fromPoint[tria_point3(i)];
  } 
  // actual point values assigned
  //--------------------------------------
#ifdef  __VERBAL2
  cout << endl << "     done.  free the auxiliary structures..."; cout.flush();
#endif
  //--------------------------------------
  // free the auxiliary structures  
  free(in.segmentlist);
  free(in.pointlist);
  free(fromPoint);
  //--------------------------------------


//--------------------------------------
// start  __DEBUG_MODE
#ifdef  __DEBUG_MODE
  cout << endl << "      Checking neighbor list...."; cout.flush();
  for (LInt j=0; j < noOfTria; j++) {
#ifdef  __VERBAL2
	cout << endl << "neighborList: " 
	  "(" <<  tria_point1(j) << ", " << tria_point2(j) << ", " << tria_point3(j) << ") -> (" 
		 << neighborTria(j, 0) << ", " << neighborTria(j, 1) << ", " << neighborTria(j, 2) << ")"; cout.flush(); 
#endif
	for (int k = 0; k < 3; k++){
	  LInt p1 = tria_point1(j);
	  LInt p2 = tria_point2(j);
	  if (k == 1)
		p2 = tria_point3(j);
	  if (k == 0)
		p1 = tria_point3(j);

	  LInt nextTria = neighborTria(j, k);
	  LInt inPoints = 0;

	  if (nextTria < 0){
		// check the triangle connectivity
		findWorkSegm(p1, p2);
		continue;
	  }

	  if (tria_point1(nextTria) == p1 ||tria_point1(nextTria) == p2)
		inPoints++;
	  if (tria_point2(nextTria) == p1 ||tria_point2(nextTria) == p2)
		inPoints++;
	  if (tria_point3(nextTria) == p1 ||tria_point3(nextTria) == p2)
		inPoints++;

	  if (inPoints != 2){
		cout << endl << " -- Error in triaNeighborList. : "; 
		cout << k << "->" << nextTria << " =  (" 
			 << tria_point1(nextTria) << ", " << tria_point2(nextTria) << ", " 
			 << tria_point3(nextTria) << ")." << endl;  cout.flush();
		exit(1);
	  }
	}
  }
#endif 
//end of __DEBUG_MODE
//--------------------------------------

//--------------------------------------
#ifdef  __VERBAL2
  cout << endl << "      Triangulation done."; cout.flush();
#endif
//--------------------------------------
  return 0;
}//triangulateDomain()


//---------------------------------------------------------------------------*
//     END OF THE DELAUNAY TRIANGULATION                                     *
//---------------------------------------------------------------------------*




//---------------------------------------------------------------------------*
//                                                                           *
//     START OF THE MADD ROUTINES                                            *
//                                                                           *
//---------------------------------------------------------------------------*

//---------------------------------------------------------------------------
// the Metis function descriptions
typedef int idxtype;
extern "C" {
  void METIS_PartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
}

extern "C" {
  void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
}

extern "C" {
  void METIS_PartGraphVKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
}

extern "C" {
  void METIS_WPartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *,
								 int *, int *, int *, float *, int *, int *, idxtype *);
}
// end of the Metis function descriprion
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// variable definitions
#define  MAX_IMBALANCE_NORMAL 0.7   // this is accepted provided that the
#define  MAX_IMBALANCE_DIFF 0.15    // imbalance difference is not greater than this one
#define  MAX_COEF 64                // returned when a separator is not accepted

DFloat  MAX_IMBALANCE_LIMIT = 0.75; // the maximum allowed imbalance
int  MAX_SMOOTH_DEPTH = 1;          // the smooth depth (no of recursions)
int  SMOOTHING_POINTS = 11;         // the number of examined points for smoothing (5+5+1)
int MAX_MADD_NODE_WEIGHT = 1000;    // the Metis max node weight
int MAX_MADD_EDGE_WEIGHT = 1600;    // the Metis max edge weight

DFloat theBasicCoeff[3] = {1.0, 0.65, 0.75}; // weight coefficients for not basic, basic, middle point
DFloat theSmoothCoeff[3] = {1.0, 0.5, 0.7}; // smoothing coefficients for not basic, basic, middle point

LLint noOfNodes;    // the number of nodes in the initial graph 
LLint* nodeGoTo;    // the collapse list of the nodes
DFloat* nodeArea;   // the area for each node
LInt * nodeMetisAA; // the pointer to the metis node
DFloat totalNodeArea;// the toatal area of the subdomain
  
LLint noOfEdges;   // the number of edges in the initial graph 
LLint* graphEdge;  // the list of edges
DFloat* edgeLength; // the lentgh for each edge
LLint GRAPH_EDGE_SIZE; // the allocated size of the edges structure

// Metis lists
idxtype *theNodesIndex, *theEdgesIndex, *theNodesWeight, *theEdgesWeight, *theNodesColour;   
LInt metisActiveNodes = 0;

// connectivity lists used for checking the Metis compontents' connectivity
int *connMetisColour;
uchar *nodeVisited;

// Arrays for keeping track of the new separator points	 
LInt seperatorPoint1[400];
LInt seperatorPoint2[400];
LInt sepTriangleIndex, noOfSeparatorTriangles;
// end of variable definitions
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// shortcuts
#define node(j, i) (3*j + i) 
#define edgeNode1(j) graphEdge[((j) << 1)] 
#define edgeNode2(j) graphEdge[((j) << 1)+1] 
#define metisColour(someNode) (theNodesColour[nodeMetisAA[getEndNode(someNode)]])
#define metisWeight(someNode) (theNodesWeight[nodeMetisAA[getEndNode(someNode)]])
#define connColour(someNode) (connMetisColour[nodeMetisAA[getEndNode(someNode)]])
#define addSeparatorPoint1(SepIndex, SepPoint) (seperatorPoint1[SepIndex] = SepPoint)
#define addSeparatorPoint2(SepIndex, SepPoint) (seperatorPoint2[SepIndex] = SepPoint)
// end of shortcuts
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------
//  allocate the initial graph nodes
//---------------------------------------------------------------------------
void
madd::allocateGraphNodes()
{
  noOfNodes = 3 * noOfTria; 
  nodeGoTo = (LLint*)malloc(noOfNodes*sizeof(LLint)+1);
  nodeArea = (DFloat*)malloc(noOfNodes*sizeof(DFloat)+1);
  // nodeWeight = (DFloat*)malloc(noOfNodes*sizeof(DFloat)+1);
  nodeMetisAA = (LInt*)malloc(noOfNodes*sizeof(LInt)+1);

  if (nodeGoTo == NULL || nodeArea == NULL || nodeMetisAA == NULL)
	InsufficientMem((char *)"madd node list");

  for (LLint i =0; i < noOfNodes; i++){
	nodeGoTo[i] = -1;
	nodeMetisAA[i] = -1;
	nodeArea[i] = 0;
  }
}//allocateGraphNodes()


//---------------------------------------------------------------------------
//  allocate the initial graph edges
//---------------------------------------------------------------------------
void
madd::allocateGraphEdges()
{
  noOfEdges = 0;
  GRAPH_EDGE_SIZE = 3 * noOfTria + 100; 
  graphEdge = (LLint*)malloc(GRAPH_EDGE_SIZE*2*sizeof(LLint)+1);
  edgeLength = (DFloat*)malloc(GRAPH_EDGE_SIZE*sizeof(DFloat)+1);

  if (graphEdge == NULL || edgeLength == NULL)
	InsufficientMem((char *)"madd edges list");
}//allocateGraphEdges()


//---------------------------------------------------------------------------
//  reallocate the initial graph edges
//---------------------------------------------------------------------------
void
madd::reallocateGraphEdges()
{
  GRAPH_EDGE_SIZE = noOfEdges + 200; 
  graphEdge = (LLint*)realloc(graphEdge, GRAPH_EDGE_SIZE*2*sizeof(LLint)+1);
  edgeLength = (DFloat*)realloc(edgeLength, GRAPH_EDGE_SIZE*sizeof(DFloat)+1);

  if (graphEdge == NULL || edgeLength == NULL)
	InsufficientMem((char *)"madd edges list");

}//reallocateGraphEdges()


//---------------------------------------------------------------------------
//  returns the end node pointed from the startNode
//  it collapses the route for efficiency
//---------------------------------------------------------------------------
LLint
madd::getEndNode(LLint startNode)
{

  // find end node
  LLint currNode = startNode;
  while (currNode < noOfNodes && nodeGoTo[currNode] >= 0){
	currNode = nodeGoTo[currNode]; 
	//	cout << endl << currNode; cout.flush();
  }
  LLint finalNode = currNode; 
  if (currNode >= noOfNodes)
	ErrorExit((char *)" Node out of link");
  //	cout << endl << "collapse the links..."; cout.flush();
  //collapse the links
  currNode = startNode;
  LLint nextNode; 
  while (currNode != finalNode){
	nextNode = 	nodeGoTo[currNode];
	nodeGoTo[currNode] = finalNode;
	currNode = nextNode;
  } 
  //	cout << endl << "Get node Done."; cout.flush();
  return finalNode;
}//getEndNode()


//---------------------------------------------------------------------------
//  returns the node in linkTriangle that is adjacent to toTriangle
//---------------------------------------------------------------------------
LLint
madd::getLinkNode(LInt linkTriangle, LInt toTriangle)
{
  LLint addNo = -1;
  if (neighborTria(linkTriangle, 0) == toTriangle)
	addNo = 0;
  else if (neighborTria(linkTriangle, 1) == toTriangle)
	addNo = 1;
  else if (neighborTria(linkTriangle, 2) == toTriangle)
	addNo = 2;
  
  if (addNo >= 0)
	return (3*linkTriangle + addNo);
  else{
	ErrorExit((char *)" Wrong link in Triangles while creating MADD graph");
  }
  return -1;
}//getLinkNode()


//---------------------------------------------------------------------------
//  collapses the node fromNode to toNode
//---------------------------------------------------------------------------
void
madd::collapseNodeTo(LLint fromNode, LLint toNode)
{
  LLint node1 = getEndNode(fromNode);
  LLint node2 = getEndNode(toNode);
  if (node1 != node2){
	nodeGoTo[node1] = node2;
	nodeArea[node1] +=  nodeArea[node2];
	nodeArea[node2] = 0;
  }
}//collapseNodeTo


//---------------------------------------------------------------------------
// insert an edge between node1, node2, with weight ...weight
//---------------------------------------------------------------------------
void
madd::insertGraphEdge(LLint node1, LLint node2, DFloat weight) 
{
  LLint endnode1 =  getEndNode(node1);
  LLint endnode2 =  getEndNode(node2);

  if (endnode1 != endnode2){
	edgeNode1(noOfEdges) = endnode1;  
	edgeNode2(noOfEdges) = endnode2;
	edgeLength[noOfEdges] = weight; 
	noOfEdges++;
  }
}//insertGraphEdge()


//---------------------------------------------------------------------------
// the main partiton routine
// returns 0 if succesful
//---------------------------------------------------------------------------
int
madd::MADDpartition()
{
  createMADDGraph();
  // create the Metis Strucur form the MADD graph
  createMetisStructure();

  // partition the Metis graph
  partitionGraph();

  // insert the separators
  int failures = 0;
  //  cout << "\n MetisActiveNodes  = " << metisActiveNodes;
  // tru to insert the seprators
  while ((!insertSeparators()) && (failures < 5)){
#ifdef __VERBAL2
	cout << "\n\n We have a failure!\n"; cout.flush();
#endif
	// if it failed, repartition and try again
	createMetisStructure();
	partitionGraph();
	failures++;
  }// while

  if (failures >= 5){
	ErrorExit((char *)"\n Separators meet at the same point. Unable to partiton!");
  }

  // the geometry has been partitioned 
  // clear the Metis lists
  free(theEdgesIndex);
  theEdgesIndex= NULL;
  free(theNodesColour);
  theNodesColour = NULL;
  // clear the madd graph lists
  free(nodeGoTo);
  nodeGoTo = NULL;
  free(nodeArea);
  nodeArea = NULL;
  free(nodeMetisAA);
  nodeMetisAA = NULL;
  free(graphEdge);
  graphEdge = NULL;
  free(edgeLength);
  edgeLength = NULL;
  madSepEnd = noOfSegms;

#ifdef __VERBAL2
  cout << endl << endl << "== MADD is done ==" << endl; cout.flush();
#endif
  return 0;
}//MADDpartition()



//---------------------------------------------------------------------------
// creates the MADD graph for the Delaunay triangulation
// returns 0 if succesful
//---------------------------------------------------------------------------
int
madd::createMADDGraph()
{
  LInt exp1, exp2, exop;
  // LInt segm;
  //  LInt endp1, endp2, opp;
  LInt p1, p2, p3;                  // the three vertices of a triangle
  DFloat *point1, *point2, *point3; // the three vertices of a triangle
  DFloat *expoint1, *expoint2;      // two boundary points
  DFloat center[2];                 // the center of a triangle
  DFloat angle1, angle2, angle3;    // three angles
  DFloat eAngle1, eAngle2;          // two boundary angles
  DFloat sepCoeff;                  // the coefficient of the separator
  LLint linkTria;                   // the adjacent triangle
  LLint linkNode;                   // the adjacent node
  LLint node1, node2, node3;        // the three nodes of the triangle
  DFloat radius, area;              // the radius and the area of the triangle

  DFloat edgelng = 0;               // total edge length is 0
  totalNodeArea = 0;                // total domain area is 0

  // allocate initial graph structure
  allocateGraphNodes();
  allocateGraphEdges();


#ifdef  __VERBAL
  cout << endl << "     Create MADD Graph..."; cout.flush();
#endif
  //----------------------------------------
  //  create the MADD graph 
  //----------------------------------------
  for (LInt thisTria=0; thisTria < noOfTria; thisTria++) {
	// the trhee triangle vertices
	p1 = tria_point1(thisTria); 
	p2 = tria_point2(thisTria); 
	p3 = tria_point3(thisTria);
	point1 = point(p1);
	point2 = point(p2);
	point3 = point(p3);
	// the three nodes of this triangle
	node1 = node(thisTria, 0);
	node2 = node(thisTria, 1);
	node3 = node(thisTria, 2);
	// its circumcenter, and radius 
	circumCenter(point1, point2, point3, center);
	radius = distance(point1, center);
	// the three internal angles, from the circumcenter
	angle1 = angle(point2, center, point3);
	angle2 = angle(point1, center, point3);
	angle3 = angle(point1, center, point2);
	// the total triangle are
	totalNodeArea += (area = areaOfTriangle(point1, point2, point3));
#ifdef  __VERBAL2
	cout << endl << "     check total quality..."; cout.flush();
#endif
	bool accepted = true; // assume the triangle is good
	// count external edges
	int extEdges = 0;
	if (triaNeighborList[node1] < 0)
	  extEdges++;
	if (triaNeighborList[node2] < 0)
	  extEdges++;
	if (triaNeighborList[node3] < 0)
	  extEdges++;

	if (!isnormal(x(center)) || !isnormal(y(center)))
	  accepted = false; // center is not valid
	if ( (angle1+angle2+angle3) < TWO_PI)
	  accepted = false; // center is outside the triangle
	if ( extEdges > 1)
	  accepted = false; //  two external edges
	if (angle1 >= PI_MINUS || angle2 >= PI_MINUS || angle3 >= PI_MINUS) 
	  accepted = false; //  center is on the edges

	if ( !accepted){
#ifdef  __VERBAL2
	  cout << " not accepted."; cout.flush();
	  cout << endl << "      collapse nodes..."; cout.flush();
#endif
	  collapseNodeTo(node2, node1);   
	  collapseNodeTo(node3, node1);
	  nodeArea[node1] +=  area;
	}

#ifdef  __VERBAL2
	cout << endl << "      try to insert inner graph edges...";    cout.flush();
#endif
	if (noOfEdges + 8 >= GRAPH_EDGE_SIZE)
	  reallocateGraphEdges();


	if (accepted){
#ifdef  __VERBAL2
	  cout << " total quality accepted."; cout.flush();
#endif
	  // assign the areas to the three nodes
	  nodeArea[node1] += areaOfTriangle(point2, center, point3) ;
	  nodeArea[node2] += areaOfTriangle(point1, center, point3) ;
	  nodeArea[node3] += areaOfTriangle(point1, center, point2) ;

	  // check if radius is too close to an edge
	  if (angle1 > M_PI - 0.4)
		collapseNodeTo(node1, node2);
	  if (angle2 > M_PI - 0.4)
		collapseNodeTo(node2, node3);
	  if (angle3 > M_PI - 0.4)
		collapseNodeTo(node3, node1);

	  // check if angles are too small
	  if (angle1  < PHI)
		collapseNodeTo(node1, node2);
	  if (angle2  < PHI)
		collapseNodeTo(node2, node3);
	  if (angle3  < PHI)
		collapseNodeTo(node3, node1);
	
	  if ((sepCoeff =  basicCoef1(center, p3)) < MAX_COEF){
	    // inner separator 1 is accepted
		edgelng = radius * sepCoeff;
		if (isnan(edgelng)){
		  cout << endl << " edgelng == NAN,  radius = " << radius << ",  sepCoeff = " << sepCoeff; cout.flush();
		  ErrorExit((char *)"Stop execution");
		} 
		insertGraphEdge(node1, node2, edgelng);
  	    // check to see how inner separator looks
		//		NewSegm(newpp, p3, 0);
	  }
	  else {
	    // inner separator 1 is not accepted
		collapseNodeTo(node2, node1);   		
	  }

	  if ((sepCoeff =  basicCoef1(center, p1)) < MAX_COEF){
	    // inner separator 2 is accepted
		edgelng = radius * sepCoeff; 
		if (isnan(edgelng)){
		  cout << endl << " edgelng == NAN,  radius = " << radius << ",  sepCoeff = " << sepCoeff; cout.flush();
		  ErrorExit((char *)"Stop execution");
		} 
		insertGraphEdge(node2, node3, edgelng);
  	    // check to see how inner separator looks
		//	NewSegm(newpp, p1, 0);
	  }
	  else {
	    // inner separator 2 is not accepted
		collapseNodeTo(node3, node2);   		
	  }

	  if ((sepCoeff =  basicCoef1(center, p2)) < MAX_COEF){
	    // inner separator 1 is accepted
		edgelng = radius * sepCoeff; 
		if (isnan(edgelng)){
		  cout << endl << " edgelng == NAN,  radius = " << radius << ",  sepCoeff = " << sepCoeff; cout.flush();
		  ErrorExit((char *)"Stop execution");
		} 
		insertGraphEdge(node1, node3, edgelng);
  	    // check to see how inner separator look
		// NewSegm(newpp, p2, 0);
	  }
	  else {
	    // inner separator 3 is not accepted
		collapseNodeTo(node3, node1);   		
	  }

	}// if accepted

#ifdef  __VERBAL2
	cout << "done."; cout.flush();
#endif

#ifdef __DEBUG_MODE
	// check the triangle connectivity
    if (triaNeighborList[node1] < 0)
	  findWorkSegm(p2, p3);
    if (triaNeighborList[node2] < 0)
	  findWorkSegm(p1, p3);
    if (triaNeighborList[node3] < 0)
	  findWorkSegm(p1, p2);
#endif

#ifdef  __VERBAL2
	cout << endl << "      insert outer edges...";    cout.flush();
#endif
	if ((linkTria = triaNeighborList[node1]) > thisTria){
	  if ((sepCoeff =  basicCoef2(p2, p3)) < MAX_COEF){
	    // outer separator 1 is accepted
		edgelng = distance(point2, point3) * sepCoeff;
		// check numerical validity for edgelng
		if (isnan(edgelng)){
		  cout << endl << " edgelng == NAN,  distance = " << distance(point2, point3) 
			   << ",  sepCoeff = " << sepCoeff; cout.flush();
		  ErrorExit((char *)"Stop execution");
		} 
		insertGraphEdge(node1, getLinkNode(linkTria, thisTria), edgelng);
	  }
	  else
	    // outer separator 1 is not accepted
		collapseNodeTo(node1, getLinkNode(linkTria, thisTria));   			  
	}
	
	if ((linkTria = triaNeighborList[node2]) > thisTria){
	  if ((sepCoeff =  basicCoef2(p1, p3)) < MAX_COEF){
	    // outer separator 2 is accepted
		edgelng = distance(point1, point3) * sepCoeff;
		// check numerical validity for edgelng
		if (isnan(edgelng)){
		  cout << endl << " edgelng == NAN,  distance = " << distance(point1, point3) 
			   << ",  sepCoeff = " << sepCoeff; cout.flush();
		  ErrorExit((char *)"Stop execution");
		} 
		insertGraphEdge(node2, getLinkNode(linkTria, thisTria), edgelng);
	  }
	  else
	    // outer separator 2 is not accepted
		collapseNodeTo(node2, getLinkNode(linkTria, thisTria));  
	}  	

	if ((linkTria = triaNeighborList[node3]) > thisTria){
	  if ((sepCoeff =  basicCoef2(p1, p2)) < MAX_COEF){
	    // outer separator 3 is accepted
		edgelng = distance(point1, point2) * sepCoeff;
		// check numerical validity for edgelng
		if (isnan(edgelng)){
		  cout << endl << " edgelng == NAN,  distance = " << distance(point1, point2) 
			   << ",  sepCoeff = " << sepCoeff; cout.flush();
		  ErrorExit((char *)"Stop execution");
		} 
		insertGraphEdge(node3, getLinkNode(linkTria, thisTria), edgelng);
	  }
	  else
	    // outer separator 3 is not accepted
		collapseNodeTo(node3, getLinkNode(linkTria, thisTria));  
   	}

	
  }// for all triangles 
  //----------------------------------------
  //  the MADD graph is created
  //----------------------------------------
  return 0;
}//createMADDGraph()




//-------------------------------------------------------------------------------
// creates the Metis graph structure from the MADD graph
// returns 0 if succesful
//-------------------------------------------------------------------------------
int
madd::createMetisStructure()
{
#ifdef  __VERBAL2
  cout << endl << "Nodes  : " << noOfNodes << "    Edges  : " << noOfEdges; cout.flush();
  cout << endl << "      Creating Metis structure"; cout.flush();
#endif 
  // clear the Metis lists
  free(theNodesIndex);
  theNodesIndex = NULL;
  free(theEdgesIndex);
  theEdgesIndex = NULL;
  free(theNodesWeight);
  theNodesWeight = NULL;
  free(theEdgesWeight);
  theEdgesWeight = NULL;
  free(theNodesColour);
  theNodesColour = NULL;



#ifdef  __VERBAL2
  cout << endl << "      Check Active Nodes..."; cout.flush();
#endif
  //----------------------------------------
  //  get active nodes
  //  store a pointer to the Metis nodes
  //----------------------------------------
  DFloat maxNodeArea = 0; // the maximum node area
  metisActiveNodes = 0;
  for (LLint thisNode = 0; thisNode < noOfNodes; thisNode++){
	LInt endNode = thisNode;
	if (nodeGoTo[thisNode] < 0){
	  // node is active
	  nodeMetisAA[thisNode] = metisActiveNodes++; //  store a pointer to the Metis nodes
	}
	else{
	  // node is not active
	  nodeMetisAA[thisNode] = -1;
	  // add the area to the End Node and collapse the link path
	  endNode = getEndNode(thisNode);
	  nodeArea[endNode] += nodeArea[thisNode];	 
	  nodeArea[thisNode] = 0;
	}
	if (maxNodeArea < nodeArea[endNode])
	  maxNodeArea = nodeArea[endNode]; // this is the max node area so far
  }// for

#ifdef  __VERBAL2
  cout << endl << "      Check Active edges..."; cout.flush();
#endif

  //----------------------------------------
  //  Allocate the degrees of the active nodes
  int* nodeDegree = (int*)malloc(metisActiveNodes*sizeof(int)+1);  
  if (nodeDegree == NULL)
	InsufficientMem((char *)"node degree list");
  //----------------------------------------
  //  get active edges,
  //  assign the degree of the active nodes
  //----------------------------------------
  for (int i=0; i < metisActiveNodes; i++)
	nodeDegree[i] = 0;

  DFloat totalEdgeLength = 0; // total edge length is 0
  DFloat maxEdgeLength = 0;   // maximum edge length is 0
  LInt activeEdges = 0;       // the active edges

  for (LLint thisEdge = 0; thisEdge < noOfEdges; thisEdge++){
	// go through all MADD edges
#ifdef  __VERBAL3
	cout << endl << "Checking edge: " << thisEdge 
		 << "    node1:" <<  edgeNode1(thisEdge)
		 << "    node2:" <<  edgeNode2(thisEdge); cout.flush();
#endif
	// get the two end nodes of the edge
	LLint node1 = getEndNode(edgeNode1(thisEdge));
	LLint node2 = getEndNode(edgeNode2(thisEdge));

	if (node1 != node2){
	  // they are not the same
	  // we have an active edge
	  activeEdges++;
	  totalEdgeLength += edgeLength[thisEdge];
#ifdef  __VERBAL3
	  cout << endl << "edgeLength[thisEdge]: " << edgeLength[thisEdge] 
	  	   << "   totalEdgeLength: " << totalEdgeLength; cout.flush();
#endif
	  if (isnan(totalEdgeLength))
		ErrorExit((char *)"isnan(totalEdgeLength)");
	  // update the degree of the two nodes
	  nodeDegree[nodeMetisAA[node1]]++;
	  nodeDegree[nodeMetisAA[node2]]++;

	  if (maxEdgeLength < edgeLength[thisEdge])
		maxEdgeLength = edgeLength[thisEdge]; // the new maximum edge length
	}// if active edge
  }//for all edges
  //----------------------------------------
  //   active nodes and edges located
  //----------------------------------------

  if (activeEdges < 1)
	ErrorExit((char *)"Not sufficient refinement");

  //----------------------------------------
  //  Allocate the Metis structures
  // allocate the start position of the connectivity list for Metis 
  if ((theNodesIndex = (idxtype *)malloc((metisActiveNodes+2) * sizeof(idxtype))) == NULL)
	InsufficientMem((char *)"Creating Graph for Metis");
  // allocate the node weigth list for Metis 
  if ((theNodesWeight = (idxtype *)malloc((metisActiveNodes+1) * sizeof(idxtype))) == NULL)
	InsufficientMem((char *)"Writing Graph");
  // allocate Metis node colour list
  if ((theNodesColour = (idxtype *)malloc((metisActiveNodes+1) * sizeof(idxtype))) == NULL)
	InsufficientMem((char *)"MADD Colour Allocating");
  // allocate the Metis edges
  if ((theEdgesIndex = (idxtype *)malloc(( activeEdges*2+2) * sizeof(idxtype))) == NULL)
	InsufficientMem((char *)"theEdgesIndex");
  if ((theEdgesWeight = (idxtype *)malloc((activeEdges*2+2) * sizeof(idxtype))) == NULL)
	InsufficientMem((char *)"theEdgesWeight");
  //  the structures are allocated 
  //----------------------------------------

#ifdef  __VERBAL2
  cout  << endl << " Active Nodes = " << metisActiveNodes << "  Total Area:" <<  totalNodeArea; cout.flush();
  cout << endl << "  Active Edges = " << activeEdges << "    Total Length = " << totalEdgeLength ; cout.flush();
  //  DFloat areaToLength = sqrt(totalNodeArea) / totalEdgeLength;
  cout << endl 
	   <<  ",  MAX_MADD_NODE_WEIGHT: " << MAX_MADD_NODE_WEIGHT
	   <<  ",  maxNodeArea: " << maxNodeArea
	   <<  ",  maxEdgeLength: " << maxEdgeLength; cout.flush();
#endif

  //--------------------------------------
  // coefficent definition for computing the metis weights
  DFloat nodeCoef = MAX_MADD_NODE_WEIGHT /  maxNodeArea; 
  // the edge should have relativly high weight for metis
  DFloat edgeCoef =  max(sqrt( nodeCoef * 1.4),  MAX_MADD_EDGE_WEIGHT / maxEdgeLength);
  //--------------------------------------

#ifdef  __VERBAL2
  cout << endl << "      Create the start position list and the weight list for Metis..."; cout.flush();
#endif
  //--------------------------------------
  // assign to the nodes index the start position of the connectivity   
  idxtype startEdge = 0;
  theNodesIndex[0] = 0;
  for (int i=0; i < metisActiveNodes;){
	startEdge += nodeDegree[i];
	theNodesIndex[(++i)] = startEdge;
  }
  //  theNodesIndex is ready
  //--------------------------------------

#ifdef  __VERBAL2
  cout << endl << "      create Metis node weights..."; cout.flush();
#endif
  //--------------------------------------
  // assign the Metis nodes weights   
  for (LLint thisNode = 0; thisNode < noOfNodes; thisNode++){
	int metisAA;
	if ((metisAA = nodeMetisAA[thisNode]) >= 0){ // if node active
	  theNodesWeight[metisAA] = (int)( nodeArea[thisNode] *  nodeCoef);
#ifdef  __VERBAL3
	  cout << endl << "Node weight : " << nodeArea[thisNode]  
		   << " -> " << theNodesWeight[metisAA]; cout.flush();
#endif
	}// if node active
  }// for all nodes
  // Metis nodes weights are ready   
  //--------------------------------------
						 
#ifdef  __VERBAL2
  cout << endl << "      create Metis edges..."; cout.flush();
#endif
  //--------------------------------------
  // create the Metis edges
  // go through all the MADD edges
  for (LLint thisEdge = 0; thisEdge < noOfEdges; thisEdge++){
#ifdef  __VERBAL3
	cout << endl << "Checking edge: " << thisEdge 
		 << "    node1:" <<  edgeNode1(thisEdge)
		 << "    node2:" <<  edgeNode2(thisEdge); cout.flush();
#endif
	// the two end nodes of this edge
	LLint node1 = getEndNode(edgeNode1(thisEdge));
	LLint node2 = getEndNode(edgeNode2(thisEdge));
#ifdef  __VERBAL3
	cout << endl << "  Link edge: " 
		 << "    node1:" <<  node1
		 << "    node2:" <<  node2; cout.flush();
#endif
	if (node1 != node2){
	  //--------------------------------------
	  // we have an active edge
	  // insert it to Metis edges
	  idxtype metisNode1 = nodeMetisAA[node1];
	  idxtype metisNode2 = nodeMetisAA[node2];
	  // find the position to be inserted
	  idxtype pos1 = theNodesIndex[metisNode1 + 1] - nodeDegree[metisNode1];
	  idxtype pos2 = theNodesIndex[metisNode2 + 1] - nodeDegree[metisNode2];
	  nodeDegree[metisNode1]--;
	  nodeDegree[metisNode2]--;
#ifdef  __VERBAL3
	  cout << endl << "  Metis " 
		   << "    node1:" <<  metisNode1
		   << "    node2:" <<  metisNode2; cout.flush();
	  cout << endl << "    pos1:" <<  pos1
		   << "    pos2:" <<  pos2; cout.flush();
	  cout << endl << "  Degree " 
		   << "    node1:" <<  nodeDegree[metisNode1]
		   << "    node2:" <<  nodeDegree[metisNode2]; cout.flush();
#endif
	  // insert the two nodes
	  theEdgesIndex[pos1] = metisNode2;
	  theEdgesIndex[pos2] = metisNode1;
	  theEdgesWeight[pos1] = theEdgesWeight[pos2] = (int)(edgeLength[thisEdge] * edgeCoef);
	  // Metis edge created
	  //--------------------------------------
#ifdef  __VERBAL3
	  cout << endl << "The edge weight  " << edgeLength[thisEdge] << " * " 
		   << edgeCoef << "  = " << (int)(edgeLength[thisEdge] * edgeCoef); cout.flush();
	  //	  cout << endl << "Edge weight : " <<  theEdgesWeight[pos1]; cout.flush();
	  if ( theEdgesWeight[pos1] <= 0)
		cout << endl << "   !! theEdgesWeight <= 0 !!" << endl; cout.flush();
	  //		ErrorExit((char *)"theEdgesWeight is <= 0 ");
	  cout << "done "; cout.flush();
#endif
	}// if active edge
  }// for all MADD edges
  // the Metis edges are created
  //--------------------------------------

  // The Metis structure is ready
  // get rid of the auxiliary structures  
  free(nodeDegree);
  nodeDegree = NULL;

  return 0;
}// createMetieStructure()





//---------------------------------------------------------------------------
// partition the Metis graph using ...Metis
// restores the connectivity and returns two connected components
// returns true if succesful
//---------------------------------------------------------------------------
bool
madd::partitionGraph()
{
  int opts[5] = {0, 0, 0, 0, 0}; // default options for MEtis
  // int opts[5] = {1, 2, 1, 2, 0}; // for METIS_PartGraphKway
  // int opts[5] = {1, 3, 1, 3, 0}; // for METIS_PartGraphVKway
  // int opts[5] = {1, 2, 1, 1, 0}; // 
  // float partWeights[3] = {0.55, 0.45, 0.0}; // for Metis balance, not used
  int theWeightsFlag = 3;
  int startFrom = 0;
  int theParts = 2;
  int edgecut;
  int AA = metisActiveNodes;

#ifdef __VERBAL2
  cout << "calling Metis..."; cout.flush();
#endif
  METIS_PartGraphRecursive(&AA, theNodesIndex, theEdgesIndex, 
						   theNodesWeight, theEdgesWeight, 
						   &theWeightsFlag, &startFrom, &theParts, opts, &edgecut,
						   theNodesColour);
  // other partition options		
  /*     
		 METIS_PartGraphVKway(&AA, theNodesIndex, theEdgesIndex, 
		 theNodesWeight, theEdgesWeight, 
		 &theWeightsFlag, &startFrom, &theParts, opts, &edgecut,
		 theNodesColour);
  */
  /*
	METIS_PartGraphKway(&AA, theNodesIndex, theEdgesIndex, 
	theNodesWeight, theEdgesWeight, 
	&theWeightsFlag, &startFrom, &theParts, opts, &edgecut,
	theNodesColour);
  */
							  
  /*  
	  METIS_WPartGraphRecursive(&AA, theNodesIndex, theEdgesIndex,
	  theNodesWeight, theEdgesWeight,
	  &theWeightsFlag, &startFrom,
	  &theParts, partWeights,
	  opts, &edgecut,
	  theNodesColour);
  */
  // end of other partition options		

  if (edgecut < 1)
	cout << endl << "Note: No edge was cut during the graph partition." << endl; cout.flush();
  //	ErrorExit((char *)"Metis  to partition.");
#ifdef __VERBAL2
  cout << "done."; cout.flush();
#endif

  //--------------------------------------
  // connectivity check
  if ((connMetisColour = (int *)malloc((metisActiveNodes+1) * sizeof(*connMetisColour))) == NULL)
	InsufficientMem((char *)"MAD Colour Allocating");
  if ((nodeVisited = (uchar *)malloc((metisActiveNodes+1) * sizeof(uchar))) == NULL)
	InsufficientMem((char *)"MAD Colour Allocating");
  // while we have more than two components, keep glueing
  while(checkConnectivity() > 2){}  
  // free connectivity lists
  free(connMetisColour);
  free(nodeVisited);
  // connectivity done
  //--------------------------------------

  // free the useless arrays
  free(theNodesWeight);
  theNodesWeight = NULL;
  free(theEdgesWeight);
  theEdgesWeight = NULL;
  free(theNodesIndex);
  theNodesIndex = NULL;

  //--------------------------------------
  // the MEtis graph is partitioned 
  //--------------------------------------
  return true;
}//partitionGraph()



//---------------------------------------------------------------------------
// Check the components connectivity for the partitioned  MEtis graph 
// If we have > 2 components glue one of them and return the remaining 
// number of components
// By invoking repeatedly we will get two connencted components.
//---------------------------------------------------------------------------
int
madd::checkConnectivity()
{
#ifdef __VERBAL2
  cout << endl << " Check connectivity..."; cout.flush();
#endif
#define markColor 1000   // a no color indicator

  // initialize connectivity
  LInt node1;
  for (int i=0; i < metisActiveNodes; i++){
	connMetisColour[i] = markColor;
	nodeVisited[i] = 0;
  }

#ifdef __VERBAL2
  cout << endl << " initial connectivity coloring done.  Color the connected components..."; cout.flush();
#endif

  //--------------------------------------
  // color the components, 
  // each with a different color
  //--------------------------------------
  int currColor = 0;
  bool graphIsColoured = false;
  while (!graphIsColoured){
	// there are still nodes to be colored
	LInt startNode;
	graphIsColoured = true;
	// find the first not colored node
	for (startNode = 0; startNode < metisActiveNodes && connMetisColour[startNode] != markColor; startNode++);
	if (startNode >= metisActiveNodes) 
	  continue; // we've reached the end

	//--------------------------------------
	// we have a new component
	connMetisColour[startNode] = currColor;
	int curNodesColour = theNodesColour[startNode];
	bool changed = true;	
	while (changed){
	  // there are stils parts of this component
	  changed = false;
	  // go through the all the nodes 
	  for (int i=startNode; i < metisActiveNodes; i++){
		if (nodeVisited[i] == 0){
		  graphIsColoured = false;  
		  if (connMetisColour[i] != markColor){
			// this is marked not visited node, and is colored.
			// mark the linked nodes
			nodeVisited[i] = 1; 
			int endLink =  theNodesIndex[i+1];
			// go through all the links of i node
			for (int l = theNodesIndex[i]; l < endLink; l++){
			  node1 = theEdgesIndex[l];
			  if (theNodesColour[node1] == curNodesColour && connMetisColour[node1] == markColor){
				// we have a new link to color
				connMetisColour[node1] = currColor;
				changed = true;
			  }
			}// for all links 
		  }// if!= markColor
		}//if !visited
	  }// for all nodes
	}// while changed
	// end of new component
	//--------------------------------------
   
	currColor++; // get next color, for next component
  }// while !graphIsColoured
  //--------------------------------------
  // the components are colored  
  // each with a different color
  //--------------------------------------

#ifdef __VERBAL2
  cout << "done." 
	   << endl << " check connectivity areas..."; cout.flush();
#endif
  //--------------------------------------
  // find how many components we have 
  // and the area of the components
  int endColor = currColor; // the number of components
  if (endColor <= 2){
	// there are only two pieces
#ifdef __VERBAL2
	cout << endl << " we have only two pieces"; cout.flush();
#endif
	return endColor;
  }
#ifdef __VERBAL2
  cout << endl << " we have a connectivity problem..."; cout.flush();
#endif

  //--------------------------------------
  // we have more than two components
  // create the component areas
  // DFloat connArea[currColor+1];
  // int colorOfArea[currColor+1];
  DFloat *connArea =  (DFloat*)malloc((endColor+1)*sizeof(DFloat) + 1);
  int *colorOfArea = (int*)malloc((endColor+1)*sizeof(int) + 1);
  if ( colorOfArea == NULL || connArea == NULL)
	InsufficientMem((char *)"connArea");

  for (int i=0; i < endColor; i++) {
	connArea[i]=0;
  }
  for (int i=0; i < endColor; i++) {
	connArea[i]=0;
  }
  // compute the size of the components
  // and find the original color of the pieces
  for (LLint metisAA = 0; metisAA < metisActiveNodes; metisAA++){
	connArea[connMetisColour[metisAA]] += theNodesWeight[metisAA];
	colorOfArea[connMetisColour[metisAA]] = theNodesColour[metisAA];
  }
  // find the size for the two subdomains
  // and the number of pieces for each
  DFloat areaPerColor[3] = {0, 0, 0};
  int piecesPerColor[3] = {0, 0, 0};
  for (int i=0; i < endColor; i++) {
	areaPerColor[colorOfArea[i]] += connArea[i];
	piecesPerColor[colorOfArea[i]]++;
  }
  //--------------------------------------
  // glue only one component	
  // get which subdomain to work on
  // get the bigger, unless it is one piece
  int maxColor = (areaPerColor[0] < areaPerColor[1])? 1: 0;
  if (piecesPerColor[maxColor] <= 1){
	maxColor = (maxColor == 1)? 0: 1;
  }
  // maxColor is our subdomain
  // get thme smallest component of this subdomain
  int color1, color2, color3;
  DFloat minArea = -1 ;
  for (int i=0; i < endColor; i++){
#ifdef __VERBAL2
	cout << endl << "Part " << i << "   area " << connArea[i] << "   color " <<  colorOfArea[i]; cout.flush();
#endif
	if ( colorOfArea[i] == maxColor && ((connArea[i] < minArea) || (minArea < 0))){
	  minArea =  connArea[i];
	  color3 = i;
	}// if
  }// for
#ifdef __VERBAL2
  cout << endl << " min Part " << color3 << "  min area " << minArea; cout.flush();
#endif
  // the smallest component is in color3
  // color3 is the component we will change
  // find the first color3 node
  for (node1=0; node1 < metisActiveNodes && connMetisColour[node1] != color3; node1++);
  // get the change color
  color2 = theNodesColour[node1];
  color1 =  (color2 == 0) ? 1 : 0;
#ifdef __VERBAL2
  cout << endl << " the metis min color: " << color2 << "  will change to " << color1; cout.flush();
#endif
  // change the color of this component
  while (node1 < metisActiveNodes){
	if (connMetisColour[node1] == color3)
	  theNodesColour[node1] = color1;
	node1++;
  }
  // component changed
  //--------------------------------------
  free(connArea);
  free(colorOfArea);

  return (endColor-1); // remaining components are endColor-1
}//checkConnectivity()


//-------------------------------------------------------------------------------
// inserts the separators in the geometry
// returns true if succesful
//-------------------------------------------------------------------------------
bool
madd::insertSeparators()
{
  //--------------------------------------
  // the following set of statements get the colors of the triangle
  // and checks the color matches
#define get_colors(thisTria)\
	p1 = tria_point1(thisTria);\
	p2 = tria_point2(thisTria);\
	p3 = tria_point3(thisTria);\
	node1 = node(thisTria, 0);\
	node2 = node(thisTria, 1);\
	node3 = node(thisTria, 2);\
	linkcolor1 = color1 =  metisColour(node1);\
	linkcolor2 = color2 =  metisColour(node2);\
	linkcolor3 = color3 =  metisColour(node3);\
	if ((linkTria1 = triaNeighborList[node1]) > thisTria){\
	  linknode1 = getLinkNode(linkTria1, thisTria);\
	  linkcolor1 =  metisColour(linknode1);\
	}\
	if ((linkTria2 = triaNeighborList[node2]) > thisTria){\
	  linknode2 = getLinkNode(linkTria2, thisTria);\
	  linkcolor2 =  metisColour(linknode2);\
	}\
	if ((linkTria3 = triaNeighborList[node3]) > thisTria){\
	  linknode3 = getLinkNode(linkTria3, thisTria);\
	  linkcolor3 =  metisColour(linknode3);\
	}\
	match1 = (color1 == linkcolor1);\
	match2 = (color2 == linkcolor2);\
	match3 = (color3 == linkcolor3);\
	sepIsColor1 = (color1 != color2) && (color1 != color3);\
	sepIsColor2 = (color2 != color1) && (color2 != color3);\
	sepIsColor3 = (color3 != color1) && (color3 != color2);\
	inSep = (sepIsColor1 || sepIsColor2 || sepIsColor3);
  // end of define
  //--------------------------------------

  // variables
  LInt p1, p2, p3;                  // the three vertices of a triangle
  DFloat *point1, *point2, *point3; // the three vertices of a triangle
  DFloat center[2];                 // the center of a triangle
  LLint node1, node2, node3;        // the thre nodes
  LLint linknode1, linknode2, linknode3; // the three adjacent nodes
  LLint linkTria, linkTria1, linkTria2, linkTria3; // the three adjacent triangles
  LLint linkNode; 
  int color1, color2, color3, color;  // store the colors
  int linkcolor1, linkcolor2, linkcolor3; // store the adjacent colors
  bool match1,  match2, match3;     // if the colors match with adjacent triangles 
  bool sepIsColor1, sepIsColor2, sepIsColor3; // locate the inner separator
  bool inSep;                       // if inner separator

  // the triangles that include separators
  LInt seperatorTriangles[400]; 
  noOfSeparatorTriangles = 0;   
#ifdef __VERBAL2
  cout << endl << "Reading the graph partition..."; cout.flush();
#endif
  //--------------------------------------
  // initialize 
  madSepStart = noOfSegms;
  LInt madPointStart = noOfPoints;
  for (int i=0; i < 16; i++) {
	partArea[i]=0;
	// partWeight[i]=0;
	edgesPerColour[i] = 0;
  }
  for (LInt i = 0; i < workNoOfSegms; i++){
	segm_colour(i) = -1;
  }
  // end initialize 
  //--------------------------------------

#ifdef __VERBAL2
  cout << endl << "color all edges first and compute part areas..."; cout.flush();
#endif

  //--------------------------------------
  //  locate the separator triangles
  //  check quality and return false if not acceptable
  //  color the edges and compute the subdomain areas
  //--------------------------------------
  bool failed = false; // for check if this is a valid partition
  for (LInt thisTria=0; thisTria < noOfTria; thisTria++) {
	// for all the triangle
	get_colors(thisTria); // get the colors
	// colour the work edges
	if (triaNeighborList[node1] < 0){
	  segm_colour(findWorkSegm(p2, p3)) = color1;
	  edgesPerColour[color1]++;
	}
	if (triaNeighborList[node2] < 0){
	  segm_colour(findWorkSegm(p1, p3)) = color2;
	  edgesPerColour[color2]++;
	}
	if (triaNeighborList[node3] < 0){
	  segm_colour(findWorkSegm(p1, p2)) = color3;
	  edgesPerColour[color3]++;
	}
	// add the part areas
	partArea[color1] += nodeArea[node1];
	partArea[color2] += nodeArea[node2];
	partArea[color3] += nodeArea[node3];
#ifdef __VERBAL3
	cout << "\n color edges for triangle " << thisTria << ": " 
		 << p1 << ", " << p2 << ", " << p3 << ". " ; cout.flush();
	cout << " colors: " 
		 << color1 << ", " << color2 << ", " << color3 << ". " ; cout.flush();
	cout << "   link colors: " 
		 << linkTria1 << ":" << linkcolor1 << ", " << linkTria2 << ":" << linkcolor2 << ", " 
		 << linkTria3 << ":" << linkcolor3 << ". " ; cout.flush();
#endif

	//--------------------------------------
    // check if we have separators in this triangle
	if (inSep || !match1 || !match2 || !match3){
	  // we have a separator
	  // add triangle to the list
	  seperatorTriangles[noOfSeparatorTriangles] = thisTria;
	  // add the separator points
	  if (sepIsColor1 || !match1){
		addSeparatorPoint1(noOfSeparatorTriangles, p2);
		addSeparatorPoint2(noOfSeparatorTriangles, p3);
	  }
	  if (sepIsColor2 || !match2){
		addSeparatorPoint1(noOfSeparatorTriangles, p1);
		addSeparatorPoint2(noOfSeparatorTriangles, p3);
	  }
	  if (sepIsColor3 || !match3){
		addSeparatorPoint1(noOfSeparatorTriangles, p1);
		addSeparatorPoint2(noOfSeparatorTriangles, p2);
	  }
	  noOfSeparatorTriangles++;
	  // separator trinagle is added
	  //--------------------------------------
#ifdef __DEBUG_MODE
	  if ((color1 != color2 || color1 != color3 || color2 != color3)
		  && !inSep){
		cout << "\n We have a HALF inner seprator!\n";
	  }   
#endif
	  //--------------------------------------
	  // check if two separators meet at the same point
	  //--------------------------------------
	  if (!match1 && !match2){
#ifdef __VERBAL3
		cout << "\n Collapse nodes: " << getEndNode(node1) << "  " <<  getEndNode(linknode1)
			 << "  that is points: " << p2 << ", " << p3; cout.flush();
#endif
		collapseNodeTo(linknode1, node1);   
		//	collapseNodeTo(linknode2, node2);   
		failed = true;;
	  }
	  if (!match1 && !match3){
#ifdef __VERBAL3
		cout << "\n Collapse nodes: " << getEndNode(node3) << "  " <<  getEndNode(linknode3)
			 << "  that is points: " << p1 << ", " << p2; cout.flush();
#endif
		//	collapseNodeTo(linknode1, node1);   
		collapseNodeTo(linknode3, node3);   
		failed = true;;
	  }
	  if (!match2 && !match3){
#ifdef __VERBAL3
	    cout << "\n Collapse nodes: " << getEndNode(node2) << "  " <<  getEndNode(linknode2)
			 << "  that is points: " << p1 << ", " << p3; cout.flush();
#endif
		collapseNodeTo(node2, linknode2);   
		//	collapseNodeTo(node3, linknode3);   
		failed = true;
	  }
	  if (!match1 && inSep){
#ifdef __VERBAL3
		cout << "\n Collapse nodes: " << getEndNode(node1) << "  " <<  getEndNode(linknode1); cout.flush();
#endif
		collapseNodeTo(linknode1, node1);   
		failed = true;;
	  }
	  if (!match2 && inSep){
#ifdef __VERBAL3
		cout << "\n Collapse nodes: " << getEndNode(node2) << "  " <<  getEndNode(linknode2); cout.flush();
#endif
		collapseNodeTo(linknode2, node2);   
		failed = true;;
	  }
	  if (!match3 && inSep){
#ifdef __VERBAL3
	    cout << "\n Collapse nodes: " << getEndNode(node3) << "  " <<  getEndNode(linknode3); cout.flush();
#endif
		collapseNodeTo(node3, linknode3);   
		failed = true;
	  }
	  //--------------------------------------
	  // end of check if two separators meet 
	  //--------------------------------------

	} // if this was a separator tringle

  }// for all triangles
  //--------------------------------------
  //  end of locating the separator triangles
  //--------------------------------------


  if (failed) // we have two seprators meeting at the same point
	return false;

  //--------------------------------------
  // still there are cases where seprators meet at the same point
  // check and correct
  //--------------------------------------
  for (int sepTria = 0; sepTria < noOfSeparatorTriangles; sepTria++) {
	// go through all the separator triangles
	// get the two boundary points of the separator
	p1 = seperatorPoint1[sepTria];
	p2 = seperatorPoint2[sepTria];
	// check if any of the next triangles have the same separator points
	for (int i = sepTria+1; i < noOfSeparatorTriangles; i++){
	  if (p1 == seperatorPoint1[i] || p2 == seperatorPoint1[i]
		  || p1 == seperatorPoint2[i] || p2 == seperatorPoint2[i]){
		// we have a dupolicate separator point
		// check if it is on inner seprator
		LInt thisTria = seperatorTriangles[i];
		get_colors(thisTria);
		if (inSep){
		  // we will get rid of this inner seperator
		  if (sepIsColor1)
			collapseNodeTo(node1, node2);   
		  if (sepIsColor2)
			collapseNodeTo(node2, node3);   
		  if (sepIsColor3)
			collapseNodeTo(node3, node1);   
		  return false; 
		}

		// check if in sepTria is inner separator
		thisTria = seperatorTriangles[sepTria];
		get_colors(thisTria);
		if (inSep){
		  // we will get rid of this inner seperator
		  if (sepIsColor1)
			collapseNodeTo(node1, node2);   
		  if (sepIsColor2)
			collapseNodeTo(node2, node3);   
		  if (sepIsColor3)
			collapseNodeTo(node3, node1);   
		  return false; 
		}

		// here we have two outer non-adjucent separators at the same point
		// destroy one of them
		if (!match1 ){
		  collapseNodeTo(linknode1, node1);   
		  return false;
		}
		if (!match2 ){
		  collapseNodeTo(linknode2, node2);   
		  return false;
		}
		if (!match3 ){
		  collapseNodeTo(linknode3, node3);   
		  return false;
		}
		
 
	  }//if we have a dupolicate separator point

	}// for the next triangles
  
  }// for all the seprator triangles
  //--------------------------------------
  // end of checking if two separators meet at the same point
  //--------------------------------------
	  
  // if we don't have serators something is wrong 
  if (noOfSeparatorTriangles < 1){
	ErrorExit((char *)"No Separators could be found.");
  }

#ifdef __VERBAL2
  cout << endl << "inserting the separators... No Of triangles:" << noOfSeparatorTriangles; cout.flush();
#endif
  //--------------------------------------
  //  we have the separator triangles
  //  the quality is accepted
  //  the domain edges are colored
  //  insert separators
  //--------------------------------------
  for (int sepTria = 0; sepTria < noOfSeparatorTriangles; sepTria++) {
	// go through all the  separator triangles
	LInt thisTria = seperatorTriangles[sepTria];
	sepTriangleIndex = sepTria; //let the rest of the methods know which is the current triangle
	get_colors(thisTria);       // get the colors
#ifdef __VERBAL3
	cout << endl << " checking separators for triangle " << thisTria << ": "
		 << p1 << ", " << p2 << ", " << p3 << ". " ; cout.flush();
	cout << " colors:" 
		 << color1 << ", " << color2 << ", " << color3 << ". " ; cout.flush();
	cout << endl << "check if there is an inner separator..."; cout.flush();
#endif
	// check if there is an inner separator
	if (inSep){
#ifdef __VERBAL3
	  cout << endl << "we have an innner separator..." 
		   << p1 << ", " << p2 << ", " << p3 << ". " ; cout.flush();
#endif
	  //--------------------------------------
	  // we have an innner separator
	  point1 = point(p1);
	  point2 = point(p2);
	  point3 = point(p3);
	  circumCenter(point1, point2, point3, center);
	  // find the two partial separators
	  LInt inp1, inp2;
	  int sepcolor1;
	  if (sepIsColor1){
		inp1 = p2; inp2 = p3; sepcolor1 = color1;
	  }
	  if (sepIsColor2){
		inp1 = p1; inp2 = p3; sepcolor1 = color2;
	  }
	  if (sepIsColor3){
		inp1 = p1; inp2 = p2; sepcolor1 = color3;
	  }
	  // insert it
	  insertInnerSeparator(inp1, inp2, center, sepcolor1);
	} // inner separator done
	//--------------------------------------

#ifdef __VERBAL3
	cout << endl << "check if there is an outer separator..."; cout.flush();
#endif
	//--------------------------------------
	// check if there is an outer separator
	if (!match1){
	  // we have an outer separator
	  insertEdgeSep(p2, p3);
	}
	if (!match2){
	  // we have an outer separator
	  insertEdgeSep(p1, p3);
	}
	if (!match3){
	  // we have an outer separator
	  insertEdgeSep(p1, p2);
	}
	// outer separator done
	//--------------------------------------

  }// for all separator triangles
  //--------------------------------------
  //  separators are inserted 
  //--------------------------------------

  return true;
}//insertSeparators() 

 


//-------------------------------------------------------------------------------
// insert an inner separator to the geometry
//-------------------------------------------------------------------------------
void
madd::insertInnerSeparator(LInt inp1,        // the first boundary point
						   LInt inp2,        // the second boundary point
						   DFloat *center,   // the center
						   int sepcolor1)    // the color of the included node (1 - 2 nodes is the partition)
{
#ifdef __VERBAL3
  cout << endl << "Checking inner separator  between points " << inp1 << ", " << inp2;
  cout.flush();
#endif

  //DFloat outnewArea[4], innewArea[4]; // keep track of the subdomain areas
  DFloat outnewArea[2], innewArea[2]; // keep track of the subdomain areas
  // the points for smoothing (check for inner and outer separators)
  LInt outop1 =  inp1;               
  LInt outop2 =  inp2;
  LInt inop1 =  inp1;
  LInt inop2 =  inp2;
  // the smoothing center
  DFloat inopCenter[2];
  x(inopCenter) = x(center); 
  y(inopCenter) = y(center);
  // this will hold the boundary edges that will require coclor switching
  LInt outchangeEdgesColour[25], inchangeEdgesColour[25];

  // get the opposite color
  int sepcolor2 = (sepcolor1 == 0) ? 1: 0;
  // get the area of the included node 
  DFloat area =  areaOfTriangle(point(inp1), center, point(inp2));
  // the initial inner areas 
  innewArea[0] = partArea[0];
  innewArea[1] = partArea[1];
  // the initial outer areas 
  outnewArea[sepcolor2] = partArea[sepcolor2] + area;
  outnewArea[sepcolor1] = partArea[sepcolor1] - area;

  //--------------------------------------
  //  start the smoothing procedure
  //  fint the best outer and inner separator
  DFloat outcoef = findOptSep(&outop1, &outop2, outnewArea, outchangeEdgesColour);
  DFloat incoef = findOptSep(&inop1, &inop2, inopCenter, innewArea, inchangeEdgesColour);
#ifdef __VERBAL3
  cout << endl << "inner separator between points " << inop1 << ", " << inop2
	   << " has smooth coeff:" << incoef;
  cout << endl << "outer separator between points " << outop1 << ", " << outop2 
	   << " has smooth coeff:" << outcoef; cout.flush();
#endif
  // if the new areas are not balanced reject them
  if (!isBalanced(outnewArea))
	outcoef = MAX_COEF; // the outer separator is not accpted due to imbalance diffenece 
  if (!isBalanced(innewArea))
	incoef = MAX_COEF; // the outer separator is not accpted due to imbalance diffenece 

  
  if ((incoef >= MAX_COEF && outcoef >= MAX_COEF) // if we have bad separators
	  ||  MAX_SMOOTH_DEPTH <= 0){                 // or no smoothing flag
#ifdef __VERBAL3
	cout << endl << "we have a bad separator or MAX_SMOOTH_DEPTH <= 0.  Inserting  inner separator..." ; cout.flush();
#endif
	//--------------------------------------
	// insert the initial inner separator 
	// compute the weight of the center and insert it to the geometry
	int centerWeight = (int)((pointWeight[inop1] + pointWeight[inop2]) * INTERPOLATION_WEIGHT_COEF / 2.0); 
	LInt newp = NewPoint(x(center), y(center), centerWeight); 
	// insert the two partial separators
	NewSegm(newp, inp1);
	NewSegm(newp, inp2);
	// the two boundary points have been used 
	addSeparatorPoint1(sepTriangleIndex, inp1);
	addSeparatorPoint2(sepTriangleIndex, inp2);
#ifdef __VERBAL3
	cout << "done." ; cout.flush();
#endif
	// initial inner separator inserted 
	//--------------------------------------
	return;
  }
  
  //--------------------------------------
  // insert the better inner or outer separator 
  if (outcoef > incoef){
	//--------------------------------------
	// the new inner separator is better
	// insert it
#ifdef __VERBAL3
	cout << endl << " proceed with the inner separator..." ; cout.flush();
#endif
	// compute the weight of the center and insert it to the geometry
	int centerWeight = (int)((pointWeight[inop1] + pointWeight[inop2]) * INTERPOLATION_WEIGHT_COEF / 2.0); 
	LInt newp = NewPoint(x(inopCenter), y(inopCenter), centerWeight); 
	// insert the two partial separators
	NewSegm(newp, inop1);
	NewSegm(newp, inop2);
	// the two boundary points have been used 
	addSeparatorPoint1(sepTriangleIndex, inop1);
	addSeparatorPoint2(sepTriangleIndex, inop2);
	// assign the new areas to the subdomains
	LInt segm1;
	int color1, color2;
	partArea[0] =  innewArea[0];
	partArea[1] =  innewArea[1];
	LInt endEdge = inchangeEdgesColour[0];
	for (int i = 1; i < endEdge; i++){
	  segm1 = inchangeEdgesColour[i];
	  color1 = segm_colour(segm1);
	  color2 = (color1 == 0) ? 1: 0;
	  segm_colour(segm1) = color2;
	  edgesPerColour[color2]++;
	  edgesPerColour[color1]--;
	}
	// the new inner separator is inserted
	//--------------------------------------
#ifdef __VERBAL3
	cout << endl << "  inner separator inserted" ; cout.flush();
#endif
	// smooth it further if needed
	//	insertInnerSeparator(inop1, inop2, inopCenter, int sepcolor1, int depth);	
	return;
  }

#ifdef __VERBAL3
  cout << endl << " proceed with the the outer separator ... "; cout.flush(); 
#endif
  //--------------------------------------
  // the outer separator is better
  //  assign the new areas and edge colours
  LInt segm1;
  int color1, color2;
  partArea[0] =  outnewArea[0];
  partArea[1] =  outnewArea[1];
  LInt endEdge = outchangeEdgesColour[0];
  for (int i = 1; i < endEdge; i++){
	segm1 = outchangeEdgesColour[i];
	color1 = segm_colour(segm1);
	color2 = (color1 == 0) ? 1: 0;
	segm_colour(segm1) = color2;
	edgesPerColour[color2]++;
	edgesPerColour[color1]--;
  }
  // smooth it further if needed
  insertOptEdgeSep(outop1, outop2, 1);
  // the outer separator is inserted
  //--------------------------------------

#ifdef __VERBAL3
  cout << "done." ; cout.flush();
#endif
  return;
}//insertInnerSeparator()


//-------------------------------------------------------------------------------
// insert the outer separator to the geometry
//-------------------------------------------------------------------------------
void
madd::insertEdgeSep(LInt p1,         // the first boundary point 
					LInt p2)         // the second boundary point
{
#ifdef __VERBAL3
  cout << endl << endl << " We have an outer separator... " << p1 <<"(" << (int)basicPoint[p1] << "), " << p2
	   << "(" << (int)basicPoint[p2] << "),"; cout.flush(); 
#endif
  // if we have smoothing
  if (MAX_SMOOTH_DEPTH <= 0){
	// find the best separator
	// and insert it
	insertOptEdgeSep(p1, p2, 0);
	return;
  }
  // else just insert and return
  insertOptEdgeSep(p1, p2, 0);
  return;
}//insertEdgeSep()


//-------------------------------------------------------------------------------
// insert the optimal outer separator to the geometry
//-------------------------------------------------------------------------------
void
madd::insertOptEdgeSep(LInt p1,         // the first boundary point
					   LInt p2,         // the second boundary point
					   int depth)       // the current smoothing depth
{
#ifdef __VERBAL3
  cout << endl << endl << " Smoothing the outer separator... " << p1 <<"(" << (int)basicPoint[p1] << "), " << p2
	   << "(" << (int)basicPoint[p2] << "),"  
	   << "  depth:" << depth << endl ; cout.flush();
#endif

  // check if we have reached the max smoothing depth
  if (depth >= MAX_SMOOTH_DEPTH){
#ifdef __VERBAL3
	cout << "MAX Depth reached. Inserting the separator and stop smoothing .";
#endif
	//--------------------------------------
	// just insert the current seperator
	NewSegm(p1, p2);
	addSeparatorPoint1(sepTriangleIndex, p1);
	addSeparatorPoint2(sepTriangleIndex, p2);
	return;
	//--------------------------------------
  }

  // continue the smoothing
  // the new areas
  DFloat area, newArea[2]; 
  newArea[0] = partArea[0];
  newArea[1] = partArea[1];
  // the two optimal boundary points
  LInt op1 =  p1;
  LInt op2 =  p2;
  // to keep track of the edges that will switch color
  LInt changeEdgesColour[25];

  //--------------------------------------
  // find the optimal outer seperator
  DFloat refCoef = findOptSep(&op1, &op2, newArea, changeEdgesColour);
  if ((op1 == p1 && op2 == p2) || (refCoef >= MAX_COEF)){
	//--------------------------------------
	// the current separator is better
	// no further smoothing
	// insert separator and return
#ifdef __VERBAL3
	cout << "Local Optimal reached. Inserting the separator and stop smoothing.";
#endif
	NewSegm(p1, p2);
	addSeparatorPoint1(sepTriangleIndex, p1);
	addSeparatorPoint2(sepTriangleIndex, p2);
	return;
	//--------------------------------------
  }  
		
  //--------------------------------------
  // check the new imbalance 
  if (!isBalanced(newArea)){
#ifdef __VERBAL3
	cout << endl << "The new separator is imbalanced"; cout.flush();
#endif
	//--------------------------------------
	// The new separator is imbalanced
	// just insert the current separator
	NewSegm(p1, p2);
	addSeparatorPoint1(sepTriangleIndex, p1);
	addSeparatorPoint2(sepTriangleIndex, p2);
	return;
	//--------------------------------------
  }

#ifdef __VERBAL3
  cout << endl << endl << "A better separator was found: " << op1 <<"(" << (int)basicPoint[op1] << "), " << op2
	   << "(" << (int)basicPoint[op2] << "),"  
	   << "  with smooth coeff:" << refCoef << endl; cout.flush();
#endif
  //--------------------------------------
  //  the new balance is accepted
  //  assign new areas and edge colours
  LInt segm1;
  int color1, color2;
  partArea[0] =  newArea[0];
  partArea[1] =  newArea[1];
  LInt endEdge = changeEdgesColour[0];
  for (int i = 1; i < endEdge; i++){
	segm1 = changeEdgesColour[i];
	color1 = segm_colour(segm1);
	color2 = (color1 == 0) ? 1: 0;
	segm_colour(segm1) = color2;
	edgesPerColour[color2]++;
	edgesPerColour[color1]--;
  }

  // smooth it further if needed
  insertOptEdgeSep(op1, op2, depth+1);
  //--------------------------------------
}//insertEdgeSep()



//-------------------------------------------------------------------------------
// find the optimal inner separator 
// returns the best coefficient 
// and the new separator in the parameters
//-------------------------------------------------------------------------------
DFloat
madd::findOptSep(LInt *op1,                // the first boundary point, initial, optimal
				 LInt *op2,                // the second boundary point, initial, optimal
				 DFloat *center,           // the center, initial, optimal
				 DFloat *newArea,          // the new areas of the subdomains
				 LInt *changeEdgesColour)  // the list of the edges to change color
{

  LInt segm1, segm2; // just two segments
  // get the smoothing interval
  int middle =  SMOOTHING_POINTS / 2; // the middle point index
  int end =  SMOOTHING_POINTS;        // the end  point index
  // the smoothing adjacent points
  LInt nextp1[34];
  LInt nextp2[34];
  // no edges will change color so far
  changeEdgesColour[0] = 0;
  // start with the initial separator
  LInt p1 = *op1;   
  LInt p2 = *op2;
#ifdef __VERBAL3
  cout << endl << endl << " Finding inner optimal separator for " << p1 <<"(" << (int)basicPoint[p1] << "), " << p2
	   << "(" << (int)basicPoint[p2] << ")" << endl ; cout.flush();
#endif
  // get the adjacent points 
  nextPointsOfPoint(p1, nextp1, end);
  nextPointsOfPoint(p2, nextp2, end);

#ifdef __VERBAL3
  cout << endl << "Computing the coefficients... "; cout.flush();
#endif
  //--------------------------------------
  // evaluate the adjacent points' separators
  //--------------------------------------
  DFloat coef = 0;
  // the initial points are our starting optimal
  LInt opt1 = middle;
  LInt opt2 = middle;
  DFloat optCoef =  smoothingCoef1(nextp1[middle], center, nextp2[middle]);
  //--------------------------------------
  // check the right-right points
  // if coef >  MAX_COEF we have an intersection or used vertex, stop
  for (int i = middle; i < end && coef <= MAX_COEF; i++){
	for (int j = middle+1; j < end && coef <= MAX_COEF; j++){
	  coef =  smoothingCoef1(nextp1[i], center, nextp2[j]);
	  if ( coef < optCoef){
		optCoef = coef;
		opt1 = i;
		opt2 = j;
	  }// if
	}// for
  }//for
  //
  //--------------------------------------
  // check the right-left points
  // if coef >  MAX_COEF we have an intersection or used vertex, stop
  for (int i = middle; i < end && coef <= MAX_COEF; i++){
	for (int j = middle-1; j >= 0 && coef <= MAX_COEF; j--){
	  coef =  smoothingCoef1(nextp1[i], center, nextp2[j]);
	  if ( coef < optCoef){
		optCoef = coef;
		opt1 = i;
		opt2 = j;
	  }// if
	}//for
  }//for
  //--------------------------------------
  // check the left-right points
  // if coef >  MAX_COEF we have an intersection or used vertex, stop
  for (int i = middle-1; i >= 0 && coef <= MAX_COEF; i--){
	for (int j = middle+1; j < end && coef <= MAX_COEF; j++){
	  coef =  smoothingCoef1(nextp1[i], center, nextp2[j]);
	  if ( coef < optCoef){
		optCoef = coef;
		opt1 = i;
		opt2 = j;
	  }// if
	}//for
  }//for
  //--------------------------------------
  // check the left-left points
  // if coef >  MAX_COEF we have an intersection or used vertex, stop
  for (int i = middle-1; i >= 0 && coef <= MAX_COEF; i--){
	for (int j = middle-1; j >= 0 && coef <= MAX_COEF; j--){
	  coef =  smoothingCoef1(nextp1[i], center, nextp2[j]);
	  if ( coef < optCoef){
		optCoef = coef;
		opt1 = i;
		opt2 = j;
	  }// if
	}//for
  }//for
#ifdef __VERBAL3
  cout << "done "; cout.flush();
#endif  
  //--------------------------------------
  // adjacent points have been evaluated
  //--------------------------------------
  // get the best points
  *op1 = nextp1[opt1];
  *op2 = nextp2[opt2];
  if (*op1 == p1 && *op2 == p2){
	// the initial points are the best
	return optCoef;
  }  

  //--------------------------------------
  // we have a bettter basic candidate
#ifdef __VERBAL3
  cout << endl << "We found a better inner separator during smoothing :"
	   << nextp1[opt1] << "(" << (int)basicPoint[nextp1[opt1]] << "),"  
	   << nextp2[opt2] << "(" << (int)basicPoint[nextp2[opt2]] << "),"; cout.flush();
#endif  
  // assign the new areas 
  // and the edges that will need color change 
  DFloat area; // the new areas
  uchar color1, color2; // two colors
  LInt np1, np2; // two next points
  int addto;

  // get the changes for point opt1
  if (opt1 < middle)
	addto = -1; // otp1 is on the left
  else
	addto = 1;  // otp1 is on the right
  LInt edgeAA = 1;
  for (int i = middle; i != opt1; i += addto){
	np1 = nextp1[i];
	np2 = nextp1[i+addto];
	segm1 = findWorkSegm(np1, np2);
	changeEdgesColour[edgeAA++] = segm1; // this segment will change color
	// compute the new area
	color1 = segm_colour(segm1);
	color2 = (color1 == 0) ? 1: 0;
	//    cout << endl << i << ". Computing area 1..."; cout.flush();
	area =  areaOfTriangle(center, point(np1), point(np2));
	newArea[color2] += area;
	newArea[color1] -= area;
  }//for

  // get the changes for point opt2
  if (opt2 < middle)
	addto = -1; // otp2 is on the left
  else
	addto = 1; // otp2 is on the right
  for (int i = middle; i != opt2; i += addto){
	np1 = nextp2[i];
	np2 = nextp2[i+addto];
	segm1 = findWorkSegm(np1, np2);
	changeEdgesColour[edgeAA++] = segm1; // this segment will change color
	// compute the new area
	color1 = segm_colour(segm1);
	color2 = (color1 == 0) ? 1: 0;
	//    cout << endl << i << ". Computing area 2..."; cout.flush();
	area =  areaOfTriangle(center, point(np1), point(np2));
	newArea[color2] += area;
	newArea[color1] -= area;
  }//for

  // set the number of edegs to change color
  changeEdgesColour[0] = edgeAA;
  // bettter inner separator is done
  //--------------------------------------
  return optCoef;
}//findOptSep()



//-------------------------------------------------------------------------------
// find the optimal outer separator 
// returns the best coefficient 
// and the new separator in the parameters
//-------------------------------------------------------------------------------
DFloat
madd::findOptSep(LInt *op1,                // the first boundary point, initial, optimal
				 LInt *op2,                // the second boundary point, initial, optimal
				 DFloat *newArea,          // the new areas of the subdomains
				 LInt *changeEdgesColour)  // the list of the edges to change color
{
#ifdef __VERBAL3
  cout << endl << endl << " Finding outer optimal separator for " << p1 <<"(" << (int)basicPoint[p1] << "), " << p2
	   << "(" << (int)basicPoint[p2] << ")" << endl ; cout.flush();
#endif
  LInt segm1, segm2;
  // get the smoothing interval
  int middle =  SMOOTHING_POINTS / 2; // the middle point index
  int end =  SMOOTHING_POINTS;        // the end  point index
  // the smoothing adjacent points
  LInt nextp1[34];
  LInt nextp2[34];
  // the initial points
  LInt p1 = *op1;
  LInt p2 = *op2;
  // no edges change color so far
  changeEdgesColour[0] = 0;

  // get the adjacent points
  nextPointsOfPoint(p1, nextp1, end);
  nextPointsOfPoint(p2, nextp2, end);

  //--------------------------------------
  // evaluate the adjacent point separators
  //--------------------------------------
  DFloat coef = 0;
  LInt opt1 = middle;
  LInt opt2 = middle;
  DFloat optCoef = smoothingCoef2(nextp1[middle], nextp2[middle]);

#ifdef __VERBAL3
  cout << endl << "Computing the coefficients... "; cout.flush();
#endif
  //--------------------------------------
  // check the right-right points
  // if coef >  MAX_COEF we have an intersection or used vertex, stop
  for (int i = middle; i < end && coef <= MAX_COEF; i++){
	for (int j = middle+1; j < end && coef <= MAX_COEF; j++){
	  coef =  smoothingCoef2(nextp1[i], nextp2[j]);
	  if ( coef < optCoef){
		optCoef = coef;
		opt1 = i;
		opt2 = j;
	  }// if
	}//for
  }//for
  //--------------------------------------
  // check the right-left points
  // if coef >  MAX_COEF we have an intersection or used vertex, stop
  for (int i = middle; i < end && coef <= MAX_COEF; i++){
	for (int j = middle-1; j >= 0 && coef <= MAX_COEF; j--){
	  coef =  smoothingCoef2(nextp1[i], nextp2[j]);
	  if ( coef < optCoef){
		optCoef = coef;
		opt1 = i;
		opt2 = j;
	  }// if
	}//for
  }//for
  //--------------------------------------
  // check the left-right points
  // if coef >  MAX_COEF we have an intersection or used vertex, stop
  for (int i = middle-1; i >= 0 && coef <= MAX_COEF; i--){
	for (int j = middle+1; j < end && coef <= MAX_COEF; j++){
	  coef =  smoothingCoef2(nextp1[i], nextp2[j]);
	  if ( coef < optCoef){
		optCoef = coef;
		opt1 = i;
		opt2 = j;
	  }// if
	}//for
  }//for
  //--------------------------------------
  // check the left-left points
  // if coef >  MAX_COEF we have an intersection or used vertex, stop
  for (int i = middle-1; i >= 0 && coef <= MAX_COEF; i--){
	for (int j = middle-1; j >= 0 && coef <= MAX_COEF; j--){
	  coef =  smoothingCoef2(nextp1[i], nextp2[j]);
	  if ( coef < optCoef){
		optCoef = coef;
		opt1 = i;
		opt2 = j;
	  }// if
	}//for
  }//for
#ifdef __VERBAL3
  cout << "done "; cout.flush();
#endif  
  //--------------------------------------
  // adjacent points have been evaluated
  //--------------------------------------
  // get the best points
  *op1 = nextp1[opt1];
  *op2 = nextp2[opt2];
  if (*op1 == p1 && *op2 == p2){
	// the initial points are the best
	// just return
	return optCoef;
  }  
#ifdef __VERBAL3
  cout << endl << "We found a better separator during smoothing :"
	   << nextp1[opt1] << "(" << (int)basicPoint[nextp1[opt1]] << "),"  
	   << nextp2[opt2] << "(" << (int)basicPoint[nextp2[opt2]] << "),"
	   << "  with smooth coeff:" << optCoef << endl; cout.flush();
#endif  
  //--------------------------------------
  // we have a bettter basic candidate
  // assign the new areas 
  // and the edges that will need color change 
  DFloat area; // the new areas
  uchar color1, color2; // two colors
  LInt bp1, bp2, np1, np2; // two base points, two next points
  int addto;

  // get the changes for point opt1
  bp1 = p1; bp2 = p2;
  if (opt1 < middle)
	addto = -1;  // opt1 is on the left
  else
	addto = 1;   // opt1 is on the right
   LInt edgeAA = 1;
  for (int i = middle; i != opt1; i += addto){
	np1 = nextp1[i];
	np2 = nextp1[i+addto];
	segm1 = findWorkSegm(np1, np2);
	// this segment will change color
	changeEdgesColour[edgeAA++] = segm1; 
	// compute the new area
	color1 = segm_colour(segm1);
	color2 = (color1 == 0) ? 1: 0;
	area =  areaOfTriangle(point(bp2), point(np1), point(np2));
	newArea[color2] += area;
	newArea[color1] -= area;
  }

  // get the changes for point opt2
  bp1 = nextp1[opt1];
  if (opt2 < middle)
	addto = -1;  // opt2 is on the left
  else
	addto = 1;   // opt2 is on the right
  for (int i = middle; i != opt2; i += addto){
	np1 = nextp2[i];
	np2 = nextp2[i+addto];
	segm1 = findWorkSegm(np1, np2);
	// this segment will change color
	changeEdgesColour[edgeAA++] = segm1;
	// compute the new area
	color1 = segm_colour(segm1);
	color2 = (color1 == 0) ? 1: 0;
	area =  areaOfTriangle(point(bp1), point(np1), point(np2));
	newArea[color2] += area;
	newArea[color1] -= area;
  }
  // assign the number of edges to change color
  changeEdgesColour[0] = edgeAA;
  // bettter outer separator is done
  //--------------------------------------
  return optCoef;
}//findOptSep()


//-------------------------------------------------------------------------------
// check if the new area is acceptably balanced
//-------------------------------------------------------------------------------
bool
madd::isBalanced(DFloat *newArea)
{
  DFloat totalArea =  newArea[1] + newArea[0];
  // the new balance
  DFloat balance = (newArea[0] < newArea[1])? 
	(newArea[1] / totalArea) : (newArea[0] / totalArea);
  // the original balance
  DFloat obalance = (partArea[0] < partArea[1])? 
	(partArea[1] / totalArea) : (partArea[0] / totalArea);

  // if the new balance is better is acceptable
  if (obalance >= balance)
	return true;

  // compute the area difference
  DFloat area = (newArea[0] > partArea[0])? (newArea[0] -  partArea[0]) : (partArea[0] - newArea[0]);
#ifdef __VERBAL3
  cout << endl << "Balance parameters:    balance: " << balance << ",  imbalance diffrenece: " 
	   << area / totalArea; cout.flush();
#endif
  // if the differnce is big and the new imbalance is > normal balance
  if (( (area / totalArea) >  MAX_IMBALANCE_DIFF && balance >  MAX_IMBALANCE_NORMAL)
  // or the new balance is > max limit
	  ||  (balance >  MAX_IMBALANCE_LIMIT)) {
	// it is not acceptable
	return false;
  }
  // otherwise is good
  return true;
}//isBalanced


//-------------------------------------------------------------------------------
// compute the basic coefficient ((used in the MADD graph generation)
// for an inner partial separator
//-------------------------------------------------------------------------------
DFloat
madd::basicCoef1(DFloat* center, // the circumcenter
				 LInt p1)        // the boundary point
{
  // get the two angles, and the minimum
  DFloat angle1 = angle(center, point(p1), point(leftOfPoint(p1)));
  DFloat angle2 = angle(center, point(p1), point(rightOfPoint(p1)));
  angle1 = min(angle1, angle2);
  // if the angle is not acceptable return the MAX_COEF
  if (angle1 < PHI) 
	return  MAX_COEF;
  //--------------------------------------
  // compute the coefficient
  // pi/2 is the max for our formula
  if (angle1 > PI_OVER_TWO)
	angle1 = PI_OVER_TWO;
  // get the type of the external point
  DFloat coeff = theBasicCoeff[basicPoint[p1]];
  // evaluate formula  
  coeff *= 1 / (angle1 - PHI + 1.0);
  // we will prefer a straight line, if possible
  coeff *= 1.15;
  // coefficient is computed
  //--------------------------------------
  return coeff;
}// basicCoef1


//-------------------------------------------------------------------------------
// compute the basic coefficient ((used in the MADD graph generation)
// for an outer separator
//-------------------------------------------------------------------------------
DFloat
madd::basicCoef2(LInt p1,  // the first  boundary point 
				 LInt p2)  // the second  boundary point 
{
  // get the minimum of the four external angles
  DFloat minangle = minExtAngle(p1, p2);
  // if the angle is not acceptable return the MAX_COEF
  if (minangle < PHI)
	return MAX_COEF;
  //--------------------------------------
  // compute the coefficient
  // pi/2 is the max for our formula
  if (minangle > PI_OVER_TWO)
	minangle = PI_OVER_TWO;
  // get the types of the external point
  DFloat coeff = theBasicCoeff[basicPoint[p1]];
  coeff *= theBasicCoeff[basicPoint[p2]];
  // evaluate formula  
  coeff *= 1 / (minangle - PHI + 1.0);
#ifdef  __VERBAL4
  cout << endl << " the Min Angle is : " << minangle * 180 / M_PI 
	   << " it contributes: " <<  1 / (minangle - PHI + 1.0)
	   << " to the coeff: " << coeff;
#endif
  // coefficient is computed
  //--------------------------------------
  return coeff;
}//basicCoef2()


//-------------------------------------------------------------------------------
// compute the smoothing coefficient ((used in the insertion of separators)
// for an inner separator
//-------------------------------------------------------------------------------
DFloat
madd::smoothingCoef1(LInt p1,          // the first  boundary point 
					 DFloat* center,   // the circumcenter
					 LInt p2)          // the second boundary point
{
  // if the points are in use by another separator reject them
  if (separatorPointUsed(p1) || separatorPointUsed(p2)){
#ifdef __VERBAL3
	cout <<"\n We have a duplicate separator: " << p1 << ", " << p2; cout.flush();
#endif
	return (MAX_COEF+1); //stop the search to this direction
  }
  // if it intersects the boundary, reject it
  if (intersectsBoundary(center, p1))
	return  (MAX_COEF+1);  //and stop the search to this direction
  if (intersectsBoundary(center, p2))
	return  (MAX_COEF+1);  //and stop the search to this direction
  // get the 5 angles
  DFloat angle1 = angle(center, point(p1), point(leftOfPoint(p1)));
  DFloat angle2 = angle(center, point(p1), point(rightOfPoint(p1)));
  DFloat angle3 = angle(center, point(p2), point(leftOfPoint(p2)));
  DFloat angle4 = angle(center, point(p2), point(rightOfPoint(p2)));
  DFloat angle5 = angle(point(p1), center, point(p2));
  // get the minimum angle
  angle1 = min(angle1, angle2);
  DFloat minangle = min(angle1, angle2);
  if (minangle > angle3)
	minangle = angle3; 
  if (minangle > angle4)
	minangle = angle4; 
  if (minangle > angle5)
	minangle = angle5; 
  // if the angle is not acceptable return the MAX_COEF
  if (minangle < PHI)
	return  MAX_COEF;
  //--------------------------------------
  // compute the coefficient
  // pi/2 is the max for our formula
  if (minangle > PI_OVER_TWO)
	minangle = PI_OVER_TWO;
  // get the types of the external point
  DFloat coeff = theSmoothCoeff[basicPoint[p1]];
  coeff *= theSmoothCoeff[basicPoint[p2]];
  // evaluate formula  
  coeff *= 1 / (minangle - PHI + 1.0);
  coeff *= (2 * distance(point(p1), center) / totalBoundaryLength);
  // we will prefer a straight line, if possible
  coeff *= 1.18;
  // coefficient is computed
  //--------------------------------------
  return coeff;
}//smoothingCoef1()


//-------------------------------------------------------------------------------
// compute the smoothing coefficient ((used in the insertion of separators)
// for an outer separator
//-------------------------------------------------------------------------------
DFloat
madd::smoothingCoef2(LInt p1, LInt p2)
{
#ifdef  __VERBAL4
  cout << endl << endl << " Smoothing coefficient of " << p1 << ", " << p2 
	   << ". The basic flags are " << (int)basicPoint[p1] 
	   << ",  " <<  (int)basicPoint[p2]; cout.flush();
#endif
  // if the points are in use by another separator reject them
  if (separatorPointUsed(p1) || separatorPointUsed(p2)){
#ifdef __VERBAL3
	cout <<"\n We have a duplicate separator: " << p1 << ", " << p2; cout.flush();
#endif
	return (MAX_COEF+1);  //stop the search to this direction
  }
  // if it intersects the boundary, reject it
  if (intersectsBoundary(p1, p2))
	return (MAX_COEF+1);  // and stop the search to this direction

  // get the minimum angle of the four angles
  DFloat minangle = minExtAngle(p1, p2);
#ifdef  __VERBAL4
  cout << endl << "min angle : " << minangle; cout.flush();
#endif
  // if the angle is not acceptable return the MAX_COEF
  if (minangle < PHI)
	return MAX_COEF;
  //--------------------------------------
  // compute the coefficient
  // pi/2 is the max for our formula
  if (minangle > PI_OVER_TWO)
	minangle = PI_OVER_TWO;
  // get the types of the external point
  DFloat coeff = theSmoothCoeff[basicPoint[p1]];
  coeff *= theSmoothCoeff[basicPoint[p2]];
#ifdef  __VERBAL4
  cout << endl  << " giving a coeeficient of " << coeff; cout.flush();
#endif
  // evaluate formula  
   coeff *= 1 / (minangle - PHI + 1.0);
   coeff *= (distance(point(p1), point(p2)) / totalBoundaryLength);
#ifdef  __VERBAL4
  cout << endl << " the Min Angle  of :" << p1 << ", " << p2 << " is : " << minangle * 180 / M_PI 
	   << " it contributes: " <<  1 / (minangle - PHI + 1.0)
	   << " to the final coeff: " << coeff << endl; cout.flush();
#endif
  // coefficient is computed
  //--------------------------------------
    return coeff;
}//smoothingCoef2()


//---------------------------------------------------------------------------
// returns true if the point is used by another separator
//---------------------------------------------------------------------------
bool
madd::separatorPointUsed(LInt thePoint)
{
  // go through the first boundary point
  for (int i=0; i < noOfSeparatorTriangles; i++){
	if ((seperatorPoint1[i] == thePoint)  &&  (i != sepTriangleIndex))
	  return true;
  }
  // go through the first boundary point
  for (int i=0; i < noOfSeparatorTriangles; i++){
	if ((seperatorPoint2[i] == thePoint)  &&  (i != sepTriangleIndex))
	  return true;
  }
  // not found, we can use it
  return false;
}//separatorPointUsed()


//---------------------------------------------------------------------------*
//                                                                           *
//     END OF THE MADD ROUTINES                                              *
//                                                                           *
//---------------------------------------------------------------------------*





				
  
  

//---------------------------------------------------------------------------*
//                                                                           *
//     START OF THE SUBDOMAIN ROUTINES                                       *
//                                                                           *
//---------------------------------------------------------------------------*


//---------------------------------------------------------------------------
// initialize theSubDomain as our work subdomain
//---------------------------------------------------------------------------
void
madd::initSubdomain(int theSubDomain)
{
  assignToWorkSubdomain(theSubDomain);
  createPointToWorkSegm();
}//initSubdomain()


//---------------------------------------------------------------------------
// make theSubDomain the workSubdomain
//---------------------------------------------------------------------------
void
madd::assignToWorkSubdomain(int theSubDomain)
{
  workSubdomain = theSubDomain;
  workSegm = subdomainSegments[theSubDomain];
  workNoOfSegms = noOfSubSegms[theSubDomain];
}//assignToWorkSubdomain()


//---------------------------------------------------------------------------
// create pointers from pointS to the work segments
// (used in identifying adjacency)
//---------------------------------------------------------------------------
void
madd::createPointToWorkSegm()
{ 
  LInt segm;
  LInt point1, point2;

  // clear all pointers
   for (LInt i = 0; i < workNoOfSegms; i++){
	 LInt segm = workSegm[i];
	 toWorkSegm1[segm_point1(segm)] = -1;
	 toWorkSegm1[segm_point2(segm)] = -1;
	 toWorkSegm2[segm_point1(segm)] = -1;
	 toWorkSegm2[segm_point2(segm)] = -1;
   }

#ifdef __VERBAL2
   cout << endl << endl << "Creating toWorkSegm list"; cout. flush();
#endif
   // assign the work segments to the points
   for (LInt i = 0; i < workNoOfSegms; i++){
	 segm =workSegm[i];
	 point1 = segm_point1(segm);
	 point2 = segm_point2(segm);
	 //     cout << endl << i << ":  " << point1 << ",  " << point2; cout. flush();
	 if (toWorkSegm1[point1] < 0)
	   toWorkSegm1[point1] = i;
	 else
	   toWorkSegm2[point1] = i;

	 if (toWorkSegm1[point2] < 0)
	   toWorkSegm1[point2] = i;
	 else
	   toWorkSegm2[point2] = i;
   }// for all work segments
#ifdef __VERBAL2
   cout << endl << "done  " << endl << endl; cout. flush();
#endif
}//createPointToWorkSegm()


//---------------------------------------------------------------------------
// put into the workHolesList the holes of theSubDomain
//---------------------------------------------------------------------------
void
madd::createWorkHoles(int theSubDomain)
{
	 uchar* subHole = subdomainHoles[theSubDomain];
	 workNoOfHoles = 0;
	 for (int i=0; i < noOfHoles; i++){
	   if (subHole[i] > 0){
		 // the hole is included
		 workHolesList[workNoOfHoles << 1] = hole_x(i);
		 workHolesList[(workNoOfHoles << 1) + 1] = hole_y(i);
		 workNoOfHoles++;
	   }//if
	 }// for all holes
}//createWorkHoles()


//---------------------------------------------------------------------------
// create the initial subdomain 0
//---------------------------------------------------------------------------
void
madd::createInitialSubdomain()
{
  noOfSubdomains = 1;
  // assign all the segments
  noOfSubSegms[0] = noOfSegms;
  subdomainSegments[0] = (LInt*)malloc((noOfSubSegms[0])*sizeof(LInt) + 1);
  LInt* tmpSubSegment =  subdomainSegments[0];
  for (LInt i = 0; i < noOfSubSegms[0]; i++)
	tmpSubSegment[i] = i;

  // get subdomain values
  subdomainArea[0] =  subdomainMinArea[0] =  subdomainLength[0]  = totalBoundaryLength;    
  subdomainWeight[0] = subdomainWeight[1] = 0;

  // check the holes of the first domain
  // initially include all holes
  for (int i=0; i < noOfHoles; i++){
	 parentSubdomainHoles[i] = 1;
  }
  findSubdomainHoles(0);

  // initialize background nodes
  subdomainBackgroundStart[0] = 0;
  subdomainBackgroundEnd[0] = backgroundListSize;
  if (POINT_WEIGHT_COEF == 0 && BACKGROUND_WEIGHTS)
	subdomainBackgroundEnd[0] = 0;

  // get the weight from the background nodes
  getBackgroundWeights(0, 1);

  //  DFloat minv, maxv, totv;
  //  evaluateGrid(0, &minv, &maxv, &totv);
#ifdef __DEBUG_MODE
  writePolyFile((char *)"init.poly", -1, 'n', 'y', 'n', -1, -1);
  writePolyFile((char *)"init0.poly", 0, 'n', 'n', 'n', -1, -1);
#endif
}//createInitialSubdomain()


//---------------------------------------------------------------------------
// create  a new set of points and segments from the work segments (of some subdomain)
//---------------------------------------------------------------------------
LInt
madd::createPointsSegms(DFloat **newPointAddress, // to the new point list
						LInt **newSegmAddress,    // to the new segment list
						LInt **fromPointAddress)  // to the reference back to the original points
{  
  DFloat *newPoint; // the new point list
  LInt *newSegm, *fromPoint; // the new segment list, reference back to the original points

   // allocate memory
   LInt estPoints = workNoOfSegms + 6;
   newPoint = (DFloat*)malloc(estPoints*2*sizeof(DFloat)+1);
   fromPoint = (LInt*)malloc(estPoints*sizeof(LInt)+1);
   newSegm = (LInt*)malloc(workNoOfSegms*2*sizeof(LInt)+1);
   LInt *toPoint = (LInt*)malloc(noOfPoints*sizeof(LInt)+1); // temporary reference to the new points
   if (newPoint == NULL || newSegm == NULL || fromPoint == NULL || toPoint == NULL){
	   InsufficientMem((char *)"createWorkPoints");
   }
   //clear the temporary reference to the new points
   for (LInt i=0; i < noOfPoints; i++){
	 toPoint[i] = -1;
   }

   // assign values to the new points and segments
   LInt p;
   LInt pointNo = 0;
   for (LInt i=0; i < workNoOfSegms; i++){
	 LInt segm = workSegm[i];
	 p = segm_point1(segm);
	 if (toPoint[p] < 0){
	   // create new point and reference
	   toPoint[p] = pointNo;
	   fromPoint[pointNo] = p;
	   newPoint[(pointNo << 1)] = point_x(p);
	   newPoint[(pointNo << 1)+1] = point_y(p);
	   pointNo++;
	 }//if
	 newSegm[(i << 1)] =  toPoint[p]; // this is the first point of the new segment
	 p = segm_point2(segm);
	 if (toPoint[p] < 0){
	   // create new point and reference
	   toPoint[p] = pointNo;
	   fromPoint[pointNo] = p;
	   newPoint[(pointNo << 1)] = point_x(p);
	   newPoint[(pointNo << 1)+1] = point_y(p);
	   pointNo++;
	 }//if
	 newSegm[(i << 1)+1] = toPoint[p];// this is the second point of the new segment

	 // check if we have enough memory
	 if (pointNo >= estPoints - 2){
	   estPoints += 60;
	   newPoint = (DFloat*)realloc(newPoint, estPoints*2*sizeof(DFloat)+1);
	   fromPoint = (LInt*)realloc(fromPoint, estPoints*sizeof(LInt)+1);
	   if (newPoint == NULL ||  fromPoint == NULL ){
		 InsufficientMem((char *)"createWorkPoints");
	   }//if
	 }//if

   }// for all workSegms

   free(toPoint);
   *newPointAddress = newPoint;
   *newSegmAddress = newSegm;
   *fromPointAddress = fromPoint;
   return pointNo;   
}//createPointsSegms()

	 
//---------------------------------------------------------------------------*
//                                                                           *
//     END OF THE SUBDOMAIN ROUTINES                                         *
//                                                                           *
//---------------------------------------------------------------------------*




//---------------------------------------------------------------------------*
//                                                                           *
//     START OF THE MEMORY ALLOCATION ROUTINES & CONSTRUCTORS                *
//                                                                           *
//---------------------------------------------------------------------------*


//---------------------------------------------------------------------------
// the constructor
//---------------------------------------------------------------------------
madd::madd()
{ 
  setAllZero();
}


//---------------------------------------------------------------------------
// the destructor
//---------------------------------------------------------------------------
madd::~madd()
{ // destructor
  clearAll();
}

//---------------------------------------------------------------------------
// print an error message and exit
//---------------------------------------------------------------------------
void   
madd::ErrorExit(char message[])
{                            
  cout << endl << endl << "ERROR: "
   << message << ". Exit Program." << endl; cout.flush();
  exit(1);
}// ErrorExit                           


//---------------------------------------------------------------------------
// error message Insufficient memory, exit
//---------------------------------------------------------------------------
void   
madd::InsufficientMem(char where[])
{  
  cout << endl << endl << "ERROR: Insufficient memory in "
	   << where << ". Exit Program." << endl; cout.flush();
  exit(1);
}//InsufficientMem


//---------------------------------------------------------------------------
// clear the partition structures
//---------------------------------------------------------------------------
void   
madd::clearPartition()
{
  clearTriangleList();  // clear the triangulation
}//clearPartition()


//---------------------------------------------------------------------------
// clear all structures
//---------------------------------------------------------------------------
void
madd::clearAll()
{ 
  // clear main structures
  free(thePointsList);
  free(theSegmsList);
  free(segmColourList);
  free(theHolesList);
  // clear subdomains
  free(noOfSubSegms);
  for (int i = 0; i < noOfSubdomains; i++){
	free(subdomainSegments[i]);
  }
  free(subdomainArea);
  free(subdomainLength);
  // clear the gradation structures
  free(backgroundNodeList); 
  free(backgroundAreaList);
  free(subdomainBackgroundStart);
  free(subdomainBackgroundEnd);
  // clear triangulation
  clearTriangleList(); 
  // reset everybody
  setAllZero();  
}//clearAll()



//---------------------------------------------------------------------------
// reset all variables
//---------------------------------------------------------------------------
void
madd::setAllZero()
{
  // reset basic structures
  segmColourList = NULL;
  theTriaPointList = thePointsList = theHolesList = NULL;
  theTriaList = theSegmsList = NULL;
  triaNeighborList = NULL;
  parentSubdomainHoles = NULL;
  toWorkSegm1 = toWorkSegm2 = NULL;

  noOfPoints = noOfSegms = noOfHoles = 0;
  TRIA_SIZE = 0;
  THE_SIZE = 0;
  POINTLIST_SIZE = SEGMLIST_SIZE = HOLELIST_SIZE = 0;

  // reset subdomain variables
  subdomainArea = subdomainLength = NULL;
  noOfSubdomains = totalNoOfSubdomains = 0;
  noOfSubSegms = NULL;
  subdomainSegments = NULL;
  subdomainBackgroundStart = subdomainBackgroundEnd = NULL;
  // reset gradation variables
  backgroundNodeList = NULL;
  backgroundAreaList = NULL;
  backgroundListSize = 0;
  libraryHandle = NULL;
  weightFunction = NULL;
  gridFunction = NULL;
  areaFunction = NULL;
}//setAllZero()  


//---------------------------------------------------------------------------
// clear the trianguleation
//---------------------------------------------------------------------------
void
madd::clearTriangleList()
{ 
#ifdef __VERBAL3
  cout << endl << "clear triangles..."; cout.flush();
  cout << endl << "free theTriaPointList:" << theTriaPointList; cout.flush();
#endif
  free(theTriaPointList);
#ifdef __VERBAL3
  cout << endl << "free theTriaList:" << theTriaList; cout.flush();
#endif
  free(theTriaList);
  free(triaNeighborList);

  theTriaPointList = NULL;
  theTriaList = NULL;
  triaNeighborList = NULL;
  noOfTria = 0;
#ifdef __VERBAL3
  cout << endl << "done. "; cout.flush();
#endif
}//clearTriangleList()


//---------------------------------------------------------------------------
// allocate the subdomains structures
//---------------------------------------------------------------------------
void
madd::allocateSubdomainStr()
{
  noOfSubSegms = (LInt*)malloc((totalNoOfSubdomains+1)*sizeof(LInt) + 1);
  subdomainSegments = (LInt**)malloc((totalNoOfSubdomains+1)*sizeof(LInt*) + 1);
  subdomainHoles = (uchar**)malloc((totalNoOfSubdomains+1)*sizeof(uchar*) + 1);
  subdomainArea = (DFloat*)malloc((totalNoOfSubdomains+1)*sizeof(DFloat) + 1);
  subdomainLength = (DFloat*)malloc((totalNoOfSubdomains+1)*sizeof(DFloat) + 1);
  subdomainMinArea = (DFloat*)malloc((totalNoOfSubdomains+1)*sizeof(DFloat) + 1);
  subdomainWeight = (DFloat*)malloc((totalNoOfSubdomains+1)*sizeof(DFloat) + 1);
  subdomainBackgroundStart = (unsigned int*)malloc((totalNoOfSubdomains+1)*sizeof(unsigned int) + 1);
  subdomainBackgroundEnd = (unsigned int*)malloc((totalNoOfSubdomains+1)*sizeof(unsigned int) + 1);
  if (noOfSubSegms == NULL || subdomainSegments == NULL 
	  || subdomainArea == NULL || subdomainLength == NULL 
	  || subdomainHoles == NULL || subdomainWeight == NULL || subdomainMinArea == NULL
	  || subdomainBackgroundStart == NULL || subdomainBackgroundEnd == NULL)
	InsufficientMem((char *)"Initialize Subdomains");
  // create all subdomain hole flags
  for (int i=0; i < totalNoOfSubdomains; i++){
	if ( (subdomainHoles[i] = (uchar*)malloc(noOfHoles*sizeof(uchar)+1)) == NULL)
      InsufficientMem((char *)"Initialize subdomainHoles");
  }
  if ( (parentSubdomainHoles = (uchar*)malloc(noOfHoles*sizeof(uchar)+1)) == NULL){
	InsufficientMem((char *)"Initialize subdomainHoles");
  }
}//allocateSubdomainStr()



//---------------------------------------------------------------------------
// allocate the point structures
//---------------------------------------------------------------------------
void
madd::allocatePoints(LInt size)
{

  size += 600;  // keep 600 places for future use

  if (size <= POINTLIST_SIZE) // no need to reallocate
	return;

#ifdef __VERBAL2
  cout << endl << "Allocating  " << size << " points..." ; cout.flush();
#endif
#ifdef __VERBAL2
  cout << endl << " basicPoint: " << (void*)basicPoint; cout.flush();
#endif

  if (thePointsList == NULL){
	// this is the initial allocation
#ifdef __VERBAL2
	cout << endl << "Allocating thePointsList: " << " -> "; cout.flush();
#endif
	thePointsList = (DFloat*)malloc(size*2*sizeof(DFloat)+1);
#ifdef __VERBAL2
	cout << thePointsList << endl << "Allocating toWorkSegm1: "; cout.flush();
#endif
	toWorkSegm1 = (LInt*)malloc(size*sizeof(LInt)+1);
#ifdef __VERBAL2
	cout << toWorkSegm1 << endl << "Allocating toWorkSegm2: "; cout.flush();
#endif
	toWorkSegm2 = (LInt*)malloc(size*sizeof(LInt)+1);
#ifdef __VERBAL2
	cout << toWorkSegm2 << endl << "Allocating basicPoint: "; cout.flush();
#endif
	basicPoint = (uchar*)malloc(size*sizeof(uchar)+1);
#ifdef __VERBAL2
	cout << (void*)basicPoint << endl << "Allocating pointWeight: "; cout.flush();
#endif
	pointWeight = (float*)malloc(size*sizeof(float)+1);
#ifdef __VERBAL2
	cout <<  (void*)pointWeight << endl << "Allocating finished" ; cout.flush();
#endif
  }
  else{  // reallocate
#ifdef __VERBAL2
	cout << endl << "Reallocating thePointsList: " << thePointsList << " -> "; cout.flush();
#endif
	thePointsList = (DFloat*)realloc(thePointsList, size*2*sizeof(DFloat)+1);
#ifdef __VERBAL2
	cout << thePointsList << endl << "Reallocating toWorkSegm1: " << toWorkSegm1 << " -> "; cout.flush();
#endif
	toWorkSegm1 = (LInt*)realloc(toWorkSegm1, size*sizeof(LInt)+1);
#ifdef __VERBAL2
	cout << toWorkSegm1 << endl << "Reallocating toWorkSegm2: " << toWorkSegm2 << " -> "; cout.flush();
#endif
	toWorkSegm2 = (LInt*)realloc(toWorkSegm2, size*sizeof(LInt)+1);
#ifdef __VERBAL2
	cout << toWorkSegm2 << endl << "Reallocating basicPoint: " << (void*)basicPoint << " -> "; cout.flush();
#endif
	basicPoint = (uchar*)realloc(basicPoint, size*sizeof(uchar)+1);
#ifdef __VERBAL2
	cout << (void*)basicPoint << endl << "Reallocating pointWeight: " <<  (void*)pointWeight << " -> "; cout.flush();
#endif
	pointWeight = (float*)realloc(pointWeight, size*sizeof(float)+1);
#ifdef __VERBAL2
	cout <<  (void*)pointWeight << endl << "Reallocating finished" ; cout.flush();
#endif
  }

  if (thePointsList == NULL || toWorkSegm1 == NULL || toWorkSegm2 == NULL 
	  || basicPoint == NULL || pointWeight == NULL){
	InsufficientMem((char *)"PointList");
	exit(1);
  }
#ifdef __VERBAL2
  cout << endl << " basicPoint: " << (void*)basicPoint << "   pointWeight: " << (void*)pointWeight; cout.flush();
  cout << endl << "done."; cout.flush();
#endif

  POINTLIST_SIZE = size;
}//allocatePoints()


//---------------------------------------------------------------------------
// allocate the background grid structure
//---------------------------------------------------------------------------
void
madd::allocateBackgroundNodes(unsigned int theSize)
{

  free(backgroundNodeList);
  free(backgroundAreaList);
  backgroundListSize = theSize;

  backgroundNodeList = (DFloat*)malloc(theSize * 2 * sizeof(DFloat)+1);
  backgroundAreaList = (float*)malloc(theSize * sizeof(float)+1);

  if (backgroundAreaList == NULL || backgroundNodeList == NULL){
	InsufficientMem((char *)"backgroundAreaList");
	exit(1);
  }
}//allocateBackgroundNodes()


//---------------------------------------------------------------------------
// allocate the holes structure
//---------------------------------------------------------------------------
void
madd::allocateHoles(LInt size)
{

  size += 2;

  if (size <= HOLELIST_SIZE)
	return;

  if (theHolesList == NULL){
	theHolesList = (DFloat*)malloc(size*2*sizeof(DFloat)+1);
	workHolesList = (DFloat*)malloc(size*2*sizeof(DFloat)+1);
	//	 noOfholeSegm = (LInt*)malloc(size*sizeof(LInt)+1);
	//	 holeColourList = (uchar*)malloc(size*sizeof(uchar)+1);
  }
  else{
	theHolesList = (DFloat*)realloc(theHolesList, size*2*sizeof(DFloat)+1);
	workHolesList = (DFloat*)realloc(workHolesList, size*2*sizeof(DFloat)+1);
	// noOfholeSegm = (LInt*)realloc(noOfholeSegm, size*sizeof(LInt)+1);
	//	 holeColourList = (uchar*)realloc(holeColourList, size*sizeof(uchar)+1);
  }
  if (theHolesList == NULL || workHolesList == NULL){ 
	// || holeColourList == NULL){
	InsufficientMem((char *)"theHolesList");
	exit(1);
  }

  HOLELIST_SIZE = size;
}//allocateHoles()


//---------------------------------------------------------------------------
// allocate the segemnts structure
//---------------------------------------------------------------------------
void
madd::allocateSegms(LInt size)
{

  size += 600;  // keep 600 places for future use

  if (size <=  SEGMLIST_SIZE)
	return;

#ifdef __VERBAL2
  cout << endl << "Allocating  " << size << " segments..." ; cout.flush();
#endif

  if (theSegmsList == NULL){
	//	 theSegmsList = (LInt*)malloc(size*2*sizeof(LInt)+1);
	theSegmsList = (LInt*)malloc(size * 2*sizeof(LInt)+1);
	//	 segmentIncludesHole = (LInt*)malloc(size*sizeof(LInt)+1);
	segmColourList = (uchar*)malloc(size*sizeof(uchar)+1);
	//	 separatorFlagList = (uchar*)malloc(size*sizeof(uchar)+1);
  }
  else{
#ifdef __VERBAL2
	cout << endl << "Reallocating theSegmsList: " << theSegmsList << " -> "; cout.flush();
#endif
	theSegmsList = (LInt*)realloc(theSegmsList, (size*2*sizeof(LInt)+1));
	//	 segmentIncludesHole = (LInt*)realloc(segmentIncludesHole, size*sizeof(LInt)+1);
#ifdef __VERBAL2
	cout << theSegmsList << endl << "Reallocating segmColourList..." << segmColourList << " -> " ; cout.flush();
#endif
	segmColourList = (uchar*)realloc(segmColourList, (size*sizeof(uchar)+1));
#ifdef __VERBAL2
	//	 cout << segmColourList << endl << "Reallocating separatorFlagList..." << separatorFlagList << " -> " ; cout.flush();
#endif
	//	 separatorFlagList = (uchar*)realloc(separatorFlagList, (size*sizeof(uchar)+1));
#ifdef __VERBAL2
	cout << endl << "Reallocating finished" ; cout.flush();
#endif
  }

  if (theSegmsList == NULL || segmColourList == NULL){
	InsufficientMem((char *)"theSegmsList");
	exit(1);
  }

#ifdef __VERBAL2
  cout << "done."; cout.flush();
#endif
  SEGMLIST_SIZE = size;
}//allocateSegms


//---------------------------------------------------------------------------*
//                                                                           *
//     END OF THE MEMORY ALLOCATION ROUTINES                                 *
//                                                                           *
//---------------------------------------------------------------------------*



//---------------------------------------------------------------------------*
//                                                                           *
//     INSERT POINTS, EDGES AND HOLES TO THE STRUCTURE                       *
//                                                                           *
//---------------------------------------------------------------------------*


//---------------------------------------------------------------------------
//  insert the segment between point1, point2
//  return the index of this segment
//---------------------------------------------------------------------------
LInt   
madd::NewSegm(LInt point1, LInt point2)
{
   if (noOfSegms >= SEGMLIST_SIZE)
      allocateSegms(noOfSegms);
  
   segm_point1(noOfSegms) = point1;
   segm_point2(noOfSegms) = point2;

   // the points become basic
   basicPoint[point1] =  1;
   basicPoint[point2] =  1;

   return (++noOfSegms)-1;
}


//---------------------------------------------------------------------------
//  insert a point at (x, y) with weight =  weight
//  returns the index of this point
//---------------------------------------------------------------------------
LInt   
madd::NewPoint(DFloat x, DFloat y, float weight)
{
   if (noOfPoints >= POINTLIST_SIZE)
      allocatePoints(noOfPoints);
  
   point_x(noOfPoints) = x;
   point_y(noOfPoints) = y;
   basicPoint[noOfPoints] = 0;
   pointWeight[noOfPoints] = weight;
   return (++noOfPoints)-1;
}


//---------------------------------------------------------------------------
//  insert a point at (x, y) with weight = 0
//  returns the index of this point
//---------------------------------------------------------------------------
LInt   
madd::NewPoint(DFloat x, DFloat y)
{
   if (noOfPoints >= POINTLIST_SIZE)
      allocatePoints(noOfPoints);
  
   point_x(noOfPoints) = x;
   point_y(noOfPoints) = y;
   basicPoint[noOfPoints] = 0;
   pointWeight[noOfPoints] = 0;
   return (++noOfPoints)-1;
}


//---------------------------------------------------------------------------
//  insert a hole at (x, y) 
//  returns the number of holes so far
//---------------------------------------------------------------------------
LInt    
madd::NewHole(DFloat x, DFloat y)
{
   if (noOfHoles >= HOLELIST_SIZE)
      allocateHoles(noOfHoles + 5);

    hole_x(noOfHoles) = x;
    hole_y(noOfHoles) = y;
    return (noOfHoles++);
}


//---------------------------------------------------------------------------*
//                                                                           *
//     END OF INSERT POINTS, EDGES AND HOLES TO THE STRUCTURE                *
//                                                                           *
//---------------------------------------------------------------------------*



//---------------------------------------------------------------------------*
//                                                                           *
//     START OF THE GRADING CONTROL ROUTINES                                 *
//                                                                           *
//---------------------------------------------------------------------------*



//---------------------------------------------------------------------------
// evaluate the subdomain over a uniform background gird
// using the gridFunction
// The grid is creted on the fly
//---------------------------------------------------------------------------
int
madd::evaluateGrid(int theSubdomain,   // the subdomain to be evaluated
				   DFloat *minValue,   // return the min value
				   DFloat *maxValue,   // return the max value
				   DFloat *meanValue,  // return the mean value
				   DFloat *totalValue) // return the sum
{
#define boundaryGrid_x(i) boundaryGridPoints[((i) << 1)]    // the x cordinate of the boundaryGridPoint
#define boundaryGrid_y(i) boundaryGridPoints[((i) << 1)+1]  // the y cordinate of the boundaryGridPoint

  //----------------------------------------
  //  create the grid bounds from the boundary 
  //  of  theSubdomain
  //----------------------------------------
  DFloat PRECISION_CHECK = GRID_STEP / 10000;  // for precision check if two points are too close 
                                               // with respect to the grid step
  DFloat *boundaryGridPoints;   // it will hold the boundaryGridPoints
  char *sideFlag;               // will hold on which side (1=left, 2=right, 3=middle) of the segment 
                                // is the boundary grid point 
  int noOfboundaryGridPoints = 0;
  int totalboundaryGridPoints = 0;
  LInt segm, p1, p2;  // a segment and the two points
  
  // make theSubdomain the work subdomain
  assignToWorkSubdomain(theSubdomain); 
  // compute the number of boundaryGridPoints
  for (LInt theWorkSegm = 0; theWorkSegm < workNoOfSegms; theWorkSegm++){
	segm =  workSegm[theWorkSegm];
	DFloat x1 = point_x(segm_point1(segm));
	DFloat x2 = point_x(segm_point2(segm));
	totalboundaryGridPoints += ((int)(fabs(x2 - x1) / GRID_STEP)) + 1;
  }
  // allocate the boundaryGridPoints lists
  boundaryGridPoints = (DFloat*)malloc((totalboundaryGridPoints+100)*2*sizeof(DFloat));
  sideFlag = (char*)malloc((totalboundaryGridPoints+100)*sizeof(char));
  if (boundaryGridPoints == NULL || sideFlag == NULL){
	InsufficientMem((char *)"boundaryGridPoints");
	exit(1);
  }

  // go through the boundary and find the 
  for (LInt theWorkSegm = 0; theWorkSegm < workNoOfSegms; theWorkSegm++){
	segm =  workSegm[theWorkSegm];

	// get the segment points in  x-y order
	if (point_x(segm_point1(segm)) < point_x(segm_point2(segm))){
	  p1 = segm_point1(segm);
	  p2 = segm_point2(segm);
	}
	else{
	  if (point_x(segm_point1(segm)) > point_x(segm_point2(segm))){
		p1 = segm_point2(segm);
		p2 = segm_point1(segm);
	  }
	  else
		// the x are equal
		if (point_y(segm_point1(segm)) < point_y(segm_point2(segm))){
		  p1 = segm_point1(segm);
		  p2 = segm_point2(segm);
		}
		else{
		  p1 = segm_point2(segm);
		  p2 = segm_point1(segm);
		}
	}

	// find line equation for p1, p2
	DFloat a, b, c;
	findLineEquation(point(p1), point(p2), &a, &b, &c);

	DFloat x1 = point_x(p1);
	DFloat y1 = point_y(p1);
	DFloat x2 = point_x(p2);
	DFloat y2 = point_y(p2);
	// get the x of the starting boundary grid point
	DFloat runx = ((long long)(x1 / GRID_STEP)) * GRID_STEP; 
	// make sure it's inside
	if (x1 - runx > PRECISION_CHECK)
	  runx += GRID_STEP;

#ifdef  __VERBAL2
	cout << endl << " Finding grid points for segm: (" << x1 << ", " << y1 << "),  (" << x2 << ", " << y2
		 << ").   Step is " << GRID_STEP << ",   starting x is " << runx << endl; cout.flush();
#endif
 
	if (runx - x2 > PRECISION_CHECK){ // no grid points on this chap
	  continue;
	}

	if (fabs(x1 - x2) <= PRECISION_CHECK){
	  // the grid is on this line (interesting)
	  // we will get the grid from the adjacent segments
	  continue;
	}

	// get start position
	int startPosition = 0;
	// get grid start points on current segment
	while (runx - x2 < PRECISION_CHECK){
	  DFloat runy = - (((a * runx) + c) / b); // this is a point on the segment, so we should get reasonable result
	  // find insert position in the array (sort-merge)
	  while ((startPosition < noOfboundaryGridPoints) 
			 && (runx - boundaryGrid_x(startPosition) > PRECISION_CHECK)) 
		startPosition++;
#ifdef  __VERBAL2
	  cout << endl << "in StartPosition " << startPosition << " we have (" 
		   << boundaryGrid_x(startPosition)  << ", " << boundaryGrid_y(startPosition) << ")"; 
#endif
	  while ((startPosition < noOfboundaryGridPoints) 
			 && (fabs(boundaryGrid_x(startPosition) - runx) <= PRECISION_CHECK ) 
			 && (runy - boundaryGrid_y(startPosition) > PRECISION_CHECK))
		startPosition++;
#ifdef  __VERBAL2
	  cout << endl << "in StartPosition " << startPosition << " we have (" 
		   << boundaryGrid_x(startPosition)  << ", " << boundaryGrid_y(startPosition) << ")"; 
	  //	  if (boundaryGrid_x(startPosition) == runx)
	  //		cout << "The x is the same";
	  // else
	  //	cout << "The x is NOT the same";
#endif

	  // check the side
	  char thisSide;
	  if (fabs(runx - x1) < PRECISION_CHECK)
		thisSide = 1; //it's on the left
	  else
		if (fabs(runx - x2) < PRECISION_CHECK)
		  thisSide = 2;//it's on the right
		else
		  thisSide = 3;//it's in the middle

	  if ((fabs(boundaryGrid_x(startPosition) - runx) < PRECISION_CHECK) 
		  && (fabs(boundaryGrid_y(startPosition) - runy) < PRECISION_CHECK)){
		//----------------------------------------
		// this is the same point
		if (thisSide == sideFlag[startPosition]){
		  // they are on the same side, get rid of it.
		  for (int gridpoint = startPosition+1; gridpoint < noOfboundaryGridPoints; gridpoint++){
			boundaryGrid_x(gridpoint-1) =  boundaryGrid_x(gridpoint);
			boundaryGrid_y(gridpoint-1) =  boundaryGrid_y(gridpoint);
		  }
		  for (int gridpoint = startPosition+1; gridpoint < noOfboundaryGridPoints; gridpoint++){
			sideFlag[gridpoint-1] = sideFlag[gridpoint]; 
		  }	  
		  noOfboundaryGridPoints--;		  
		}
		runx += GRID_STEP;
		continue; 
		//----------------------------------------
	  }// if the same

	  // shift the arrays
	  for (int gridpoint = noOfboundaryGridPoints; gridpoint > startPosition; gridpoint--){
		boundaryGrid_x(gridpoint) =  boundaryGrid_x(gridpoint-1);
		boundaryGrid_y(gridpoint) =  boundaryGrid_y(gridpoint-1);
	  }
	  for (int gridpoint = noOfboundaryGridPoints; gridpoint > startPosition; gridpoint--){
		sideFlag[gridpoint] = sideFlag[gridpoint-1]; 
	  }	  
	  noOfboundaryGridPoints++;

	  // insert grid point
#ifdef  __VERBAL2
	  cout << endl << " inserting grid point (" << runx << ", " << runy << "),  "; 
	  cout << endl << " startPosition is " << startPosition;
	  if (startPosition > 0)
		cout << endl << " previous point is  (" <<  boundaryGrid_x(startPosition-1) <<", " 
			 <<  boundaryGrid_y(startPosition-1) << ")."; 
	  if (startPosition+1 < noOfboundaryGridPoints)
		cout << endl << " next point is  (" <<  boundaryGrid_x(startPosition+1) <<", " 
			 <<  boundaryGrid_y(startPosition+1) << ")."; 
	  cout.flush();
#endif
	  boundaryGrid_x(startPosition) = runx;
	  boundaryGrid_y(startPosition) = runy;
	  sideFlag[startPosition] = thisSide;
		
	  startPosition++;
	  runx += GRID_STEP;
	}// for grid points on segment

  }//for all work segments
  //----------------------------------------
  //  the boundary grid is created 
  //----------------------------------------

#ifdef  __VERBAL2
  cout << endl << " The boundary grid list is created." << endl; cout.flush();
#endif
  // we have the grid boundary array sorted on x-y
  // scan and evaluate function
  DFloat theMin, theMax, theTotal;
  theTotal = 0;
  // get initial min-max
  if (gridFunction != NULL)
	theMin = theMax = gridFunction(boundaryGrid_x(0), boundaryGrid_y(0));
#ifdef __VERBAL2
  cout << endl << "-> the grid value at (" << boundaryGrid_x(0) << ", " << boundaryGrid_y(0)
	   << ") is " << theMin << endl;
#endif

  //----------------------------------------
  // go through all boundary grid points
  // get the vertical lines of the grid
  // and evalutate the contained grid points
  //----------------------------------------
  int boundpoint = 0;
  while( boundpoint < noOfboundaryGridPoints-1){
	DFloat runx =  boundaryGrid_x(boundpoint);
	DFloat runy;
#ifdef  __VERBAL2
	cout << endl << " evaluating grid line at (" << boundaryGrid_x(boundpoint) << ", " << boundaryGrid_y(boundpoint) << "),  (" 
		 <<  boundaryGrid_x(boundpoint+1) << ", " << boundaryGrid_y(boundpoint+1) << ")"; cout.flush();
#endif
	// check if the two consequent points are on the same vertical line
	if (fabs(runx - boundaryGrid_x(boundpoint+1)) > PRECISION_CHECK){
	  // it does not pair! Get the next.
	  boundpoint++;
	  continue;
	}

	DFloat functionValue; 
	// shortcut for getting the min and the max
#define updateGridValue() 	{ theMin = min(theMin, functionValue);\
	                          theMax = max(theMax, functionValue);}
    // get the start - end y coordinates
	DFloat y1 = boundaryGrid_y(boundpoint);
	DFloat y2 = boundaryGrid_y(boundpoint+1);

	// evaluate the grid-boundary inersection 
	runy = y1;
	functionValue =  gridFunction(runx, runy);
	updateGridValue();
	runy = y2;
	functionValue =  gridFunction(runx, runy);
	updateGridValue();

	//----------------------------------------
	// go through the grid points of the vertical line
	runy = ((long long)(y1 / GRID_STEP)) * GRID_STEP;
	if (runy <= y1)
	  runy += GRID_STEP;
#ifdef  __VERBAL2
	cout << endl << "   y starts at " << runy; cout.flush();
#endif
	while (runy < y2){
	  functionValue =  gridFunction(runx, runy);
	  updateGridValue();
#ifdef  __VERBAL2
	  cout << endl << " value at grid point (" << runx << ", " << runy << "),  is " << functionValue; cout.flush();
#endif
#ifdef __WRITE_GRID
	  // insert the grid points for test
	  NewPoint(runx, runy);
#endif
	  runy += GRID_STEP;
	}// while
	// go through the grid points done
	//----------------------------------------
	  
	boundpoint += 2; // get next pair
				
  }// for all bounadary grid points 
  //----------------------------------------
  //  the grid points have been evalutated
  //----------------------------------------

  // return the values  
  *minValue = theMin;
  *maxValue = theMax;
  *totalValue = theTotal;
  free(boundaryGridPoints);
  free(sideFlag);
  // done
  return 0;
}// evaluateGrid()


//---------------------------------------------------------------------------
//  assigns the subdomain weights and min-areas from the background grid
//  to the two created subdomains from the parent domain
//---------------------------------------------------------------------------
void
madd::getBackgroundWeights(int parentSubdomain,  // thye parent domain, also the first subdomain
						   int subdomain2)       // the second subdomain
{
  // the  first subdomain
  int subdomain1 = parentSubdomain;
  // The start-end indeces of the background points for the subdomains
  unsigned int startWeight = subdomainBackgroundStart[subdomain1];
  unsigned int endWeight = subdomainBackgroundEnd[subdomain1];
  unsigned int startWeight1 = startWeight;
  unsigned int endWeight1 =  endWeight;
  // the two min-areas, start from the initial values
  DFloat min1 = subdomainMinArea[subdomain1];
  DFloat min2 = subdomainMinArea[subdomain2];
  DFloat max1 = 0;
  DFloat max2 = 0;
  DFloat total1 = 0;
  DFloat total2 = 0;
#ifdef  __VERBAL2
  cout << endl << "------------------------------------------";
  cout << endl << "getBackgroundWeights....";
  cout << endl << "startWeight: " << startWeight 
	   << ",   endWeight: " << endWeight;
  cout	<< endl << "Initial Weights: " << subdomainWeight[subdomain1] << ",  " << subdomainWeight[subdomain2];
#endif

  //----------------------------------------
  // Evaluate the background grid
  // and assign to the two subdomains
  // the included  background grid points
  //----------------------------------------
  if (startWeight < endWeight){
	// there are background points
	startWeight1 = startWeight;
	endWeight1 =  endWeight;
	// get work domain from subdomain1
	initSubdomain(subdomain1);
	int count=1;
	// go through the included background points
	// and locate in which subdomain they belong
	for (unsigned int i = startWeight1; i < endWeight1;){
	  // the weight or area of the node
	  DFloat thisWeight =  weightOfNode(i);
	  // find the  subdomain
	  if (containsPoint(backgroundNode(i))){
		// background node included in subdomain1
		// get the weight and the area 
		//	subdomainWeight[subdomain1] += thisWeight;
		min1 = min(min1, thisWeight);
		max1 = max(max1, thisWeight);
		total1 += thisWeight;
		// coninue to the next background node
		i++;
	  }
	  else{
		// background node included in subdomain2
		// get the weight and the area 
		//		subdomainWeight[subdomain2] += thisWeight;
		min2 = min(min2, thisWeight);
		max2 = max(max2, thisWeight);
		total2 += thisWeight;
		//----------------------------------------
		// move this background node to the end
		endWeight1--;
		DFloat mx = background_x(endWeight1);
		DFloat my = background_y(endWeight1);
		background_x(endWeight1) = 	background_x(i);
		background_y(endWeight1) = 	background_y(i);
		background_x(i) = mx;
		background_y(i) = my;
		backgroundAreaList[i] = backgroundAreaList[endWeight1];
		backgroundAreaList[endWeight1] = thisWeight;
		//----------------------------------------
	  }//else
	}// for all included  background points
  }// if
  //----------------------------------------
  // background grid is evaluated
  //----------------------------------------
  
  // assign the min areas
  subdomainMinArea[subdomain1] =  min1;
  subdomainMinArea[subdomain2] =  min2;
  subdomainWeight[subdomain1] =  total1;
  subdomainWeight[subdomain2] =  total2;
  subdomainWeight[subdomain1] =  max1;
  subdomainWeight[subdomain2] =  max2;
  // assign the start-end indeces of the
  // background grid points for the
  // two subdomains
  subdomainBackgroundStart[subdomain1] = startWeight;
  subdomainBackgroundEnd[subdomain1] = endWeight1;
  subdomainBackgroundStart[subdomain2] = endWeight1;
  subdomainBackgroundEnd[subdomain2] = endWeight;

  //----------------------------------------
  // background grid is done
  //----------------------------------------
#ifdef  __VERBAL2
  cout << endl << "Subdomain " << subdomain1 << ": "
	   << "startWeight: " << subdomainBackgroundStart[subdomain1]   
	   << ",  endWeight: " << subdomainBackgroundEnd[subdomain1]   
	   << ",  Weight: " << subdomainWeight[subdomain1];   
  cout << endl << "Subdomain " << subdomain2 << ": "
	   << "startWeight: " << subdomainBackgroundStart[subdomain2]   
	   << ",  endWeight: " << subdomainBackgroundEnd[subdomain2]   
	   << ",  Weight: " << subdomainWeight[subdomain2];   
  cout << endl << "------------------------------------------"; cout.flush();
#endif
  return;
}//evaluateGrid()

//---------------------------------------------------------------------------*
//                                                                           *
//     END OF THE GRADING CONTROL ROUTINES                                   *
//                                                                           *
//---------------------------------------------------------------------------*




//---------------------------------------------------------------------------*
//                                                                           *
//     START OF THE GEOMETRY ROUTINES                                        *
//                                                                           *
//---------------------------------------------------------------------------*


// short cuts
// returns the y coordinate for y on the line defined be (x1, y1), (x2, y2)
#define find_y(x1, y1, x2, y2, x) (x2 == x1)? x1: ((y1 - y2) * (x2 - x) / (x2 - x1) + y2)


//---------------------------------------------------------------------------
// check if a partial inner separator intersects (internally) the boundary
//---------------------------------------------------------------------------
bool
madd::intersectsBoundary(DFloat *center,  // the circumcenter
						 LInt checkPoint) // the boundary point
{
  LInt segm, p1, p2;
  // go through all the subdomain segments
  // and check if they intersect the partial separator
  for (LInt i = 0; i < workNoOfSegms; i++){
	segm =  workSegm[i];
	p1 = segm_point1(segm);
	p2 = segm_point2(segm);
	if ( p1 != checkPoint && p2 != checkPoint){
	  // the segment does not contain the checkPoint
	  if (segmentsIntersect(center, point(checkPoint), point(p1), point(p2)) >= 0)
		// the two segments intersect
		return true;
	}
  }//for

  // done, no inersection
  return false;
}//intersectsBoundary()


//---------------------------------------------------------------------------
// check if a partial outer separator intersects (internally) the boundary
//---------------------------------------------------------------------------
bool
madd::intersectsBoundary(LInt checkPoint1,  // the first boundary point
						 LInt checkPoint2)  // the second boundary point
{
  LInt segm, p1, p2;
  // go through all the subdomain segments
  // and check if they intersect the partial separator
  for (LInt i = 0; i < workNoOfSegms; i++){
	segm =  workSegm[i];
	p1 = segm_point1(segm);
	p2 = segm_point2(segm);
	if ( p1 != checkPoint1 && p2 != checkPoint1 && p1 != checkPoint2 && p2 != checkPoint2 ){
	  // the segment is not adjacent to the partial outer separator
	  if (segmentsIntersect(point(checkPoint1), point(checkPoint2), point(p1), point(p2)) >= 0)
		// the two segments intersect
		return true;
	}
  }//for
   
  // done, no inersection
  return false;
}//intersectsBoundary()


//---------------------------------------------------------------------------
//  return the work segemnt defined by two points
//---------------------------------------------------------------------------
LInt
madd::findWorkSegm(LInt p1,    // the first point 
				   LInt p2)    // the second point 
{
  // get the first segment from p1
  LInt segm = toWorkSegm1[p1];
  //  cout << endl << "finding toWorkSegm of: " << p1 << ",  " << p2; cout.flush();
  // cout << endl << "the work segments: " << toWorkSegm1[p1]  << ",  " << toWorkSegm2[p1]
  //	   << "  and   " << toWorkSegm1[p2]  << ",  " << toWorkSegm2[p2]; cout.flush();

  // check if segm works 
  if (segm == toWorkSegm1[p2] || segm == toWorkSegm2[p2])
	return segm;

  // get the second segment from p1
  segm = toWorkSegm2[p1];
  // check if segm works 
  if (segm == toWorkSegm1[p2] || segm == toWorkSegm2[p2])
	return segm;

  // the connectivity list didn't find it
  // go through the whole worklist
  for (LInt i = 0; i < workNoOfSegms; i++){
	segm = workSegm[i];
	if ((segm_point1(segm) == p1 || segm_point2(segm) == p1)
		&& (segm_point1(segm) == p2 || segm_point2(segm) == p2))
	  return i;
  }

  // there is something wrong here, no segment found
  cout << endl << "Cannot find segment : " << p1 << ", " << p2; cout.flush();
  ErrorExit((char *)"findWorkSegm didn't find the segment");
  return -1;
}//findWorkSegm()


//---------------------------------------------------------------------------
//  return the point "left" of a point
//---------------------------------------------------------------------------
LInt
madd::leftOfPoint(LInt point)
{
  // get the "left" segment 
  LInt segm = workSegm[toWorkSegm1[point]];
  // if segm < 0 there is no adjacent point
  if (segm < 0)
	return point;

  // find the other point of segm
  // and return it
  if (segm_point1(segm) != point)
	return segm_point1(segm);

  return segm_point2(segm);
}//leftOfPoint()


//---------------------------------------------------------------------------
//  return the point "right" of a point
//---------------------------------------------------------------------------
LInt
madd::rightOfPoint(LInt point)
{
  // get the "right" segment 
  LInt segm = workSegm[toWorkSegm2[point]];
  // if segm < 0 there is no "right" point
  if (segm < 0)
	return point;

  // find the other point of segm
  // and return it
  if (segm_point1(segm) != point)
	return segm_point1(segm);

  return segm_point2(segm);
}//rightOfPoint()


//---------------------------------------------------------------------------
//  return the next point of the sequence point1, point2
//---------------------------------------------------------------------------
LInt
madd::nextOfPoints(LInt point1,   // the first point 
				   LInt point2)   // the second point 
{
  LInt nextPoint;

  //----------------------------------------
  // get the first segment of point2
  LInt segm = workSegm[toWorkSegm1[point2]];
  // if segm < 0 no adjacent point
  if (segm < 0)
	return point2;

  // suppose next point is on the "left"
  nextPoint = segm_point1(segm);
  // if it's not the first two, it's the next
  if (nextPoint != point1 && nextPoint != point2 )
	return nextPoint;

  // suppose next point is on the "right"
  nextPoint = segm_point2(segm);
  // if it's not the first two, it's the next
  if (nextPoint != point1 && nextPoint != point2 )
	return nextPoint;

  //----------------------------------------
  // get the second segment of point2
  segm = workSegm[toWorkSegm2[point2]];
  // if segm < 0 no adjacent point
  if (segm < 0)
	return point2;

  // suppose next point is on the "left"
  nextPoint = segm_point1(segm);
  if (nextPoint != point1 && nextPoint != point2 )
	return nextPoint;

  // suppose next point is on the "right"
  nextPoint = segm_point2(segm);
  if (nextPoint != point1 && nextPoint != point2 )
	return nextPoint;

  // we shouldn't be here
  return point2;
}//nextOfPoints()	


//---------------------------------------------------------------------------
//  return the adjacent numberOfPoints-1 points of a point
//---------------------------------------------------------------------------
void
madd::nextPointsOfPoint(LInt point,          // the point for which we want the adjacent
						LInt *nextPoint,     // a list of the ajacent points in order
                                             // in the middle is the given point
						int numberOfPoints)  // the total number of points
{
  // initialize the points = point
  for (int i=0; i < numberOfPoints; i++)
	nextPoint[i] = point;
  // get the middle index
  int middle = (numberOfPoints / 2);
  // get the first "right" point
  LInt rightPoint = nextPoint[middle+1] = rightOfPoint(point);
  // get the first "left" point
  LInt leftPoint = nextPoint[middle-1] = leftOfPoint(point);
  //----------------------------------------
  // from the middle outwards compute the next points in the sequence
  //----------------------------------------
  for (int i = 2; i <= middle; i++){
	if (leftPoint == rightPoint){
	  //----------------------------------------
	  // we are in a circle, stop here
	  while(i <= middle){
		nextPoint[middle-i] = nextPoint[middle+i] = leftPoint;
		i++;
	  }// while
	  //----------------------------------------
	}//if
	else{
	  // no circle so far
	  // get the next points left (down), right (up)
	  leftPoint = nextPoint[middle-i] = nextOfPoints(nextPoint[middle-i+2], leftPoint);
	  rightPoint = nextPoint[middle+i] = nextOfPoints(nextPoint[middle+i-2], rightPoint);
	}
  }// for
  //----------------------------------------
  // the adjacent points have been computed
  //----------------------------------------
  /*
	cout << endl << "Adjacent points: "; 
	for (int i = 0; i < numberOfPoints; i++){
	cout  << (nextPoint[i]) << ", ";
	}
	cout.flush();
  */
}//nextPointsOfPoint()


//---------------------------------------------------------------------------
//  return the segment "left" of a segment
//---------------------------------------------------------------------------
LInt
madd::leftOfSegm(LInt aSegm)
{
  // get the segment from the work segment
  LInt curSegm = workSegm[aSegm];

  // get the first point
  LInt point = segm_point1(curSegm);
  // get the two segments of  the first point
  LInt segm1 = toWorkSegm1[point];
  LInt segm2 = toWorkSegm2[point];

  // return the adjacent segment
  if (segm1 != aSegm && segm1 >= 0)
	return segm1;
  if (segm2 != aSegm && segm2 >= 0)
	return segm2;

  // we shouldn't be here
  ErrorExit((char *)"Left segment was not found");
  return -1;
}//leftOfSegm()


//---------------------------------------------------------------------------
//  return the segment "right" of a segment
//---------------------------------------------------------------------------
LInt
madd::rightOfSegm(LInt aSegm)
{
  // get the segment from the work segment
  LInt curSegm = workSegm[aSegm];

  // get the second point
  LInt point = segm_point2(curSegm);
  // get the two segments of  the second point
  LInt segm1 = toWorkSegm1[point];
  LInt segm2 = toWorkSegm2[point];

  // return the adjacent segment
  if (segm1 != aSegm && segm1 >= 0)
	return segm1;
  if (segm2 != aSegm && segm2 >= 0)
	return segm2;

  // we shouldn't be here
  ErrorExit((char *)"Right segment was not found");
  return -1;
}//rightOfSegm()


//---------------------------------------------------------------------------
//  returns true if the work subdomain contains aPoint
//---------------------------------------------------------------------------
bool
madd::containsPoint(DFloat *aPoint)
{	 
// define for this function a precision limit 
#define PRECISION_TOLERANCE 0.000001
  // the two coordinates of the point
  register DFloat pointx = x(aPoint);
  register DFloat pointy = y(aPoint);
  DFloat x1, x2, y1, y2, y;

  // counter switch for the number of segments above aPoint
  register char included = 0;  // if == 0 we have even number of segments
                               // if == 1 we have even number of segments 

#ifdef __VERBAL3
  int xpass = 0;
  cout << endl << "Includes check for point (" << pointx << ", " << pointy << ")." << endl; cout.flush();
#endif
  //----------------------------------------
  // count how many subdomain segments are above aPoint
  // if it's even, aPoint is outside
  // if it's odd, aPoint is inside
  //----------------------------------------
  // go through all the segments 
  for (LInt w = 0; w < workNoOfSegms; w++){
	LInt segm = workSegm[w];
	// get the coordinates of the two end points
	x1 = point_x(segm_point1(segm));
	y1 = point_y(segm_point1(segm));
	x2 = point_x(segm_point2(segm));
	y2 = point_y(segm_point2(segm));
	// put the two points in order
	if (x1 > x2){
	  register DFloat temp = x1;
	  x1 = x2;
	  x2 = temp;
	  temp = y1;
	  y1 = y2;
	  y2 = temp;
	}

#ifdef __VERBAL3
	cout << endl << " checking segment (" << x1 << ", " << y1 << ") - (" << x2 << ", " << y2 << "): "; cout.flush();
#endif
	// check if aPoint is too close to this segment	
	if (fabs(x1 - pointx) < PRECISION_TOLERANCE && fabs(y1 - pointy) < PRECISION_TOLERANCE)
	  return true; // accepted being very close
	if (fabs(x2 - pointx) < PRECISION_TOLERANCE && fabs(y2 - pointy) < PRECISION_TOLERANCE)
	  return true; // accepted being very close

	// check if the segment's x projection includes the 
	// x projection of aPoint 
	if ((x1 - pointx) > PRECISION_TOLERANCE || (pointx > x2))
	  continue;  // out of x-interval, not above
	if ((x2 - pointx) < PRECISION_TOLERANCE)
	  continue; // ignore the count of segments at the far right point

#ifdef __VERBAL3
	cout << "   x-passed = " << ++xpass; cout.flush();
	cout << endl << ++xpass << "   " << x1 << "  " << y1;
	cout << endl << ++xpass << "   " << x2 << "  " << y2; cout.flush();
#endif
	//----------------------------------------
    // aPoint is in the x-interval of the segment
	// check if the segment is above
	// get the projection of aPoint to the segment
	y = find_y(x1, y1, x2, y2, pointx);
	if (fabs(y - pointy) < PRECISION_TOLERANCE)
	  return true; // accepted as too close to segment
    // check if the projection os above
	if (y > pointy){
	  // we have a hit, change counting flag
	  included ^= 1;
	}
	//----------------------------------------
#ifdef __VERBAL3
	cout << " included = " << (int)included; cout.flush();
#endif
  }// for all segments
  // count how many subdomain segments are above aPoint is done
  // if included == 0 we have even number of segments
  // if included == 1 we have odd number of segments
  //----------------------------------------

#ifdef __VERBAL3
  cout << endl << " included = " << (int)included << endl; cout.flush();
#endif

  if (included == 0){
	// we have even number of segments above, not included
	return false;
  }

  // we have odd number of segments above, included
  return true;
  //----------------------------------------
#undef PRECISION_TOLERANCE
}//containsPoint()


//---------------------------------------------------------------------------
//  assigns the holes to the subdomain 
//  the parent holes are in parentSubdomainHoles[]
//---------------------------------------------------------------------------
void
madd::findSubdomainHoles(LInt subdomain)
{
#ifdef __VERBAL2
  cout << endl << "define subdomain holes for ..."<< subdomain; cout.flush();
#endif
  // get subdomain in the work structure
  initSubdomain(subdomain);
  // the holes list of subdomain
  uchar* workSubHole = subdomainHoles[subdomain];

  //----------------------------------------
  // go through the holes and see if subdomain includes them	
  for (int i=0; i < noOfHoles; i++){
	if (parentSubdomainHoles[i] > 0)
	  // this was a parent hole, check it
	  workSubHole[i] = checkThisHole(i); 
	else
	  // this was not a parent hole, reject it
	  workSubHole[i] = 0; // no hole 
  }// for all holes
  //----------------------------------------

}//findSubdomainHoles()



//---------------------------------------------------------------------------
//  an interface for checking if a hole belongs to the work subdomain
//  if not included return 0
//  if included return 1
//---------------------------------------------------------------------------
uchar
madd::checkThisHole(LInt theHole)
{	 
#ifdef __VERBAL2
  cout << endl << "Checking Hole: " << theHole << "..."; cout.flush();
#endif

  return containsHole(hole(theHole));
}//checkThisHole()


//---------------------------------------------------------------------------
//  check if a hole belongs to the work subdomain
//  if not included return 0
//  if included return 1
//---------------------------------------------------------------------------
uchar
madd::containsHole(DFloat *aPoint)
{	
  // the coordinates of the hole point
  DFloat holex = x(aPoint);
  DFloat holey = y(aPoint);
  // the  hole point
  DFloat* thisHole = aPoint;

  DFloat x1, x2, y1, y2, y;
  LInt segm;

  //----------------------------------------
  // find the closest segment on the top
  // of the   hole point
  //----------------------------------------
  LInt topSegm = -1;    // initially no top segment
  DFloat mindist = -1;  // initially the min distance is -1
  // go through all the work segments
  for (LInt w = 0; w < workNoOfSegms; w++){
	LInt j = workSegm[w];
	// get the coordinates of the end points
	x1 = point_x(segm_point1(j));
	y1 = point_y(segm_point1(j));
	x2 = point_x(segm_point2(j));
	y2 = point_y(segm_point2(j));
	// check if the segment's x range covers the hole point
	if ((x1 - holex) * (x2 - holex) <= 0){
	  // we have a hit
	  y = find_y(x1, y1, x2, y2, holex);
	  // get the y-height
	  DFloat dist = (y - holey);
	  //		   cout << endl << " we have a top seg hit." << endl 
	  //		<< "   x1: " << x1 << "   x2: " << x2 << "    holex: " << holex << endl  
	  //	<< "   y1: " << y1 << "   y2: " << y2 << "    holey: " << holey << endl  
	  //	<< "   dist: " << dist << "   mindist " << mindist << endl; cout.flush();  

	  // if the y-height is >= 0 and is closer (or th first) 
	  // grab the segment
	  if ((dist >= 0) && (mindist < 0 || mindist > dist)){
		//	 cout << "       Accepted." << endl; cout.flush();  
		mindist = dist;
		topSegm = w;
	  }//if
	}//if
  }// for all the work segments
  //----------------------------------------
  // the closest segment on the top is found
  //----------------------------------------

  // if top segement not found, retrun 0 (not contained)
  if (topSegm < 0)
	return 0;

#ifdef __VERBAL2
  cout << endl << "   Top segment found: " << topSegm << ". Checking connectivity, angles:" << endl; cout.flush();
#endif

  //----------------------------------------
  // find the rest of the segments and compute 
  // the total angle
  //----------------------------------------
  segm = workSegm[topSegm];
  DFloat totalAngle = sangle(point(segm_point1(segm)), aPoint, point(segm_point2(segm)));
  LInt nextSegm = rightOfSegm(topSegm);
  LInt previousSegm = topSegm;
  LInt currentSegm = nextSegm;
  LInt counter = 1;
 
  //  while we haven't gone a full circle get the angle for next segment
  while (counter < noOfSegms && nextSegm != topSegm){
#ifdef __VERBAL2
	cout << "    " << nextSegm << ": " << totalAngle << ", "; cout.flush(); 
#endif
	// shift to the current segment
	currentSegm = nextSegm;
    segm = workSegm[currentSegm];
	// find the orientation
	if ( rightOfSegm(currentSegm) != previousSegm){
	  // it's on the right
	  // compute angle
	  totalAngle += sangle(point(segm_point1(segm)), aPoint, point(segm_point2(segm)));
	  // get next segment
	  nextSegm = rightOfSegm(currentSegm);
	}
	else{
	  // it's on the left
	  // compute angle
	  totalAngle += sangle(point(segm_point2(segm)), aPoint, point(segm_point1(segm)));
	  // get next segment
	  nextSegm = leftOfSegm(currentSegm);
	}//else
	// shift the current segment
	previousSegm = currentSegm;
	counter++; // the counter will force exit in a DEAD Loop
  }// while
  //----------------------------------------
  // the total angle is computed
  //----------------------------------------
#ifdef __VERBAL2
  cout << endl << "End of connectivity ";
  cout << "   noOfSegs: " << noOfSegms << "   counter: " << counter 
	   << "   angle: " << totalAngle << "   2 PI : " << TWO_PI_MINUS << endl; cout.flush();  
#endif
  // if the angle is less than 2*pi, the hole is not inside
  if (counter > noOfSegms || fabs(totalAngle) < TWO_PI_MINUS){
	// no hole
	return 0;
  }

  // the angle is 2*pi, the hole is inside
  return 1;
}//containsHole()


//---------------------------------------------------------------------------
//  refine all segments according to a distance dist
//---------------------------------------------------------------------------
void  
madd::breakAllSegm(DFloat dist,  uchar basicFlag)
{
  LInt endSegm = noOfSegms;
  for (LInt i=0; i < endSegm; i++) {
	breakSegm(i, dist, basicFlag);
  }
  return;
}//breakAllSegm()


//---------------------------------------------------------------------------
//  refine all segments into number of pieces
//---------------------------------------------------------------------------
void  
madd::breakAllSegm(int pieces,  uchar basicFlag)
{
  LInt endSegm = noOfSegms;
  for (LInt i=0; i < endSegm; i++) {
	breakSegm(i, pieces, basicFlag);
  }
  return;
}//breakAllSegm()


//---------------------------------------------------------------------------
//  refine one segment into subsegmants of size ...size
//  in weighted gradation, it takes into account the weights and the gradation factor
//---------------------------------------------------------------------------
void 
madd::breakSegm(LInt theSegm, DFloat size,  uchar basicFlag)
{
  if (AREA_WEIGHT_COEF <= 0)
	return;

  LInt n1, n2;
  DFloat dist;
  // get the two points
  n1 = segm_point1(theSegm);
  n2 = segm_point2(theSegm);
  // get the weights
  DFloat weight1 = pointWeight[n1];
  DFloat weight2 = pointWeight[n2];

  // compute the refine coefficient
  DFloat refCoef;
  if ((weight1 + weight2) <= 0)
	refCoef = 1;
  else
	refCoef = 1 / ( ((weight1 + weight2) * POINT_WEIGHT_COEF / (2*maxPointWeight)) + 1); 

  // compute the weighted length 
  size *= refCoef * (1 /  (DFloat)AREA_WEIGHT_COEF);

  // if the distance is small do not break 
  if ((dist = distance(point(n1), point(n2))) <= (size * 1.1))
	return;

  //  cout << endl << "Breaking segment: " << theSegm << " (" << n1 << ", " << n2
  // << "),   with distance: " << dist << "  into size: " << size; cout.flush();

  // compute the number of subsegments
  int pieces = (int)((dist / size) + 0.8);

  // break theSegm into pieces 
  if (pieces > 1)
	breakSegm(theSegm, pieces,  basicFlag);

  // done
  return;
}//breakSegm()



//---------------------------------------------------------------------------
//  refine the inserted set of separators
//---------------------------------------------------------------------------
int 
madd::refineThisSepLength()
{
  DFloat size = refineLength; 
  // go hrough all the new separators
  for (LInt j=madSepStart; j < madSepEnd; j++) {
	breakSegm(j, size, 0);
  }
  madSepEnd = noOfSegms; 
  return 0;
}


//---------------------------------------------------------------------------
//  refine one segment into subsegments
//  in the weighted case a graded refinement will be produced
//---------------------------------------------------------------------------
void 
madd::breakSegm(LInt theSegm,     // the segment to be refined
				int pieces,       // the number of pieces 
				uchar basicFlag)  // if > 0, this will be assigned to the new points
                                  // else 0 will be assigned, except to the middle point
{
  if (pieces < 2) {
	return;
  }

  // we want even number of pieces
  if ((pieces & 1) == 1)
	pieces++; // or pieces-- ?

#ifdef __VERBAL2
  cout << endl << "Breaking  segment " << theSegm << " into " << pieces << " pieces..."; cout.flush();
#endif
  LInt n1, n2;
  // the middle point index
  int middlePoint = (pieces)/2;
  // the end points of the segment
  n1 = segm_point1(theSegm);
  n2 = segm_point2(theSegm);
  // the weights of the end points
  DFloat weight1 = pointWeight[n1];
  DFloat weight2 = pointWeight[n2];

#ifdef __VERBAL3
  cout << endl << "points and weights found: " << n1 << "(" << weight1 << "),  " << n2 << "(" << weight2 << ")"; cout.flush();
#endif

  // check if reallocation is needed
  if (noOfPoints + pieces + 5 >=  POINTLIST_SIZE)
	allocatePoints(noOfPoints + pieces);
  if (noOfSegms + pieces + 5 >= SEGMLIST_SIZE)
	allocateSegms(noOfSegms + pieces);

  //----------------------------------------
  // create the new points and segments
  //----------------------------------------

  // the end point coordinates
  register DFloat  x1 = point_x(n1);
  register DFloat  y1 = point_y(n1);
  register DFloat  x2 = point_x(n2);
  register DFloat  y2 = point_y(n2);
  // the refining coefficient for the garded case
  register DFloat breakCoef;
  // the next point location
  register DFloat* nextPoint = point(noOfPoints);
  nextPoint--;
  // the next segment location
  register LInt* nextSegm = segm(noOfSegms);
  nextSegm--;
  // the next basic flag location
  register uchar* nextBasicPoint = basicPoint+(noOfPoints-1);
  // the next point weight location
  register float* nextPointWeight = pointWeight+(noOfPoints-1);
  // the left point of the subsegment
  register LInt leftPoint = n1;
  // the right point of the subsegment
  register LInt rightPoint = noOfPoints;
  // the new point coordinates
  DFloat newX;
  DFloat newY;
#ifdef __VERBAL2
  cout << endl << "Inserting points and segms..."; cout.flush();
#endif
  for (register int i = 1; i < pieces; i++) {
	//  assign the cordinates of the new point
	*(++nextPoint) = newX = (x2-x1)*(breakCoef = (DFloat)i/(DFloat)pieces) + x1;
	*(++nextPoint) = newY = (y2-y1)*breakCoef + y1;
	// assign the weight  of the new point
	pointWeight[rightPoint] = (uint)(((weight2 - weight1) * breakCoef + weight1) * INTERPOLATION_WEIGHT_COEF) ;
	if (weightFunction != NULL)
	  pointWeight[rightPoint] += weightFunction(newX, newY);
	// assign the basicFlag of the new point
	if (basicFlag != 0)
	  *(++nextBasicPoint) = basicFlag;
	else
	  // if it's in the middle is 2, otherwise 0
	  *(++nextBasicPoint) = (i == middlePoint)? 2: 0;
    // create the subsegment	
	*(++nextSegm) = leftPoint;
	leftPoint = *(++nextSegm) = rightPoint++;
  }//for 

  // assign the new list sizes   
  noOfPoints += (pieces-1);
  noOfSegms += (pieces-1);

  // replace the old segment
  segm_point1(theSegm) = noOfPoints-1;
  segm_point2(theSegm) = n2;

  //----------------------------------------
  // the new points and segments have been created
  //----------------------------------------
#ifdef __VERBAL2
  cout << "done. "; cout.flush();
  cout << endl << pieces << " segments created."; cout.flush();
#endif
  return;
}//breakSegm()


//---------------------------------------------------------------------------
//  compute the domain lengths
//---------------------------------------------------------------------------
void
madd::computeDomainLengths()
{
  DFloat *p1 = point(segm_point1(0));
  DFloat *p2;
  // the minimum x and the maximum x coordinate
  DFloat  minx, maxx;
  minx = maxx = x(p1);
  // the minimum y and the maximum y coordinate
  DFloat  miny, maxy;
  miny = maxy = y(p1);
  // the total boundary length
  DFloat totalLength = 0;

  // go through all the boundary segments and compute:
  for (LInt i=0; i < noOfSegms; i++) {
	p1 = point(segm_point1(i));
	p2 = point(segm_point2(i));
	// the min x coordinate	
	minx = min(min(x(p1), x(p2)), minx); 
	// the max x coordinate	
	maxx = max(max(x(p1), x(p2)), maxx); 
	// the min y coordinate	
	miny = min(min(y(p1), y(p2)), miny); 
	// the max y coordinate	
	maxy = max(max(y(p1), y(p2)), maxy); 
	// the total boundary length
	totalLength += distance(p1, p2);
  }
  // assign the values to the global variables
  totalBoundaryLength = totalLength;
  domainMinX = minx;
  domainMaxX = maxx;
  domainMinY = miny;
  domainMaxY = maxy;

  return;
}


//---------------------------------------------------------------------------
//  return the total length of the separator of  subdomain
//---------------------------------------------------------------------------
double
madd::separatorLengthOfSubdomain(int subdomain)
{
  DFloat *p1, *p2;
  LInt segm;

  // the total separator length
  DFloat totalLength = 0;
  // make subdomain the work subdomain 
  assignToWorkSubdomain(subdomain);
  // go through all the subdomain segments
  for (LInt i = 0; i < workNoOfSegms; i++){
	// if the segmenmt is a seprator add the length
	if (segm = workSegm[i] >= noOfExtSegms){
	  p1 = point(segm_point1(segm));
	  p2 = point(segm_point2(segm));
	  totalLength += distance(p1, p2);
	}//if
  }//for all the subdomain segments
  return totalLength;
}//separatorLengthOfSubdomain()


//---------------------------------------------------------------------------
//  return the total area of the subdomain
//---------------------------------------------------------------------------
DFloat
madd::areaOfSubdomain(int theSubDomain)
{
  initSubdomain(theSubDomain);
  // create the work holes
  createWorkHoles(theSubDomain);
  // get the Delauany triangulation
  triangulateDomain(NULL);
  DFloat totalArea =0;
  for (LInt thisTria=0; thisTria < noOfTria; thisTria++) {
	// the trhee triangle vertices
	DFloat* point1 = point(tria_point1(thisTria));
	DFloat* point2 = point(tria_point2(thisTria));
	DFloat* point3 = point(tria_point3(thisTria));
	// add the total triangle area
	totalArea += areaOfTriangle(point1, point2, point3);
  }
  return totalArea;
}//subdomainArea()


//---------------------------------------------------------------------------
//  return the maximum weight of the boundary points
//---------------------------------------------------------------------------
float
madd::computeMaxPointWeight()
{
  float maxWeight = 0;
  if (weightFunction != NULL){
	for (LInt i=0; i < noOfPoints; i++) {
	  pointWeight[i] += weightFunction(point_x(i), point_y(i));
	  if (pointWeight[i] > maxWeight)
		maxWeight = pointWeight[i];
	}
  }
  else{
	for (LInt i=0; i < noOfPoints; i++) {
	  if (pointWeight[i] > maxWeight)
		maxWeight = pointWeight[i];
	}
  }

  return maxWeight;
}//computeMaxPointWeight()


//---------------------------------------------------------------------------
// returns true if the angles from  point1 to the boundPoint are >= PHI
//---------------------------------------------------------------------------
bool
madd::goodExtAngle(DFloat* point1, LInt boundPoint)
{
  // check the first angle
  if  (angle(point1, point(boundPoint), point(leftOfPoint(boundPoint))) < PHI)
	return false;
  // check the second angle
  if  (angle(point1, point(boundPoint), point(rightOfPoint(boundPoint))) < PHI)
	return false;

  // angles accepted
  return true;
}//goodExtAngle()


//---------------------------------------------------------------------------
// returns the minimum of the four angles formed by 
// the separator between two boundary points 
//---------------------------------------------------------------------------
DFloat
madd::minExtAngle(LInt boundPoint1,  // the first boundary point
				  LInt boundPoint2)  // the second boundary point
{

  DFloat minangle, angle1, angle2, angle3, angle4;
  //  cout << "Computing external angles of points: " 
  //	   << leftOfPoint(boundPoint1) << ", " << boundPoint1 << ", " << rightOfPoint(boundPoint1)
  //	   << "  and  " << leftOfPoint(boundPoint2) << ", " << boundPoint2 << ", " << rightOfPoint(boundPoint2);
  // cout.flush();   

  // the four angles
  angle1 = angle(point(boundPoint1), point(boundPoint2), point(leftOfPoint(boundPoint2)));
  angle2 = angle(point(boundPoint1), point(boundPoint2), point(rightOfPoint(boundPoint2)));
  angle3 = angle(point(boundPoint2), point(boundPoint1), point(leftOfPoint(boundPoint1)));
  angle4 = angle(point(boundPoint2), point(boundPoint1), point(rightOfPoint(boundPoint1)));
  /*
	cout << endl << " The four external angles: " << angle1 *180 / M_PI << ",  " << angle2 *180 / M_PI
	<< ",  " << angle3 *180 / M_PI  << ",  " << angle4 *180 / M_PI;
  */

  // get the minimum
  minangle = min(angle1, angle2);
  minangle = min(minangle, angle3);
  minangle = min(minangle, angle4);

  return minangle;
}//minExtAngle()


//---------------------------------------------------------------------------*
//                                                                           *
//    BASIC GEOMETRY AND MATH ROUTINES                                       *
//                                                                           *
//---------------------------------------------------------------------------*


//---------------------------------------------------------------------------
// returns the distance between two points 
//---------------------------------------------------------------------------
DFloat 
madd::distance(DFloat* p1, DFloat* p2)
{
   return sqrt(sqr(x(p1) - x(p2)) + sqr(y(p1) - y(p2)));
}


//---------------------------------------------------------------------------
// returns the inner product of the vectors (p1-p0) * (p2 - p0) 
//---------------------------------------------------------------------------
DFloat 
madd::innProduct(DFloat* p1, DFloat* p0, DFloat* p2)
{
      return (x(p1) - x(p0)) * (x(p2) - x(p0)) +  (y(p1) - y(p0)) * (y(p2) - y(p0));
}

   
//---------------------------------------------------------------------------
// returns the area of the triangle p1, p0, p2 
//---------------------------------------------------------------------------
DFloat 
madd::areaOfTriangle(DFloat* p1, DFloat* p0, DFloat* p2)
{
   return fabs((x(p1) - x(p0)) * (y(p2) - y(p0)) - (y(p1) - y(p0)) * (x(p2) - x(p0)));
}
   
   
//---------------------------------------------------------------------------
// returns the signed area of the triangle p1, p0, p2 
//---------------------------------------------------------------------------
DFloat 
madd::sArea(DFloat* p1, DFloat* p0, DFloat* p2)
{
   return ((x(p1) - x(p0)) * (y(p2) - y(p0)) - (y(p1) - y(p0)) * (x(p2) - x(p0)));
}


//---------------------------------------------------------------------------
// returns the cos of the angle p1,p0,p2 
//---------------------------------------------------------------------------
DFloat 
madd::cos(DFloat* p1, DFloat* p0, DFloat* p2)
{
  DFloat denom = distance(p1, p0) * distance(p2, p0);

  if (denom == 0)
	return 1; // the angle is zero

  DFloat thecos = innProduct(p1, p0, p2) / denom;

  if (thecos > 1.0) return 1.0;
  if (thecos < -1.0) return -1.0;

  return thecos;
}


//---------------------------------------------------------------------------
// returns the angle formed by p1,p0,p2 
//---------------------------------------------------------------------------
DFloat 
madd::angle(DFloat* p1, DFloat* p0, DFloat* p2)
{
      DFloat thecos = cos(p1, p0, p2);
      return acos(thecos); 
}


//---------------------------------------------------------------------------
// returns the signed angle formed by p1,p0,p2 
//---------------------------------------------------------------------------
DFloat 
madd::sangle(DFloat* p1, DFloat* p0, DFloat* p2)
{
      DFloat thecos = cos(p1, p0, p2);
	  DFloat area = sArea(p1, p0, p2);
	  DFloat theangle;
      if (area >= 0)
		theangle = acos(thecos); 
      else
		theangle = -acos(thecos);
	  /*	  
	  cout << endl << "The angle of (" << x(p1) << ", " << y(p1)
		   << "),  (" << x(p0) << ", " << y(p0)
		   << "),  (" << x(p2) << ", " << y(p2)
		   << ")   is: " << theangle; cout.flush();
	  */
	  return theangle;
}


//---------------------------------------------------------------------------
//  computes the circumCenter of the triangle p1, p2, p3 
//---------------------------------------------------------------------------
void
madd::circumCenter(DFloat* p1, DFloat* p2, DFloat* p3,  // the three points of the triangle
				   DFloat* result)                      // will return the circumcenter
{
  DFloat a_x, a_y, b_x, b_y, lamda;

  a_x = x(p2) - x(p1);
  a_y = y(p2) - y(p1);
  b_x = x(p3) - x(p1);
  b_y = y(p3) - y(p1);

  lamda = ((b_x * b_x + b_y * b_y) - (b_x * a_x + b_y * a_y)) * 0.5;
  lamda /= (b_y * a_x - b_x *a_y);

  x(result) = x(p1) + a_x * 0.5 - lamda * a_y;
  y(result) = y(p1) + a_y * 0.5 + lamda * a_x;

  return;
}


//---------------------------------------------------------------------------
//  computes the middle of two points p1, p2 
//---------------------------------------------------------------------------
void 
madd::middle(DFloat* p1, DFloat* p2,    // the two points
			 DFloat* middle)            // will return middle 
{
  x(middle) = (x(p1) + x(p2)) * 0.5;
  y(middle) = (y(p1) + y(p2)) * 0.5;
}


//---------------------------------------------------------------------------
//  checks if two segments defined by p1, p2 and q1, q2 intersect 
//  returns -1 if they dont instersect
//           0 if they are adjacent
//           1 if they intersect
//---------------------------------------------------------------------------
int 
madd::segmentsIntersect(DFloat* p1, DFloat* p2, // the first segment
						DFloat* q1, DFloat* q2) // the second segment
{
  DFloat a1, b1, c1;
  findLineEquation(p1, p2, &a1, &b1, &c1);

  DFloat a2, b2, c2;
  findLineEquation(q1, q2, &a2, &b2, &c2);

  DFloat result11 = a2 * x(p1) + b2 * y(p1) + c2;
  DFloat result12 = a2 * x(p2) + b2 * y(p2) + c2;
  DFloat result1 = result11 * result12;
  if ( result1 > 0)
	return -1;

  DFloat result21 = a1 * x(q1) + b1 * y(q1) + c1;
  DFloat result22 = a1 * x(q2) + b1 * y(q2) + c1;
  DFloat result2 = result21 * result22;
  if ( result2 > 0)
	return -1;

  // check if they are adjacent
  if (result1 * result2 == 0)
	return 0;

  // the segments intersect
  return 1;

}

//---------------------------------------------------------------------------
//  computes the intersection of the lines defined by p1, p2 and q1, q2  
//---------------------------------------------------------------------------
int 
madd::intersectionOf(DFloat* p1, DFloat* p2, // the first two points
					 DFloat* q1, DFloat* q2, // the second two points
					 DFloat* intersecion)    // the intersection point
{
  DFloat a1, b1, c1;
  findLineEquation(p1, p2, &a1, &b1, &c1);

  DFloat a2, b2, c2;
  findLineEquation(q1, q2, &a2, &b2, &c2);

  DFloat inter_x, inter_y;
  if (a1 !=0)
	inter_y = (c1 * a2 - a1 * c2) / (a1 * b2 - b1 * a2);
  else
	inter_y = - (c1 / b1);

  if (a1 != 0)
	inter_x = (-c1 - b1 * inter_y) / a1;
  else
	inter_x = (-c2 - b2 * inter_y) / a2;

  intersecion[0] =inter_x; 
  intersecion[1] =inter_y; 

  return 0;

}


//---------------------------------------------------------------------------
//  computes the equation ax+by+c=0 of the line 
//  defined by two points p1, p2
//---------------------------------------------------------------------------
void 
madd::findLineEquation(DFloat* p1,  DFloat* p2,          // the two points
					   DFloat* a, DFloat* b, DFloat* c)  // the coefficients of the computed equation
{
  DFloat a1, b1, c1;

  if (x(p1) == x(p2)){
	b1 = 0;
	a1 = 1.0;
  }
  else{
	b1 = 1.0;
	a1 = (y(p2) - y(p1)) / (x(p1) - x(p2));
  }

  c1 = - a1 * x(p1) - b1 * y(p1);

  *a = a1;
  *b = b1;
  *c = c1;
  return;
}


//---------------------------------------------------------------------------*
//                                                                           *
//     END OF THE GEOMETRY ROUTINES                                          *
//                                                                           *
//---------------------------------------------------------------------------*





//---------------------------------------------------------------------------*
//                                                                           *
//     FILE IO ROUTINES                                                      *
//                                                                           *
//---------------------------------------------------------------------------*



//---------------------------------------------------------------------------
// writes the domain and all the subdomains to poly files
// returns 0 if succesful
//---------------------------------------------------------------------------
int   
madd::writePolyFileAll(char* fileName)  
{
  char name[32];
  // write the domain
  sprintf(name, "%s.poly", fileName);  
  writePolyFile(name, -1, 'n', 'y', 'n', -1);
  // write the subdomains
  for (int i=0; i < noOfSubdomains; i++){
	sprintf(name, "%s-%d.poly", fileName, i);
	writePolyFile(name, i, 'n', 'n', 'y', -1);
  } 
  return 0;
}//writePolyFileAll()


//---------------------------------------------------------------------------
// writes the whole domain or one subdomain to a poly file
// returns 0 if succesful
//---------------------------------------------------------------------------
int   
madd::writePolyFile(char* name,               // the name of the file
					int subdomain,            // the subdomain, if >= 0 one subdomain
					                          //                if = -1 the decomposition
                                              //                if < -1 the initial domain
					char includeTriangles,    // 'y' = include the triangulation
					char includeAllHoles,     // 'y' = include all holes
					char includeAllPoints,    // 'y' = include all points
					int hole)                 // if >= 0, write only this hole
{
  if (VERBAL > 0){
	cout << endl << "Start writing file : " << name << endl; cout.flush();
  }
  if (noOfPoints >  POINTLIST_SIZE || noOfSegms >  SEGMLIST_SIZE){
	cout << endl << " ERROR: The pointers are out of the lists. Exiting." << endl; cout.flush();
	exit(1);
  }

  DFloat* newPoint;
  LInt  *fromPoint, *newSeg;
  LInt noPoints, noSegm, cnt;

  // open the output file
  ofstream oFile;
  oFile.open(name);
  if (!oFile.good()) {
	cerr << endl << "ERROR IN " << name << endl;
	return -1;
  }
  oFile.precision(8);
  oFile << "# subdomain is " << subdomain;
  //oFile << ",   triangles inluded = " << includeTriangles << endl;
  noPoints = noSegm = 0;

  //----------------------------------------
  // initialize the (sub)domain to write
  if (subdomain >= 0){
	// write a subdomain
	initSubdomain(subdomain);
	noPoints = createPointsSegms(&newPoint, &newSeg, &fromPoint);
	noSegm = workNoOfSegms;
  }
  else{
	if (subdomain == -1){
	  // write the decomposition
	  noPoints = noOfPoints;
	  newPoint = point(0);
	  newSeg = segm(0);
	  fromPoint = NULL;
	  noSegm = noOfSegms;
	  includeAllHoles = 'y';
	}
	else{
	  // write the only the external boundary
	  noPoints = noOfExtPoints;
	  noSegm = noOfExtSegms;
	  newPoint = point(0);
	  newSeg = segm(0);
	  fromPoint = NULL;
	  includeAllHoles = 'y';
	}//else
  }//esle

  //----------------------------------------
  // write the points
  if (VERBAL > 1){
	cout << "     Start writing points... " << endl; cout.flush();
  }
  oFile << "# the points" << endl;
  cnt = 0;
  if (subdomain >= 0){
	// write the subdomain points
	if  (includeAllPoints == 'y'){
	  unsigned int startWeight = subdomainBackgroundStart[subdomain];
	  unsigned int endWeight = subdomainBackgroundEnd[subdomain];
	  oFile << noPoints+(endWeight-startWeight) << "    2  1  0" << endl;
	  for (LInt i=0; i < noPoints; i++) {
		oFile << (++cnt) << "    "  << newPoint[2*i] << "  " <<  newPoint[2*i+1] 
			  << "   " << pointWeight[fromPoint[i]] << endl;
	  }
	  for (LInt i=startWeight; i < endWeight; i++) {
		oFile << (++cnt) << "    "  << background_x(i) << "  " <<   background_y(i) 
			  << "   " << backgroundAreaList[i] << endl;
	  }
	}
	else{
	  oFile << noPoints << "    2  1  0" << endl;
	  for (LInt i=0; i < noPoints; i++) {
		oFile << (++cnt) << "    "  << newPoint[2*i] << "  " <<  newPoint[2*i+1] 
			  << "   " << pointWeight[fromPoint[i]] << endl;
	  }
	}
  }
  else{
	// write all points
	oFile << noPoints << "    2  1  0" << endl;
	for (LInt i=0; i < noPoints; i++) {
	  oFile << (++cnt) << "    "  << newPoint[2*i] << "  " <<  newPoint[2*i+1] 
			<< "   " << pointWeight[i] << endl;
	}
  }

  //----------------------------------------
  // write the segments
  if (VERBAL > 1){
	cout << "     Start writing segments... " << endl; cout.flush();
  }
  cnt = 0;
  oFile << "# the segments" << endl;
  if (includeTriangles == 'y') {
	oFile << noSegm+(noOfTria*3)  << "   0" << endl;
  }
  else
	oFile << noSegm  << "   0" << endl;

  for (LInt i=0; i < noSegm; i++) {
	oFile << (++cnt) << "    "  << newSeg[2*i]+1 << "  " 
		  <<   newSeg[2*i+1]+1 << "   " << endl;
  }

  //----------------------------------------
  // write the holes
  if (VERBAL > 1){
	cout << "     Start writing holes... " << endl; cout.flush();
  }
  oFile << "# the holes" << endl;
  if (includeAllHoles == 'y') {
	// write all holes
	oFile << noOfHoles << endl;
	for (LInt i=0; i < noOfHoles; i++) {
	  oFile << i+1 << "    "  << hole_x(i) << "  " <<  hole_y(i)  << endl;
	}//for
  }
  else {
	if (hole >=0){
	  // write only this hole
	  oFile << 1 << endl;
	  oFile << 1 << "    "  << hole_x(hole) << "  " << hole_y(hole) << endl;
	}
	else{
	  // write the holes of this subdomain
	  createWorkHoles(subdomain);
	  oFile << workNoOfHoles << endl;
	  LInt j = 0;
	  for (LInt i=0; i < workNoOfHoles; i++) {
		oFile << i+1 << "    "  << workHolesList[2*i] << "  " <<  workHolesList[2*i+1]
		  //<< "   0" 
			  << endl;
	  }//for
	}//else
  }//else

  //----------------------------------------
  // check if the file is good and close
  if (!oFile.good()) {
	cerr << endl << "ERROR IN " << name << endl;
	return -1;
  }
  oFile.flush();
  oFile.close();
  // file has been created
  //----------------------------------------
   
  // free the auxiliary strucures
  if (subdomain >= 0){
	free(newPoint);
	free(fromPoint);
	free(newSeg);
  }

  if (VERBAL > 0){
	cout << endl << "End writting file : " << name << endl; cout.flush();
  }

  // file write done
  //----------------------------------------
  return 0;
}//writePolyFile()


//---------------------------------------------------------------------------
// writes the decomposition in points, segments, subdomains -segments
// returns 0 if succesful
//---------------------------------------------------------------------------
int   
madd::writeSubdomainFile(char* name)               // the name of the file
{
  // open the output file
  ofstream oFile;
  oFile.open(name);
  if (!oFile.good()) {
	cerr << endl << "ERROR IN " << name << endl;
	return -1;
  }
  oFile.precision(8);

  //----------------------------------------
  // write the decomposition

  //----------------------------------------
  // write the points
  if (VERBAL > 1){
	cout << "     Start writing points... " << endl; cout.flush();
  }
  oFile << noOfPoints << endl;
  oFile.flush();
  // write all points
  for (LInt i=0; i < noOfPoints; i++) {
	oFile << (i+1) << "   "  << point_x(i) << "  " <<  point_y(i) << endl;
  }

  //----------------------------------------
  // write the segments
  if (VERBAL > 1){
	cout << "     Start writing segments... " << endl; cout.flush();
  }
  oFile << noOfSegms  <<  endl;

  for (LInt i=0; i < noOfSegms; i++) {
	oFile << (i+1) << "    "  << segm_point1(i)+1 << "  " << segm_point2(i)+1 << endl;
  }

  //----------------------------------------
  // write the subdomains
  oFile << noOfSubdomains << endl;
  for (int i=0; i < noOfSubdomains; i++){
	assignToWorkSubdomain(i);
	oFile << (i+1) << "  " << workNoOfSegms << endl;
	// go through all the subdomain segments
	for (LInt j = 0; j < workNoOfSegms; j++){
	  // if the segmenmt is a seprator add the length
	  oFile << workSegm[j]+1 << "  ";
	}//for all the subdomain segments
	oFile << endl;
  }

  //----------------------------------------
  // write the holes of the subdomains
  for (int i=0; i < noOfSubdomains; i++){
	assignToWorkSubdomain(i);
	createWorkHoles(i);
	oFile << (i+1) << "  " << workNoOfHoles << endl;
	for (LInt j=0; j < workNoOfHoles; j++) {
	  oFile << workHolesList[2*j] << "  " <<  workHolesList[2*j+1] << endl;
	}//for all the subdomain segments
  }


  //----------------------------------------
  // check if the file is good and close
  if (!oFile.good()) {
	cerr << endl << "ERROR IN " << name << endl;
	return -1;
  }
  oFile.flush();
  oFile.close();
  // file has been created
  //----------------------------------------

  if (VERBAL > 0){
	cout << endl << "End writting file : " << name << endl; cout.flush();
  }

  // file write done
  //----------------------------------------
  return 0;
}//writeSubdomainFile()




//----------------------------------------
// these are for the read file routines
FILE* iFile;
LInt in1;
DFloat inputs[6];
//----------------------------------------

//---------------------------------------------------------------------------
// reads a poly file
// returns 0 if succesful
//---------------------------------------------------------------------------
int   
madd::readPolyFile( char* name)
{
  if (VERBAL > 0){
	cout << endl << "Start reading file : " << name << endl; cout.flush();
  }

  // open the file
  iFile = fopen(name, "r");
  if (iFile == NULL) {
	cerr << endl << "ERROR open " << name << endl;
	fclose(iFile);  
	return -1;
  }

  //----------------------------------------
  // read the points
  // read the head
  if (readFileLine() < 2) {
	cout << endl << "ERROR reading points head " << name << endl; cout.flush();
	fclose(iFile);  
	return -1;
  }
  LInt pointno = in1;            // the number of poinst to read
  LInt zeroPoint = noOfPoints;   // the structure start position of the new points
  LInt inZeroPoint = 0;          // start from 0 in the file
  int pointAtr = (int)inputs[1]; // the number of attributes
  int readSize = (pointAtr > 0) ? 4: 3; // get the size of the line
  inputs[2] = 0;
  if (VERBAL > 1){
	cout << endl << "Reading " << in1 << " points..."; cout.flush();
  }
  allocatePoints(noOfPoints + pointno);
  for (LInt i = 0; i < pointno; i++) {
	if (readFileLine() < readSize) {
	  cout << endl << "ERROR reading points in " << name << endl; cout.flush();
	  fclose(iFile);  
	  return -1;
	}
	if (i == 0)
	  inZeroPoint = in1;

	NewPoint(inputs[0], inputs[1], (float)inputs[2]);
	//	  cout << inputs[2] << endl;
  }

  //----------------------------------------
  // read the segments
  if (readFileLine() < 1) {
	fclose(iFile);  
	return 1;
  }
  LInt segmno = in1;
  if (VERBAL > 1){
	cout << endl << "reading " << in1 << " segments..."; cout.flush();
  }
  allocateSegms(noOfSegms + segmno);
  for (LInt i = 0; i < segmno; i++) {
	if (readFileLine() < 3) {
	  cerr << endl << "ERROR  reading segments in " << name << endl;
	  fclose(iFile);  
	  return -1;
	}
	//	  cout << endl << i << ": " << inputs[0] << ", " <<  inputs[1]; cout.flush();
	//      NewSegm(LInt(inputs[0])+zeroPoint-1, LInt(inputs[1])+zeroPoint-1, 1, (uchar)inputs[2]);
	NewSegm(LInt(inputs[0])+zeroPoint-inZeroPoint, LInt(inputs[1])+zeroPoint-inZeroPoint);
  }

  //----------------------------------------
  // read the holes
  if (readFileLine() < 1) {
	fclose(iFile);  
	return 1;
  }
  LInt holeno = in1;
  if (VERBAL > 1){
	cout << endl << "reading " << in1 << " holes..."; cout.flush();
  }
  allocateHoles(noOfHoles + holeno);
  for (LInt i = 0; i < holeno; i++) {
	if (readFileLine() < 3) {
	  cerr << endl << "ERROR reading holes in " << name << endl;
	  fclose(iFile);  
	  return -1;
	}
	NewHole(inputs[0], inputs[1]);
  }

  // close
  fclose(iFile);
  if (VERBAL > 0){
	cout << endl << "End reading file : " << name << endl; cout.flush();
  }
  // file is read
  //----------------------------------------
  return 0;
}//readPolyFile()



//---------------------------------------------------------------------------
// reads a background grid from a poly file
// returns 0 if succesful
//---------------------------------------------------------------------------
int   
madd::readBackgroundFile(char* name)
{
  //   clearAll();
  if (VERBAL > 0){
	cout << endl << "reading background grid file : " << name << endl; cout.flush();
  }
  // open the file
  iFile = fopen(name, "r");
  if (iFile == NULL) {
	cout << endl << "ERROR IN " << name << endl;
	fclose(iFile);  
	return -1;
  }

  // get number of points
  if (readFileLine() < 1) {
	cout << endl << "ERROR IN " << name << endl;
	fclose(iFile);  return -1;
  }
  LInt pointno = in1;
  int readSize = 4;

  // allocate the background structure
  allocateBackgroundNodes(pointno);

  //----------------------------------------
  // read points
  for (LInt i = 0; i < pointno; i++) {
	if (readFileLine() < 4) {
	  cout << endl << "ERROR IN " << name << endl;
	  fclose(iFile);  
	  return -1;
	}
	background_x(i) = inputs[0];
	background_y(i) = inputs[1];
	weightOfNode(i) = (float)inputs[2];
	//	  cout << inputs[2] << endl;
  }//for

  fclose(iFile);
  // file is read
  //----------------------------------------
  if (VERBAL > 0){
	cout << endl << "End reading weights file : " << name << endl; cout.flush();
  }
  return 0;
}//readBackgroundFile()


//---------------------------------------------------------------------------
// reads a line of a poly file
// returns the number of parameters read
//---------------------------------------------------------------------------
int 
madd::readFileLine()
{
  char inBuffer[180];
  int retValue;

  // clear inputs
  for (int i=0; i < 6; i++)
	inputs[i] = 0;
  // read until the line is not a comment and not EOF
  do{
	if (fgets(inBuffer, 180, iFile) == NULL)
	  return -1;
	//	cout << inBuffer << endl;
  }while(inBuffer[0] == '#');
	
  retValue = sscanf(inBuffer, "%ld%lf%lf%lf%lf%lf%lf", &in1, inputs, inputs+1, inputs+2, 
					  inputs+3, inputs+4, inputs+5);
  // return the number of read parameters
  return retValue;
}//readFileLine()


//---------------------------------------------------------------------------*
//                                                                           *
//     END OF FILE IO ROUTINES                                               *
//                                                                           *
//---------------------------------------------------------------------------*



//---------------------------------------------------------------------------*
//                                                                           *
//     START OF THE INTERFACE UTILITIES                                      *
//                                                                           *
//---------------------------------------------------------------------------*



//---------------------------------------------------------------------------
//  computes the statistics for the quality of the decomposition
//---------------------------------------------------------------------------
void 
madd::getSepStatistics(double* anglDistribution,     // returns the angle distribution 
					   double* minAngle,             // returns the min angle
					   double* meanAngle,            // returns the mean angle
					   double* phi,                  // returns the lower bound angle 
					   long int* finalSegments,      // returns the number of final segments
					   double* minWeightedArea,      // returns the min area (according to expectation)
					   double* maxWeightedArea,      // returns the max area (according to expectation)
					   double* meanWeightedArea,     // returns the mean area (according to expectation) 
					   double* minCommunication,     // returns the min sep/area ratio (according to expectation)
					   double* maxCommunication,     // returns the max sep/area ratio (according to expectation) 
					   double* meanCommunication)    // returns the mean sep/area ratio (according to expectation)
{
  // the angle distribution
  LInt theAngleDistribution[20];
  for (int i=0; i < 19; i++)
	theAngleDistribution[i]= 0;
  double totalSumOfAngles;
  double theMinAngle = M_PI;
  LInt totalNoOfAngles = 0;

  // segments containing the same point
  LInt adjSegms[20];
  // the sorted list of angles that are formed on one point
  double theAngle[191];

  //----------------------------------------
  // for all points compute the angle statistics
  //----------------------------------------
  for (LInt thisPoint = 0; thisPoint < noOfPoints; thisPoint++){
#ifdef __DEBUG_MODE
	cout << endl << "Checking point: " << thisPoint;	cout.flush();
#endif
	int noOfadjSegms = 0;
	int noOfAngles = 0;
	int externSegms = 0;
	//----------------------------------------
	// find adjacent segments to point
	for (LInt segm = 0; segm < noOfSegms; segm++){
	  if (segm_point1(segm) == thisPoint || segm_point2(segm) == thisPoint){
		adjSegms[noOfadjSegms] = segm;
		noOfadjSegms++;
		if (segm < noOfExtSegms)
		  externSegms++;
	  }// if
	}// for
#ifdef __DEBUG_MODE
	cout << "   noOfadjSegms:" << noOfadjSegms;
	cout << "   externSegms:" << externSegms;
	cout.flush();
#endif

	if (externSegms == noOfadjSegms)
	  continue; // all segments are external

	//----------------------------------------
	// separators exist on this point
	// compute the angles
	// go through all pairs of segments and compute the angles
	// store the angles into a sorted list
	//----------------------------------------
	for (int i=0; i < noOfadjSegms - 1; i++){
#ifdef __DEBUG_MODE
	  cout << endl << " Inside first loop, i: " << i; cout.flush();
#endif
	  // get the other point in point1
	  LInt point1 = (segm_point1(adjSegms[i]) == thisPoint) ? segm_point2(adjSegms[i]): segm_point1(adjSegms[i]);
	  // check if this is a boundary segment
	  int outExtern = (adjSegms[i] < noOfExtSegms)? 1: 0;
	  // get the rest of the segments
	  for (int j=i+1; j < noOfadjSegms; j++){
		// find the number of boundary segments 
	    int isExtern = outExtern + ((adjSegms[j] < noOfExtSegms)? 1: 0);
#ifdef __DEBUG_MODE
		cout << endl << " Inside second loop, j: " << j << "  isExtern: " << isExtern; cout.flush();
#endif
		// check if they are both boundary segments  
		if (isExtern > 1)
		  continue; // these are two  boundary segments, ignore them
		// get the other point in point1
		LInt point2 = (segm_point1(adjSegms[j]) == thisPoint) ? segm_point2(adjSegms[j]): segm_point1(adjSegms[j]); 

#ifdef __DEBUG_MODE
		cout << endl << " Computing angle of " << point1 << ", " << thisPoint << ", " << point2 << endl; cout.flush();
#endif
		// compute the angle
		double thisAngle =  angle(point(point1), point(thisPoint), point(point2));

		//----------------------------------------
		// if the two segments are not colinear
		// or one is boundary segment insert the angle
		if (thisAngle < PI_MINUS || isExtern > 0){ // they are not on the same line or one is external
#ifdef __DEBUG_MODE
		  cout << endl << " Adding angle " << thisAngle << "..." << endl; cout.flush();
#endif
		  // add the angle to the sorted list
		  int a;
		  // find the position
		  for (a=0; a < noOfAngles && thisAngle > theAngle[a]; a++);
		  // shift the list
		  for (int b = noOfAngles-1; b >= a; b--){
			theAngle[b+1] = theAngle[b];
		  }//for
		  // insert angle
		  theAngle[a] = thisAngle;
		  noOfAngles++;
		  //  angle inserted
		  //----------------------------------------
		}// if 

	  }// for all next segments
	}//for all adjacent segments

	//----------------------------------------
	// we have a sorted list of angles
	// check the fisrt noOfadjSegms-1 angles
	//----------------------------------------
#ifdef __DEBUG_MODE
	cout << "   noOfAngles:" << noOfAngles;
	cout.flush();
#endif
	if (noOfAngles < 1)
	  continue; // no angles (they are on the same line), ignore them

	//----------------------------------------
	// get the sum, min of the angles
	double thisSum = 0;
	int activeAngles = noOfadjSegms - 1;
#ifdef __DEBUG_MODE
	cout << "   activeAngles:" << activeAngles;
	cout.flush();
#endif

	for (int a=0; a < activeAngles; a++){
	  // get the sum
	  thisSum += theAngle[a];
	  int angleGroup = (int)((theAngle[a] * 180) / (10 * M_PI));
	  if (angleGroup >= 0 && angleGroup <= 19){
		// get the distribution sum
		theAngleDistribution[angleGroup]++;
	  }//if
	}//for

	if (externSegms < 1 && fabs(thisSum - M_PI) > 0.01){
	  // all angles count, a total of 2xPI
	  // add the remaining angle
	  int angleGroup = (int)(((TWO_PI - thisSum) * 180) / (10 * M_PI));
	  if (angleGroup >= 0 && angleGroup <= 19){
		// get the distribution sum
		theAngleDistribution[angleGroup]++;
	  }
	  activeAngles++;
	  // the sum is 2*pi
	  thisSum = TWO_PI;
	}
	// add to the statistics
	totalNoOfAngles += activeAngles;
	totalSumOfAngles += thisSum;
	// get the min angle
	if (theMinAngle > theAngle[0])
	  theMinAngle = theAngle[0];
	//----------------------------------------
	// done
	//----------------------------------------
#ifdef __DEBUG_MODE
	cout << "   activeAngles:" << activeAngles;
	cout << "   thisSum:" << thisSum;
	cout << "   theMinAngle:" << theMinAngle;
	cout << "   totalNoOfAngles:" << totalNoOfAngles;
	cout << "   totalSumOfAngles:" << totalSumOfAngles;
	cout.flush();
#endif
  }// for all points

  //----------------------------------------
  // angle statistics have been compouted
  // compute balance and communication estimators
  //----------------------------------------
  double minWArea, maxWArea, meanWArea; //balance estimators
  if (BACKGROUND_AREAS || areaFunction != NULL)
	minWArea = maxWArea = meanWArea = subdomainArea[0] / (subdomainMinArea[0] * DEC_AREA_RATIO);
  else
	minWArea = maxWArea = meanWArea = subdomainArea[0] / (double)(subdomainWeight[0] + 1.0);
  double minSArea, maxSArea, meanSArea; //communication estimators
  minSArea = maxSArea = meanSArea = separatorLengthOfSubdomain(0) / sqrt(subdomainArea[0]);
  double totalDomainArea=0;
  for (int i=0; i < noOfSubdomains; i++){
	totalDomainArea += subdomainArea[i];
  }
  //----------------------------------------
  // go throug all the subdomains and
  // compute balance and communication estimators
  //----------------------------------------
  for (int i=1; i < noOfSubdomains; i++){
	
	//----------------------------------------
	// compute balance estimator
	double thisWArea;
	if (BACKGROUND_AREAS || areaFunction != NULL)
	  thisWArea = subdomainArea[i] / (subdomainMinArea[i] * DEC_AREA_RATIO);
	else
	  thisWArea = subdomainArea[i] / totalDomainArea;
	if (thisWArea > maxWArea)
	  maxWArea = thisWArea;
	if (thisWArea < minWArea)
	  minWArea = thisWArea;
	meanWArea += thisWArea;

	//----------------------------------------
	// compute communication estimator
	double thisSArea = separatorLengthOfSubdomain(i) / sqrt(subdomainArea[i]);
	if (thisSArea > maxSArea)
	  maxSArea = thisSArea;
	if (thisSArea < minSArea)
	  minSArea = thisSArea;
	meanSArea += thisSArea;
  } // for all subdomains
  //----------------------------------------
  //  balance and communication estimators have been computed
  //----------------------------------------
 
  meanWArea /= noOfSubdomains;
  meanSArea /= noOfSubdomains;

  //----------------------------------------
  // fill the return values 
  //----------------------------------------
  *minAngle = (theMinAngle * 180) / M_PI;
  *meanAngle =  ((totalSumOfAngles / totalNoOfAngles) * 180) / M_PI;
  *phi = (PHI * 180) / M_PI;
  *finalSegments = noOfSegms;

  *minWeightedArea = minWArea;
  *maxWeightedArea = maxWArea;
  *meanWeightedArea = meanWArea; 
  *minCommunication = minSArea;
  *maxCommunication = maxSArea;
  *meanCommunication = meanSArea;

  totalNoOfAngles = 0;
  for (int i=0; i < 19; i++){
	totalNoOfAngles += theAngleDistribution[i];
  }
  for (int i=0; i < 19; i++){
	//	cout << i << ". " << theAngleDistribution[i] << endl;
	anglDistribution[i] = (double)theAngleDistribution[i] / totalNoOfAngles;
  }

  //----------------------------------------
  // statistics done
  //----------------------------------------
  return;
  
#undef __DEBUG_MODE
}



//---------------------------------------------------------------------------
//  find and write a partitioning straight line for each subdomain
//---------------------------------------------------------------------------
int   
madd::writeSubpartitions(char* name)
{
  // open the output file
  ofstream oFile;
  oFile.open(name);
  if (!oFile.good()) {
	cerr << endl << "ERROR IN " << name << endl;
	return -1;
  }
  //oFile.precision(8);
  //----------------------------------------
  // write the partitions of the subdomains
  oFile << noOfSubdomains << endl;
  for (int i=0; i < noOfSubdomains; i++){
	// write the subdomain we work on
	oFile << i+1;  
	// init the values
	initSubdomain(i);
	LInt startPoint, segmStart, segmEnd, p1, p2;
	startPoint = segmStart = p1 = segm_point1(workSegm[0]);
	segmEnd = p2 =  segm_point2(workSegm[0]);
	DFloat vector_x, vector_y, ref_x, ref_y;
	DFloat ref2_x, ref2_y;
	vector_x = vector_y = 0;
	ref_x = ref_y = 0;
	DFloat totalLength = 0;	
	
	//----------------------------------------
	// go through the segments in oriented way
	do{
	  // get the (almost) co-linaer segments
	  do{
		segmEnd = p2;
		LInt p = nextOfPoints(p1, p2);
		p1 = p2;
		p2 = p;		
	  }while((p1 != startPoint) &&
			 angle(point(segmStart), point(segmEnd), point(p2)) > (M_PI - 0.2));
	  //compute length, middle and direction vector 
	  DFloat length = distance(point(segmStart), point(segmEnd));
	  vector_x += (point_x(segmEnd) - point_x(segmStart)) * sqr(length) * length;
	  vector_y += (point_y(segmEnd) - point_y(segmStart)) * sqr(length) * length;
	  ref_x += ((point_x(segmEnd) + point_x(segmStart))) * sqrt(length);
	  ref_y += ((point_y(segmEnd) + point_y(segmStart))) * sqrt(length);
	  totalLength += sqrt(length);

	  segmStart = p1;
	}while(p1 != startPoint);
    //----------------------------------------

	ref_x /= (totalLength * 2);
	ref_y /= (totalLength * 2);
	//vector_x /= totalLength;
	//vector_y /= totalLength;
	DFloat vector2_x = vector_y;
	DFloat vector2_y = -vector_x;

	// write the first line that goes through ref in the direction of vector
	// as ax+by+c=0 line 
	oFile << "   " << vector_y;  // write the a
	oFile << " " << -vector_x; // write the b
	oFile << " " << (ref_y * vector_x - ref_x * vector_y); // write the c
	// write the second line that goes through ref in the orthogonal direction of vector
	// as ax+by+c=0 line 
	oFile << "   " << vector2_y;  // write the a
	oFile << " " << -vector2_x; // write the b
	oFile << " " << (ref_y * vector2_x - ref_x * vector2_y) << endl; // write the c

#ifdef __DEBUG_MODE
	// insert the segment into the geometry to see how it looks like
	p1 = NewPoint(ref_x, ref_y);
	p2 = NewPoint(ref_x+ 0.0005 * vector_x, ref_y + 0.0005 * vector_y);
	NewSegm(p1, p2);
	p2 = NewPoint(ref_x+ 0.0005 * vector2_x, ref_y + 0.0005 * vector2_y);
	NewSegm(p1, p2);
#endif
#undef __DEBUG_MODE
  }// for all subdomains


  //----------------------------------------
  // check if the file is good and close
  if (!oFile.good()) {
	cerr << endl << "ERROR IN " << name << endl;
	return -1;
  }
  oFile.flush();
  oFile.close();
  // file has been created
  //----------------------------------------

  if (VERBAL > 0){
	cout << endl << "End writting file : " << name << endl; cout.flush();
  }

  // file write done
  //----------------------------------------
  return 0;
}//writeSubpartitions



//---------------------------------------------------------------------------
//  interface for reading a weighted background grid
//---------------------------------------------------------------------------
int   
madd::readBackgroundWeights(char* name)
{
  BACKGROUND_WEIGHTS = true;
  BACKGROUND_AREAS = false;
  return readBackgroundFile(name);
}


//---------------------------------------------------------------------------
//  interface for reading an area background grid
//---------------------------------------------------------------------------
int   
madd::readBackgroundAreas(char* name)
{
  BACKGROUND_WEIGHTS = false;
  BACKGROUND_AREAS = true;
  return readBackgroundFile(name);
}


//---------------------------------------------------------------------------
//  SET PARAMETER FUNCTIONS
//  THEY ALL RETURN 0 IF SUCCESFULL
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
//  sets the weight function
//  the function is loaded dynamically
//---------------------------------------------------------------------------
int
madd::setWeightFunction(char *libraryName,   // the library name
						char *functionName)  // the function
{
  // get the library
  libraryHandle = dlopen(libraryName, RTLD_LAZY);
  if (libraryHandle == NULL){
	cout << endl << "Error opening library " << libraryName << ".  " << dlerror() << endl;
	cout.flush();
	return 1;
  }
     
  // get the function
  weightFunction = (double (*)(double, double))dlsym(libraryHandle, functionName);
  char* errorMsg = dlerror();;
  if (weightFunction == NULL || errorMsg != NULL){
	cout << endl << "Error locating function " << functionName <<".  " << errorMsg << endl;
	cout.flush();
  	return 2;	
  }

  // this is the function to evaluate on the grid
  gridFunction = weightFunction;
  //  cout << endl << "weightFunction(0,0) = " << (*weightFunction)(0.0, 0.0);
  
  return 0;
} 


//---------------------------------------------------------------------------
//  sets the area function
//  the function is loaded dynamically
//---------------------------------------------------------------------------
int
madd::setAreaFunction(char *libraryName, char *functionName)
{
  // get the library
  libraryHandle = dlopen(libraryName, RTLD_LAZY);
  if (libraryHandle == NULL){
	cout << endl << "Error opening library " << libraryName << ".  " << dlerror() << endl;
	cout.flush();
	return 1;
  }
     
  // get the function
  areaFunction = (double (*)(double, double))dlsym(libraryHandle, functionName);
  char* errorMsg = dlerror();;
  if (areaFunction == NULL || errorMsg != NULL){
	cout << endl << "Error locating function " << functionName <<".  " << errorMsg << endl;
	cout.flush();
  	return 2;	
  }

  // this is the function to evaluate on the grid
  gridFunction = areaFunction; 
  // cout << endl << "areaFunction(0,0) = " << (*areaFunction)(0.0, 0.0);
  // cout << endl << "areaFunction(1,1) = " << (*areaFunction)(1.0, 1.0) << endl;
  
  return 0;
} 


//---------------------------------------------------------------------------
//  sets minimum acceptable angle formed during the decomposition
//---------------------------------------------------------------------------
int
madd::setPhi(double phiValue)
{

  
  if (phiValue > 0.4)
	return 1;
  if (phiValue < 0.2)
	return -1;

  PHI = M_PI * phiValue;

  return 0;
}


//---------------------------------------------------------------------------
//  sets the min area ratio for decomposing a subdomain
//  when an area function is given or an area background grid is given 
//  the subdomain to be decomposed must have ratio area/excpected_area
//  greater than this ratio
//  if the min area ratio is 1, then all the created subdomain will have
//  area <= excpected_area 
//---------------------------------------------------------------------------
int
madd::setDecAreaRatio(double ratioValue)
{
  if (ratioValue <= 0)
	return -1;

  DEC_AREA_RATIO = ratioValue;

  return 0;
}


//---------------------------------------------------------------------------
//  sets the step length for the uniform background grid
//  this grid will be created automatically on the fly 
//---------------------------------------------------------------------------
int
madd::setGridStep(double gridValue)
{
 
  if (gridValue <= 0)
	return -1;

  GRID_STEP = gridValue;

  return 0;
}


//---------------------------------------------------------------------------
//  sets the maximum number of the subdomains
//  this is for the cases we use an area criterio
//  if the number of created subdomains reaches this number, 
//  the decomposition will stop
//---------------------------------------------------------------------------
int
madd::setMaxNoOfSubdomaind(int maxValue)
{
 
  if (maxValue < 0)
	return -1;
 
  MAX_NO_OF_SUBDOMAINS = maxValue;

  return 0;
}


//---------------------------------------------------------------------------
//  sets the verbal level
//  0: no messages
//  1: basic messages
//  2: some details
//  > 2: more details
//---------------------------------------------------------------------------
int
madd::setVerbal(int verbValue)
{
 
  if (verbValue < 0)
	return -1;
 
  VERBAL = verbValue;

  return 0;
}


//---------------------------------------------------------------------------
//  sets the uniform refining and decomposition factor
//  the larger the parameter the more refining we will obtain
//  and, in the cases of graded decompositions, a more smooth decomposition
//
//  By default the value is 2 for graded decompositions, 3 for unioform
//  Use a value of 2.1-4, if more refinement is required
//---------------------------------------------------------------------------
int
madd::setUniformRefineLevel(double uniformValue)
{

  if (uniformValue < 0.0)
	return -1;
  if (uniformValue > 10.0)
	return 1;

  AREA_WEIGHT_COEF = uniformValue;

  return 0;
}


//---------------------------------------------------------------------------
//  sets the graded refining and decomposition factor
//  the larger the parameter, the more gradation we will obtain
//
//  By default the value is 0. 
//  In most cases a value = uniformValue gives good results
//---------------------------------------------------------------------------
int
madd::setAdaptiveRefineLevel(double adaptiveValue)
{

  if (adaptiveValue < 0.0)
	return -1;
  if (adaptiveValue > 10.0)
	return 1;

  POINT_WEIGHT_COEF = adaptiveValue;

  return 0;
}


//---------------------------------------------------------------------------
//  sets the maximum acceptble area imbalance during the madd partition procedure
//  default is 0.7
//  a value 0.7 - 0.8 works fine.
//---------------------------------------------------------------------------
int
madd::setMaxImbalanceLevel(double imbalanceValue)
{

  if (imbalanceValue < 0.6)
	return -1;
  if (imbalanceValue > 0.9)
	return 1;

  MAX_IMBALANCE_LIMIT = imbalanceValue;

  return 0;
}


//---------------------------------------------------------------------------
// the madd partiton procedure smoothens the separators repeatedly
// this value defines the number of the smoothing iterations
// default is = 4
// a zero value will result no smoothing during the insertion of 
// the separators
//---------------------------------------------------------------------------
int
madd::setMaxSmoothLevel(int smoothLevel)
{
  if (smoothLevel > 10)
	return 1;
  if (smoothLevel < 0)
	return -1;

  MAX_SMOOTH_DEPTH = smoothLevel;

  return 0;
}


//---------------------------------------------------------------------------
// the madd partiton procedure smoothens the separators by examining
// the adjacent points left and right of the current separator
// this function sets the number of points, left and right, that
// the smoothing procedure will examine.
// For example, a value of 3 wiil cause the algorithm to examine 
// 3 points left and 3 points right
// Default is 5
//---------------------------------------------------------------------------
int
madd::setSmoothVertices(int smoothVertices)
{
  if (smoothVertices > 10)
	return 1;
  if (smoothVertices < 0)
	return -1;

  SMOOTHING_POINTS = smoothVertices * 2 + 1;

  return 0;
}


//---------------------------------------------------------------------------
// sets the interpolation coefficient for computing the weigts of the inserted
// pointe from the weights of the existing segment end points.
// a Value of 0 will result no interpolation
// a value of 1 will cause a (weighted) linear interpolation
// Default: 0
//---------------------------------------------------------------------------
int
madd::setInterpolationWeightLevel(double weightValue)
{

  if (weightValue < 0.0)
	return -1;
  if (weightValue > 1.0)
	return 1;

  INTERPOLATION_WEIGHT_COEF = weightValue;

  return 0;
}


//---------------------------------------------------------------------------
// sets the maximum weight of the Metis graph nodes
// Default : 1000
//---------------------------------------------------------------------------
int
madd::setMetiseMaxNode(int weightValue)
{

  if (weightValue < 10)
	return -1;
  if (weightValue > 10000)
	return 1;

  MAX_MADD_NODE_WEIGHT = weightValue;

  return 0;
}


//---------------------------------------------------------------------------
// sets the maximum weight of the Metis graph edges
// Default : 1600
//---------------------------------------------------------------------------
int
madd::setMetiseMaxEdge(int weightValue)
{

  if (weightValue < 10)
	return -1;
  if (weightValue > 10000)
	return 1;

  MAX_MADD_EDGE_WEIGHT = weightValue;

  return 0;
}


//---------------------------------------------------------------------------
//  sets the weight of all point to the given value
//---------------------------------------------------------------------------
void
madd::weightAllPoints(int weight)
{
  for (LInt i=0; i < noOfPoints; i++){
	pointWeight[i] = weight;
  }

}


//---------------------------------------------------------------------------
//  sets the weights of one subdomain to the given value
//---------------------------------------------------------------------------
void
madd::weightSubdomainPoints(int weight, int theSubDomain)
{
  initSubdomain(theSubDomain);
  LInt segm, p1, p2;
  for (LInt i = 0; i < workNoOfSegms; i++){
	segm =  workSegm[i];
	p1 = segm_point1(segm);
	p2 = segm_point2(segm);
	pointWeight[p1] = pointWeight[p2] = weight;  
  }
}



//---------------------------------------------------------------------------*
//                                                                           *
//     END OF THE INTERFACE UTILITIES                                        *
//                                                                           *
//---------------------------------------------------------------------------*

