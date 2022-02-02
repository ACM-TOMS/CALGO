


#ifndef __DEF_MADD
#define __DEF_MADD

//---------------------------------------------------------------------------
// define the data types
typedef double DFloat; 
typedef int LInt;
typedef unsigned int UInt;
typedef unsigned int uint;
typedef int LLint; 
typedef char uchar;
//---------------------------------------------------------------------------

class madd {
public:                    

//---------------------------------------------------------------------------
// the constructor
  madd();
//---------------------------------------------------------------------------
// the destructor
  ~madd();
//---------------------------------------------------------------------------
// the main decomposition routine
// input the number of subdomains to be created
// if subDomains <= 0 use the area criterion
// returns the number of the created subdomains
  int decompose(int subDomains);
//---------------------------------------------------------------------------
//  computes the statistics for the quality of the decomposition
 void getSepStatistics(double* anglDistribution,     // returns the angle distribution 
					   double* minAngle,             // returns the min angle
					   double* meanAngle,            // returns the mean angle
					   double* phi,                  // returns the lower bound angle 
					   long int* finalSegments,      // returns the number of final segments
					   double* minWeightedArea,      // returns the min area (according to expectation)
					   double* maxWeightedArea,      // returns the max area (according to expectation)
					   double* meanWeightedArea,     // returns the mean area (according to expectation) 
					   double* minCommunication,     // returns the min sep/area ratio (according to expectation)
					   double* maxCommunication,     // returns the max sep/area ratio (according to expectation) 
					   double* meanCommunication);    // returns the mean sep/area ratio (according to expectation)
//---------------------------------------------------------------------------
//  interface for reading a weighted background grid from a file
//  returns 0 if succesful
 int readBackgroundWeights(char* name);
//---------------------------------------------------------------------------
//  interface for reading an area background grid from a file
//  returns 0 if succesful
 int readBackgroundAreas(char* name);
//---------------------------------------------------------------------------
// writes the domain and all the subdomains to poly files
// returns 0 if succesful
  int  writePolyFileAll(char* fileName);
//---------------------------------------------------------------------------
// writes the whole domain or one subdomain to a poly file
// returns 0 if succesful
  int writePolyFile(char* name,               // the name of the file
					int subdomain,            // the subdomain, if >= 0 one subdomain
					                          //                if = -1 the decomposition
                                              //                if < -1 the initial domain
					char includeTriangles,    // 'y' = include the triangulation
					char includeAllHoles,     // 'y' = include all holes
					char includeAllPoints,    // 'y' = include all points
					int hole) ;                // if >= 0, write only this hole
//---------------------------------------------------------------------------
// writes the decomposition in points, segments, subdomains -segments
// returns 0 if succesful
//---------------------------------------------------------------------------
   int  writeSubdomainFile(char* name);        // the name of the file
//---------------------------------------------------------------------------
// reads a poly file
// returns 0 if succesful
  int readPolyFile( char* name);
//---------------------------------------------------------------------------
//  find and write a partitioning straight line for each subdomain
//---------------------------------------------------------------------------
  int  writeSubpartitions(char* name);
//---------------------------------------------------------------------------
//  SET PARAMETER FUNCTIONS
//  THEY ALL RETURN 0 IF SUCCESFULL
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  sets the weight function
//  the function is loaded dynamically
  int setWeightFunction(char *libraryName,   // the library name
						char *functionName); // the function
//---------------------------------------------------------------------------
//  sets the area function
//  the function is loaded dynamically
  int setAreaFunction(char *libraryName,     // the library name
					  char *functionName);   // the function
//---------------------------------------------------------------------------
//  sets minimum acceptable angle formed during the decomposition
  int setPhi(double phiValue);
//---------------------------------------------------------------------------
//  sets the min area ratio for decomposing a subdomain
//  when an area function is given or an area background gri is given 
//  the subdomain to be decomposed must have ratio area/excpected_area
//  greater than this ratio
//  if the min area ratio is 1, then all the created subdomain will have
//  area <= excpected_area 
  int setDecAreaRatio(double ratioValue);
//---------------------------------------------------------------------------
//  sets the step length for the uniform background grid
//  this grid will be created automatically on the fly 
  int setGridStep(double gridValue);
//---------------------------------------------------------------------------
//  sets the maximum number of the subdomains
//  this is for the cases we use an area criterio
//  if the number of created subdomains reaches this number, 
//  the decomposition will stop
  int setMaxNoOfSubdomaind(int maxValue);
//---------------------------------------------------------------------------
//  sets the verbal level
//  0: no messages, 1: basic messages, 2: some details,  > 2: more details
  int setVerbal(int verbValue);
//---------------------------------------------------------------------------
//  sets the uniform refining and decomposition factor
//  the larger the parameter the more refining we will obtain
//  and, in the cases of graded decompositions, a more smooth decomposition
//
//  By default the value is 2 for graded decompositions, 3 for uniform
//  Use a value of 2.1-4, if more refinement is required
  int setUniformRefineLevel(double uniformValue);
//---------------------------------------------------------------------------
//  sets the graded refining and decomposition factor
//  the larger the parameter, the more gradation we will obtain
//
//  By default the value is 0. 
//  In most cases a value = uniformValue gives good results
  int setAdaptiveRefineLevel(double adaptiveValue);
//---------------------------------------------------------------------------
//  sets the maximum acceptble area imbalance during the madd partition procedure
//  default is 0.7
//  a value 0.7 - 0.8 works fine.
  int setMaxImbalanceLevel(double imbalanceValue);
//---------------------------------------------------------------------------
// the madd partiton procedure smoothens the separators repeatedly
// this value defines the number of the smoothing iterations
// default is = 4
// a zero value will result no smoothing during the insertion of 
// the separators
  int setMaxSmoothLevel(int smoothLevel);
//---------------------------------------------------------------------------
// the madd partiton procedure smoothens the separators by examining
// the adjacent points left and right of the current separator
// this function sets the number of points, left and right, that
// the smoothing procedure will examine.
// For example, a value of 3 wiil cause the algorithm to examine 
// 3 points left and 3 points right
// Default is 5
  int setSmoothVertices(int smoothVertices);
//---------------------------------------------------------------------------
// sets the interpolation coefficient for computing the weigts of the inserted
// pointe from the weights of the existing segment end points.
// a Value of 0 will result no interpolation
// a value of 1 will cause a (weighted) linear interpolation
// Default: 0
  int setInterpolationWeightLevel(double weightValue);
//---------------------------------------------------------------------------
// sets the maximum weight of the Metis graph nodes
// Default : 1000
  int setMetiseMaxNode(int weightValue);
//---------------------------------------------------------------------------
// sets the maximum weight of the Metis graph edges
// Default : 1600
  int setMetiseMaxEdge(int weightValue);
//---------------------------------------------------------------------------
//  sets the weight of all point to the given value
  void weightAllPoints(int weight);
//---------------------------------------------------------------------------
//  sets the weights of one subdomain to the given value
  void weightSubdomainPoints(int weight, int theSubDomain);

  int cylinder(); // temporary here


private:
//---------------------------------------------------------------------------
//  private variables
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// the weight fuctions
   void *libraryHandle;
   double (*weightFunction)(double x, double y);
   double (*areaFunction)(double x, double y);
   double (*gridFunction)(double x, double y);
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// the basic pslg sturcture and variables
   LInt noOfPoints, noOfSegms;
   LInt noOfExtPoints, noOfExtSegms;
   LInt THE_SIZE, TRIA_SIZE;
   LInt POINTLIST_SIZE, SEGMLIST_SIZE;
   LInt HOLELIST_SIZE ;
   LInt noOfHoles, workNoOfHoles;
   DFloat *thePointsList;
   uchar  *basicPoint;
   float  *pointWeight;
   LInt  *theSegmsList;
   DFloat *theHolesList;
   DFloat totalBoundaryLength;
   DFloat domainMinX, domainMaxX, domainMinY, domainMaxY; 
   float maxPointWeight;
    DFloat refineLength;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// subdomain variables
   uchar** subdomainHoles;
   uchar* parentSubdomainHoles;
   int noOfSubdomains, totalNoOfSubdomains, workSubdomain;
   LInt* noOfSubSegms;
   LInt *toWorkSegm1, *toWorkSegm2;
   LInt** subdomainSegments;
   uchar *segmColourList;
   DFloat *subdomainArea;
   DFloat *subdomainWeight;
   DFloat *subdomainMinArea;
   DFloat *subdomainLength;
   unsigned int *subdomainBackgroundStart;
   unsigned int *subdomainBackgroundEnd;

//---------------------------------------------------------------------------
// the work structure
   DFloat *workHolesList;
   LInt* workSegm;
   LInt workNoOfSegms;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// the background weight list (unstructured)
   DFloat *backgroundNodeList; 
   float *backgroundAreaList;
   unsigned int backgroundListSize;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// the triangulation structure
   //DFloat *theTriaArea;
   LInt noOfTria;
   LInt *triaNeighborList;
   DFloat *theTriaPointList;
   LInt  *theTriaList;
   LInt noOfTriangles, trianglePoints;
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------
// madd variables
   LInt madSepStart, madSepEnd;
   DFloat partArea[16];
   LInt pointsPerColour[17];
   LInt edgesPerColour[17];
   DFloat partLength[16];

  
//---------------------------------------------------------------------------*
//   MAIN ROUTINES                                                           *
//---------------------------------------------------------------------------*
//---------------------------------------------------------------------------
// the main partiton routine
// it partitions theSubDomain into two subdomains
// returns 0 if successful
  int partition(int theSubDomain);
//---------------------------------------------------------------------------*
// the Delaunay triangulation routine
// triangulates the work subdomain using Shewchuck's Triangle
// gets a string of flags passed to Triangle
  int triangulateDomain(char *inFlags);


//---------------------------------------------------------------------------*
//     GRADING CONTROL ROUTINES                                              *
//---------------------------------------------------------------------------*
//---------------------------------------------------------------------------
// evaluate the subdomain over a uniform background gird
// using the gridFunction
// The grid is creted on the fly
  int evaluateGrid(int theSubdomain,   // the subdomain to be evaluated
				   DFloat *minValue,   // return the min value
				   DFloat *maxValue,   // return the max value
				   DFloat *meanValue,  // return the mean value
				   DFloat *totalValue);// return the sum
//---------------------------------------------------------------------------
//  assigns the subdomain weights and min-areas from the background grid
//  to the two created subdomains from the parent domain
  void getBackgroundWeights(int parentSubdomain,  // thye parent domain, also the first subdomain
						   int subdomain2) ;      // the second subdomain


//---------------------------------------------------------------------------*
//   SUBDOMAIN ROUTINES                                                      *
//---------------------------------------------------------------------------*
//---------------------------------------------------------------------------
// initialize theSubDomain as our work subdomain
  void initSubdomain(int theSubDomain);
//---------------------------------------------------------------------------
// make theSubDomain the workSubdomain
  void assignToWorkSubdomain(int theSubDomain);
//---------------------------------------------------------------------------
// create pointers from pointS to the work segments
// (used in identifying adjacency)
  void createPointToWorkSegm();
//---------------------------------------------------------------------------
// put into the workHolesList the holes of theSubDomain
  void createWorkHoles(int theSubDomain);
//---------------------------------------------------------------------------
// create the initial subdomain 0
  void createInitialSubdomain();
//---------------------------------------------------------------------------
// create  a new set of points and segments from the work segments (of some subdomain)
  LInt createPointsSegms(DFloat **newPointAddress, // to the new point list
						LInt **newSegmAddress,    // to the new segment list
						LInt **fromPointAddress);  // to the reference back to the original points


//---------------------------------------------------------------------------*
//    THE MADD ROUTINES                                                      *
//---------------------------------------------------------------------------*
//---------------------------------------------------------------------------
//  allocate the initial graph nodes
  void allocateGraphNodes();
//---------------------------------------------------------------------------
//  allocate the initial graph edges
  void allocateGraphEdges();
//---------------------------------------------------------------------------
//  reallocate the initial graph edges
  void reallocateGraphEdges();
//---------------------------------------------------------------------------
//  returns the end node pointed from the startNode
//  it collapses the route for efficiency
  LLint getEndNode(LLint startNode);
//---------------------------------------------------------------------------
//  returns the node in linkTriangle that is adjacent to toTriangle
  LLint getLinkNode(LInt linkTriangle, LInt toTriangle);
//---------------------------------------------------------------------------
//  collapses the node fromNode to toNode
  void collapseNodeTo(LLint fromNode, LLint toNode);
//---------------------------------------------------------------------------
// insert an edge between node1, node2, with weight ...weight
  void insertGraphEdge(LLint node1, LLint node2, DFloat weight); 
//---------------------------------------------------------------------------
// the main partiton routine
// returns 0 if succesful
  int MADDpartition();
//---------------------------------------------------------------------------
// creates the MADD graph for the Delaunay triangulation
// returns 0 if succesful
  int createMADDGraph();
//-------------------------------------------------------------------------------
// creates the Metis graph structure from the MADD graph
// returns 0 if succesful
  int createMetisStructure();
//---------------------------------------------------------------------------
// partition the MEtis graph 
// returns true if succesful
  bool partitionGraph();
//---------------------------------------------------------------------------
// Check the components connectivity for the partitioned  MEtis graph 
// If we have > 2 components glue one of them and return the remaining 
// number of components
// By invoking repeatedly we will get two connencted components.
  int checkConnectivity();
//-------------------------------------------------------------------------------
// inserts the separators in the geometry
// returns true if succesful
  bool insertSeparators();
//-------------------------------------------------------------------------------
// insert an inner separator to the geometry
  void insertInnerSeparator(LInt inp1,        // the first boundary point
						   LInt inp2,        // the second boundary point
						   DFloat *center,   // the center
						   int sepcolor1);    // the color of the included node (1 - 2 nodes is the partition)
//-------------------------------------------------------------------------------
// insert the outer separator to the geometry
  void insertEdgeSep(LInt p1,         // the first boundary point 
					LInt p2);         // the second boundary point
//-------------------------------------------------------------------------------
// insert the optimal outer separator to the geometry
  void insertOptEdgeSep(LInt p1,         // the first boundary point
					    LInt p2,         // the second boundary point
					   int depth);       // the current smoothing depth
//-------------------------------------------------------------------------------
// find the optimal inner separator 
// returns the best coefficient 
// and the new separator in the parameters
  DFloat findOptSep(LInt *op1,                // the first boundary point, initial, optimal
			  	   LInt *op2,                // the second boundary point, initial, optimal
				   DFloat *center,           // the center, initial, optimal
				   DFloat *newArea,          // the new areas of the subdomains
				   LInt *changeEdgesColour);  // the list of the edges to change color
//-------------------------------------------------------------------------------
// find the optimal outer separator 
// returns the best coefficient 
// and the new separator in the parameters
  DFloat findOptSep(LInt *op1,                // the first boundary point, initial, optimal
				   LInt *op2,                 // the second boundary point, initial, optimal
				   DFloat *newArea,           // the new areas of the subdomains
				   LInt *changeEdgesColour);  // the list of the edges to change color
//-------------------------------------------------------------------------------
// check if the new area is acceptably balanced
  bool isBalanced(DFloat *newArea);
//-------------------------------------------------------------------------------
// compute the basic coefficient ((used in the MADD graph generation)
// for an inner partial separator
  DFloat basicCoef1(DFloat* center, // the circumcenter
				   LInt p1);        // the boundary point
//-------------------------------------------------------------------------------
// compute the basic coefficient ((used in the MADD graph generation)
// for an outer separator
  DFloat basicCoef2(LInt p1,  // the first  boundary point 
				   LInt p2);  // the second  boundary point 

//-------------------------------------------------------------------------------
// compute the smoothing coefficient ((used in the insertion of separators)
// for an inner separator
  DFloat smoothingCoef1(LInt p1,          // the first  boundary point 
					   DFloat* center,    // the circumcenter
					   LInt p2);          // the second boundary point
//-------------------------------------------------------------------------------
// compute the smoothing coefficient ((used in the insertion of separators)
// for an outer separator
  DFloat smoothingCoef2(LInt p1, LInt p2);
//---------------------------------------------------------------------------
// returns true if the point is used by another separator
  bool separatorPointUsed(LInt thePoint);



//---------------------------------------------------------------------------*
//     FILE IO ROUTINES                                                      *
//---------------------------------------------------------------------------*
//---------------------------------------------------------------------------
// reads a background grid from a poly file
// returns 0 if succesful
  int readBackgroundFile(char* name);
//---------------------------------------------------------------------------
// reads a line of a poly file
// returns the number of parameters read
//---------------------------------------------------------------------------
  int readFileLine();
  

//---------------------------------------------------------------------------*
//     MEMORY ALLOCATION & CLEAR ROUTINES                                    *
//---------------------------------------------------------------------------*
//---------------------------------------------------------------------------
// clear the partition structures
  void  clearPartition();
//---------------------------------------------------------------------------
// clear all structures
  void clearAll();
//---------------------------------------------------------------------------
// reset all variables
  void setAllZero();
//---------------------------------------------------------------------------
// clear the trianguleation
  void clearTriangleList();
//---------------------------------------------------------------------------
// allocate the subdomains structures
  void allocateSubdomainStr();
//---------------------------------------------------------------------------
// allocate the point structures
  void allocatePoints(LInt size);
//---------------------------------------------------------------------------
// allocate the background grid structure
  void allocateBackgroundNodes(unsigned int theSize);
//---------------------------------------------------------------------------
// allocate the holes structure
  void allocateHoles(LInt size);
//---------------------------------------------------------------------------
// allocate the segments structure
  void allocateSegms(LInt size);
//---------------------------------------------------------------------------
// print an error message and exit
  void  ErrorExit(char message[]);
//---------------------------------------------------------------------------
// error message Insufficient memory, exit
  void InsufficientMem(char where[]);



//---------------------------------------------------------------------------*
//      GEOMETRY ROUTINES                                                    *
//---------------------------------------------------------------------------*
//---------------------------------------------------------------------------
// check if a partial inner separator intersects (internally) the boundary
bool intersectsBoundary(DFloat *center,  // the circumcenter
						 LInt checkPoint); // the boundary point
//---------------------------------------------------------------------------
// check if a partial outer separator intersects (internally) the boundary
bool intersectsBoundary(LInt checkPoint1,  // the first boundary point
						 LInt checkPoint2);  // the second boundary point
//---------------------------------------------------------------------------
//  return the work segemnt defined by two points
//---------------------------------------------------------------------------
  LInt findWorkSegm(LInt p1,    // the first point 
				   LInt p2);    // the second point 
//---------------------------------------------------------------------------
//  return the point "left" of a point
 LInt leftOfPoint(LInt point);
//---------------------------------------------------------------------------
//  return the point "right" of a point
 LInt rightOfPoint(LInt point);
//---------------------------------------------------------------------------
//  return the next point of the sequence point1, point2
  LInt nextOfPoints(LInt point1,   // the first point 
				   LInt point2);   // the second point 
//---------------------------------------------------------------------------
//  return the adjacent numberOfPoints-1 points of a point
  void nextPointsOfPoint(LInt point,          // the point for which we want the adjacent
				 		 LInt *nextPoint,     // a list of the ajacent points in order
                                              // in the middle is the given point
						int numberOfPoints);  // the total number of points
//---------------------------------------------------------------------------
//  return the segment "left" of a segment
  LInt leftOfSegm(LInt aSegm);
//---------------------------------------------------------------------------
//  return the segment "right" of a segment
  LInt rightOfSegm(LInt aSegm);
//---------------------------------------------------------------------------
//  returns true if the work subdomain contains aPoint
  bool containsPoint(DFloat *aPoint);
//---------------------------------------------------------------------------
//  assigns the holes to the subdomain 
//  the parent holes are in parentSubdomainHoles[]
  void findSubdomainHoles(LInt subdomain);
//---------------------------------------------------------------------------
//  an interface for checking if a hole belongs to the work subdomain
//  if not included return 0
//  if included return 1
  uchar checkThisHole(LInt theHole);
//---------------------------------------------------------------------------
//  check if a hole belongs to the work subdomain
//  if not included return 0
//  if included return 1
  uchar containsHole(DFloat *aPoint);
//---------------------------------------------------------------------------
//  refine all segments according to a distance dist
  void breakAllSegm(DFloat dist,  uchar basicFlag);
//---------------------------------------------------------------------------
//  refine all segments into a number of pieces
  void breakAllSegm(int pieces,  uchar basicFlag);
//---------------------------------------------------------------------------
//  refine one segment into subsegmants of size ...size
//  in weighted gradation, it takes into account the weights and the gradation factor
  void breakSegm(LInt theSegm, DFloat size,  uchar basicFlag);
//---------------------------------------------------------------------------
//  refine the inserted set of separators
 int refineThisSepLength();
//---------------------------------------------------------------------------
//  refine one segment into subsegments
//  in the weighted case a graded refinement will be produced
  void breakSegm(LInt theSegm,     // the segment to be refined
			  	 int pieces,       // the number of pieces 
				uchar basicFlag);  // if > 0, this will be assigned to the new points
                                   // else 0 will be assigned, except to the middle point
//---------------------------------------------------------------------------
//  compute the domain lengths
  void computeDomainLengths();
//---------------------------------------------------------------------------
//  return the total area of the subdomain
  DFloat areaOfSubdomain(int subdomain);
//---------------------------------------------------------------------------
//  return the total length of the separators of  subdomain
  double separatorLengthOfSubdomain(int subdomain);
//---------------------------------------------------------------------------
//  return the maximum weight of the boundary points
  float computeMaxPointWeight();
//---------------------------------------------------------------------------
// returns true if the angles from  point1 to the boundPoint are >= PHI
  bool goodExtAngle(DFloat* point1, LInt boundPoint);
//---------------------------------------------------------------------------
// returns the minimum of the four angles formed by 
// the separator between two boundary points 
  DFloat minExtAngle(LInt boundPoint1,  // the first boundary point
			 	   LInt boundPoint2);  // the second boundary point


//---------------------------------------------------------------------------*
//    BASIC GEOMETRY AND MATH ROUTINES                                       *
//---------------------------------------------------------------------------*
//---------------------------------------------------------------------------
// returns the distance between two points 
  DFloat distance(DFloat* p1, DFloat* p2);
//---------------------------------------------------------------------------
// returns the inner product of the vectors (p1-p0) * (p2 - p0) 
  DFloat innProduct(DFloat* p1, DFloat* p0, DFloat* p2);
//---------------------------------------------------------------------------
// returns the area of the triangle p1, p0, p2 
  DFloat areaOfTriangle(DFloat* p1, DFloat* p0, DFloat* p2);
//---------------------------------------------------------------------------
// returns the signed area of the triangle p1, p0, p2 
  DFloat sArea(DFloat* p1, DFloat* p0, DFloat* p2);
//---------------------------------------------------------------------------
// returns the cos of the angle p1,p0,p2 
  DFloat cos(DFloat* p1, DFloat* p0, DFloat* p2);
//---------------------------------------------------------------------------
// returns the angle formed by p1,p0,p2 
  DFloat angle(DFloat* p1, DFloat* p0, DFloat* p2);
//---------------------------------------------------------------------------
// returns the signed angle formed by p1,p0,p2 
  DFloat sangle(DFloat* p1, DFloat* p0, DFloat* p2);
//---------------------------------------------------------------------------
//  computes the circumCenter of the triangle p1, p2, p3 
void circumCenter(DFloat* p1, DFloat* p2, DFloat* p3,  // the three points of the triangle
				   DFloat* result);                   // will return the circumcenter
//---------------------------------------------------------------------------
//  computes the middle of two points p1, p2 
  void middle(DFloat* p1, DFloat* p2,    // the two points
			 DFloat* middle);            // will return middle 
//---------------------------------------------------------------------------
//  checks if two segments defined by p1, p2 and q1, q2 intersect 
//  returns -1 if they dont instersect
//           0 if they are adjacent
//           1 if they intersect
int segmentsIntersect(DFloat* p1, DFloat* p2,    // the first segment
						DFloat* q1, DFloat* q2); // the second segment
//---------------------------------------------------------------------------
//  computes the intersection of the lines defined by p1, p2 and q1, q2  
//---------------------------------------------------------------------------
int intersectionOf(DFloat* p1, DFloat* p2, // the first two points
					 DFloat* q1, DFloat* q2, // the second two points
					 DFloat* intersecion);    // the intersection point
//---------------------------------------------------------------------------
//  computes the equation ax+by+c=0 of the line 
//  defined by two points p1, p2
void findLineEquation(DFloat* p1,  DFloat* p2,           // the two points
					   DFloat* a, DFloat* b, DFloat* c);  // the coefficients of the computed equation


//---------------------------------------------------------------------------*
//     INSERT POINTS, EDGES AND HOLES TO THE STRUCTURE                       *
//---------------------------------------------------------------------------*
//---------------------------------------------------------------------------
//  insert the segment between point1, point2
//  return the index of this segment
  LInt NewSegm(LInt point1, LInt point2);
//---------------------------------------------------------------------------
//  insert a point at (x, y) with weight =  weight
//  returns the index of this point
  LInt NewPoint(DFloat x, DFloat y, float weight);
//---------------------------------------------------------------------------
//  insert a point at (x, y) with weight = 0
//  returns the index of this point
  LInt NewPoint(DFloat x, DFloat y);
//---------------------------------------------------------------------------
//  insert a hole at (x, y) 
//  returns the number of holes so far
  LInt NewHole(DFloat x, DFloat y);

};

#endif
//---------------------------------------------------------------------------*





