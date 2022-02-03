/*********************************************************************/
/*                                                                   */
/* Andrey Chernikov, January 2006                                    */
/*                                                                   */
/* cdt.cc                                                            */
/*                                                                   */
/* This file is a utility program for constructing                   */
/*   initial constrained Delaunay triangulations.                    */
/*   Uses Shewchuk's Triangle.                                       */
/*                                                                   */
/*********************************************************************/

#include "pcdmio.h"

extern "C" {
#define REAL double
#define VOID int
#include "triangle.h"
}

//==================================================================================================
void CDT (double* points, int numPoints, 
          int* segments, int numSegments, 
          double* holes, int numHoles,
          int** tris, int* numTris)
{
    triangulateio                                   *in, *out;
    
    in = (triangulateio*) calloc(1, sizeof(triangulateio));
    out = (triangulateio*) calloc(1, sizeof(triangulateio));
    
    in->numberofpoints = numPoints;
    in->pointlist = points;
    
    in->numberofsegments = numSegments;
    in->segmentlist = segments;
    
    in->numberofholes = numHoles;
    in->holelist = holes;
    
    triangulate("pzNYYQ", in, out, 0);
    
    *numTris = out->numberoftriangles;
    *tris = out->trianglelist;
    
    free(in);
    free(out);
}

//==================================================================================================
int main (int argc, char* argv[])
{
    PcdmIO      pcdmio;
    FILE        *F;
    int         i, j, d;
    int         *thisRegionTriangleList, *thisRegionPointIdList;
        
    pcdmio.oneFlag = 0;
    pcdmio.read(argv[1], 0);
    pcdmio.pullApart();
    
    F = fopen(argv[2], "w");
    assert(F != 0);
    
    pcdmio.regionNumberOfTriangles = (int*) malloc(pcdmio.numberOfRegions * sizeof(int));
    pcdmio.regionLocalTriangleLists = (int**) malloc(pcdmio.numberOfRegions * sizeof(int*));
    for (i = 0; i < pcdmio.numberOfRegions; i++) {
        CDT(pcdmio.regionLocalPointLists[i], pcdmio.regionNumberOfPoints[i], 
            pcdmio.regionLocalSegmentLists[i], pcdmio.regionNumberOfSegments[i], 
            pcdmio.regionHoleLists[i], pcdmio.regionNumberOfHoles[i],
            &(pcdmio.regionLocalTriangleLists[i]), &(pcdmio.regionNumberOfTriangles[i]));
        thisRegionTriangleList = pcdmio.regionLocalTriangleLists[i];
        thisRegionPointIdList = pcdmio.regionLocalPointIdLists[i];
        fprintf(F, "%d %d\n", i + 1, pcdmio.regionNumberOfTriangles[i]);
        for (j = 0; j < pcdmio.regionNumberOfTriangles[i]; j++) {
            fprintf(F, "%d", j + 1);
            for (d = 0; d < 3; d++)
                fprintf(F, " %d", thisRegionPointIdList[thisRegionTriangleList[j * 3 + d]] + 1);
            fprintf(F, "\n");
        }
    }
    
    fclose(F);
}
