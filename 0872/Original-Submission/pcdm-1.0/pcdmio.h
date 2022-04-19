/*********************************************************************/
/*                                                                   */
/* Copyright (C) 2006  Andrey Nikolayevich Chernikov                 */
/*                                                                   */
/* pcdmio.h                                                          */
/*                                                                   */
/* PCDM I/O interface class declaration                              */
/*                                                                   */
/*********************************************************************/

#ifndef __PCDMIO_H__
#define __PCDMIO_H__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "defs.h"

//==================================================================================================
class PcdmIO
{
        public:
            
    int     oneFlag;
    int     numberOfPoints;
    double  *pointList;
    int     numberOfSegments;
    int     *segmentList;
    int     numberOfRegions;
    
    int     *regionNumberOfPoints;
    int     *regionNumberOfSegments;
    int     *regionNumberOfTriangles;
    int     *regionNumberOfHoles;
    int     **regionGlobalSegmentLists;
    int     **regionGlobalTriangleLists;
    
    int     *segmentToRegions;
    
    double  **regionLocalPointLists;
    int     **regionLocalPointIdLists;
    int     **regionLocalSegmentLists;
    int     **regionLocalSegmentIdLists;
    double  **regionLocalHoleLists;
    int     **regionLocalTriangleLists;
    
    int     ***regionSegmentInnerPointList;
    
    double  **regionHoleLists;
        
    double  areaBound;
    double  angleBound;
    int     aggregation;
    int     consequtiveDistribution;
    int     pollFrequency;
    char    *polyDir;
    char    *statDir;
    
    //----------------------------------------------------------------------------------------------
    PcdmIO () : 
            oneFlag(0),
            numberOfPoints(0),
            pointList(0),
            numberOfSegments(0),
            segmentList(0),
            numberOfRegions(0),
            
            regionNumberOfPoints(0),
            regionNumberOfSegments(0),
            regionNumberOfTriangles(0),
            regionNumberOfHoles(0),
            regionGlobalSegmentLists(0),
            regionGlobalTriangleLists(0),
            
            segmentToRegions(0),
            
            regionLocalPointLists(0),
            regionLocalPointIdLists(0),
            regionLocalSegmentLists(0),
            regionLocalSegmentIdLists(0),
            regionLocalTriangleLists(0),
            
            regionSegmentInnerPointList(0),
            
            regionHoleLists(0),
            
            areaBound(-1), 
            angleBound(-1), 
            aggregation(-1),
            consequtiveDistribution(-1),
            pollFrequency(-1),
            polyDir(0),
            statDir(0)
    {}
    
    //----------------------------------------------------------------------------------------------
    void freeVectors (int** v, int n)
    {
        if (v) {
            for (int i = 0; i < n; i++)
                FREE(v[i]);
            FREE(v);
        }
    }
    
    //----------------------------------------------------------------------------------------------
    void freeVectors (double** v, int n)
    {
        if (v) {
            for (int i = 0; i < n; i++)
                FREE(v[i]);
            FREE(v);
        }
    }
    
    //----------------------------------------------------------------------------------------------
    void clear ()
    {
        FREE(pointList);
        FREE(segmentList);
        FREE(regionNumberOfPoints);
        FREE(regionNumberOfSegments);
        FREE(regionNumberOfTriangles);
        FREE(regionNumberOfHoles);
        freeVectors(regionGlobalSegmentLists, numberOfRegions);
        freeVectors(regionGlobalTriangleLists, numberOfRegions);
        FREE(segmentToRegions);
        freeVectors(regionLocalPointLists, numberOfRegions);
        freeVectors(regionLocalPointIdLists, numberOfRegions);
        freeVectors(regionLocalSegmentLists, numberOfRegions);
        freeVectors(regionLocalTriangleLists, numberOfRegions);
        freeVectors(regionHoleLists, numberOfRegions);
    }
    
    //----------------------------------------------------------------------------------------------
    ~PcdmIO ()
    {
        clear();
    }

    //----------------------------------------------------------------------------------------------
    void read (char* inputRegionFileName, char* inputTriangleFileName);

    //----------------------------------------------------------------------------------------------
    void pullApart ();
};

#endif
