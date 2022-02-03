/*********************************************************************/
/*                                                                   */
/* Copyright (C) 2006  Andrey Nikolayevich Chernikov                 */
/*                                                                   */
/* pcdmio.cc                                                         */
/*                                                                   */
/* PCDM I/O interface class implementation                           */
/*                                                                   */
/*********************************************************************/

#include "pcdmio.h"

//==================================================================================================
void PcdmIO::read (char* inputRegionFileName, char* inputTriangleFileName)
{
    FILE        *F;
    int         i, j, d, n, counter, idxIncrement, thisRegionNumberOfSegments;
    
    idxIncrement = ((oneFlag == 0) ? (-1) : (0));
    
    F = fopen(inputRegionFileName, "r");
    assert (F != 0);
    
    fscanf(F, "%d", &(numberOfPoints));
    MALLOC(pointList, double*, numberOfPoints * 2 * sizeof(double));
    
    for (i = 0; i < numberOfPoints; i++)
        fscanf(F, "%*d %lf %lf", &(pointList[i * 2 + 0]), 
                                 &(pointList[i * 2 + 1]));
        
    fscanf(F, "%d", &(numberOfSegments));
    MALLOC(segmentList, int*, numberOfSegments * 2 * sizeof(int));
    for (i = 0; i < numberOfSegments; i++) {
        fscanf(F, "%d", &counter);
        assert(counter == i + 1);
        for (d = 0; d < 2; d++) {
            fscanf(F, "%d", &n);
            segmentList[i * 2 + d] = n + idxIncrement;
        }
    }
    
    fscanf(F, "%d", &(numberOfRegions));
    MALLOC(regionNumberOfSegments, int*, numberOfRegions * sizeof(int));
    MALLOC(regionGlobalSegmentLists, int**, numberOfRegions * sizeof(int*));
    for (i = 0; i < numberOfRegions; i++) {
        fscanf(F, "%d %d", &counter, &(thisRegionNumberOfSegments));
        assert(counter == i + 1);
        regionNumberOfSegments[i] = thisRegionNumberOfSegments;
        MALLOC(regionGlobalSegmentLists[i], int*, thisRegionNumberOfSegments * sizeof(int));
        for (j = 0; j < thisRegionNumberOfSegments; j++) {
            fscanf(F, "%d", &n);
            regionGlobalSegmentLists[i][j] = n + idxIncrement;
        }
    }
    
    MALLOC(regionNumberOfHoles, int*, numberOfRegions * sizeof(int));
    MALLOC(regionHoleLists, double**, numberOfRegions * sizeof(double*));
    for (i = 0; i < numberOfRegions; i++) {
        fscanf(F, "%d %d", &counter, &(regionNumberOfHoles[i]));
        assert(counter == i + 1);
        MALLOC(regionHoleLists[i], double*, regionNumberOfHoles[i] * 2 * sizeof(double));
        for (j = 0; j < regionNumberOfHoles[i]; j++) {
            fscanf(F, "%lf %lf", &(regionHoleLists[i][j * 2 + 0]), 
                                 &(regionHoleLists[i][j * 2 + 1]));
        }
    }
    
    fclose(F);
    
    if (inputTriangleFileName) {
        F = fopen(inputTriangleFileName, "r");
        assert(F != 0);
        MALLOC(regionNumberOfTriangles, int*, numberOfRegions * sizeof(int));
        MALLOC(regionGlobalTriangleLists, int**, numberOfRegions * sizeof(int*));
        for (i = 0; i < numberOfRegions; i++) {
            fscanf(F, "%d %d", &counter, &(regionNumberOfTriangles[i]));
            assert(counter == i + 1);
            MALLOC(regionGlobalTriangleLists[i], int*, 
                   regionNumberOfTriangles[i] * 3 * sizeof(int));
            for (j = 0; j < regionNumberOfTriangles[i]; j++) {
                fscanf(F, "%d", &counter);
                assert(counter == j + 1);
                for (d = 0; d < 3; d++) {
                    fscanf(F, "%d", &n);
                    regionGlobalTriangleLists[i][j * 3 + d] = n + idxIncrement;
                }
            }
        }
        fclose(F);
    }
}

//==================================================================================================
void PcdmIO::pullApart ()
{
    int         i, j, n, d, p, s, idxIncrement;
    int         thisRegionNumberOfSegments, thisRegionNumberOfPoints, thisRegionNumberOfTriangles;
    int         *globalPointMarkers;
    int         *thisRegionGlobalSegments, *thisRegionLocalSegments;
    int         *thisRegionGlobalTriangles, *thisRegionLocalTriangles;
    double      *thisRegionPoints;
    int         *thisRegionPointIdList;
    
    idxIncrement = ((oneFlag == 0) ? (0) : (-1));
    
    MALLOC(globalPointMarkers, int*, numberOfPoints * sizeof(int));
    for (i = 0; i < numberOfPoints; i++)
        globalPointMarkers[i] = -1;
    
    MALLOC(regionNumberOfPoints, int*, numberOfRegions * sizeof(int));
    MALLOC(regionLocalPointLists, double**, numberOfRegions * sizeof(double*));
    MALLOC(regionLocalPointIdLists, int**, numberOfRegions * sizeof(int*));
    MALLOC(regionLocalSegmentLists, int**, numberOfRegions * sizeof(int*));
    MALLOC(segmentToRegions, int*, numberOfSegments * 2 * sizeof(int));
    if (regionNumberOfTriangles != 0)
        MALLOC(regionLocalTriangleLists, int**, numberOfRegions * sizeof(int*));
    
    for (i = 0; i < numberOfSegments * 2; i++)
        segmentToRegions[i] = -10;
    
    for (i = 0; i < numberOfRegions; i++) {
        thisRegionNumberOfSegments = regionNumberOfSegments[i];
        thisRegionGlobalSegments = regionGlobalSegmentLists[i];
        
        // ** count the local points **
        n = 0;
        for (j = 0; j < thisRegionNumberOfSegments; j++) {
            for (d = 0; d < 2; d++) {
                p = segmentList[(thisRegionGlobalSegments[j] + idxIncrement) * 2 + d] + idxIncrement;
                if (globalPointMarkers[p] < 0)
                    globalPointMarkers[p] = n++;
            }
        }
        
        // ** allocate the memory for the local points **
        thisRegionNumberOfPoints = n;
        regionNumberOfPoints[i] = thisRegionNumberOfPoints;
        MALLOC(thisRegionPoints, double*, thisRegionNumberOfPoints * 2 * sizeof(double));
        regionLocalPointLists[i] = thisRegionPoints;
        MALLOC(thisRegionPointIdList, int*, thisRegionNumberOfPoints * sizeof(int));
        regionLocalPointIdLists[i] = thisRegionPointIdList;
        
        // ** clear the local point markers **
        for (j = 0; j < thisRegionNumberOfSegments; j++)
            for (d = 0; d < 2; d++) {
                p = segmentList[(thisRegionGlobalSegments[j] + idxIncrement) * 2 + d] + idxIncrement;
                globalPointMarkers[p] = -1;
            }
        
        // ** copy the local point coordinates **
        n = 0;
        for (j = 0; j < thisRegionNumberOfSegments; j++) {
            for (d = 0; d < 2; d++) {
                p = segmentList[(thisRegionGlobalSegments[j] + idxIncrement) * 2 + d] + idxIncrement;
                if (globalPointMarkers[p] < 0) {
                    globalPointMarkers[p] = n;
                    thisRegionPoints[n * 2 + 0] = pointList[p * 2 + 0];
                    thisRegionPoints[n * 2 + 1] = pointList[p * 2 + 1];
                    thisRegionPointIdList[n] = p - idxIncrement;
                    n++;
                }
            }
        }
        
        MALLOC(thisRegionLocalSegments, int*, thisRegionNumberOfSegments * 2 * sizeof(int));
        regionLocalSegmentLists[i] = thisRegionLocalSegments;
        
        // ** construct the list of local segments as the indexes to local points **
        for (j = 0; j < thisRegionNumberOfSegments; j++)
            for (d = 0; d < 2; d++) {
                p = segmentList[(thisRegionGlobalSegments[j] + idxIncrement) * 2 + d] + idxIncrement;
                thisRegionLocalSegments[j * 2 + d] = globalPointMarkers[p] - idxIncrement;
            }
        
        // ** construct the list of local triangles as the indexes to local points **
        if (regionNumberOfTriangles != 0) {
            assert(regionGlobalTriangleLists != 0);
            thisRegionNumberOfTriangles = regionNumberOfTriangles[i];
            thisRegionGlobalTriangles = regionGlobalTriangleLists[i];
            MALLOC(thisRegionLocalTriangles, int*, thisRegionNumberOfTriangles * 3 * sizeof(int));
            regionLocalTriangleLists[i] = thisRegionLocalTriangles;
            for (j = 0; j < thisRegionNumberOfTriangles; j++)
                for (d = 0; d < 3; d++)
                    thisRegionLocalTriangles[j * 3 + d] = 
                            globalPointMarkers[thisRegionGlobalTriangles[j * 3 + d] + idxIncrement] 
                            - idxIncrement;
        }
        
        // ** clear the local point markers **
        for (j = 0; j < thisRegionNumberOfSegments; j++)
            for (d = 0; d < 2; d++) {
                p = segmentList[(thisRegionGlobalSegments[j] + idxIncrement) * 2 + d] + idxIncrement;
                globalPointMarkers[p] = -1;
            }
            
        for (j = 0; j < thisRegionNumberOfSegments; j++) {
            s = thisRegionGlobalSegments[j] + idxIncrement;
            if (segmentToRegions[s * 2 + 0] < 0)
                segmentToRegions[s * 2 + 0] = i - idxIncrement;
            else {
                assert(segmentToRegions[s * 2 + 1] < 0);
                segmentToRegions[s * 2 + 1] = i - idxIncrement;
            }
        }
    }
    
    FREE(globalPointMarkers);
}

