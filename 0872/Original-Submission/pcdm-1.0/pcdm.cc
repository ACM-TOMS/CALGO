/*********************************************************************/
/*                                                                   */
/* Parallel Constrained Delaunay Meshing (PCDM)                      */
/* Copyright (C) 2006  Andrey Nikolayevich Chernikov                 */
/*                                                                   */
/* This program is free software; you can redistribute it and/or     */
/* modify it under the terms of the GNU General Public License       */
/* as published by the Free Software Foundation; either version 2    */
/* of the License, or (at your option) any later version.            */
/*                                                                   */
/* This program is distributed in the hope that it will be useful,   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     */
/* GNU General Public License for more details.                      */
/*                                                                   */
/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software       */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA     */
/* 02110-1301, USA.                                                  */
/*                                                                   */
/* Author's contact address:                                         */
/*   Andrey N. Chernikov                                             */
/*   Department of Computer Science                                  */
/*   College of William and Mary                                     */
/*   PO Box 8795,                                                    */
/*   Williamsburg, VA 23187, USA                                     */
/*   phone: 757-221-3436                                             */
/*   email: ancher@cs.wm.edu                                         */
/*                                                                   */
/* Supervisor's contact address:                                     */
/*   Nikos P. Chrisochoides                                          */
/*   Department of Computer Science                                  */
/*   College of William and Mary                                     */
/*   PO Box 8795,                                                    */
/*   Williamsburg, VA 23187, USA                                     */
/*   phone: 757-221-3466                                             */
/*   email: nikos@cs.wm.edu                                          */
/*                                                                   */
/*********************************************************************/

//************************ INCLUDES ****************************************************************
#include "pcdm.h"
#ifdef MEM_CHECK
#include <mcheck.h>
#endif
#ifndef SEQUENTIAL
#include <mpi.h>
#endif
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <new>
#include <set>
#include <map>
#include <vector>
#include <deque>
#include <algorithm>
using namespace std;

extern "C" {
#define REAL double
void exactinit();
REAL orient2d(REAL* pa, REAL* pb, REAL* pc);
REAL incircle(REAL* pa, REAL* pb, REAL* pc, REAL* pd);
}

extern "C" {
    
typedef int idxtype;    

void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, 
        int *, int *, int *, int *, int *, idxtype *); 

void METIS_PartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, 
        int *, int *, int *, int *, int *, idxtype *); 
}   
                 
//************************ FORWARD DECLARATIONS ****************************************************

// ** classes **
class Region;

// ** functions **
void request (int tgt, char* buffer, int bytes);
void mpiPoll ();
void processRegion (REGION* reg);
void msgHandler (int src, void* buf, int size);

//************************ GLOBAL VARIABLES ********************************************************

const int NEXT3 [3] = {1, 2, 0};
const int PREV3 [3] = {2, 0, 1};

double                      timer0 [NUM_TIMES];
timeval                     timer1 [NUM_TIMES];
timeval                     timer2;

int                         numProcs = 1, myProc = 0;

map <int, REGION*>          gidToRegion;
map <int, int>              regionToProc;

int                         dummyZero = 0;

// ** PARAMETERS **

int                         argC = 0;
char                        **argV = 0;
double                      angleBound = RADIANS(20.7);
double                      cosAngleBound = cos(angleBound);
double                      areaBound = -1;
char                        *polyDir = 0;
char                        *statDir = 0;
char                        *inputRegionFileName = 0, *inputTriangleFileName = 0;
bool                        consequtiveDistribution = false;
UINT                        msgAggregation = 512;
UINT                        pollFrequency = 512;
vector <int>                procCount;
vector <double>             procWork;

// ** TERMINATION **

bool                        allFinished = false;
bool                        stateActive = true;
bool                        ownFinishedToken = false;
bool                        processingSplitMsg = false;
int                         myColor = 0;

// ** STATISTICS **

int                         numSplitSent = 0, numSplitRecv = 0, numSplitBytes = 0, numInitBytes = 0;
int                         numPointNumReceived = 0;
int                         originalNumPoints = 0, originalNumRegions = 0;
bool                        computeHistogram = false;
        
int                         procedure = -1;

//************************ USAGE *******************************************************************

void printHelp ()
{
    fprintf(stdout, "%s", 
"                                                                               \n\
pcdm <input_region_file_name> <input_triangle_file_name>                        \n\
     -a <double> [-q <double>] [-m <int>] [-o <directory>]                      \n\
     [-s <directory>] [-c] [-f <int>] [-h]                                      \n\
     [{-c} or {-W <num_procs_1> <memory_1> ... <num_procs_k> <memory_k>}]       \n\
                                                                                \n\
    input_region_file_name --- contains the boundaries of the subdomains,       \n\
                               required;                                        \n\
    input_triangle_file_name --- contains the initial constrained Delaunay      \n\
                                 triangulations of the subdomains, required;    \n\
                                                                                \n\
    -a --- the upper bound on the desired triangle area, required;              \n\
    -q --- the lower bound on the desired minimal angle in degrees;             \n\
           the default is 20.1 degrees;                                         \n\
    -m --- maximum number of `split' messages that can be accumulated for a     \n\
           particular subdomain before sending them; 512 by default;            \n\
    -o --- the directory to save the resulting files;                           \n\
           if this parameter is omitted, the files are not saved;               \n\
    -s --- the directory to save the statistics files;                          \n\
           if this parameter is omitted, the files are not saved;               \n\
    -Q --- compute the angle histogram; OFF by default;                         \n\
    -f --- poll frequency, i.e. number of points inserted between calls         \n\
           to poll(); 512 by default;                                           \n\
    -h --- print help and exit;                                                 \n\
    -c --- use consequtive mapping of regions onto processors instead of calling\n\
           MeTiS (e.g. regions 1...i1 -> processor 0, i1+1...i2 ->processor 1,  \n\
           and so on); OFF by default;                                          \n\
    -W --- if given, must be the last parameter; assuming the processes are     \n\
           taken from k subclusters, with memory_i (i=1,...,k) GB of memory     \n\
           each, the subdomains are distributed proportional to the amount of   \n\
           memory available for each process.                                   \n\
           (num_procs_1 + ... + num_procs_k = total number of processes)        \n\
");
    fflush(stdout);
}


//************************ GEOMETRIC ROUTINES ******************************************************

//==================================================================================================
double dist (double* A, double* B)
{
    double dx = A[0] - B[0];
    double dy = A[1] - B[1];
    return sqrt(dx * dx + dy * dy);
}


//==================================================================================================
double angle (double* A, double* B, double* C)
{
    double AminusBx = A[0] - B[0];
    double AminusBy = A[1] - B[1];
    double CminusBx = C[0] - B[0];
    double CminusBy = C[1] - B[1];
    double lenAminusB2 = AminusBx * AminusBx + AminusBy * AminusBy;
    double lenCminusB2 = CminusBx * CminusBx + CminusBy * CminusBy;
    double t = (AminusBx * CminusBx + AminusBy * CminusBy) / sqrt (lenAminusB2 * lenCminusB2);
    if (t >= 1)
        return 0;
    if (t <= -1)
        return PI;
    return acos (t);
}

        
//==================================================================================================
void circumcenter (double* O, double* D, double* A, double* ccenter)
{
    double DminusOx = D[0] - O[0];
    double DminusOy = D[1] - O[1];
    double AminusOx = A[0] - O[0];
    double AminusOy = A[1] - O[1];
    double len2DminusO = DminusOx * DminusOx + DminusOy * DminusOy;
    double len2AminusO = AminusOx * AminusOx + AminusOy * AminusOy;
    double d = 0.5 / (DminusOx * AminusOy - AminusOx * DminusOy);
    ccenter[0] = O[0] - (DminusOy * len2AminusO - AminusOy * len2DminusO) * d;
    ccenter[1] = O[1] + (DminusOx * len2AminusO - AminusOx * len2DminusO) * d;
}


//==================================================================================================
double signedArea (double *A, double *B, double *C)
{
    double a = 0.5 * (  A[0] * (B[1] - C[1]) 
                      + B[0] * (C[1] - A[1]) 
                      + C[0] * (A[1] - B[1]));
    return a;
}


//************************ CLASS DECLARATIONS ******************************************************

//==================================================================================================
template <class T>
class OrderedPair
{
        public:
                
    T n0, n1;
        
    //----------------------------------------------------------------------------------------------
    OrderedPair (T n0_, T n1_)
    {
        if (n0_ < n1_) {
            n0 = n0_;
            n1 = n1_;
        }
        else {
            n0 = n1_;
            n1 = n0_;
        }
    }
    
    //----------------------------------------------------------------------------------------------
    bool operator<(const OrderedPair& other) const
    {
        if (n0 == other.n0)
            return (n1 < other.n1);
        return (n0 < other.n0);
    }
};

//==================================================================================================
class Rational2
{
        private:
                
    UINT num;
    UINT pow;   
    
        public:
                
    //----------------------------------------------------------------------------------------------
    void normalize ()
    {
        while ((num) && (pow)) {
            if (num & ONE) 
                break;
            else {
                num = (num >> 1);
                pow--;
            }
        } 
    }
                            
    //----------------------------------------------------------------------------------------------
    Rational2 (): num (0), pow (0) {}
        
    //----------------------------------------------------------------------------------------------
    Rational2 (UINT n, UINT p) : num (n), pow (p) { normalize (); }
        
    //----------------------------------------------------------------------------------------------
    Rational2 (char* buf) 
    { 
        memcpy(&num, buf, sizeof(UINT));
        memcpy(&pow, &(buf[sizeof(UINT)]), sizeof(UINT));
    }

    //----------------------------------------------------------------------------------------------
    int pack (char* buf) 
    { 
        memcpy (buf, &num, sizeof (UINT));
        memcpy (&(buf [sizeof (UINT)]), &pow, sizeof (UINT));
        return sizeof(Rational2);
    }
        
    //----------------------------------------------------------------------------------------------
    void set (UINT n, UINT p)
    {
        num = n;
        pow = p;
        normalize ();
    }
        
    //----------------------------------------------------------------------------------------------
    void set (UINT n)
    {
        num = n;
        pow = 0;
    }
        
    //----------------------------------------------------------------------------------------------
    void set (const Rational2& r)
    {
        num = r.num;
        pow = r.pow;
    }
    
    //----------------------------------------------------------------------------------------------
    UINT numerator () const  { return num; }
    
    //----------------------------------------------------------------------------------------------
    UINT power()      const  { return pow; }
        
    //----------------------------------------------------------------------------------------------
    double asDouble ()
    {
        double x = ldexp ((double)num, - (int)pow);
        assert ((x >= 0) && (x <= 1));
        return x;
    }
    
    //----------------------------------------------------------------------------------------------
    bool operator< (const Rational2& r) const
    {
        if (pow < r.pow)
            return (num * (ONE << (r.pow - pow)) < r.num);
        
        if (pow > r.pow)
            return (num < r.num * (ONE << (pow - r.pow)));
        
        return (num < r.num);
    }
    
    //----------------------------------------------------------------------------------------------
    bool operator== (const Rational2& r) const
    {
        return ((num == r.num) && (pow == r.pow));
    }
};

//--------------------------------------------------------------------------------------------------
inline Rational2 mid (const Rational2& r0, const Rational2& r1)
{
    if (r0.power () < r1.power ())
        return Rational2 (r0.numerator () * (ONE << (r1.power () - r0.power ())) 
                          + r1.numerator (), r1.power () + 1);
    
    if (r0.power () > r1.power ())
        return Rational2 (r1.numerator () * (ONE << (r0.power () - r1.power ())) 
                          + r0.numerator (), r0.power () + 1);
    
    return Rational2 (r0.numerator () + r1.numerator (), r0.power () + 1);
} 

//==================================================================================================
#define OBJECTS_PER_BLOCK 128
class Pool
{
        public:
                
    int          objSize;
    int          blockSize;
    char         *currentBlock;
    char         *firstBlock;
    int          numUsed;
    char         *nextBlock;
        
    //----------------------------------------------------------------------------------------------
    void init ()
    {
        currentBlock = 0;
        firstBlock = 0;
        numUsed = 0;
        nextBlock = 0;
    }    
    
    //----------------------------------------------------------------------------------------------
    Pool (int objSize_) : 
        objSize(objSize_),
        blockSize(objSize_ * OBJECTS_PER_BLOCK)
    {
        init();
    }
        
    //----------------------------------------------------------------------------------------------
    void clear ()
    {
        if (firstBlock != 0) {
            while (true) {
                if (firstBlock == currentBlock) {
                    FREE(firstBlock);
                    break;
                }
                else {
                    memcpy(&nextBlock, &(firstBlock[blockSize]), sizeof(char*));
                    FREE(firstBlock);
                    firstBlock = nextBlock;
                }
            }
        }
        init();
    }
    
    //----------------------------------------------------------------------------------------------
    ~Pool ()
    {
        clear();
    }
    
    //----------------------------------------------------------------------------------------------
    void* allocate ()
    {
        if ((currentBlock == 0) || (numUsed >= OBJECTS_PER_BLOCK)) {
            CALLOC(nextBlock, char*, blockSize + sizeof(char*));
            if (currentBlock == 0)
                firstBlock = nextBlock;
            else
                memcpy(&(currentBlock[blockSize]), &nextBlock, sizeof(char*));
            currentBlock = nextBlock;
            numUsed = 0;
        }
        return (void*) &(currentBlock[(numUsed++) * objSize]); 
    } 
};

//==================================================================================================
class Point
{
        public:
                
    double              coord[2];
    int                 gid;        // > 0 -- initial global id
                                    // = 0 -- internal point
                                    // < 0 -- edge global id
    int                 lid;
    Rational2           frac;
    
    //----------------------------------------------------------------------------------------------
    void init (double x_, double y_, int gid_)
    {
        coord[0] = x_;
        coord[1] = y_;
        gid = gid_;
    }
    
    //----------------------------------------------------------------------------------------------
    void init (double* coord_, int gid_)
    {
        coord[0] = coord_[0];
        coord[1] = coord_[1];
        gid = gid_;
    }
    
    //----------------------------------------------------------------------------------------------
    void init (double x_, double y_, int gid_, Rational2 frac_)
    {
        init (x_, y_, gid_);
        frac.set (frac_);
    }
    
    //----------------------------------------------------------------------------------------------
    void init (double* coord_, int gid_, Rational2 frac_)
    {
        init (coord_[0], coord_[1], gid_);
        frac.set (frac_);
    }
            
    //----------------------------------------------------------------------------------------------
    void init (double* coord_)
    {
        init (coord_[0], coord_[1], 0);
    }
    
    //----------------------------------------------------------------------------------------------
    Point (double* coord_, int gid_)
    {
        init (coord_[0], coord_[1], gid_);
    }
    
    //----------------------------------------------------------------------------------------------
    Point (double* coord_)
    {
        init (coord_[0], coord_[1], 0);
    }
    
    //----------------------------------------------------------------------------------------------
    Point (double x_, double y_, int gid_)
    {
        init (x_, y_, gid_);
    }
    
    //----------------------------------------------------------------------------------------------
    Point (double x_, double y_)
    {
        init (x_, y_, 0);
    }
    
    //----------------------------------------------------------------------------------------------
    Point (double* coord_, int gid_, Rational2 frac_)
    {
        init (coord_[0], coord_[1], gid_);
        frac.set (frac_);
    }
    
    //----------------------------------------------------------------------------------------------
    Point (double x_, double y_, int gid_, Rational2 frac_)
    {
        init (x_, y_, gid_);
        frac.set (frac_);
    }
};

//==================================================================================================
class TriangleQuality
{
        public:
                
    int                 qIdx;

    //----------------------------------------------------------------------------------------------
    TriangleQuality () : qIdx(0)
    { }
    
    //----------------------------------------------------------------------------------------------
    void computeQuality (double *A, double *B, double *C)
    {
        double  angle, cosA, a;
        int     minAngleIdx;
        
        cosMinAngle(A, B, C, &cosA, &minAngleIdx);
        setMinAngleIdx (minAngleIdx);
        if (cosA >= cosAngleBound) {
            if (cosA >= 1)
                angle = 0;
            else {
                if(cosA <= -1)
                    angle = PI;
                else
                    angle = acos(cosA);
            }
            setBad((int)(angle / angleBound * Q_COUNT));
        }
        else {
            a = area(A, B, C);
            assert(a > 0);
            if (a > areaBound)
                setBig();
        }
    }
    
    //----------------------------------------------------------------------------------------------
    double angleABC (double AminusBx, double AminusBy, double CminusBx, double CminusBy,
                     double lenAminusB2, double lenCminusB2)
    {
        double cosine = (AminusBx * CminusBx + AminusBy * CminusBy) / 
                            sqrt(lenAminusB2 * lenCminusB2);
        if (cosine >= 1)
            return 0;
        else {
            if (cosine <= -1)
                return PI;
            else
                return acos(cosine);
        }
    }
    
    //----------------------------------------------------------------------------------------------
    double cosAngleABC (double AminusBx, double AminusBy, double CminusBx, double CminusBy,
                        double lenAminusB2, double lenCminusB2)
    {
        return ((AminusBx * CminusBx + AminusBy * CminusBy) / sqrt(lenAminusB2 * lenCminusB2));
    }
    
    //----------------------------------------------------------------------------------------------
    double minAngle (double *A, double *B, double *C)
    {
        double AminusBx = A[0] - B[0];
        double AminusBy = A[1] - B[1];
        double CminusBx = C[0] - B[0];
        double CminusBy = C[1] - B[1];
        double AminusCx = A[0] - C[0];
        double AminusCy = A[1] - C[1];
        
        double lenAminusB2 = AminusBx * AminusBx + AminusBy * AminusBy;
        double lenCminusB2 = CminusBx * CminusBx + CminusBy * CminusBy;
        double lenAminusC2 = AminusCx * AminusCx + AminusCy * AminusCy;
        
        if (lenAminusB2 < lenCminusB2) {
            if (lenAminusB2 < lenAminusC2)
                // min angle is C
                return angleABC(AminusCx, AminusCy, -CminusBx, -CminusBy, lenAminusC2, lenCminusB2);
            else 
                // min angle is B
                return angleABC(AminusBx, AminusBy, CminusBx, CminusBy, lenAminusB2, lenCminusB2);
        }
        else {
            if (lenAminusC2 < lenCminusB2) 
                // min angle is B
                return angleABC(AminusBx, AminusBy, CminusBx, CminusBy, lenAminusB2, lenCminusB2);
            else 
                // min angle is A
                return angleABC(AminusBx, AminusBy, AminusCx, AminusCy, lenAminusB2, lenAminusC2);
        }
        
        assert(0);
    }
    
    //----------------------------------------------------------------------------------------------
    void cosMinAngle (double *A, double *B, double *C, double *cosA, int *minAngleIdx)
    {
        double AminusBx = A[0] - B[0];
        double AminusBy = A[1] - B[1];
        double CminusBx = C[0] - B[0];
        double CminusBy = C[1] - B[1];
        double AminusCx = A[0] - C[0];
        double AminusCy = A[1] - C[1];
        
        double lenAminusB2 = AminusBx * AminusBx + AminusBy * AminusBy;
        double lenCminusB2 = CminusBx * CminusBx + CminusBy * CminusBy;
        double lenAminusC2 = AminusCx * AminusCx + AminusCy * AminusCy;
        
        if (lenAminusB2 < lenCminusB2) {
            if (lenAminusB2 < lenAminusC2) {
                // min angle is C
                *cosA = cosAngleABC(AminusCx, AminusCy, -CminusBx, -CminusBy, 
                                    lenAminusC2, lenCminusB2);
                *minAngleIdx = 2;
            }
            else {
                // min angle is B
                *cosA = cosAngleABC(AminusBx, AminusBy, CminusBx, CminusBy, 
                                    lenAminusB2, lenCminusB2);
                *minAngleIdx = 1;
            }
        }
        else {
            if (lenAminusC2 < lenCminusB2) {
                // min angle is B
                *cosA = cosAngleABC(AminusBx, AminusBy, CminusBx, CminusBy, 
                                    lenAminusB2, lenCminusB2);
                *minAngleIdx = 1;
            }
            else {
                // min angle is A
                *cosA = cosAngleABC(AminusBx, AminusBy, AminusCx, AminusCy, 
                                    lenAminusB2, lenAminusC2);
                *minAngleIdx = 0;
            }
        }
    }
    
    //----------------------------------------------------------------------------------------------
    double area (double *A, double *B, double *C)
    {
        return signedArea(A, B, C);
    }
    
    //----------------------------------------------------------------------------------------------
    bool isBad () const
    {
        return ((qIdx & (ONE << LOG_Q_COUNT)) != 0);
    }
    
    int howBad () const
    {
        return (qIdx & (Q_COUNT - 1));
    }
    
    void setBad (int i)
    {
        qIdx |= ((ONE << LOG_Q_COUNT) | i);
    }
    
    //----------------------------------------------------------------------------------------------
    bool isBig () const
    {
        return ((qIdx & (ONE << (LOG_Q_COUNT + 1))) != 0);
    }
    
    void setBig ()
    {
        qIdx |= (ONE << (LOG_Q_COUNT + 1));
    }
    
    //----------------------------------------------------------------------------------------------
    bool isDeleted () const
    {
        return ((qIdx & (ONE << (LOG_Q_COUNT + 4))) != 0);
    }
    
    void setDeleted ()
    {
        qIdx |= (ONE << (LOG_Q_COUNT + 4));
    }
    
    void setNotDeleted ()
    {
        qIdx &= ~(ONE << (LOG_Q_COUNT + 4));
    }
    
    //----------------------------------------------------------------------------------------------
    int getMinAngleIdx () const
    {
        if ((qIdx & (ONE << (LOG_Q_COUNT + 5))) != 0)
            return 0;
        if ((qIdx & (ONE << (LOG_Q_COUNT + 6))) != 0)
            return 1;
        if ((qIdx & (ONE << (LOG_Q_COUNT + 7))) != 0)
            return 2;
        assert(0);
    }
    
    //----------------------------------------------------------------------------------------------
    void setMinAngleIdx (int idx)
    {
        if (idx == 0) {
            qIdx |= (ONE << (LOG_Q_COUNT + 5));
            return;
        }
        if (idx == 1) {
            qIdx |= (ONE << (LOG_Q_COUNT + 6));
            return;
        }
        if (idx == 2) {
            qIdx |= (ONE << (LOG_Q_COUNT + 7));
            return;
        }
        assert(0);
    }
};

//==================================================================================================
class Triangle: public TriangleQuality
{
    
        public:
                
    Point*              points[3];
    Triangle*           neis[3];
    
    //----------------------------------------------------------------------------------------------
    void init (Point* p0, Point* p1, Point* p2, int qIdx_)
    {
        assert ((p0 != 0) && (p1 != 0) && (p2 != 0));
        points [0] = p0;
        points [1] = p1;
        points [2] = p2;
        neis   [0] = 0;
        neis   [1] = 0;
        neis   [2] = 0;
        qIdx       = qIdx_;
    }
    
    //----------------------------------------------------------------------------------------------
    void init (Point* p0, Point* p1, Point* p2, int qIdx_, bool checkQuality)
    {
        init(p0, p1, p2, qIdx_);

        if (checkQuality) {
            computeQuality(p0->coord, p1->coord, p2->coord);
        }
    }
    
    //----------------------------------------------------------------------------------------------
    Triangle (Point* p0, Point* p1, Point* p2, int qIdx_)
    {
        init(p0, p1, p2, qIdx_);
    }
    
    //----------------------------------------------------------------------------------------------
    Triangle (Point* p0, Point* p1, Point* p2, int qIdx_, bool checkQuality)
    {
        init(p0, p1, p2, qIdx_, checkQuality);
    }
    
    //--------------------------------------------------------------------------------------------------
    int insideTest (double* coord)
    {
        for (int i = 0; i < 3; i++) {
            if (orient2d(points[PREV3[i]]->coord, points[i]->coord, coord) < 0)
                return NEXT3[i];
        }
        return -1;
    }

    //----------------------------------------------------------------------------------------------
    double computeArea ()
    {
        return area(points[0]->coord, points[1]->coord, points[2]->coord);
    }
    
    //----------------------------------------------------------------------------------------------
    void angles (double* a)
    {
        a[0] = angle (points[2]->coord, points[0]->coord, points[1]->coord);
        a[1] = angle (points[0]->coord, points[1]->coord, points[2]->coord);
        a[2] = PI - a[0] - a[1];
    }
    
    //----------------------------------------------------------------------------------------------
    void ccenter (double* coord)
    {
        circumcenter (points[0]->coord, points[1]->coord, points[2]->coord, coord);
    }
    
    //----------------------------------------------------------------------------------------------
    void setNei (Point* p0, Point* p1, Triangle* nei)
    {
        for (int i = 0; i < 3; i++)
            if ((points[i] == p0) && (points[NEXT3[i]] == p1)) {
                neis[PREV3[i]] = nei;
                return;
            }
        assert (0);
    }
    
    //----------------------------------------------------------------------------------------------
    Triangle* getNei (Point* p0, Point* p1)
    {
        for (int i = 0; i < 3; i++)
            if ((points[i] == p0) && (points[NEXT3[i]] == p1))
                return neis[PREV3[i]];
        assert (0);
    }
    
    //----------------------------------------------------------------------------------------------
    int pointIndex (Point* point)
    {
        for (int i = 0; i < 3; i++)
            if (points[i] == point)
                return i;
        assert (0);
    }
};

//==================================================================================================
struct ltpt
{
    bool operator()(const Point* p0, const Point* p1) const
    {
        return (p0->frac < p1->frac);
    }
};

//--------------------------------------------------------------------------------------------------
class Segment
{
        public:
    
    typedef set <Point*, ltpt>                      PointContainer;
    typedef map < OrderedPair<Point*>, Triangle* >  PointsToTriangleMap;
                
    int                                     gid;
    Point                                   *point0, *point1;
    int                                     rid0, rid1;
    PointContainer                          points;
    PointsToTriangleMap                     pointsToTri;  
          
    //----------------------------------------------------------------------------------------------
    void init (int gid_, Point* point0_, Point* point1_, int rid0_, int rid1_)
    {
        gid    = gid_;
        point0 = point0_;
        point1 = point1_;
        rid0   = rid0_;
        rid1   = rid1_;
    }
     
    //----------------------------------------------------------------------------------------------
    Segment (int gid_, Point* point0_, Point* point1_, int rid0_, int rid1_)
    {
        init(gid_, point0_, point1_, rid0_, rid1_);
    }

    //----------------------------------------------------------------------------------------------
    Segment (int gid_, Point* point0_, Point* point1_)
    {
        init(gid_, point0_, point1_, -1, -1);
    }
    
    //----------------------------------------------------------------------------------------------
    Segment* copy ()
    {
        Segment* segment;
        
        NEW (segment, Segment (gid, point0, point1, rid0, rid1));
        return segment;
    }
    
    //----------------------------------------------------------------------------------------------
    void addRegion (int rid)
    {
        if (rid0 < 0) {
            rid0 = rid;
        }
        else {
            if (rid1 < 0) {
                rid1 = rid;
            }
            else
                assert (0);
        }
    }
    
    //----------------------------------------------------------------------------------------------
    bool insertPoint (Point* point)
    {
        return (points.insert(point)).second;
    }
    
    //----------------------------------------------------------------------------------------------
    int otherRegion (int rid)
    {
        if (rid == rid0) {
            return rid1;
        }
        if (rid == rid1) {
            return rid0;
        }
        assert (0);
    }
    
    //----------------------------------------------------------------------------------------------
    double length ()
    {
        return dist(point0->coord, point1->coord);
    }
    
    //----------------------------------------------------------------------------------------------
    void assignPointNumbers (int * nextPointNumber)
    {
        PointContainer::iterator    pI;
        
        FOR(pI, points) {
            assert((*pI)->gid < 0);
            (*pI)->gid = (*nextPointNumber)++;
        }
    }
};

//==================================================================================================
class CavityItem
{
        public:
                
    Triangle    *tri, *nei;
    int         triSide;
    Point       *p0, *p1;
    CavityItem  *next;
    
    //----------------------------------------------------------------------------------------------
    void init (Triangle* tri_, Triangle* nei_, int triSide_)
    {
        tri           = tri_;
        nei           = nei_;
        triSide       = triSide_;
        p0            = tri->points[NEXT3[triSide]];
        p1            = tri->points[PREV3[triSide]];
        next          = 0;
    }
};

//==================================================================================================
class Cavity
{
        public:
                
    Triangle                        *trisSelected[MAX_CAVITY_SIZE], *trisRejected[MAX_CAVITY_SIZE];
    CavityItem                      items[MAX_CAVITY_SIZE];
    int                             numSelected, numRejected, numItems;
    
    //----------------------------------------------------------------------------------------------
    Cavity () : 
        numSelected(0), numRejected(0), numItems(0)
    {  }
    
    //----------------------------------------------------------------------------------------------
    void clear ()
    {
        numSelected  = 0;
        numRejected  = 0;
        numItems     = 0;
    }    
    
    //----------------------------------------------------------------------------------------------
    void setSelected (Triangle* tri)
    {
        assert(numSelected < MAX_CAVITY_SIZE);
        trisSelected[numSelected++] = tri;
    }
    
    //----------------------------------------------------------------------------------------------
    void setRejected (Triangle* tri)
    {
        assert(numRejected < MAX_CAVITY_SIZE);
        trisRejected[numRejected++] = tri;
    }
    
    //----------------------------------------------------------------------------------------------
    bool isSelected (Triangle* tri)
    {
        for (int i = 0; i < numSelected; i++) {
            if (trisSelected[i] == tri)
                return true;
        }
        return false;
    }
    
    //----------------------------------------------------------------------------------------------
    bool isRejected (Triangle* tri)
    {
        for (int i = 0; i < numRejected; i++) {
            if (trisRejected[i] == tri)
                return true;
        }
        return false;
    }
            
    //----------------------------------------------------------------------------------------------
    void construct (Triangle* base, double* xy)
    {
        int                 i, j, qSize = 0;
        Triangle            *tri, *nei;
        static Triangle     *Q[MAX_CAVITY_SIZE];
        CavityItem          *ci, *cj;
                
        setSelected(base);
        Q[qSize++] = base;
        while (qSize > 0) {
            tri = Q[--qSize];
            for (i = 0; i < 3; i++) {
                nei = tri->neis[i];
                if (nei) {
                    if (!isSelected(nei)) {
                        if (isRejected(nei)) {
                            addItem(tri, nei, i);
                        }
                        else {
                            if (incircle(nei->points[0]->coord,
                                         nei->points[1]->coord, 
                                         nei->points[2]->coord, xy) >= 0) {
                                setSelected(nei);
                                
                                assert(qSize < MAX_CAVITY_SIZE);
                                Q[qSize++] = nei;
                            }
                            else {
                                setRejected(nei);
                                addItem(tri, nei, i);
                            }
                        }
                    }
                }
                else {
                    addItem(tri, nei, i);
                }
            }
        }
        
        for (i = 0; i < numItems - 1; i++) {
            ci = &(items[i]);
            for (j = i + 1; j < numItems; j++) {
                cj = &(items[j]);
                if (ci->p0 == cj->p1) {
                    cj->next = ci;
                }
                else
                    if (ci->p1 == cj->p0) {
                        ci->next = cj;
                    }
            }
        }
    }
    
    //----------------------------------------------------------------------------------------------
    void addItem (Triangle* tri_, Triangle* nei_, int triSide_)
    {
        assert(numItems < MAX_CAVITY_SIZE);
        items[numItems++].init(tri_, nei_, triSide_);
    }
};

//==================================================================================================
struct gti
{
    bool operator()(const int i0, const int i1) const
    {
        return (i0 > i1);
    }
};

//--------------------------------------------------------------------------------------------------
class Region
{
        public:
    
    typedef deque <Point*>                      PointContainer;
    typedef Pool                                PointStorage;
    typedef deque <Triangle*>                   TriangleContainer;
    typedef Pool                                TriangleStorage;
    typedef pair<double, double>                HoleType;
    typedef vector <HoleType>                   HoleContainer;
    typedef map <int, Segment*>                 GidSegmentMap;
    typedef map <OrderedPair<int>, Segment*>    IdpairSegmentMap;
    typedef vector <Point*>                     SplitPointContainer;
    typedef map <int, SplitPointContainer*>     RegionPointsMap;
                            
        private:
                
    TriangleContainer                           bigTris;
    TriangleContainer                           badTris[Q_COUNT];
    set <int>                                   minQ;
    int                                         numTris;
    Cavity                                      cavity;
    Point                                       *splitSegmentPoint;
    
        public:
 
    int                                         gid;
    double                                      area;
    PointContainer                              points;
    PointStorage                                pointStorage;
    TriangleContainer                           tris;
    TriangleStorage                             triStorage;
    TriangleContainer                           tris0;
    HoleContainer                               holes;
    GidSegmentMap                               gidToSegment;
    IdpairSegmentMap                            idpairToSegment;
    RegionPointsMap                             regionToPoints;
    bool                                        refinementFinished;
    
        public:
                
    //----------------------------------------------------------------------------------------------
    int getGid () const { return gid; }
    
    //----------------------------------------------------------------------------------------------
    double getArea () const { return area; }
                
    //----------------------------------------------------------------------------------------------
    void computeArea () 
    {
        TriangleContainer::iterator     tI;
        Triangle                        *tri;
        
        area = 0;
        FOR (tI, tris) {
            tri = *tI;
            if (!(tri->isDeleted())) {
                area += tri->computeArea();
            }
        }
    }
                
    //----------------------------------------------------------------------------------------------
    PointContainer & getPoints () { return points; }
                
    //----------------------------------------------------------------------------------------------
    TriangleContainer & getTriangles () { return tris; }
                
    //----------------------------------------------------------------------------------------------
    GidSegmentMap & getGidSegmentMap () { return gidToSegment; }
                
    //----------------------------------------------------------------------------------------------
    bool isRefinementFinished () const { return refinementFinished; }
                
    //----------------------------------------------------------------------------------------------
    void setRefinementFinished (bool v) { refinementFinished = v; }
                
    //----------------------------------------------------------------------------------------------
    void init (int gid_)
    {
        gid     = gid_;
        numTris = 0;
        refinementFinished = false;
        splitSegmentPoint = 0;
    }
    
    //----------------------------------------------------------------------------------------------
    Region (int gid_) : pointStorage(sizeof(Point)), triStorage(sizeof(Triangle))
    {
        init (gid_);
    }
    
    //----------------------------------------------------------------------------------------------
    void clear ()
    {
        deque <Triangle*> :: iterator       tI;
        deque <Point*> :: iterator          pI;
        map <int, Segment*> :: iterator     gsI;
        int                                 i;
        RegionPointsMap :: iterator         rI;
        vector <Point*>                     *q;
        
        pointStorage.clear();
        points.clear();
        
        triStorage.clear();
        tris.clear();
        tris0.clear();
        numTris = 0;
        
        FOR (gsI, gidToSegment) 
            delete (*gsI).second;
        gidToSegment.clear();
        idpairToSegment.clear();
        
        bigTris.clear();
        for (i = 0; i < Q_COUNT; i++)
            badTris[i].clear();
        minQ.clear();
        
        FOR (rI, regionToPoints) {
            q = (*rI).second;
            assert(q->size() == 0);
            delete q;
        }
        regionToPoints.clear();
    }
    
    //----------------------------------------------------------------------------------------------
    ~Region ()
    {
        clear();
    }
    
    //----------------------------------------------------------------------------------------------
    bool check ()
    {
        int                                 i, n;
        Triangle                            *tri, *nei;
        deque <Triangle*> :: iterator       tI;
        deque <Point*> :: iterator          pI;
        
        FOR (pI, points) {
            (*pI)->lid = 0;
        }
        
        n = 0;
        FOR (tI, tris) {
            tri = *tI;
            if (! tri->isDeleted ()) {
                n++;
                for (i = 0; i < 3; i++) {
                    assert (tri->points[i] != 0);
                    tri->points[i]->lid++;
                }
            }
        }
        assert (numTris == n);
        
        n = 0;
        FOR (pI, points) {
            n += (*pI)->lid;
        }
        if (3 * numTris != n) {
            assert (0);
        }
        
        FOR (tI, tris) {
            tri = *tI;
            if (! tri->isDeleted ()) {
                for (i = 0; i < 3; i++) {
                    nei = tri->neis[i];
                    if (nei != 0) 
                        assert (nei->getNei (tri->points[PREV3[i]], tri->points[NEXT3[i]]) == tri);
                }
            }
        }
        
        return true;
    }
    
    //----------------------------------------------------------------------------------------------
    bool addPoint (Point* point)
    {
        PUSH_BACK (points, point);
        if (point->gid < 0) {
            return (gidToSegment[- point->gid]->insertPoint(point));
        }
        return true;
    }
    
    //----------------------------------------------------------------------------------------------
    Segment* getSegment (Point* p0, Point* p1)
    {
        IdpairSegmentMap :: iterator                isI;
        Segment*                                    segment;
        
        if ((p0->gid != 0) && (p1->gid != 0)) {
            if (p0->gid > 0) {
                if (p1->gid > 0) {
                    isI = idpairToSegment.find (OrderedPair<int>(p0->gid, p1->gid));
                    if (isI != idpairToSegment.end ())
                        return (*isI).second;
                }
                else {
                    segment = gidToSegment [- p1->gid];
                    if ((segment->point0->gid == p0->gid) || (segment->point1->gid == p0->gid))
                        return segment;
                }
            }
            else {
                if (p1->gid > 0) {
                    segment = gidToSegment [- p0->gid];
                    if ((segment->point0->gid == p1->gid) || (segment->point1->gid == p1->gid))
                        return segment;
                }
                else {
                    if (p0->gid == p1->gid)
                        return gidToSegment [- p0->gid];
                }
            }
        }
        return 0;
    }
    
    //----------------------------------------------------------------------------------------------
    Triangle* addTriangle (Point* p0_, Point* p1_, Point* p2_, int qIdx_, bool checkQuality_)
    {
        int                                                     i, minAngleIdx;
        Segment*                                                segments[3];
        Point                                                   *p0, *p1;
        pair <Segment::PointsToTriangleMap::iterator, bool>     insertResult;
        Triangle                                                *tri, *oldTri;
        
        if (tris0.empty ()) {
            tri = (Triangle*) triStorage.allocate();
            tri->init(p0_, p1_, p2_, qIdx_, checkQuality_);
            PUSH_BACK(tris, tri);
        }
        else {
            tri = tris0.back ();
            tris0.pop_back();
            tri->init(p0_, p1_, p2_, qIdx_, checkQuality_);
        }
                
        for (i = 0; i < 3; i++) {
            p0 = tri->points[NEXT3[i]];
            p1 = tri->points[PREV3[i]];
            segments[i] = getSegment(p0, p1);
            if (segments[i] != 0) {
                insertResult = segments[i]->pointsToTri.insert(
                        pair<OrderedPair<Point*>, Triangle*>(OrderedPair<Point*>(p0, p1), tri));
                assert(insertResult.second);
            }
        }
        
        minAngleIdx = tri->getMinAngleIdx();
        if ((tri->isBad ()) && 
            ((segments[NEXT3[minAngleIdx]] == 0) || (segments[PREV3[minAngleIdx]] == 0))) {
                i = tri->howBad ();
                PUSH_FRONT (badTris[i], tri);
                minQ.insert (i);
            }
        else {
            if (tri->isBig ()) {
                PUSH_FRONT (bigTris, tri);
            }
        }
        
        numTris++;
        
        return tri;
    }
        
    //----------------------------------------------------------------------------------------------
    void deleteTriangle (Triangle* tri)
    {
        int                                     i, n;
        Segment*                                segment;
        Point                                   *p0, *p1;
        
        assert (tri != 0);
        assert (! tri->isDeleted ());
        tri->setDeleted ();
        p1 = tri->points[2];
        for (i = 0; i < 3; i++) {
            p0 = tri->points [i];
            segment = getSegment(p0, p1);
            if (segment != 0) {
                n = segment->pointsToTri.erase (OrderedPair<Point*>(p0, p1));
                if (n != 1) {
                    assert (0);
                }
            }
            p1 = p0;
        }
        
        if ((! tri->isBig ()) && (! tri->isBad ()))
            PUSH_BACK (tris0, tri);
        
        numTris--;
    }

    //----------------------------------------------------------------------------------------------
    Triangle* topQ ()
    {
        Triangle*   tri = 0;
        int         i;
        
        while ((tri == 0) && (! minQ.empty ())) {
            i = *(minQ.begin ());
            tri = badTris[i].back ();
            if (tri->isDeleted ()) {
                PUSH_BACK (tris0, tri);
                badTris[i].pop_back ();
                if (badTris[i].empty ())
                    minQ.erase (minQ.begin ());
                tri = 0;
            }
        };
            
        while ((tri == 0) && (! bigTris.empty ())) {
            tri = bigTris.back ();
            if (tri->isDeleted ()) {
                PUSH_BACK (tris0, tri);
                bigTris.pop_back ();
                tri = 0;
            }
        }
        
        
        return tri;
    }
    
    //----------------------------------------------------------------------------------------------
    void addSegment (Segment* segment)
    {
        gidToSegment[segment->gid] = segment;
        idpairToSegment[OrderedPair<int>(segment->point0->gid, segment->point1->gid)] = segment;
    }
    
    //----------------------------------------------------------------------------------------------
    void triangulateCavity (Point* point, void* segmentPoint0, void* segmentPoint1)
    {
        int                                     i, numNewTris = 0;
        static Triangle                         *newTris[MAX_CAVITY_SIZE];
        CavityItem                              *ci, *first;
        Triangle                                *tri, *nei;
        
        for (i = 0; i < cavity.numSelected; i++)
            deleteTriangle(cavity.trisSelected[i]);
                    
        first = &(cavity.items[0]);
        ci = first;
        do {
            if ((ci->p0 != segmentPoint0) || (ci->p1 != segmentPoint1)) {
                tri = addTriangle(point, ci->p0, ci->p1, 0, true);
                nei = ci->nei;
                if (nei != 0) {
                    nei->setNei(ci->p1, ci->p0, tri);
                    tri->neis[0] = nei;
                }
                assert(numNewTris < MAX_CAVITY_SIZE);
                newTris[numNewTris++] = tri;
            }
            else {
                assert(numNewTris < MAX_CAVITY_SIZE);
                newTris[numNewTris++] = 0;
            }
            ci = ci->next;
        } while (ci != first);
            
        nei = newTris[numNewTris - 1];
        for (i = 0; i < numNewTris; i++) {
            tri = newTris[i];
            if ((tri != 0) && (nei != 0)) {
                tri->neis[2] = nei;
                nei->neis[1] = tri;
            }
            nei = tri;
        }
    }
    
    //----------------------------------------------------------------------------------------------
    void splitSegmentSrc (Triangle* tri, int side)
    {
        int                             neiId;
        Point                           *p0, *p1, *point;
        Segment*                        segment;
        Rational2                       frac, frac0, frac1;
        double                          f, coord[2];
        bool                            res;
        
        procedure = 1;
        
        p0 = tri->points[NEXT3[side]];
        p1 = tri->points[PREV3[side]];
        
        segment = getSegment (p0, p1);
        assert (segment != 0);

        if (p0 == segment->point0) {
            frac0.set (0);
            if (p1 == segment->point1)
                frac1.set (1);
            else 
                frac1.set (p1->frac);
        }
        else {
            if (p0 == segment->point1) {
                frac0.set (1);
                if (p1 == segment->point0)
                    frac1.set (0);
                else
                    frac1.set (p1->frac);
            }
            else {
                frac0.set (p0->frac);
                if (p1 == segment->point0)
                    frac1.set (0);
                else {
                    if (p1 == segment->point1)
                        frac1.set (1);
                    else
                        frac1.set (p1->frac);
                }
            }
        }
        
        frac = mid(frac0, frac1);
        f = frac.asDouble();
        coord[0] = segment->point0->coord[0] + 
                   (segment->point1->coord[0] - segment->point0->coord[0]) * f;
        coord[1] = segment->point0->coord[1] + 
                   (segment->point1->coord[1] - segment->point0->coord[1]) * f;
        point = (Point*) pointStorage.allocate();
        point->init(coord, - segment->gid, frac);
        res = addPoint(point);
        assert(res);
        
        cavity.construct(tri, point->coord);
        triangulateCavity(point, p0, p1);
        cavity.clear();
        
        neiId = segment->otherRegion(gid);
        if (neiId >= 0) 
            splitMessage(neiId, point);
    }
    
    //----------------------------------------------------------------------------------------------
    void flushQueue (int neiId, vector <Point*> *q)
    {
        int                             size, curSize, segmentGid;
        int                             msg[MSG_HEADER_INTS];
        char*                           buf;
        vector <Point*> :: iterator     pI;
        map <int, int> :: iterator      rpI;
        
        msg[0] = MSG_SPLIT;
        msg[1] = neiId;
        msg[2] = q->size();
        msg[3] = myProc;
        
        size = MSG_HEADER_SIZE + q->size() * (sizeof(int) + sizeof(Rational2));
        MALLOC(buf, char*, size);
        
        memcpy(buf, &msg, MSG_HEADER_SIZE);
        curSize = MSG_HEADER_SIZE;
        FOR (pI, *q) {
            segmentGid = - (*pI)->gid;
            memcpy (&(buf [curSize]), &segmentGid, sizeof (int));
            curSize += sizeof (int);
            curSize += (*pI)->frac.pack (&(buf [curSize]));
        }
        q->clear();

        rpI = regionToProc.find (neiId);
        request(((rpI == regionToProc.end ()) ? 0 : (*rpI).second), buf, size);
        numSplitSent++;  // EWD Rule 0
        numSplitBytes += size;
        
        FREE(buf);
    }
    
    //----------------------------------------------------------------------------------------------
    void cleanAllQueues ()
    {
        RegionPointsMap :: iterator                      rI;
        vector <Point*>                                  *q;
        int                                              neiId;
        
        FOR (rI, regionToPoints) {
            neiId = (*rI).first;
            q = (*rI).second;
            if (q->size() > 0) {
                flushQueue(neiId, q);
            }
        }
    }
    
    //----------------------------------------------------------------------------------------------
    void splitMessage (int neiId, Point* point)
    {
        vector <Point*>                             *q;
        RegionPointsMap :: iterator                 rI;
     
        rI = regionToPoints.find(neiId);
        if (rI == regionToPoints.end ()) {
            NEW(q, SplitPointContainer);
            regionToPoints[neiId] = q;
        }
        else {
            q = (*rI).second;
        }
     
        PUSH_BACK(*q, point);
        if (q->size() >= msgAggregation)
            flushQueue(neiId, q);
    }
    
    //----------------------------------------------------------------------------------------------
    void insertPointOnSegment (Point* point, Segment* segment, Point* p0, Point* p1)
    {
        pair < set <Point*, ltpt> :: iterator, bool >   insertResult;
        Segment::PointsToTriangleMap :: iterator        ptI;  
        Point                                           *tmp;
        double                                          f;
        Triangle*                                       baseTri;
        int                                             idx;
        
        insertResult = segment->points.insert(point);
        assert(insertResult.second);
        
        ptI = segment->pointsToTri.find(OrderedPair<Point*>(p0, p1));
        assert(ptI != segment->pointsToTri.end());
        
        baseTri = (*ptI).second;
        idx = baseTri->pointIndex(p0);
        if (baseTri->points[NEXT3[idx]] != p1) {
            assert(baseTri->points[PREV3[idx]] == p1);
            tmp = p0;
            p0 = p1;
            p1 = tmp;
        }
        
        f = point->frac.asDouble ();
        point->coord[0] = segment->point0->coord[0] + 
                          (segment->point1->coord[0] - segment->point0->coord[0]) * f;
        point->coord[1] = segment->point0->coord[1] + 
                          (segment->point1->coord[1] - segment->point0->coord[1]) * f;
        addPoint(point);
        
        cavity.construct(baseTri, point->coord);
        triangulateCavity(point, p0, p1);
        cavity.clear();
    }
    
    //----------------------------------------------------------------------------------------------
    void splitSegmentsDst (int count, char* buf)
    {
        int                                         i, j, segmentGid, curSize = 0, iterCount;
        Segment                                     *segment;
        map <int, Segment*> :: iterator             gsI;
        set <Point*, ltpt> :: iterator              pI;
        Point                                       *point, *p0, *p1;
        Rational2                                   frac, history[16], history0[16], history1[16];
        
        procedure = 2;
        
        for (i = 0; i < count; i++) {
            memcpy(&segmentGid, &(buf[curSize]), sizeof(int));
            curSize += sizeof(int);
            
            gsI = gidToSegment.find(segmentGid);
            assert(gsI != gidToSegment.end());
            segment = (*gsI).second;
            
            if (splitSegmentPoint == 0)
                splitSegmentPoint = (Point*) pointStorage.allocate();
            splitSegmentPoint->init(0, 0, -segmentGid, Rational2(&(buf[curSize])));
            curSize += sizeof(Rational2);
            
            if (segment->points.find(splitSegmentPoint) == segment->points.end()) {
                
                pI = segment->points.upper_bound(splitSegmentPoint);
                if (pI == segment->points.end()) {
                    p1 = segment->point1;
                    if (segment->points.size() == 0) {
                        p0 = segment->point0;
                        frac.set(1, 1);
                    }
                    else {
                        p0 = *(--pI);
                        frac.set(p0->frac.numerator() + (ONE << p0->frac.power()), 
                                 p0->frac.power() + 1);
                    }
                }
                else {
                    p1 = *(pI);
                    if (pI == segment->points.begin()) {
                        p0 = segment->point0;
                        frac.set(p1->frac.numerator(), p1->frac.power() + 1);
                    }
                    else {
                        p0 = *(--pI);
                        frac = mid(p0->frac, p1->frac);
                    }
                }
                
                iterCount = 0;
                while (true) {
                    if (frac == splitSegmentPoint->frac) {
                        insertPointOnSegment(splitSegmentPoint, segment, p0, p1);
                        break;
                    }
                    else { 
                        point = (Point*) pointStorage.allocate();
                        point->init(0, 0, -segmentGid, frac);
                        insertPointOnSegment(point, segment, p0, p1);
                        if (frac < splitSegmentPoint->frac) {
                            p0 = point;
                            if (p1 == segment->point1)
                                frac.set(frac.numerator() + (ONE << frac.power()), 
                                         frac.power() + 1);
                            else    
                                frac = mid(frac, p1->frac);
                        }
                        else {
                            p1 = point;
                            if (p0 == segment->point0)
                                frac.set(frac.numerator(), frac.power() + 1);
                            else
                                frac = mid(p0->frac, frac);
                        }
                    }
                }
                
                splitSegmentPoint = 0;
            }
        }
        
        processRegion(this);
    }
    
    //----------------------------------------------------------------------------------------------
    int splitTriangle (Triangle*& tri, double* coord)
    {
        CavityItem*                             ci;
        Point*                                  point;
        int                                     i, side;
        
        procedure = 0;
        
        cavity.construct(tri, coord);
        
        for (i = 0; i < cavity.numItems; i++) {
            ci = &(cavity.items[i]);
            if (ci->nei == 0) {
                if (angle(ci->p0->coord, coord, ci->p1->coord) >= 2 * PI / 3) {
                    tri = ci->tri;
                    side = ci->triSide;
                    cavity.clear();
                    return side;
                }
            }
        }
        
        point = (Point*) pointStorage.allocate();
        point->init(coord);
        addPoint(point);
        triangulateCavity(point, &dummyZero, &dummyZero);
        cavity.clear();
        return -1;
    }
    
    //----------------------------------------------------------------------------------------------
    bool refine (int count)
    {
        BEGIN_TIMING (TIME_REFINE_PCDM);
        
        int                                 side, i = 0;
        bool                                notdone = true;
        Triangle                            *tri, *nei;
        double                              coord[2];

        while ((i++ < count) && (notdone)) {
            tri = topQ();
            if (tri) {
                tri->ccenter(coord);
                while (true) {
                    side = tri->insideTest(coord);
                    if (side >= 0) {
                        // if outside
                        nei = tri->neis[side];
                        if (nei != 0) {
                            // if the neighbor exists
                            tri = nei;
                        }
                        else {
                            // the point lies across the boundary
                            // ** split the corresponding segment **
                            splitSegmentSrc(tri, side);
                            break;
                        }
                    }
                    else {
                        // the point lies inside a triangle, not necessarily the one we began with
                        side = splitTriangle(tri, coord);
                        if (side >= 0) {
                            splitSegmentSrc(tri, side);
                        }
                        break;
                    }
                }
            }
            else
                notdone = false;
        }
        
        END_TIMING (TIME_REFINE_PCDM);
        
        return notdone;
    }
    
    //----------------------------------------------------------------------------------------------
    void copyTo (void* dst, char* buf, int size, int * curSize)
    {
        memcpy (dst, &(buf[*curSize]), size); 
        (*curSize) += size;
    }
    
    //----------------------------------------------------------------------------------------------
    Region (char* buf, int & size, int src) : 
            pointStorage(sizeof(Point)), triStorage(sizeof(Triangle))
    {
        int                                 i, j, k, n, idx, id, lid[7], curSize = 0;
        Point*                              point;
        Segment*                            segment;
        vector <Point*>                     newPoints;
        vector <Point*> :: iterator         pI;
        double                              xy[2];
        map <int, Segment*> :: iterator     gsI;
        Triangle                            **newTris, *tri;
        int                                 *triNeis;

        // ** gid **
        
        copyTo(&id, buf, sizeof(int), &curSize);
        init (id);
        
        // ** refinementFinished **
        
        copyTo(&refinementFinished, buf, sizeof(bool), &curSize);
        
        // ** area **
        
        copyTo(&area, buf, sizeof(double), &curSize);
        
        // ** points **
        
        copyTo(&n, buf, sizeof(int), &curSize);
        
        for (i = 0; i < n; i++) {
            copyTo(&id, buf, sizeof(int), &curSize);
            
            copyTo(xy, buf, 2 * sizeof(double), &curSize);
            
            point = (Point*) pointStorage.allocate();
            point->init(xy[0], xy[1], id, Rational2 (&(buf[curSize])));
            curSize += sizeof(Rational2);
            
            PUSH_BACK (newPoints, point);
        }
        
        // ** segments **
        
        copyTo(&n, buf, sizeof(int), &curSize);
        
        for (i = 0; i < n; i++) {
            copyTo(lid, buf, 5 * sizeof (int), &curSize);
            NEW (segment, Segment(lid[2], newPoints[lid[0]], newPoints[lid[1]], lid[3], lid[4]));
            addSegment(segment);
        }

        // ** triangles **
        
        copyTo(&n, buf, sizeof(int), &curSize);
        
        MALLOC(newTris, Triangle**, n * sizeof(Triangle*));
        MALLOC(triNeis, int*, n * 3 * sizeof(int));
        
        k = 0;
        for (i = 0; i < n; i++) {
            copyTo(lid, buf, 7 * sizeof(int), &curSize);
            tri = addTriangle(newPoints[lid[0]], newPoints[lid[1]], newPoints[lid[2]], 
                              lid[3], false);
            newTris[i] = tri;
            for (j = 0; j < 3; j++)
                triNeis[k++] = lid[4 + j];
        }
        
        k = 0;
        for (i = 0; i < n; i++) {
            tri = newTris[i];
            for (j = 0; j < 3; j++) {
                idx = triNeis[k++];
                if (idx >= 0)
                    tri->neis[j] = newTris[idx];
            }
        }
                        
        FREE(triNeis);
        FREE(newTris);
        
        // ** holes **
        
        copyTo(&n, buf, sizeof(int), &curSize);
        
        for (i = 0; i < n; i++) {
            copyTo(xy, buf, 2 * sizeof(double), &curSize);
            PUSH_BACK(holes, HoleType(xy[0], xy[1]));
        }

        // ** **
                            
        FOR (pI, newPoints)
            addPoint (*pI);
        
        size = curSize;
          
    }
    
    //----------------------------------------------------------------------------------------------
    int packedSize ()
    {
        return   sizeof (int) + sizeof (bool) + sizeof (double)
               + sizeof (int) + points.size () * (sizeof(int) + 2 * sizeof(double) + 
                                                  sizeof(Rational2))
               + sizeof (int) + gidToSegment.size () * (5 * sizeof (int))
               + sizeof (int) + numTris * 7 * sizeof (int)
               + sizeof (int) + holes.size() * 2 * sizeof(double);
    }
    
    //----------------------------------------------------------------------------------------------
    void copyFrom (void* src, char* buf, int size, int * curSize)
    {
        memcpy (&(buf[*curSize]), src, size); 
        *curSize += size;
    }
        
    //----------------------------------------------------------------------------------------------
    int pack (char* buf, int tgt)
    {
        int                                 i, idx, n, curSize = 0;
        deque <Triangle*> :: iterator       tI;
        deque <Point*> :: iterator          pI;
        map <int, Segment*> :: iterator     gsI;
        HoleContainer :: iterator           hI;
        Triangle                            *tri, *nei;
        Segment                             *segment;
        map <Triangle*, int>                triToNumber;
        
        // ** gid **
        
        copyFrom (&gid, buf, sizeof(int), &curSize);
        
        // ** refinementFinished **
        
        copyFrom (&refinementFinished, buf, sizeof(bool), &curSize);
        
        // ** area **
        
        copyFrom (&area, buf, sizeof(double), &curSize);
        
        // ** points **
        
        n = points.size ();
        copyFrom (&n, buf, sizeof(int), &curSize);
        
        i = 0;
        FOR (pI, points) {
            copyFrom (&((*pI)->gid), buf, sizeof(int), &curSize);
            
            copyFrom ((*pI)->coord, buf, 2 * sizeof(double), &curSize);
            
            curSize += (*pI)->frac.pack (&(buf[curSize]));
            
            (*pI)->lid = i++;
        }
        
        // ** segments **
        
        n = gidToSegment.size ();
        copyFrom (&n, buf, sizeof(int), &curSize);
        
        FOR (gsI, gidToSegment) {
            segment = (*gsI).second;
            copyFrom (&(segment->point0->lid), buf, sizeof(int), &curSize);
            copyFrom (&(segment->point1->lid), buf, sizeof(int), &curSize);
            copyFrom (&(segment->gid), buf, sizeof(int), &curSize);
            copyFrom (&(segment->rid0), buf, sizeof(int), &curSize);
            copyFrom (&(segment->rid1), buf, sizeof(int), &curSize);
        }
            
        // ** triangles **
        
        copyFrom (&numTris, buf, sizeof(int), &curSize);
        
        n = 0;
        FOR (tI, tris) {
            tri = *tI;
            if (! tri->isDeleted ()) {
                triToNumber[tri] = n++;
            }
        }
        
        n = 0;    
        FOR (tI, tris) {
            tri = *tI;
            if (! tri->isDeleted ()) {
                n++;
                for (i = 0; i < 3; i++)
                    copyFrom (&(tri->points[i]->lid), buf, sizeof(int), &curSize);
                copyFrom (&(tri->qIdx), buf, sizeof(int), &curSize);
                for (i = 0; i < 3; i++) {
                    nei = tri->neis[i];
                    idx = ((nei) ? (triToNumber[nei]) : (-1));
                    copyFrom (&idx, buf, sizeof(int), &curSize);
                }
            }
        }
        assert(n == numTris);
        
        // ** holes **
        
        n = holes.size ();
        copyFrom (&n, buf, sizeof(int), &curSize);
        
        FOR (hI, holes) {
            copyFrom (&((*hI).first),  buf, sizeof(double), &curSize);
            copyFrom (&((*hI).second), buf, sizeof(double), &curSize);
        }
        
        assert(curSize == packedSize());
        
        return curSize;
    }
    
    //----------------------------------------------------------------------------------------------
    void Region::findTriNeis ()
    {
        int                                             i;
        REGION::PointContainer::iterator                pI;
        REGION::TriangleContainer::iterator             tI;
        Point                                           *p0, *p1, *point;
        set <Triangle*>                                 res, *set0, *set1;
        Triangle*                                       tri;
        map < Point*, set <Triangle*> * >               ptToTris;
        map < Point*, set <Triangle*> * > :: iterator   ptI;
        set <Triangle*>                                 *triSet;
        pair <set <Triangle*> :: iterator, bool>        insertResult;
                
        FOR (tI, tris) {
            tri = *tI;
            if (!(tri->isDeleted()))
                for (i = 0; i < 3; i++) {
                    point = tri->points[i];
                    ptI = ptToTris.find(point);
                    if (ptI == ptToTris.end()) {
                        NEW(triSet, set <Triangle*>);
                        ptToTris[point] = triSet;
                    }
                    else
                        triSet = (*ptI).second;
                    insertResult = triSet->insert(tri);
                    assert(insertResult.second);
                }
        }
        
        FOR (tI, tris) {
            tri = *tI;
            if (!(tri->isDeleted()))
                for (i = 0; i < 3; i++) {
                    p0 = tri->points[i];
                    p1 = tri->points[NEXT3[i]];
                    set0 = ptToTris[p0];
                    set1 = ptToTris[p1];
                    res.clear();
                    insert_iterator < set <Triangle*> > ii (res, res.begin());
                    set_intersection(set0->begin(), set0->end(), set1->begin(), set1->end(), ii);
                    if (res.size() == 2) {
                        res.erase(tri);
                        tri->neis[PREV3[i]] = *(res.begin());
                    }
                    assert(res.size() == 1);
                }
        }
        
        FOR (ptI, ptToTris) 
            delete (*ptI).second;
    }

    //----------------------------------------------------------------------------------------------
    void writePoly (char* customPolyDir, char* prefix)
    {
        int                                             i, n;
        char                                            fname[255];
        char                                            *myPolyDir;
        FILE*                                           F;
        PointContainer :: iterator                      pI;
        TriangleContainer :: iterator                   tI;
        map <int, Segment*> :: iterator                 gsI;
        Triangle                                        *tri;
        Segment                                         *segment;
        
        myPolyDir = ((customPolyDir == 0) ? (polyDir) : (customPolyDir));

        if (myPolyDir) {
            
            sprintf(fname, "%s/%s%d.node", myPolyDir, prefix, gid);
            F = fopen (fname, "w");
            assert(F);
            fprintf (F, "%d 2 0 1\n", points.size ());
            i = 0;
            FOR (pI, points) {
                fprintf (F, "%d %lf %lf %d\n", ++i, (*pI)->coord[0], (*pI)->coord[1], (*pI)->gid);
                (*pI)->lid = i;
            }
            fclose (F);
            
            sprintf (fname, "%s/%s%d.ele", myPolyDir, prefix, gid);
            F = fopen (fname, "w");
            fprintf (F, "%d 3 0 \n", numTris);
            i = 0;
            FOR (tI, tris) {
                tri = *tI;
                if (! tri->isDeleted ())
                    fprintf (F, "%d %d %d %d\n", ++i, tri->points[0]->lid, tri->points[1]->lid,
                             tri->points[2]->lid);
            }
            fclose (F);
            
            sprintf (fname, "%s/%s%d.poly", myPolyDir, prefix, gid);
            F = fopen (fname, "w");
            FOR (gsI, gidToSegment) {
                segment = (*gsI).second;
                segment->point0->lid = -1;
                segment->point1->lid = -1;
            }
            n = 0;
            FOR (gsI, gidToSegment) {
                segment = (*gsI).second;
                if (segment->point0->lid < 0)
                    segment->point0->lid = ++n;
                if (segment->point1->lid < 0)
                    segment->point1->lid = ++n;
            }
            FOR (gsI, gidToSegment) {
                segment = (*gsI).second;
                segment->point0->lid = -1;
                segment->point1->lid = -1;
            }
            fprintf (F, "%d 2 0 0\n", n);
            n = 0;
            FOR (gsI, gidToSegment) {
                segment = (*gsI).second;
                if (segment->point0->lid < 0) {
                    segment->point0->lid = ++n;
                    fprintf (F, "%d %lf %lf\n", n, segment->point0->coord[0], 
                                                   segment->point0->coord[1]);
                }
                if (segment->point1->lid < 0) {
                    segment->point1->lid = ++n;
                    fprintf (F, "%d %lf %lf\n", n, segment->point1->coord[0], 
                                                   segment->point1->coord[1]);
                }
            }
            fprintf (F, "%d 0\n", gidToSegment.size());
            i = 0;
            FOR (gsI, gidToSegment) 
                fprintf (F, "%d %d %d\n", ++i, (*gsI).second->point0->lid, 
                                               (*gsI).second->point1->lid);
            fprintf (F, "%d\n", holes.size());
            for (i = 0; i < (int)holes.size(); i++)
                fprintf (F, "%d %lf %lf\n", i + 1, holes[i].first, holes[i].second);
            fclose (F);
    
        }
    }
 
    //----------------------------------------------------------------------------------------------
    void writePoly ()
    {
        writePoly (0, "");
    }
    
    //----------------------------------------------------------------------------------------------
    int numberOfTriangles () const
    {
        return numTris;
    } 
    
    //----------------------------------------------------------------------------------------------
    int numberOfPoints () const
    {
        return points.size();
    } 
    
    //----------------------------------------------------------------------------------------------
    int numberOfSegments () const
    {
        return gidToSegment.size();
    } 
    
    //----------------------------------------------------------------------------------------------
    int numberOfPointsToMark ()
    {
        GidSegmentMap::iterator     gsI;
        Segment                     *segment;
        int                         myGid, neiGid, n = 0;
        
        myGid = getGid();
        FOR(gsI, gidToSegment) {
            segment = (*gsI).second;
            neiGid = segment->otherRegion(myGid);
            if ((neiGid < 0) || (myGid < neiGid))
                n += segment->points.size();
        }
        return n;
    } 
};

//==================================================================================================
int packRegion (char* buf, REGION* reg, int tgt)
{
    map <int, REGION*> :: iterator          grI;
    int                                     size;
    
    grI = gidToRegion.find (reg->getGid());
    assert (grI != gidToRegion.end ());
    
    reg->cleanAllQueues ();
    size = reg->pack (buf, tgt);
    gidToRegion.erase (grI);
    regionToProc [reg->getGid()] = tgt;
    delete reg;
    
    return size;
}

//==================================================================================================
REGION* unPackRegion (char* buf, int& size, int src)
{
    REGION* reg;
    NEW (reg, REGION (buf, size, src));
    gidToRegion [reg->getGid()] = reg;
    regionToProc [reg->getGid()] = myProc;
    return reg;
}

//==================================================================================================
int packRegionsOneProc (char* buf, vector <REGION*> part, int tgt)
{
    vector <REGION*> :: iterator            rI;
    int                                     curSize = MSG_HEADER_SIZE;
    int                                     msg[MSG_HEADER_INTS];
    
    msg[0] = MSG_OBJ_DATA;
    msg[1] = part.size ();
    memcpy (buf, &msg, MSG_HEADER_SIZE);
    FOR (rI, part) {
        curSize += packRegion (&(buf[curSize]), *rI, tgt);
    }
    
    return curSize;
}

//==================================================================================================
int unPackRegionsOneProc (int numRegions, char* buf, int src)
{
    int                                     i, size, curSize = MSG_HEADER_SIZE;
    
    for (i = 0; i < numRegions; i++) {
        unPackRegion (&(buf[curSize]), size, src);
        curSize += size;
    }
    
    return curSize;
}

//==================================================================================================
int packRegionsAllProc (char* buf, vector < vector <REGION*> *> parts)
{
    vector < vector <REGION*> *> :: iterator    pI;
    int                                         i = 0, curSize = 0;
    
    assert ((int) parts.size () == numProcs);
    FOR (pI, parts) {
        curSize += packRegionsOneProc (&(buf[curSize]), **pI, i++);
    }
    
    return curSize;
}

//==================================================================================================
void processRegion (REGION* reg)
{
    while (reg->refine(pollFrequency)) {
#ifndef SEQUENTIAL
        mpiPoll ();
#endif
    }
        
    if (!reg->isRefinementFinished()) {
        reg->setRefinementFinished(true);
    }
}


//*********************** COMMUNICATION ROUTINES ***************************************************

#ifndef SEQUENTIAL

//==================================================================================================
void initMPI (int argc, char* argv[])
{
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank (MPI_COMM_WORLD, &myProc);
}

//==================================================================================================
void finalizeMPI ()
{
    MPI_Finalize();		
}

//==================================================================================================
void barrier ()
{
    BEGIN_TIMING (TIME_SYNC);
    MPI_Barrier (MPI_COMM_WORLD);
    END_TIMING (TIME_SYNC);
}

//==================================================================================================
void bcast (void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
    BEGIN_TIMING (TIME_COMM);
    MPI_Bcast (buffer, count, datatype, root, comm);
    END_TIMING (TIME_COMM);
}

//==================================================================================================
int scan (void* sendbuf,void* recvbuf,int count,
          MPI_Datatype datatype,MPI_Op op,MPI_Comm comm)
{
    int result;
    
    BEGIN_TIMING (TIME_COMM);
    result = MPI_Scan (sendbuf, recvbuf, count, datatype, op, comm);
    END_TIMING (TIME_COMM);
    
    return result;
}

//==================================================================================================
void scatter (void* sendbuf, int sendcount, MPI_Datatype sendtype,
              void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    BEGIN_TIMING (TIME_COMM);
    MPI_Scatter (sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    END_TIMING (TIME_COMM);
}
        
//==================================================================================================
void scatterv (void* sendbuf, int* sendcounts, int* displs, MPI_Datatype sendtype,
               void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    BEGIN_TIMING (TIME_COMM);
    MPI_Scatterv (sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
    END_TIMING (TIME_COMM);
}

//==================================================================================================
void allReduce (void* sendbuf, void* recvbuf, int count, 
                 MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    BEGIN_TIMING (TIME_COMM);
    MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);
    END_TIMING (TIME_COMM);
}


//==================================================================================================
void allGather (void* sendbuf, int sendcount, MPI_Datatype sendtype,
                 void* recvbuf, int recvcount, MPI_Datatype recvtype,
                 MPI_Comm comm)
{
    BEGIN_TIMING (TIME_COMM);
    MPI_Allgather (sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
    END_TIMING (TIME_COMM);
}

//==================================================================================================
void allGatherV (void* sendbuf, int sendcount, 
                 MPI_Datatype sendtype,
                 void* recvbuf, int* recvcounts, int* displs,
                 MPI_Datatype recvtype, MPI_Comm comm)
{
    BEGIN_TIMING (TIME_COMM);
    MPI_Allgatherv (sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
    END_TIMING (TIME_COMM);
}

//==================================================================================================
void allToAll (void* sendbuf, int sendcount, MPI_Datatype sendtype,
               void* recvbuf, int recvcount, MPI_Datatype recvtype,
               MPI_Comm comm)
{
    BEGIN_TIMING (TIME_COMM);
    MPI_Alltoall (sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
    END_TIMING (TIME_COMM);
}

//==================================================================================================
void allToAllV (void* sendbuf, int* sendcounts, int* sdispls, 
                MPI_Datatype sendtype,
                void* recvbuf, int* recvcounts, int* rdispls,
                MPI_Datatype recvtype, MPI_Comm comm)
{
    BEGIN_TIMING (TIME_COMM);
    MPI_Alltoallv (sendbuf, sendcounts, sdispls, sendtype,
                   recvbuf, recvcounts, rdispls, recvtype, comm);
    END_TIMING (TIME_COMM);
}

//==================================================================================================
void mpiPoll ()
{
    int         flag, src, len, tag;
    MPI_Status  status;
    char*       buf;

    while (true) {
        
        BEGIN_TIMING (TIME_MPI_PROBE);
        MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        END_TIMING (TIME_MPI_PROBE);
        
        if (flag) {
            
            tag = status.MPI_TAG;
            src = status.MPI_SOURCE;
            
            BEGIN_TIMING (TIME_COMM);
            MPI_Get_count (&status, MPI_BYTE, &len);
            END_TIMING (TIME_COMM);
            
            MALLOC (buf, char*, len);
            
            BEGIN_TIMING (TIME_COMM);
            MPI_Recv (buf, len, MPI_BYTE, src, tag, MPI_COMM_WORLD, &status);
            END_TIMING (TIME_COMM);
            
            msgHandler (src, buf, len);
            FREE (buf);
        }
        else
            break;
    };
}

#endif

//==================================================================================================
void terminationProbe (int msgCountReceived, int colorReceived)
{
    int                                 msgSend[MSG_HEADER_INTS];
    
    if (myProc == 0) {
        if ((colorReceived == 0) && (!stateActive) &&
            (msgCountReceived + numSplitSent - numSplitRecv == 0)) {
            msgSend[0] = MSG_TERMINATE;                             // termination
        }
        else {      
                                                        // EWD Rule 1 and Rule 5
            msgSend[0] = MSG_FINISHED;
            msgSend[1] = 0;
            msgSend[2] = 0;                                         // EWD Rule 6
            myColor = 0;                                            // EWD Rule 6
        }
    }
    else {
            msgSend[0] = MSG_FINISHED;
            msgSend[1] = msgCountReceived + numSplitSent - numSplitRecv;   // EWD Rule 2
            if (myColor == 1) {
                msgSend[2] = 1;                                     // EWD Rule 4
                myColor = 0;                                        // EWD Rule 7
            }
            else
                msgSend[2] = colorReceived;                         // EWD Rule 4
    }
    request((myProc - 1 + numProcs) % numProcs, (char*)&msgSend, MSG_HEADER_SIZE);
}

//==================================================================================================
void msgHandler (int src, void* buf, int size)
{
    int                                 s;
    int                                 msgRecv[MSG_HEADER_INTS], msgSend[MSG_HEADER_INTS];
    map <int, REGION*> :: iterator      grI, grJ;
    map <int, Segment*>  :: iterator    gsI;
    map <int, int> :: iterator          rpI;
    
    memcpy (&msgRecv, buf, MSG_HEADER_SIZE);
    
    switch (msgRecv[0]) {
        
        case MSG_SPLIT:
            grI = gidToRegion.find(msgRecv[1]);
            assert(grI != gidToRegion.end());
            if (grI != gidToRegion.end()) {
                processingSplitMsg = true;
                myColor = 1;                            // EWD Rule 3
                numSplitRecv++;                         // EWD Rule 0
                if (msgRecv[3] != src) {
                    msgSend[0] = MSG_UPDATE;
                    msgSend[1] = msgRecv[1];
                    request(msgRecv[3], (char*)&msgSend, MSG_HEADER_SIZE);
                }
                (*grI).second->splitSegmentsDst(msgRecv[2], &(((char*)buf)[MSG_HEADER_SIZE]));
                processingSplitMsg = false;
                if (ownFinishedToken) {
                    ownFinishedToken = false;
                    terminationProbe(0, 0);             // EWD Rule 2
                }
            }
            else {
                rpI = regionToProc.find(msgRecv[1]);
                assert (rpI != regionToProc.end ());
                request((*rpI).second, (char*)buf, size);
            }
            break;
            
        case MSG_UPDATE:
            regionToProc[msgRecv[1]] = src;
            break;
            
        case MSG_OBJ_DATA:
            s = unPackRegionsOneProc(msgRecv[1], (char*)buf, src);
            assert (s == size);
            break;
            
        case MSG_FINISHED:
            if ((processingSplitMsg) || (stateActive))
                ownFinishedToken = true;
            else
                terminationProbe(msgRecv[1], msgRecv[2]);
            break;
            
        case MSG_TERMINATE:
            allFinished = true;
            if (myProc != 0) {
                msgSend[0] = MSG_TERMINATE;
                request((myProc - 1 + numProcs) % numProcs, (char*)&msgSend, MSG_HEADER_SIZE);
            }
            break;
        
        case MSG_POINT_NUM:
            grI = gidToRegion.find(msgRecv[1]);
            assert(grI != gidToRegion.end());
            gsI = (*grI).second->gidToSegment.find(msgRecv[2]);
            assert(gsI != (*grI).second->gidToSegment.end());
            (*gsI).second->assignPointNumbers(&(msgRecv[3]));
            numPointNumReceived++;
            break;
            
        default:
            assert(0);
    }
}

//==================================================================================================
void request (int tgt, char* buffer, int bytes)
{
#ifndef SEQUENTIAL
    mpiPoll();
    BEGIN_TIMING(TIME_COMM);
    MPI_Send(buffer, bytes, MPI_BYTE, tgt, 0, MPI_COMM_WORLD);
    END_TIMING(TIME_COMM);
#else
    msgHandler(myProc, buffer, bytes);
#endif
}

//==================================================================================================
void computeMapping (vector < vector <REGION*> *> & partitions)
{
    if (consequtiveDistribution) {
        int                                     i, n = gidToRegion.size();
        int                                     curCluster = -1, clusterAvailProcs = 0;
        double                                  totalWork = 0;
        map <int, REGION*> :: iterator          grI;
        double                                  procWeight, cumProcWeightSought = 0;
        int                                     cumProcWeightAssigned = 0;
        
        for (i = 0; (unsigned)i < procWork.size(); i++)
            totalWork += procCount[i] * procWork[i];
        
        grI = gidToRegion.begin();
        for (i = 0; i < numProcs; i++) {
            while ((clusterAvailProcs == 0) && (curCluster < (int) procCount.size())) {
                curCluster++;
                clusterAvailProcs = procCount[curCluster];
            }
            assert((unsigned)curCluster < procCount.size());
            clusterAvailProcs--;
                
            procWeight = procWork[curCluster] / totalWork * n;
            cumProcWeightSought += procWeight;
            while ((cumProcWeightAssigned <= cumProcWeightSought) && (grI != gidToRegion.end())) {
                PUSH_BACK(*(partitions[i]), ((*grI).second));
                cumProcWeightAssigned++;
                grI++;
            }
        }
        assert(grI == gidToRegion.end());
    }
    else {
    
        map <int, REGION*> :: iterator                  grI;
        REGION :: GidSegmentMap :: iterator             gsI;
        REGION                                          *region;
        Segment                                         *segment;
        int                                             i, i0, i1, rid0, rid1, lid0, lid1, w;
        int                                             *xadj, *adjncy, *vwgt, *adjwgt, *iadj;
        map <OrderedPair<int>, double>                  neisToWgt;
        map <OrderedPair<int>, double> :: iterator      nwJ;
        map <int, double>                               neiToWgt;
        map <int, double> :: iterator                   nwI;
        map <int, int>                                  gidToLid;
        vector <REGION*>                                lidToRegion;
        double                                          maxVWgt = 0, maxAdjWgt = 0, vCoef, adjCoef;
        int                                             numVerts, numEdges;
        int                                             wgtflag = 3, numflag = 0, nparts = numProcs;
        int                                             options[5] = {0, 0, 0, 0, 0};
        int                                             edgecut;
        int                                             *part;
        
        
        i = 0;
        FOR (grI, gidToRegion) {
            gidToLid [(*grI).first] = i++;
            PUSH_BACK (lidToRegion, (*grI).second);
        }
 
        numVerts = gidToRegion.size ();
        MALLOC (xadj, int*, (numVerts + 1) * sizeof (int));
 
        // ** find neisToWgt, xadj, maxVWgt, maxAdjWgt **
        xadj [0] = 0;
        i = 0;
        FOR (grI, gidToRegion) {
            rid0 = (*grI).first;
            region = (*grI).second;
            region->computeArea();
            neiToWgt.clear ();
            FOR (gsI, region->getGidSegmentMap()) {
                segment = (*gsI).second;
                rid1 = segment->otherRegion (rid0);
                if (rid1 >= 0) {
                    nwI = neiToWgt.find (rid1);
                    if (nwI == neiToWgt.end ())
                        neiToWgt [rid1] = segment->length ();
                    else
                        (*nwI).second += segment->length ();
                }
            }
            i++;
            xadj [i] = xadj [i - 1] + neiToWgt.size ();
            if (region->getArea() > maxVWgt)
                maxVWgt = region->getArea();
            FOR (nwI, neiToWgt) {
                OrderedPair<int> ip (rid0, (*nwI).first);
                if (neisToWgt.find (ip) == neisToWgt.end ()) {
                    neisToWgt [ip] = (*nwI).second;
                    if ((*nwI).second > maxAdjWgt)
                        maxAdjWgt = (*nwI).second;
                }
            }
        }
        numEdges = neisToWgt.size ();
        vCoef   = ((double)1e7) / maxVWgt / numVerts;
        adjCoef = ((double)1e7) / maxAdjWgt / numEdges;
 
        MALLOC (adjncy, int*, 2 * numEdges * sizeof (int));
        MALLOC (vwgt, int*, numVerts * sizeof (int));
        MALLOC (adjwgt, int*, 2 * numEdges * sizeof (int));
        CALLOC (iadj, int*, (numVerts + 1) * sizeof (int));
        CALLOC (part, int*, numVerts * sizeof (int));
 
        // ** find vwgt **
        FOR (grI, gidToRegion)
            vwgt [gidToLid [(*grI).first]] = (int) (vCoef * (*grI).second->getArea());
 
        // ** find adjncy, adjwgt **
        FOR (nwJ, neisToWgt) {
            lid0 = gidToLid [(*nwJ).first.n0];
            lid1 = gidToLid [(*nwJ).first.n1];
            w = (int) (adjCoef * (*nwJ).second);
            i0 = xadj [lid0] + (iadj [lid0] ++);
            i1 = xadj [lid1] + (iadj [lid1] ++);
            adjncy [i0] = lid1;
            adjncy [i1] = lid0;
            adjwgt [i0] = w;
            adjwgt [i1] = w;
        }
 
        BEGIN_TIMING (TIME_METIS);
        if (nparts <= 8)
            METIS_PartGraphRecursive (&numVerts, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                                      &nparts, options, &edgecut, part);
        else
            METIS_PartGraphKway (&numVerts, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                                 &nparts, options, &edgecut, part);
        END_TIMING (TIME_METIS);
        
        for (i = 0; i < numVerts; i++) {
            PUSH_BACK(*(partitions[part[i]]), lidToRegion[i]);
        }
        
        FREE (xadj);
        FREE (adjncy);
        FREE (vwgt);
        FREE (adjwgt);
        FREE (iadj);
        FREE (part);
    }
}

//==================================================================================================
void parseCmdLine (int argc, char* argv[])
{
    int                                         i = 3, j;
    double                                      f;
    
    if (argc <= 4) {
        printHelp();
	    exit(0);
    }
    
    while (i < argc) {
        if (strcmp (argv[i], "-a") == 0) {
            areaBound = atof (argv[i + 1]);
            i += 2;
        }
        if (strcmp (argv[i], "-q") == 0) {
            angleBound = RADIANS(atof(argv[i + 1]));
            cosAngleBound = cos(angleBound);
            i += 2;
        } else
        if (strcmp (argv[i], "-m") == 0) {
            msgAggregation = atoi (argv[i + 1]);
            i += 2;
        } else
        if (strcmp (argv[i], "-o") == 0) {
            polyDir = argv[i + 1];
            i += 2;
        } else
        if (strcmp (argv[i], "-s") == 0) {
            statDir = argv[i + 1];
            i += 2;
        } else
        if (strcmp (argv[i], "-Q") == 0) {
            computeHistogram = true;
            i++;
        } else
        if (strcmp (argv[i], "-c") == 0) {
            consequtiveDistribution = true;
            i++;
        } else
        if (strcmp (argv[i], "-f") == 0) {
            pollFrequency = atoi(argv[i + 1]);
            i += 2;
        } else 
	    if ((strcmp (argv[i], "-h") == 0) || (strcmp (argv[i], "--help") == 0)){
	        printHelp();
	        exit(0);
	    } else 
	    if (strcmp (argv[i], "-W") == 0){
            consequtiveDistribution = true;
	        j = i + 1;
            while (j < argc) {
                PUSH_BACK(procCount, atoi(argv[j++]));
                sscanf(argv[j++], "%lf", &f);
                PUSH_BACK(procWork, f);
            }
	        i = j;
	    } else 
            assert(0);
    }
    assert(areaBound > 0);
    inputRegionFileName = argv[1];
    inputTriangleFileName = argv[2];
    assert(inputRegionFileName);
    assert(inputTriangleFileName);
}

//==================================================================================================
void loadData (
#ifdef PCDM_LIBRARY
               PcdmIO * pcdmio
#endif
              )
{
    int                                         i, j, d;
    Point                                       *point, *pt[3];
    vector <Point*>                             pointVector;
    vector <REGION*>                            regionVector;
    vector <Point*> :: iterator                 pI;
    vector < vector <REGION*> *>                part;
    vector < vector <REGION*> *> :: iterator    prI;
    vector <REGION*> :: iterator                rI;
    Pool                                        pointStorage(sizeof(Point));
    map <int, int> :: iterator                  rpI;
    REGION                                      *region;
    Segment                                     *segment;
    map <int, REGION*> :: iterator              grI;
    int                                         idxIncrement, idIncrement;
#ifndef PCDM_LIBRARY
    PcdmIO                                      *pcdmio;
#endif
#ifndef SEQUENTIAL
    int                                         recvcount, sendcount = 0;
    int                                         *sendcounts = 0, *displs = 0;
    char                                        *sendbuf = 0, *recvbuf = 0;
    int                                         s, n, size, *regionMap;
    vector <REGION*>                            *regionVectorPtr;
    int                                         k;
#endif
    double                                      *thisRegionPointList;
    int                                         *thisRegionLocalSegmentList;
    int                                         *thisRegionGlobalSegmentList;
    int                                         *thisRegionTriangleList;
    int                                         *thisRegionPointIdList;
    double                                      *thisRegionHoleList;
    
    if (myProc == 0) {
        
        BEGIN_TIMING (TIME_IO_READ);
        
#ifndef PCDM_LIBRARY
        NEW(pcdmio, PcdmIO);
        pcdmio->oneFlag = 0;
        pcdmio->read(inputRegionFileName, inputTriangleFileName);
#endif
        idxIncrement = ((pcdmio->oneFlag == 0) ? (0) : (-1));
        idIncrement = ((pcdmio->oneFlag == 0) ? (1) : (0));
        originalNumPoints  = pcdmio->numberOfPoints;
        originalNumRegions = pcdmio->numberOfRegions;
        
        pcdmio->pullApart();
        
        for (i = 0; i < pcdmio->numberOfRegions; i++) {
            NEW(region, REGION(i + 1));
            gidToRegion[region->getGid()] = region;
            PUSH_BACK(regionVector, region);
            
            thisRegionPointList = pcdmio->regionLocalPointLists[i];
            thisRegionPointIdList = pcdmio->regionLocalPointIdLists[i];
            for (j = 0; j < pcdmio->regionNumberOfPoints[i]; j++) {
                point = (Point*) region->pointStorage.allocate();
                point->init(&(thisRegionPointList[j * 2]), thisRegionPointIdList[j] + idIncrement);
                region->addPoint(point);
                PUSH_BACK(pointVector, point);
            }
            
            thisRegionLocalSegmentList = pcdmio->regionLocalSegmentLists[i];
            thisRegionGlobalSegmentList = pcdmio->regionGlobalSegmentLists[i];
            for (j = 0; j < pcdmio->regionNumberOfSegments[i]; j++) {
                for (d = 0; d < 2; d++) {
                    pt[d] = pointVector[thisRegionLocalSegmentList[j * 2 + d] + idxIncrement];
                }
                s = thisRegionGlobalSegmentList[j] + idxIncrement;
                NEW(segment, Segment(s + 1, pt[0], pt[1]));
                region->addSegment(segment);
                segment->addRegion(pcdmio->segmentToRegions[s * 2 + 0] + idIncrement);
                segment->addRegion(pcdmio->segmentToRegions[s * 2 + 1] + idIncrement);
            }
        
            thisRegionTriangleList = pcdmio->regionLocalTriangleLists[i];// !!!
            for (j = 0; j < pcdmio->regionNumberOfTriangles[i]; j++) {
                for (d = 0; d < 3; d++) {
                    pt[d] = pointVector[thisRegionTriangleList[j * 3 + d] + idxIncrement];
                }
                region->addTriangle(pt[0], pt[1], pt[2], 0, true);
            }
            
            thisRegionHoleList = pcdmio->regionHoleLists[i];
            for (j = 0; j < pcdmio->regionNumberOfHoles[i]; j++) {
                PUSH_BACK(region->holes, REGION::HoleType(thisRegionHoleList[j * 2 + 0],
                                                          thisRegionHoleList[j * 2 + 1]));
            }
            
            region->findTriNeis();
            pointVector.clear();
        }
        
#ifndef PCDM_LIBRARY
        delete pcdmio;
#endif
        
        END_TIMING (TIME_IO_READ);
    }
    
    if (numProcs == 1) 
        return;
    
#ifndef SEQUENTIAL
    
    bcast(&originalNumPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
    bcast(&originalNumRegions, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    MALLOC(regionMap, int*, (originalNumRegions + numProcs) * sizeof(int));
    
    if (myProc == 0) {
        
        // ** compute the mapping **
        
        for (i = 0; i < numProcs; i++) {
            NEW(regionVectorPtr, vector <REGION*>);
            PUSH_BACK(part, regionVectorPtr);
        }
                
        computeMapping(part);
        
        MALLOC(sendcounts, int*, numProcs * sizeof (int));
        for (i = 0; i < numProcs; i++)
            sendcounts [i] = MSG_HEADER_SIZE;
        sendcount = numProcs * MSG_HEADER_SIZE;
        k = 0;
        for (i = 0; i < numProcs; i++) {
            regionMap[k++] = part[i]->size();
            FOR (rI, *(part[i])) {
                region = (*rI);
                regionMap[k++] = region->gid;
                s = region->packedSize ();
                sendcounts[i] += s;
                sendcount += s;
            }
        }
    }
    
    bcast(regionMap, originalNumRegions + numProcs, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (myProc != 0) {
        k = 0;
        for (i = 0; i < numProcs; i++) {
            n = regionMap[k++];
            for (j = 0; j < n; j++)
                regionToProc[regionMap[k++]] = i;
        }
    }
    
    FREE(regionMap);
    
    if (myProc == 0) {
        
        MALLOC (displs, int*, numProcs * sizeof (int));
        displs [0] = 0;
        for (i = 1; i < numProcs; i++)
            displs [i] = displs [i-1] + sendcounts [i-1];
        
        MALLOC (sendbuf, char*, sendcount);
        size = packRegionsAllProc (sendbuf, part);
        assert((size == sendcount));
        
        FOR (prI, part)
            delete *prI;
    }

    numInitBytes += sizeof(int) + sendcount;
    
    scatter(sendcounts, 1, MPI_INT, &recvcount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    MALLOC(recvbuf, char*, recvcount);
    
    scatterv(sendbuf, sendcounts, displs, MPI_BYTE, 
             recvbuf, recvcount, MPI_BYTE, 0, MPI_COMM_WORLD);
    
    if (myProc == 0) {
        FREE(sendbuf);
        FREE(displs);
        FREE(sendcounts);
    }
    
    msgHandler(0, recvbuf, recvcount);
    
    FREE(recvbuf);
    
#endif
}

//==================================================================================================
#define NUM_STATS 23

void statistics (bool elementwise)
{
    int                                         i, j, arg[NUM_STATS], res[NUM_STATS];
    double                                      angles[3], minAngle;
    deque <Triangle*> :: iterator               tI;
    FILE*                                       statFile;
    char                                        statFileName[255];
    time_t                                      curTime;
    map <int, REGION*> :: iterator              grI;
    REGION                                      *region;
    
    if (statDir) {
        
        sprintf (statFileName, "%s/%d", statDir, myProc);
        if (elementwise)
            statFile = fopen (statFileName, "w");
        else
            statFile = fopen (statFileName, "a");
        assert (statFile != 0);
        
        if (elementwise) {
            
            time (&curTime);
            fprintf(statFile, "================================================================\n");
            fprintf(statFile, "# %s", ctime(&curTime));
            fprintf(statFile, "# generated with ");
            for (i = 0; i < argC; i++)
                fprintf(statFile, " %s", argV[i]);
            fprintf(statFile, "\n");
            fprintf(statFile, "----------------------------------------------------------------\n");
        
            minAngle = PI;
            for (i = 0; i < NUM_STATS; i++)
                arg[i] = 0;
            
            FOR (grI, gidToRegion) {
                region = (*grI).second;
                arg[18] += region->numberOfTriangles ();
                if (computeHistogram) {
                    FOR (tI, region->getTriangles()) {
                        if (! (*tI)->isDeleted ()) {
                            (*tI)->angles (angles);
                            for (i = 0; i < 3; i++) {
                                j = (int) (angles[i] / RADIANS(10));
                                if (j < 0) j = 0; else
                                if (j > 17) j = 17;
                                arg[j] += 1;
                                if (angles[i] < minAngle)
                                    minAngle = angles[i];
                            }
                        }
                    }
                }
            }
            arg[19] = numSplitSent;
            arg[20] = gidToRegion.size();
            arg[21] = numInitBytes;
            arg[22] = numSplitBytes;
    
#ifndef SEQUENTIAL
            allReduce (arg, res, NUM_STATS, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
            memcpy(res, arg, NUM_STATS * sizeof(int));
#endif
            fprintf(statFile, "Global number of triangles = %d\n", res[18]);
            fprintf(statFile, "Local number of triangles = %d\n",  arg[18]);
            fprintf(statFile, "Global number of subdomains = %d\n",res[20]);
            fprintf(statFile, "Local number of subdomains = %d\n", arg[20]);
                
            fprintf(statFile, "Global number of split messages = %d\n", res[19]);
            fprintf(statFile, "Global size of initialization messages = %d\n", res[21]);
            fprintf(statFile, "Global size of split messages = %d\n", res[22]);
            fprintf(statFile, "----------------------------------------------------------------\n");
            fprintf(statFile, "\tAngle histogram:\n");
            for (i = 0; i < 18; i++)
                fprintf (statFile, "%d -- %d:\t%d\n", i*10, (i+1)*10, res[i]);
            fprintf(statFile, "----------------------------------------------------------------\n");
        }
        else {
            fprintf(statFile, "Time Total = %lf\n",                timer0[TIME_TOTAL]);
            fprintf(statFile, "  Time Metis = %lf\n",              timer0[TIME_METIS]);
            fprintf(statFile, "  Time Read Data = %lf\n",          timer0[TIME_IO_READ]);
            fprintf(statFile, "  Time Refine PCDM = %lf\n",        timer0[TIME_REFINE_PCDM]);
            fprintf(statFile, "  Time Write Data = %lf\n",         timer0[TIME_IO_WRITE]);
            fprintf(statFile, "  Time Communication = %lf\n",      timer0[TIME_COMM]);
            fprintf(statFile, "  Time Synchronization = %lf\n",    timer0[TIME_SYNC]);
            fprintf(statFile, "  Time Probe = %lf\n",              timer0[TIME_MPI_PROBE]);
            fprintf(statFile, "  Time Clean = %lf\n",              timer0[TIME_FREE]);
            fprintf(statFile, "----------------------------------------------------------------\n");
            fprintf(statFile, "  Time Stat = %lf\n" ,              timer0[TIME_STATS]);
        }
        fclose (statFile);
    }
}

//==================================================================================================
void pcdmInit(int argc, char* argv[])
{
    int i;
    
#ifdef MEM_CHECK
    mtrace();
#endif
    
    for (i = 0; i < NUM_TIMES; i++)
        timer0[i] = 0;
    
    BEGIN_TIMING (TIME_TOTAL);
        
    exactinit ();

    argC = argc;
    argV = argv;
        
#ifndef SEQUENTIAL    
    initMPI(argc, argv);
#endif
}

//==================================================================================================
void refineAll ()
{
    map <int, REGION*> :: iterator              grI;
    REGION::GidSegmentMap::iterator             gsI;
    map <int, int> :: iterator                  rpI;
    Segment::PointContainer::iterator           pI;
    REGION::PointContainer::iterator            ptI;
    REGION                                      *region;

    int                                         numPointNumExpected = 0, savedMsgAggregation;
    int                                         locNumPointsToMark = 0, cumNumPointsToMark = 0;
    int                                         nextPointNumberLocal = 0;
    int                                         myGid, neiGid;
    Segment                                     *segment;
    int                                         msg[MSG_HEADER_INTS];

    allFinished = false;
    ownFinishedToken = false;
    stateActive = true;
    myColor = 0;
    processingSplitMsg = false;
    
    FOR (grI, gidToRegion) {
        region = (*grI).second;
        processRegion(region);
    }
    
    savedMsgAggregation = msgAggregation;
    msgAggregation = 1;
    FOR (grI, gidToRegion) 
        (*grI).second->cleanAllQueues();
    
    stateActive = false;
    if ((myProc == 0) || (ownFinishedToken)) {
        ownFinishedToken = false;
        terminationProbe(1, 1);             
    }
    
#ifndef SEQUENTIAL
    while (!allFinished)
        mpiPoll();
#endif

    msgAggregation = savedMsgAggregation;
    
    barrier();
    
    locNumPointsToMark = 0;
    FOR(grI, gidToRegion)
        locNumPointsToMark += (*grI).second->numberOfPointsToMark();
#ifndef SEQUENTIAL
    scan(&locNumPointsToMark, &cumNumPointsToMark, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    numPointNumReceived = 0;
#else
    cumNumPointsToMark = locNumPointsToMark;
#endif
    nextPointNumberLocal = originalNumPoints + cumNumPointsToMark - locNumPointsToMark + 1;

    FOR (grI, gidToRegion) {
        region = (*grI).second;
        myGid = region->getGid();
        FOR(gsI, region->gidToSegment) {
            segment = (*gsI).second;
            if (segment->points.size() > 0) {
                neiGid = segment->otherRegion(myGid);
                assert(neiGid != myGid);
                if (neiGid > 0) {
                    if (myGid < neiGid) {
                        msg[0] = MSG_POINT_NUM;
                        msg[1] = neiGid;
                        msg[2] = (*gsI).first;
                        msg[3] = nextPointNumberLocal;
                        request(regionToProc[neiGid], (char*)&msg, MSG_HEADER_SIZE);
                        segment->assignPointNumbers(&nextPointNumberLocal);
                    }
                    else
                        numPointNumExpected++;
                }
                else
                    segment->assignPointNumbers(&nextPointNumberLocal);
            }
        }
    }
    
#ifndef SEQUENTIAL
    while (numPointNumReceived < numPointNumExpected) {
        mpiPoll();
    }
#endif
}

//==================================================================================================
void pcdmFinalize ()
{
#ifndef SEQUENTIAL    
    finalizeMPI();
#endif
    
    END_TIMING(TIME_TOTAL);
}

//==================================================================================================
//      MAIN
//==================================================================================================
#ifdef PCDM_LIBRARY
void pcdmRefine (struct PcdmIO * pcdmio)
#else
int main (int argc, char* argv[])
#endif
{
    map <int, REGION*> :: iterator                  grI;
    REGION                                          *region;
#ifdef PCDM_LIBRARY
    int                                             i, j, k, d, nr, idxIncrement;
    Triangle                                        *tri;
    Segment                                         *segment;
    Point                                           *point;
    REGION::PointContainer :: iterator              pI;
    REGION::GidSegmentMap :: iterator               gsI;
    REGION::TriangleContainer :: iterator           tI;
    Segment::PointContainer :: iterator             psI;
    double                                          *thisRegionPointList;
    double                                          *thisRegionHoleList;
    int                                             *thisRegionTriangleList;
    int                                             *thisRegionPointIdList;
    int                                             *thisRegionSegmentList;
    int                                             *thisRegionSegmentIdList;
    int                                             *thisSegmentInnerPointList;
    int                                             **thisRegionSegmentInnerPointList;
#endif
        
    try {
#ifdef PCDM_LIBRARY
        assert(pcdmio->areaBound > 0);
        areaBound = pcdmio->areaBound;
        if (pcdmio->angleBound > 0) {
            angleBound = pcdmio->angleBound;
            cosAngleBound = cos(angleBound);
        }
        if (pcdmio->aggregation > 0)
            msgAggregation = pcdmio->aggregation;
        if (pcdmio->consequtiveDistribution > 0)
            consequtiveDistribution = (pcdmio->consequtiveDistribution != 0);
        if (pcdmio->pollFrequency > 0)
            pollFrequency = pcdmio->pollFrequency;
        polyDir = pcdmio->polyDir;
        statDir = pcdmio->statDir;
#else
        pcdmInit(argc, argv);
        parseCmdLine(argc, argv);
#endif
        
        loadData(
#ifdef PCDM_LIBRARY
                 pcdmio
#endif
                );
        
#ifndef SEQUENTIAL
        barrier();
#endif
        
        refineAll();
        
#ifndef SEQUENTIAL
        barrier();
#endif
    
        BEGIN_TIMING(TIME_STATS);    
        statistics(true);
        END_TIMING(TIME_STATS);    
        
#ifdef PCDM_LIBRARY
        idxIncrement = ((pcdmio->oneFlag == 0) ? (0) : (1));
        pcdmio->clear();
        pcdmio->numberOfRegions = gidToRegion.size();
        nr = pcdmio->numberOfRegions;
        MALLOC(pcdmio->regionNumberOfPoints,        int*,     nr * sizeof(int));
        MALLOC(pcdmio->regionLocalPointLists,       double**, nr * sizeof(double*));
        MALLOC(pcdmio->regionLocalPointIdLists,     int**,    nr * sizeof(int*));
         
        MALLOC(pcdmio->regionNumberOfSegments,      int*,     nr * sizeof(int));
        MALLOC(pcdmio->regionLocalSegmentLists,     int**,    nr * sizeof(int*));
        MALLOC(pcdmio->regionLocalSegmentIdLists,   int**,    nr * sizeof(int*));
        MALLOC(pcdmio->regionSegmentInnerPointList, int***,   nr * sizeof(int**));
        
        MALLOC(pcdmio->regionNumberOfTriangles,     int*,     nr * sizeof(int));
        MALLOC(pcdmio->regionLocalTriangleLists,    int**,    nr * sizeof(int*));
        
        MALLOC(pcdmio->regionNumberOfHoles,         int*,     nr * sizeof(int));
        MALLOC(pcdmio->regionHoleLists,             double**, nr * sizeof(double*));
        i = 0;
#endif
        FOR (grI, gidToRegion) {
            
            BEGIN_TIMING(TIME_IO_WRITE);
        
            region = (*grI).second;
#ifdef PCDM_LIBRARY
            
                // ** points **
            pcdmio->regionNumberOfPoints[i] = region->numberOfPoints();
            MALLOC(thisRegionPointList, double*, 
                   pcdmio->regionNumberOfPoints[i] * 2 * sizeof(double));
            pcdmio->regionLocalPointLists[i] = thisRegionPointList;
            MALLOC(thisRegionPointIdList, int*, 
                   pcdmio->regionNumberOfPoints[i] * sizeof(int));
            pcdmio->regionLocalPointIdLists[i] = thisRegionPointIdList;
            j = 0;
            
            FOR (pI, region->points) {
                point = (*pI);
                thisRegionPointList[j * 2 + 0] = point->coord[0];
                thisRegionPointList[j * 2 + 1] = point->coord[1];
                thisRegionPointIdList[j] = point->gid;
                point->lid = j++;
            }
            
                // ** segments **
            pcdmio->regionNumberOfSegments[i] = region->numberOfSegments();
            MALLOC(thisRegionSegmentList, int*, 
                   pcdmio->regionNumberOfSegments[i] * 2 * sizeof(int));
            pcdmio->regionLocalSegmentLists[i] = thisRegionSegmentList;
            MALLOC(thisRegionSegmentIdList, int*, 
                   pcdmio->regionNumberOfSegments[i] * sizeof(int));
            pcdmio->regionLocalSegmentIdLists[i] = thisRegionSegmentIdList;
            MALLOC(thisRegionSegmentInnerPointList, int**, 
                   pcdmio->regionNumberOfSegments[i] * sizeof(int));
            pcdmio->regionSegmentInnerPointList[i] = thisRegionSegmentInnerPointList;
            j = 0;
            
            FOR (gsI, region->gidToSegment) {
                segment = (*gsI).second;
                thisRegionSegmentList[j * 2 + 0] = segment->point0->lid + idxIncrement;
                thisRegionSegmentList[j * 2 + 1] = segment->point1->lid + idxIncrement;
                thisRegionSegmentIdList[j] = segment->gid;
                MALLOC(thisSegmentInnerPointList, int*, segment->points.size() * 3 * sizeof(int));
                thisRegionSegmentInnerPointList[j] = thisSegmentInnerPointList;
                k = 0;
                FOR(psI, segment->points) {
                    point = *psI;
                    thisSegmentInnerPointList[k++] = point->lid + idxIncrement;
                    thisSegmentInnerPointList[k++] = point->frac.numerator();
                    thisSegmentInnerPointList[k++] = point->frac.power();
                }
                j++;
            }
            
                // ** triangles **
            pcdmio->regionNumberOfTriangles[i] = region->numberOfTriangles();
            MALLOC(thisRegionTriangleList, int*, 
                    pcdmio->regionNumberOfTriangles[i] * 3 * sizeof(int));
            pcdmio->regionLocalTriangleLists[i] = thisRegionTriangleList;
            k = 0;
            FOR (tI, region->tris) {
                tri = *tI;
                if (!(tri->isDeleted()))
                    for (d = 0; d < 3; d++)
                        thisRegionTriangleList[k++] = tri->points[d]->lid + idxIncrement;
            }
            
                // ** holes **
            pcdmio->regionNumberOfHoles[i] = region->holes.size();
            MALLOC(thisRegionHoleList, double*, 
                    pcdmio->regionNumberOfHoles[i] * 2 * sizeof(double));
            pcdmio->regionHoleLists[i] = thisRegionHoleList;
            for (j = 0; j < (int) region->holes.size(); j++) {
                thisRegionHoleList[j * 2 + 0] = region->holes[i].first;
                thisRegionHoleList[j * 2 + 1] = region->holes[i].second;
            }
            
            i++;
#endif
            
            region->writePoly(0, "");
            
            END_TIMING(TIME_IO_WRITE);
            
            BEGIN_TIMING(TIME_FREE);    
            delete region;
            END_TIMING(TIME_FREE);
        }
        gidToRegion.clear();
        regionToProc.clear();
        
        barrier();
        
#ifndef PCDM_LIBRARY
        pcdmFinalize();
#endif

        statistics(false);
        
    }
    catch (std::exception& e) {             
        fprintf(stderr, "%s\n", e.what());
        fflush(stderr);                    
        exit(1);                         
    }                                       
    
#ifndef PCDM_LIBRARY
    return 0;   
#endif
}
