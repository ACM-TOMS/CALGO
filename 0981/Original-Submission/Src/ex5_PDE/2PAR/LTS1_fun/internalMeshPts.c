/** internalMeshPts

    Return a (col-wise) matrix XY of size (n^2, 2) where n = NXYval-2:
            XY(row,col) = XY[col*n^2 + row]
    The two columns of XY contain the abscissas and ordinates of the
    internal mesh points.
 */

#include <stdlib.h>
#include <stdio.h>

double *internalMeshPts(unsigned int NXYval, double XYmin, double XYmax)
{
    unsigned int n = NXYval-2;
    double h = 1.0/(NXYval-1); /* mesh step */

    double *XY = (double*)malloc(n*n*2*sizeof(double));
    if ( XY == NULL )
    {   fprintf(stderr, "\n***   ERROR IN internalMeshPts: DYNAMIC ALLOCATION OF XY IS FAILED. ***\n");
        exit(1);
    }

    unsigned int i,j;
    for (j=0; j<n; j++)     /* y(j) */
        for (i=0; i<n; i++) /* x(i) */
        {
            XY[i*n+j]     = XYmin+h*(j+1); /* XY(row,1) */
            XY[(i+n)*n+j] = XYmin+h*(i+1); /* XY(row,2) */
        }

    return XY;
}
