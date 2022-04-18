/*
	Hilbert.c
	
	Generate the points of the Hilbert curve by recursion.
	
	Greg Breinholt, ETHZ, 1997.
	
	email: 	breinholt@computer.org
	
	Compiles with standard ANSI C.
	
*/




#include <stdio.h>
#include <math.h>

void Hilbert(int x, int y, int lg, int i1, int i2);

int resolution = -1;




void main(void) {
	int width = -1;
	float p = 0.0;
	
	while (1) { 												/* Repeat forever */
		while (width < 2) {  										/* Check valid width */
			printf("Enter the Hilbert curve width.\n");
			scanf("%d", &width);

			if (width<0){
			  exit(0);
			}
			
			p = (log10(width)/ log10(2));								/* Check width is result of 2^m*/
			if (p != ((int)p)) {								
				printf("Curve width must be >= 2, and the result of 2^m (m = 1,2,3,4...).\n\n");
				width = -1;
			}
		}
		
		while (resolution < 0) {										/* Check valid resolution */
			printf("Enter the plotting resolution (1 = points only, 3 = aesthetic curve plotting).\n");
			scanf("%d", &resolution);
		}
		
		printf("\n");
		Hilbert(0,0,width,0,0);										/* Start recursion */
		printf("\n\n");
		
		width = -1; 
		resolution = -1;
	}		
}


void Hilbert(int x, int y, int lg, int i1, int i2){
			/* 	 x:	initial x co-ordinate, 
				 y:	initial y co-ordinate, 
				lg:	curve width (2^m), 
				i1:	starting point of the unit shape, 
				i2:	end point of the unit shape.  
			*/

	if (lg==1) {													/* Unit shape reached */
			printf("%d%c%d\n",x*resolution,',',y*resolution);							/* Output co-ordinates to console */
			return;												/* Exit recusion*/
	}
		
	lg >>= 1;														/* Divide by 2 */
	Hilbert(x + i1*lg, 	y + i1*lg, 	lg,	 i1,	1-i2);
	Hilbert(x + i2*lg, 	y + (1- i2)*lg,	lg,	 i1,	i2);
	Hilbert(x + (1-i1)*lg, 	y + (1-i1)*lg,	lg,	 i1,	i2);
	Hilbert(x + (1-i2)*lg, 	y + i2*lg, 	lg,	 1-i1,	i2);
}



/* Note
			Use: 	printf("%d%c%d\n",x * resolution,',',y * resolution);
					to output the points needed to plot a line curve.
		
			Use: 	LineTo(x * resolution, y * resolution); 
					to draw the curve.
					
			Use resolution 	= 1 to give just the points, 
							= 3 to space the points aesthetically,
*/
