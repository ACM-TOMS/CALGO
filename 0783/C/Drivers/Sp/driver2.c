/*  Driver routine for Pcp2Nurb
 *  
 *	generates a twisted cube from a planar cut polyhedron (pcp)
 *      to illustrate surface of medium complexity
 *
 *      HINT: if confused read the ACM TOMS paper
 *
 *	     Pcp2Nurb: smooth free-form surfacing
 *	     with linearly-trimmed bicubic B-splines
 *				by Jorg Peters
 *
 *  input: file with (72) pcp subnets
 *  output: Open Inventor NurbsSurfaces
 *
 *  jorg@cs.purdue.edu Jan 1996
 *  IRIX 5.3
 *		COMPILATION and EXECUTION sequence
 *  cc -g -c nurb_iv2.c; cc -g -c pcp2nurb.c;
 *  cc -g -o nurb_iv2 nurb_iv2.o pcp2nurb.o -lm
 *  nurb_iv2 < in2 > out; ivview out;
 */

#define PTS	6	/* number of B-spline control points per row */
#define DEG	3	/* degree of B-spline patch */
#define PFL	4	/* number of points in profile (trimming) curve */
#define KTS	DEG+PTS+1	/* number of knots per row */
#define DIM	3	
#define SUB	9	/* size of subnet */

void  Pcp2Nurb( float Ci[][DIM], int valence[],
		float ctl_pt[][DIM], float knots[]);

main()
{
/*--------------------------VARS------------------------------------*/

    /* nine-point subnet */
    float	Pcp[SUB][DIM];  
    int		valence[2];

    /* profile (trimming) curve and control points of NURBS */
    float	profil[4][2] = {
		    {0,3},
		    {-3,0},
		    {0,-3},
		    {3,0} };
    float	ctl_pt[PTS*PTS][DIM], knots[KTS];

    int		i,j, no_subnets;
    void draw_patch( float profil[][2], float ctl_pt[][DIM], float knots[]);

/*------------------------------------------------------------------*/
    /* write the Inventor file */
    printf("#Inventor V2.0 ascii\n\n");
    printf("Separator {\n");

    scanf("%d", &no_subnets);
    for (i=0; i< no_subnets; i++) {
	scanf("%d %d",&valence[0], &valence[1]);
	for (j=0; j< SUB; j++) 
	    scanf("%f %f %f",&Pcp[j][0], &Pcp[j][1], &Pcp[j][2]);

	/* transform subset Pcp of planar cut polyhedron to a B-spline patch */
	Pcp2Nurb(Pcp,valence, ctl_pt, knots);

	draw_patch(profil,ctl_pt,knots);
    }
    printf("}\n");
/*------------------------------------------------------------------*/
}

void draw_patch( float	profil[][2], float ctl_pt[][DIM], float knots[])
{
    int	i;
    printf("Separator {\n");
	printf("Material { ambientColor %f %f %f }\n", 
	    0.8, 0.1, 0.0);
	printf("ProfileCoordinate2 {\n");
	printf("\tpoint  [\n");
	for (i=0; i< PFL-1; i++) 
	    printf("\t\t%2.1f %2.1f ,\n", profil[i][0], profil[i][1]);
	i = PFL-1;
	printf("\t\t%2.1f %2.1f\n", profil[i][0], profil[i][1]);
	printf("\t]\n");
	printf("}\n");

	printf("NurbsProfile {\n");
	printf("\tknotVector[\n\t 0,");
	for (i=0; i<= PFL; i++) 
	    printf(" %d,",i);
	i = PFL;
	printf(" %d]\n",i);
	printf("\tindex[\n\t");
	for (i=0; i< PFL; i++) 
	    printf(" %d,",i);
	i = 0;
	printf(" %d]\n",i);
	printf("}\n");

	printf("Coordinate3 {\n");
	printf("\tpoint  [\n");
	for (i=0; i< PTS*PTS-1; i++) 
	    printf("\t\t%f %f %f,\n", ctl_pt[i][0], ctl_pt[i][1], ctl_pt[i][2]);
	i = PTS*PTS-1;
	printf("\t\t%f %f %f\n", ctl_pt[i][0], ctl_pt[i][1], ctl_pt[i][2]);
	printf("\t]\n");
	printf("}\n");

	printf("NurbsSurface {\n");
	printf("\tnumUControlPoints\t %d\n", PTS);
	printf("\tnumVControlPoints\t %d\n", PTS);
	printf("\tuKnotVector [\n\t\t");
	for (i=0; i< KTS-1; i++) 
	    printf(" %2.1f,",knots[i]);
	i = KTS-1;
	printf(" %2.1f]\n",knots[i]);
	printf("\n\tvKnotVector [\n\t\t");
	for (i=0; i< KTS-1; i++) 
	    printf(" %2.1f,",knots[i]);
	i = KTS-1;
	printf(" %2.1f]\n",knots[i]);
	printf("}\n");
    printf("}\n");
}
