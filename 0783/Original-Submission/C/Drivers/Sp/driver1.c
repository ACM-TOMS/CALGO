/*  Driver routine for Pcp2Nurb
 *	generates three Inventor NurbsSurface s
 *      that smooth the corner of a cube with vertices +/-2 
 *	using an adjustable blend ratio with default 0.5.
 *	A fourth biquadratic patch is C1 attached.
 *
 *      HINT: if confused read the ACM TOMS paper
 *
 *	     Pcp2Nurb: smooth free-form surfacing
 *	     with linearly-trimmed bicubic B-splines
 *				by Jorg Peters
 *
 *
 *  jorg@cs.purdue.edu Jan 1996
 *  IRIX 5.3
 *		COMPILATION and EXECUTION sequence
 *  cc -g -c nurb_iv1.c; cc -g -c pcp2nurb.c;
 *  cc -g -o nurb_iv1 nurb_iv1.o pcp2nurb.o -lm
 *  nurb_iv1 0.5 > out; ivview out;
 *  (try also other values in [0,1]: 0.35 0.65 0.0 1.0 --
 *   10 and -10 are not  legal values -- so try them, too)
 */

#define PTS	6	/* number of B-spline control points per row */
#define DEG	3	/* degree of B-spline patch */
#define PFL	4	/* number of points in profile (trimming) curve */
#define KTS	DEG+PTS+1	/* number of knots per row */
#define DIM	3	

const float	m_pi = 3.141592654;
void  Pcp2Nurb( float Ci[][DIM], int valence[],
		float ctl_pt[][DIM], float knots[]);

main( int argc, char *argv[])
{
/*--------------------------VARS------------------------------------*/
    /* profile (trimming) curve and control points of NURBS */
    float	profil[4][2] = {
		    {0,3},
		    {-3,0},
		    {0,-3},
		    {3,0} };
    float	ctl_pt[PTS*PTS][DIM], knots[KTS];

    /*  9 (control and vertex) points of a planar-cut polyhedron
     *          arranged in the order (cf. Figure 5 of the accompanying paper)
     *      8 7 6
     *      1 0 5    4,8 are vertex points
     *      2 3 4
     */
    float	Pcp[9][DIM] = {
		      {  1, 1, 2},	
		      { -1, 1, 2},
		      { -1, 2, 1},
		      {  1, 2, 1},
		      {  1.333333, 1.333333, 1.333333},
			  /* V_1 = avg of 1,1,2; 1,2,1; 2,1,1 */
		      {  2, 1, 1},
		      {  2,-1, 1},
		      {  1,-1, 2},
		      {  0, 0, 2}};
    int		valence[2] = {3,4};	/* V_1 has 3 neighbors */

    /* control points for an abutting biquadratic patch:
     * reuse Pcp[3], Pcp[0], Pcp[7], Pcp[1], Pcp[2], and
     * specify 4 additional vertices: */
    float       Cbi2[4][DIM] = {
			{ -1, -1, 2},	
			{ -3, -1, 2}, { -3, 1, 2}, { -3, 2, 1}};
    void draw_rot_patch( float axis[], float angle, float acol[]);
    void draw_patch( float profil[][2], float ctl_pt[][DIM], float knots[]);
    void draw_bi2( float ctl_pt[][DIM], float more_ctl_pt[][DIM]);
    void draw_cube();
    float	axis[3], acol[3], 
		blend_ratio = 0.5;
    double	atof();
    int		i;

/*------------------------------------------------------------------*/
    if (argc > 1)
	blend_ratio = (float) atof(argv[1]);		
    /* to minimize the amount of code describing the planar cut polyhedron
       I take advantage of symmetry:     */
    Pcp[0][0] = Pcp[0][1] = Pcp[3][0] = Pcp[3][2] = Pcp[5][1] = Pcp[5][2] = 
	2*(1-blend_ratio)+ 0*blend_ratio;
    /* adjust centroids in case the blend ratio has been changed */
    for (i=0; i< 3; i++) {
	Pcp[4][i] = (Pcp[0][i]+Pcp[3][i]+Pcp[5][i])/3;
	Pcp[8][i] = (Pcp[1][i]+Pcp[0][i]+Pcp[7][i]+Cbi2[0][i])/4;
    }


    /* transform subset Pcp of planar cut polyhedron to a B-spline patch */
    Pcp2Nurb(Pcp,valence, ctl_pt, knots);

    /* write the Inventor file */
    printf("#Inventor V2.0 ascii\n\n");
    printf("Separator {\n");
	draw_cube();
	draw_patch(profil,ctl_pt,knots);
	for (i=0; i<3; i++)
	    axis[i] = 1.0;
	acol[0]= 0.1;
	acol[1]= 0.5;
	acol[2]= 0.7;
	draw_rot_patch(axis,2*m_pi/3, acol);
	acol[0]= 0.7;
	acol[1]= 0.5;
	acol[2]= 0.1;
	draw_rot_patch(axis,4*m_pi/3, acol);
	draw_bi2(Pcp, Cbi2);
    printf("}\n");
/*------------------------------------------------------------------*/
}

void draw_rot_patch( float axis[], float angle, float acol[])
{
    printf("Separator {\n");
	printf("Rotation { rotation %f %f %f %f }\n", 
	    axis[0], axis[1], axis[2], angle);
	printf("Material { ambientColor %f %f %f }\n", 
	    acol[0], acol[1], acol[2]);
	printf("USE PATCH\n"); 
    printf("}\n");
}

void draw_patch( float	profil[][2], float ctl_pt[][DIM], float knots[])
{
    int	i;
    printf("DEF PATCH \n");
    printf("Separator {\n");
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

void draw_bi2( float ctl_pt[][DIM], float more_ctl_pt[][DIM])
{
    int	i,j,kts,dg;

    dg = 2;  /* surface degree */
    printf("Separator {\n");

	printf("Material { ambientColor %f %f %f }\n", 
	    1.0, 0.0, 0.0);
	printf("Coordinate3 {\n");
	printf("\tpoint  [\n");
	i = 7;
	printf("\t\t%f %f %f,\n", ctl_pt[i][0], ctl_pt[i][1], ctl_pt[i][2]);
	i = 0;
	printf("\t\t%f %f %f,\n", ctl_pt[i][0], ctl_pt[i][1], ctl_pt[i][2]);
	i = 3;
	printf("\t\t%f %f %f,\n", ctl_pt[i][0], ctl_pt[i][1], ctl_pt[i][2]);
	j = 0;
	printf("\t\t%f %f %f,\n", 
	    more_ctl_pt[j][0], more_ctl_pt[j][1], more_ctl_pt[j][2]);
	i = 1;
	printf("\t\t%f %f %f,\n", ctl_pt[i][0], ctl_pt[i][1], ctl_pt[i][2]);
	i = 2;
	printf("\t\t%f %f %f,\n", ctl_pt[i][0], ctl_pt[i][1], ctl_pt[i][2]);
	j = 1;
	printf("\t\t%f %f %f,\n", 
	    more_ctl_pt[j][0], more_ctl_pt[j][1], more_ctl_pt[j][2]);
	j = 2;
	printf("\t\t%f %f %f,\n", 
	    more_ctl_pt[j][0], more_ctl_pt[j][1], more_ctl_pt[j][2]);
	j = 3;
	printf("\t\t%f %f %f\n", 
	    more_ctl_pt[j][0], more_ctl_pt[j][1], more_ctl_pt[j][2]);
	printf("\t]\n");
	printf("}\n");

	kts = dg+3+1; /* number of knots in row = dg + coeffs + 1 */
	printf("NurbsSurface {\n");
	printf("\tnumUControlPoints\t %d\n", dg+1);
	printf("\tnumVControlPoints\t %d\n", dg+1);
	printf("\tuKnotVector [\n\t\t");
	for (i=0; i< kts-1; i++) 
	    printf(" %2.1f,",1.0*i);
	i = kts-1;
	printf(" %2.1f]\n",1.0*i);
	printf("\n\tvKnotVector [\n\t\t");
	for (i=0; i< kts-1; i++) 
	    printf(" %2.1f,",1.0*i);
	i = kts-1;
	printf(" %2.1f]\n",1.0*i);
	printf("}\n");
    printf("}\n");
}

/* write Inventor node for cube */
void draw_cube()
{
    printf("Separator {\n");
	printf("DrawStyle {\n");
	    printf("\tstyle LINES \n");
	printf("}\n");
	printf("Transform {\n");
	    printf("\tscaleFactor %f %f %f\n",2.0,2.0,2.0);
	printf("}\n");
	printf("Cube {\n");
	printf("\t width 2 \n");
	printf("\t height 2 \n");
	printf("\t depth 2 \n");
	printf("}\n");
    printf("}\n");
}
