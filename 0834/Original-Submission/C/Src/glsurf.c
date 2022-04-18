
/*
 *  glsurf.c:  Surface plotting program using OpenGL.
 *
 *  R. Renka
 *  renka@cs.unt.edu
 *  11/29/03
 *
 *  This program reads in a data set defining a surface
 *  (and, optionally, another data set defining a space
 *  curve assumed to lie in the surface), and displays
 *  the surface (and curve) from an interactively chosen
 *  perspective with various user-specified options.
 *
 *  Compile and link command for Linux with Mesa:
 *
 *     gcc -O3 -o glsurf glsurf.c -L/usr/X11R6/lib -lm -lX11 
 *         -lXmu -lXi -lXext -lglut -lMesaGL -lMesaGLU 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>

struct rect {    /* Rectangle in color selection window */
   int x0, y0;   /* Window coordinates of lower left corner */
   int w, h;     /* Width and height in window coordinates */
};

/*
 * Static array dimensions:
 */
#define IW1 32     /* Width of image1 (must be a power of 2) */
#define IW2 64     /* Width of image2 (must be a power of 2) */
#define IW3 32     /* Width of image3 (must be a power of 2) */
#define NC 36      /* Number of contour intervals and
                        associated color table entries */
#define LTITLE 80  /* Length of title string */
#define NH 180     /* Number of hue values in hstable */
#define NS 140     /* Number of saturation values in hstable */
#define NV 130     /* Number of intensity values for color
                        selection */

/*
 *  Global variables:
 *
 *  aspect = Flag with value True iff the aspect ratio of
 *           the viewing volume (intersected with the 
 *           projection plane) is to be preserved.
 *  bbox = Flag with value True iff the bounding box is
 *         to be displayed.
 *  bottom_color = Material color of back faces (bottom of
 *                 surface).
 *  colors = Array dimensioned nc by 3 containing color
 *           values (red, green, and blue components)
 *           stored in one-to-one correspondence with
 *           contours.
 *  contours = Array of length nc containing an increasing
 *             sequence of function values defining contour
 *             levels.  Refer to cFillTriangle.
 *  cplot = Flag with value true iff a (flat) contour plot
 *          is to be created.
 *  cpoints = Space curve vertex array.
 *  cs_color = Currently selected color (RGB) in color
 *             selection window.
 *  curve_color = Color for drawing space curve.
 *  cutaway = Fraction of the surface to be cut away by
 *            clippling against the near plane.
 *  decal = Flag with value True if the image1 (color-filled
 *          contours) colors are to replace the surface
 *          color (both top_color and bottom_color);  value
 *          False if the image1 colors are to be modulated
 *          by top_color and bottom_color.
 *  dmin,dmax = Initial values of zmin and zmax (before any
 *              scaling):  used for labels on bounding box.
 *  font_number = Font number for call to displayString.
 *  have_vnormals = Flag with value True iff vertex normals
 *                  vnormals are included in the input data
 *                  set.
 *  hide_curve = Flag with value True iff depth testing is
 *               to be enabled while drawing a space curve.
 *  hstable = Array of RGB color components associated with
 *            (hue,saturation) (column,row) pairs, V = 1.
 *  image1 = Array dimensioned IW1 by 3 containing a texture
 *           image used to produce a color-filled contour
 *           plot of z values (pasted onto the surface).
 *  image2 = Array dimensioned IW2 by 3 containing a texture
 *           image used to produce contour lines (represent-
 *           ing z values) on the surface.
 *  image3 = Array dimensioned IW3 by 3 containing a texture
 *           image used to produce contour lines (represent-
 *           ing mean curvature) on the surface.
 *  indc = Array of length nv containing color/contour
 *         indices (0 to nc-1) in one-to-one correspondence
 *         with vertices.
 *  indices = Array of vertex indices defining triangles.
 *  lbdown = Flag with value True iff the left mouse button
 *           is down.
 *  light_angle = Rotation angle (in degrees) associated
 *                with light_axis defining the location of
 *                a point light source; altered by 'x' and
 *                'y' keys.
 *  light_axis = Axis of rotation defining the location of
 *               a point light source; altered by 'x' and
 *               'y' keys.
 *  lighting = Flag with value true iff an illumination
 *             model is to be applied.
 *  line_width = Line width for drawing the space curve if
 *               np > 0.
 *  nc = Number of contour intervals.
 *  np = Number of points (vertices) defining a space curve
 *       or 0 if no space curve file is specified.
 *  nt = Number of triangles.
 *  nv = Number of vertices in the surface.
 *  orient_angle = Rotation angle (in degrees) associated
 *                 with orient_axis defining the viewer
 *                 orientation; altered by x, y, and z keys.
 *  orient_axis = Axis of rotation defining the viewer
 *                orientation (pitch, yaw, and roll);
 *                altered by x, y, and z keys.
 *  perspective = Flag with value False for orthographic
 *                projection, or True for perspective
 *                projection.
 *  polyfill = Flag with value False for drawing a wireframe
 *             mesh, or True for drawing filled triangles.
 *  r = Radius of a bounding sphere centered at (xc,yc,zc).
 *  rhstable = Rectangle containing hstable.
 *  rvtable = Rectangle containing a table of intensities.
 *  rcolor = Rectangle containing the current color.
 *  rh,rs,rv,rr,rg,rb = Rectangles containing text strings
 *                      with HSV and RGB values.
 *  raccept = Rectangle representing 'Accept' button.
 *  rcancel = Rectangle representing 'Cancel' button.
 *  shininess = Shininess factor (exponent in specular
 *              lighting term).
 *  smooth = Flag with value False for flat shading, or
 *           True for smooth shading.
 *  sp_color = Specular color.
 *  tex_flag = Flag containing a texture name, or 0 if
 *             no texture is to be applied.
 *  title = Optional character string written on the plot.
 *  title_pos = Raster position of (lower left corner of)
 *              title.
 *  top_color = Material color of front faces.
 *  trareas = Array of triangle areas.
 *  trnormals = Array of triangle normals used to compute
 *              vnormals and/or for display with flat
 *              shading.
 *  vareas = Array of vertex areas (one-third of the sum of
 *           areas of the triangles containing a vertex).
 *  vcrvs = Array of vertex mean curvatures normalized to
 *          [0,1] for use as 1-D texture coordinates
 *  vertices = Vertex array.
 *  view_angle = Rotation angle (in degrees) associated
 *               with view_axis defining the viewer
 *               location; altered by arrow keys.
 *  view_axis = Axis of rotation defining the viewer
 *              location; altered by arrow keys.
 *  vnormals = Array of vertex normals.
 *  vsf = Scale factor for r defining the size of the
 *        (square) viewing volume (intersected with the
 *        projection plane).
 *  win_height,win_width = height and width of window 1:
 *                         stored by reshape and used by
 *                         functions mouse and specialKey.
 *  win3_height,win3_width = height and width of window 3:
 *                           stored by cckReshape and used
 *                           by function cckDisplay.
 *  xs,ys = Starting coordinates (relative to the upper left
 *          corner of the window) of a mouse move (triggered
 *          by a left button press) used to rotate about the
 *          x or y axis.
 *  xc,yc,zc = Coordinates of the center c of the bounding
 *             box associated with the vertices.
 *  xmin,xmax,ymin,ymax,zmin,zmax = Coordinates of the
 *                                  bounding box.
 *  zcon = Height of the (flat) contour plot above the x-y
 *         plane (below the plane if negative).
 */
static GLboolean aspect = GL_TRUE;
static GLboolean bbox = GL_TRUE;
static GLfloat bottom_color[4] = {0.1, 0.5, 0.8, 1.0};
static GLdouble contours[NC];
static GLboolean cplot = GL_FALSE;
static GLdouble *cpoints;
static GLfloat cs_color[3] = {0.78, 0.57, 0.11};
static GLfloat curve_color[3] = {1.0, 0.0, 0.0};  /* red */
static GLdouble cutaway = 0.0;
static GLboolean decal = GL_TRUE;
static GLdouble dmin, dmax;
static unsigned int font_number = 3;
static GLboolean have_vnormals = GL_FALSE;
static GLboolean hide_curve = GL_FALSE;
static GLubyte hstable[NS][NH][3];
static GLuint *indc, *indices;
static GLboolean lbdown = GL_FALSE;
static GLdouble light_angle = 0.0;
static GLdouble light_axis[3] = {1.0, 0.0, 0.0};
static GLboolean lighting = GL_TRUE;
static GLfloat line_width = 2.0;
static int nc = NC;
static int np = 0;
static int nt, nv;
static GLdouble orient_angle = 0.0;
static GLdouble orient_axis[3] = {1.0, 0.0, 0.0};
static GLboolean perspective = GL_FALSE;
static GLboolean polyfill = GL_TRUE;
static GLdouble r;
struct rect rhstable = {25, 235, NH, NS};
struct rect rvtable = {257, 240, 18, NV};
struct rect rcolor = {25, 105, 60, 60};
struct rect rh = {128, 173, 49, 25};
struct rect rr = {223, 173, 49, 25};
struct rect rs = {128, 121, 49, 25};
struct rect rg = {223, 121, 49, 25};
struct rect rv = {128, 68, 49, 25};
struct rect rb = {223, 68, 49, 25};
struct rect raccept = {50, 25, 90, 25};
struct rect rcancel = {160, 25, 90, 25};
static GLfloat shininess = 32.0;
static GLboolean smooth = GL_TRUE;
static GLfloat sp_color[4] = {1.0, 1.0, 1.0, 1.0};  /* white */
static int tex_flag = 0;
static char title[LTITLE] = "";
static GLfloat title_pos[3] = {0., 0., 0.};
static GLfloat top_color[4] = {0.78, 0.57, 0.11, 1.0};
static GLdouble *trareas, *trnormals, *vareas, *vcrvs,
                *vertices, *vnormals;
static GLdouble view_angle = -90.0;
static GLdouble view_axis[3] = {1.0, 0.0, 0.0};
static GLdouble vsf = 1.0;
static GLfloat win_height, win_width;
static GLfloat win3_height, win3_width;
static int xs, ys;
static GLdouble xc, yc, zc;
static GLdouble xmin, xmax, ymin, ymax, zmin, zmax;
static GLdouble zcon;

/*
 *  Color table and texture images:
 */
static GLfloat colors[NC][3] = {
          {0.0, 0.0, .75}, {.75, 0.0, .75}, {0.0, .75, 0.0},
          {.75, .75, 0.0}, {0.0, .75, .75}, {.75, 0.0, 0.0},
          {0.0, 0.0, 0.8}, {0.8, 0.0, 0.8}, {0.0, 0.8, 0.0},
          {0.8, 0.8, .75}, {0.0, 0.8, 0.8}, {0.8, 0.0, 0.0},
          {0.0, 0.0, .85}, {.85, 0.0, .85}, {0.0, .85, 0.0},
          {.85, .85, 0.0}, {0.0, .85, .85}, {.85, 0.0, 0.0},
          {0.0, 0.0, 0.9}, {0.9, 0.0, 0.9}, {0.0, 0.9, 0.0},
          {0.9, 0.9, 0.0}, {0.0, 0.9, 0.9}, {0.9, 0.0, 0.0},
          {0.0, 0.0, .95}, {.95, 0.0, .95}, {0.0, .95, 0.0},
          {.95, .95, 0.0}, {0.0, .95, .95}, {.95, 0.0, 0.0},
          {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {0.0, 1.0, 0.0},
          {1.0, 1.0, 0.0}, {0.0, 1.0, 1.0}, {1.0, 0.0, 0.0}};
/*
 *  image1 is applied in decal mode (replaces the surface
 *  color) if decal = True, or modulate mode (each component
 *  of the texture color is scaled by the corresponding
 *  component of the fragment color) otherwise.
 */
static GLubyte image1[IW1][3] = {{100, 0, 0},  {0, 0, 100},
               {100, 54, 0},  {100, 0, 100}, {0, 100, 0},
               {54, 0, 100},  {100, 100, 0}, {0, 100, 100},
               {137, 0, 0},   {0, 0, 137},   {137, 73, 0},
               {137, 0, 137}, {0, 137, 0},   {73, 0, 137},
               {137, 137, 0}, {0, 137, 137}, {187, 0, 0},
               {0, 0, 187},   {187, 100, 0}, {187, 0, 187},
               {0, 187, 0},   {100, 0, 187}, {187, 187, 0},
               {0, 187, 187}, {255, 0, 0},   {0, 0, 255},
               {255, 137, 0}, {255, 0, 255}, {0, 255, 0},
               {137, 0, 255}, {255, 255, 0}, {0, 255, 255}};
/*
 *  image2 is applied in blend mode:  the second column (of
 *  color triples) corresponds to a contour line (with color
 *  specified in a call to glTexEnvfv); the remaining
 *  columns are initialized to zeros leaving the fragment
 *  color (surface color) unaltered.
 */
static GLubyte image2[IW2/8][8][3] = {
                              {{0, 0, 0}, {255, 255, 255}},
                              {{0, 0, 0}, {255, 255, 255}},
                              {{0, 0, 0}, {255, 255, 255}},
                              {{0, 0, 0}, {255, 255, 255}},
                              {{0, 0, 0}, {255, 255, 255}},
                              {{0, 0, 0}, {255, 255, 255}},
                              {{0, 0, 0}, {255, 255, 255}},
                              {{0, 0, 0}, {255, 255, 255}}};
/*
 *  image3 (mean curvature) is applied in decal mode
 *  (replaces the surface color).
 */
static GLubyte image3[IW3][3] = {{0, 0, 255},  {15, 0, 255},
               {31, 0, 255},  {47, 0, 255}, {63, 0, 255},
               {79, 0, 255},  {95, 0, 255}, {111, 0, 255},
               {127, 0, 255},  {143, 0, 255}, {159, 0, 255},
               {175, 0, 255},  {191, 0, 255}, {207, 0, 255},
               {223, 0, 255},  {239, 0, 255}, {255, 0, 255},
               {255, 0, 239},  {255, 0, 223}, {255, 0, 207},
               {255, 0, 191},  {255, 0, 175}, {255, 0, 159},
               {255, 0, 143},  {255, 0, 127}, {255, 0, 111},
               {255, 0, 95},  {255, 0, 79}, {255, 0, 63},
               {255, 0, 47},  {255, 0, 31}, {255, 0, 15}};

/*
 *  Function prototypes
 */
static void benchmark(void);
static void cckDisplay(void);
static void cckKey(unsigned char key, int x, int y);
static void cckReshape(int w, int h);
static void cFillTriangle(int i1, int i2, int i3,
                          GLdouble zval);
static void computeNormals(void);
static void csDisplay(void);
static void csInit(void);
static void csKey(unsigned char key, int x, int y);
static void csMotion(int x, int y);
static void csMouse(int button, int state, int x, int y);
static void csReshape(int w, int h);
static void display(void);
static int displayString(GLfloat pos[3], char *string,
                         unsigned int fontno);
static void drawRect(GLuint x0, GLuint y0, GLuint w, GLuint h,
                     GLfloat brcolor[3], GLfloat tlcolor[3]);
static float fmax(float x, float y, float z);
static int getLin(char *s, int len);
static void getRasterPos(int button, int state, int x, int y);
static void hsv2rgb(float h, float s, float v, float *r,
                    float *g, float*b);
static void init1(void);
static void inputData(int argc, char *argv[]);
static void key(unsigned char key, int x, int y);
static void makeMenu(void);
static void menu(int item);
static void motion(int x, int y);
static void mouse(int button, int state, int x, int y);
static void normalizeVC(void);
static void reshape(int w, int h);
static void rgb2hsv(float r, float g, float b, float *h,
                    float *s, float *v);
static void scaleData(GLdouble sf);
static int screenDump(int reversebw);
static void specialKey(int key, int x, int y);
static void storeC(void);
static void updateR (unsigned char c, GLdouble d,
                     GLdouble *a, GLdouble *u);


/*
 *----------------------------------------------------------
 *
 *  benchmark:  Compute and print the numbers of frames and
 *              triangles displayed per second (for the
 *              current perspective) averaged over 5 seconds.
 *
 *  glsurf function called:  display
 *  OpenGL functions called:  glutGet
 *
 *----------------------------------------------------------
 */
static void benchmark(void)
{
   int start_time, end_time;
   int ndraws;
   double seconds, fps, tps;

   printf("Benchmarking...\n");

   ndraws = 0;
   start_time = glutGet(GLUT_ELAPSED_TIME);
   do {
      display();
      ndraws++;
      end_time = glutGet(GLUT_ELAPSED_TIME);
   } while (end_time - start_time < 5000);   /* 5 seconds */

   seconds = (double) (end_time - start_time)/1000.0;
   tps = nt*ndraws/seconds;
   fps = (double) ndraws/seconds;
   printf("Result:  triangles/sec: %g  frames/sec: %g\n",
          tps, fps);
   return;
}


/*
 *----------------------------------------------------------
 *
 *  cckDisplay:  Contour color key window display callback.
 *               Partitions the window into a vertically
 *               stacked set of color-filled rectangles with
 *               colors from image1.
 *
 *  glsurf function called:  DisplayString
 *  OpenGL functions called:  glColor3fv, glColor3ubv,
 *                            glFlush, glRectf
 *
 *----------------------------------------------------------
 */
static void cckDisplay(void)
{
   GLfloat black[3] = {0.0, 0.0, 0.0};   /* Label colors */
   GLfloat white[3] = {1.0, 1.0, 1.0};
   GLfloat dy;                /* Rectangle height */
   unsigned int fontno = 5;   /* Font for labels */
   int k;                     /* Rectangle/color index */
   GLfloat pos[3];            /* Label position */
   char text[11];             /* Label text */
   GLfloat x1,y1;    /* Coordinates of lower left corner */
   GLfloat x2,y2;    /* Coordinates of upper right corner */

   x1 = 0.0;
   x2 = win3_width;
   y1 = 0.0;
   dy = win3_height/(GLfloat) IW1;
   for (k = 0; k < IW1; k++) {
      y2 = y1 + dy;
      glColor3ubv(&image1[k][0]);
      glRectf(x1, y1, x2, y2);
      y1 = y2;
   }
/*
 *  Label lowest, middle, and highest rectangles.
 */
   glColor3fv(white);
   pos[0] = 3.0;
   pos[1] = 3.0;
   pos[2] = 0.0;
   sprintf(text, "%-10.4g", zmin);
   displayString(pos, text, fontno);

   pos[0] = 3.0;
   pos[1] = win3_height/2.0-dy+3.0;
   pos[2] = 0.0;
   sprintf(text, "%-10.4g", (zmin+zmax)/2.0);
   displayString(pos, text, fontno);

   glColor3fv(black);
   pos[1] = win3_height-dy+3.0;
   sprintf(text, "%-10.4g", zmax);
   displayString(pos, text, fontno);

   glFlush();
   return;
}


/*
 *----------------------------------------------------------
 *
 *  cckKey:  Keyboard callback associated with the contour
 *           color key window:
 *
 *           Escape ==> Cancel selection window.
 *
 *  The mouse coordinates (x,y) are not used.
 *
 *  glsurf function called:  None
 *  OpenGL functions called:  glutHideWindow, glutSetWindow
 *
 *----------------------------------------------------------
 */
static void cckKey(unsigned char key, int x, int y)
{
   if (key == 27) {                 /* Escape key pressed */
      glutHideWindow();
      glutSetWindow(1);
   }
   return;
}


/*
 *----------------------------------------------------------
 *
 *  cckReshape:  Set the viewport, projection matrix, and
 *               modelview matrix for the contour color key
 *               window (window 3):  two-dimensional render-
 *               ing.
 *
 *  This function is called before the first call to
 *  cckDisplay.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  glLoadIdentity, glMatrixMode,
 *                            glScalef, glTranslatef,
 *                            gluOrtho2D, glViewport
 *
 *----------------------------------------------------------
 */
static void cckReshape(int w, int h)
{
   win3_width = (GLfloat) w;
   win3_height = (GLfloat) h;
/*
 *  Set the viewport to the entire window (the default).
 */
   glViewport(0, 0, (GLint) w, (GLint) h);
/*
 *  Set up an orthogonal projection for data space equal to
 *  viewport (identity).
 */
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluOrtho2D(0.0, (GLfloat) w, 0.0, (GLfloat) h);
/*
 *  Set the modelview matrix to the identity.
 */
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();

   return;
}



/*
 *----------------------------------------------------------
 *
 *  cFillTriangle:
 *
 *  Given the vertices of a triangle along with vertex data
 *  values and color indices, and a sequence of contour
 *  levels with associated color values, this function
 *  creates a color-filled contour plot of the linear
 *  interpolant of the data values on the triangle;  i.e.,
 *  the triangle is partitioned into a set of polygons, each
 *  associated with a constant color value (range of data
 *  values), and the polygons are rendered with color-fill.
 *
 *  Note that this function calls only the OpenGL functions
 *  glBegin, glEnd, glColor3fv, and glVertex3fv.  Thus, an
 *  appropriate drawing context must be established before
 *  the first call to this function:  create a window and
 *  set up the viewing and projection matrices.  Also, 
 *  lighting should not be enabled when this function is
 *  called.
 *
 *  Input parameters:
 *
 *      i1,i2,i3 = Indices (into vertices and indc) of a
 *                 counterclockwise-ordered sequence of
 *                 vertices defining a triangle.  It is
 *                 assumed that the indices are in the
 *                 range 0 to nv-1 and distinct.
 *
 *      zval = Height of the contour plot relative to the
 *             x-y plane.  The third components of the 
 *             vertices are replaced by zval in the calls to
 *             glVertex3fv.  This allows the contour plot to
 *             be drawn above or below a surface plot.
 *
 *  Global data:
 *
 *      nv = Number of vertices.
 *
 *      vertices = Array dimensioned nv by 3 containing
 *                 vertex coordinates (first two columns)
 *                 and data values (third column).
 *
 *      nc = Number of contour levels (intervals).
 *
 *      contours = Array of length nc containing an increas-
 *                 ing sequence of function values defining
 *                 the contour levels:  contours[ic] is the
 *                 maximum function value associated with
 *                 contour level ic for ic in the range 0
 *                 to nc-1.  Thus, to partition [zmin,zmax]
 *                 uniformly into nc subintervals, set
 *                 contours[ic] to zmin + (ic+1)*dz for
 *                 dz = (zmax-zmin)/nc.
 *
 *      colors = Array dimensioned nc by 3 containing color
 *               values (red, green, and blue components)
 *               stored in one-to-one correspondence with
 *               contours.
 *
 *      indc = Array of length nv containing color/contour
 *             indices (0 to nc-1) in one-to-one correspon-
 *             dence with vertices.  If contour values are
 *             uniformly distributed in [zmin,zmax], then
 *             indc[i] = max{nc-1, (vertices[i][2]-zmin)/dz},
 *             for dz = (zmax-zmin)/nc, i = 0 to nv-1.
 *
 *  If two vertices of the triangle have identical data
 *  values, the triangle is not filled.
 *
 *  No parameters are altered by this function.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  glBegin, glEnd, glColor3fv,
 *                            glVertex3dv
 *
 *----------------------------------------------------------
 */
static void cFillTriangle(int i1, int i2, int i3,
                          GLdouble zval)
{
   int i, ic, ic1, ic2, ic3, j1, j2, j3, n;
   GLdouble s, x1, x2, x3, xp, xq, y1, y2, y3, yp, yq;
   GLdouble z1, z2, z3, zp;
   GLdouble vert[5][3];
/* 
 *  Store the vertex coordinates, data values, and color
 *  indices in local variables.  The vertices are cyclically
 *  permuted so that ic1 (the color index of the first
 *  vertex) is the smallest.
 */
   j1 = i1;
   j2 = i2;
   j3 = i3;
   ic1 = indc[j1];
   ic2 = indc[j2];
   ic3 = indc[j3];
   if (ic1 > ic2  ||  ic1 > ic3) {
      if (ic2 <= ic3) {
         j1 = i2;
         j2 = i3;
         j3 = i1;
      } else {
         j1 = i3;
         j2 = i1;
         j3 = i2;
      }
      ic1 = indc[j1];
      ic2 = indc[j2];
      ic3 = indc[j3];
   }

   x1 = vertices[3*j1];
   x2 = vertices[3*j2];
   x3 = vertices[3*j3];
   y1 = vertices[3*j1+1];
   y2 = vertices[3*j2+1];
   y3 = vertices[3*j3+1];
   z1 = vertices[3*j1+2];
   z2 = vertices[3*j2+2];
   z3 = vertices[3*j3+2];
/*
 *  The triangle is partitioned into a set of polygons
 *  defined by intersecting the triangle with half-spaces
 *  associated with contour/color indices (ranges of z-
 *  values).  The contour lines (and hence the polygons)
 *  are ordered by increasing indices.
 *
 *  n is the number of vertices in the current polygon.
 *  vert contains the ccw-ordered sequence of vertices.
 *
 *  Initialize the first polygon with p1 = (x1,y1).
 */
   n = 1;
   vert[0][0] = x1;
   vert[0][1] = y1;
   vert[0][2] = zval;
/*
 *  Loop on contour lines p-q such that p=(xp,yp) lies on
 *  edge p1-p2.
 */
   for (ic = ic1; ic < ic2; ic++) {
      zp = contours[ic];
      if (z1 == z2) return;
      s = (zp-z1)/(z2-z1);
      if (s < 0.0) s = 0.0;
      if (s > 1.0) s = 1.0;
      xp = x1 + s*(x2-x1);
      yp = y1 + s*(y2-y1);
      if (ic < ic3) {
/*
 *  q lies on edge p1-p3.
 */
         if (z1 == z3) return;
         s = (zp-z1)/(z3-z1);
         if (s < 0.0) s = 0.0;
         if (s > 1.0) s = 1.0;
         xq = x1 + s*(x3-x1);
         yq = y1 + s*(y3-y1);
      } else {
/*
 *  q lies on edge p3-p2.
 */
         if (z2 == z3) return;
         s = (zp-z3)/(z2-z3);
         if (s < 0.0) s = 0.0;
         if (s > 1.0) s = 1.0;
         xq = x3 + s*(x2-x3);
         yq = y3 + s*(y2-y3);
      }
/*
 *  Append p and q to the polygon.
 */
      vert[n][0] = xp;
      vert[n][1] = yp;
      vert[n][2] = zval;
      vert[n+1][0] = xq;
      vert[n+1][1] = yq;
      vert[n+1][2] = zval;
      n = n + 2;
      if (ic == ic3) {
/*
 *  Append p3 to the polygon.
 */
         vert[n][0] = x3;
         vert[n][1] = y3;
         vert[n][2] = zval;
         n = n + 1;
      }
/*
 *  Set the current color to colors[ic], and render
 *  the polygon.
 */
      glBegin(GL_POLYGON);
         glColor3fv((GLfloat *) &colors[ic]);
         for (i = 0; i < n; i++) {
            glVertex3dv(&vert[i][0]);
         }
      glEnd();
/*
 *  Initialize the next polygon with q and p.
 */
      n = 2;
      vert[0][0] = xq;
      vert[0][1] = yq;
      vert[0][2] = zval;
      vert[1][0] = xp;
      vert[1][1] = yp;
      vert[1][2] = zval;
   }
/*
 *  Append p2 to the polygon.  This completes the last
 *  polygon if p3 has already been included (ic3 < ic2).
 */
   vert[n][0] = x2;
   vert[n][1] = y2;
   vert[n][2] = zval;
   n = n + 1;
   ic = ic2;
   if (ic3 >= ic2) {
/*
 *  Loop on contour lines p-q such that p lies on edge p2-p3
 *  (and q lies on edge p3-p1).
 */
      for (ic = ic2; ic < ic3; ic++) {
         zp = contours[ic];
         if (z2 == z3) return;
         s = (zp-z2)/(z3-z2);
         if (s < 0.0) s = 0.0;
         if (s > 1.0) s = 1.0;
         xp = x2 + s*(x3-x2);
         yp = y2 + s*(y3-y2);
         if (z1 == z3) return;
         s = (zp-z1)/(z3-z1);
         if (s < 0.0) s = 0.0;
         if (s > 1.0) s = 1.0;
         xq = x1 + s*(x3-x1);
         yq = y1 + s*(y3-y1);
/*
 *  Append p and q.
 */
         vert[n][0] = xp;
         vert[n][1] = yp;
         vert[n][2] = zval;
         vert[n+1][0] = xq;
         vert[n+1][1] = yq;
         vert[n+1][2] = zval;
         n = n + 2;
/*
 *  Set the current color to colors[ic], and render
 *  the polygon.
 */
         glBegin(GL_POLYGON);
            glColor3fv((GLfloat *) &colors[ic]);
            for (i = 0; i < n; i++) {
               glVertex3dv(&vert[i][0]);
            }
         glEnd();
/*
 *  Initialize a new polygon with q and p.
 */
         n = 2;
         vert[0][0] = xq;
         vert[0][1] = yq;
         vert[0][2] = zval;
         vert[1][0] = xp;
         vert[1][1] = yp;
         vert[1][2] = zval;
      }
/*
 *  Append p3 to the polygon and set ic = ic3.
 */
      vert[n][0] = x3;
      vert[n][1] = y3;
      vert[n][2] = zval;
      n = n + 1;
      ic = ic3;
   }
/*
 *  Set the current color to colors[ic], and render
 *  the last polygon.
 */
   glBegin(GL_POLYGON);
      glColor3fv((GLfloat *) &colors[ic]);
      for (i = 0; i < n; i++) {
         glVertex3dv(&vert[i][0]);
      }
   glEnd();

   return;
}


/*
 *----------------------------------------------------------
 *
 *  computeNormals:  Compute triangle normals trnormals as
 *                   normalized cross products of directed
 *                   edges and, unless have_vnormals = True,
 *                   compute vertex normals vnormals as
 *                   normalized angle-weighted averages of
 *                   the triangle normals shared by each
 *                   vertex.  Triangle areas (trareas),
 *                   vertex areas (vareas), and approxima-
 *                   tions to mean curvatures at the
 *                   vertices (vcrvs) are also computed.
 *
 *  A vertex that is not contained in any triangle has its
 *  normal vector, area, and curvature set to zeros.
 *
 *  The program is terminated with exit code 1 if duplicate
 *  vertices are encountered.  A warning message is printed
 *  (but the program is not necessarily terminated) if a
 *  triangle has area zero.
 *
 *  glsurf functions called:  NormalizeVC
 *  OpenGL functions called:  None
 *
 *----------------------------------------------------------
 */
static void computeNormals(void) 
{
/*
 *  Local variables:
 *
 *  i,i1,i2,i3 =  Vertex indices
 *  k =           Triangle index
 *  a1,a2,a3 =    Interior angles of a triangle or their
 *                  cotangents
 *  dx1,dy1,dz1 = Components of edge from second to third
 *                  vertex
 *  dx2,dy2,dz2 = Components of edge from third to first
 *                  vertex
 *  dx3,dy3,dz3 = Components of edge from first to second
 *                  vertex
 *  s12,s13,s23,s21,s31,s32 = Terms used to compute vcrvs
 *  vn1,vn2,vn3 = Euclidean vector norms of the sides of a
 *                  triangle
 *  vnorm =       Euclidean vector norm or scale factor used
 *                  for normalization
 *  x,y,z =       Components of a vertex or triangle normal
 *
 */
   int i, k;
   GLuint i1, i2, i3;
   GLdouble a1, a2, a3, dx1, dx2, dx3, dy1, dy2, dy3, dz1,
            dz2, dz3, s12, s13, s23, s21, s31, s32,
            vn1, vn2, vn3, vnorm, x, y, z;
/*
 *  Compute triangle normals and areas.
 */
   for (k = 0; k < nt; k++) {
      i1 = indices[3*k];
      i2 = indices[3*k+1];
      i3 = indices[3*k+2];
      dx1 = vertices[3*i3] - vertices[3*i2];
      dy1 = vertices[3*i3+1] - vertices[3*i2+1];
      dz1 = vertices[3*i3+2] - vertices[3*i2+2];
      dx2 = vertices[3*i1] - vertices[3*i3];
      dy2 = vertices[3*i1+1] - vertices[3*i3+1];
      dz2 = vertices[3*i1+2] - vertices[3*i3+2];
      x = dy1*dz2 - dy2*dz1;
      y = dz1*dx2 - dz2*dx1;
      z = dx1*dy2 - dx2*dy1;
      vnorm = sqrt(x*x + y*y + z*z);
      if (vnorm == 0.0) {
         printf("Triangle %d has area equal to zero.\n", k);
         trnormals[3*k] = 0.0;
         trnormals[3*k+1] = 0.0;
         trnormals[3*k+2] = 0.0;
      } else {
         trnormals[3*k] = x/vnorm;
         trnormals[3*k+1] = y/vnorm;
         trnormals[3*k+2] = z/vnorm;
      }
      trareas[k] = 0.5*vnorm;
   }
   if (!have_vnormals) {
/*
 *  Compute unnormalized vertex normals by accumulating
 *  contributions from the triangles:  initialize normals
 *  to zero vectors and make one pass through the triangles.
 */
      for (i = 0; i < nv; i++) {
         vnormals[3*i] = 0.0;
         vnormals[3*i+1] = 0.0;
         vnormals[3*i+2] = 0.0;
      }
      for (k = 0; k < nt; k++) {
         i1 = indices[3*k];
         i2 = indices[3*k+1];
         i3 = indices[3*k+2];
         dx1 = vertices[3*i3] - vertices[3*i2];
         dy1 = vertices[3*i3+1] - vertices[3*i2+1];
         dz1 = vertices[3*i3+2] - vertices[3*i2+2];
         dx2 = vertices[3*i1] - vertices[3*i3];
         dy2 = vertices[3*i1+1] - vertices[3*i3+1];
         dz2 = vertices[3*i1+2] - vertices[3*i3+2];
         dx3 = vertices[3*i2] - vertices[3*i1];
         dy3 = vertices[3*i2+1] - vertices[3*i1+1];
         dz3 = vertices[3*i2+2] - vertices[3*i1+2];
         vn1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
         vn2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
         vn3 = dx3*dx3 + dy3*dy3 + dz3*dz3;
         if (vn1 == 0.0 || vn2 == 0.0 || vn3 == 0.0) {
            printf("Triangle %d has duplicate vertices.\n", k);
            exit(1);
         }
         vn1 = sqrt(vn1);
         vn2 = sqrt(vn2);
         vn3 = sqrt(vn3);
         a1 = acos(-(dx3*dx2+dy3*dy2+dz3*dz2)/(vn3*vn2));
         a2 = acos(-(dx1*dx3+dy1*dy3+dz1*dz3)/(vn1*vn3));
         a3 = acos(-(dx2*dx1+dy2*dy1+dz2*dz1)/(vn2*vn1));
         x = trnormals[3*k];
         y = trnormals[3*k+1];
         z = trnormals[3*k+2];
         vnormals[3*i1] += a1*x;      /* x-components */
         vnormals[3*i2] += a2*x;
         vnormals[3*i3] += a3*x;
         vnormals[3*i1+1] += a1*y;    /* y-components */
         vnormals[3*i2+1] += a2*y;
         vnormals[3*i3+1] += a3*y;
         vnormals[3*i1+2] += a1*z;    /* z-components */
         vnormals[3*i2+2] += a2*z;
         vnormals[3*i3+2] += a3*z;
      }
/*
 *  Normalize normals.
 */
      for (i = 0; i < nv; i++) {
         x = vnormals[3*i];
         y = vnormals[3*i+1];
         z = vnormals[3*i+2];
         vnorm = sqrt(x*x + y*y + z*z);
         if (vnorm != 0.0) {
            vnormals[3*i] /= vnorm;
            vnormals[3*i+1] /= vnorm;
            vnormals[3*i+2] /= vnorm;
         }
         else printf("Vertex %d has a zero normal vector.\n", i);
      }
   }
/*
 *  Compute approximations to vertex mean curvature values:
 *  H(v) = .5*<N(v),L(v)>, where v is a vertex, H denotes
 *  mean curvature, N(v) is the unit normal vector at v, and
 *  L is an approximation to the Laplace-Beltrami operator:
 *  L(p) = (1.5/A)*Sum[(cot a(i) + cot b(i))*(u(i)-v)],
 *  where the sum is over all neighbors u(i) of v (elements
 *  of the 1-ring), a(i) = Angle(p,q(i-1),q(i)), b(i) =
 *  Angle(p,q(i+1),q(i)), and A = denotes the sum of the areas
 *  of the triangles that contain p.
 *
 *  The vertex curvatures are initialized to zeros, the
 *  triangle contributions to all vertices are accumulated
 *  in a single pass through the triangles, and the curva-
 *  tures are normalized in a final loop on vertices.
 */
   for (i = 0; i < nv; i++) {
      vareas[i] = 0.0;
      vcrvs[i] = 0.0;
   }
   for (k = 0; k < nt; k++) {
      i1 = indices[3*k];
      i2 = indices[3*k+1];
      i3 = indices[3*k+2];
      dx1 = vertices[3*i3] - vertices[3*i2];
      dy1 = vertices[3*i3+1] - vertices[3*i2+1];
      dz1 = vertices[3*i3+2] - vertices[3*i2+2];
      dx2 = vertices[3*i1] - vertices[3*i3];
      dy2 = vertices[3*i1+1] - vertices[3*i3+1];
      dz2 = vertices[3*i1+2] - vertices[3*i3+2];
      dx3 = vertices[3*i2] - vertices[3*i1];
      dy3 = vertices[3*i2+1] - vertices[3*i1+1];
      dz3 = vertices[3*i2+2] - vertices[3*i1+2];
      vnorm = trareas[k];
      a1 = -0.5*(dx3*dx2+dy3*dy2+dz3*dz2)/vnorm;
      a2 = -0.5*(dx1*dx3+dy1*dy3+dz1*dz3)/vnorm;
      a3 = -0.5*(dx2*dx1+dy2*dy1+dz2*dz1)/vnorm;
      s12 = vnormals[3*i1]*dx2+vnormals[3*i1+1]*dy2+
            vnormals[3*i1+2]*dz2;
      s13 = vnormals[3*i1]*dx3+vnormals[3*i1+1]*dy3+
            vnormals[3*i1+2]*dz3;
      s23 = vnormals[3*i2]*dx3+vnormals[3*i2+1]*dy3+
            vnormals[3*i2+2]*dz3;
      s21 = vnormals[3*i2]*dx1+vnormals[3*i2+1]*dy1+
            vnormals[3*i2+2]*dz1;
      s31 = vnormals[3*i3]*dx1+vnormals[3*i3+1]*dy1+
            vnormals[3*i3+2]*dz1;
      s32 = vnormals[3*i3]*dx2+vnormals[3*i3+1]*dy2+
            vnormals[3*i3+2]*dz2;
      vcrvs[i1] += a3*s13 - a2*s12;
      vcrvs[i2] += a1*s21 - a3*s23;
      vcrvs[i3] += a2*s32 - a1*s31;
      vnorm /= 3.0;
      vareas[i1] += vnorm;
      vareas[i2] += vnorm;
      vareas[i3] += vnorm;
   }
/*
 *  Scale curvatures values.
 */
   for (i = 0; i < nv; i++) {
      vnorm = 0.0;
      if (vareas[i] != 0.0) vnorm = 0.25/vareas[i];
      vcrvs[i] *= vnorm;
   }
/*
 *  Normalize curvatures to [0,1] for texture mapping.
 */
   normalizeVC();
   return;
}


/*
 *----------------------------------------------------------
 *
 *  csDisplay:  Color selection window display callback.
 *
 *  glsurf functions called:  displayString, drawRect,
 *                            rgb2hsv
 *  OpenGL functions called:  glBegin, glClear, glColor3fv,
 *                            glColor3ub, glDrawPixels,
 *                            glEnd, glRasterPos2i, glRecti,
 *                            glutSwapBuffers, glVertex2i
 *
 *----------------------------------------------------------
 */
static void csDisplay(void)
{
   GLfloat black[3] = {0.0, 0.0, 0.0};   /* Border colors */
   GLfloat white[3] = {1.0, 1.0, 1.0};
   int fontno = 6;    /* Font number for strings (18-point
                           proportional spaced Helvetica) */
   float h, s, v;     /* HSV components of currently
                           selected color */
   int j;             /* Column index for vtable */
   int jv;            /* Row index for vtable */
   GLfloat pos[3];    /* Text string position */
   GLubyte r, g, b;   /* RGB components of currently
                           selected color (cs_color) */
   char text[8];      /* Text string */
   float vt;          /* Intensity value indexed by jv
                           divided by v */
   int x0, y0;        /* Coordinates of the lower left
                           corner of a rectangle */
/*
 * Table of intensities associated with the current
 * hstable entry:
 */
   GLubyte vtable[NV][18][3];
/*
 *  Clear the color buffer to the background color.
 */
   glClear(GL_COLOR_BUFFER_BIT);
/*
 *  Copy hstable to the window and draw its outline.
 */
   x0 = rhstable.x0;
   y0 = rhstable.y0;
   glRasterPos2i(x0, y0);
   glDrawPixels(rhstable.w, rhstable.h, GL_RGB,
                GL_UNSIGNED_BYTE, hstable);
   drawRect(x0, y0, rhstable.w, rhstable.h, white,
            black);
/*
 *  Create vtable, copy it to the frame buffer, and draw its
 *  outline.
 */
   r = (GLubyte) (255.0*cs_color[0] + 0.5);
   g = (GLubyte) (255.0*cs_color[1] + 0.5);
   b = (GLubyte) (255.0*cs_color[2] + 0.5);
   rgb2hsv(cs_color[0], cs_color[1], cs_color[2], &h, &s, &v);

   for (jv = 0; jv < NV; jv++) {
      vt = ((float) jv)/((float) (NV-1));
      if (v != 0.0) vt = vt/v;
      for (j = 0; j < rvtable.w; j++) {
         vtable[jv][j][0] = (GLubyte) (vt*r + 0.5);
         vtable[jv][j][1] = (GLubyte) (vt*g + 0.5);
         vtable[jv][j][2] = (GLubyte) (vt*b + 0.5);
      }
   }
   x0 = rvtable.x0;
   y0 = rvtable.y0;
   glRasterPos2i(x0, y0);
   glDrawPixels((GLint) rvtable.w, (GLint) rvtable.h, GL_RGB,
                GL_UNSIGNED_BYTE, vtable);
   drawRect(x0, y0, rvtable.w, rvtable.h, white,
            black);
/*
 *  Draw a triangular pointer (to the current intensity
 *  value) on the right edge of the vtable rectangle.
 */
   jv = (int) (((float) (NV-1))*v);
   x0 += rvtable.w+2;
   y0 += jv;
   glColor3fv(black);
   glBegin(GL_TRIANGLES);
      glVertex2i(x0, y0);
      glVertex2i(x0+5, y0-5);
      glVertex2i(x0+5, y0+5);
   glEnd();
/*
 *  Draw a rectangle filled with the currently selected
 *  color.
 */
   glColor3ub(r, g, b);
   x0 = rcolor.x0;
   y0 = rcolor.y0;
   glRecti(x0, y0, x0+rcolor.w-1, y0+rcolor.h-1);
   drawRect(x0, y0, rcolor.w, rcolor.h, white, black);
/*
 *  Draw rectangles and text strings specifying the
 *  currently selected color components.
 */
   glColor3fv(black);
   x0 = rh.x0;
   y0 = rh.y0;
   glRecti(x0, y0, x0+rh.w-1, y0+rh.h-1);
   drawRect(x0, y0, rh.w, rh.h, white, black);
   pos[0] = 108.0;
   pos[1] = 180.0;
   pos[2] = 0.0;
   sprintf(text, "H:  %3u", (int) h);
   glColor3fv(white);
   displayString(pos, text, fontno);

   glColor3fv(black);
   x0 = rr.x0;
   y0 = rr.y0;
   glRecti(x0, y0, x0+rr.w-1, y0+rr.h-1);
   drawRect(x0, y0, rr.w, rr.h, white, black);
   pos[0] = 203.0;
   sprintf(text, "R:  %3u", r);
   glColor3fv(white);
   displayString(pos, text, fontno);

   glColor3fv(black);
   x0 = rs.x0;
   y0 = rs.y0;
   glRecti(x0, y0, x0+rs.w-1, y0+rs.h-1);
   drawRect(x0, y0, rs.w, rs.h, white, black);
   pos[0] = 108.0;
   pos[1] = 128.0;
   sprintf(text, "S:  %3u", (int) (255.0*s + 0.5));
   glColor3fv(white);
   displayString(pos, text, fontno);

   glColor3fv(black);
   x0 = rg.x0;
   y0 = rg.y0;
   glRecti(x0, y0, x0+rg.w-1, y0+rg.h-1);
   drawRect(x0, y0, rg.w, rg.h, white, black);
   pos[0] = 203.0;
   sprintf(text, "G:  %3u", g);
   glColor3fv(white);
   displayString(pos, text, fontno);

   glColor3fv(black);
   x0 = rv.x0;
   y0 = rv.y0;
   glRecti(x0, y0, x0+rv.w-1, y0+rv.h-1);
   drawRect(x0, y0, rv.w, rv.h, white, black);
   pos[0] = 108.0;
   pos[1] = 75.0;
   sprintf(text, "V:  %3u", (int) (255.0*v + 0.5));
   glColor3fv(white);
   displayString(pos, text, fontno);

   glColor3fv(black);
   x0 = rb.x0;
   y0 = rb.y0;
   glRecti(x0, y0, x0+rb.w-1, y0+rb.h-1);
   drawRect(x0, y0, rb.w, rb.h, white, black);
   pos[0] = 203.0;
   sprintf(text, "B:  %3u", b);
   glColor3fv(white);
   displayString(pos, text, fontno);
/*
 * Draw buttons at the bottom of the window.
 */
   x0 = raccept.x0;
   y0 = raccept.y0;
   drawRect(x0, y0, raccept.w, raccept.h, black, white);
   pos[0] = (float) (x0+7);
   pos[1] = (float) (y0+7);
   sprintf(text, "ACCEPT");
   glColor3fv(white);
   displayString(pos, text, fontno);

   x0 = rcancel.x0;
   y0 = rcancel.y0;
   drawRect(x0, y0, rcancel.w, rcancel.h, black, white);
   pos[0] = (float) (x0+7);
   pos[1] = (float) (y0+7);
   sprintf(text, "CANCEL");
   glColor3fv(white);
   displayString(pos, text, fontno);

   glutSwapBuffers();
   return;
}


/*
 *----------------------------------------------------------
 *
 *  csInit:  Color selection window initialization:
 *           construct the Hue/Saturation table hstable
 *           (using v = 1).
 *
 *  glsurf function called:  hsv2rgb
 *  OpenGL function called:  glPixelStorei
 *
 *----------------------------------------------------------
 */
static void csInit(void)
{
   int nh = NH;     /* Number of columns (hues) */
   int ns = NS;     /* Number of rows (saturation values) */
   int jh, js;      /* Column and row indices */
   float h, s, v;   /* HSV color components */
   float r, g, b;   /* RGB color components */

   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
   v = 1.0;
   for (js = 0; js < ns; js++) {
      s = ((float) js)/((float) (ns-1));
      for (jh = 0; jh < nh; jh++) {
         h = 360.0*((float) jh)/((float) (nh-1));
         hsv2rgb(h, s, v, &r, &g, &b);
         hstable[js][jh][0] = (GLubyte) (255.0*r + 0.5);
         hstable[js][jh][1] = (GLubyte) (255.0*g + 0.5);
         hstable[js][jh][2] = (GLubyte) (255.0*b + 0.5);
      }
   }
   return;
}


/*
 *----------------------------------------------------------
 *
 *  csKey:  Keyboard callback associated with the color
 *          selection window:
 *          Enter ==> Accept current color,
 *          Escape ==> Cancel selection window.
 *
 *  The mouse coordinates (x,y) are not used.
 *
 *  glsurf function called:  csMouse
 *  OpenGL functions called:  glutHideWindow, glutSetWindow
 *
 *----------------------------------------------------------
 */
static void csKey(unsigned char key, int x, int y)
{
   if (key == 10  ||  key == 13) {  /* Enter key pressed */
      csMouse(GLUT_LEFT_BUTTON, GLUT_UP,
              raccept.x0+1, 398-raccept.y0);
   } else {
   if (key == 27) {                 /* Escape key pressed */
      glutHideWindow();
      glutSetWindow(1);
   }}
   return;
}


/*
 *----------------------------------------------------------
 *
 *  csMotion:  Callback triggered by mouse motion in the
 *             color selection window with one or more mouse
 *             buttons pressed.
 *
 *  This function updates cs_color if the mouse position
 *  is inside either of the rectangles rhstable or rvtable.
 *
 *  glsurf function called:  fmax
 *  OpenGL function called:  glutPostRedisplay
 *
 *----------------------------------------------------------
 */
static void csMotion(int x, int y)
{
   int jh, js;       /* New hstable indices for (x,y) */
   int jv;           /* New intensity index for y */
   GLubyte r, g, b;  /* RGB components of new hstable entry */
   float v;          /* New intensity value divided by v0 */
   float v0;         /* Current or previous intensity value */
   int x0, y0;       /* Lower left corner of a rectangle */
   int x1, y1;       /* Upper right corner of a rectangle */

   y = 399 - y;  /* Map glut y-coordinate to viewport */

   x0 = rhstable.x0;
   y0 = rhstable.y0;
   x1 = x0 + rhstable.w - 1;
   y1 = y0 + rhstable.h - 1;
   if (x >= x0  &&  x <= x1  &&  y >= y0  &&  y <= y1) {
      jh = x - x0;
      js = y - y0;
      v0 = fmax(cs_color[0], cs_color[1], cs_color[2]);
      r = hstable[js][jh][0]; /* Current color with v = 1 */
      g = hstable[js][jh][1];
      b = hstable[js][jh][2];
      cs_color[0] = v0*((float) r)/255.0;
      cs_color[1] = v0*((float) g)/255.0;
      cs_color[2] = v0*((float) b)/255.0;
      glutPostRedisplay();
      return;
   }

   x0 = rvtable.x0;
   y0 = rvtable.y0;
   x1 = x0 + rvtable.w - 1;
   y1 = y0 + rvtable.h - 1;
   if (x >= x0  &&  x <= x1  &&  y >= y0  &&  y <= y1) {
      jv = y - y0;
      v0 = fmax(cs_color[0], cs_color[1], cs_color[2]);
      v = ((float) jv)/((float) (NV-1));
      if (v0 != 0.0) {
         v = v/v0;
         cs_color[0] *= v;
         cs_color[1] *= v;
         cs_color[2] *= v;
      } else {
         cs_color[0] = v;
         cs_color[1] = v;
         cs_color[2] = v;
      } 
      glutPostRedisplay();
   }
   return;
}


/*
 *----------------------------------------------------------
 *
 *  csMouse:  Callback triggered by a mouse button press or
 *            release with the mouse position in the color
 *            selection window.
 *
 *  If the left or right mouse button is pressed with the
 *  mouse position in one of the six color specification
 *  rectangles, the appropriate cs_color component is
 *  incremented (left button) or decremented (right button).
 *
 *  If the left mouse button is pressed with the mouse
 *  position inside the 'Accept' or 'Cancel' rectangle, the
 *  appearance of the rectangle is altered.  If the left
 *  button is released inside one of these rectangles, the
 *  window state is changed to hidden (after updating
 *  top_color to cs_color in the case of 'Accept').
 *
 *  If none of the above events occur, csMotion is called so
 *  that a button press in the rhtable or rvtable rectangle
 *  changes the current color.
 *
 *  glsurf functions called:  csMotion, drawRect, hsv2rgb,
 *                            rgb2hsv
 *  OpenGL functions called:  glMaterialfv, glutHideWindow,
 *                            glutPostRedisplay,
 *                            glutSetWindow
 *
 *----------------------------------------------------------
 */
static void csMouse(int button, int state, int x, int y)
{
   GLfloat black[3] = {0.0, 0.0, 0.0};   /* Border colors */
   GLfloat white[3] = {1.0, 1.0, 1.0};
   float h, s, v;     /* HSV components of currently
                           selected color */
   float hstep = 1.0;       /* Increment step for h */
   float step = 1.0/255.0;  /* Increment step for s, v,
                                 r, g, and b */
   int x0, y0;   /* Lower left corner of a rectangle */
   int x1, y1;   /* Upper right corner of a rectangle */

   y = 399 - y;  /* Map glut y-coordinate to viewport */

   if (state == GLUT_DOWN && button != GLUT_MIDDLE_BUTTON) {
/*
 *  Left or right button press.
 */
      rgb2hsv(cs_color[0], cs_color[1], cs_color[2], &h, &s, &v);

      x0 = rh.x0;
      y0 = rh.y0;
      x1 = x0 + rh.w - 1;
      y1 = y0 + rh.h - 1;
      if (x > x0  &&  x < x1  &&  y > y0  &&  y < y1) {
         if (button == GLUT_LEFT_BUTTON) {
            h += hstep;
            if (h >= 360.0) h -= 360.0;
         } else {
            h -= hstep;
            if (h < 0.0) h += 360.0;
         }
         hsv2rgb(h, s, v, &cs_color[0], &cs_color[1], &cs_color[2]);
         glutPostRedisplay();
         return;
      }
      x0 = rs.x0;
      y0 = rs.y0;
      x1 = x0 + rs.w - 1;
      y1 = y0 + rs.h - 1;
      if (x > x0  &&  x < x1  &&  y > y0  &&  y < y1) {
         if (button == GLUT_LEFT_BUTTON) {
            s += step;
            if (s > 1.0) s = 1.0;
         } else {
            s -= step;
            if (s < 0.0) s = 0.0;
         }
         hsv2rgb(h, s, v, &cs_color[0], &cs_color[1], &cs_color[2]);
         glutPostRedisplay();
         return;
      }
      x0 = rv.x0;
      y0 = rv.y0;
      x1 = x0 + rv.w - 1;
      y1 = y0 + rv.h - 1;
      if (x > x0  &&  x < x1  &&  y > y0  &&  y < y1) {
         if (button == GLUT_LEFT_BUTTON) {
            v += step;
            if (v > 1.0) v = 1.0;
         } else {
            v -= step;
            if (v < 0.0) v = 0.0;
         }
         hsv2rgb(h, s, v, &cs_color[0], &cs_color[1], &cs_color[2]);
         glutPostRedisplay();
         return;
      }
      x0 = rr.x0;
      y0 = rr.y0;
      x1 = x0 + rr.w - 1;
      y1 = y0 + rr.h - 1;
      if (x > x0  &&  x < x1  &&  y > y0  &&  y < y1) {
         if (button == GLUT_LEFT_BUTTON) {
            cs_color[0] += step;
            if (cs_color[0] > 1.0) cs_color[0] = 1.0;
         } else {
            cs_color[0] -= step;
            if (cs_color[0] < 0.0) cs_color[0] = 0.0;
         }
         glutPostRedisplay();
         return;
      }
      x0 = rg.x0;
      y0 = rg.y0;
      x1 = x0 + rg.w - 1;
      y1 = y0 + rg.h - 1;
      if (x > x0  &&  x < x1  &&  y > y0  &&  y < y1) {
         if (button == GLUT_LEFT_BUTTON) {
            cs_color[1] += step;
            if (cs_color[1] > 1.0) cs_color[1] = 1.0;
         } else {
            cs_color[1] -= step;
            if (cs_color[1] < 0.0) cs_color[1] = 0.0;
         }
         glutPostRedisplay();
         return;
      }
      x0 = rb.x0;
      y0 = rb.y0;
      x1 = x0 + rb.w - 1;
      y1 = y0 + rb.h - 1;
      if (x > x0  &&  x < x1  &&  y > y0  &&  y < y1) {
         if (button == GLUT_LEFT_BUTTON) {
            cs_color[2] += step;
            if (cs_color[2] > 1.0) cs_color[2] = 1.0;
         } else {
            cs_color[2] -= step;
            if (cs_color[2] < 0.0) cs_color[2] = 0.0;
         }
         glutPostRedisplay();
         return;
      }
   }

   if (button != GLUT_LEFT_BUTTON) {
      csMotion(x, 399-y);
      return;
   }
/*
 *  Left button press or release.
 */
   x0 = raccept.x0;
   y0 = raccept.y0;
   x1 = x0 + raccept.w - 1;
   y1 = y0 + raccept.h - 1;
   if (x > x0  &&  x < x1  &&  y > y0  &&  y < y1) {
/*
 *  'Accept' button.
 */
      if (state == GLUT_DOWN) {
         drawRect(x0, y0, raccept.w, raccept.h, white,
                  black);
         glutSwapBuffers();
      } else {                 /* 'Accept' button released */
         top_color[0] = cs_color[0];
         top_color[1] = cs_color[1];
         top_color[2] = cs_color[2];
         glutHideWindow();
         glutSetWindow(1);
         glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,
                      top_color);
      }
      return;
   }

   x0 = rcancel.x0;
   y0 = rcancel.y0;
   x1 = x0 + rcancel.w - 1;
   y1 = y0 + rcancel.h - 1;
   if (x > x0  &&  x < x1  &&  y > y0  &&  y < y1) {
/*
 *  'Cancel' button.
 */
      if (state == GLUT_DOWN) {
         drawRect(x0, y0, rcancel.w, rcancel.h, white,
                  black);
         glutSwapBuffers();
      } else {                /* 'Cancel' button released */
         glutHideWindow();
         glutSetWindow(1);
      }
      return;
   }
   csMotion(x, 399-y);
   return;
}


/*
 *----------------------------------------------------------
 *
 *  csReshape:  Set the viewport, projection matrix, and
 *              modelview matrix for the color selection
 *              window (window 2):  two-dimensional render-
 *              ing with a 300 by 400 viewport.
 *
 *  This function is called before the first call to
 *  csDisplay.
 *
 *  If the window is too small, a request for a larger
 *  window is executed.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  glLoadIdentity, glMatrixMode,
 *                            glScalef, glTranslatef,
 *                            gluOrtho2D, glutReshapeWindow,
 *                            glViewport
 *
 *----------------------------------------------------------
 */
static void csReshape(int w, int h)
{
   if (w < 300  ||  h < 400) {     /* Test for sufficient */
      glutReshapeWindow(300, 400); /*   window extent.    */
      return;
   }
/*
 *  Set the viewport to the entire window (the default).
 */
   glViewport(0, 0, (GLint) w, (GLint) h);
/*
 *  Set up an orthogonal projection for data space equal to
 *  viewport (identity).
 */
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluOrtho2D(0.0, (GLfloat) w, 0.0, (GLfloat) h);
/*
 *  Set the modelview matrix to the identity.
 */
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();

   return;
}


/*
 *----------------------------------------------------------
 *
 *  display:  Clear the color and depth buffers, construct
 *            the modelview matrix, and draw the surface,
 *            optionally including a contour plot, space
 *            curve, title, and the bounding box.
 *
 *  glsurf functions called:  cFillTriangle, displayString
 *  OpenGL functions called:  glBegin, glBindTexture,
 *                            glClear, glColor3f,
 *                            glColor3fv, glDisable,
 *                            glDisableClientState,
 *                            glDrawElements, glEnable,
 *                            glEnableClientState, glEnd,
 *                            glLineWidth, glLoadIdentity,
 *                            glMatrixMode, glNormal3dv,
 *                            glPolygonOffset, glRotated,
 *                            glTranslated, glutSwapBuffers,
 *                            glVertex3d, glVertex3dv
 *
 *----------------------------------------------------------
 */
static void display(void)
{
   int k;   /* Triangle index or space curve vertex index */
/*
 *  Clear the color and depth buffers.
 */
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
/*
 *  Set up the modelview matrix:
 *    Translate by -c, rotate about view_axis,
 *    translate by -3*r*e3, and rotate about orient_axis.
 */
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glRotated(orient_angle, orient_axis[0], orient_axis[1],
             orient_axis[2]);
   glTranslated(0.0, 0.0, -3.0*r);
   glRotated(view_angle, view_axis[0], view_axis[1],
             view_axis[2]);
   glTranslated(-xc, -yc, -zc);
/*
 *  Make texture object tex_flag active, or disable texture
 *  mapping if tex_flag = 0.
 */
   if (tex_flag) {
      glEnable(GL_TEXTURE_1D);
      glBindTexture(GL_TEXTURE_1D, tex_flag);
      if (tex_flag <= 2) {
/*
 *  Automatically generate a texture s-coordinate (for each
 *  vertex) based on distance from the x-y plane (z value
 *  in object coordinates).
 */
         glDisableClientState(GL_TEXTURE_COORD_ARRAY);
         glEnable(GL_TEXTURE_GEN_S);
      } else {
/*
 *  Take the source of the texture data (for glDrawElements)
 *  to be the array specified in init1.
 */
         glDisable(GL_TEXTURE_GEN_S);
         glEnableClientState(GL_TEXTURE_COORD_ARRAY);
      }
   } else {
      glDisable(GL_TEXTURE_1D);
   }
   if (!polyfill) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
   if (smooth) {
/*
 *  Draw the surface using vertex normals.
 */
      glDrawElements(GL_TRIANGLES, 3*nt, GL_UNSIGNED_INT,
                     indices);
   } else {
/*
 *  Draw the surface using triangle normals.
 */
      glBegin(GL_TRIANGLES);
      for (k = 0; k < nt; k++){
         glNormal3dv(&trnormals[3*k]);
         glVertex3dv(&vertices[3*indices[3*k]]);
         glVertex3dv(&vertices[3*indices[3*k+1]]);
         glVertex3dv(&vertices[3*indices[3*k+2]]);
      }
      glEnd();
   }
   glDisable(GL_TEXTURE_1D);
   if (!polyfill) {
/*
 *  Wireframe mode:  restore polygon fill mode, enable
 *  polygon offsetting, remove hidden lines by drawing
 *  the surface again with the background color, and
 *  disable offsetting.
 */
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1.0, 1.0);
      glDisable(GL_LIGHTING);
      glColor3f(0.0, 0.0, 0.0);
      glDrawElements(GL_TRIANGLES, 3*nt, GL_UNSIGNED_INT,
                     indices);
      if (lighting) glEnable(GL_LIGHTING);
      glDisable(GL_POLYGON_OFFSET_FILL);
   }

   if (cplot) {
/*
 * Display a contour plot
 */
      glDisable(GL_LIGHTING);
      for (k = 0; k < nt; k++) {
         cFillTriangle(indices[3*k], indices[3*k+1],
                       indices[3*k+2], zcon);
      }
      if (lighting) glEnable(GL_LIGHTING);
   }

   if (np > 0) {
/*
 *  Draw a space curve.
 */
      glDisable(GL_LIGHTING);
      if (!hide_curve) {
         glDisable(GL_DEPTH_TEST);
      }
      glLineWidth(line_width);
      glColor3fv(curve_color);
      glBegin(GL_LINE_STRIP);
      for (k = 0; k < np; k++) {
         glVertex3dv(&cpoints[3*k]);
      }
      glEnd();
      glLineWidth(1.0);
      if (lighting) glEnable(GL_LIGHTING);
      glEnable(GL_DEPTH_TEST);
   }

   if (bbox) {
/*
 * Display the bounding box
 */
      unsigned int fontno = 5;    /* Font for labels */
      GLfloat pos[3];             /* Label position */
      char text[30];              /* Label text */

      glDisable(GL_LIGHTING);
      glColor3f(0.7, 0.7, 0.7);   /* grey */
      glBegin(GL_LINE_LOOP);
        glVertex3d(xmin, ymin, zmin);
        glVertex3d(xmax, ymin, zmin);
        glVertex3d(xmax, ymax, zmin);
        glVertex3d(xmin, ymax, zmin);
      glEnd();
      glBegin(GL_LINE_LOOP);
        glVertex3d(xmin, ymin, zmax);
        glVertex3d(xmax, ymin, zmax);
        glVertex3d(xmax, ymax, zmax);
        glVertex3d(xmin, ymax, zmax);
      glEnd();
      glBegin(GL_LINES);
        glVertex3d(xmin, ymin, zmin);
        glVertex3d(xmin, ymin, zmax);
        glVertex3d(xmax, ymin, zmin);
        glVertex3d(xmax, ymin, zmax);
        glVertex3d(xmax, ymax, zmin);
        glVertex3d(xmax, ymax, zmax);
        glVertex3d(xmin, ymax, zmin);
        glVertex3d(xmin, ymax, zmax);
      glEnd();

      pos[0] = (GLfloat) xmin;       /* Label two corners */
      pos[1] = (GLfloat) ymin;
      pos[2] = (GLfloat) (zmin + 0.04*(zmax-zmin));
      sprintf(text, "(%3.2f, %3.2f, %3.2f)", xmin, ymin, dmin);
      glColor3f(1.0, 1.0, 1.0);   /* white */
      displayString(pos, text, fontno);
      pos[0] = (GLfloat) xmax;
      pos[1] = (GLfloat) ymax;
      pos[2] = (GLfloat) (zmax + 0.04*(zmax-zmin));
      sprintf(text, "(%3.2f, %3.2f, %3.2f)", xmax, ymax, dmax);
      displayString(pos, text, fontno);
      if (lighting) glEnable(GL_LIGHTING);
   }

   if (*title != '\0') {
/*
 * Display the title string
 */
      glDisable(GL_DEPTH_TEST);
      glDisable(GL_LIGHTING);
      displayString(title_pos, title, font_number);
      if (lighting) glEnable(GL_LIGHTING);
      glEnable(GL_DEPTH_TEST);
   }

   glutSwapBuffers();

   return;
}


/*
 *----------------------------------------------------------
 *
 *  displayString:  Writes string to the screen with lower
 *                  left corner at raster position pos
 *                  using font number fontno (defined in
 *                  the table below).
 *
 *  Note that the raster position is actually a location in
 *  the user coordinate space; i.e., the current modelview
 *  and projection matrices are applied to pos in order to
 *  obtain the screen coordinates of the first character.
 *  If those coordinates are not in the view volume, nothing
 *  is written.  Otherwise the string is written but clipped
 *  against the window boundaries.
 *
 *  The string is written in the currently selected color.
 *
 *  If fontno is invalid, an error message is written to
 *  the standard output unit and the return value is 1.
 *  Otherwise the return value is 0.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  glLineWidth, glPopMatrix,
 *                            glPushMatrix, glRasterPos3fv,
 *                            glScalef, glTranslatef,
 *                            glutBitmapCharacter,
 *                            glutStrokeCharacter
 *
 *----------------------------------------------------------
 */
static int displayString(GLfloat pos[3], char *string,
                         unsigned int fontno)
{
   unsigned int i;            /* string index */
   unsigned int nfonts = 9;   /* Number of fonts */
   GLfloat sf;                /* Scale factor for stroked fonts */

   static void *font[9] = {
      GLUT_BITMAP_8_BY_13,         /* 8 x 13 fixed width */
      GLUT_BITMAP_9_BY_15,         /* 9 x 15 fixed width */  
      GLUT_BITMAP_TIMES_ROMAN_10,  /* 10-point proportional spaced */
      GLUT_BITMAP_TIMES_ROMAN_24,  /* 24-point proportional spaced */
      GLUT_BITMAP_HELVETICA_10,    /* 10-point proportional spaced */
      GLUT_BITMAP_HELVETICA_12,    /* 12-point proportional spaced */
      GLUT_BITMAP_HELVETICA_18,    /* 18-point proportional spaced */
      GLUT_STROKE_ROMAN,           /* proportional spaced stroke */
      GLUT_STROKE_MONO_ROMAN       /* mono-spaced Roman simplex */
   };

   if (fontno >= nfonts) {
      printf("displayString:  Invalid font number = %u\n",
             fontno);
      return 1;
   }
   if (fontno < 7 ) {
      glRasterPos3fv(pos);
      i = 0;
      while (string[i] != '\0') {
         glutBitmapCharacter(font[fontno], string[i]);
         i++;
      }
   } else {
      sf = 0.0005*r;
      glLineWidth(3.0);
      glPushMatrix();
      glTranslatef(pos[0], pos[1], pos[2]);
      glScalef(sf, sf, sf);
      i = 0;
      while (string[i] != '\0') {
         glutStrokeCharacter(font[fontno], string[i]);
         i++;
      }
      glPopMatrix();
      glLineWidth(1.0);
   }
   return 0;
}


/*
 *----------------------------------------------------------
 *
 *  drawRect:  Draws the outline of an axis-alligned rectan-
 *             gle in the z = 0 plane with lower left corner
 *             at (x0,y0), width w, and height h.  The origin
 *             of the coordinate system is at the lower left.
 *
 *  The bottom (lower) and right edges are drawn with color
 *  brcolor; the top and left edges are drawn with color
 *  tlcolor, where brcolor and tlcolor are arrays of length
 *  three containing red, green, and blue components normal-
 *  ized to [0,1].  Taking brcolor and tlcolor to be black
 *  (0,0,0) and white (1,1,1), respectively, gives the
 *  appearance of a button sticking out of the screen.
 *  Reversing those colors makes the button appear to extend
 *  inward as if it were pressed.
 *
 *  The outline is drawn with lines two pixels wide.  If the
 *  rectangle is to be filled, it should be done before the
 *  outline is drawn (or only the interior should be filled)
 *  in order to preserve the outline colors.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  glBegin, glColor3fv, glEnd,
 *                            glLineWidth, glVertex2i
 *
 *----------------------------------------------------------
 */
static void drawRect(GLuint x0, GLuint y0, GLuint w, GLuint h,
                     GLfloat brcolor[3], GLfloat tlcolor[3])
{
   glLineWidth(2.0);
   glBegin(GL_LINES);
     glColor3fv(brcolor);
     glVertex2i(x0, y0);            /* Bottom */
     glVertex2i(x0+w-1, y0);

     glVertex2i(x0+w-1, y0);        /* Right */
     glVertex2i(x0+w-1, y0+h-1);

     glColor3fv(tlcolor);
     glVertex2i(x0+w-1, y0+h-1);    /* Top */
     glVertex2i(x0, y0+h-1);

     glVertex2i(x0, y0+h-1);        /* Left */
     glVertex2i(x0, y0);
   glEnd();
   glLineWidth(1.0);
   return;
}


/*
 *----------------------------------------------------------
 *
 *  fmax:  Returns max{x, y, z}.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  None
 *
 *----------------------------------------------------------
 */
static float fmax(float x, float y, float z)
{
   float fm;

   if (x >= y) {
      fm = z;
      if (x >= z) fm = x;
   } else {
      fm = y;
      if (x < z) fm = (y >= z) ? y : z;
   }
   return fm;
}


/*
 *----------------------------------------------------------
 *
 *  getLin:  Read a line from stdin into s.
 *
 *  On input:  s is a char array of length at least len.
 *             Input is terminated by end of file, a new-
 *             line character, or when len-1 characters
 *             have been read.
 *
 *  On output:  s contains a (null-terminated) string of
 *              length L for return value L <= len.  The
 *              newline character is not included.  No
 *              characters are left in the input stream.
 *
 *  This function uses the constant EOF defined in <stdio.h>.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  None
 *
 *----------------------------------------------------------
 */
static int getLin(char *s, int len)
{
   int c, i;

   i = 0;
   c = getchar();
   while (i < len-1  &&  c != EOF  &&  c != '\n') {
      s[i] = c;
      i++;
      c = getchar();
   }
   s[i] = '\0';
   while (c != EOF && c != '\n') {
      c = getchar();                /* flush input stream */
   }
   return i+1;
}


/*
 *----------------------------------------------------------
 *
 *  getRasterPos:  Callback triggered by a mouse button
 *                 press or release with the mouse position
 *                 in the current window (except for buttons
 *                 attached to menus).
 *
 *  Set title_pos to the world coordinates (half way between
 *  the near and far clipping planes) corresponding to the
 *  mouse position when the left button is pressed.  The
 *  mappings defined by the modelview, projection, and
 *  viewport matrices must be inverted.
 * 
 *  glsurf functions called:  None
 *  OpenGL functions called:  glGetDoublev, glGetIntegerv,
 *                            glutGet, glutMouseFunc,
 *                            gluUnProject
 *
 *----------------------------------------------------------
 */
static void getRasterPos(int button, int state, int x,
                         int y)
{
   GLdouble mvmatrix[16];    /* Modelview matrix */
   GLdouble projmatrix[16];  /* Projection matrix */
   GLint viewport[4];
   GLdouble wx, wy, wz;      /* Computed world coordinates */
   GLdouble yr;              /* OpenGL y coordinate */

   if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
      glGetIntegerv(GL_VIEWPORT, viewport);
      glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
      glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);
/*
 *  Compute yr = h-1-y for window height h.
 */
      yr = glutGet(GLUT_WINDOW_HEIGHT) - 1 - (GLint) y;
      gluUnProject((GLdouble) x, (GLdouble) yr, 0.5, mvmatrix,
                   projmatrix, viewport, &wx, &wy, &wz);
      title_pos[0] = (GLfloat) wx;
      title_pos[1] = (GLfloat) wy;
      title_pos[2] = (GLfloat) wz;
      glutMouseFunc(mouse);  /* Restore standard callback */
      glutPostRedisplay();
   }
   return;
}


/*
 *----------------------------------------------------------
 *
 *  hsv2rgb:  Converts a color specified in the HSV system
 *            to the RGB system.
 *
 *  On input:  h = Hue in the range 0 to 360 (degrees),
 *                 where red, green, and blue correspond
 *                 to h = 0, h = 120, and h = 240, respec-
 *                 tively.
 *
 *             s = Saturation value in [0,1], where s = 0
 *                 implies greyscale (h is irrelevant and
 *                 r = g = b = v), and s = 1 implies a
 *                 pure hue (at least one of r, g, and b
 *                 is 0).
 *
 *             v = Value of intensity in [0,1], where
 *                 v = 0 implies black, and v = 1 implies
 *                 the maximum brightness.
 *
 *  A multiple of 360 is added to h if necessary, and the
 *  values of s and v are clamped to the interval [0,1].
 *
 *  On output:  r,g,b = Red, green, and blue color compo-
 *                      nents in the range [0,1].
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  None
 *
 *----------------------------------------------------------
 */
static void hsv2rgb(float h, float s, float v, float *r,
                    float *g, float*b)
{
   int i;         /* Integer part of h/60 */
   float f;       /* Fractional part of h/60 */
   float p,q,t;   /* Scaled down intensities */

   while (h < 0.0) {   /* Adjust h if necessary */
      h += 360.0;
   }
   while (h >= 360.0) {
      h -= 360.0;
   }
   if (s < 0.0) {      /* Clamp s to [0,1] */
      s = 0.0;
   } else {
   if (s > 1.0) {
      s = 1.0;
   }}

   if (v < 0.0) {      /* Clamp v to [0,1] */
      v = 0.0;
   } else {
   if (v > 1.0) {
      v = 1.0;
   }}

   h /= 60.0;     /* Map h to [0,6), and compute i and f */
   i = (int) h;
   f = h - (float) i;

   p = v*(1.0 - f*s);        /* Compute p, q, and t */
   q = v*(1.0 - (1.0-f)*s);
   t = v*(1.0 - s);
   switch (i) {
   case 0:
      *r = v;  *g = q;  *b = t;
      break;
   case 1:
      *r = p;  *g = v;  *b = t;
      break;
   case 2:
      *r = t;  *g = v;  *b = q;
      break;
   case 3:
      *r = t;  *g = p;  *b = v;
      break;
   case 4:
      *r = q;  *g = t;  *b = v;
      break;
   case 5:
      *r = v;  *g = t;  *b = p;
      break;
   default:
      break;
   }
   return;
}


/*
 *----------------------------------------------------------
 *
 *  init1:  Initialization of window 1.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  glBindTexture, glClearColor,
 *                            glEnable, glEnableClientState,
 *                            glGenTextures, glLightModeli,
 *                            glLightfv, glMaterialf,
 *                            glMaterialfv, glNormalPointer,
 *                            glPixelStorei,
 *                            glTexCoordPointer, glTexEnvf,
 *                            glTexGenfv, glTexGeni,
 *                            glTexImage1D, glTexParameteri,
 *                            glVertexPointer
 *
 *----------------------------------------------------------
 */
static void init1(void)
{
   GLfloat ambient_light[4] = {0.4, 0.4, 0.4, 1.0};
   GLfloat light_position[4] = {0.0, 0.0, 1.0, 0.0};
   GLuint tex_names[2];  /* Texture names (1 and 2) */
   GLfloat xyplane[4] = {0.0, 0.0, 1.0, 0.0};
/*
 *  Enable depth test and lighting.
 */
   glEnable(GL_DEPTH_TEST);
   if (lighting) glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
/*
 *  Set background, material, and light colors.  (Material
 *    colors and shininess are global variables.)
 */
   glClearColor (0.0, 0.0, 0.0, 0.0);
   
   glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,
                top_color);
   glMaterialfv(GL_FRONT, GL_SPECULAR, sp_color);
   glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE,
                bottom_color);
   glMaterialfv(GL_BACK, GL_SPECULAR, sp_color);
   glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS,
                shininess);

   glLightfv(GL_LIGHT0, GL_POSITION, light_position);
   glLightfv(GL_LIGHT0, GL_AMBIENT, ambient_light); 

/*
 *  Activate client-side vertex and normal arrays, and
 *  specify the sources of the array data for use by
 *  glDrawElements.
 */
   glEnableClientState(GL_VERTEX_ARRAY);
   glEnableClientState(GL_NORMAL_ARRAY);
   glVertexPointer(3, GL_DOUBLE, 0, vertices);
   glNormalPointer(GL_DOUBLE, 0, vnormals);
   glTexCoordPointer(1, GL_DOUBLE, 0, vcrvs);

/*
 *  Initialize for 1-D texture mapping (used for contour
 *  plots).  Set pixel packing alignment, and generate
 *  names.
 */
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
   glGenTextures(3, tex_names);
/*
 *  Texture image1:  Create texture object and assign a
 *  name, change the default texture properties of the
 *  object, and define the texture.
 */
   glBindTexture(GL_TEXTURE_1D, tex_names[0]);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S,
                   GL_REPEAT);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER,
                   GL_NEAREST);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER,
                   GL_NEAREST);
   glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, IW1, 0, GL_RGB,
                GL_UNSIGNED_BYTE, image1);
   glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
             GL_DECAL);
/*
 *  Texture image2.
 */
   glBindTexture(GL_TEXTURE_1D, tex_names[1]);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S,
                   GL_REPEAT);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER,
                   GL_NEAREST);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER,
                   GL_NEAREST);
   glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, IW2, 0, GL_RGB,
                GL_UNSIGNED_BYTE, image2);
/*
 *  Texture image3.
 */
   glBindTexture(GL_TEXTURE_1D, tex_names[2]);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S,
                   GL_REPEAT);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER,
                   GL_NEAREST);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER,
                   GL_NEAREST);
   glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, IW3, 0, GL_RGB,
                GL_UNSIGNED_BYTE, image3);
   glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
             GL_DECAL);
/*
 *  Automatically generate a texture s-coordinate (for each
 *  vertex) based on distance from the x-y plane (z value
 *  in object coordinates).
 */
   glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
   glTexGenfv(GL_S, GL_OBJECT_PLANE, xyplane);

   return;
}


/*
 *----------------------------------------------------------
 *
 *  inputData:  Allocate storage and read nv, vertices, nt,
 *              indices, and (optionally) vnormals from a
 *              user-specified file defining a surface.  If
 *              a second filename is specified on the
 *              command line, parameters np and cpoints
 *              defining a space curve are also read in.
 *
 *  On input:  argc and argv are the parameters input to
 *             main.
 *
 *  On output:  Unless the program is terminated with exit
 *              code 1 (and an error message), global
 *              variables nv, nt, vertices, and indices are
 *              initialized with the input values.  Also, if
 *              normals are included in the data set, then
 *              array vnormals is filled and have_vnormals
 *              is set to True.  Storage is also allocated
 *              for indc, trareas, trnormals, vareas, and
 *              vcrvs.  If a second data set is specified np
 *              and cpoints are also stored.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  None
 *
 *----------------------------------------------------------
 */
static void inputData(int argc, char *argv[])
{
   int i;
   FILE *fp;
   char *filemsg1 = "\n\n"
     "  This program must be called with a path to a file\n"
     "  containing the following records:\n\n"
     "    nv = number of vertices (at least 3)\n"
     "    x y z = coordinates of the first vertex\n"
     "    x y z = coordinates of the second vertex\n"
     "    . . .\n"
     "    x y z = coordinates of the last vertex\n";
   char *filemsg2 =
     "    nt = number of triangles (at least 1)\n"
     "    i1 i2 i3 = CCW-ordered vertex indices of first triangle\n"
     "    i1 i2 i3 = CCW-ordered vertex indices of second triangle\n"
     "    . . .\n"
     "    i1 i2 i3 = CCW-ordered vertex indices of last triangle\n\n"
     "  (Indices are in the range 1 to nv.)\n"
     "  The triangle list may optionally be followed by nv and a\n"
     "  sequence of (x,y,z) coordinates defining vertex normals.\n";
   char *readerror = "  Error reading file.\n";
   char *memerror = "  Unable to allocate sufficient memory.\n";

/*
 *  Open file with pointer fp, and read nv.
 */
   if (argc < 2) {
      printf("%s", filemsg1);
      printf("%s", filemsg2);
      exit(1);
   }
   fp = fopen(argv[1], "r");
   if (fp == NULL) {
     printf("Cannot open file %s.\n", argv[1]);
     exit(1);
   }
   if (fscanf(fp, "%d", &nv) != 1) {    
      printf("%s", readerror);
      exit(1);
   }
   if (nv < 3) {
      printf("Invalid value:  nv = %d.\n", nv);
      exit(1);
   }
/*
 *  Allocate storage for vertices, vnormals, indc, vareas,
 *  and vcrvs.
 */
   vertices = (GLdouble *) calloc (3*nv, sizeof(GLdouble));
   vnormals = (GLdouble *) calloc (3*nv, sizeof(GLdouble));
   indc = (GLuint *) calloc (nv, sizeof(GLuint));
   vareas = (GLdouble *) calloc (nv, sizeof(GLdouble));
   vcrvs = (GLdouble *) calloc (nv, sizeof(GLdouble));
   if (vertices == NULL || vnormals == NULL ||
       indc == NULL || vareas == NULL || vcrvs == NULL) {
      printf("%s", memerror);
      exit(1);
   }
/*
 *  Read vertices.
 */
   for (i = 0; i < nv; i++) {
      if (fscanf(fp, "%lf%lf%lf", vertices+3*i, vertices+3*i+1,
                                  vertices+3*i+2) != 3) {
         printf("%s", readerror);
         exit(1);
      }
   }
/*
 *  Read nt.
 */
   if (fscanf(fp, "%d", &nt) != 1) {    
      printf("%s", readerror);
      exit(1);
   }
   if (nt < 1) {
      printf("Invalid value:  nt = %d.\n", nt);
      exit(1);
   }
/*
 *  Allocate storage for indices, trnormals, and trareas,
 *  and read indices.
 */
   indices = (GLuint *) calloc (3*nt, sizeof(GLuint));
   trnormals = (GLdouble *) calloc (3*nt, sizeof(GLdouble));
   trareas = (GLdouble *) calloc (nt, sizeof(GLdouble));
   if (indices == NULL || trnormals == NULL || trareas == NULL) {
      printf("%s", memerror);
      exit(1);
   }
   for (i = 0; i < nt; i++) {
      if (fscanf(fp, "%u%u%u", indices+3*i, indices+3*i+1,
                               indices+3*i+2) != 3) {   
         printf("%s", readerror);
         exit(1);
      }
      if (indices[3*i  ] < 1 || indices[3*i  ] > (GLuint) nv ||
          indices[3*i+1] < 1 || indices[3*i+1] > (GLuint) nv ||
          indices[3*i+2] < 1 || indices[3*i+2] > (GLuint) nv) {
         printf("The i-th index triple is invalid for "
                "i = %d.\n", i);
         exit(1);
      }
/*
 *  Decrement indices.
 */
      indices[3*i] -= 1;
      indices[3*i+1] -= 1;
      indices[3*i+2] -= 1;
   }
/*
 *  Test for the inclusion of vertex normals.
 */
   if (fscanf(fp, "%d", &i) == 1) {
      if (i != nv) {
         printf("Invalid value:  nv = %d, i = %d.\n", nv, i);
         exit(1);
      }
/*
 *  Read vnormals.
 */
      have_vnormals = GL_TRUE;
      for (i = 0; i < nv; i++) {
         if (fscanf(fp, "%lf%lf%lf", vnormals+3*i, vnormals+3*i+1,
                                     vnormals+3*i+2) != 3) {
            printf("%s", readerror);
            exit(1);
         }
      }
   }
/*
 *  Test for a curve file.
 */
   if (argc <= 2) return;
/*
 *   Open the file, read np, and test for errors.
 */
   fp = fopen(argv[2], "r");
   if (fp == NULL) {
     printf("Cannot open file %s.\n", argv[2]);
     exit(1);
   }
   if (fscanf(fp, "%d", &np) != 1) {
      printf("%s", readerror);
      exit(1);
   }
   if (np < 1) {
      printf("Invalid value:  np = %d.\n", np);
      exit(1);
   }
/*
 *  Allocate storage for cpoints.
 */
   cpoints = (GLdouble *) calloc (3*np, sizeof(GLdouble));
   if (cpoints == NULL) {
      printf("%s", memerror);
      exit(1);
   }
/*
 *  Read cpoints.
 */
   for (i = 0; i < np; i++) {
      if (fscanf(fp, "%lf%lf%lf", cpoints+3*i, cpoints+3*i+1,
                                  cpoints+3*i+2) != 3) {
         printf("%s", readerror);
         exit(1);
      }
   }
   return;
}


/*
 *----------------------------------------------------------
 *
 *  key:  Associate ASCII keys with menu entries.
 *        The mouse coordinates (x,y) are not used.
 *
 *  glsurf function called:  menu
 *  OpenGL functions called:  None
 *
 *----------------------------------------------------------
 */
static void key(unsigned char key, int x, int y)
{
   menu((int) key);
   return;
}


/*
 *----------------------------------------------------------
 *
 *  makeMenu:  Create a menu (for window 1) attached to the
 *             right button.
 *
 *  glsurf function called:  menu
 *  OpenGL functions called:  glutAddMenuEntry,
 *                            glutAddSubMenu,
 *                            glutAttachMenu,
 *                            glutCreateMenu, glutSetMenu
 *
 *----------------------------------------------------------
 */
static void makeMenu(void)
{
   glutCreateMenu(menu);

   glutAddMenuEntry("", 0);
   glutAddMenuEntry(" b:  Toggle display of the bounding Box", 'b');
   glutAddMenuEntry(" C:  Toggle display of Contour plot", 'C');
   glutAddMenuEntry(" c:  Select top (front face) material Color", 'c');
   glutAddMenuEntry(" f:  Toggle wiremesh/polygon-Fill mode", 'f');
   glutAddMenuEntry(" K:  Increase cutaway fraction", 'K');
   glutAddMenuEntry(" k:  Decrease cutaway fraction", 'k');
   glutAddMenuEntry(" m:  Swap top and bottom Material colors", 'm');
   glutAddMenuEntry(" r/Home:  Restore defaults", 'r');
   glutAddMenuEntry(" s:  Toggle flat/smooth Shading", 's');
   glutAddMenuEntry(" t:  Cycle texture maps (contouring)", 't');
   glutAddMenuEntry(" u:  Toggle decal/modulate mode", 'u');
   glutAddMenuEntry(" w:  Copy the Window to a Postscript file", 'w');
   glutAddMenuEntry(" W:  Same as w with black and white reversed", 'W');
   glutAddMenuEntry(" X:  Decrease pitch", 'X');
   glutAddMenuEntry(" x:  Increase pitch", 'x');
   glutAddMenuEntry(" Y:  Decrease yaw", 'Y');
   glutAddMenuEntry(" y:  Increase yaw", 'y');
   glutAddMenuEntry(" Z:  Decrease roll", 'Z');
   glutAddMenuEntry(" z:  Increase roll", 'z');
   glutAddMenuEntry(" PgUp:  Zoom in (middle mouse button)", 0);
   glutAddMenuEntry(" PgDn:  Zoom out (shifted middle mouse button)", 0);
   glutAddMenuEntry(" <esc>:  Terminate the program", 27);

   glutCreateMenu(menu);
   glutAddMenuEntry(" a:  Toggle Aspect ratio preservation mode",'a');
   glutAddMenuEntry(" B:  Compute a Benchmark (timings)", 'B');
   glutAddMenuEntry(" D:  Double Data values (z-values)", 'D');
   glutAddMenuEntry(" d:  Scale down Data values (by half)", 'd');
   glutAddMenuEntry(" E:  Double shininess Exponent", 'E');
   glutAddMenuEntry(" e:  Scale down shininess Exponent (by half)", 'e');
   glutAddMenuEntry(" h:  Toggle hidden line removal for space curve", 'h');
   glutAddMenuEntry(" I:  Rotate light source clockwise about x axis", 'I');
   glutAddMenuEntry(" i:  Rotate light source CCW about the x axis", 'i');
   glutAddMenuEntry(" J:  Rotate light source clockwise about y axis", 'J');
   glutAddMenuEntry(" j:  Rotate light source CCW about the y axis", 'j');
   glutAddMenuEntry(" L:  Toggle illumination", 'L');
   glutAddMenuEntry(" p:  Toggle orthographic/perspective Projection", 'p');

   glutCreateMenu(menu);
   glutAddMenuEntry(" T:  Replace Text (keyboard string entry)", 'T');
   glutAddMenuEntry(" l:  Change text string Location (left button)", 'l');
   glutAddMenuEntry(" 0:  Set 8 by 13 fixed width font", '0');
   glutAddMenuEntry(" 1:  Set 9 by 15 fixed width font", '1');
   glutAddMenuEntry(" 2:  Set 10-point prop. spaced Times Roman", '2');
   glutAddMenuEntry(" 3:  Set 24-point prop. spaced Times Roman", '3');
   glutAddMenuEntry(" 4:  Set 10-point prop. spaced Helvetica", '4');
   glutAddMenuEntry(" 5:  Set 12-point prop. spaced Helvetica", '5');
   glutAddMenuEntry(" 6:  Set 18-point prop. spaced Helvetica", '6');
   glutAddMenuEntry(" 7:  Set prop. spaced Roman stroke font", '7');
   glutAddMenuEntry(" 8:  Set mono-spaced Roman stroke font", '8');

   glutSetMenu(1);
   glutAddSubMenu(" Additional Options -->", 2);
   glutAddSubMenu(" Text -->", 3);

   glutAttachMenu(GLUT_RIGHT_BUTTON);
   return;
}


/*
 *----------------------------------------------------------
 *
 *  menu:  Respond to menu items selected with right button.
 *
 *  glsurf functions called:  benchmark, cckDisplay,
 *                            computeNormals, getLin,
 *                            reshape, scaleData,
 *                            screenDump, storeC, updateR
 *  OpenGL functions called:  glEnable, glDisable, glLightfv,
 *                            glLoadIdentity, glMaterialf,
 *                            glMaterialfv, glPolygonMode,
 *                            glPushMatrix, glPopMatrix,
 *                            glRotated, glShadeModel,
 *                            glTexEnvf, glutGet,
 *                            glutHideWindow, glutMouseFunc,
 *                            glutPostRedisplay,
 *                            glutSetWindow, glutShowWindow
 *
 *----------------------------------------------------------
 */
static void menu(int item)
{
/*
 *  Local variables:
 *
 *  black = Contour line color
 *  i = Loop index
 *  light_position = Initial (unrotated) position of light
 *                   source in eye coordinates
 *  xstep = Rotation angle (in degrees) for for altering
 *          pitch or yaw
 *  zstep = Rotation angle (in degrees) for for altering roll
 */
   GLfloat black[4] = {0.0, 0.0, 0.0, 1.0};
   unsigned int i;  /* Loop index */
   GLfloat light_position[4] = {0.0, 0.0, 1.0, 0.0};
   GLdouble xstep = 0.15;
   GLdouble zstep = 5.00;

   switch (item) {
/*
 *  Menus 1 and 2:
 */
   case 'a':   /* Toggle aspect ratio mode */
      aspect = !aspect;
      reshape(glutGet(GLUT_WINDOW_WIDTH),
              glutGet(GLUT_WINDOW_HEIGHT));
      break;

   case 'B':   /* Compute a benchmark */
      benchmark();
      break;

   case 'b':   /* Toggle bbox */
      bbox = !bbox;
      break;

   case 'C':   /* Toggle contour plotting */
      cplot = !cplot;
      glutSetWindow(3);  /* Make window 3 current */
      if (cplot) {
         glutShowWindow();  /* Make window 3 visible */
      } else {
         if (tex_flag != 1) glutHideWindow();  /* Hide window 3 */
      }
      glutSetWindow(1);  /* Make window 1 current */
      break;

   case 'c':   /* Select (front face) material color */
      glutSetWindow(2);  /* Make window 2 current */
      glutShowWindow();  /* Make window 2 visible */
      break;

   case 'D':   /* Scale data values up */
      scaleData(2.0);
      have_vnormals = GL_FALSE;
      computeNormals();
      storeC();
      reshape(glutGet(GLUT_WINDOW_WIDTH),
              glutGet(GLUT_WINDOW_HEIGHT));
      if (cplot || tex_flag == 1) {
         glutSetWindow(3);  /* Make window 3 current */
         cckDisplay();
         glutSetWindow(1);  /* Make window 1 current */
      }
      break;

   case 'd':   /* Scale data values down */
      scaleData(0.5);
      have_vnormals = GL_FALSE;
      computeNormals();
      storeC();
      reshape(glutGet(GLUT_WINDOW_WIDTH),
              glutGet(GLUT_WINDOW_HEIGHT));
      if (cplot || tex_flag == 1) {
         glutSetWindow(3);  /* Make window 3 current */
         cckDisplay();
         glutSetWindow(1);  /* Make window 1 current */
      }
      break;

   case 'E':   /* Scale shininess exponent up */
      shininess *= 2.0;
      if (shininess > 128.0) shininess = 128.0;
      glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS,
                   shininess);
      break;

   case 'e':   /* Scale shininess exponent down */
      shininess *= 0.5;
      glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS,
                   shininess);
      break;

   case 'f':   /* Toggle polygon fill mode */
      polyfill = !polyfill;
      break;

   case 'h':   /* Toggle hide_curve */
      hide_curve = !hide_curve;
      break;

   case 'I':   /* Rotate the light source clockwise about the x axis */
      glPushMatrix();
      glLoadIdentity();
      updateR ('x', -3.0, &light_angle, light_axis);
      glRotated(light_angle, light_axis[0], light_axis[1],
                light_axis[2]);
      glLightfv(GL_LIGHT0, GL_POSITION, light_position);
      glPopMatrix();
      break;

   case 'i':   /* Rotate the light source CCW about the x axis */
      glPushMatrix();
      glLoadIdentity();
      updateR ('x', 3.0, &light_angle, light_axis);
      glRotated(light_angle, light_axis[0], light_axis[1],
                light_axis[2]);
      glLightfv(GL_LIGHT0, GL_POSITION, light_position);
      glPopMatrix();
      break;

   case 'J':   /* Rotate the light source clockwise about the y axis */
      glPushMatrix();
      glLoadIdentity();
      updateR ('y', -3.0, &light_angle, light_axis);
      glRotated(light_angle, light_axis[0], light_axis[1],
                light_axis[2]);
      glLightfv(GL_LIGHT0, GL_POSITION, light_position);
      glPopMatrix();
      break;

   case 'j':   /* Rotate the light source CCW about the y axis */
      glPushMatrix();
      glLoadIdentity();
      updateR ('y', 3.0, &light_angle, light_axis);
      glRotated(light_angle, light_axis[0], light_axis[1],
                light_axis[2]);
      glLightfv(GL_LIGHT0, GL_POSITION, light_position);
      glPopMatrix();
      break;

   case 'K':   /* Increase cutaway fraction */
      cutaway += 0.05;
      if (cutaway > 1.0) cutaway = 1.0;
      reshape(glutGet(GLUT_WINDOW_WIDTH),
              glutGet(GLUT_WINDOW_HEIGHT));
      break;

   case 'k':   /* Decrease cutaway fraction */
      cutaway -= 0.05;
      if (cutaway < 0.0) cutaway = 0.0;
      reshape(glutGet(GLUT_WINDOW_WIDTH),
              glutGet(GLUT_WINDOW_HEIGHT));
      break;

   case 'L':   /* Toggle illumination */
      lighting = !lighting;
      if (lighting) {
         glEnable(GL_LIGHTING);
      } else {
         glDisable(GL_LIGHTING);
      }
      break;

   case 'm':   /* Swap top and bottom material colors */
      for (i = 0; i < 4; i++) {
         GLfloat t;
         t = top_color[i];
         top_color[i] = bottom_color[i];
         bottom_color[i] = t;
      }
      glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,
                   top_color);
      glMaterialfv(GL_FRONT, GL_SPECULAR, sp_color);
      glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE,
                   bottom_color);
      glMaterialfv(GL_BACK, GL_SPECULAR, sp_color);
      break;

   case 'p':   /* Toggle projection type */
      perspective = !perspective;
      reshape(glutGet(GLUT_WINDOW_WIDTH),
              glutGet(GLUT_WINDOW_HEIGHT));
      break;

   case 'r':   /* Restore defaults */
      aspect = GL_TRUE;
      bbox = GL_TRUE;
      cutaway = 0.0;
      decal = GL_TRUE;
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
                GL_DECAL);
      hide_curve = GL_FALSE;
      light_angle = 0.0;
      light_axis[0] = 1.0;
      light_axis[1] = 0.0;
      light_axis[2] = 0.0;
      lighting = GL_TRUE;
      orient_angle = 0.0;
      orient_axis[0] = 1.0;
      orient_axis[1] = 0.0;
      orient_axis[2] = 0.0;
      perspective = GL_FALSE;
      polyfill = GL_TRUE;
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      shininess = 32.0;
      smooth = GL_TRUE;
      glShadeModel(GL_SMOOTH);
      tex_flag = 0;
      glutSetWindow(3);  /* Make window 3 current */
      glutHideWindow();  /* Hide window 3 */
      glutSetWindow(1);  /* Make window 1 current */
      view_angle = -90.0;
      view_axis[0] = 1.0;
      view_axis[1] = 0.0;
      view_axis[2] = 0.0;
      vsf = 1.0;
      glLoadIdentity();
      glRotated(light_angle, light_axis[0], light_axis[1],
                light_axis[2]);
      glLightfv(GL_LIGHT0, GL_POSITION, light_position);
      reshape(glutGet(GLUT_WINDOW_WIDTH),
              glutGet(GLUT_WINDOW_HEIGHT));
      break;

   case 's':   /* Toggle shading method */
      smooth = !smooth;
      if (smooth) {
         glShadeModel(GL_SMOOTH);
      } else {
         glShadeModel(GL_FLAT);
      }
      break;

   case 't':   /* Switch texture maps */
      tex_flag += 1;
      if (tex_flag > 3) tex_flag = 0;
      if (tex_flag == 1) {
         if (decal) {
            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
                      GL_DECAL);
         } else {
            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
                      GL_MODULATE);
         }
      }
      if (tex_flag == 2) {
         glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
                   GL_BLEND);
         glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR,
                    black);
      }
      if (tex_flag == 3) {
         glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
                   GL_DECAL);
      }
      glutSetWindow(3);  /* Make window 3 current */
      if (tex_flag == 1) {
         glutShowWindow();  /* Make window 3 visible */
      } else {
         if (!cplot) glutHideWindow();  /* Hide window 3 */
      }
      glutSetWindow(1);  /* Make window 1 current */
      break;

   case 'u':   /* Toggle decal/modulate mode */
      decal = !decal;
      if (tex_flag == 1) {
         if (decal) {
            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
                      GL_DECAL);
         } else {
            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
                      GL_MODULATE);
         }
      }
      break;

   case 'W':   /* Copy the window to a Postscript file */
      screenDump(1);   /* Reverse black and white */
      break;
      
   case 'w':   /* Copy the window to a Postscript file */
      screenDump(0);
      break;

   case 'X':   /* Decrease pitch */
      updateR ('x', xstep, &orient_angle, orient_axis);
      break;

   case 'x':   /* Increase pitch */
      updateR ('x', -xstep, &orient_angle, orient_axis);
      break;

   case 'Y':   /* Decrease yaw */
      updateR ('y', -xstep, &orient_angle, orient_axis);
      break;

   case 'y':   /* Increase yaw */
      updateR ('y', xstep, &orient_angle, orient_axis);
      break;

   case 'Z':   /* Decrease roll */
      updateR ('z', zstep, &orient_angle, orient_axis);
      break;

   case 'z':   /* Increase roll */
      updateR ('z', -zstep, &orient_angle, orient_axis);
      break;

   case 27:   /* Terminate program */
      exit(0);
      break;

/*
 *  Menu 3:
 */
   case 'T':   /* Input a new character string for title */
      getLin(title, LTITLE);
      break;

   case 'l':   /* Get a new title location title_pos */
      glutMouseFunc(getRasterPos);
      break;
 
   case '0':  case '1':  case '2':  case '3':  case '4':
   case '5':  case '6':  case '7':  case '8':
      font_number = item - '0';
      break;

   default:
      break;
   }
   glutPostRedisplay();
   return;
}


/*
 *----------------------------------------------------------
 *
 *  motion:  Callback triggered by mouse motion in the
 *           primary window with one or more mouse buttons
 *           pressed.
 *
 *  If the left mouse button is down (lbdown = True),
 *  rotate about the y axis by x-xs if |x-xs| > |y-ys|, or
 *  rotate about the x axis by y-ys if |x-xs| <= |y-ys|, and
 *  update (xs,ys) to the current mouse location.  This
 *  allows the user to drag the surface around one of the
 *  axes at the rate of one degree of rotation per pixel
 *  of mouse motion.
 *
 *  glsurf function called:  updateR
 *  OpenGL function called:  glutPostRedisplay
 * 
 *----------------------------------------------------------
 */
static void motion(int x, int y)
{
   int dx, dy;

   if (lbdown) {
      dx = x - xs;
      dy = y - ys;
      if (abs(dx) <= abs(dy)) {
/*
 *  Rotate about the x-axis by dy.
 */
         updateR ('x', (GLdouble) dy, &view_angle,
                  view_axis);
      } else {
/*
 *  Rotate about the y-axis by dx.
 */
         updateR ('y', (GLdouble) dx, &view_angle,
                  view_axis);
      }
      xs = x;
      ys = y;
      glutPostRedisplay();
   }
}


/*
 *----------------------------------------------------------
 *
 *  mouse:  Callback triggered by a mouse button press or
 *          release with the mouse position in the primary
 *          window (except for buttons attached to menus).
 *
 *  Initialize (xs,ys) and set lbdown to True on left button
 *  press, or clear lbdown to False on left button release.
 * 
 *  Zoom in or zoom out (if a shift key is down) on middle
 *  button press.
 *
 *  glsurf function called:  specialKey
 *  OpenGL functions called:  glutGet, glutGetModifiers,
 *                            glutPostRedisplay
 * 
 *----------------------------------------------------------
 */
static void mouse(int button, int state, int x, int y)
{
   if (button == GLUT_LEFT_BUTTON) {
      if (state == GLUT_DOWN) {
/*
 *  Left button press:  Initialize (xs,ys).
 */
         xs = x;
         ys = y;
         lbdown = GL_TRUE;
      } else {
         lbdown = GL_FALSE;
      }
      return;
   }

   if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
/*
 *  Middle button press:  Zoom in or out.
 */
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
         specialKey(GLUT_KEY_PAGE_DOWN, x, y);
      } else {
         specialKey(GLUT_KEY_PAGE_UP, x, y);
      }
   }
   return;
}


/*
 *----------------------------------------------------------
 *
 *  normalizeVC:  Normalize vertex curvatures vcrv:  values
 *                are clamped to [vm-t*vsd,vm+t*vsd], where
 *                vm is the mean, vsd is the standard devi-
 *                ation, and t = threshold.  The resulting
 *                values are then mapped to [0,1] for
 *                texture mapping.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  None
 *
 *----------------------------------------------------------
 */
static void normalizeVC(void)
{
   int i;                 /* Index for vcrvs */
   GLdouble t = 1.5;      /* Threshold */
   GLdouble vcmin,vcmax;  /* Minimum and maximum curvature
 *                           values before normalization */
   GLdouble vm,vsd;       /* Mean and standard deviation */
/*
 *  Compute vm and vsd.
 */
   vm = 0.0;
   for (i = 0; i < nv; i++) {
      vm += vcrvs[i];
   }
   vm /= (GLdouble) nv;

   vsd = 0.0;
   for (i = 0; i < nv; i++) {
      vsd += (vcrvs[i]-vm)*(vcrvs[i]-vm);
   }
   vsd = sqrt(vsd/(GLdouble) nv);
/*
 *  Clamp vertex curvature values to [vm-t*vsd,vm+t*vsd].
 */
   for (i = 0; i < nv; i++) {
      if (vcrvs[i] < vm-t*vsd) vcrvs[i] = vm-t*vsd;
      if (vcrvs[i] > vm+t*vsd) vcrvs[i] = vm+t*vsd;
   }
/*
 *  Compute the range of clamped values.
 */
   vcmin = 0.0;
   vcmax = 0.0;
   for (i = 0; i < nv; i++) {
      if (vcrvs[i] < vcmin) vcmin = vcrvs[i];
      if (vcrvs[i] > vcmax) vcmax = vcrvs[i];
   }
/*
 *  Normalize the values to [0,1] for texture mapping.
 */
   for (i = 0; i < nv; i++) {
      vcrvs[i] = (vcrvs[i]-vcmin)/(vcmax-vcmin);
   }
   return;
}


/*
 *----------------------------------------------------------
 *
 *  reshape:  Store window dimensions, set the viewport,
 *            construct a projection, and initialize the
 *            modelview matrix.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  glFrustum, glLoadIdentity,
 *                            glMatrixMode, glOrtho,
 *                            glTranslated, glViewport
 *
 *----------------------------------------------------------
 */
static void reshape(int w, int h)
{
   win_width = (GLfloat) w;
   win_height = (GLfloat) h;
   if (aspect) { 
/*
 *  Use a square viewport (assuming square pixels) centered
 *  in the window.
 */
      if (w <= h) {
         glViewport(0, (GLint) (h-w)/2, (GLint) w, (GLint) w);
      } else {
         glViewport((GLint) (w-h)/2, 0, (GLint) h, (GLint) h);
      }
   } else {
/*
 *  Set the viewport to the entire window.
 */
      glViewport(0, 0, (GLint) w, (GLint) h);
   }

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   if (perspective) {
/*
 *  Perspective projection.
 */
      glFrustum(-vsf*r, vsf*r, -vsf*r, vsf*r,
                2.0*r*(1.0+cutaway), 4.0*r);
   } else {
/*
 *  Orthographic projection.
 */
      glOrtho(-vsf*r, vsf*r, -vsf*r, vsf*r,
              2.0*r*(1.0+cutaway), 4.0*r);
   }

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();

   return;
}


/*
 *----------------------------------------------------------
 *
 *  rgb2hsv:  Converts a color specified in the RGB system
 *            to the HSV system.
 *
 *  On input:  r,g,b = Red, green, and blue color compo-
 *                     nents in the range [0,1].  Values are
 *                     clamped if necessary.
 *
 *  On output:  h = Hue in the range 0 to 360 (degrees),
 *                  where red, green, and blue correspond
 *                  to h = 0, h = 120, and h = 240, respec-
 *                  tively.  If s = 0, h is arbitrarily 
 *                  chosen to be 0.
 *
 *              s = Saturation value in [0,1], with s = 0
 *                  corresponding to greyscale (r = g = b),
 *                  and s = 1 corresponding to a pure hue
 *                  (at least one of r, g, and b is 0).
 *
 *              v = Value of intensity in [0,1], with v = 0
 *                  corresponding to black (r = g = b = 0),
 *                  and v = 1 corresponding to maximum
 *                  brightness (r = 1 or g = 1 or b = 1).
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  None
 *
 *----------------------------------------------------------
 */
static void rgb2hsv(float r, float g, float b, float *h,
                    float *s, float *v)
{
   float cmax, cmin;   /* Maximum and minimum components */
   float delta;        /* cmax - cmin */

   if (r >= g) {       /* Compute cmax and cmin */
      if (r >= b) {
         cmax = r;
         cmin = (g <= b) ? g : b;
      } else {
         cmax = b;
         cmin = g;
      }
   } else {            /* r < g */
      if (r >= b) {
         cmax = g;
         cmin = b;
      } else {      
         cmax = (g >= b) ? g : b;
         cmin = r;
      }
   }
   delta = cmax - cmin;
   *v = cmax;
   if (delta == 0.0) {  /* r=g=b  ==>  h=s=0 */
      *s = 0.0;
      *h = 0.0;
   } else {
      *s = delta/cmax;
      if (r == cmax)
         *h = (g-b)/delta;   /* h is between magenta and yellow */
      else if (g == cmax)
         *h = 2.0 + (b-r)/delta;  /* h between yellow and cyan */
      else
         *h = 4.0 + (r-g)/delta;  /* h between cyan and magenta */

      *h *= 60.0;                 /* Convert h to degrees */
      if (*h < 0.0) *h += 360.0;
   }
   return;
}


/*
 *----------------------------------------------------------
 *
 *  scaleData:  Scale z values (3rd components of vertices)
 *              by sf.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  None
 *
 *----------------------------------------------------------
 */
static void scaleData(GLdouble sf)
{
   int i;             /* Index for vertices */

   for (i = 0; i < nv; i++) {
      vertices[3*i+2] *= sf;
   }
   return;
}


/*
 *----------------------------------------------------------
 *
 *  screenDump:  Create a level-2 Encapsulated Postscript
 *               file containing a 24-bit RGB color image of
 *               the current window.
 *
 *  The file name is "screenX.ps", where X is initially 'A'
 *  and is incremented on each call to this function.
 *
 *  The return value is 0 on success, or 1 if an error is
 *  encountered in memory allocation, file creation, or
 *  writing to the file, in which case a message is written
 *  to standard output.
 *
 *  If reversebw is not 0, black and white are reversed in
 *  the image file.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  glPixelStorei, glReadBuffer,
 *                            glReadPixels, glutGet
 *
 *----------------------------------------------------------
 */
static int screenDump(int reversebw)
{
   static char filename[] = "screen@.ps";         /* File name */
   static char hexdigit[] = "0123456789ABCDEF";   /* Conversion table */

   char fbuffer[73];    /* Output file buffer */
   FILE *fp;            /* File pointer */
   float a;             /* Aspect ratio of the window
                             (assuming square pixels) */
   GLint w,h;           /* Width and height in pixels */
   GLubyte *image;      /* Image buffer read from frame buffer */
   int k;               /* Index for fbuffer */
   int line_size;       /* Image line size in bytes */
   int lstr;            /* File image row length in bytes */
   int nch;             /* Number of characters written,
                             or error flag  */
   int pw,ph;           /* Width and height in points */
   int px1,py1;         /* Lower left corner coordinates
                             in points */
   int px2,py2;         /* Upper right corner coordinates
                             in points */
   long i;              /* Index for image buffer */
   long image_size;     /* Image size in bytes */
/*
 *  Get window dimensions in pixels, and compute line_size
 *  and image_size.
 */
   w = glutGet(GLUT_WINDOW_WIDTH);
   h = glutGet(GLUT_WINDOW_HEIGHT);
   if (!glutGet(GLUT_WINDOW_RGBA)) {
      printf("screenDump:  Error:  the window must be in "
             "RGB mode\n");
      return 1;
   }
   line_size = w*3;            /* Three bytes per pixel */
   image_size = (long) line_size*h;
/*
 *  Allocate storage for the image.
 */
   image = (GLubyte *) malloc(image_size);
   if (image == NULL) {
      printf("screenDump:  Error allocating image "
             "memory of size %ld\n", image_size);
      return 1;
   }
/*
 *  Open an ASCII file.
 */
   filename[6]++;
   fp = fopen(filename, "w");
   if (fp == NULL) {
      printf("screenDump:  Error opening file\n");
      free(image);
      return 1;
   }
/*
 *  Compute corner coordinates of the Postscript display
 *  area and bounding box.  The coordinates, specified in
 *  user space units (points, at 72 points/inch), are
 *  chosen to leave at least (1/2)-inch margins on an 8.5
 *  by 11 inch page (with origin at the lower left), to
 *  preserve the aspect ratio a, and to center the image.
 */
   a = (float) w / (float) h;
   pw = 540;
   ph = (int) (((float) pw / a) + 0.5);
   if (ph <= 720) {
      px1 = 36;
      px2 = 576;
      py1 = 396 - ph/2;
      py2 = py1 + ph;
   } else {
      ph = 720;
      py1 = 36;
      py2 = 756;
      pw = (int) (720.0*a + 0.5);
      px1 = 306 - pw/2;
      px2 = px1 + pw;
   }
/*
 *  Write header comments and string definition.  A string
 *  is associated with a scanline (image row) and thus has
 *  length equal to the number of bytes per scanline (with
 *  24-bit pixels).
 */
   lstr = 3*w;
   nch = fprintf(fp, "%%!PS-Adobe-2.0 EPSF-1.2\n"
                     "%%%%BoundingBox: %d %d %d %d\n"
                     "%%%%Title:  Screen Dump\n"
                     "%%%%Creator:  glsurf\n"
                     "%%%%EndComments\n"
                     "/picstr %d string def\n",
                        px1, py1, px2, py2, lstr);
   if (nch < 0) {
      printf("screenDump:  Write error");
      free(image);
      fclose(fp);
      return 1;
   }
/*
 *  Write commands to map the unit square in user space to
 *  the display area, and output the colorimage command:
 *  width, height, bits/component, image matrix, read-string
 *  procedure, flag specifying a single source, and number
 *  of components.  The image matrix (six numbers associated
 *  with homogeneous coordinates) maps the unit square in
 *  user space to the source image.
 */
   nch = fprintf(fp, "%d %d translate\n"
                     "%d %d scale\n"
                     "%d %d %d\n"
                     "[%d %d %d %d %d %d]\n"
                     "{currentfile picstr readhexstring pop} "
                     "false 3 colorimage\n",
                       px1, py1, pw, ph,
                       w, h, 8,
                       w, 0, 0, h, 0, 0);
   if (nch < 0) {
      printf("screenDump:  Write error");
      free(image);
      fclose(fp);
      return 1;
   }
/*
 *  Read the three color planes from the front frame buffer
 *  into the image buffer.
 */
   glPixelStorei(GL_PACK_ALIGNMENT, 1);
   glReadBuffer(GL_FRONT);
   glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, image);
/*
 *  Write out the image data as a string of hex digits, two
 *  per pixel, left-to-right within bottom-to-top.
 */
   k = 0;
   fbuffer[72] = '\0';
   for (i = 0; i < image_size; i += 3) {
      unsigned int r, g, b;   /* Color components of a pixel */
      r = (unsigned int) *(image+i);
      g = (unsigned int) *(image+i+1);
      b = (unsigned int) *(image+i+2);
      if (reversebw) {
         if (r == 0  &&  g == 0  &&  b == 0) {
            r = 0xff;  g = 0xff;  b = 0xff;   /* Black to white */
         } else
         if (r == 0xff  &&  g == 0xff  &&  b == 0xff) {
            r = 0;  g = 0;  b = 0;            /* White to black */
         }
      }
      fbuffer[k] = hexdigit[r >> 4];      /* High hex digit of r */
      fbuffer[k+1] = hexdigit[r & 0x0f];  /* Low hex digit of r */
      fbuffer[k+2] = hexdigit[g >> 4];
      fbuffer[k+3] = hexdigit[g & 0x0f];
      fbuffer[k+4] = hexdigit[b >> 4];
      fbuffer[k+5] = hexdigit[b & 0x0f];
      k = k + 6;
      if (k == 72) {
         fprintf(fp, "%s\n", fbuffer);
         k = 0;
      }
   }
   if (k > 0) {
      fbuffer[k] = '\0';
      fprintf(fp, "%s\n", fbuffer);
   }
/*
 *  Write the showpage command and end-of-file indicator.
 */
   nch = fprintf(fp, "showpage\n"
                     "%%%%EOF\n");
/*
 *  HP's interpreters require a one-byte End-of-Postscript-
 *  Job indicator (to eliminate a timeout error message):
 *  ASCII 4.
 */
   nch = fprintf(fp, "\x04");
   if (nch < 0) {
      printf("screenDump:  Write error");
      free(image);
      fclose(fp);
      return 1;
   }
   free(image);
   fclose(fp);
   return 0;
}


/*
 *----------------------------------------------------------
 *
 *  specialKey:  Rotate viewer location about x or y axis
 *               in response to arrow keys, zoom in or out
 *               on Page Up or Page Down keys, or restore
 *               defaults if Home key is pressed.
 *
 *  glsurf functions called:  menu, reshape, updateR
 *  OpenGL functions called:  glutGetModifiers,
 *                            glutPostRedisplay
 *
 *----------------------------------------------------------
 */
static void specialKey(int key, int x, int y)
{
   GLdouble step;     /* Rotation angle in degrees */

   if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
      step = 3.0;
   } else {
      step = 15.0;
   }
   switch (key) {
      case GLUT_KEY_UP:
         updateR ('x', step, &view_angle, view_axis);
         break;
      case GLUT_KEY_DOWN:
         updateR ('x', -step, &view_angle, view_axis);
         break;
      case GLUT_KEY_LEFT:
         updateR ('y', step, &view_angle, view_axis);
         break;
      case GLUT_KEY_RIGHT:
         updateR ('y', -step, &view_angle, view_axis);
         break;

      case GLUT_KEY_PAGE_UP:     /* Zoom in */
         vsf *= 0.8;
         reshape(win_width, win_height);
         break;

      case GLUT_KEY_PAGE_DOWN:   /* Zoom out */
         vsf *= 1.25;
         reshape(win_width, win_height);
         break;

      case GLUT_KEY_HOME:        /* Restore defaults */
         menu((int) 'r');
         return;

      default:
         break;
   }
   glutPostRedisplay();
 
   return;
}


/*
 *----------------------------------------------------------
 *
 *  storeC:  Compute the corner coordinates (xmin,ymin,zmin)
 *           and (xmax,ymax,zmax), center c = (xc,yc,zc),
 *           and radius r of the bounding box, the contour
 *           intervals contours, the vertex contour/color
 *           indices indc, and the height zcon.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  None
 *
 *----------------------------------------------------------
 */
static void storeC(void)
{
   int i;                /* Index for vertices and indc*/
   int ic;               /* Index for contours */
   GLdouble x, y, z;     /* Components of a vertex */
   GLdouble dx, dy, dz;  /* Component interval widths */

   xmin = vertices[0];
   xmax = xmin;
   ymin = vertices[1];
   ymax = ymin;
   zmin = vertices[2];
   zmax = zmin;

   for (i = 1; i < nv; i++) {
      x = vertices[3*i];
      y = vertices[3*i+1];
      z = vertices[3*i+2];
      if (x < xmin) xmin = x;
      if (x > xmax) xmax = x;
      if (y < ymin) ymin = y;
      if (y > ymax) ymax = y;
      if (z < zmin) zmin = z;
      if (z > zmax) zmax = z;
   }

   xc = (xmin+xmax)/2.0;
   yc = (ymin+ymax)/2.0;
   zc = (zmin+zmax)/2.0;
 
   dx = xmax - xmin;
   dy = ymax - ymin;
   dz = zmax - zmin;
   r = sqrt(dx*dx + dy*dy + dz*dz)/2.0;

/*
 *  Compute the data subinterval width dz associated with
 *  each contour interval, and store an increasing sequence
 *  of contour values uniformly distributed in [zmin,zmax].
 */
   dz /= (GLdouble) nc;
   contours[0] = zmin + dz;
   for (ic = 1; ic < nc; ic++) {
      contours[ic] = contours[ic-1] + dz;
   }

/*
 *  Compute and store vertex contour/color indices indc.
 */
   for (i = 0; i < nv; i++) {
      ic = (GLuint) ((vertices[3*i+2] - zmin)/dz);
      if (ic >= nc) ic = nc - 1;
      indc[i] = ic;
   }

   zcon = zmax;
   return;
}


/*
 *----------------------------------------------------------
 *
 *  updateR:  Update an angle a and axis of rotation (unit
 *            vector) u to reflect an incremental rotation
 *            through angle d about one of the coordinate
 *            axes; i.e., if R denotes the rotation operator
 *            associated with angle a and axis u, and Q
 *            denotes the rotation through angle d about the
 *            coordinate axis, then R is effectively
 *            replaced by Q*R.
 *
 *  On input:  c = 'x', 'y', or 'z', specifying the coordi-
 *             nate axis, d is the incremental rotation
 *             angle in degrees, a is the address of the
 *             angle (in degrees) to be updated, and u is an
 *             array of length 3 containing the unit direc-
 *             tion vector to be updated.
 *
 *  On output:  Unless c is not one of the valid characters,
 *              a and u are updated with the new angle and
 *              rotation axis, respectively.
 *
 *  glsurf functions called:  None
 *  OpenGL functions called:  None
 *
 *----------------------------------------------------------
 */
static void updateR (unsigned char c, GLdouble d,
                     GLdouble *a, GLdouble *u)
{
   GLdouble ad2;   /* a/2 in radians */
   GLdouble ca;    /* cos(a/2) */
   GLdouble cd;    /* cos(d/2) */
   GLdouble cfac;  /* Degrees-to-radians conversion factor */
   GLdouble un;    /* Normalization factor for u */
   GLdouble dd2;   /* d/2 in radians */
   GLdouble q[4];  /* Quaternion equivalent of (a,u) */
   GLdouble r[4];  /* Quaternion equivalent of (a,u) for
                        the updated values */
   GLdouble sa;    /* sin(a/2) */
   GLdouble sd;    /* sin(d/2) */
/*
 *  Convert (a,u) to the equivalent quaternion
 *  q = (cos(a/2),sin(a/2)*u)
 */
   cfac = atan(1.0)/45.0;
   ad2 = cfac*(*a)/2.0;
   ca = cos(ad2);
   sa = sin(ad2);
   q[0] = ca;
   q[1] = sa*u[0];
   q[2] = sa*u[1];
   q[3] = sa*u[2];
/*
 *  Convert (d,v) to the equivalent quaterion
 *  p = (cos(d/2),sin(d/2)*v) for v = e1, e2, or e3.
 */
   dd2 = cfac*d/2.0;
   cd = cos(dd2);
   sd = sin(dd2);
/*
 *  Compute the product r = p*q.
 */
   switch (c) {
      case 'x':
         r[0] = cd*q[0] - sd*q[1];
         r[1] = cd*q[1] + sd*q[0];
         r[2] = cd*q[2] - sd*q[3];
         r[3] = cd*q[3] + sd*q[2];
         break;

      case 'y':
         r[0] = cd*q[0] - sd*q[2];
         r[1] = cd*q[1] + sd*q[3];
         r[2] = cd*q[2] + sd*q[0];
         r[3] = cd*q[3] - sd*q[1];
         break;

      case 'z':
         r[0] = cd*q[0] - sd*q[3];
         r[1] = cd*q[1] - sd*q[2];
         r[2] = cd*q[2] + sd*q[1];
         r[3] = cd*q[3] + sd*q[0];
         break;

      default:
         return;
   }
/*
 *  Convert r to (a,u), ensuring that |u| = 1.
 */
   *a = 2.0*acos(r[0])/cfac;
   un = sqrt(r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);
   if (un > 0) {
      u[0] = r[1]/un;
      u[1] = r[2]/un;
      u[2] = r[3]/un;
   }
   return;
}


/*
 *----------------------------------------------------------
 *
 *  main:  Input data, compute normals, create windows,
 *         create a menu, register callback functions, and
 *         enter a loop.
 *
 *  glsurf functions called:  computeNormals, csInit, init1, 
 *                            inputData, makeMenu, storeC
 *  OpenGL functions called:  glClearColor, glShadeModel,
 *                            glutCreateWindow,
 *                            glutDisplayFunc,
 *                            glutHideWindow, glutInit,
 *                            glutInitDisplayMode,
 *                            glutInitWindowPosition,
 *                            glutInitWindowSize,
 *                            glutKeyboardFunc,
 *                            glutMainLoop, glutMotionFunc,
 *                            glutMouseFunc,
 *                            glutReshapeFunc, 
 *                            glutSetWindow, glutSpecialFunc
 *
 *----------------------------------------------------------
 */
int main(int argc, char *argv[])
{
   inputData(argc, argv);   /* Input data */
   computeNormals();        /* Compute surface normals */
   storeC();                /* Compute center, radius,
                                 contours, indc, and zcon */
   dmin = zmin;             /* Save initial range of data */
   dmax = zmax;             /*   values (for labels).     */

   glutInit(&argc, argv);   /* Initialize Glut library */
/*
 *  Create the primary window (surface display):  window 1.
 */
   glutInitWindowPosition(100, 50);
   glutInitWindowSize(500, 500);    /* Size in pixels */
   glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
   glutCreateWindow(argv[1]);       /* Create window 1 */

   makeMenu();                 /* Create a menu attached
                                    to the right button */

   glutKeyboardFunc(key);      /* Register callbacks */
   glutSpecialFunc(specialKey);
   glutMotionFunc(motion);
   glutMouseFunc(mouse);
   glutReshapeFunc(reshape); 
   glutDisplayFunc(display);

   init1();                    /* One-time initialization */
/*
 *  Create the second window (color selection):  window 2.
 */
   glutInitWindowPosition(650, 250);
   glutInitWindowSize(300, 400);         /* Size in pixels */
   glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
   glutCreateWindow("Color Selection");  /* Create window 2 */
   glClearColor (0.5, 0.5, 0.5, 1.0);    /* Grey background */
   glShadeModel(GL_FLAT);                /* Default is smooth */

   glutKeyboardFunc(csKey);      /* Register callbacks */
   glutMouseFunc(csMouse);
   glutMotionFunc(csMotion);
   glutReshapeFunc(csReshape);
   glutDisplayFunc(csDisplay);
   csInit();                     /* Create hstable */
   glutHideWindow();             /* Initially hidden */
/*
 *  Create the third window (contour color key):  window 3.
 */
   glutInitWindowPosition(605, 50);
   glutInitWindowSize(100, 500);           /* Size in pixels */
   glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
   glutCreateWindow("Contour Color Key");  /* Create window 3 */
   glShadeModel(GL_FLAT);                /* Default is smooth */
   glutKeyboardFunc(cckKey);             /* Register callbacks */
   glutReshapeFunc(cckReshape);
   glutDisplayFunc(cckDisplay);
   glutHideWindow();             /* Initially hidden */

   glutSetWindow(1);             /* Make window 1 current */
   glutMainLoop();
   return 0;
}
