/****************************************************************

Generate a random instance for solving for the conics that meets
4 lines and 2 points (1 point being [-1,0,0]).  

This creates the parameter homotopy file used by Bertini as 
well as the polynomial system used by alphaCertified.

September 26, 2010
Jonathan D. Hauenstein   

****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void createInputFiles(char *inputName, char *certifyName);
void createInputFile(char *name, mpq_t **A, mpq_t **B);
void createCertifyFile(char *name, mpq_t **A, mpq_t **B);
void createRandomRat(char **str);

int main()
{ // seed the random number generator
  srand((unsigned int) time(NULL));

  // construct input file
  createInputFiles("Conic2.bertini", "Conic2.poly");

  return 0;
}

void createInputFiles(char *inputName, char *certifyName)
{ // create the input files
  int i, j;
  char *str = NULL;

  // allocate and initialize memory
  mpq_t **A = (mpq_t **)malloc(5 * sizeof(mpq_t *)), **B = (mpq_t **)malloc(5 * sizeof(mpq_t *));
  for (i = 0; i < 5; i++)
  {
    A[i] = (mpq_t *)malloc(3 * sizeof(mpq_t));
    B[i] = (mpq_t *)malloc(3 * sizeof(mpq_t));
    for (j = 0; j < 3; j++)
    { // initialize
      mpq_init(A[i][j]);
      mpq_init(B[i][j]);

      // create random numbers
      createRandomRat(&str);
      mpq_set_str(A[i][j], str, 10);
      mpq_canonicalize(A[i][j]);
      createRandomRat(&str);
      mpq_set_str(B[i][j], str, 10);
      mpq_canonicalize(B[i][j]);
    }
  }

  // create the input file for Bertini and the input file for alphaCertified
  createInputFile(inputName, A, B);
  createCertifyFile(certifyName, A, B);

  // clear memory
  free(str);
  for (i = 0; i < 5; i++)
  {
    for (j = 0; j < 3; j++)
    {
      mpq_clear(A[i][j]);
      mpq_clear(B[i][j]);
    }
    free(A[i]);
    free(B[i]);
  }
  free(A);
  free(B);

  return;
}

void createInputFile(char *name, mpq_t **A, mpq_t **B)
{
  int j, k;
  mpf_t temp;
  FILE *OUT = fopen(name, "w");

  // print random numbers using 300 digits
  mpf_init2(temp, 1024);

  // initial configuration data
  fprintf(OUT, "CONFIG\nUSERHOMOTOPY: 1;\nPRINTPATHMODULUS: 100;\nMAXNORM: 1e100;\nSECURITYLEVEL: 1;\nTRACKTOLBEFOREEG: 1e-7;\nTRACKTOLDURINGEG: 1e-7;\nFINALTOL: 1e-9;\nSHARPENDIGITS: 75;\nCONDNUMTHRESHOLD: 1e100;\nMAXSTEPSIZE: 0.02;\nCOEFFBOUND: 25;\nDEGREEBOUND: 3;\nSAMPLEFACTOR: 0.25;\nEND;\n");
  // input declarations
  fprintf(OUT, "INPUT\nvariable t1,t2,t3,t4,a,b,c,e,f,g;\nfunction f1,f2,f3,f4,f5,f6,f7,f8,f9,f10;\nconstant p111,v111,p112,v112,p113,v113,p121,v121,p122,v122,p123,v123,p131,v131,p132,v132,p133,v133,p141,v141,p142,v142,p143,v143;\nconstant x15,y15,z15,x05,y05,z05;\nconstant p011,v011,p012,v012,p013,v013,p021,v021,p022,v022,p023,v023,p031,v031,p032,v032,p033,v033,p041,v041,p042,v042,p043,v043;\npathvariable t;\nparameter s;\ns = t;\none_minus_s = 1 - s;\n");
  // starting contant values
  fprintf(OUT, "p111 = -.43430128096558950900885065493639558553695678710938 + I*-.90076767112926969804931331964326091110706329345703;\nv111 = .84169583859478780407670228669303469359874725341797 + I*-.53995195646669991162980295484885573387145996093750;\np112 = .85029736124152752174687464048474794253706932067871e-1 + I*-.99637841404491345187466322386171668767929077148438;\nv112 = .39522265332853717678673888258344959467649459838867 + I*-.91858535493221893375448416918516159057617187500000;\np113 = .63143726363749130836566791913355700671672821044922 + I*.77542696760558782465011518070241436362266540527344;\nv113 = .98071111717300918364514927816344425082206726074219 + I*-.19546279608475009004209255181194748729467391967773;\np121 = -.83797690901853982836655632127076387405395507812500 + I*-.54570568986563960933722228219266980886459350585938;\nv121 = .76532301540561686881147807071101851761341094970703 + I*-.64364639522835342955175974566373042762279510498047;\np122 = .99314670367952484486551156805944629013538360595703 + I*.11687439826794455977054809636683785356581211090088;\nv122 = .76185598631517148504599390435032546520233154296875 + I*-.64774644430960592877966064406791701912879943847656;\np123 = -.14420489238870623505128776287165237590670585632324 + I*-.98954785079406926140421774107380770146846771240234;\nv123 = .68196836727835075109283025085460394620895385742188 + I*-.73138166919311053959518176270648837089538574218750;\np131 = .65481401335674405128628450256655924022197723388672 + I*.75579005544637445357381011490360833704471588134766;\nv131 = .66366404128464795419972688250709325075149536132812 + I*-.74803077497234626846989158366341143846511840820312;\np132 = -.38515023867666903223394569977244827896356582641602 + I*-.92285388531842083015277466984116472303867340087891;\nv132 = .97485914781713112198247017659014090895652770996094 + I*-.22282199603552815814211385259113740175962448120117;\np133 = .79050723801252975597719796496676281094551086425781 + I*-.61245269747940656035467554829665459692478179931641;\nv133 = -.55452096927409122439911470792139880359172821044922 + I*-.83216975109368296337208903423743322491645812988281;\np141 = -.70191820130700333102424792741658166050910949707031 + I*.71225756484149838065889071003766730427742004394531;\nv141 = .84778390295906769225098287279251962900161743164062 + I*.53034182739388957550374925631331279873847961425781;\np142 = -.70911523877004656313260966271627694368362426757812 + I*.70509260253111427640959618656779639422893524169922;\nv142 = .63357068188932053054429616167908534407615661621094 + I*-.77368481376481812450407460346468724310398101806641;\np143 = .13372318114258330279398023776593618094921112060547 + I*-.99101872375102373347033335448941215872764587402344;\nv143 = -.85242882424413080055813907165429554879665374755859 + I*-.52284328397500590135393849777756258845329284667969;\n");
  fprintf(OUT, "x15 = -.76602514768138119105600480907014571130275726318359 + I*-.64281060439270765183294997768825851380825042724609;\ny15 = -.95293578205969198258173946669558063149452209472656 + I*-.30317222047919084593203820077178534120321273803711;\nz15 = -.74542607526671844642152109372545965015888214111328 + I*.66658830346208197692448038651491515338420867919922;\n");

  // new constant values
  for (j = 1; j <= 4; j++)
    for (k = 1; k <= 3; k++)
    { // print values
      fprintf(OUT, "p0%d%d = ", j, k);
      mpf_set_q(temp, A[j-1][k-1]);
      mpf_out_str(OUT, 10, 0, temp);
      fprintf(OUT, ";\nv0%d%d = ", j, k);
      mpf_set_q(temp, B[j-1][k-1]);
      mpf_out_str(OUT, 10, 0, temp);
      fprintf(OUT, ";\n");
    }
  for (k = 0; k < 3; k++)
  {
    fprintf(OUT, "%s05 = ", k == 0 ? "x" : (k == 1 ? "y" : "z"));
    mpf_set_q(temp, A[4][k]);
    mpf_out_str(OUT, 10, 0, temp);
    fprintf(OUT, ";\n");
  }

  // print the homotopy
  fprintf(OUT, "x1 = (p111*s + one_minus_s*p011) + t1*(v111*s + one_minus_s*v011);\ny1 = (p112*s + one_minus_s*p012) + t1*(v112*s + one_minus_s*v012);\nx2 = (p121*s + one_minus_s*p021) + t2*(v121*s + one_minus_s*v021);\ny2 = (p122*s + one_minus_s*p022) + t2*(v122*s + one_minus_s*v022);\nx3 = (p131*s + one_minus_s*p031) + t3*(v131*s + one_minus_s*v031);\ny3 = (p132*s + one_minus_s*p032) + t3*(v132*s + one_minus_s*v032);\nx4 = (p141*s + one_minus_s*p041) + t4*(v141*s + one_minus_s*v041);\ny4 = (p142*s + one_minus_s*p042) + t4*(v142*s + one_minus_s*v042);\nx5 = x15*s + one_minus_s*x05;\ny5 = y15*s + one_minus_s*y05;\nf1 = (z15*s + z05*one_minus_s) - (f*(x5 + 1) + g*y5);\nf2 = a*x5^2 + b*x5*y5 + c*y5^2 + (a+1)*x5 + e*y5 + 1;\nf3 = (p113*s + one_minus_s*p013 + t1*(v113*s + one_minus_s*v013)) - (f*(x1 + 1) + g*y1);\nf4 = (p123*s + one_minus_s*p023 + t2*(v123*s + one_minus_s*v023)) - (f*(x2 + 1) + g*y2);\nf5 = (p133*s + one_minus_s*p033 + t3*(v133*s + one_minus_s*v033)) - (f*(x3 + 1) + g*y3);\nf6 = (p143*s + one_minus_s*p043 + t4*(v143*s + one_minus_s*v043)) - (f*(x4 + 1) + g*y4);\nf7 = a*x1^2 + b*x1*y1 + c*y1^2 + (a+1)*x1 + e*y1 + 1;\nf8 = a*x2^2 + b*x2*y2 + c*y2^2 + (a+1)*x2 + e*y2 + 1;\nf9 = a*x3^2 + b*x3*y3 + c*y3^2 + (a+1)*x3 + e*y3 + 1;\nf10= a*x4^2 + b*x4*y4 + c*y4^2 + (a+1)*x4 + e*y4 + 1;\nEND;\n");

  // close the file
  fclose(OUT);

  mpf_clear(temp);

  return;
}

void createCertifyFile(char *name, mpq_t **A, mpq_t **B)
{
  int j, k;
  mpq_t temp1, temp2, one;
  FILE *OUT = fopen(name, "w");
  mpq_init(temp1);  mpq_init(temp2);  mpq_init(one);
  mpq_set_ui(one, 1, 1);

  fprintf(OUT, "10 10\n\n");

  fprintf(OUT, "3\n");
  fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 "); mpq_out_str(OUT, 10, A[4][2]); fprintf(OUT, " 0\n");
  fprintf(OUT, "0 0 0 0 0 0 0 0 1 0 "); mpq_add(temp1, A[4][0], one); mpq_neg(temp1, temp1); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
  fprintf(OUT, "0 0 0 0 0 0 0 0 0 1 "); mpq_neg(temp1, A[4][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");

  fprintf(OUT, "5\n");
  fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 "); mpq_add(temp1, A[4][0], one); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
  fprintf(OUT, "0 0 0 0 1 0 0 0 0 0 "); mpq_add(temp1, A[4][0], one); mpq_mul(temp1, temp1, A[4][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
  fprintf(OUT, "0 0 0 0 0 1 0 0 0 0 "); mpq_mul(temp1, A[4][0], A[4][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
  fprintf(OUT, "0 0 0 0 0 0 1 0 0 0 "); mpq_mul(temp1, A[4][1], A[4][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
  fprintf(OUT, "0 0 0 0 0 0 0 1 0 0 "); mpq_out_str(OUT, 10, A[4][1]); fprintf(OUT, " 0\n");

  for (j = 0; j < 4; j++)
  {
    fprintf(OUT, "6\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 "); mpq_out_str(OUT, 10, A[j][2]); fprintf(OUT, " 0\n");
    for (k = 0; k < 4; k++) fprintf(OUT, "%d ", k == j);
    fprintf(OUT, "0 0 0 0 0 0 "); mpq_out_str(OUT, 10, B[j][2]); fprintf(OUT, " 0\n");
    for (k = 0; k < 4; k++) fprintf(OUT, "%d ", k == j);
    fprintf(OUT, "0 0 0 0 1 0 "); mpq_neg(temp1, B[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 4; k++) fprintf(OUT, "%d ", k == j);
    fprintf(OUT, "0 0 0 0 0 1 "); mpq_neg(temp1, B[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 1 0 "); mpq_add(temp1, A[j][0], one); mpq_neg(temp1, temp1); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 1 "); mpq_neg(temp1, A[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");

    fprintf(OUT, "13\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 "); mpq_add(temp1, A[j][0], one); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 4; k++) fprintf(OUT, "%d ", k == j);
    fprintf(OUT, "1 0 0 0 0 0 "); mpq_mul(temp1, A[j][0], B[j][0]); mpq_add(temp1, temp1, temp1); mpq_add(temp1, temp1, B[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 4; k++) fprintf(OUT, "%d ", k == j);
    fprintf(OUT, "0 1 0 0 0 0 "); mpq_mul(temp1, A[j][0], B[j][1]); mpq_mul(temp2, B[j][0], A[j][1]); mpq_add(temp1, temp1, temp2); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 4; k++) fprintf(OUT, "%d ", 2*(k == j));
    fprintf(OUT, "0 1 0 0 0 0 "); mpq_mul(temp1, B[j][0], B[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 4; k++) fprintf(OUT, "%d ", k == j);
    fprintf(OUT, "0 0 1 0 0 0 "); mpq_mul(temp1, A[j][1], B[j][1]); mpq_add(temp1, temp1, temp1); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 4; k++) fprintf(OUT, "%d ", 2*(k == j));
    fprintf(OUT, "1 0 0 0 0 0 "); mpq_mul(temp1, B[j][0], B[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 4; k++) fprintf(OUT, "%d ", 2*(k == j));
    fprintf(OUT, "0 0 1 0 0 0 "); mpq_mul(temp1, B[j][1], B[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 4; k++) fprintf(OUT, "%d ", k == j);
    fprintf(OUT, "0 0 0 1 0 0 "); mpq_out_str(OUT, 10, B[j][1]); fprintf(OUT, " 0\n");
    for (k = 0; k < 4; k++) fprintf(OUT, "%d ", k == j);
    fprintf(OUT, "0 0 0 0 0 0 "); mpq_out_str(OUT, 10, B[j][0]); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 1 0 0 0 0 0 "); mpq_mul(temp1, A[j][0], A[j][0]); mpq_add(temp1, temp1, A[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 1 0 0 0 0 "); mpq_mul(temp1, A[j][0], A[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 1 0 0 0 "); mpq_mul(temp1, A[j][1], A[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 1 0 0 "); mpq_out_str(OUT, 10, A[j][1]); fprintf(OUT, " 0\n");
  }

  fclose(OUT);
  mpq_clear(temp1);
  mpq_clear(temp2);
  mpq_clear(one);

  return;
}

void createRandomRat(char **str)
{ // create random rational number
  int i, curr = 0;

  // allocate
  *str = realloc(*str, 25 * sizeof(char));

  // + or -
  if (rand() % 2)
  {
    (*str)[curr] = '-';
    curr++;
  }

  for (i = 0; i < 10; i++)
  {
    (*str)[curr] = 48 + (rand() % 10);
    curr++;
  }
  (*str)[curr] = '/';
  curr++;
  for (i = 0; i < 10; i++)
  {
    (*str)[curr] = 48 + (rand() % 10);
    curr++;
  }
  (*str)[curr] = '\0';

  return;
}




