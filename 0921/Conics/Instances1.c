/****************************************************************

Generate a random instance for solving for the conics that meets
6 lines and the point [-1,0,0].

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
  createInputFiles("Conic1.bertini", "Conic1.poly");

  return 0;
}

void createInputFiles(char *inputName, char *certifyName)
{ // create the input files
  int i, j;
  char *str = NULL;

  // allocate and initialize memory
  mpq_t **A = (mpq_t **)malloc(6 * sizeof(mpq_t *)), **B = (mpq_t **)malloc(6 * sizeof(mpq_t *));
  for (i = 0; i < 6; i++)
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
  for (i = 0; i < 6; i++)
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
  fprintf(OUT, "INPUT\nvariable t1,t2,t3,t4,t5,t6,a,b,c,e,f,g;\nfunction f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12;\nconstant p111,v111,p112,v112,p113,v113,p121,v121,p122,v122,p123,v123,p131,v131,p132,v132,p133,v133,p141,v141,p142,v142,p143,v143,p151,v151,p152,v152,p153,v153,p161,v161,p162,v162,p163,v163;\n");
  fprintf(OUT, "constant p011,v011,p012,v012,p013,v013,p021,v021,p022,v022,p023,v023,p031,v031,p032,v032,p033,v033,p041,v041,p042,v042,p043,v043,p051,v051,p052,v052,p053,v053,p061,v061,p062,v062,p063,v063;\npathvariable t;\nparameter s;\ns = t;\none_minus_s = 1-s;\n");
  // old contant values
  fprintf(OUT, "p111 = -.43430128096558950900885065493639558553695678710938 + I*-.90076767112926969804931331964326091110706329345703;\nv111 = .84169583859478780407670228669303469359874725341797 + I*-.53995195646669991162980295484885573387145996093750;\np112 = .85029736124152752174687464048474794253706932067871e-1 + I*-.99637841404491345187466322386171668767929077148438;\nv112 = .39522265332853717678673888258344959467649459838867 + I*-.91858535493221893375448416918516159057617187500000;\np113 = .63143726363749130836566791913355700671672821044922 + I*.77542696760558782465011518070241436362266540527344;\nv113 = .98071111717300918364514927816344425082206726074219 + I*-.19546279608475009004209255181194748729467391967773;\np121 = -.83797690901853982836655632127076387405395507812500 + I*-.54570568986563960933722228219266980886459350585938;\nv121 = .76532301540561686881147807071101851761341094970703 + I*-.64364639522835342955175974566373042762279510498047;\np122 = .99314670367952484486551156805944629013538360595703 + I*.11687439826794455977054809636683785356581211090088;\nv122 = .76185598631517148504599390435032546520233154296875 + I*-.64774644430960592877966064406791701912879943847656;\np123 = -.14420489238870623505128776287165237590670585632324 + I*-.98954785079406926140421774107380770146846771240234;\nv123 = .68196836727835075109283025085460394620895385742188 + I*-.73138166919311053959518176270648837089538574218750;\np131 = .65481401335674405128628450256655924022197723388672 + I*.75579005544637445357381011490360833704471588134766;\nv131 = .66366404128464795419972688250709325075149536132812 + I*-.74803077497234626846989158366341143846511840820312;\np132 = -.38515023867666903223394569977244827896356582641602 + I*-.92285388531842083015277466984116472303867340087891;\nv132 = .97485914781713112198247017659014090895652770996094 + I*-.22282199603552815814211385259113740175962448120117;\np133 = .79050723801252975597719796496676281094551086425781 + I*-.61245269747940656035467554829665459692478179931641;\nv133 = -.55452096927409122439911470792139880359172821044922 + I*-.83216975109368296337208903423743322491645812988281;\np141 = -.70191820130700333102424792741658166050910949707031 + I*.71225756484149838065889071003766730427742004394531;\nv141 = .84778390295906769225098287279251962900161743164062 + I*.53034182739388957550374925631331279873847961425781;\np142 = -.70911523877004656313260966271627694368362426757812 + I*.70509260253111427640959618656779639422893524169922;\nv142 = .63357068188932053054429616167908534407615661621094 + I*-.77368481376481812450407460346468724310398101806641;\np143 = .13372318114258330279398023776593618094921112060547 + I*-.99101872375102373347033335448941215872764587402344;\nv143 = -.85242882424413080055813907165429554879665374755859 + I*-.52284328397500590135393849777756258845329284667969;\n");
;
  fprintf(OUT, "p151 = -.76602514768138119105600480907014571130275726318359 + I*-.64281060439270765183294997768825851380825042724609;\nv151 = -.95293578205969198258173946669558063149452209472656 + I*-.30317222047919084593203820077178534120321273803711;\np152 = -.74542607526671844642152109372545965015888214111328 + I*.66658830346208197692448038651491515338420867919922;\nv152 = .99978912190881030763733861022046767175197601318359 + I*-.20535620584969226792848573381888854783028364181519e-1;\np153 = -.66082270020815792044288627948844805359840393066406e-1 + I*-.99781417788528947721005124549265019595623016357422;\nv153 = .95051891242748354216018924489617347717285156250000 + I*-.31066669779310723820131556749402079731225967407227;\np161 = -.81121955050442173806857226736610755324363708496094 + I*.58474168731107589724871331782196648418903350830078;\nv161 = -.39258761698952898688830259743554051965475082397461 + I*-.91971460953193662213323023024713620543479919433594;\np162 = -.23163329922465047960855599740170873701572418212891 + I*-.97280317366377022647583316938835196197032928466797;\nv162 = -.63982647014457916245788737796829082071781158447266 + I*.76851941296386783175620394104043953120708465576172;\np163 = .98667790950555145368383591630845330655574798583984 + I*.16268590256612514011180792294908314943313598632812;\nv163 = -.85655771575413686136357682698871940374374389648438 + I*-.51605123736122882061039263135171495378017425537109;\n");

  // new constant values
  for (j = 1; j <= 6; j++)
  {
    for (k = 1; k <= 3; k++)
    { // generate random real numbers in [-1,1] and print them
      fprintf(OUT, "p0%d%d = ", j, k);
      mpf_set_q(temp, A[j-1][k-1]);
      mpf_out_str(OUT, 10, 0, temp);
      fprintf(OUT, ";\nv0%d%d = ", j, k);
      mpf_set_q(temp, B[j-1][k-1]);
      mpf_out_str(OUT, 10, 0, temp);
      fprintf(OUT, ";\n");
    }
  }

  // setup homotopy between old and new constant values
  fprintf(OUT, "x1 = (s*p111 + one_minus_s*p011) + t1*(s*v111 + one_minus_s*v011);\ny1 = (s*p112 + one_minus_s*p012) + t1*(s*v112 + one_minus_s*v012);\nz1 = (s*p113 + one_minus_s*p013) + t1*(s*v113 + one_minus_s*v013);\nx2 = (s*p121 + one_minus_s*p021) + t2*(s*v121 + one_minus_s*v021);\ny2 = (s*p122 + one_minus_s*p022) + t2*(s*v122 + one_minus_s*v022);\nz2 = (s*p123 + one_minus_s*p023) + t2*(s*v123 + one_minus_s*v023);\nx3 = (s*p131 + one_minus_s*p031) + t3*(s*v131 + one_minus_s*v031);\ny3 = (s*p132 + one_minus_s*p032) + t3*(s*v132 + one_minus_s*v032);\nz3 = (s*p133 + one_minus_s*p033) + t3*(s*v133 + one_minus_s*v033);\nx4 = (s*p141 + one_minus_s*p041) + t4*(s*v141 + one_minus_s*v041);\ny4 = (s*p142 + one_minus_s*p042) + t4*(s*v142 + one_minus_s*v042);\nz4 = (s*p143 + one_minus_s*p043) + t4*(s*v143 + one_minus_s*v043);\nx5 = (s*p151 + one_minus_s*p051) + t5*(s*v151 + one_minus_s*v051);\ny5 = (s*p152 + one_minus_s*p052) + t5*(s*v152 + one_minus_s*v052);\nz5 = (s*p153 + one_minus_s*p053) + t5*(s*v153 + one_minus_s*v053);\nx6 = (s*p161 + one_minus_s*p061) + t6*(s*v161 + one_minus_s*v061);\ny6 = (s*p162 + one_minus_s*p062) + t6*(s*v162 + one_minus_s*v062);\nz6 = (s*p163 + one_minus_s*p063) + t6*(s*v163 + one_minus_s*v063);\n");
  // print the functions
  fprintf(OUT, "f1 = f*(x1+1) + g*y1 - z1;\nf2 = f*(x2+1) + g*y2 - z2;\nf3 = f*(x3+1) + g*y3 - z3;\nf4 = f*(x4+1) + g*y4 - z4;\nf5 = f*(x5+1) + g*y5 - z5;\nf6 = f*(x6+1) + g*y6 - z6;\n");
  fprintf(OUT, "f7 = a*x1^2 + b*x1*y1 + c*y1^2 + (a+1)*x1 + e*y1 + 1;\nf8 = a*x2^2 + b*x2*y2 + c*y2^2 + (a+1)*x2 + e*y2 + 1;\nf9 = a*x3^2 + b*x3*y3 + c*y3^2 + (a+1)*x3 + e*y3 + 1;\nf10= a*x4^2 + b*x4*y4 + c*y4^2 + (a+1)*x4 + e*y4 + 1;\nf11= a*x5^2 + b*x5*y5 + c*y5^2 + (a+1)*x5 + e*y5 + 1;\nf12= a*x6^2 + b*x6*y6 + c*y6^2 + (a+1)*x6 + e*y6 + 1;\n");

  // close the file
  fprintf(OUT, "END;");
  fclose(OUT);

  mpf_clear(temp);

  return;
}

void createCertifyFile(char *name, mpq_t **p, mpq_t **v)
{
  int j, k;
  mpq_t temp1, temp2;
  FILE *OUT = fopen(name, "w");
  mpq_init(temp1);  mpq_init(temp2); 

  fprintf(OUT, "12 12\n\n");

  for (j = 0; j < 6; j++)
  {
    fprintf(OUT, "6\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 0 0 "); mpq_out_str(OUT, 10, p[j][2]); fprintf(OUT, " 0\n");
    for (k = 0; k < 6; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 0 0 0 1 "); mpq_neg(temp1, v[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 6; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 0 0 1 0 "); mpq_neg(temp1, v[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 6; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 0 0 0 0 "); mpq_out_str(OUT, 10, v[j][2]); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 1 0 "); mpq_set_ui(temp1, 1, 1); mpq_add(temp1, p[j][0], temp1); mpq_neg(temp1, temp1); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 0 1 "); mpq_neg(temp1, p[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");

    fprintf(OUT, "13\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 0 0 "); mpq_set_ui(temp1, 1, 1); mpq_add(temp1, p[j][0], temp1); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 6; k++) fprintf(OUT, "%d ", 2*(j == k));
    fprintf(OUT, "1 0 0 0 0 0 "); mpq_mul(temp1, v[j][0], v[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 6; k++) fprintf(OUT, "%d ", 2*(j == k));
    fprintf(OUT, "0 0 1 0 0 0 "); mpq_mul(temp1, v[j][1], v[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 6; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 0 1 0 0 "); mpq_out_str(OUT, 10, v[j][1]); fprintf(OUT, " 0\n");
    for (k = 0; k < 6; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "1 0 0 0 0 0 "); mpq_mul(temp1, p[j][0], v[j][0]); mpq_add(temp1, temp1, temp1); mpq_add(temp1, temp1, v[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 6; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 1 0 0 0 0 "); mpq_mul(temp1, p[j][0], v[j][1]); mpq_mul(temp2, v[j][0], p[j][1]); mpq_add(temp1, temp1, temp2); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 6; k++) fprintf(OUT, "%d ", 2*(j == k));
    fprintf(OUT, "0 1 0 0 0 0 "); mpq_mul(temp1, v[j][0], v[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 6; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 1 0 0 0 "); mpq_mul(temp1, p[j][1], v[j][1]); mpq_add(temp1, temp1, temp1); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 6; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 0 0 0 0 "); mpq_out_str(OUT, 10, v[j][0]); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 1 0 0 0 0 0 "); mpq_mul(temp1, p[j][0], p[j][0]); mpq_add(temp1, temp1, p[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 1 0 0 0 0 "); mpq_mul(temp1, p[j][0], p[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 1 0 0 0 "); mpq_mul(temp1, p[j][1], p[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 1 0 0 "); mpq_out_str(OUT, 10, p[j][1]); fprintf(OUT, " 0\n");
  }

  fclose(OUT);
  mpq_clear(temp1);
  mpq_clear(temp2);

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


