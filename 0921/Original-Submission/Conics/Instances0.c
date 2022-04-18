/****************************************************************

Generate a random instance for solving for the conics that meets
8 lines.

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
  createInputFiles("Conic0.bertini", "Conic0.poly");

  return 0;
}

void createInputFiles(char *inputName, char *certifyName)
{ // create the input files
  int i, j;
  char *str = NULL;

  // allocate and initialize memory
  mpq_t **A = (mpq_t **)malloc(8 * sizeof(mpq_t *)), **B = (mpq_t **)malloc(8 * sizeof(mpq_t *));
  for (i = 0; i < 8; i++)
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
  for (i = 0; i < 8; i++)
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
  fprintf(OUT, "INPUT\nvariable t1,t2,t3,t4,t5,t6,t7,t8,a,b,c,d,e,f,g,h;\nfunction f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16;\nconstant p111,v111,p112,v112,p113,v113,p121,v121,p122,v122,p123,v123,p131,v131,p132,v132,p133,v133,p141,v141,p142,v142,p143,v143,p151,v151,p152,v152,p153,v153,p161,v161,p162,v162,p163,v163,p171,v171,p172,v172,p173,v173,p181,v181,p182,v182,p183,v183;\nconstant p011,v011,p012,v012,p013,v013,p021,v021,p022,v022,p023,v023,p031,v031,p032,v032,p033,v033,p041,v041,p042,v042,p043,v043,p051,v051,p052,v052,p053,v053,p061,v061,p062,v062,p063,v063,p071,v071,p072,v072,p073,v073,p081,v081,p082,v082,p083,v083;\npathvariable t;\nparameter s;\ns = t;\n");
  // old contant values
  fprintf(OUT, "p111 = -.43430128096558950900885065493639558553695678710938 + I*-.90076767112926969804931331964326091110706329345703;\nv111 = .84169583859478780407670228669303469359874725341797 + I*-.53995195646669991162980295484885573387145996093750;\np112 = .85029736124152752174687464048474794253706932067871e-1 + I*-.99637841404491345187466322386171668767929077148438;\nv112 = .39522265332853717678673888258344959467649459838867 + I*-.91858535493221893375448416918516159057617187500000;\np113 = .63143726363749130836566791913355700671672821044922 + I*.77542696760558782465011518070241436362266540527344;\nv113 = .98071111717300918364514927816344425082206726074219 + I*-.19546279608475009004209255181194748729467391967773;\np121 = -.83797690901853982836655632127076387405395507812500 + I*-.54570568986563960933722228219266980886459350585938;\nv121 = .76532301540561686881147807071101851761341094970703 + I*-.64364639522835342955175974566373042762279510498047;\np122 = .99314670367952484486551156805944629013538360595703 + I*.11687439826794455977054809636683785356581211090088;\nv122 = .76185598631517148504599390435032546520233154296875 + I*-.64774644430960592877966064406791701912879943847656;\np123 = -.14420489238870623505128776287165237590670585632324 + I*-.98954785079406926140421774107380770146846771240234;\nv123 = .68196836727835075109283025085460394620895385742188 + I*-.73138166919311053959518176270648837089538574218750;\np131 = .65481401335674405128628450256655924022197723388672 + I*.75579005544637445357381011490360833704471588134766;\nv131 = .66366404128464795419972688250709325075149536132812 + I*-.74803077497234626846989158366341143846511840820312;\np132 = -.38515023867666903223394569977244827896356582641602 + I*-.92285388531842083015277466984116472303867340087891;\nv132 = .97485914781713112198247017659014090895652770996094 + I*-.22282199603552815814211385259113740175962448120117;\np133 = .79050723801252975597719796496676281094551086425781 + I*-.61245269747940656035467554829665459692478179931641;\nv133 = -.55452096927409122439911470792139880359172821044922 + I*-.83216975109368296337208903423743322491645812988281;\np141 = -.70191820130700333102424792741658166050910949707031 + I*.71225756484149838065889071003766730427742004394531;\nv141 = .84778390295906769225098287279251962900161743164062 + I*.53034182739388957550374925631331279873847961425781;\np142 = -.70911523877004656313260966271627694368362426757812 + I*.70509260253111427640959618656779639422893524169922;\nv142 = .63357068188932053054429616167908534407615661621094 + I*-.77368481376481812450407460346468724310398101806641;\np143 = .13372318114258330279398023776593618094921112060547 + I*-.99101872375102373347033335448941215872764587402344;\nv143 = -.85242882424413080055813907165429554879665374755859 + I*-.52284328397500590135393849777756258845329284667969;\n");
  fprintf(OUT, "p151 = -.76602514768138119105600480907014571130275726318359 + I*-.64281060439270765183294997768825851380825042724609;\nv151 = -.95293578205969198258173946669558063149452209472656 + I*-.30317222047919084593203820077178534120321273803711;\np152 = -.74542607526671844642152109372545965015888214111328 + I*.66658830346208197692448038651491515338420867919922;\nv152 = .99978912190881030763733861022046767175197601318359 + I*-.20535620584969226792848573381888854783028364181519e-1;\np153 = -.66082270020815792044288627948844805359840393066406e-1 + I*-.99781417788528947721005124549265019595623016357422;\nv153 = .95051891242748354216018924489617347717285156250000 + I*-.31066669779310723820131556749402079731225967407227;\np161 = -.81121955050442173806857226736610755324363708496094 + I*.58474168731107589724871331782196648418903350830078;\nv161 = -.39258761698952898688830259743554051965475082397461 + I*-.91971460953193662213323023024713620543479919433594;\np162 = -.23163329922465047960855599740170873701572418212891 + I*-.97280317366377022647583316938835196197032928466797;\nv162 = -.63982647014457916245788737796829082071781158447266 + I*.76851941296386783175620394104043953120708465576172;\np163 = .98667790950555145368383591630845330655574798583984 + I*.16268590256612514011180792294908314943313598632812;\nv163 = -.85655771575413686136357682698871940374374389648438 + I*-.51605123736122882061039263135171495378017425537109;\np171 = -.41578287056849805303215816820738837122917175292969 + I*.90946391052191810633331670032930560410022735595703;\nv171 = -.72753629796411400931788193702232092618942260742188 + I*-.68606919122248299913735536392778158187866210937500;\np172 = -.91175548765879310675330771118751727044582366943359 + I*.41073340590226697921494292131683323532342910766602;\nv172 = .84318935492661795105817645890056155622005462646484 + I*.53761669592604177658756725577404722571372985839844;\np173 = -.72215398704856292866338662861380726099014282226562 + I*.69173233189570093593090405192924663424491882324219;\nv173 = -.63975520834612931153628778702113777399063110351562 + I*.76857873597569736912049620514153502881526947021484;\np181 = -.85731236922212716855540293181547895073890686035156 + I*.51479656329344625209643027119454927742481231689453;\nv181 = -.92344419129425636683095035550650209188461303710938 + I*-.38373275279144059002334188335225917398929595947266;\np182 = .40905418411576688342279339849483221769332885742188 + I*.91251009553724071743374679499538615345954895019531;\nv182 = -.69829061679580328725336357820197008550167083740234 + I*.71581437153422433627980581150040961802005767822266;\np183 = .99885446575784220257787637820001691579818725585938 + I*-.47851397425944866559355261870223330333828926086426e-1;\nv183 = -.76960576044917805571543567566550336778163909912109 + I*-.63851936030432321356187230776413343846797943115234;\n");

  // new constant values
  for (j = 1; j <= 8; j++)
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

  // print the functions
  fprintf(OUT, "one_minus_s = 1-s;\n");
  fprintf(OUT, "x1 = (p011 + t1*v011)*one_minus_s + s*(p111 + t1*v111);\ny1 = (p012 + t1*v012)*one_minus_s + s*(p112 + t1*v112);\nz1 = (p013 + t1*v013)*one_minus_s + s*(p113 + t1*v113);\nx2 = (p021 + t2*v021)*one_minus_s + s*(p121 + t2*v121);\ny2 = (p022 + t2*v022)*one_minus_s + s*(p122 + t2*v122);\nz2 = (p023 + t2*v023)*one_minus_s + s*(p123 + t2*v123);\nx3 = (p031 + t3*v031)*one_minus_s + s*(p131 + t3*v131);\ny3 = (p032 + t3*v032)*one_minus_s + s*(p132 + t3*v132);\nz3 = (p033 + t3*v033)*one_minus_s + s*(p133 + t3*v133);\nx4 = (p041 + t4*v041)*one_minus_s + s*(p141 + t4*v141);\ny4 = (p042 + t4*v042)*one_minus_s + s*(p142 + t4*v142);\nz4 = (p043 + t4*v043)*one_minus_s + s*(p143 + t4*v143);\nx5 = (p051 + t5*v051)*one_minus_s + s*(p151 + t5*v151);\ny5 = (p052 + t5*v052)*one_minus_s + s*(p152 + t5*v152);\nz5 = (p053 + t5*v053)*one_minus_s + s*(p153 + t5*v153);\nx6 = (p061 + t6*v061)*one_minus_s + s*(p161 + t6*v161);\ny6 = (p062 + t6*v062)*one_minus_s + s*(p162 + t6*v162);\nz6 = (p063 + t6*v063)*one_minus_s + s*(p163 + t6*v163);\nx7 = (p071 + t7*v071)*one_minus_s + s*(p171 + t7*v171);\ny7 = (p072 + t7*v072)*one_minus_s + s*(p172 + t7*v172);\nz7 = (p073 + t7*v073)*one_minus_s + s*(p173 + t7*v173);\nx8 = (p081 + t8*v081)*one_minus_s + s*(p181 + t8*v181);\ny8 = (p082 + t8*v082)*one_minus_s + s*(p182 + t8*v182);\nz8 = (p083 + t8*v083)*one_minus_s + s*(p183 + t8*v183);\n");
  fprintf(OUT, "f1 = z1 - (f*x1 + g*y1 + h);\nf2 = z2 - (f*x2 + g*y2 + h);\nf3 = z3 - (f*x3 + g*y3 + h);\nf4 = z4 - (f*x4 + g*y4 + h);\nf5 = z5 - (f*x5 + g*y5 + h);\nf6 = z6 - (f*x6 + g*y6 + h);\nf7 = z7 - (f*x7 + g*y7 + h);\nf8 = z8 - (f*x8 + g*y8 + h);\n");
  fprintf(OUT, "f9 = (a*x1 + b*y1 + d)*x1 + (c*y1 + e)*y1 + 1;\nf10= (a*x2 + b*y2 + d)*x2 + (c*y2 + e)*y2 + 1;\nf11= (a*x3 + b*y3 + d)*x3 + (c*y3 + e)*y3 + 1;\nf12= (a*x4 + b*y4 + d)*x4 + (c*y4 + e)*y4 + 1;\nf13= (a*x5 + b*y5 + d)*x5 + (c*y5 + e)*y5 + 1;\nf14= (a*x6 + b*y6 + d)*x6 + (c*y6 + e)*y6 + 1;\nf15= (a*x7 + b*y7 + d)*x7 + (c*y7 + e)*y7 + 1;\nf16= (a*x8 + b*y8 + d)*x8 + (c*y8 + e)*y8 + 1;\n");

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

  fprintf(OUT, "16 16\n\n");

  for (j = 0; j < 8; j++)
  {
    fprintf(OUT, "7\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "); mpq_out_str(OUT, 10, p[j][2]); fprintf(OUT, " 0\n");
    for (k = 0; k < 8; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 0 0 0 0 0 0 "); mpq_out_str(OUT, 10, v[j][2]); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 "); mpq_neg(temp1, p[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 "); mpq_neg(temp1, p[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 "); mpq_set_ui(temp1, 1, 1); mpq_neg(temp1, temp1); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 8; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 0 0 0 0 1 0 "); mpq_neg(temp1, v[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 8; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 0 0 0 1 0 0 "); mpq_neg(temp1, v[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");

    fprintf(OUT, "14\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "); mpq_set_ui(temp1, 1, 1); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 "); mpq_mul(temp1, p[j][0], p[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 "); mpq_mul(temp1, p[j][0], p[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 "); mpq_mul(temp1, p[j][1], p[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 "); mpq_out_str(OUT, 10, p[j][0]); fprintf(OUT, " 0\n");
    fprintf(OUT, "0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 "); mpq_out_str(OUT, 10, p[j][1]); fprintf(OUT, " 0\n");
    for (k = 0; k < 8; k++) fprintf(OUT, "%d ", 2*(j == k));
    fprintf(OUT, "1 0 0 0 0 0 0 0 "); mpq_mul(temp1, v[j][0], v[j][0]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 8; k++) fprintf(OUT, "%d ", 2*(j == k));
    fprintf(OUT, "0 0 1 0 0 0 0 0 "); mpq_mul(temp1, v[j][1], v[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 8; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 0 1 0 0 0 0 "); mpq_out_str(OUT, 10, v[j][0]); fprintf(OUT, " 0\n");
    for (k = 0; k < 8; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 0 0 1 0 0 0 "); mpq_out_str(OUT, 10, v[j][1]); fprintf(OUT, " 0\n");
    for (k = 0; k < 8; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "1 0 0 0 0 0 0 0 "); mpq_mul(temp1, p[j][0], v[j][0]); mpq_add(temp1, temp1, temp1); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 8; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 1 0 0 0 0 0 0 "); mpq_mul(temp1, p[j][0], v[j][1]); mpq_mul(temp2, p[j][1], v[j][0]); mpq_add(temp1, temp1, temp2); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 8; k++) fprintf(OUT, "%d ", 2*(j == k));
    fprintf(OUT, "0 1 0 0 0 0 0 0 "); mpq_mul(temp1, v[j][0], v[j][1]); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
    for (k = 0; k < 8; k++) fprintf(OUT, "%d ", j == k);
    fprintf(OUT, "0 0 1 0 0 0 0 0 "); mpq_mul(temp1, p[j][1], v[j][1]); mpq_add(temp1, temp1, temp1); mpq_out_str(OUT, 10, temp1); fprintf(OUT, " 0\n");
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


