
#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>
#include "sdpa_mpicopy.h"
using namespace sdpa;

#define DEFAULT_PARAMETER_FILE ((char *)"./param.sdpa")

static void message(char* argv0)
{
  cout << endl;
  cout << "*** Please assign data file and output file.***" << endl;
  cout << endl;
  cout << "---- option type 1 ------------" << endl;
  cout << argv0 <<" DataFile OutputFile [InitialPtFile]"
    " [-pt parameters] [-dimacs] [-numThreads numThreads]"<< endl;
  cout << "parameters = 0 default, 1 fast (unstable),"
    " 2 slow (stable)" << endl;
  cout << "  -dimacs : printout dimacs information incurring additional computation cost " << endl;
  cout << "  -numThreads: Number of pthreads for internal computation" << endl;
  cout << "example1-1: " << argv0
       << " example1.dat example1.result" << endl;
  cout << "example1-2: " << argv0
       << " example1.dat-s example1.result" << endl;
  cout << "example1-3: " << argv0
       << " example1.dat example1.result example1.ini" << endl;
  cout << "example1-4: " << argv0
       << " example1.dat example1.result -pt 2" << endl;
  cout << "example1-5: " << argv0
       << " example1.dat example1.result -dimacs" << endl;
  cout << "example1-6: " << argv0
       << " example1.dat example1.result -numThreads 4" << endl;

  cout << endl;
  cout << "---- option type 2 ------------" << endl;
  cout << argv0 << " [option filename]+ " << endl;
  cout << "  -dd : data dense :: -ds : data sparse     " << endl;
  cout << "  -id : init dense :: -is : init sparse     " << endl;
  cout << "  -o  : output     :: -p  : parameter       " << endl;
  cout << "  -pt : parameters , 0 default, 1 fast (unstable)" << endl;
  cout << "                     2 slow (stable)         " << endl;
  cout << "  -dimacs : printout dimacs information incurring additional computation cost " << endl;
  cout << "  -numThreads: Number of pthreads for internal computation" << endl;
  cout << "example2-1: " << argv0
       << " -o example1.result -dd example1.dat" << endl;
  cout << "example2-2: " << argv0
       << " -ds example1.dat-s -o example1.result "
       << "-p param.sdpa" << endl;
  cout << "example2-3: " << argv0
       << " -ds example1.dat-s -o example1.result "
       << "-pt 2" << endl;
  cout << "example2-4: " << argv0
       << " -ds example1.dat-s -o example1.result "
       << "-dimacs" << endl;
  cout << "example2-5: " << argv0
       << " -ds example1.dat-s -o example1.result "
       << "-numThreads 4" << endl;
}

static void argumentAnalysis(SDPA& Problem1,
			     int argc, char** argv,
			     char*& inputFileName,
			     char*& resultFileName,
			     char*& initFileName,
			     char*& paramFileName,
			     SDPA::SparseType& isInputSparse,
			     SDPA::SparseType& isInitSparse,
			     SDPA::ParameterType& parameterType,
			     bool& isDimacs, int& numThreads)
{
  if (argv[1][0] == '-') {
    // rsdpa argument
    
    for (int index = 0; index < argc; ++index) {
      char* target = argv[index];
      if (strcmp(target,"-dd")==0 && index+1 < argc) {
	inputFileName = argv[index+1];
	isInputSparse = SDPA::DENSE;
	index++;
	continue;
      }
      if (strcmp(target,"-ds")==0 && index+1 < argc) {
	inputFileName = argv[index+1];
	isInputSparse = SDPA::SPARSE;
	continue;
      }
      if (strcmp(target,"-id")==0 && index+1 < argc) {
	initFileName = argv[index+1];
	isInitSparse = SDPA::DENSE;
	continue;
      }
      if (strcmp(target,"-is")==0 && index+1 < argc) {
	initFileName = argv[index+1];
	isInitSparse = SDPA::SPARSE;
	index++;
	continue;
      }
      if (strcmp(target,"-o")==0 && index+1 < argc) {
	resultFileName = argv[index+1];
	index++;
	continue;
      }
      if (strcmp(target,"-p")==0 && index+1 < argc) {
	paramFileName = argv[index+1];
	index++;
	continue;
      }
      if (strcmp(target,"-k")==0 && index+1 < argc) {
	double KAPPA = atof(argv[index+1]);
	Problem1.setKappa(KAPPA);
	if (MpiSt::iam == 0) {
	  rMessage("Kappa = " << KAPPA);
	}
	index++;
	continue;
      }
      if (strcmp(target,"-dimacs")==0) {
	isDimacs = true;
	continue;
      }
      if (strcmp(target,"-pt")==0 && index+1 < argc) {
	int tmp = atoi(argv[index+1]);
	switch (tmp) {
	case 0:
	  parameterType = SDPA::PARAMETER_DEFAULT;
	  break;
	case 1:
	  parameterType = SDPA::PARAMETER_UNSTABLE_BUT_FAST;
	  break;
	case 2:
	  parameterType = SDPA::PARAMETER_STABLE_BUT_SLOW;
	  break;
	default:
	  parameterType = SDPA::PARAMETER_DEFAULT;
	}
	index++;
	paramFileName = NULL;
	continue;
      }
      if (strcmp(target,"-numThreads")==0 && index+1 < argc) {
	numThreads = atoi(argv[index+1]);
	index++;
	continue;
      }
    }
  }
  else { // SDPA argument
    inputFileName = argv[1];
    int len = strlen(inputFileName);
    if (inputFileName[len-1] == 's'
	&& inputFileName[len-2] == '-') {
      isInputSparse = SDPA::SPARSE;
    }
	
    resultFileName = argv[2];

    paramFileName = DEFAULT_PARAMETER_FILE;

    for (int index=3; index<argc; ++index) {
      if (strcmp(argv[index],"-dimacs")==0) {
	isDimacs = true;
      }
      else if (strcmp(argv[index],"-numThreads")==0 && index+1 < argc) {
	numThreads = atoi(argv[index+1]);
	++index;
      }
      else if (strcmp(argv[index],"-pt")==0 && index+1 < argc) {
	int tmp = atoi(argv[index+1]);
	switch (tmp) {
	case 0:
	  parameterType = SDPA::PARAMETER_DEFAULT;
	  break;
	case 1:
	  parameterType = SDPA::PARAMETER_UNSTABLE_BUT_FAST;
	  break;
	case 2:
	  parameterType = SDPA::PARAMETER_STABLE_BUT_SLOW;
	  break;
	default:
	  parameterType = SDPA::PARAMETER_DEFAULT;
	}
	index++;
	paramFileName = NULL;
      } // end of "-pt"
      else {
	initFileName = argv[index];
	int len = strlen(initFileName);
	if (initFileName[len-1] == 's'
	    && initFileName[len-2] == '-') {
	  isInitSparse = SDPA::SPARSE;
	}
      }
    } // end of 'for'
  }
  
}

 
int main(int argc, char** argv)
{
  int info = 0;
  MpiSt::initialize(argc,argv);
  // MpiSt::display();
  TimeStart(ALL_START1);
  time_t ltime;
  time(&ltime);
  char string_time[1024];
  strcpy(string_time,ctime(&ltime));
  string_time[strlen(string_time)-1]='\0';
  if (MpiSt::iam == 0) {
    // cout << "let me see your ..." << endl;
    fprintf(stdout,"SDPA start at [%s]\n",string_time);
  }
  if (argc == 1) {
    if (MpiSt::iam == 0) {
      message(argv[0]);
    }
    MpiSt::finalize();
    return 0;
  }
  
  char* inputFileName  = NULL;
  char* resultFileName = NULL;
  char* initFileName   = NULL;
  char* paramFileName  = NULL;

  SDPA::SparseType isInputSparse = SDPA::DENSE;
  SDPA::SparseType isInitSparse  = SDPA::DENSE;
  SDPA::ParameterType parameterType = SDPA::PARAMETER_DEFAULT;

  SDPA Problem1;
  bool isDimacs = false;
  int  numThreads = 0; // 0 means automatic computation

  argumentAnalysis(Problem1, argc, argv,
		   inputFileName, resultFileName, initFileName,
		   paramFileName, isInputSparse, isInitSparse,
		   parameterType, isDimacs, numThreads);
  if (inputFileName == NULL || resultFileName == NULL) {
    if (MpiSt::iam == 0) {
      message(argv[0]);
    }
    MpiSt::finalize();
    return 0;
  }

  if (MpiSt::iam == 0) {
    Problem1.setDisplay(stdout);
  }
  else {
    Problem1.setDisplay(NULL);
    // for debug
    // Problem1.setDisplay(stdout);
  }
    
  FILE* fpresult = NULL;
  if (MpiSt::iam == 0) {
    if ((fpresult = fopen(resultFileName,"w")) == NULL) {
      rMessage("Cannot Open Result File : " << resultFileName);
      info = -1;
    }
    else {
      fprintf(fpresult,"SDPA start at [%s]\n",string_time);
      Problem1.setResultFile(fpresult);
    }
  }
  else {
    #if 1
    fpresult = NULL;
    Problem1.setResultFile(NULL);
    #else // for debug
    char RESULT_FILE_NAME[1024];
    rHere();
    sprintf(RESULT_FILE_NAME,"%s.%03d.result",resultFileName,MpiSt::iam);
    if ((fpresult = fopen(RESULT_FILE_NAME,"w")) == NULL) {
      rError("Cannot Open Result File : " << RESULT_FILE_NAME);
    }
    fprintf(fpresult,"SDPA start at [%s]\n",string_time);
    Problem1.setResultFile(fpresult);
    #endif
  }
  MpiCopy::allSendRecieveI(1,&info);
  if (info != 0) {
    Problem1.terminate();
    MpiSt::finalize();
    return 0;
  }
  
  if (paramFileName == NULL) {
    if (MpiSt::iam == 0) {
      if (parameterType == SDPA::PARAMETER_DEFAULT) {
	fprintf(stdout  ,"set   is DEFAULT\n");
	fprintf(fpresult,"set   is DEFAULT\n");
      }
      else if (parameterType == SDPA::PARAMETER_UNSTABLE_BUT_FAST) {
	fprintf(stdout  ,"set   is UNSTABLE_BUT_FAST\n");
	fprintf(fpresult,"set   is UNSTABLE_BUT_FAST\n");
      }
      else if (parameterType == SDPA::PARAMETER_STABLE_BUT_SLOW) {
	fprintf(stdout  ,"set   is STABLE_BUT_SLOW\n");
	fprintf(fpresult,"set   is STABLE_BUT_SLOW\n");
      }
    }
    Problem1.setParameterType(parameterType);
  }

  if (paramFileName) {
    if (MpiSt::iam == 0) {
      fprintf(stdout  ,"param is %s\n", paramFileName);
      info = Problem1.readParameter(paramFileName,fpresult);
    }
    MpiCopy::allSendRecieveI(1,&info);
    if (info != 0) {
      Problem1.terminate();
      MpiSt::finalize();
      return 0;
    }
    Problem1.mpiCopyParameter();
    // Problem1.param.display(stdout);
  }

  if (MpiSt::iam == 0) {
    fprintf(stdout  ,"data  is %s", inputFileName);
    if (isInputSparse == SDPA::SPARSE) {
      fprintf(stdout  ," : sparse\n");
    }
    else {
      fprintf(stdout  ," : dense\n");
    }
    info = Problem1.readInput(inputFileName, fpresult,
			      isInputSparse);
  }
  MpiCopy::allSendRecieveI(1,&info);
  if (info != 0) {
    Problem1.terminate();
    MpiSt::finalize();
    return 0;
  }
  Problem1.mpiCopyInput();
  
  // Data Copy
  
  if (initFileName) {
    Problem1.setInitPoint(true);
    if (MpiSt::iam == 0) {
      fprintf(stdout  ,"init  is %s", initFileName);
      if (isInitSparse == SDPA::SPARSE) {
	fprintf(stdout  ," : sparse\n");
      }
      else {
	fprintf(stdout  ," : dense\n");
      }
      info = Problem1.readInit(initFileName, fpresult, isInitSparse);
    }
    MpiCopy::allSendRecieveI(1,&info);
    if (info != 0) {
      Problem1.terminate();
      MpiSt::finalize();
      return 0;
    }
    Problem1.mpiCopyInit();
  }

  // MpiSt::barrier(); rHere(); MpiSt::barrier();
  if (MpiSt::iam == 0) {
    fprintf(stdout  ,"out   is %s\n", resultFileName);
    fprintf(fpresult,"out    is %s\n", resultFileName);
    if (isDimacs) {
      fprintf(stdout  ,"Dimacs information will be computed after the iteration.\n");
      fprintf(fpresult,"Dimacs information will be computed after the iteration.\n");
    }
  }

  Problem1.setNumThreads(numThreads);
  Problem1.initializeSolve();
  // MpiSt::barrier(); rHere(); MpiSt::barrier();
  Problem1.solve();

  if (MpiSt::iam==0 && isDimacs == true) {
    double dimacs_error[7];
    Problem1.getDimacsError(dimacs_error);
    Problem1.printDimacsError(dimacs_error,
			      Problem1.getParameterPrintInformation(),
			      stdout);
    Problem1.printDimacsError(dimacs_error,
			      Problem1.getParameterPrintInformation(),
			      fpresult);
  }
  
  Problem1.terminate();
  // MpiSt::barrier(); rHere(); MpiSt::barrier();

  TimeEnd(ALL_END1);
  double all_time = TimeCal(ALL_START1,ALL_END1);
  if (MpiSt::iam == 0) {
    time(&ltime);
    strcpy(string_time,ctime(&ltime));
    string_time[strlen(string_time)-1]='\0';
    fprintf(stdout  ,"SDPA end at [%s]\n",string_time);
    fprintf(stdout  ,"ALL TIME = %.6lf\n", all_time);
  }
  if (fpresult != NULL) {
    fprintf(fpresult,"SDPA end at [%s]\n",string_time);
    fprintf(fpresult,"ALL TIME = %.6lf\n", all_time);
    fclose(fpresult);
  }
  MpiSt::barrier();
  MpiSt::finalize();
  return 0;
}

