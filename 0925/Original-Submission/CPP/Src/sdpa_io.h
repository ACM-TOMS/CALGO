#ifndef __sdpa_io_h__
#define __sdpa_io_h__

#define lengthOfString 256

#include "sdpa_block.h"
#include "sdpa_parts.h"
  
namespace sdpa {

class IO
{
public:
  static void read(FILE* fpData, FILE* fpout, int& m, char* str);
  static void read(FILE* fpData, int& nBlock);
  static void read(FILE* fpData, BlockStruct& bs);
  static void read(FILE* fpData, Vector& b);
  static void read(FILE* fpData, DenseLinearSpace& xMat,
		   Vector& yVec, DenseLinearSpace& zMat,
		   BlockStruct& bs, bool inputSparse);
  static void read(FILE* fpData, int m,
		   BlockStruct& bs,
		   InputData& inputData, bool isDataSparse);


  // 2008/02/27 kazuhide nakata   
  // without LP_ANonZeroCount
  static void setBlockStruct(FILE* fpData, InputData& inputData,
			     int m, BlockStruct& bs,
                             long position, bool isDataSparse);
  
  // 2008/02/27 kazuhide nakata   
  // without LP_ANonZeroCount
  static void setElement(FILE* fpData, InputData& inputData, int m,
			 BlockStruct& bs,
                         long position, bool isDataSparse);

  static void printHeader(FILE* fpout, FILE* Display);

  static void printOneIteration(int pIteration,
				AverageComplementarity& mu,
				RatioInitResCurrentRes& theta,
				SolveInfo& solveInfo,
				StepLength& alpha,
				DirectionParameter& beta,
				FILE* fpout,
				FILE* Display);

  static void printLastInfo(int pIteration,
			    AverageComplementarity& mu,
			    RatioInitResCurrentRes& theta,
			    SolveInfo& solveInfo,
			    StepLength& alpha,
			    DirectionParameter& beta,
			    Residuals& currentRes,
			    Phase & phase,
			    Solutions& currentPt,
			    InputData& inputData,
                            WorkVariables& work,
			    double cputime,
			    ComputeTime& com,
			    Parameter& param,
			    FILE* fpout,
			    FILE* Display,
			    bool printTime = true);

  static void computeDimacs(double* dimacs_error,
			    SolveInfo& solveInfo,
			    Residuals& currentRes,
			    Solutions& currentPt,
			    InputData& inputData,
			    WorkVariables& work);
  
  static void printDimacs(double* dimacs_error,
			  char* printFormat,
			  FILE* fpout);

  static void printSolution(BlockStruct& bs, Solutions& currentPt,
			    Parameter& param, FILE* fpout);

};

} // end of namespace 'sdpa'

#endif // __sdpa_io_h__
