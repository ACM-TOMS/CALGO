#ifndef VBF_tools__H
#define VBF_tools__H

#include <set>

NTL_CLIENT

long IsDigit(long a)
{
   if (a > NTL_MAX_INT || a < NTL_MIN_INT)
      return 0;

   int b = (int) a;

   if (isdigit(b))
      return 1;
   else 
      return 0;
}

long IsValid(long a)
{

   if (a == '+' || a == '\n')
      return 1;

   int b = (int) a;

   if (isdigit(b) || isspace(b))
      return 1;
   else 
      return 0;
}

void Tokenize(const NTL_SNS string& str,
              NTL_SNS vector<NTL_SNS string>& tokens,
              const NTL_SNS string& delimiters = "+")
{
    // Skip delimiters at beginning.
    NTL_SNS string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    NTL_SNS string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (NTL_SNS string::npos != pos || NTL_SNS string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void PickMin(long& min, NTL_SNS set<long>& S)
{
    NTL_SNS set<long>::iterator it;

    min = -1;

    for (it=S.begin(); it!=S.end(); it++)
    {
	if (min == -1) 
	{
	   min = *it;
	} else {
 	   if (*it < min) min = *it;
        }
    }
}

long Factorial(long val)
{
  long Result = 1;

  for(long i = 2; i <= val; ++i)
  {
    Result *= i;
  }

  return Result;
}

long Combination(long N, long R)
{
   return (Factorial(N) / (Factorial(N-R) * Factorial(R)));
}

#endif
