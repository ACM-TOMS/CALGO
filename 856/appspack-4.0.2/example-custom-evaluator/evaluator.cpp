// $Id: evaluator.cpp,v 1.1.2.1 2005/06/29 00:01:32 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/example-custom-evaluator/Attic/evaluator.cpp,v $ 

#include "evaluator.hpp"

CustomEvaluator::CustomEvaluator(): cnt(0)
{
}

CustomEvaluator::~CustomEvaluator()
{
  print();
}

void CustomEvaluator::operator()(int tag_in, const APPSPACK::Vector& x_in, 
				 bool& isF_out, double& f_out, string& msg_out)
{
  if (constraint(x_in))
  {
    isF_out = true;
    f_out = feval(x_in);
    msg_out = "Success";
  }
  else
  {
    isF_out = false;
    msg_out = "Constraint Violation";
  }
  cnt ++;
  return;
}

void CustomEvaluator::print() const
{
  cout << "Number of evaluations is: " << cnt << endl;
}

double CustomEvaluator::feval(const APPSPACK::Vector& x)
{
  double f = 0;

  for (int i = 0; i < x.size(); i ++)
    f += (i + 1) * x[i] * x[i];

  return(f);
} 

bool CustomEvaluator::constraint(const APPSPACK::Vector& x)
{
  double tmp = 0;

  for (int i = 0; i < x.size(); i ++)
    tmp += x[i] * x[i];

  return (tmp >= 1.0);
} 
