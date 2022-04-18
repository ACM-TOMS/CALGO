// $Id: evaluator.hpp,v 1.1.2.1 2005/06/29 00:01:32 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/example-custom-evaluator/Attic/evaluator.hpp,v $ 

#include "APPSPACK_Evaluator_Interface.hpp"

class CustomEvaluator : public APPSPACK::Evaluator::Interface
{

public:

  // Constructor - can take any arguments that are appropriate. In
  // this case, we don't have any.
  CustomEvaluator();

  // Destructor
  ~CustomEvaluator();

  // Evaluates the function at x, returning the function value or any
  // error messages. It must fill in isF_out, f_out, and msg_out. For
  // successful evaluations, msg_out should be something like
  // "Success".
  void operator()(int tag_in, const APPSPACK::Vector& x_in, 
		  bool& isF_out, double& f_out, string& msg_out);

  // Prints out information about the evaluator
  void print() const;

private:
  
  // Return true is the constraint is satisfied, false otherwise.
  bool constraint(const APPSPACK::Vector& x);

  // Evaluate function at x
  double feval(const APPSPACK::Vector& x);

  // Number of calls to the evaluator
  int cnt;

};
