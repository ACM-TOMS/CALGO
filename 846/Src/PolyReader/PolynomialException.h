//polyException.h
//Author: Xing Li
//Data:  02/20/2000

#ifndef INCLUDE_POLYNOMIALEXCEPTION_H
#define INCLUDE_POLYNOMIALEXCEPTION_H
 
#include <string>

class PolynomialException {
public:
    enum errorType { ILLEGAL_SYMBOL, WRONG_PLACE_SYMBOL,  UNPAIRED_PAREN,
                NON_CONST_DINOMINATOR, NEGATIVE_EXPONENT};

    PolynomialException(errorType e, const std::string& s="");
   ~PolynomialException();

    PolynomialException& operator=(const PolynomialException&);
    
    std::string toString();

private:
    PolynomialException();

    errorType m_error;
    std::string m_symbol;
};

#endif //INCLUDE_POLYNOMIALEXCEPTION_H

