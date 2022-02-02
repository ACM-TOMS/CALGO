//PolynomialException.cpp
//Author: Xing Li
//Data:  02/20/2000

#include "PolynomialException.h"
 

PolynomialException::PolynomialException(errorType e, const std::string& s):
m_error(e), m_symbol(s)
{
}


PolynomialException:: ~PolynomialException()
{
}
  
PolynomialException& PolynomialException::operator=(const PolynomialException& e)
{
    m_error = e.m_error;
    m_symbol = e.m_symbol; 

    return *this;
}

std::string PolynomialException ::toString()
{
        switch( m_error)
        {
        case ILLEGAL_SYMBOL:
            return std::string("parse error:  illegal symbol \' ") + m_symbol + std::string("\'");

        case WRONG_PLACE_SYMBOL:
            return std::string("parse error: illegal usage of \'") + m_symbol+ std::string("\'");

        case UNPAIRED_PAREN:
            return std::string("parse error: unpaired parenthesis \'") + m_symbol+ std::string("\'");
           
        case NON_CONST_DINOMINATOR:
            return std::string("Non constant term cannot be in dinominator.\n");

        case NEGATIVE_EXPONENT:
            return std::string("Exponent cannot be negative.\n");
        }

        return std::string("It is impossible!");
}




