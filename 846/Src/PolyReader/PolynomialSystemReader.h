//PolynomialSystemReader.h
//Author: Xing Li
//Date: 02/24/2000


#ifndef INCLUDE_POLYNOMIALSYSTEMREADER_H
#define INCLUDE_POLYNOMIALSYSTEMREADER_H

#include "ExpandedPolynomial.h"

#include <string>
#include <vector>
#include <complex>
#include <iostream>


class PolynomialSystemReader
{
 public: 
  PolynomialSystemReader(istream&);
  void print();
 
  void getSupportSize(int *ss) const;
  void getCoefficient(complex<double> *coef) const;
  void getSupports(int** sp) const;
  //void createSupportCoefficientFiles(const char*, const char *) const;
  void GetSupport(int* SptSize, int** Spt) const;
  void GetDimensions(int &nVar,int &nPts) const;

  int dimension() const
    { return m_variable.size(); }

  int equationNumber() const
    { return m_poly.size(); }

  const vector<string>& variable() const 
  { return m_variable; }

 private:
  PolynomialSystemReader();
  PolynomialSystemReader(const PolynomialSystemReader&);
  PolynomialSystemReader& operator=(const PolynomialSystemReader&);

 private:
  vector<ExpandedPolynomial> m_poly;
  vector<string> m_variable;
};

#endif //INCLUDE_PolynomialSystemReader_H
