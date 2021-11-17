//ExpandedPolynomial.h
//Author: Xing Li
//Date: 03/02/2000

#ifndef INCLUDE_EXPANDEDPOLYNOMIAL_H
#define INCLUDE_EXPANDEDPOLYNOMIAL_H

#include<complex>
#include<utility>
#include<vector>
#include<string>

using namespace std;

class ExpandedPolynomial
{
 public:

  struct Pterm;
  typedef complex<double> ValueType;
  typedef std::pair<int, int> VarDegree; // <variable, degree>
  typedef std::vector<VarDegree> VarDegrees;
  typedef vector<Pterm> VPterm;

  struct Pterm 
  {
    Pterm(ValueType);
    Pterm(ValueType, const VarDegrees&);
    Pterm(const Pterm& pt);
    ~Pterm();

    Pterm& operator=(const Pterm& pt);
    bool operator<(const Pterm& p) const;
    bool operator>(const Pterm& p) const;
    bool operator==(const Pterm& p) const;
	bool operator!=(const Pterm& p) const;
    Pterm& operator*=(const Pterm& pt);
    Pterm operator*(const Pterm& pt);
	void print() const;

    ValueType m_coefficient;
    VarDegrees m_degrees;
  };

  ExpandedPolynomial(); 
  ExpandedPolynomial(const Pterm& p);
  ~ExpandedPolynomial();
  ExpandedPolynomial& operator+=(const Pterm& p);
  ExpandedPolynomial operator+(const Pterm& p);
  ExpandedPolynomial& operator+=(const ExpandedPolynomial& sp);
  ExpandedPolynomial operator+(const ExpandedPolynomial& sp);
  ExpandedPolynomial& operator-=(const Pterm& p);
  ExpandedPolynomial operator-(const Pterm& p);
  ExpandedPolynomial& operator-=(const ExpandedPolynomial& sp);
  ExpandedPolynomial operator-(const ExpandedPolynomial& sp);
  ExpandedPolynomial& operator*=(const Pterm& p);
  ExpandedPolynomial operator*(const Pterm& p);
  ExpandedPolynomial& operator*=(const ExpandedPolynomial& sp);
  ExpandedPolynomial operator*(const ExpandedPolynomial& sp);
 ExpandedPolynomial& operator/=(const ExpandedPolynomial& sp);
  ExpandedPolynomial operator/(const ExpandedPolynomial& sp);
  ExpandedPolynomial& operator^=(int);
  ExpandedPolynomial operator^(int);

  void setVariable(vector<string>);
  void print() const;
  void getSupport(int**) const;
  void getCoefficient(complex<double>* ) const;
 

  int size() const { return m_terms.size(); }
  const vector<string>& variable() const  { return m_variables;}

 private:
  VPterm m_terms;
  vector<string> m_variables;
};

#endif //INCLUDE_EXPANDEDPOLYNOMIAL_H

