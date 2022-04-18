//ExpandedPolynomial.cpp
//Author: Xing Li
//Date: 03/02/2000

#include "ExpandedPolynomial.h"
#include "PolynomialException.h"

#include <iostream>
#include <algorithm>
#include <cstring>

//Defintion of ExpandedPolynomial::Pterm
ExpandedPolynomial::Pterm::Pterm(ValueType v ): 
m_coefficient(v), m_degrees()
{
}

ExpandedPolynomial::Pterm::Pterm(ValueType v, const VarDegrees& d): 
m_coefficient(v), m_degrees(d)
{
}
      
ExpandedPolynomial::Pterm::Pterm(const Pterm& pt): 
m_coefficient(pt.m_coefficient), m_degrees(pt.m_degrees) 
{
}

ExpandedPolynomial::Pterm::~Pterm()
{
}

ExpandedPolynomial::Pterm& 
ExpandedPolynomial::Pterm::operator=(const Pterm& pt)
{
  if(this != &pt)
  {
    m_coefficient = pt.m_coefficient;
    m_degrees = pt.m_degrees;
  }
  return *this;
}

bool ExpandedPolynomial::Pterm::operator<(const Pterm& p) const
{
  VarDegrees::const_iterator it1 = m_degrees.begin();
  VarDegrees::const_iterator it2 = p.m_degrees.begin();
  for(;it1 != m_degrees.end() && it2 != p.m_degrees.end() 
	&& it1->first == it2->first && it1->second == it2->second ; it1++, it2++)
    ;
  
  return ( (it1 != m_degrees.end() && it2 != p.m_degrees.end() 
       && (it1->first<it2->first || (it1->first == it2->first && it1->second < it2->second) ))
      || (it1 == m_degrees.end() && it2 != p.m_degrees.end() ) );
}

bool ExpandedPolynomial::Pterm::operator>(const Pterm& p) const
{
  VarDegrees::const_iterator it1 = m_degrees.begin();
  VarDegrees::const_iterator it2 = p.m_degrees.begin();
  for(;it1 != m_degrees.end() && it2 != p.m_degrees.end() 
	&& it1->first == it2->first && it1->second == it2->second ; it1++, it2++)
    ;
  return ( (it1 != m_degrees.end() && it2 != p.m_degrees.end() 
	    && (it1->first>it2->first || (it1->first==it2->first && it1->second >it2->second) ))
	   ||(it1 != m_degrees.end() && it2 == p.m_degrees.end() ) );
}

bool ExpandedPolynomial::Pterm::operator==(const Pterm& p) const
{
  if ( m_degrees.size() != p.m_degrees.size() )
    return false;
  
  VarDegrees::const_iterator it1 = m_degrees.begin();
  VarDegrees::const_iterator it2 = p.m_degrees.begin();
  for(  ;it1 != m_degrees.end() && it2 != p.m_degrees.end() 
	  && it1->first == it2->first && it1->second == it2->second; it1++, it2++)
    ;
  
  return (it1 == m_degrees.end() && it2 == p.m_degrees.end() );
}

bool ExpandedPolynomial::Pterm::operator!=(const Pterm& p) const
{
  return !operator==(p);
}


ExpandedPolynomial::Pterm& ExpandedPolynomial::Pterm::operator*=(const Pterm& pt) 
{
  m_coefficient *= pt.m_coefficient;
  for(VarDegrees::const_iterator it2 = pt.m_degrees.begin();
      it2 != pt.m_degrees.end(); it2++)
  {
    VarDegrees::iterator it1 = m_degrees.begin();
    for( ;it1 != m_degrees.end() && it1->first < it2->first; ++it1)
      ;
    if(it1 == m_degrees.end() || it1->first > it2->first )
      m_degrees.insert(it1, *it2 );
    else
      it1->second += it2->second;
  }
  
  return *this;
}

ExpandedPolynomial::Pterm ExpandedPolynomial::Pterm::operator*(const Pterm& pt)
{
  return Pterm(*this) *= pt;
}

void ExpandedPolynomial::Pterm::print() const
{
  std::cout<<'+';
  std::cout<<m_coefficient;
  
  for(VarDegrees::const_iterator itr = m_degrees.begin();
      itr != m_degrees.end(); itr++)
    {
      std::cout<<'x'<<(*itr).first;
      if( (*itr).second > 1 )
	std::cout<<'^'<<(*itr).second;
    }
}




//Definition of ExpandedPolynomial
ExpandedPolynomial::ExpandedPolynomial():
m_terms(), m_variables()
{
}

ExpandedPolynomial::ExpandedPolynomial(const Pterm& p): 
m_terms(1,p), m_variables()
{ 
}

ExpandedPolynomial::~ExpandedPolynomial()
{
}

ExpandedPolynomial& ExpandedPolynomial::operator+=(const Pterm& p)
{
  VPterm::iterator itr = m_terms.begin();
  for( ;itr != m_terms.end() && *itr != p; ++itr )
  ;
  
  if( itr == m_terms.end() )
  m_terms.insert(itr, p);
  else
  {
    itr->m_coefficient += p.m_coefficient;
    if( abs( itr->m_coefficient ) < 1.0e-10 )
    m_terms.erase(itr);
  }
  
  return *this;
}

ExpandedPolynomial ExpandedPolynomial::operator+(const Pterm& p)
{
  return ExpandedPolynomial(*this) += p;
}

ExpandedPolynomial& ExpandedPolynomial::operator+=(const ExpandedPolynomial& sp)
{
  if( this != &sp)
  for(VPterm::const_iterator itr = sp.m_terms.begin();
      itr != sp.m_terms.end(); ++itr)
  operator+=(*itr);
  else
  for(VPterm::iterator itr = m_terms.begin(); itr != m_terms.end(); ++itr)
  itr->m_coefficient *= 2;
  
  return *this;
}

ExpandedPolynomial ExpandedPolynomial::operator+(const ExpandedPolynomial& sp)
{
  return ExpandedPolynomial(*this) += sp;
}    
    
ExpandedPolynomial& ExpandedPolynomial::operator-=(const Pterm& p)
{
  VPterm::iterator itr = m_terms.begin();
  for( ; itr != m_terms.end() && (*itr) != p; ++itr )
  ;
  if( itr == m_terms.end() )
  m_terms.insert(itr, p)->m_coefficient *= -1;
  else
  {
    itr->m_coefficient -= p.m_coefficient;
    if( abs( itr->m_coefficient ) < 1.0e-10 )
    m_terms.erase(itr);
  }
  
  return *this;
}

ExpandedPolynomial ExpandedPolynomial::operator-(const Pterm& p)
{
  return ExpandedPolynomial(*this) -= p;
}

ExpandedPolynomial& ExpandedPolynomial::operator-=(const ExpandedPolynomial& sp)
{
  if( this != & sp )
  for(VPterm::const_iterator itr = sp.m_terms.begin();
      itr != sp.m_terms.end(); ++itr)
  operator -=(*itr);
  else
  m_terms.clear();
  
  return *this;
}

ExpandedPolynomial ExpandedPolynomial::operator-(const ExpandedPolynomial& sp)
{
  return ExpandedPolynomial(*this) -= sp;
}    
    
ExpandedPolynomial& ExpandedPolynomial::operator*=(const Pterm& p)
{
  for(VPterm::iterator itr = m_terms.begin(); itr != m_terms.end(); ++itr )
  (*itr) *= p;
  return *this;
}

ExpandedPolynomial ExpandedPolynomial::operator*(const Pterm& p)
{
  return ExpandedPolynomial(*this) *= p;
}

ExpandedPolynomial& ExpandedPolynomial::operator*=(const ExpandedPolynomial& sp)
{
  if( this != &sp)
  {
    ExpandedPolynomial temp(*this);
    VPterm::const_iterator itr = sp.m_terms.begin();
    
    operator*=( *itr);
    ++itr;
    for( ; itr != sp.m_terms.end(); ++itr)
    operator+=( temp * (*itr) );
  }
  else
  {
    ExpandedPolynomial temp(*this);
    VPterm::const_iterator itr = temp.m_terms.begin();
    
    operator*=( *itr);
    ++itr;
    for( ; itr != temp.m_terms.end(); ++itr)
    operator+=( temp * (*itr) );
  }
  
  return *this;
}

ExpandedPolynomial ExpandedPolynomial::operator*(const ExpandedPolynomial& sp)
{
  return ExpandedPolynomial(*this) *= sp;
}

ExpandedPolynomial& ExpandedPolynomial::operator/=(const ExpandedPolynomial& sp)
{
    if(sp.m_terms.size() > 1 ||  sp.m_terms.begin()->m_degrees.size() > 0 )
        throw PolynomialException( PolynomialException::NON_CONST_DINOMINATOR);
    
    for(VPterm::iterator itr = m_terms.begin(); itr != m_terms.end(); ++itr)
        itr->m_coefficient /= sp.m_terms.begin()->m_coefficient;
  
  return *this;
}

ExpandedPolynomial ExpandedPolynomial::operator/(const ExpandedPolynomial& sp)
{
  return ExpandedPolynomial(*this) /= sp;
}
    
ExpandedPolynomial& ExpandedPolynomial::operator^=(int deg)
{
  if( deg < 0 )
      throw PolynomialException(PolynomialException::NEGATIVE_EXPONENT);
  
  if( deg == 0 )
  {
    m_terms.clear();
    m_terms.push_back(Pterm(1.0) );
  }
  else
  {
    ExpandedPolynomial temp(*this);
    for(int i =0; i < deg - 1; ++i)
    operator*=(temp);
  }
  return *this;
}

ExpandedPolynomial ExpandedPolynomial::operator^(int deg)
{
  return ExpandedPolynomial(*this)^=deg;
}

void ExpandedPolynomial::setVariable(vector<string> newv)
{
  if(m_variables.empty() )
    m_variables = newv;

//    copy(m_variables.begin(), m_variables.end(), ostream_iterator<string>(cout, " ") );
  
  sort(newv.begin(), newv.end() );
  for(vector<string>::iterator it1 = m_variables.begin(); it1 != m_variables.end(); ++it1)
    { 
      vector<string>::iterator it2 = newv.begin();
      for( ; it2 != newv.end() && *it1 > *it2 ; ++it2 )
	;
      if( it2 == newv.end() || *it2 != *it1 )
	newv.insert(it2, *it1);
    }
  
  for(int i = m_variables.size() -1; i >= 0; --i)
    {
      int q = 0;
      for(vector<string>::iterator itv = newv.begin();  m_variables[i] != *itv; ++itv, ++q)
	;
      for(VPterm::iterator it = m_terms.begin(); it != m_terms.end(); ++it)
	{
	  VarDegrees::iterator itd = it->m_degrees.begin();
	  while( itd != it->m_degrees.end() && itd->first < i)
	    ++itd;
	  if(itd != it->m_degrees.end() && itd->first == i )
	    itd->first = q;
	}
//        cout<<'('<<i<<"->"<<q<<")   ";
    }
  m_variables = newv;
//    cout<<"________";
//    copy(m_variables.begin(), m_variables.end(), ostream_iterator<string>(cout, " ") );
//    cout<<endl;
} 


void ExpandedPolynomial::print() const
{
  for(VPterm::const_iterator itr = m_terms.begin(); itr != m_terms.end(); ++itr)
  {
    if( imag(itr->m_coefficient) > - 1.0e-15 && imag(itr->m_coefficient) < 1.0e-15 )
      {
	cout<< ( (real(itr->m_coefficient) >= 0)? " +" :"") ;
	cout<<  real(itr->m_coefficient );
      }
    else
      cout<< itr->m_coefficient;

    for(VarDegrees::const_iterator it = itr->m_degrees.begin(); it != itr->m_degrees.end(); ++it)
      {
	cout<<m_variables[it->first];
	if( it->second > 1 )
	  cout<<'^'<<it->second;
      }
  }
}

void ExpandedPolynomial::getSupport(int** sp) const
{
  int** p = sp;
  for(VPterm::const_iterator itr = m_terms.begin(); itr != m_terms.end(); ++itr, ++p)
    {
      memset( *p, 0, m_variables.size()*sizeof(int) );
      for(VarDegrees::const_iterator it = itr->m_degrees.begin(); 
	  it != itr->m_degrees.end(); ++it)
	(*p)[it->first] = it->second;
    }
}

void ExpandedPolynomial::getCoefficient(complex<double>* coef) const
{
  complex<double>* p = coef;
  for(VPterm::const_iterator itr = m_terms.begin(); itr != m_terms.end(); ++itr, ++p)
    *p = itr->m_coefficient; 
}


