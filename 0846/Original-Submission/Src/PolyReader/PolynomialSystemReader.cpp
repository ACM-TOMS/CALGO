//PolynomialSystemReader.cpp
//Author: Xing Li
//Date: 02/24/2000


#include <algorithm>
#include <numeric>
#include <fstream>
#include "PolynomialSystemReader.h"
#include "GenericPolynomial.h"

#ifdef _MSC_VER
#define for if(0){} else for
#endif

PolynomialSystemReader::PolynomialSystemReader(istream& is)
{
  char c;
  while( (c = is.get() ) != EOF && ( c == ' ' || c == '\n' ) )
    ;
  
  if( c == EOF || c != '{' )
    {
      cout<<"expect a { as the begining of the polynomial system.\n" ;
      exit(0);
    }
  
  vector<string> vs;
  for(;;)
    {
      string str;
      while( (c = is.get() ) != EOF && c != ';'  && c != '}')
	if( c != '\n' )
	  str.append(&c, &c +1);
      
      if( c == EOF )
	{
	  cout<<"The polynomial system is expected to end with }\n";
	  exit(-1);
	}
      
      int count = 0;
      for(string::iterator it = str.begin(); it != str.end(); ++it)
	if( *it != ' ' )
	  count++;
      
      if(count < 3 )
        {
	  if( c == '}' && count == 0 )
	    break;
	  
	  
	  cout<<str<< c <<" unexpected.\n";
	  exit(1);
        }
      
      //       cout<<str<<endl;
      m_poly.push_back( GenericPolynomial( str ).expand() );
      
      if( c == '}'  )
	break;
      
    }  
  //     for( vector<ExpandedPolynomial>::const_iterator it = m_poly.begin(); it != m_poly.end(); ++it)
  //       it->print(),cout<<endl;
  
  vector<ExpandedPolynomial>::iterator  beg = m_poly.begin();
  for( vector<ExpandedPolynomial>::const_iterator it = beg + 1; it != m_poly.end(); ++it)
    beg->setVariable( it->variable() );
  
  for( vector<ExpandedPolynomial>::iterator it = beg + 1; it != m_poly.end(); ++it)
    it->setVariable( beg->variable() );
  
  m_variable = beg->variable();
  //     for(vector<string>::const_iterator it = m_variable.begin(); it != m_variable.end(); ++it)
  //       cout<< *it <<"  ";
  //     cout<<endl;
}

void PolynomialSystemReader::print()
{
  for(vector<ExpandedPolynomial>::const_iterator it = m_poly.begin(); it != m_poly.end(); ++it)
    it->print(), cout<<';'<<endl;
}

void PolynomialSystemReader::getSupportSize(int *ss) const 
{
  for(unsigned int i = 0; i < m_poly.size(); ++i)
    ss[i] = m_poly[i].size();
}

void PolynomialSystemReader::getCoefficient(complex<double> *coef) const 
{
  complex<double> * p = coef;
  for(vector<ExpandedPolynomial>::const_iterator it = m_poly.begin(); 
      it != m_poly.end(); ++it)
    {
      it->getCoefficient( p );
      p += it->size();
    }
}


void PolynomialSystemReader::getSupports(int** sp) const 
{
  int** pp = sp;
  for(vector<ExpandedPolynomial>::const_iterator it = m_poly.begin(); it != m_poly.end(); ++it)
    {
      it->getSupport( pp );
      pp += it->size();
    }
}

/*
void PolynomialSystemReader::createSupportCoefficientFiles(const char* dsupp, const char* dcoef) const
{
  typedef int* Pint;
  typedef Pint* PPint;
  typedef complex<double> Dcomplex;
  typedef Dcomplex* PDcomplex;
  typedef vector<int> Vint;
  typedef vector<Vint> VVint;
  typedef vector<double> Vdouble;

  int dim = m_variable.size();
  Pint sptSize = new int[dim];
  getSupportSize( sptSize );
  
  int totalPoints = std::accumulate(sptSize, sptSize + dim, 0);
  Pint sptMem = new int[ dim *(totalPoints + dim) ];
  PDcomplex coef = new complex<double>[totalPoints  ];
  PPint spt = new int*[totalPoints+dim];
  spt[0] = sptMem;
  for(int i = 1; i < totalPoints + dim; ++i)
      spt[i] = spt[i-1] + dim;
 
  getCoefficient( coef );
  getSupports( spt );

  std::ofstream dataCoef(dcoef), dataSupp(dsupp);
  dataSupp << dim << "  " << dim <<endl;
  for(int i = 0; i < dim; ++i) {
          dataSupp << 1 << "  ";
  }
  dataSupp << endl;
  for(int i = 0, n = 1; i < dim; ++i) {
          dataSupp << n << "   " ;
          dataSupp << ( n+=sptSize[i]-1 ) <<endl; ++n;
  }

  for(int i = 0; i < totalPoints; ++i) {
	dataCoef << real( coef[i] )<< endl;
	for(int j = 0; j < dim; ++j)
		dataSupp << spt[i][j] <<"  ";

	dataSupp <<endl;
}  */

void PolynomialSystemReader::GetDimensions(int &nVar,int &nPts) const
{
  // added by Tangan Gao

  typedef int* Pint;
  typedef Pint* PPint;

  nVar = m_variable.size();
  Pint sptSize = new int[nVar];
  getSupportSize( sptSize );
  
  nPts = std::accumulate(sptSize, sptSize + nVar, 0);

  delete [] sptSize;
}

void PolynomialSystemReader::GetSupport(int* SptSize, int** Spt) const
{
  // added by Tangan Gao

  typedef int* Pint;
  typedef Pint* PPint;

  int dim = m_variable.size();
  Pint sptSize = new int[dim];
  getSupportSize( sptSize );
  
  int totalPoints = std::accumulate(sptSize, sptSize + dim, 0);
  Pint sptMem = new int[ dim *(totalPoints + dim) ];
  PPint spt = new int*[totalPoints+dim];
  spt[0] = sptMem;
  for(int i = 1; i < totalPoints + dim; ++i)
      spt[i] = spt[i-1] + dim;
 
  getSupports( spt );

  for (int i=0; i<dim; i++) SptSize[i] = sptSize[i];

  for (int i=0; i<totalPoints; i++)
  {
      for (int j=0; j<dim; j++) Spt[i][j] =  spt[i][j];
  }

  delete []sptMem;
  delete []spt;
  delete []sptSize;
}

