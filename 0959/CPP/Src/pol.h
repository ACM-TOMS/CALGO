
#ifndef NTL_POL__H
#define NTL_POL__H

#include <NTL/vec_GF2.h>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

NTL_CLIENT

class pol {
public:

   string  _pol__rep;
   vec_GF2 _vec__rep;	

   pol() { _pol__rep = ""; clear(_vec__rep); }

   pol(const pol& a) { _pol__rep = ""; clear(_vec__rep); *this = a; }

   ~pol() { }

   pol& operator=(const pol& a) { _pol__rep = a._pol__rep; _vec__rep = a._vec__rep; return *this; }

   pol& operator=(string a) { _pol__rep = a; return *this; }

   pol& operator=(vec_GF2 a) { _vec__rep = a; return *this; }

};

string long2string(const long& number)
{
   ostringstream oss;
   oss << number;
   return oss.str();
}

string to_pol(const vec_GF2& a)
{
   vec_GF2 v;
   long spacen = a.length();
   long n = logtwo(spacen);
   long i, j, last;
   string s;
   v.SetLength(n);
         
   last = 0;
   i = spacen-1;
   while (last == 0)
   {
      if (a[i] == 1) last = i;
      i--;
   }
      	    
   for (i = 0; i < last; i++)
   {
      if (a[i] == 1)
      {
         v = to_vecGF2(i, n);
         if (IsZero(v))
         {
            s += "1";
         }
         else
         {
      	    for (j = 0; j < n; j++)
      	    {
      	       if (v[j] == 1)
      	       {
                  s += "x";
		  s += long2string(j+1);
               }
            }
         }
         s += "+";
      }
   }
         
   v = to_vecGF2(last, n);
   for (i = 0; i < n; i++)
   {
      if (v[i] == 1)
      {
         s += "x";
	 s += long2string(i+1);
      }
   }

   return s;
}


long len_vec(string& str)
{
   vector<string> monoms;
   unsigned long i, j, n = 0;

   Tokenize(str, monoms,"+");

   for (i = 0; i < monoms.size(); i++)
   {
      vector<string> indexes;
      Tokenize(monoms[i], indexes, "x");

      for (j = 0; j < indexes.size(); j++)
      {
	 istringstream Sindex(indexes[j]);
         unsigned long index;

	 Sindex >> index;
	 if (index > n) n = index;
      }
   }
   
   return n;
}


vec_GF2 to_vec(string& str, long n)
{
   vector<string> monoms;
   unsigned long i, j, spacen, coef;
   vec_GF2 v, anf;

   Tokenize(str, monoms,"+");
   v.SetLength(n);
   spacen = (1 << n);
   anf.SetLength(spacen);

   for (i = 0; i < monoms.size(); i++)
   {
      vector<string> indexes;

      if (monoms[i] == "1")
      {
	 anf[0] = 1;
      }
      else
      {
         Tokenize(monoms[i], indexes, "x");
         clear(v);
         for (j = 0; j < indexes.size(); j++)
         {
	    istringstream Sindex(indexes[j]);
            long index;

	    Sindex >> index;
	    v(index) = 1;
         }
         coef = conv_long(v);
         anf[coef] = 1;
      }
   }
   
   return anf;
}
  
  
string rep(pol& a)
{ 

   if (a._pol__rep == "") {
      a._pol__rep = to_pol(a._vec__rep);
   }

   return a._pol__rep;
}


vec_GF2 vec(pol& a)
{ 
   long n;

   if (a._vec__rep.length() <= 0) {
      n = len_vec(a._pol__rep);
      a._vec__rep = to_vec(a._pol__rep, n);
   }

   return a._vec__rep;
}


NTL_SNS ostream& operator<<(NTL_SNS ostream& s, pol a)
{
   string str;
   
   str = rep(a);
   s << str;

   return s;
}


NTL_SNS istream& operator>>(NTL_SNS istream& s, pol& x)
{
   string str;
 
   s >> str;
   x = str;
   
   return s;
}

#endif

