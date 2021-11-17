
#ifndef NTL_vec_pol__H
#define NTL_vec_pol__H

#include "pol.h"
#include <NTL/mat_GF2.h>

NTL_CLIENT

class vec_pol {  
public:
   vector<string> _pol__rep;  
   mat_GF2 _mat__rep;
    
   vec_pol() { _pol__rep.clear(); clear(_mat__rep); }  
   vec_pol(const vec_pol& a) { _pol__rep.clear(); clear(_mat__rep); *this = a; }     
   ~vec_pol() { }  
    
   long length() const { return _pol__rep.size(); }  
   
   string& operator[](long l__i)   
   {  
      return _pol__rep[l__i];  
   }  
  
   const string& operator[](long l__i) const 
   {   
      return _pol__rep[l__i];  
   }  
    
   string& operator()(long l__i) { return (*this)[l__i-1]; }  
   const string& operator()(long l__i) const { return (*this)[l__i-1]; } 
   
   vec_pol& operator=(const vec_pol& a) { _pol__rep = a._pol__rep; _mat__rep = a._mat__rep; return *this; }
         
   vec_pol& operator=(vector<string> a) { _pol__rep = a; return *this; }

   vec_pol& operator=(mat_GF2 a) { _mat__rep = a; return *this; }

};  

vector<string> to_vec_pol(const mat_GF2& A)
{
   long m = A.NumCols(); 
   long i;
   mat_GF2 ATr;
   vector<string> polynoms;

   ATr = transpose(A);

   for (i = 0; i < m; i++)
   {      	
      string str;
   
      str = to_pol(ATr[i]);
      polynoms.push_back(str);
   }

   return polynoms;
}


mat_GF2 to_mat(vector<string>& str)
{
   long i, temp, spacen, n = 0, m = str.size();
   vec_GF2 v, w;
   mat_GF2 A, ATr;

   for (i = 0; i < m; i++)
   {  
      temp = len_vec(str[i]);
      if (temp > n) n = temp;
   }  

   spacen = (1 << n);
   ATr.SetDims(m, spacen);
   w.SetLength(spacen);

   for (i = 0; i < m; i++)
   {  
      clear(w);
      w = to_vec(str[i],n);
      ATr[i] = w;
   }  
   A = transpose(ATr);

   return A;
}


vector<string> rep(vec_pol& a)
{

   if (a._pol__rep.size() <= 0) {
      a._pol__rep = to_vec_pol(a._mat__rep);
   }

   return a._pol__rep; 
}


mat_GF2 mat(vec_pol& a)
{

   if (a._mat__rep.NumRows() <= 0) {
      a._mat__rep = to_mat(a._pol__rep);
   }

   return a._mat__rep; 
}


NTL_SNS istream& operator>>(NTL_SNS istream& s, vec_pol& a)
{
   vector<string> v;
   string str;

   while (!s.eof()) {
      getline(s,str);
      if (str != "")
      {
         v.push_back(str);
      }       
   }
   
   a = v;
   
   return s;
} 
  
NTL_SNS ostream& operator<<(NTL_SNS ostream& s, vec_pol& x)
{
   vector<string> str;
   unsigned long i;

   str = rep(x);

   for (i = 0; i < str.size(); i++)
   {  
      s << str[i] << endl;
   }

   return s;

} 

#endif
