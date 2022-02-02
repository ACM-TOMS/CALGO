#ifndef ERROR_H
#define ERROR_H

#include <string>

//exception struct
//===============       
struct Error{
  std::string fError;
  Error( std::string s ): fError( s ) {};
};


//exception struct
//===============       
struct Warning{
  std::string fWarning;
  Warning( std::string s ): fWarning( s ) {};
};


#endif

