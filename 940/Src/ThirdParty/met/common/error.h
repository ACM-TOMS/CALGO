//-*- C++ -*-

#ifndef ERROR_H_
#define ERROR_H_

#include <iostream>
#include <iomanip>

// global error function
static void gerror(const char* msg="") {
    std::cerr << msg << std::endl;
  exit(1);
}

#endif // ERROR_H_
