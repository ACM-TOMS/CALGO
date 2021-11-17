#ifndef ATTRIBUTE_H
#define ATTRIBUTE_H

#include <string>

using namespace std;

class Attribute{
  public:
    Attribute(const string&);
    Attribute(const string& n, const string& v);
    Attribute(const Attribute&);

    Attribute& operator=(const Attribute&);
        
    string name;
    string value;
};

#endif



