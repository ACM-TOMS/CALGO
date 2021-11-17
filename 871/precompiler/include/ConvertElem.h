#ifndef CONVERT_ELEM_H
#define CONVERT_ELEM_H

#include "MatchElem.h"
#include "XMLParser.h"
#include "error.h"

#include <vector>
using namespace std;

enum stypes {
              cSource,
              cTarget,
              cAssignment,
              cError
            };

class ConvertElem{
  public:
    //constructors/destructor
    //=======================
    ConvertElem();
    ConvertElem(matchtypes m){ mtype=m; }
    ConvertElem(const ConvertElem& m);
    ~ConvertElem(){}

    //operators
    //=========
    ConvertElem& operator=(const ConvertElem& m);
  
    //members
    //=======

    void setOperation(XMLNode::iterator);
    void setSource(XMLNode::iterator);

    bool isError(){ return mtype==mError; }
    bool isEmpty(){ return mtype==mEmpty; }


    matchtypes mtype;
    string operation,name;
    string operationArg;
    string sourceName; //name of source tag
    stypes sourceType;  //can be source,target or assignment
    MatchElem target;
    bool hasOperation;
    

};


#endif

