#ifndef MATCH_ELEM_H
#define MATCH_ELEM_H

#include <string>
using namespace std;


enum matchtypes {mKeyword,mToken,mConvert,mEmpty,mError};

class MatchElem{

  public:
    //constructor/destructor
    //======================
    
    MatchElem();
    MatchElem(matchtypes);
    MatchElem(const string&);
    MatchElem(const string&,const string&);
    MatchElem(const string&,const string&,matchtypes);
    MatchElem(const MatchElem&);    
    ~MatchElem(){}
    
    
    //members
    //=======
    void setKeyword(const string&);

    bool isKeyword() const {return mtype==mKeyword;}
    bool isToken() const {return mtype==mToken;}
    bool isEmpty() const {return mtype==mEmpty;}
    bool isError() const {return mtype==mError;}
        
    //operator
    //========
    MatchElem& operator=(const MatchElem&);
    
    string name;
    string keyword;
    string token;
    matchtypes mtype;
  
  private:
    void init();
};

#endif

