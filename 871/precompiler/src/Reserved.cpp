/*=============================================================================
author        :Walter Schreppers
filename      :Reserved.cpp
created       :/
modified      :18/04/2000
version       :2
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#include "Reserved.h"




/*-----------------------------------------------------------------------------
name        :Reserved
description :constructor
parameters  :/
return      :Reserved
exceptions  :/
algorithm   :add all C++ reserved words to the private vector<string>resWords 
-----------------------------------------------------------------------------*/
Reserved::Reserved(){
	push_back("and");
  push_back("and_eq");
  push_back("asm");
  push_back("auto");
  push_back("bitand");
  push_back("bitor");
  push_back("break");
  push_back("case");
  push_back("catch");
  push_back("class");
  push_back("const");
  push_back("const_cast");
  push_back("continue");
  push_back("default");
  push_back("delete");
  push_back("do");
  push_back("dynamic_cast");
  push_back("else");
  push_back("enum");
  push_back("explicit");
  push_back("export");
  push_back("extern");
  push_back("false");
  push_back("for");
  push_back("friend");
  push_back("goto");
  push_back("if");
  push_back("inline");
  push_back("mutable");
  push_back("namespace");
  push_back("new");
  push_back("not");
  push_back("not_eq");
  push_back("operator");
  push_back("or");
  push_back("or_eq");
  push_back("private");
  push_back("protected");
  push_back("public");
  push_back("register");
  push_back("reinterpret_cast");
  push_back("return");
  push_back("short");
  push_back("signed");
  push_back("sizeof");
  push_back("static");
  push_back("static_cast");
  push_back("struct");
  push_back("switch");
  push_back("template");
  push_back("this");
  push_back("throw");
  push_back("true");
  push_back("try");
  push_back("typedef");
  push_back("typeid");
  push_back("typename");
  push_back("union");
  push_back("unsigned");
  push_back("using");
  push_back("virtual");
  //push_back("void"); regard void as a basic type like char
  push_back("volatile");
  push_back("while");
  push_back("xor");
  push_back("xor_eq");
}

/*-----------------------------------------------------------------------------
name        :Reserved
description :constructor
parameters  :char *
return      :Reserved
exceptions  :/
algorithm   :open a file containing the reserved words.
             this function is not needed when using the default constructor
-----------------------------------------------------------------------------*/
Reserved::Reserved(char* fileName){
	ifstream input(fileName);
	string str;
	while(!input.eof()){
		input>>str;
		if(!input.eof()) push_back(str);
	}
}


/*-----------------------------------------------------------------------------
name        :printAll
description :
parameters  :char *
return      :void
exceptions  :/
algorithm   :this functio was made to easily implement the default constructor
             it prints all the reserved words with code around it so that it may
             be inserted into the default constructor
-----------------------------------------------------------------------------*/
void Reserved::printAll(){
  for(vector<string>::iterator p=begin();p<end();p++){
    cout<<"  push_back(\""<<*p<<"\");"<<endl;
  }
}


/*-----------------------------------------------------------------------------
name        :isReserved
description :return true if the given char* is a reserved word
parameters  :char *
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool Reserved::isReserved(char* cptr){
	vector<string>::iterator p=find(begin(),end(),stripSpaces(string(cptr)));
	return p!=end();
}
		
		
		

/*-----------------------------------------------------------------------------
name        :isReserved
description :return true if the given string is a reserved word
parameters  :string
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool Reserved::isReserved(string s){
	vector<string>::iterator p=find(begin(),end(),stripSpaces(s));
	return p!=end();
}



/*-----------------------------------------------------------------------------
name        :stripSpaces
description :
description :strip characters of the string so that only the name is left over
parameters  :const string& s
return      :/
exceptions  :/
algorithm   :strip the string from spaces,tabs , ')', '(' , '*', '&' 
-----------------------------------------------------------------------------*/
string Reserved::stripSpaces(const string& s){
  string strRet="";
  for( unsigned int i=0; i < s.length(); i++ ){
    if( (s[i]!=' ') &&
        (s[i]!='(') &&
        (s[i]!=')') &&
        (s[i]!='*') &&
        (s[i]!='&') &&
        (s[i]!='\t') )
    {
      strRet+=s[i];
    }
  }
  return strRet;  
}


		



/*
  this main function was used to create the default constructor by streaming its output to a file
*/

/*
  int main(){
	  Reserved resresWords("my-reserved");
    resprintAll();
  }
*/
