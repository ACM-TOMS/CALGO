/*=============================================================================
author        :Walter Schreppers
filename      :ConfigFile.cpp
created       :/
modified      :15/05/2000
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#include "ConfigFile.h"

/*-----------------------------------------------------------------------------
name        :ConfigFile
description :Constructor
parameters  :char*
return      :/
exceptions  :/
algorithm   :use loadFile to load a configuration file
-----------------------------------------------------------------------------*/
ConfigFile::ConfigFile(char* filename){
  loadFile(filename);
}


/*-----------------------------------------------------------------------------
name        :ConfigFile
description :Constructor
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
ConfigFile::ConfigFile(){
  this->clear();
}



/*-----------------------------------------------------------------------------
name        :loadFile
description :load a config file
parameters  :char*
return      :/
exceptions  :/
algorithm   :read a file with filename line per line and use addLine to insert
             it
-----------------------------------------------------------------------------*/
void ConfigFile::loadFile(char* filename){
  this->clear();
  
  /* === old style ===
  ifstream from(filename);
  if(from){                   //file exists if it doesn't exist nothing will be
                              //added to the vector with addLine and locate will
                              //always return false 
    string strLine;
    while(!from.eof()){
      getline(from,strLine);
      if(!from.eof()) addLine(strLine);    
    }
  }
  */
  
  ifstream from(filename);
  if(from){
    XMLParser xml;
    xml.setIncludeWhites(false); //skip white space textnodes 
    xml.parse(&from);
    readTree(xml.tree);
  }
  
}


/*-----------------------------------------------------------------------------
name        :readTree
description :get rules from conversion tree. Use addElem member to add the
             ConfigElem's to this vector for doing locate's
parameters  :XMLNode*
return      :/
exceptions  :/
algorithm   : traverse tree and get variable and function tags to be skipped
-----------------------------------------------------------------------------*/
void ConfigFile::readTree(XMLNode* tree){
  for(XMLNode::iterator i=tree->begin();i!=tree->end();++i){
    if( (*i)->nodeType() == XMLNode::ElementNode ){
      ElementNode* e=(ElementNode*) *i; //it was an element so cast it
      if( e->tagName() == "skip" ){
        if( e->hasChildren() ){
          addElement( e );
        }
        else{
          cerr<<"Warning: "<<e->tagName()<<" tag "<<e->attribute( "name", "" )<<" does not contain child elements."<<endl;
        }
      }//is skip
      else{
        cerr<<"Warning: "<<e->tagName()<<" is an invalid tagname in this context."<<endl;
      }
    }//is element
    else{
      cerr<<"Warning: invalid tag in skip configuration file"<<endl;
    }
  }//for

}



/*-----------------------------------------------------------------------------
name        :addElement
description :add a ConfigElem to the vector of this class
parameters  :ElementNode*
return      :/
exceptions  :/
algorithm   : check for tags variable and function. Read their name attribute
              and use this to add a ConfigElem.
              If there is a parse error then this rule is ignored
-----------------------------------------------------------------------------*/
void ConfigFile::addElement(ElementNode* e){

  string func="", var="";
  bool add=false;
  
  for(XMLNode::iterator i=e->begin(); i!=e->end(); ++i){
    if( (*i)->nodeType() == XMLNode::ElementNode ){
      ElementNode* skip=(ElementNode*) *i;
      if( skip->tagName() == "function" ){
        add=true;
        func = skip->attribute( "name", "/" );
      }
      else if( skip->tagName() == "variable" ){
        add=true;
        var = skip->attribute( "name", "/" );
      }
      else{
        cerr<<"Warning : skip tag "<< e->attribute( "name", "" ) <<" contains invalid element, rule ignored."<<endl;
        add=false;
      }
    }//is element
  }//for
  
  if( add ){
    ConfigElem cElem( func, var );
    this->push_back( cElem );        
  }

}

		
/*-----------------------------------------------------------------------------
name        :addLine
description :this is depricated, used for old skip file format
parameters  :char*
return      :/
exceptions  :/
algorithm   : read the inputfile line per line and add elements if the
                line contains '::' the string left of '::' is a procedure name
                the string to the right is a variable name.
              
              if the line doesn't contain '::' skip the line
              
              also skip lines with a '#' (can be used for comments in 
                the configfile)
-----------------------------------------------------------------------------*/
void ConfigFile::addLine(string s){
    int pos=s.find("::");
    //strip comments
    if(s.find("#")!=string::npos){ //check if line contains '#' 
      int posComment=s.find("#");
      s=s.substr(0,posComment); //strip everything after '#'
    }
        
    if(s.find("::")!=string::npos){  //check if line contains '::' and add ConfigElem if true
      string proc=s.substr(0,pos);
      string  var=s.substr(pos+2,s.length()-pos-2);
      ConfigElem e(proc,var);
      this->push_back(e);
    }
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
string ConfigFile::stripSpaces(const string& s){
  string strRet="";
  for(unsigned int i=0;i<s.length();i++){
    if((s[i]!=' ')&&(s[i]!='(')&&(s[i]!=')')
      &&(s[i]!='*')&&(s[i]!='&')&&(s[i]!='\t')){
      strRet+=s[i];
    }
  }
  return strRet;  
}
		

/*-----------------------------------------------------------------------------
name        :locate
description :
parameters  :const string& proc,const string var
return      :/
exceptions  :/
algorithm   : 
              for c++ class style member functions we strip classname::
              ex. ConfigFile::stripSpaces becomes stripSpaces
              same for variable names ...
              
              we return true in the following cases
              - if in the configfile contained procedurename::
                lookup will return true for every variable
                with proc==procedurename
              
              - if the configfile contained ::variablename
                lookup will return true for every variable which
                name is variablename regardless of the procedurename proc
                
              - if the configfile contained procedurename::variablename
                lookup returns true only if proc==procedurename and
                var==variable name
              
              In other cases false is returned.
              We will use this function in pre.lex to skip convertion of
              certain variables.
-----------------------------------------------------------------------------*/
bool ConfigFile::locate(const string& proc,const string& var){
  bool found=0;
  
  /*
    we strip spaces and brakkets from the strings proc and var
    because pre.lex sometimes adds leading or trailing spaces
    this has no influence on locate algorithm
    since spaces/brakkets are not allowed in variable names of C/C++
  */
  string strProc=stripSpaces(proc);
  if(strProc.find("::")!=string::npos) strProc.erase(0,strProc.find("::")+2); //erase classname::
  
  string strVar=stripSpaces(var);
  if(strVar.find("::")!=string::npos) strVar.erase(0,strVar.find("::")+2);
    
  if(this->size()>0){                 //if vector empty locate must return false
    vector<ConfigElem>::iterator p=this->begin();
    while((p!=this->end())&&(found==0)){
      if(p->getProcedure().length()==0){
        if(p->getVar()==strVar){
          found=1;
        }
      }
      else if(p->getVar().length()==0){
        if(p->getProcedure()==strProc){
          found=1;
        }
      }
      else if((p->getProcedure()==strProc)&&(p->getVar()==strVar)){
        found=1;
      }
      ++p;
    }
  }
  return found;
}



/*
 // some test code


int main(){
  ConfigFile skip("test.conf");
  cout<<"constructed ConfigFile skip\n";
  if(skip.locate("blah","varB")){
    cout<<"found varB"<<endl;
  }
  else{
    cout<<"not found"<<endl;
  }
  
  cout<<endl<<endl<<"printing the entire file"<<endl;
  for(ConfigFile::iterator p=skip.begin();p!=skip.end();++p){
    cout<<"proc="<<p->getProcedure()<<" length="<<p->getProcedure().length()<<endl;
    cout<<"var ="<<p->getVar()<<      " length="<<p->getVar().length()<<endl;
  }
}    
*/


