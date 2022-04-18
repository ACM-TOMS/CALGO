  /*==================================================================*/
  /* description  : XMLParser                                         */
  /* author       : Walter Schreppers                                 */
  /* created      : 9/1/2002                                          */
  /* modified     : 11/2/2002                                         */
  /* bugs         : added extra check in stripWhiteSpaceAndQuotes to  */
  /*                resolve whitespace bug                            */
  /*==================================================================*/  


#ifndef XMLPARSE_FUNCTIONS_H
#define XMLPARSE_FUNCTIONS_H

/*============================================================================
                           PARSE FUNCTIONS 
============================================================================*/

extern void setXMLParser(XMLParser*);
extern int getXMLLineNr();
extern int getXMLColumn();
extern void resetXMLLineCol();

void XMLParser::init(){
  tagStack.clear();   /* stack for matching open/close tags */  
  emptyAllowed=true;  /* document can have only 1 root => tagStack can be empty only ones*/  
  bSemError=0;        /* set to 1 when there is a semantic error*/  
  bSynError=0;        /* set to 1 when there is a syntax error*/
  tree=0; 
  treePos=0; 
  attributes.clear(); 
  resetXMLLineCol();
}


void XMLParser::setIncludeWhites(bool b){
  bIncludeWhites=b;
}

int XMLParser::parse(ifstream* in){
  init();  
  //set new input stream for lexer
  setXMLParser(this);
  xmlLexer.switch_streams(in,0);
  bool ok=1-yyparse();
  if(tagStack.size()!=0){
    semanticError("Tag '"+*tagStack.rbegin()+"' is still open");
  }

  return ( 
           ok&&
           !bSynError && 
           !bSemError   
         );
}


int XMLParser::parse(){
  init();
  //this defaults to cin as input stream.
  setXMLParser(this);
  xmlLexer.switch_streams(0,0);
  bool ok=1-yyparse();
  if(tagStack.size()!=0){
    semanticError("Tag '"+*tagStack.rbegin()+"' is still open");
  }
  return ( 
           ok &&
           !bSynError &&
           !bSemError   
         );
}




int XMLParser::errorLine(){
  return getXMLLineNr();
}

int XMLParser::errorColumn(){
  return getXMLColumn();
}



XMLParser::~XMLParser(){
  //if(tree!=0) delete tree;
  //if(treePos!=0) delete treePos;
}

/*============================================================================
                           CONVERSION FUNCTIONS 
============================================================================*/

/*-----------------------------------------------------------------------------
name        :charpToStr
description :convert char* to string
parameters  :char*
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string XMLParser::charpToStr(char *c){
	string s="";
	while (*c!='\0'){
		s+=*c;
		c++;
	}
	return s;
}









/*============================================================================
                            ERRORS AND WARNINGS
============================================================================*/

//is defined in bison.cc
//void yyerror(char* msg){
//  cout<<"Syntax error : "<<msg<<" at line "<<getXMLLineNr()<<" , col"<<getXMLColumn()<<endl;
//}


int XMLParser::syntaxError( const string& msg ) {
  cerr<<"SYNTAX ERROR at line "   <<getXMLLineNr()<<", col "<<getXMLColumn()<<" : "<< msg <<endl;
  //exit(-1);
  bSynError=1;
  YYACCEPT;
};


int XMLParser::semanticError( const string& s ){
  cerr<<"SEMANTIC ERROR at line " <<getXMLLineNr()<<", col "<<getXMLColumn()<<" : "<< s <<endl;
  bSemError=1;
  return 0;
}

void XMLParser::warning( const string& s ) {
  cerr<<"WARNING at line "        <<getXMLLineNr()<<", col "<<getXMLColumn()<<" : "<< s <<endl;
};




/*============================================================================
                      DESTINATION FILE CONSTRUCTION
============================================================================*/

//...



/*============================================================================
                          SEMANTIC CHECKING 
============================================================================*/

void XMLParser::checkVersion(const string& version){
  string v=version.substr(1,version.length()-2);
  if(v!="1.0"){
    warning("Wrong XML version: "+v+", should be 1.0");  
  }
}



void XMLParser::newTextNode(const string& s){
  TextNode* textNode=new TextNode();
  if(bIncludeWhites){
    textNode->setData(s);
  }
  else{
    textNode->setData(stripWhiteAndQuotes(s));
  }
  if( treePos==0 ) {
    semanticError("You need to define a root tag before using a text node");
  }
  else if( tagStack.size()==0 ){
    semanticError("Text node found outside document root");
  }
  else{
    treePos->appendChild(textNode);
  }
}

void XMLParser::newWhiteTextNode(const string& s){
  TextNode* textNode=new TextNode();
  textNode->setDataLiteral(s); //no data replacements => faster

  //whitespace inside document becomes a textnode
  //outside document root it is ignored
  if( (treePos!=0) && (tagStack.size()!=0) ) {
    treePos->appendChild(textNode);
  }
}


string XMLParser::stripWhiteAndQuotes(const string& orig){
  if(orig.size()==0) return "";
  string result=orig;

  string::size_type a=0;
  while( ( (result[a]==' ') ||
          (result[a]=='\t') ||
          (result[a]=='\'') ||
          (result[a]=='"')  ||
          (result[a]=='\n') )&&
          (a<result.length()) ) a++;

  string::size_type b=result.size()-1;
  while( ( (result[b]==' ') || 
          (result[b]=='\t') ||
          (result[b]=='\'') ||
          (result[b]=='"')  ||
          (result[b]=='\n') )&&
          (b>0) ) b--;

  return result.substr(a,b-a+1);
}

void XMLParser::openTag(){
  ElementNode* newNode=new ElementNode(tagName);
  for(unsigned int i=0;i<attributes.size();i++){
    Attribute a( attributes[i].name, stripWhiteAndQuotes(attributes[i].value) );
    newNode->attributes.push_back(a);
  }
  attributes.clear();

  if(tagStack.size()==0){ //stack is empty
    if(!emptyAllowed){
      semanticError("XML document may only have 1 root, you may not open the tag '"+tagName+"' here!");
    }
    else{
      emptyAllowed=false;
      tree=newNode;
      treePos=tree;
    }
  }
  else{ //stack is not empty
    treePos->appendChild(newNode);
  }

  if(!bSingleTag){
    treePos=newNode; //go deeper into tree
    tagStack.push_back(tagName);
  }
}



void XMLParser::showTagStack(){
  cout<<"  tagStack contents:"<<endl;
  for(list<string>::iterator i=tagStack.begin();i!=tagStack.end();i++){
    cout<<"    '"<<*i<<"'"<<endl;
  }
}

void XMLParser::closeTag(const string& closeName){
  if(tagStack.size()!=0){
    list<string>::iterator s=tagStack.end(); //last element
    s--;
    if(*s!=closeName){
      semanticError("The closing tag '"+closeName+"' does not match the opening tag '"+*s+"'.");
      //showTagStack();
    }
    else{
      tagStack.pop_back();
      treePos=treePos->getParent();
    }
  }
  else{
    semanticError("The closing tag '"+closeName+"' does not have a matching open tag.");
  } 
}


#endif


