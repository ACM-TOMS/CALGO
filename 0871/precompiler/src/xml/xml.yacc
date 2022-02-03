  /*==================================================================*/
  /* description  : XMLParser                                         */
  /* author       : Walter Schreppers                                 */
  /* created      : 9/1/2002                                          */
  /* modified     : /                                                 */
  /* bugs         : /                                                 */
  /*==================================================================*/  


%name XMLParser


  /*======================= INCLUDES ===========================*/

%header{
  #include <stdio.h>
  #include <string> 
  #include <list>
  #include <iostream>
  #include <fstream>
  #include <FlexLexer.h>
  #include "xmlnode.h"
  #include "textnode.h"
  #include "elementnode.h"
  #include "pinode.h"
  #include "commentnode.h"
    
%}




  /*========================== TYPES ===========================*/
  
  %union {
    string*   str;
  }




  /*========================== TOKENS ==========================*/

  %token PROLOG_START VERSION EQ NR PROLOG_END
  %token DOCTYPE COMMENT_BEGIN COMMENT_END COMMENT_DATA
  %token OPEN_TAG_BEGIN OPEN_TAG_END ATTRIBUTE
  %token SINGLE_TAG_END
  %token PROCESS_TAG_TARGET PROCESS_TAG_DATA PROCESS_TAG_END
  %token CLOSE_TAG CLOSE_TAG_END        
  %token COMMENT TEXTNODE WHITESPACE



  /*===================== TYPE ASSOCIATIONS ====================*/
  
  %type<str> ATTRIBUTE
  %type<str> OPEN_TAG_BEGIN
  %type<str> CLOSE_TAG 
  %type<str> TEXTNODE
  %type<str> WHITESPACE
  %type<str> PROCESS_TAG_TARGET
  %type<str> PROCESS_TAG_DATA
  %type<str> process_data
  %type<str> NR
  %type<str> COMMENT_DATA
  %type<str> comment_text

  /*=========================== ORDER ==========================*/

  /* %left '*' '/' '\\' ARRAY_MULT ARRAY_DIV LEFT_ARRAY_DIV */
  %left WHITESPACE


  

/*=============Parser class members=============*/
%define MEMBERS \
  virtual ~XMLParser();                       \
  void init();                                \
  int parse(ifstream* in);                    \
  int parse();                                \
  int errorLine();                            \
  int errorColumn();                          \
  int syntaxError(const string&);             \
  int semanticError(const string&);           \
  void warning(const string&);                \
  string charpToStr(char*);                   \
  string stripWhiteAndQuotes(const string&);  \
  string stripWhites(const string&);          \
  void checkVersion(const string&);           \
  void openTag();                             \
  void singleTag();                           \
  void closeTag(const string&);               \
  void showTagStack();                        \
  void newTextNode(const string&);            \
  void newWhiteTextNode(const string&);       \
  void setIncludeWhites(bool);                \
  yyFlexLexer xmlLexer;                       \
  bool bSemError;                             \
  bool bSynError;                             \
  bool bSingleTag;                            \
  bool bIncludeWhites;                        \
  list<string> tagStack;                      \
  bool emptyAllowed;                          \
  string tagName;                             \
  XMLNode *tree,*treePos;                     \
  vector<Attribute> attributes;        

  
  

/*==============Lexer body call=================*/

%define LEX_BODY { return xmlLexer.yylex(); }




/*===== The error reporting funtion body =======*/

  /*  %define ERROR xmlError  */
%define ERROR_BODY {  } 



/*========= Constructor params + code ==========*/

%define CONSTRUCTOR_PARAM
%define CONSTRUCTOR_CODE { \
  init(); \
  bIncludeWhites=1; \
}



%start xml_doc

%%

xml_doc
  : internal_parts                  {}
  | prolog internal_parts           {}
  | doctype internal_parts          {}
  | prolog doctype internal_parts   {}
  ;



doctype
  : DOCTYPE               {}
  | WHITESPACE doctype    {}

prolog
  :PROLOG_START VERSION EQ NR PROLOG_END  {
                                            checkVersion(*$4);
                                            delete $4; //NR is string*
                                          }

  |PROLOG_START VERSION EQ NR process_data PROLOG_END  {
                                                  checkVersion(*$4);
                                                  //warning("Unused prolog data:"+*$5);
                                                  delete $4;
                                                  delete $5;
                                                }
  | WHITESPACE prolog                           {}
  ;                                          


internal_parts
  : part                  {}
  | part internal_parts   {}
  ;
  
part
  : open_tag                  {
                                openTag();
                              }
  | comment                   {}
  | CLOSE_TAG CLOSE_TAG_END   {
                                closeTag(*$1);
                                delete($1);
                              }
  | TEXTNODE                  {
                                newTextNode(*$1);
                                delete $1;
                              }
  | process_tag               {}
  | WHITESPACE                {
                                if(bIncludeWhites) newWhiteTextNode(*$1);
                                delete $1;
                              }
  ;


open_tag
  : OPEN_TAG_BEGIN open_tag_end   {
                                    tagName=string(*$1);
                                    delete $1;
                                    
                                  }
  ;
  

open_tag_end
  : OPEN_TAG_END              {
                                bSingleTag=false;
                              }
  | SINGLE_TAG_END            {
                                bSingleTag=true;
                              }
  | ATTRIBUTE open_tag_end    {
                                Attribute a(*$1);
                                attributes.push_back(a);
                                delete $1; 
                              }
  ;


process_tag 
  : PROCESS_TAG_TARGET process_data PROCESS_TAG_END     {
                                                          PINode* pNode=new PINode(*$1,*$2);
                                                          if(treePos==0){
                                                            semanticError("You need to define a root tag before using a processing instruction");
                                                          }
                                                          else{
                                                            treePos->appendChild(pNode);
                                                          }
                                                          delete $1;
                                                          delete $2;
                                                        }
  ;

process_data
  : PROCESS_TAG_DATA {
                        $$=new string(*$1);
                        delete $1;
                     }
  | process_data PROCESS_TAG_DATA  {
                                    $$=new string(*$1+*$2);
                                    delete $1;
                                    delete $2;
                                   }


comment
  : COMMENT_BEGIN comment_text COMMENT_END  {
                                              CommentNode* cNode=new CommentNode(*$2);
                                              if(treePos==0){
                                                semanticError("no comments allowed outside of a root");
                                              }
                                              else{
                                                treePos->appendChild(cNode);
                                              }
                                              delete $2;
                                            }


comment_text
  : COMMENT_DATA                  {
                                    $$=new string(*$1);
                                    delete $1;
                                  }
  | comment_text COMMENT_DATA     {
                                    $$=new string(*$1+*$2);
                                    delete $1;
                                    delete $2;
                                  }
%%


#include "XMLParseFunctions.h"
