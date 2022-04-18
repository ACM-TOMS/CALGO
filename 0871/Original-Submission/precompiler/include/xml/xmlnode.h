#ifndef XMLNODE_H
#define XMLNODE_H

#include <string>
#include <list>
#include <iostream>

using namespace std;

//BUGS: prevSibling and nextSibling sometimes segfault and are not optimal, in short, don't use em! :)


class XMLNode:public list<XMLNode*> {

  public:
    enum NodeType { Unknown=0, ElementNode = 1, TextNode = 2, 
                    CommentNode=3, DocumentNode=4, PINode=5,
                    DocumentTypeNode=6 };
                    
                    //PINode is same as ProcessingInstructionNode of QT domdocument...
                    
                    
                    /* hey man I'm not that obsessed with xml yet , I'll do these when I need em :)
                    CDATASectionNode = 4, EntityReferenceNode = 5, 
                    EntityNode = 6, 
                    DocumentNode = 9, DocumentTypeNode = 10, 
                    DocumentFragmentNode = 11, NotationNode = 12, 
                    BaseNode = 21, CharacterDataNode = 22 };
                    */
  
    //constructors/destructor
    //=======================
    XMLNode();
    XMLNode(XMLNode*); //give parent
    virtual ~XMLNode();
  
    //copy constructor
    //================
    XMLNode( const XMLNode& );
  
    //public members
    //==============
    void setParent(XMLNode* p);
    XMLNode* getParent(){return parent;}
    

    virtual void show(int indent=0);
    void showTree(XMLNode*, int indent=0) const;
    
    
    void appendChild(XMLNode*);
    void appendSibling(XMLNode*);    
    
    
    XMLNode* firstChild();
    XMLNode* lastChild();
    
    XMLNode* nextSibling();   
    XMLNode* prevSibling();
    
    virtual NodeType nodeType(){
      return Unknown;
    }
    
    XMLNode::iterator lookup(); //gives location in parent list as iterator (used by prevSibling and nextSibling)
    
    void setName(const string&);
    string getName() const;  

    bool hasChildren(){ return size()!=0; }
  
  protected:
    XMLNode* parent;
    
        
  private:
    string name;
    void destroy(XMLNode*);
    NodeType fType;
};

#endif



