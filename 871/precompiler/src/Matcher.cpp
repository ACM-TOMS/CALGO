#include "Matcher.h"

void Matcher::init(){
  sources.clear(); //stores source tags
  assigns.clear(); //stores rhs tags
  targets.clear(); //stores target tags
  
  converts.clear();  //stores convert tags
  convertMap.clear();
}

Matcher::Matcher(){
  init();
}


Matcher::~Matcher(){
  init();
}


bool Matcher::addElement(ElementNode* parent, XMLNode::iterator i){
  if( (*i)->nodeType()==XMLNode::ElementNode ){
    ElementNode* content=(ElementNode*) *i;
    string name=parent->attribute("name","");
    string val="";
    matchtypes mType=mEmpty;
    if(content->tagName()=="keyword"){
      mType=mKeyword;
    }
    else if( content->tagName()=="token" ){
      mType=mToken;
    }
    if(content->size()==1){
      XMLNode::iterator j=content->begin();
      if( (*j)->nodeType()==XMLNode::TextNode ){
        TextNode* t=(TextNode*) *j;
        val=t->data();
      }
    }
     
    if( parent->tagName() == "source" ) 
      sources.push_back( MatchElem(name, val, mType) );
    else if( parent->tagName() == "target" )
      targets.push_back( MatchElem(name, val, mType) );
    else if( parent->tagName() == "rhs" )
      assigns.push_back( MatchElem(name, val, mType) );
    else{
      Warning( "Found unsupported tag, named:" + parent->tagName() );
      return false;
    }
    return true;
  }
  else{
    return false;
  }

}


void Matcher::addEdges(const string& name, vector<MatchElem>& v, const string& target, int key){
  for(vector<MatchElem>::iterator i= v.begin(); i!=v.end();++i){
    if( i->name == name ) {
      if( i->isKeyword() ){
        graph.addEdge( i->keyword, target, key );
      }
      else if( i->isToken() ){
        graph.addEdge( i->token, target, key );
      }
    }
  }
}



void Matcher::addConvert(const ConvertElem& cElem, int key ){
  
  converts.push_back( cElem );
  convertMap.insert( cpair( key, cElem ) );
  
  
  //add the nodes to the graph
  string targetNode;
  if(cElem.target.isKeyword()){
    targetNode=cElem.target.keyword;
  }
  else {
    throw( Error("Error: target element "+cElem.target.name+" is not a keyword!"));
  }
  
  //graph.addEdge( word or token in source, assign or target, word in target, key);
  
  switch(cElem.sourceType){
    case cSource:      addEdges(cElem.sourceName, sources, targetNode,key); break;
    case cTarget:      addEdges(cElem.sourceName, targets, targetNode,key); break;
    case cAssignment:  addEdges(cElem.sourceName, assigns, targetNode,key); break;
    default: cerr<<"Warning: invalid source type in convert"<<endl;
  }
  
}



MatchElem Matcher::findMatchElem(vector<MatchElem> v, const string& name){
  MatchElem mNotFound(mError);
  
  for(vector<MatchElem>::iterator i=v.begin(); i!=v.end(); ++i){
    if( i->name==name ){
        return *i;
    }
  }
  
  return mNotFound;
}


MatchElem Matcher::getMatchElem(XMLNode::iterator i){
  if( (*i)->nodeType()==XMLNode::ElementNode){
    ElementNode* e=(ElementNode*) *i;
    if( e->tagName() == "source"){
      return findMatchElem(sources,e->attribute("name","") );
    }
    else if( e->tagName()=="target" ){
      return findMatchElem(targets, e->attribute("name","") );
    }
    else if( e->tagName()=="rhs" ){
      return findMatchElem(assigns, e->attribute("name","") );
    }
    else{
      cerr<<"found undefined tag '"<<e->tagName()<<"' in convert rule."<<endl;
      return MatchElem();
    }
  }
  else{
    cerr<<"Warning convert element contains invalid node."<<endl;
      return MatchElem();
  }
}




ConvertElem Matcher::getConvert(ElementNode* e){
  ConvertElem conv;
  if( e->hasChildren() ){
    if(e->size()==2){
      XMLNode::iterator i=e->begin();
      conv.setSource(i);
      ++i;
      conv.target=getMatchElem(i);
      
      conv.name=e->attribute( "name", "" );
      conv.hasOperation=false;
      conv.mtype=mConvert;
    }
    else if(e->size()==3){
      XMLNode::iterator i=e->begin();
      conv.setSource(i);
      ++i;
      conv.target=getMatchElem(i);
      ++i;
      conv.setOperation(i);

      conv.name=e->attribute( "name", "" );
      conv.mtype=mConvert;
    }
    else{
      cerr<<"Warning: convert tag '"<<e->attribute("name","")<<"' has wrong number of conversion elements."<<endl;
      conv.mtype=mError;
    }
  }
  else{
    cerr<<"Warning: convert tag '"<<e->attribute("name","")<<"' has no conversion elements."<<endl;
    conv.mtype=mError;
  }
  return conv;
}


Matcher::Matcher(XMLNode* tree){
  init();
  int key=0; //unique key for converts in the convertmap

  for(XMLNode::iterator i=tree->begin();i!=tree->end();++i){
    if( (*i)->nodeType() == XMLNode::ElementNode ){
      ElementNode* e=(ElementNode*) *i; //it was an element so cast it
      if(
          ( e->tagName() == "source" ) || 
          ( e->tagName() == "target" ) ||
          ( e->tagName() == "rhs")
        ){
        if( e->hasChildren() ){
          if((e->tagName()== "target") && (e->size()>1)){
            cerr<<"Error : target element contains more than 1 child!"<<endl;
            //..throw error
            exit(1);
          }
          for(XMLNode::iterator i=e->begin(); i!=e->end(); ++i){
            if( !addElement(e,i) )
              cerr<<"Warning: invalid element in "<<e->tagName()<<" tag '"<<e->attribute("name","")<<"'."<<endl;
          }

        }
        else{
          cerr<<"Warning: "<<e->tagName()<<" tag '"<<e->attribute("name","")<<"' does not contain child elements."<<endl;
        }
      }//is source element
      else if( e->tagName() == "convert" ){
        key++;
        ConvertElem conv=getConvert( e );        
        if( !conv.isError() && !conv.isEmpty() ){
          addConvert(conv,key);
        }        
        else{
          cerr<<"Warning: convert tag '"<<e->attribute("name","")<<"' has no conversion elements."<<endl;
        }

      }
    }//is element
  }//for
  
}



//use dijkstra's shortest path algo to find suitable path,
//return vector containing ConvertElem's which follow this path.
vector<ConvertElem> Matcher::convertLookup(const string& start, const string& dest){
  //for each edge along path add the corresponding
  //ConvertElem in a vector
  vector<ConvertElem> retPath;
  retPath.clear();

  
  try{
  
    //find shortest path in graph
    graph.unweighted(start);
    //graph.printPath(dest);
    vector<int> path;
    graph.getPath(dest,path);  
  
    for(vector<int>::iterator i=path.begin();i!=path.end();++i){
      cmap::const_iterator itr = convertMap.find( *i );
      if( itr!=convertMap.end() ){
        ConvertElem conv=itr->second;
        retPath.push_back( conv );
      }
    }
  
  }
  catch(Error e){
    throw Warning("Matcher::convertLookup : could not find path from "+start+" to "+dest );
  }
  
  return retPath;
}


MatchElem Matcher::sourceLookup(const string& word){
  MatchElem mNotFound(mError);
  
  for(vector<MatchElem>::iterator i=sources.begin();i!=sources.end();++i){
    if( ( i->isKeyword() ) &&
        ( word==i->keyword )    ){
        return *i;
    }
  }
  
  return mNotFound;
}

MatchElem Matcher::assignLookup(const string& word){
  MatchElem mNotFound(mError);
  
  for(vector<MatchElem>::iterator i=assigns.begin();i!=assigns.end();++i){
    if( ( i->isKeyword() ) &&
        ( word==i->keyword )    ){
        return *i;
    }
  }
  
  return mNotFound;
}


ConvertElem Matcher::convertLookup(const string& word){
  ConvertElem cNotFound(mError);
  
  for(vector<ConvertElem>::iterator i=converts.begin();i!=converts.end();++i){
  #ifdef _DEBUG_
    cout << "Matcher::convertLookup i->sourceName=" <<i->sourceName<<" word="<<word<<endl;
  #endif
    if( i->sourceName == word ){
      return *i;
    }
  }
  
  return cNotFound;
}


