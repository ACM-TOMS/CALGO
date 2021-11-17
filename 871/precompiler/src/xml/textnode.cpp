#include "textnode.h"


TextNode::TextNode(){
  fData="";
}


TextNode::TextNode(const string& s){
  fData=s;
}


TextNode::TextNode(const TextNode& a):XMLNode(a){
  if( this != &a ){
    fData=a.fData;
  }
}



TextNode& TextNode::operator=(const TextNode& a){
  fData=a.fData;
  return *this;
}


void TextNode::replaceAll(string& s, const string& findStr,const string& repStr){
  string::size_type pos=s.find(findStr);
  while(pos!=string::npos){
    s.replace(pos,findStr.length(),repStr);
    pos= s.find(findStr);
  }
}


void TextNode::setData(const string& s){
  fData=s;

  //replace &lt; and &gt; , ...

  replaceAll( fData, "&lt;",   "<"  );
  replaceAll( fData, "&gt;",   ">"  );
  replaceAll( fData, "&amp;",  "&"  );


  replaceAll( fData, "&quot;", "\"" );
  replaceAll( fData, "&quot;", "\'" );
  replaceAll( fData, "&nbsp;", " "  );
  replaceAll( fData, "&copy;", "©"  );
  replaceAll( fData, "&deg;",  "°"  );

  // ... maybe some more in the future ...

}

string TextNode::data(){
  return fData;
}


XMLNode::NodeType TextNode::nodeType(){
    return XMLNode::TextNode;
}

void TextNode::show(int indent){
  string s="";
  for(int i=0;i<indent;i++) s+="  ";
  cout<<s+fData<<endl;
}
