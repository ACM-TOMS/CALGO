#ifndef MATCHER_H
#define MATCHER_H

#include <vector>
#include <string>

#include "Graph.h"


#include "MatchElem.h"
#include "ConvertElem.h"
#include "XMLParser.h"


typedef map<int,ConvertElem> cmap;
typedef pair<int,ConvertElem> cpair;


class Matcher{
  public:
    Matcher();
    Matcher(XMLNode*);
    ~Matcher();
    
    MatchElem sourceLookup(const string&);
    MatchElem assignLookup(const string&);
    ConvertElem convertLookup(const string&);
    vector<ConvertElem> convertLookup(const string& start, const string& dest);
  
    //public members (for easy access)
    vector<MatchElem> sources;
    vector<MatchElem> targets;
    vector<MatchElem> assigns;
    vector<ConvertElem> converts;
    cmap convertMap;
  
  protected:
  
  private:
    void init();
    bool addElement(ElementNode*, XMLNode::iterator);
    void addConvert(const ConvertElem&,int);
    void addEdges(const string&, vector<MatchElem>&, const string&,int);
    ConvertElem getConvert(ElementNode*);
    MatchElem findMatchElem(vector<MatchElem> v, const string& name);
    MatchElem getMatchElem(XMLNode::iterator);

    
    Graph graph;

};

#endif
