#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <stack>
#include <list>
using namespace std;

#include <limits.h>

#include "error.h"

#define G_INFINITY INT_MAX
#define UNKNOWN -1

struct Vertex
{
    string           name;      // Vertex name
    vector<Vertex *> adj;       // Adjacent vertices
    vector<int>      edges;     // convert id's on edges
    int              dist;      // Cost
    Vertex          *path;      // Previous vertex on shortest path
    int              convertId; // convert id of path

    Vertex( const string & nm ) : name( nm )
      { reset( ); }

    void reset( )
      { dist = G_INFINITY; path = NULL; convertId=UNKNOWN;}
};

typedef map<string,Vertex *> vmap;
typedef pair<string,Vertex *> vpair;


class Graph {
  public:
    Graph( ) { }
    ~Graph( );
    
    void addEdge( const string & sourceName, const string & destName ); //no convertId ...
    void addEdge( const string & sourceName, const string & destName, int convertId); //with convertId on edges

    void unweighted( const string & startName );
    void getPath( const string& , vector<int>&) const;

    void printPath( const string & destName ) const; //show path on cout

  private:
    Vertex * getVertex( const string & vertexName );
    void printPath( const Vertex & dest ) const;
    void getPath( const Vertex & dest, vector<int>& ) const;
    void clearAll( );

    vmap vertexMap;
    vector<Vertex *> allVertices;
    
};


#endif
