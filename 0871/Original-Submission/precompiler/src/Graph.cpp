#include "Graph.h"


void Graph::addEdge( const string & sourceName, const string & destName ){
    Vertex * v = getVertex( sourceName );
    Vertex * w = getVertex( destName );
    v->adj.push_back( w );
}

void Graph::addEdge( const string & sourceName, const string & destName, int convertId){
    Vertex * v = getVertex( sourceName );
    Vertex * w = getVertex( destName );
    v->adj.push_back( w );
    v->edges.push_back( convertId );
}


void Graph::printPath( const string & destName ) const {
    vmap::const_iterator itr = vertexMap.find( destName );

    if( itr == vertexMap.end( ) )
    {
      //cerr << "Destination vertex '" + destName + "' not found" << endl;
      return;
    }

    const Vertex & w = *(*itr).second;
    if( w.dist == G_INFINITY ){
      //cerr << destName << " is unreachable";
    }
    else{
      printPath( w );
      cout << endl;
    }
}

void Graph::getPath( const string& destName , vector<int>& cpath) const {
  cpath.clear();
  vmap::const_iterator itr = vertexMap.find( destName );

    if( itr == vertexMap.end( ) )
    {
        throw Error( "Graph::getPath : Destination vertex '" + destName + "' not found" );
    }

    const Vertex & w = *(*itr).second;
    if( w.dist == G_INFINITY ){
        throw Error( "Graph::getPath : "+ destName + " is unreachable" );
    }
    else{
        getPath( w , cpath);
    }
}

void Graph::getPath( const Vertex & dest, vector<int>& cpath ) const{
    if( dest.path != NULL )
    {
        getPath( *dest.path , cpath );
        cpath.push_back( dest.convertId );
    }
}


// If vertexName is not present, add it to vertexMap
// In either case, return the Vertex
Vertex * Graph::getVertex( const string & vertexName ) {
    vmap::iterator itr = vertexMap.find( vertexName );

    if( itr == vertexMap.end( ) )
    {
        Vertex *newv = new Vertex( vertexName );
        allVertices.push_back( newv );
        vertexMap.insert( vpair( vertexName, newv ) );
        return newv;
    }
    return (*itr).second;
}

void Graph::printPath( const Vertex & dest ) const{
    if( dest.path != NULL )
    {
        printPath( *dest.path );
        cout <<"-"<<dest.convertId<<"->";
    }
    cout << dest.name;
}

void Graph::clearAll( ){
    for( unsigned int i = 0; i < allVertices.size( ); i++ )
        allVertices[ i ]->reset( );
}

Graph::~Graph( ){
    for( unsigned int i = 0; i < allVertices.size( ); i++ )
        delete allVertices[ i ];
}


void Graph::unweighted( const string & startName ){
    clearAll( );

    vmap::iterator itr = vertexMap.find( startName );

    if( itr == vertexMap.end( ) )
    {
        throw Error( "Graph::unweighted : " + startName + " is not a vertex in this graph" );
        return;
    }

    Vertex *start = (*itr).second;
    list<Vertex *> q;
    q.push_back( start ); start->dist = 0;

    while( !q.empty( ) )
    {
        Vertex *v = q.front( ); q.pop_front( );

        for( unsigned int i = 0; i < v->adj.size( ); i++ )
        {
            Vertex *w = v->adj[ i ];
            if( w->dist == G_INFINITY )
            {
                w->dist = v->dist + 1; //every weight is 1
                w->path = v;
                if( i< v->edges.size() ) w->convertId = v->edges[ i ]; //set edge name in pathEdge
                q.push_back( w );
            }
        }
    }
}

