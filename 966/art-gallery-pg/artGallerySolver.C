/*****************************************************************************
 * This code is part of Art Gallery Solver (AGSol) Package, which aims the 
 * resolution of the Art Gallery Problem With Point Guards.
 *
 * This software version (1.0.2) has been tested under and is compatible with 
 * CGAL 3.9 and GLPK 4.52.
 *
 * Authors:
 *  Davi Colli Tozoni - davi.tozoni@gmail.com
 *  Marcelo Castilho Couto - coutomarcelo@gmail.com
 *
 * AGSol Concept and Design: 
 *  Davi Colli Tozoni, Marcelo Castilho Couto, Pedro Jussieu de Rezende & Cid 
 * Carvalho de Souza.
 * 
 * Other information can be found at:
 *  http://www.ic.unicamp.br/~cid/Problem-instances/Art-Gallery/index.html
 *
 * --
 * 
 *  This program is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation, either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.
 *      
 ****************************************************************************/


#include "artGallerySolver.h"

/**
 * Loads a polygon from file.
 */
PolygonWithHolesExt ArtGallerySolver::loadPol(char* name) {
    PolygonWithHolesExt pol;

    std::ifstream file;
    file.open(name);
    PolygonWithHoles yada;
    CGAL::set_ascii_mode(file);
    file >> yada;
    pol = PolygonWithHolesExt(yada);
    file.close();

    std::string polName(name);
    _polName = polName;

    return pol;
}

/**
 * Prints the polygon in .pol (CGAL) format.
 */
void ArtGallerySolver::printPol(PolygonWithHolesExt pol) {
    cout << "[ArtGallerySolver] " << pol << std::endl;
}

/**
 * Selects the initial witnesses based on the discretization technique chosen.
 */
vector<Point> ArtGallerySolver::selectInitialDisPoints(){
    double t0 = CPUTIME(ruse);
    double t1;
    vector<Point> disPoints;

    switch(_mode) {
        case ALL_VERTICES:
            {
                for(Vertex_iterator vit_tmp = _polygon.outer_boundary().vertices_begin(); vit_tmp != _polygon.outer_boundary().vertices_end(); ++vit_tmp){
                    disPoints.push_back(*vit_tmp);
                    _witnessTest.insert(pair<Point,bool>(*vit_tmp,true));
                }

                for(PolygonWithHoles::Hole_const_iterator holeIt = _polygon.holes_begin(); holeIt != _polygon.holes_end(); ++holeIt) {
                    for(Vertex_iterator vit_tmp = (*holeIt).vertices_begin(); vit_tmp != (*holeIt).vertices_end(); ++vit_tmp){
                        disPoints.push_back(*vit_tmp);
                        _witnessTest.insert(pair<Point,bool>(*vit_tmp,true));
                    }
                } 
                break;
            }
        case CONVEX_VERTICES:
            {
                /* Identifies all convex vertices from the boundary */
                Polygon::Vertex_circulator vit_minus1 = _polygon.outer_boundary().vertices_circulator();
                Polygon::Vertex_circulator vit = _polygon.outer_boundary().vertices_circulator();
                Polygon::Vertex_circulator vit_plus1 = _polygon.outer_boundary().vertices_circulator();
                for(;(*vit_plus1) != (*vit); vit_plus1++);
                for(;(*vit_minus1) != (*vit); vit_minus1--);
                vit_plus1++;
                vit_minus1--;

                for(int i=0; i<_polygon.outer_boundary().size(); i++) {
                    if (CGAL::left_turn(*vit_minus1, *vit, *vit_plus1)) {
                        disPoints.push_back(*vit);
                        _witnessTest.insert(pair<Point,bool>(*vit,true));
                    }
                    vit_plus1++;
                    vit_minus1++;
                    vit++;
                }

                /* Identifies all reflex vertices from holes */
                for(PolygonWithHoles::Hole_const_iterator holeIt = _polygon.holes_begin(); holeIt != _polygon.holes_end(); ++holeIt) {
                    vit_minus1 = (*holeIt).vertices_circulator();
                    vit = (*holeIt).vertices_circulator();
                    vit_plus1 = (*holeIt).vertices_circulator();

                    for(;(*vit_plus1) != (*vit); vit_plus1++);
                    for(;(*vit_minus1) != (*vit); vit_minus1--);
                    vit_plus1++;
                    vit_minus1--;

                    for(int i=0; i<(*holeIt).size(); i++) {
                        if (CGAL::left_turn(*vit_minus1, *vit, *vit_plus1)) {
                            disPoints.push_back(*vit);
                            _witnessTest.insert(pair<Point,bool>(*vit,true));
                        }
                        vit_plus1++;
                        vit_minus1++;
                        vit++;
                    }
                } 
                break;
            }
        case CHWA_POINTS:
            {
                Point a, b;
                bool aIsConvex = false;
                bool bIsConvex = false;

                /* Identifies all convex vertices */
                Polygon::Vertex_circulator vit_minus1 = _polygon.outer_boundary().vertices_circulator();
                Polygon::Vertex_circulator vit = _polygon.outer_boundary().vertices_circulator();
                Polygon::Vertex_circulator vit_plus1 = _polygon.outer_boundary().vertices_circulator();
                for(;(*vit_plus1) != (*vit); vit_plus1++);
                for(;(*vit_minus1) != (*vit); vit_minus1--);
                vit_plus1++;
                vit_minus1--;

                /* Analyzes first vertex */
                a = (*vit);
                if (CGAL::left_turn(*vit_minus1, *vit, *vit_plus1)) {
                    aIsConvex = true;
                }
                vit_plus1++;
                vit_minus1++;
                vit++;

                for(int i=0; i<_polygon.outer_boundary().size(); i++) {
                    /* Analyzes new vertex */
                    b = (*vit);    
                    if (CGAL::left_turn(*vit_minus1, *vit, *vit_plus1)) {
                        bIsConvex = true;
                    }
                    else{
                        bIsConvex = false;
                    }

                    if (bIsConvex && !aIsConvex){
                        if(_witnessTest.find(b) == _witnessTest.end()){
                            disPoints.push_back(b);
                            _witnessTest.insert(pair<Point,bool>(b,true));
                        }
                    }
                    else if (aIsConvex && !bIsConvex){
                        if(_witnessTest.find(a) == _witnessTest.end()){
                            disPoints.push_back(a);
                            _witnessTest.insert(pair<Point,bool>(a,true));
                        }
                    }
                    else if (!aIsConvex && !bIsConvex){
                        Vector medium = (b - a)/2;
                        Point p_medium = Point(a.x()+medium.x(), a.y()+medium.y());
                        disPoints.push_back(p_medium);
                        _witnessTest.insert(pair<Point,bool>(p_medium,true));
                    }

                    /* Holds information about current vertex */
                    a = b;
                    aIsConvex = bIsConvex;

                    vit_plus1++;
                    vit_minus1++;
                    vit++;
                }

                /* Holes ... */
                for(PolygonWithHoles::Hole_const_iterator holeIt = _polygon.holes_begin(); holeIt != _polygon.holes_end(); ++holeIt) {
                    aIsConvex = false;
                    bIsConvex = false;

                    vit_minus1 = (*holeIt).vertices_circulator();
                    vit = (*holeIt).vertices_circulator();
                    vit_plus1 = (*holeIt).vertices_circulator();

                    for(;(*vit_plus1) != (*vit); vit_plus1++);
                    for(;(*vit_minus1) != (*vit); vit_minus1--);
                    vit_plus1++;
                    vit_minus1--;

                    /* Analyzes first vertex */
                    a = (*vit);
                    if (CGAL::left_turn(*vit_minus1, *vit, *vit_plus1)) {
                        aIsConvex = true;
                    }
                    vit_plus1++;
                    vit_minus1++;
                    vit++;

                    for(int i=0; i<(*holeIt).size(); i++) {
                        /* Analyzes new vertex */
                        b = (*vit);    
                        if (CGAL::left_turn(*vit_minus1, *vit, *vit_plus1)) {
                            bIsConvex = true;
                        }
                        else{
                            bIsConvex = false;
                        }

                        if (bIsConvex && !aIsConvex){
                            if(_witnessTest.find(b) == _witnessTest.end()){
                                disPoints.push_back(b);
                                _witnessTest.insert(pair<Point,bool>(b,true));
                            }
                        }
                        else if (aIsConvex && !bIsConvex){
                            if(_witnessTest.find(a) == _witnessTest.end()){
                                disPoints.push_back(a);
                                _witnessTest.insert(pair<Point,bool>(a,true));
                            }
                        }
                        else if (!aIsConvex && !bIsConvex){
                            Vector medium = (b - a)/2;
                            Point p_medium = Point(a.x()+medium.x(), a.y()+medium.y());
                            disPoints.push_back(p_medium);
                            _witnessTest.insert(pair<Point,bool>(p_medium,true));
                        }

                        /* Holds information about the current vertex */
                        a = b;
                        aIsConvex = bIsConvex;

                        vit_plus1++;
                        vit_minus1++;
                        vit++;
                    }

                }
                break;
            }
        case CHWA_POINTS_EXTENDED:
            {
                Point a, b;
                bool aIsConvex = false;
                bool bIsConvex = false;

                /* Identifies all convex vertices */
                Polygon::Vertex_circulator vit_minus1 = _polygon.outer_boundary().vertices_circulator();
                Polygon::Vertex_circulator vit = _polygon.outer_boundary().vertices_circulator();
                Polygon::Vertex_circulator vit_plus1 = _polygon.outer_boundary().vertices_circulator();
                for(;(*vit_plus1) != (*vit); vit_plus1++);
                for(;(*vit_minus1) != (*vit); vit_minus1--);
                vit_plus1++;
                vit_minus1--;

                /* Analyzes first vertex */
                a = (*vit);
                if (CGAL::left_turn(*vit_minus1, *vit, *vit_plus1)) {
                    aIsConvex = true;
                }
                vit_plus1++;
                vit_minus1++;
                vit++;

                for(int i=0; i<_polygon.outer_boundary().size(); i++) {
                    /* Analyzes new vertex */
                    b = (*vit);    
                    if (CGAL::left_turn(*vit_minus1, *vit, *vit_plus1)) {
                        bIsConvex = true;
                    }
                    else{
                        bIsConvex = false;
                    }

                    if ((bIsConvex && !aIsConvex) || (aIsConvex && !bIsConvex)){
                        if(_witnessTest.find(a) == _witnessTest.end()){
                            disPoints.push_back(a);
                            _witnessTest.insert(pair<Point,bool>(a,true));
                        }
                        if(_witnessTest.find(b) == _witnessTest.end()){
                            disPoints.push_back(b);
                            _witnessTest.insert(pair<Point,bool>(b,true));
                        }
                    }
                    else if (!aIsConvex && !bIsConvex){
                        Vector medium = (b - a)/2;
                        Point p_medium = Point(a.x()+medium.x(), a.y()+medium.y());
                        disPoints.push_back(p_medium);
                        _witnessTest.insert(pair<Point,bool>(p_medium,true));
                    }

                    /* Holds information about current vertex */
                    a = b;
                    aIsConvex = bIsConvex;

                    vit_plus1++;
                    vit_minus1++;
                    vit++;
                }

                /* Holes... */
                for(PolygonWithHoles::Hole_const_iterator holeIt = _polygon.holes_begin(); holeIt != _polygon.holes_end(); ++holeIt) {
                    aIsConvex = false;
                    bIsConvex = false;

                    vit_minus1 = (*holeIt).vertices_circulator();
                    vit = (*holeIt).vertices_circulator();
                    vit_plus1 = (*holeIt).vertices_circulator();

                    for(;(*vit_plus1) != (*vit); vit_plus1++);
                    for(;(*vit_minus1) != (*vit); vit_minus1--);
                    vit_plus1++;
                    vit_minus1--;

                    /* Analizes first vertex */
                    a = (*vit);
                    if (CGAL::left_turn(*vit_minus1, *vit, *vit_plus1)) {
                        aIsConvex = true;
                    }
                    vit_plus1++;
                    vit_minus1++;
                    vit++;

                    for(int i=0; i<(*holeIt).size(); i++) {
                        /* Analizes new vertex */
                        b = (*vit);    
                        if (CGAL::left_turn(*vit_minus1, *vit, *vit_plus1)) {
                            bIsConvex = true;
                        }
                        else{
                            bIsConvex = false;
                        }

                        if ((bIsConvex && !aIsConvex) || (aIsConvex && !bIsConvex)){
                            if(_witnessTest.find(a) == _witnessTest.end()){
                                disPoints.push_back(a);
                                _witnessTest.insert(pair<Point,bool>(a,true));
                            }
                            if(_witnessTest.find(b) == _witnessTest.end()){
                                disPoints.push_back(b);
                                _witnessTest.insert(pair<Point,bool>(b,true));
                            }
                        }
                        else if (!aIsConvex && !bIsConvex){
                            Vector medium = (b - a)/2;
                            Point p_medium = Point(a.x()+medium.x(), a.y()+medium.y());
                            disPoints.push_back(p_medium);
                            _witnessTest.insert(pair<Point,bool>(p_medium,true));
                        }

                        /* Holds information about current vertex */
                        a = b;
                        aIsConvex = bIsConvex;

                        vit_plus1++;
                        vit_minus1++;
                        vit++;
                    }

                }
                break;
            }
        default:
            break;
    } 

    t1 = CPUTIME(ruse);
    _initDisTime += t1 - t0;

    return disPoints;
}

/**
 * Adds each new witness and its visibility polygon in the Structure.
 */
void ArtGallerySolver::insertDisPointsToArrange(vector<Point> dis){
    double t0 = CPUTIME(ruse);
    double t1;
    PolygonExt visPol;
    tr1::unordered_map<Point, PolygonExt, PointHash>::iterator it;

    _visNewHiddenPoints.clear();

    for(int i=0; i<dis.size(); i++){
        it = _visCache.find(dis[i]);
        if(it == _visCache.end()){
            /* Computes visibility polygon */
            visPol = _polygon.getVisibility(dis[i]);

            _visCache.insert(pair<Point,PolygonExt>(dis[i],visPol));
        }
        else{
            visPol = it->second;
        }

        _hiddenPoints.push_back(dis[i]);
        _visHiddenPoints.push_back(visPol);
        _visNewHiddenPoints.push_back(visPol);
    }

    t1 = CPUTIME(ruse);
    _insertDisTime += t1 - t0;

    return;
}

/**
 * Selects guard candidates for the next iteration. The guard candidates are
 * all the vertices of the polygon and also the vertices of light AVPs obtained 
 * from the visibility arrangement.
 */
void ArtGallerySolver::selectGuardCandidates(){
    double t0 = CPUTIME(ruse);
    double t1;
    vector<Point> oldGuardCandidates;
    vector<PolygonExt> visOldGuardCandidates; 
    tr1::unordered_map<Point, int, PointHash> candidatesUsed; 
    tr1::unordered_map<Point, int, PointHash>::iterator itMap; 

    _groups.clear();

    oldGuardCandidates = _guardCandidates;
    visOldGuardCandidates = _visGuardCandidates;

    /* Cleans initial set */
    _guardCandidates.clear();
    _visGuardCandidates.clear();

    /* Adds all vertices */
    vector<int> newGroup;
    for(Vertex_iterator vIt = _polygon.outer_boundary().vertices_begin(); vIt != _polygon.outer_boundary().vertices_end(); ++vIt){
        _guardCandidates.push_back(*vIt);
        candidatesUsed.insert(pair<Point, int>(*vIt, _guardCandidates.size() - 1));    

        newGroup.push_back(_guardCandidates.size() - 1);
    } 
    _groups.push_back(newGroup);
    newGroup.clear();

    PolygonWithHoles::Hole_const_iterator holeIt = _polygon.holes_begin();
    for(; holeIt != _polygon.holes_end(); ++holeIt) {
        for(Vertex_iterator vIt = (*holeIt).vertices_begin(); vIt != (*holeIt).vertices_end(); ++vIt){
            _guardCandidates.push_back(*vIt);
            candidatesUsed.insert(pair<Point, int>(*vIt, _guardCandidates.size() - 1));    
            newGroup.push_back(_guardCandidates.size() - 1);
        }
        _groups.push_back(newGroup);
        newGroup.clear();
    } 
    
    /* Adds vertices of light AVPs */
    if(_mainLoopIterations == 1){
        _guardsGrid = new AVPLightGrid();
    }
    _guardsGrid->setPolygon(_polygon);
    _guardsGrid->setVisibilityPolygons(_visHiddenPoints);
    if(_mainLoopIterations == 1){
        _guardsGrid->makeGrid();
    }
    else{
        _guardsGrid->addToGrid(_visNewHiddenPoints);
    }
    std::vector<PolygonExt> lightAVPs = _guardsGrid->getStruct2();

    _arrangement = _guardsGrid->getStruct1();
    _lightAVPs = lightAVPs;

    int AVPsSize = lightAVPs.size();
    for (int j = 0; j < AVPsSize; j++) {
        PolygonExt guardPolygon = lightAVPs.at(j);
        for (Polygon::Vertex_iterator vIt = guardPolygon.vertices_begin(); vIt != guardPolygon.vertices_end(); vIt++) {
            itMap = candidatesUsed.find(*vIt);
            if (itMap == candidatesUsed.end()) {
                _guardCandidates.push_back(*vIt);
                candidatesUsed.insert(pair<Point, int>(*vIt, _guardCandidates.size() - 1));    
                newGroup.push_back(_guardCandidates.size() - 1);
            }
            else {
                newGroup.push_back(itMap->second);
            }
        }
        _groups.push_back(newGroup);
        newGroup.clear();
    }

    t1 = CPUTIME(ruse);
    _selectGuardCandTime += t1 - t0;

    return; 
}

/**
 * Based on the uncovered regions, selects new witnesses.
 */
vector<Point> ArtGallerySolver::selectNewDisPoints(){
    double t0 = CPUTIME(ruse);
    double t1;
    vector<Point> newDis;
    list<PolygonWithHoles> notCoveredList;
    _notCoveredFirst.polygons_with_holes(back_inserter(notCoveredList));
    list<PolygonWithHoles>::const_iterator it;

    cout << "[ArtGallerySolver] # of uncovered regions found: " << notCoveredList.size() << endl;

    PolygonWithHoles PH = _polygon;
    for (it = notCoveredList.begin(); it != notCoveredList.end(); ++it) {
        PolygonWithHoles QH = (*it);

        Polygon P = PH.outer_boundary();
        Polygon Q = QH.outer_boundary();
        for(Edge_iterator eit_P = P.edges_begin(); eit_P != P.edges_end(); ++eit_P){
            Point p_s = (*eit_P).source(); 
            Point p_t = (*eit_P).target();

            for(Edge_iterator eit_Q = Q.edges_begin(); eit_Q != Q.edges_end(); ++eit_Q){
                Point q_s = (*eit_Q).source();
                Point q_t = (*eit_Q).target();

                if((*eit_P).has_on(q_s) && (*eit_P).has_on(q_t)){
                    Vector medium = (q_t - q_s)/2;
                    Point p_medium = Point(q_s.x()+medium.x(), q_s.y()+medium.y());
                    if(_witnessTest.find(p_medium) == _witnessTest.end()){
                        newDis.push_back(p_medium);
                        _witnessTest.insert(pair<Point,bool>(p_medium,true));
                    }
                }
                else if((*eit_P).has_on(q_s)){
                    if(_witnessTest.find(q_s) == _witnessTest.end()){
                        newDis.push_back(q_s);
                        _witnessTest.insert(pair<Point,bool>(q_s,true));
                    }
                }
                else if((*eit_P).has_on(q_t)){
                    if(_witnessTest.find(q_t) == _witnessTest.end()){
                        newDis.push_back(q_t);
                        _witnessTest.insert(pair<Point,bool>(q_t,true));
                    }
                }
            }
        }

        for(PolygonWithHoles::Hole_iterator holeIt = PH.holes_begin(); holeIt != PH.holes_end(); ++holeIt) {

            P = (*holeIt);
            for(Edge_iterator eit_P = P.edges_begin(); eit_P != P.edges_end(); ++eit_P){
                Point p_s = (*eit_P).source(); 
                Point p_t = (*eit_P).target();

                Q = QH.outer_boundary();
                for(Edge_iterator eit_Q = Q.edges_begin(); eit_Q != Q.edges_end(); ++eit_Q){
                    Point q_s = (*eit_Q).source();
                    Point q_t = (*eit_Q).target();

                    if((*eit_P).has_on(q_s) && (*eit_P).has_on(q_t)){
                        Vector medium = (q_t - q_s)/2;
                        Point p_medium = Point(q_s.x()+medium.x(), q_s.y()+medium.y());
                        if(_witnessTest.find(p_medium) == _witnessTest.end()){
                            newDis.push_back(p_medium);
                            _witnessTest.insert(pair<Point,bool>(p_medium,true));
                        }
                    }
                    else if((*eit_P).has_on(q_s)){
                        if(_witnessTest.find(q_s) == _witnessTest.end()){
                            newDis.push_back(q_s);
                            _witnessTest.insert(pair<Point,bool>(q_s,true));
                        }
                    }
                    else if((*eit_P).has_on(q_t)){
                        if(_witnessTest.find(q_t) == _witnessTest.end()){
                            newDis.push_back(q_t);
                            _witnessTest.insert(pair<Point,bool>(q_t,true));
                        }
                    }
                }

                for(PolygonWithHoles::Hole_iterator qHoleIt = QH.holes_begin(); qHoleIt != QH.holes_end(); ++qHoleIt) {
                    Q = (*qHoleIt);
                    for(Edge_iterator eit_Q = Q.edges_begin(); eit_Q != Q.edges_end(); ++eit_Q){
                        Point q_s = (*eit_Q).source();
                        Point q_t = (*eit_Q).target();

                        if((*eit_P).has_on(q_s) && (*eit_P).has_on(q_t)){
                            Vector medium = (q_t - q_s)/2;
                            Point p_medium = Point(q_s.x()+medium.x(), q_s.y()+medium.y());
                            if(_witnessTest.find(p_medium) == _witnessTest.end()){
                                newDis.push_back(p_medium);
                                _witnessTest.insert(pair<Point,bool>(p_medium,true));
                            }
                        }
                        else if((*eit_P).has_on(q_s)){
                            if(_witnessTest.find(q_s) == _witnessTest.end()){
                                newDis.push_back(q_s);
                                _witnessTest.insert(pair<Point,bool>(q_s,true));
                            }
                        }
                        else if((*eit_P).has_on(q_t)){
                            if(_witnessTest.find(q_t) == _witnessTest.end()){
                                newDis.push_back(q_t);
                                _witnessTest.insert(pair<Point,bool>(q_t,true));
                            }
                        }
                    }
                }
            }
        }

        if (QH.number_of_holes() == 0) {
            Point aux = internalPoint(QH.outer_boundary());
            Point newPoint = findSimplerInteriorPoint(QH.outer_boundary(), aux);
            if(_witnessTest.find(newPoint) == _witnessTest.end()){
                newDis.push_back(newPoint);
                _witnessTest.insert(pair<Point,bool>(newPoint,true));
            }
        }
        else {
            Point aux = internalPoint(QH);
            Point newPoint = findSimplerInteriorPoint(QH, aux);
            if(_witnessTest.find(newPoint) == _witnessTest.end()){
                newDis.push_back(newPoint);
                _witnessTest.insert(pair<Point,bool>(newPoint,true));
            }
        }

    }

    t1 = CPUTIME(ruse);
    _selectNewDisTime += t1 - t0;

    return newDis;
}

/**
 * Sets the discretization technique to be used.
 */
void ArtGallerySolver::setMode(char * mode){
    _mode = stringModeToInt(mode);
    _sMode = mode;
}

/**
 * Converts mode from string to integer.
 */
int ArtGallerySolver::stringModeToInt(char * mode){
    if(!strcmp(mode, "ALL_VERTICES"))
        return ALL_VERTICES;
    else if(!strcmp(mode, "CONVEX_VERTICES"))
        return CONVEX_VERTICES;
    else if(!strcmp(mode, "CHWA_POINTS"))
        return CHWA_POINTS;
    else if(!strcmp(mode, "CHWA_POINTS_EXTENDED"))
        return CHWA_POINTS_EXTENDED;
    else 
        return -1;
}

/**
 * Gets cardinality of the current solution.
 */
int ArtGallerySolver::getCardinality(){
    return _guards.size();
}

/**
 * Adds new points to be covered during the AGPFC iterative procedure.
 */
void ArtGallerySolver::addIterationGridPoints() {
    std::vector<Point> grid;
    std::list<PolygonWithHoles> notCoveredList;

    _notCovered.polygons_with_holes(std::back_inserter(notCoveredList));

    std::cout << "[ArtGallerySolver] # of uncovered regions: " << notCoveredList.size() << std::endl;
    
    std::list<PolygonWithHoles>::iterator it;
    /* Seeks for non simple uncovered regions */
    for (it = notCoveredList.begin(); it != notCoveredList.end(); ++it) {
        Polygon tmpPol = it->outer_boundary();
        if (!tmpPol.is_simple() && it->number_of_holes() == 0) {
            vector<PolygonExt> vet = splitNonSimple(tmpPol);
            for (int i=0; i<vet.size(); i++) {
                PolygonWithHolesExt newPol = vet[i];
                notCoveredList.push_back(newPol);
            }
            it = notCoveredList.erase(it);
        }
    }
    
    
    for (it = notCoveredList.begin(); it != notCoveredList.end(); ++it) {
        RT xCentral = 0;
        RT yCentral = 0;
        Polygon tmpPol = it->outer_boundary();
        RT numVertex = RT(tmpPol.size());

        if (tmpPol.is_convex() && it->number_of_holes() == 0) {
            Polygon::Vertex_iterator vIt = tmpPol.vertices_begin();
            for (; vIt != tmpPol.vertices_end(); vIt++) {
                xCentral += (*vIt).x();
                yCentral += (*vIt).y();
            }
            xCentral /= numVertex;
            yCentral /= numVertex;
            _gridPoints.push_back(Point(xCentral,yCentral));
            grid.push_back(Point(xCentral,yCentral));
        } else {
            /* Changes the behavior of the algorithm when we find a non convex
             * uncovered area */

            /* First thing is to identify a convex vertex v and its adjacents. */
            Polygon::Vertex_circulator vIminus1 = tmpPol.vertices_circulator();
            Polygon::Vertex_circulator vI = tmpPol.vertices_circulator();
            Polygon::Vertex_circulator vIplus1 = tmpPol.vertices_circulator();
            for(;(*vIplus1) != (*vI); vIplus1++);
            for(;(*vIminus1) != (*vI); vIminus1--);
            vIplus1++;
            vIminus1--;
            Point convexVertex = *vI;

            do {
                if (CGAL::left_turn(*vIminus1, *vI, *vIplus1)) {
                    convexVertex = *vI;
                    break;
                }
                vIplus1++;
                vIminus1++;
                vI++;
            } while (*vI != convexVertex);

            Line lineAB(*vIplus1, *vIminus1);
            Polygon polAVB;
            polAVB.push_back(*vIminus1);
            polAVB.push_back(*vI);
            polAVB.push_back(*vIplus1);

            RT maxDist = RT(-1);
            Point maxVertex;

            /* For each other vertex q, if q is inside polAVB, compute distance to q (orthogonal to ab). 
             * If distance is a new maximum, save the point. */
            vI++;
            vI++;
            for (; *vI != *vIminus1; vI++) {
                if (polAVB.bounded_side(*vI) == CGAL::ON_BOUNDED_SIDE) {
                    RT dist = CGAL::squared_distance(*vI, lineAB);
                    if (dist > maxDist) {
                        maxVertex = *vI;
                        maxDist = dist;
                    }
                }
            }

            /* Verify holes */    
            for(PolygonWithHoles::Hole_const_iterator holeIt = it->holes_begin(); holeIt != it->holes_end(); ++holeIt) {
                for (Polygon::Vertex_iterator vIt = (*holeIt).vertices_begin(); vIt != (*holeIt).vertices_end(); vIt++) {
                    if (polAVB.bounded_side(*vIt) == CGAL::ON_BOUNDED_SIDE) {
                        RT dist = CGAL::squared_distance(*vIt, lineAB);
                        if (dist > maxDist) {
                            maxVertex = *vIt;
                            maxDist = dist;
                        }
                    }
                }
            }
            
            /* If no point is inside, returns midpoint of ab, or centroid of polAVP.
              Otherwise, if some point inside qv is internal, return its midpoint. */
            if (maxDist == RT(-1)) {
                xCentral = ((convexVertex).x() + (*vIminus1).x() + (*vIplus1).x()) / RT(3);
                yCentral = ((convexVertex).y() + (*vIminus1).y() + (*vIplus1).y()) / RT(3);
            } else {
                xCentral = ((convexVertex).x() + (maxVertex).x()) / RT(2);
                yCentral = ((convexVertex).y() + (maxVertex).y()) / RT(2);
            }
            grid.push_back(Point(xCentral,yCentral));
            _gridPoints.push_back(Point(xCentral,yCentral));
        }
    }
    
    std::cout << "[ArtGallerySolver] " << grid.size() << " new witnesses were chosen..." << std::endl;

    /* Constructs the matrix of visibility */
    int nRows = grid.size();
    int nCols = _guardCandidates.size();
    int qtdElem = 0;

    int previousSize = _matrix.size();
    tr1::unordered_map<pair<Point, Point>, bool>::iterator itMap;
    tr1::unordered_map<Point, PolygonExt>::iterator itVis;
    bool areVisible;

    _matrix.resize(previousSize + nRows);
    for(int i = previousSize; i < (nRows + previousSize); i++) {
        _matrix[i].resize(nCols);
    }

    for(int i = 0; i < nRows; i++) {
        Point p = grid.at(i);
        PolygonExt visPolCache = _polygon.getVisibility(p);
        _visCache.insert(pair<Point, PolygonExt>(p, visPolCache));

        for(int j = 0; j < nCols; j++) {
            pair<Point, Point> pointsPair(_guardCandidates.at(j), p);
            
            if (visPolCache.bounded_side(_guardCandidates.at(j)) == CGAL::ON_UNBOUNDED_SIDE)
                areVisible = false;
            else
                areVisible = true;

            if (areVisible) {
                _matrix[i + previousSize][j] = true;
                qtdElem++;
            } else {
                _matrix[i + previousSize][j] = false;
            }
        }
    }

    _solverBox = new PreSolver(_matrix, _groups, (double) _guards.size());

    return;
}

/**
 * Solves the current SCP instance.
 */
RT ArtGallerySolver::solveIteration() {
    double t0, t1;
    tr1::unordered_map<Point,PolygonExt>::iterator itVis;
    PolygonExt visPolCache;
    vector<int> guardsIndex; 

    t0 = CPUTIME(ruse);
    _solverBox->solve(_solverMode);
    t1 = CPUTIME(ruse);
    _scpResolTime += t1 - t0;

    _IPSolved ++;

    vector<int> opt = _solverBox->getOptimalSolution();
    _isLastOptimal = _solverBox->isOptimal();
    delete _solverBox;

    _polCovered = PolygonSet();
    PolygonSet thePol(_polygon);
    _guards.clear();
    _visGuards.clear();

    for (int k=0; k<opt.size(); k++) {
        _guards.push_back(_guardCandidates.at(opt[k]));
        guardsIndex.push_back(opt[k]);

        itVis = _visCache.find(_guardCandidates.at(opt[k]));
        if (itVis != _visCache.end()) {
            visPolCache = itVis->second;
        }
        else {
            visPolCache = _polygon.getVisibility(_guardCandidates.at(opt[k]));
            _visCache.insert(pair<Point,PolygonExt>(_guardCandidates.at(opt[k]),visPolCache));
        }
        _visGuards.push_back(visPolCache);
        _polCovered.join(visPolCache);
    }

    std::list<PolygonWithHoles> notCovered;
    thePol.difference(_polCovered);
    _notCovered = thePol;
    thePol.polygons_with_holes(std::back_inserter(notCovered));

    std::list<PolygonWithHoles>::const_iterator it;

    RT areaMiss = 0;

    for (it = notCovered.begin(); it != notCovered.end(); ++it) {
        areaMiss += it->outer_boundary().area();
    }

    std::cout << "[ArtGallerySolver] Area miss: " << areaMiss.to_double() << std::endl;

    _areaMiss = areaMiss;

    return areaMiss;
}

/**
 * Initializes the matrix used in SCP resolution.
 */
void ArtGallerySolver::initSolver(PolygonWithHolesExt pol, vector<Point> guardCandidates) {
    double t0 = CPUTIME(ruse);
    double t1;

    if (_polygon.outer_boundary().orientation() == CGAL::CLOCKWISE){
        _polygon.outer_boundary().reverse_orientation();
    }

    _remaining = RT(-1);
    _iteration = 0;
    _notCovered = PolygonSet();
    _polCovered = PolygonSet();
    _notCoveredFirst = PolygonSet();
    _guards.clear();
    _areaMiss = pol.area();

    _notCovered.clear();
    _polCovered.clear();
    _notCoveredFirst.clear();

    /* Constructs the matrix of visibility */
    std::vector<Point> grid = _hiddenPoints;
    int nRows = grid.size();
    int nCols = _guardCandidates.size();

    _matrix.resize(nRows);
    for (int i=0; i<nRows; i++) {
        _matrix[i].resize(nCols);
    }

    tr1::unordered_map<pair<Point, Point>, bool>::iterator itMap;
    tr1::unordered_map<Point,PolygonExt>::iterator itVis;
    bool areVisible;
    int reap = 0;
    int notReap = 0;

    if (_mainLoopIterations == 1) {
        for(int j = 0; j < nCols; j++) {
            Point pVis = _guardCandidates.at(j);

            for(int i =0; i < nRows; i++) {
                Point p = grid.at(i);

                pair<Point, Point> pointsPair(pVis, p);

                itVis = _visCache.find(p);
                PolygonExt visPolCache = itVis->second;

                if (visPolCache.bounded_side(pVis) == CGAL::ON_UNBOUNDED_SIDE)
                    areVisible = false;
                else
                    areVisible = true;

                _visibilityTest.insert(pair<pair<Point, Point>, bool>(pointsPair, areVisible));

                if (areVisible) {
                    _matrix[i][j] = true;
                }
                else {
                    _matrix[i][j] = false;
                }

            }
        }
    }
    else {    
        for(int j = 0; j < nCols; j++) {
            Point pVis = _guardCandidates.at(j);

            for(int i =0; i < nRows; i++) {
                Point p = grid.at(i);

                pair<Point, Point> pointsPair(pVis, p);

                itMap = _visibilityTest.find(pointsPair);

                if(itMap == _visibilityTest.end()){

                    itVis = _visCache.find(p);
                    PolygonExt visPolCache = itVis->second;

                    if (visPolCache.bounded_side(pVis) == CGAL::ON_UNBOUNDED_SIDE)
                        areVisible = false;
                    else
                        areVisible = true;

                    _visibilityTest.insert(pair<pair<Point, Point>, bool>(pointsPair, areVisible));
                    notReap++;
                }
                else{
                    areVisible = itMap->second;
                    reap++;
                }

                if (areVisible) {
                    _matrix[i][j] = true;
                }
                else {
                    _matrix[i][j] = false;
                }
            }
        }
    }

    _solverBox = new PreSolver(_matrix, _groups, _lowerBound);

    t1 = CPUTIME(ruse);
    _initSolverTime += t1 - t0;

    return;
}

/**
 * Solves an AGPFC instance using an iterative procedure.
 */
void ArtGallerySolver::runArt() {
    double t0 = CPUTIME(ruse);
    double t1;

    while ((_remaining != _zero) && (_guards.size() < _upperBound)) {
        stepArt();
        if (_iteration == 1){
            _areaMissFirst = _areaMiss;
            _notCoveredFirst.join(_notCovered);
            if (_isLastOptimal)
                _lowerBound = _guards.size();
            _guardsFirst = _guards;
            _hiddenPointsFirst = _hiddenPoints;
            _arrangementFirst = _arrangement;
            _lightAVPsFirst = _lightAVPs; 

            /* Compares new LB with old UB */
            if (_lowerBound == _upperBound || ((_mainLoopIterations+1) % AGPFC_FREQ) != 0) {
                t1 = CPUTIME(ruse);
                _runSolverTime += t1 - t0;
                return;
            }
        }
    }

    if (_iteration > _maxHorIt)
        _maxHorIt = _iteration;

    /* Sets new upper bound */
    if ((_guards.size() < _upperBound) && (_remaining == _zero)) {
        _upperBound = _guards.size();
        _solution = _guards;
    }

    _gridPoints.clear();

    t1 = CPUTIME(ruse);
    _runSolverTime += t1 - t0;

    return;
}

/**
 * Iterates in AGPFC algorithm.
 */
void ArtGallerySolver::stepArt() {
    if (_remaining != _zero) {
        ++_iteration;
        _remaining = solveIteration();

        if (_guards.size() >= _upperBound)
            return;

        if (_remaining > _zero) {
            addIterationGridPoints();
        }
    }
}

/**
 * Prints best solution found in file.
 */
void ArtGallerySolver::writeBestSol(char* idFile, std::ostream* solFile) {
    (* solFile) << idFile << endl;
    (* solFile) << _solution.size() << " ";
    for (int i=0; i<_solution.size(); i++) {
        (* solFile) << _solution[i] << " ";
    }
}

/**
 * Prints information about the solving process in file.
 */
void ArtGallerySolver::writeLog(char* idFile, std::ostream* logFile) {
    (* logFile) << "idFile polSize guards iterations initWitnSize finalWitnSize initCandSize finalCandSize maxHorizontal IPSolved initDisTime insertDisTime selectCandTime initSolverTime solverTime scpResolTime newDisTime totalTime" << endl;
    (* logFile) << idFile << " " << _polygon.size() << " " << getCardinality() << " " << _mainLoopIterations << " " << _initWitnSize << " " << _finalWitnSize << " " << _initCandSize << " " << _finalCandSize << " " << _maxHorIt << " " << _IPSolved << " " << _initDisTime << " " << _insertDisTime << " " << _selectGuardCandTime << " " << _initSolverTime << " " << _runSolverTime << " " << _scpResolTime << " " << _selectNewDisTime << " " << _procTime << endl;
}

/**
 * This function implements the full algorithm on a higher level code. The
 * algorithm follows an iterative process of obtaining lower and upper bounds
 * that lead to an optimal solution. These bounds are drawn up through the
 * resolution of instances of the aforementioned AGPW and AGPFC problems using
 * Integer Linear Programming (ILP) techniques.
 */
int ArtGallerySolver::solveProblem(PolygonWithHolesExt pol){ 
    double t0 = CPUTIME(ruse);
    double t1;
    vector<Point> dis;

    if (_solverMode == 0){
        _solverMode = 1;
    }
   
    cout << "[ArtGallerySolver] Instance: " << _polName << endl;
    cout << "[ArtGallerySolver] Parameters: " << endl;
    cout << "[ArtGallerySolver] == Witness Set Option: " << _sMode << endl;
    cout << "[ArtGallerySolver] == Solver Mode Option: " << _solverMode << endl;
    cout << endl;
   
    _visibilityTest.rehash(pol.size()*pol.size());
    _visCache.rehash(10*pol.size());

    _mainLoopIterations = 0;
    _lowerBound = 0;
    _polygon = pol;
    _upperBound = _polygon.size();

    _procTime = 0;
    _initDisTime = 0;
    _insertDisTime = 0;
    _selectGuardCandTime = 0;
    _initSolverTime = 0;
    _runSolverTime = 0;
    _scpResolTime = 0;
    _selectNewDisTime = 0;

    _IPSolved = 0;
    _maxHorIt = 0;

    if (_polygon.outer_boundary().orientation() == CGAL::CLOCKWISE){
        _polygon.outer_boundary().reverse_orientation();
    }        
    
    for(PolygonWithHoles::Hole_iterator holeIt = _polygon.holes_begin(); holeIt != _polygon.holes_end(); ++holeIt) {
        if ((*holeIt).orientation() != CGAL::CLOCKWISE){
            (*holeIt).reverse_orientation();
        }
    }

    cout << "[ArtGallerySolver] " << "Selecting initial discretization..." << endl;
    dis = selectInitialDisPoints();
    _initWitnSize = dis.size();

    cout << "[ArtGallerySolver] " << "Initiating Loop..." << endl;

    /* Main Loop */
    while(1){
        _mainLoopIterations++;
        
        if (_mainLoopIterations > MAX_ITERATIONS) {
            cout << "[ArtGallerySolver] Maximum number of iterations was exceeded!" << endl;
            exit(125);
        }

        cout << endl << "[ArtGallerySolver] " << "Initiating iteration " << _mainLoopIterations << endl;

        cout << "[ArtGallerySolver] " << "Inserting new witnesses..." << endl;
        insertDisPointsToArrange(dis);

        cout << "[ArtGallerySolver] " << "Selecting guards..." << endl;
        selectGuardCandidates();
        
        if (_mainLoopIterations == 1) {
            _finalCandSize = _guardCandidates.size();
            _initCandSize = _guardCandidates.size();
        }
        else {
            _finalCandSize = _guardCandidates.size();
        } 

        cout << "[ArtGallerySolver] " << "Initiating solver..." << endl;
        initSolver(pol, _guardCandidates);

        cout << "[ArtGallerySolver] " << "Running solver..." << endl;
        runArt();

        cout << "[ArtGallerySolver] " << "Comparing upper bound with lower bound..." << endl;
        cout << "[ArtGallerySolver] " << "Lower Bound: " << _lowerBound << endl;
        cout << "[ArtGallerySolver] " << "Upper Bound: " << _upperBound << endl;
            
        _finalWitnSize = _hiddenPointsFirst.size();
        
        if(_upperBound == _lowerBound){
            cout << "[ArtGallerySolver] " << "Found optimal solution!!" << endl;
            break;
        }

        cout << "[ArtGallerySolver] " << "Selecting new witnesses..." << endl;
        dis = selectNewDisPoints();
        
        t1 = CPUTIME(ruse); 
        _procTime = t1-t0;
        writePartialLog();

    }

    delete _guardsGrid;

    t1 = CPUTIME(ruse); 
    _procTime = t1-t0;

    printf("[ArtGallerySolver] CPU time: %f secs.\n", t1-t0);
    return _upperBound;
}

/**
 * Prints information about the partial solving process in file.
 */
void ArtGallerySolver::writePartialLog() {
    ofstream log;
    char buffer[256];
    sprintf(buffer, "%s.par", _nameLog);

    log.open(buffer);

    log << "idFile polSize LB UB iterations initWitnSize finalWitnSize initCandSize finalCandSize maxHorizontal IPSolved initDisTime insertDisTime selectCandTime initSolverTime solverTime scpResolTime newDisTime totalTime" << endl;
    log << _polName << " " << _polygon.size() << " " << _lowerBound << " " << _upperBound << " " << _mainLoopIterations << " " << _initWitnSize << " " << _finalWitnSize << " " << _initCandSize << " " << _finalCandSize << " " << _maxHorIt << " " << _IPSolved << " " << _initDisTime << " " << _insertDisTime << " " << _selectGuardCandTime << " " << _initSolverTime << " " << _runSolverTime << " " << _scpResolTime << " " << _selectNewDisTime << " " << _procTime << endl;

    log.close();
}


int verifyParameters(char * witnessMode, char * solverMode) {
    int res = 1;
    
    if (strcmp(witnessMode, "ALL_VERTICES") && strcmp(witnessMode, "CHWA_POINTS_EXTENDED") && 
            strcmp(witnessMode, "CONVEX_VERTICES") && strcmp(witnessMode, "CHWA_POINTS")) {
        cout << "The only available options for initial witness discretization are: " << endl;
        cout << "  ALL_VERTICES" << endl;
        cout << "  CONVEX_VERTICES" << endl;
        cout << "  CHWA_POINTS" << endl;
        cout << "  CHWA_POINTS_EXTENDED" << endl;

        res = 0;
    }

    if (atoi(solverMode) <= 0){
        cout << "The solver mode parameter only accepts positive integer values" << endl;
        
        res = 0;
    }

    return res;
}

/**
 * Main function.
 */
int main(int argc, char ** argv){
    ArtGallerySolver art_davi = ArtGallerySolver(); 
    char * name_pol = NULL;
    char * name_log = NULL;
    char * witness_mode = NULL;
    char * solver_mode = NULL;
    ofstream log;
    ofstream sol;
    PolygonWithHolesExt pol;

    if (argc != 5) {
        cout << "[ArtGallerySolver] Usage: ./artGallerySolver <file .pol> <log file> <witness discretization (e.g. CHWA_POINTS)> <solver mode (e.g. 2)>" << endl;
        cout << "\t<file .pol> is the file representing the polygon. See the file format in Section 4.1 of the README file." << endl;  
        cout << "\t<log file> is the file where the summary log will be outputted. See the file format in Section 4.2 of the README file." << endl;
        cout << "\t<witness discretization> consists in the technique chosen to select the initial Witness Set of our method. Current possibilities are:" << endl;
        cout << "\t\tALL_VERTICES: includes all vertices of the polygon of input" << endl;
        cout << "\t\tCONVEX_VERTICES: includes all convex vertices of the polygon of input" << endl;
        cout << "\t\tCHWA_POINTS: assembles the initial witness set for our algorithm from the midpoints of all reflex-reflex edges and all convex vertices from convex-reflex edges" << endl;
        cout << "\t\tCHWA_POINTS_EXTENDED: uses a combination of points in CHWA_POINTS with the reflex vertices of the polygon" << endl;
        cout << "\t<solver mode> is a positive integer that represents the mode used to solve SCP instances. Current possibilities are:" << endl;
        cout << "\t\t1: GLPK + Lagrangian Heuristic" << endl;
        cout << "\t\t2: GLPK" << endl;
        cout << "\t\t3: XPRESS + Lagrangian Heuristic" << endl;
        cout << "\t\t4: XPRESS" << endl;

        exit(1);
    }
    
    name_pol = argv[1];
    name_log = argv[2];
    witness_mode = argv[3];
    solver_mode= argv[4];
    
    if (!verifyParameters(witness_mode, solver_mode))
        exit(1);

    pol = art_davi.loadPol(name_pol);

    art_davi.setMode(witness_mode);
    art_davi.setSolverMode(atoi(solver_mode));
    art_davi.setLogFile(name_log);

    cout << "[ArtGallerySolver] --Initiating Iterative Solution--" << endl;

    art_davi.solveProblem(pol);

    cout << "[ArtGallerySolver] --Finalizing Iterative Solution--" << endl;

    log.open(name_log); 
    art_davi.writeLog(name_pol, &log); 
    log.close();

    sol.open((getFileName(string(name_pol)) + ".sol").c_str()); 
    art_davi.writeBestSol(name_pol, &sol); 
    sol.close(); 

    return 0;
}
