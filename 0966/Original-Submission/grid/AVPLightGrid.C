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


#include "AVPLightGrid.h"
        
/**
 * Destrutor. Frees arrangement's and observer's memory
 */
AVPLightGrid::~AVPLightGrid() {
    delete _arr;
    delete _obs;
}

/**
 * Adds edges from new visibility polygons in the arrangement.
 */
void AVPLightGrid::addToGrid(std::vector<PolygonExt> vet){
    void *lixo = malloc(1);
    struct mallinfo info;
    int mem1, mem2;
    free(lixo);
    info = mallinfo();
    mem1 = info.uordblks;
    
    double tmp1 = CPUTIME(ruseTime);
    double tmp2 = CPUTIME(ruseTime);
    double tmpArr = 0.0;
    double tmpStruct = 0.0;

    _struct1.clear();
    _struct2.clear();

    int countVis = 0;
    int visSize = vet.size();
    for (; countVis < visSize; countVis++) {
        PolygonExt visPol = vet.at(countVis);

        for (Polygon::Edge_const_iterator eIt = visPol.edges_begin(); eIt != visPol.edges_end(); eIt++) {
            _obs->setOriginalSeg((*eIt));
            CGAL::insert((*_arr), (*eIt));
        }
    }

    tmp2 = CPUTIME(ruseTime);
    tmpArr = tmp2 - tmp1; 

    _vertexNum = _arr->number_of_vertices();
    _edgeNum = _arr->number_of_edges();
    _faceNum = _arr->number_of_faces();

    Arrangement::Ccb_halfedge_const_circulator  curr, first;
    Arrangement::Face_const_iterator            fit;

    tmp1 = CPUTIME(ruseTime);
    for (fit = _arr->faces_begin(); fit != _arr->faces_end(); ++fit) {
        if (!fit->is_unbounded()) {
            bool shadowFace = true;
            curr = fit->outer_ccb();
            first = curr;
            RT xCentral;
            RT yCentral;
            RT numVertex;
            PolygonExt pol;
            bool enterGrid = false;
            do {
                if (curr->data() != -1) {
                   enterGrid = true;
                }
                shadowFace = (shadowFace & (curr->data() != 1)); //0

                numVertex += 1;
                xCentral += curr->target()->point().x();
                yCentral += curr->target()->point().y();
                pol.push_back(curr->target()->point());
                ++curr;
            } while(curr != first);
            xCentral /= numVertex;
            yCentral /= numVertex;
            if (shadowFace) {
                if (enterGrid) {
                    _grid.push_back(Point(xCentral,yCentral));
                    _struct2.push_back(pol);
                }
            } else {
                if (enterGrid) {
                    _struct1.push_back(pol);
                }
            }
        }
    }
    tmp2 = CPUTIME(ruseTime);
    tmpStruct = tmp2 - tmp1; 

    info = mallinfo();
    mem2 = info.uordblks;
    _gridMem = mem2 - mem1;

    _gridTime = tmpArr + tmpStruct;
}

/**
 * Constructs the grid.
 */
void AVPLightGrid::makeGrid() {
    void *lixo = malloc(1);
    struct mallinfo info;
    int mem1, mem2;
    free(lixo);
    info = mallinfo();
    mem1 = info.uordblks;
    
    double tmp1 = CPUTIME(ruseTime);
    double tmp2 = CPUTIME(ruseTime);
    double tmpArr = 0.0;
    double tmpStruct = 0.0;

    _arr = new Arrangement();
    _obs = new MyObserver((*_arr));

    for (Polygon::Edge_const_iterator eIt = _polygon.outer_boundary().edges_begin(); eIt != _polygon.outer_boundary().edges_end(); eIt++) {
        _obs->setOriginalSeg((*eIt));
        CGAL::insert((*_arr), (*eIt));
    }

    PolygonWithHoles::Hole_const_iterator holeIt = _polygon.holes_begin();
    for(; holeIt != _polygon.holes_end(); ++holeIt) {

        Polygon holeItPol(*holeIt);
        for (Polygon::Edge_const_iterator eIt = holeItPol.edges_begin(); eIt != holeItPol.edges_end(); eIt++) {
            _obs->setOriginalSeg((*eIt));
            CGAL::insert((*_arr), (*eIt));
        }
    }

    _obs->setStarterEdge(false);
    int countVis = 0;
    int visSize = _visPol.size();
    for (; countVis < visSize; countVis++) {
        PolygonExt visPol = _visPol.at(countVis);

        for (Polygon::Edge_const_iterator eIt = visPol.edges_begin(); eIt != visPol.edges_end(); eIt++) {
            _obs->setOriginalSeg((*eIt));
            CGAL::insert((*_arr), (*eIt));
        }
    }

    tmp2 = CPUTIME(ruseTime);
    tmpArr = tmp2 - tmp1; 

    _vertexNum = _arr->number_of_vertices();
    _edgeNum = _arr->number_of_edges();
    _faceNum = _arr->number_of_faces();

    Arrangement::Ccb_halfedge_const_circulator  curr, first;
    Arrangement::Face_const_iterator            fit;

    tmp1 = CPUTIME(ruseTime);
    for (fit = _arr->faces_begin(); fit != _arr->faces_end(); ++fit) {
        if (!fit->is_unbounded()) {
            bool shadowFace = true;
            curr = fit->outer_ccb();
            first = curr;
            RT xCentral;
            RT yCentral;
            RT numVertex;
            PolygonExt pol;
            bool enterGrid = false;
            do {
                if (curr->data() != -1) {
                   enterGrid = true;
                }
                shadowFace = (shadowFace & (curr->data() != 1)); //0

                numVertex += 1;
                xCentral += curr->target()->point().x();
                yCentral += curr->target()->point().y();
                pol.push_back(curr->target()->point());
                ++curr;
            } while(curr != first);
            xCentral /= numVertex;
            yCentral /= numVertex;
            if (shadowFace) {
                if (enterGrid) {
                    _grid.push_back(Point(xCentral,yCentral));
                    _struct2.push_back(pol);
                }
            } else {
                if (enterGrid) {
                    _struct1.push_back(pol);
                }
            }
        }
    }
    tmp2 = CPUTIME(ruseTime);
    tmpStruct = tmp2 - tmp1; 

    info = mallinfo();
    mem2 = info.uordblks;
    _gridMem = mem2 - mem1;

    _gridTime = tmpArr + tmpStruct;
}
        
