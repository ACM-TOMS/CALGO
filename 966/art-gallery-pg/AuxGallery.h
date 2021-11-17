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


#ifndef AUX_GALLERY_H
#define AUX_GALLERY_H

#include "PolygonExt.h"
#include "PolygonWithHolesExt.h"

#include <CGAL/copy_n.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_set_2.h>

#include <limits>
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>
#include <tr1/unordered_map>

typedef Polygon::Vertex_const_iterator Vertex_iterator;

using namespace std;

/**
 * Implements the hash operation of a <Point, Point> key.
 */
typedef struct
{
      long operator() (const pair<Point, Point> &k) const { 
          std::stringstream out;
          out << k.first << ":" << k.second;
          return std::tr1::hash<std::string>()(out.str());
      }
} PointPointHash;

/**
 * Implements the hash operation of a Point key.
 */
typedef struct
{
      long operator() (const Point &k) const { 
          std::stringstream out;
          out << k;
          return std::tr1::hash<std::string>()(out.str());
      }
} PointHash;

/**
 * If a hole is a non simple polygon, splits the same into simple ones.
 */
vector<PolygonExt> splitNonSimple(PolygonExt pol);

/**
 * Gets file name, whithout the whole path.
 */
string getFileName(string polFileName);

/**
 * Finds the size of the fraction x.
 */
int fractionSize(RT x);

/**
 * Finds an internal point of a Polygon With Holes.
 */
Point internalPoint(PolygonWithHoles polWHoles);

/**
 * Finds an internal point of a Polygon Without Holes.
 */
Point internalPoint(Polygon pol);

/**
 * Finds a simpler point inside the polygon.
 */
Point findSimplerInteriorPoint(PolygonWithHoles polWHoles, Point p);

/**
 * Finds a simpler point inside the polygon.
 */
Point findSimplerInteriorPoint(Polygon pol, Point p);

/**
 * Returns a truncated num/den.
 */
RT * truncateFraction(string num, string den, int size);

#endif
