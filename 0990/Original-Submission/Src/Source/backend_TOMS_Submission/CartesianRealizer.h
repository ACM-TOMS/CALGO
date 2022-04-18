/*
 This file is part of EASAL. 

 EASAL is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 EASAL is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef CARTESIANREALIZER_H_
#define CARTESIANREALIZER_H_

#include "Orientation.h"
#include "ActiveConstraintGraph.h"
#include "ConvexChart.h"

#include <vector>

/**
 * This class contains routines that computes Orientation that represents transformations
 * of the rigid molecular composite s relative to each other. It computes Cartesian realization
 * of an ACG with the parameter lengths taken from cayleyPoint ? and active constraint lengths
 * for a specific flip. It intentionally ignores the remainder of the assembly constraint system,
 * namely atom markers not in G and their constraints.
 */
class CartesianRealizer {
public:

	/**
	 * @brief computes the Orientation by leveraging partial 3-tree techniques.
	 * ActiveConstraintGraph which is a complete 3-tree is built up from a base tetrahedron
	 * by adding, at each step, a new vertex edge-connected to the face of a tetrahedron.
	 *
	 * @param flipno Determines which orientation of a realization to compute
	 * @param[out] fail True if it is non-realizable or distorted
	 */
	static Orientation* computeRealization(ActiveConstraintGraph* cgk,
			ConvexChart* descr, int flipno, bool& fail);
	static Orientation* computeRealization(ActiveConstraintGraph* cgk,
			ConvexChart* descr, int flipno, double paramconnections[][12],
			bool& fail);

private:
	/** @brief clean toA and toB members of the class */
	static void destructor();

	/** @brief Checks if the transformation from fromA to toA corrupts or not*/
	static bool isDistorted();

	/** @brief Creates a new Orientation with fromB, toB and flip initialization */
	static Orientation* getOrienation();

	/**
	 * @brief Finds Cartesian coordinates of the vertices of base(first) tetrahedron by known edge lengths
	 *
	 * The vertices are achieved from tetrahedra[0]
	 *
	 * @param mirror Determines which side of the face to put the 4th vertex.
	 */
	static bool setBaseTetra(bool mirror);

	/**
	 * @brief Sets the distance of any undefined edge connected to the vertex
	 * if the neighbor vertex is also location-wise known.
	 *
	 * @param vertex The index of vertex whose Cartesian location is recently computed
	 *
	 * previous name was findLengthsUsingTwoKnownPoints
	 */
	static void computeLengthOfTheEdgesConnectedToVertex(int vertex);

	/**
	 * @brief Finds Cartesian coordinates of the vertex that is connected to the face of a tetrahedron
	 * The indices of vertices are taken from 'tetrahedron_index' th tetrahedron
	 * Computes Cartesian coordinates of the FIRST vertex of the tetrahedron.
	 * Uses last 3 vertices as triangular face to connect.
	 *
	 * First, locates the tetrahedron (i.e. 4 vertices) on the space without using pre-computed location of the vertices of the face.
	 * Then kind of computes transformation of the face and then apply same transformation to first vertex
	 *
	 *
	 * @param tetrahedron_index The index of tetrahedron
	 * @param mirror Determines which side of the face to put the 4th vertex.
	 *
	 * @see locateTetrahedron()
	 * @see Utils::matApp()
	 *
	 * previous name was findLocation
	 */
	static bool locateVertex(int tetrahedron_index, bool mirror);

	/**
	 * @brief This method uses a set of 6 lengths and outputs the locations of the 4 vertices
	 * The vertices are achieved from tetrahedra[tetrahedron_index]
	 * The lengths of edges are achieved from edge_length array
	 *
	 *  vertex indices           edge indices
	 *	   2---3                  +-1-+   /\
	 *	   |   |                  2   3  4  5
	 *	   0---1                  +-0-+ /    \
	 *
	 *	First point is located at the origin
	 *	Second point is located on x axis
	 *	Third point is located on x-y plane
	 *	Fourth point is located on x-y-z plane
	 *
	 *  @param tetrahedron_index The index of tetrahedron
	 *	@param mirror Determines which side of the z-axis to put 4th point.
	 *	@param[out] positions The locations of the 4 vertices
	 *
	 *	@see Utils::volumeTetra()
	 *	@see Utils::lenToTetra()
	 *
	 *	previous name was mylenToTetra
	 */
	static void locateTetrahedron(int tetrahedron_index, double positions[][3],
			bool mirror);

	/** @brief Computes the translation and Euler angles needed for Jacobian computation */
	static void compute_TB_EAB();

	/** @brief Print Cartesian location of indices and the lengths of edges in between those vertices for debug purpose*/
	static void printPositionsAndConnections(vector<int> ind);

	/////////////////////////
	// Variables
	////////////////////////

	/**
	 * maximum number of atoms that can contribute to activeConstraintGraph for 2 rigid body packing case.
	 * For n=2, dof is 6. Hence there can be 6 contact pairs with each owning distinct atoms that results in 12 different atoms.
	 */
	static const int NO = 12;

	static bool partial3tree;

	/**
	 * Vector of tetrahedrons. Each tetrahedron are represented by 4 vertex indices.
	 * Each sub vector hold the indices of vertices of one tetrahedron.
	 *
	 * @see CayleyParameterization::built3tree()
	 */
	static std::vector<std::vector<int> > tetrahedra;

	static Eigen::Vector3d TB;
	static Eigen::Vector3d eaB;

	static bool debug;

	/** false if volume is negative */
	static bool realizable;
	// todo check if volume computations are accurate enough and check the error messages printed to the console

	static PointSet* helA, *helB;

	/*
	 * vector of contact + parameter "atom" indices for the vertices of contact graph
	 * (atoms ,that are not in the actual contact list but hold for parameters, also exist in this list)
	 * the unique vertex list of the contact graph
	 */
	static std::vector<int> verticesA;
	static std::vector<int> verticesB;

	static ActiveConstraintGraph* graph;

	/**
	 * For 3*3 graphs there are 3 tetrahedrons
	 * For 3*4 graphs there are 4 tetrahedrons. One of the tetrahedron is fixed. (The one where all the vertices located in one helix)
	 * For 3*5 graphs there are 5 tetrahedrons. Two of the tetrahedrons are fixed. (The ones where all the vertices located in one helix)
	 * For 3*6 graphs there are 6 tetrahedrons. Three of the tetrahedrons are fixed. (The ones where all the vertices located in one helix)
	 *
	 * there are 3 tetrahedrons hence 2^3=8 flips.
	 *  If all vertices of the tetrahedron is in the same rigid molecule, then it does not contribute to the flip set.
	 *  todo checks the effects of this situation
	 *
	 *  for 2 helices, there are only 3 tetrahedrons satisfying above and being partial 3-tree graph.
	 */
	static int flipNo;

	/** Cartesian coordinates of first three vertices of ActiveConstraintGraph in helixA */
	static double *fromA[3];

	/** Cartesian coordinates of first three vertices of ActiveConstraintGraph in helixB */
	static double *fromB[3];

	/** Transformed Cartesian coordinates of first three vertices of ActiveConstraintGraph */
	static double *toA[3], *toB[3];

	/** Cartesian coordinates of vertices in ActiveConstraintGraph */
	static double positions[NO][3];

	/** All fixed distances plus current distance values of non-edges of ActiveConstraintGraph
	 * previous name was pconnections
	 */
	static double edge_length[NO][NO];

};

#endif /* CARTESIANREALIZER_H_ */
