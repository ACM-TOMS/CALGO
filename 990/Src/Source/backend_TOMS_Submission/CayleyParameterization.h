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
#ifndef CAYLEYPARAMETERIZATION_H_
#define CAYLEYPARAMETERIZATION_H_

#include "ActiveConstraintGraph.h"
#include "PointSet.h"
#include "CayleyPoint.h"

#include <string>
#include <vector>
#include <utility> // for pair/*
 * This class chooses non-edges in an ActiveConstraintGraph that complete the graph into 3-tree.
 * Those non-edges are called the parameters. The complexity of the sampling algorithm varies
 * based on the choice of non-edges and the order in which they are fixed.
 */

class CayleyParameterization {
public:

	/**
	 * maximum number of atoms that can contribute to activeConstraintGraph for 2 rigid body packing case.
	 * For n=2, dof is 6. Hence there can be 6 contact pairs with each owning distinct atoms that results in 12 different atoms.
	 */
	static const int NO = 12;

	/**
	 * @brief Constructor with active constraint graph initialization
	 * It calls the helper methods to determine the parameters that make the graph 3-tree
	 *
	 * @param basic True keeps the parameters unordered to fasten the process, False built3tree to order parameters.
	 */
	CayleyParameterization(ActiveConstraintGraph* cgK, bool basic = true);

	/** @brief Destructor */
	virtual ~CayleyParameterization();

	/////////////////////////
	// Getter/Setter
	////////////////////////

	/**
	 * @return Non-edges of the ActiveConstraintGraph
	 * i.e. set of parameters, each represented by <helix A "vertex" index , helix B "vertex" index>
	 */
	vector<pair<int, int> > getParameters();

	/** @return True if ActiveConstraintGraph is partial 3-tree, False otherwise */
	bool is_partial3tree();

	/** @return Vector of tetrahedrons each represented by 4 vertex indices. */
	vector<vector<int> > getTetras();

	/** @return Adjacency map for dependency of parameters */
	vector<vector<int> > getUpdateList();

	/**
	 * @brief Finding bounds for each non-edge is either solving a linear inequality
	 * or solving a single quadratic of one variable.
	 *
	 * @return Set of integers corresponding to each parameter that express
	 * what inequality is needed to compute parameter range.
	 * i.e. triangular, tetrahedral inequality or short-range.
	 */
	vector<int> getBoundaryComputationWay();

	/////////////////////////
	// Helper functions
	////////////////////////

	/**
	 * @brief The parameters of an active constraint graph are selected as maximal
	 * 3-realizable (3-tree) extension by leveraging the convex parametrization theory.
	 *
	 * @detail Checks if the contacts of currentGraph is subset of any isomorphisms of the graphs in the list of complete3trees
	 * If found then the edges in that graph (except contacts) will be used as parameter set.
	 * If the contacts are not part of any pre-determined graphs i.e. it is not partial 3-tree, then assign random parameters
	 *
	 * @return Set of parameters each represented by <helix A "vertex" index , helix B "vertex" index>
	 */
	virtual vector<pair<int, int> > defineParameters();

	/**
	 * @brief An alternative method to pick the parameters for 5 dimensional region such that
	 * the range of each parameter has similar length. It aims sampling more uniformly on
	 * Cartesian space by setting Cayley parameter space as spheric as possible.
	 *
	 * @detail Set parameters such that the deviation of parameter intervals is minimal
	 * i.e. choose a parameter set that will give almost a "spheric" parameter space.
	 * note that sampling occurs from min to max values of a parameter step by step.
	 * ---
	 * for 5d case, we have 1 contact and need 5 parameters. There are 8 non-edges in contact graph. we need to eliminate 3 non-edges
	 * to find that set:
	 * order those 8 non-edges according to length of intervals. The set of non-edges in the middle will give the parameters with min deviation.
	 * So eliminate the 3 non-edges from top or down such that remaining 5 non-edges will be our parameter set with the min deviation.
	 *
	 * @return Set of parameters each represented by <helix A "vertex" index , helix B "vertex" index>
	 *
	 * todo remove this function when submitting to TOMS???
	 * todo put the flag to run this method or not into the settings.ini file
	 */
	vector<pair<int, int> > parameterMinDeviation();

	/**
	 * @brief 3-tree formed by starting with a 4-vertex complete graph and then repeatedly
	 * adding vertices in such a way that each added vertex has exactly 3 neighbors that form a clique.
	 *
	 * Finds a 4 clique as a base tetrahedra. Then complete 3-tree is
	 * built up from a base tethedra by adding, at each step, a new vertex edge-connected
	 * to the face of a tetrahedra. Store the tetrahedras in the order they are created in
	 * attribute tetrahedra.
	 *
	 * previous name was createTetrahedras
	 */
	void built3tree();

	/**
	 * @brief Search for a 4-vertex complete graph
	 * Then updates its vertex_position s to be True(known)
	 * And create a new tetrahedron with these 4 vertices (base tetrahedron)
	 * Then add this tetrahedron to tetrahedra array.
	 *
	 * @return True if able to find, False otherwise
	 */
	bool find4Clique();

	/**
	 * @brief Search for a vertex that is connected to the face of a tetrahedron
	 * Then updates its vertex_position to be True(known)
	 * And create a new tetrahedron with this vertex plus the 3 vertices of the face that it is connected to.
	 * Then add this tetrahedron to tetrahedra array.  HERE ADD THE LABEL FOR MISSING EDGE ???!!!
	 *
	 * @return The index of vertex that is just identified to be known by location wise
	 *
	 * previous name was findPointThatIsConnectedToFaceOfTetraHedra  or findPointsAcordingToKnownLocations
	 */
	int identifyRigidVertex();

	/**
	 * @brief Search for any undefined edges connected to the vertex if the neighbor vertex is also location-wise known.
	 * Then updates this edge to be True(known)
	 *
	 * @param vertex The index of vertex which has recently added to 3-tree
	 *
	 * previous name was findLengthsUsingTwoKnownPoints
	 */
	void identifyFixedEdgesConnectedToVertex(int vertex);

	/**
	 * @bried Display the information related tetrahedrons on the console
	 */
	void printTetras();

	/*
	 * @brief Order parameters to guarantee efficient sampling.
	 * Also sets the updateList and boundaryComputationWay arrays.
	 *
	 * The parameter in the first tetrahedra will be at the end, if the parameter takes place
	 * in the same tetrahedra, then max degree will be put to the end
	 *
	 * starts from first tetrahedra.
	 * finds the parameter with max degree then puts it at the end
	 * then find the next parameter with max degree, then put it before end.
	 * then go to next tetrahedra, do same staff.
	 * till paramss empty.
	 *
	 * so this algorithm orders the parameters in the reverse order they take place in the tetrahedras.
	 * if there are multiple parameters in one tetrahedra, then the parameter with higher degree will be put to the end side
	 */
	void orderParameters();

private:

	ActiveConstraintGraph* currentGraph;

	PointSet* helA;
	PointSet* helB;

	/*
	 * Set of Contacts
	 * A Contact is a pair of integer where each integer is the "index of the vertex in constraint graph".
	 * i.e. A Contact is represented by <helix A "vertex" index , helix B "vertex" index>
	 *
	 * i th vertex from helix B is represented by i+NO/2 as index number
	 */
	vector<pair<int, int> > contacts;

	/*
	 * Set of Parameters
	 * i.e. non-edges in an ActiveConstraintGraph that complete the graph into 3-tree.
	 *
	 * A Parameter is a pair of integer where each integer is the "index of the vertex in constraint graph".
	 * i.e. A Parameter is represented by <helix A "vertex" index , helix B "vertex" index>
	 *
	 * i th vertex from helix B is represented by i+NO/2 as index number
	 */
	vector<pair<int, int> > parameters;
// vector<int> *v

	/*
	 * vector of contact + parameter "atom" indices for the vertices of contact graph
	 * (atoms ,that are not in the actual contact list but hold for parameters, also exist in this list)
	 * the unique vertex list of the contact graph
	 */
	vector<int> verticesA;
	vector<int> verticesB;

	/**
	 * Vector of ordered tetrahedrons that helps to define the order of parameters
	 * that will be important later during sampling procedure. This data is later on
	 * passed to ConvexChart and CartesianRealizer to help their computations.
	 *
	 * Each tetrahedron are represented by 4 vertex indices.
	 * Each sub vector hold the indices of vertices of one tetrahedron.
	 *
	 * i th vertex from helix B is represented by i+NO/2 as index number
	 *
	 * @see built3tree()
	 */
	vector<vector<int> > tetrahedra;

	/**
	 * shows if the length of the edge between vertices of ActiveConstraintGraph is known or not
	 *
	 * The following lengths are known:
	 * - interatomic distances inside of helixA
	 * - interatomic distances inside of helixB
	 * - the distance between atoms that are in contact
	 * - the length of the edge that is used as Cayley parameter
	 * - the length of the edge between vertices that has known locations
	 *
	 * i th vertex from helix B is represented by i+NO/2 as index number
	 *
	 * previous name was p
	 */
	bool edge[NO][NO];

	/**
	 * shows if the length of the edge is known by the vertices coordinates that it is connected to.
	 *
	 * i th vertex from helix B is represented by i+NO/2 as index number
	 */
	bool fixedEdge[NO][NO];

	/**
	 * show known positions of vertices of ActiveConstraintGraph
	 *
	 * The following vertices have identified positioning:
	 * - 4 vertices of initial 4-clique that is used as base tetrahedron
	 * - vertices that have exactly 3 neighbors that form a clique.
	 * i.e. vertices that are edge-connected to the face of a tetrahedron.
	 *
	 * Note that i th vertex from helix B is represented by i+NO/2 as index number
	 *
	 * previous name was pos
	 */
	bool vertex_position[NO];

	/**
	 * adjacency map for dependency of parameters. It provides the set of
	 * parameters whose range will be updated when one parameter fixed.
	 */
	vector<vector<int> > updateList;

	/**
	 * inequalities that express the range of a parameter can be classified into
	 * either a linear or non-linear class. This variable is characterization
	 * of the parameter that tells what inequality is needed to compute parameter range.
	 * i.e. triangular or tetrahedral inequality.
	 *
	 * Set of integers corresponding to each parameter that takes values as:
	 * -2 if the parameter is actually a contact in case of short range sampling
	 * -1 if triangular inequality is used to compute the range of the parameter
	 * the index of the tetrahedron in tetrahedra if tetrahedral inequality is used to compute the range of the parameter
	 */
	vector<int> boundaryComputationWay;

	/** A flag to show whether ActiveConstraintGraph is partial 3-tree or not. */
	bool partial3tree;

};

#endif /* CAYLEYPARAMETERIZATION_H_ */
