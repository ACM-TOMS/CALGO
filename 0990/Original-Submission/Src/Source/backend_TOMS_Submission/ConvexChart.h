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
#ifndef CONVEXCHART_H_
#define CONVEXCHART_H_

#include "ActiveConstraintRegion.h"
#include "CayleyParameterization.h"

#include <vector>
#include <utility>

/*
 * This class contains the equations and inequalities of the MolecularUnit interaction.
 * It regulates the equations and defines the grid which the algorithm needs to step through.
 * This class is used to determine the chart that parameterize the regions. i.e.
 * it computes the range of parameters of ActiveConstraintGraph.
 * An exact convex chart yields feasible Cayley points for the current active constraint
 * region. The resulting Cayley configuration space is convex, before collisions or other
 * (e.g. angle) constraints are introduced. The range of parameters are computed by
 * triangle and tetrahedra inequalities.
 */

class ConvexChart {
public:
	/**
	 * maximum number of atoms that can contribute to activeConstraintGraph for 2 rigid body packing case.
	 * For n=2, dof is 6. Hence there can be 6 contact pairs with each owning distinct atoms that results in 12 different atoms.
	 */
	static const int NO = 12;

	/**
	 * @brief Constructor with active constraint graph, CayleyParameterization and PredefinedInteractions initialization
	 * And initializes the class variables and computeFixEdgeDistances.
	 *
	 * @param dense A flag for refined sampling
	 */
	ConvexChart(ActiveConstraintGraph* cgk, bool dense,
			CayleyParameterization* cparam, PredefinedInteractions *df);

	/** @brief Destructor */
	virtual ~ConvexChart();

	/**
	 * @brief Compute the lengths within the MolecularUnit by the MolecularUnit atom positions.
	 * Note that getAtomAt is only used for distances within the helix, since it gives the data as it was read rather than transformed.
	 * Lengths between Helices are defined by the radii.
	 */
	virtual void computeFixEdgeDistances();

	/**
	 * @brief It initializes the boundaries of convex chart as tight as possible.
	 *
	 * The Cayley point will be set to witness point if it is fresh sampling.
	 * Otherwise if it is continuation of an older sampling, then Cayley point will be
	 * set to the last point where it had stopped before.
	 */
	bool initializeChart(bool contin, ActiveConstraintRegion * region);

	/** @return True if sampling is completed, False otherwise */
	virtual bool doneGrid();

	/**
	 * @brief Sets parameter point to the next grid point within the computed range.
	 * It follows a direction i.e. sample along one Cayley parameter p.
	 * When it reach the end of that direction i.e. out of the range of the parameter p,
	 * then it walks one step on the next parameter p+1th while setting previous parameter
	 * (first p parameters) values to proximity of witness point.
	 *
	 * It also updates min max range for parameters from 0th till p'th  (not just p'th parameter)
	 * Range computation is required in every iteration for dependent parameters.
	 */
	virtual void stepGrid();

	/**
	 * @brief Sets parameter point to the neighbour grid point in dimension specified by 'dir'.
	 * Walk one step ahead (forward and backward) from the valid point along parameter specified by 'dir'.
	 * Allow the point to be less than edge_lengthLower since the purpose of this method is to
	 * search for a colliding configuration to start binary search in between.
	 *
	 * todo name in the paper was stepNeighbour()
	 */
	virtual bool stepAround();

	/**
	 * @brief Sets parameter point to the somewhere between current point and neighbour grid point according to binary
	 * search procedure in boundaryDetection.
	 *
	 * if stepAround causes collision, then we start binary search in between valid point and colliding point.
	 * if since we start binary search from the colliding point, first we need to walk back through the valid point.
	 *
	 * @param valid If True, walk through colliding point. If false, walk through valid point
	 */
	virtual void stepGridBinary(bool valid);

	/** @brief Reset/return back to grid point 'p' where it was before  after binary search */
	virtual void setInitialPoint(std::vector<double> p);

	/** @brief Reset the direction to -1 */
	virtual void setDir();

	/*those methods can be used from Jacobian sampling*/
	virtual void updateBoundaries();
	virtual void updateBoundariesSpecific(int x);
	virtual void updateBoundariesDependent(int x);

	/** @return Cayley point: the current location of sampling on Cayley space */
	std::vector<double> getPoint();

	/** @return Cayley point (achieved from param_length): the current location of sampling on Cayley space */
	std::vector<double> getParamPoint();

	/** @brief assigns wp to witnessPoint */
	void setWitnessPoint(CayleyPoint* wp);

	/** @return Vector of tetrahedrons each represented by 4 vertex indices. */
	std::vector<std::vector<int> > getTetras();

	/** @return Atom with index 'vertex' where i th vertex from helix B is represented by i+NO/2 as index number*/
	Point * getAtom(int vertex);

	/**
	 * @param ipair is not required to have the order of  <helixAvertex, helixBvertex> i.e. ipair can be set as 2-8 or 8-2
	 * @return -1 if ipair is not found in ivector, otherwise return the index of ipair in ivector
	 */
	int find_pair_in_vector(std::vector<std::pair<int, int> > ivector,
			std::pair<int, int> ipair);

private:
	/**
	 * @brief Computes the range of the non-edge v1-v2 through triangular inequalities.
	 * Uses all triangles that this edge is shared.
	 * the order of the parameter affects the result
	 */
	virtual void setRangeByTriangleInequality(int v1, int v2);

	/**
	 * @brief Computes the range of the non-edge v1-v2 through tetrahedral inequality using tetrahedron indexed by tetraNo
	 */
	virtual void setRangeByTetrahedralInequality(int v1, int v2, int tetraNo);

	/**
	 *  @brief Computes the range of the parameter v1-v2 in order to eliminate sampling infeasible grid points.
	 *  When there is one parameter in a tetrahedra, the range of parameter can be computed through tetrahedral inequalities.
	 *  If there are 2 unfixed parameters at the moment in a tetrahedra, then the range of one of the parameters will be computed
	 *  through triangular inequalities and fixed. And the range of other parameter will be computed
	 *  through tetrahedral inequalities. The range computation way depends on the order of parameters.
	 */
	void computeRange(int v1, int v2, int paramNo);

	/*
	 * @brief This function stores the valid values of edge_length into the param_length, param_lengthLower and param_lengthUpper.
	 * It also ensures that there is a diagonal mirror on each.
	 */
	virtual void setMinMaxParam(int v1, int v2);

	/** @return The step size computed according to Settings::Sampling::dynamicStepSizeWithin */
	double dynamicStepSize(int v1, int v2, int paramNo);

public:

	/** upper bound of the parameters range */
	double param_lengthUpper[NO][NO];

	/** lower bound of the parameters range */
	double param_lengthLower[NO][NO];

	/** current value of parameters */
	double param_length[NO][NO];

	/*
	 * Set of Parameters
	 * i.e. non-edges in an ActiveConstraintGraph that complete the graph into 3-tree.
	 *
	 * A Parameter is a pair of integer where each integer is the "index of the vertex in constraint graph".
	 * i.e. A Parameter is represented by <helix A "vertex" index , helix B "vertex" index>
	 */
	vector<pair<int, int> > parameters;

	/** Flag to show if sampling is completed or not*/
	bool gridDone;

	//----------JACOBIAN MEMBERS-------------
	/** up(1) or down(-1) from the mid point for 6 parameters */
	int direction[6];
	bool isBoundary[6]; //for jac sampling
	double mjpoint[6][6]; //middle point for jacobian sampling
	std::vector<int> independent_directions; //for jacobian sampling, independent rows on jacobian matrix
	double mj_traj_point[6][6]; //middle trajectory point for jacobian sampling

	// static const int pr = 42;//42; // 360/Settings::Sampling::gridSteps[3]; //number of partitions along 1 dimension
	// bool visited[pr][pr][pr][pr][pr];

	int start_coordinats[6];

	int end_coordinats[6];

	static const int pr1 = 54, pr2 = 54, pr3 = 16, pr4 = 38, pr5 = 5, pr6 = 38; // human20/rat20 molecule
	// static const int pr1=42, pr2=42, pr3=9, pr4=38, pr5=5, pr6=38;  //dim = 5 no grid Constraint //just set the size 'a little' (+-1) bigger than Cartesian grid boundaries to prevent rounding issues.
	// static const int pr1=60, pr2=60, pr3=60, pr4=38, pr5=4, pr6=38;   //dim = 2, 3
	//no need to set big values for rest of the pr s for low dimensions, once initialized less than pr would be enough. since it will not be increased anymore.

	bool visited[pr1][pr2][pr3][pr4][pr5][pr6];

	//----------JACOBIAN MEMBERS END-------------

	double stepSize;
	bool partial3tree;

	ActiveConstraintGraph* currentGraph;

	/*
	 * Set of Contacts
	 * A Contact is a pair of integer where each integer is the "index of the vertex in constraint graph".
	 * i.e. A Contact is represented by <helix A "vertex" index , helix B "vertex" index>
	 */
	std::vector<std::pair<int, int> > contacts;

private:

	/* The first NO/2 vertices represent helixA, The second NO/2 vertices represent helixB.
	 * lets say A = {0, ..., NO/2 -1} and  B = {NO/2, ..., NO-1}
	 * edge_length[A][A] represent internal distances of helixA
	 * edge_length[B][B] represent internal distances of helixB
	 * edge_length[A][B] is symmetric to edge_length[B][A] and represent distances between helixA and helixB
	 */
	double edge_length[NO][NO]; //fixed

	/**
	 * Vector of tetrahedrons. Each tetrahedron are represented by 4 vertex indices.
	 * Each sub vector hold the indices of vertices of one tetrahedron.
	 *
	 * @see CayleyParameterization::built3tree()
	 */
	std::vector<std::vector<int> > tetrahedra; //first is base tetrahedron, for the rest the first index will be computed by the last three indices.

	/** Step size during binary sampling, is halved repeatedly. */
	double stepSize_of_stepGridContact;

	/** The middle(witness) point where the sampling starts */
	double witnessPoint[6];  //vector<double> witnessPoint;

	/** Triangular inequality boundaries. Not tight as tetrahedral one. */
	double maxOfMax[6], minOfMin[6];

	/**
	 * Direction around a valid point in which we are searching for a contact boundary.
	 * Remember in order to find a new contact, we were walking in between a valid point and colliding point.
	 */
	int dir;

	/** binary: during binary search to find out a new contact, it determines how much we should half the step size. */
	int bin;

	/**
	 * adjacency map for dependency of parameters. It provides the set of
	 * parameters whose range will be updated when one parameter fixed.
	 */
	std::vector<std::vector<int> > updateList;

	/**
	 * inequalities that express the range of a parameter can be classified into
	 * either a linear or non-linear class. This variable is characterization
	 * of the parameter that tells what inequality is needed to compute parameter range.
	 * i.e. triangular or tetrahedral inequality.
	 *
	 * Set of integers corresponding to each parameter that takes values as:
	 * -2 if the parameter is actually a contact in case of short range sampling
	 * -1 if triangular inequality is used to compute the range of the parameter
	 * the index of the tetrahedron in tetra if tetrahedral inequality is used to compute the range of the parameter
	 */
	std::vector<int> boundaryComputationWay;

	/** Cayley point: the current location of sampling on Cayley space */
	double point[6];	//vector<double> point;

	PointSet* helA;
	PointSet* helB;
	PredefinedInteractions *df;

	/*
	 * vector of contact + parameter "atom" indices for the vertices of contact graph
	 * (atoms ,that are not in the actual contact list but hold for parameters, also exist in this list)
	 * the unique vertex list of the contact graph
	 */
	std::vector<int> verticesA;
	std::vector<int> verticesB;

	/** A flag for refined sampling */
	bool dense;

};

#endif /* CONVEXCHART_H_ */
