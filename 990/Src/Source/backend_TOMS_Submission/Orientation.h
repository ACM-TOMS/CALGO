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

/*
 Set and get FromTo code is repeated over and over again! This should just
 call the one function.

 size_t used in place of unsigned int again.

 */

#ifndef ORIENTATION_H_
#define ORIENTATION_H_

#include <vector>
#include <cstdlib>

class ActiveConstraintGraph;
class PointSet;

#include <Eigen/Core>
#include <Eigen/LU>

/**
 * The Euclidean transformation of molecular composite
 */
class Orientation {
public:
	/////////////////////////////////
	// Constructors/Destructors
	/////////////////////////////////

	/** @brief Copy Constructor */
	Orientation(Orientation* copy);
	/** @brief Constructor with from/to initilization  */
	Orientation(double *fromB[], double *toB[]);
	/** @brief Constructor with form/to initilization  */
	Orientation(double fromB[][3], double toB[][3]);
	/** @brief Destructor, no code */
	virtual ~Orientation();

	/////////////////////////////////
	// From/To 3x3 matrices
	/////////////////////////////////

	/**
	 * @param fromB The "from" 3x3 matrix for the Orientation.
	 * @param toB The "to" 3x3 matrix for the Orientation.
	 */
	void setFromTo(double fromB[][3], double toB[][3]);

	/**
	 * @param fromB The "from" 3x3 matrix to be filled with Orientation's "from".
	 * @param toB The "to" 3x3 matrix to be filled with Orientation's "to".
	 */
	void getFromTo(double fromB[][3], double toB[][3]);

	/**
	 * @brief Compares the "from" and "to" matrices of two Orientation instances.
	 *
	 * @param other The other Orientation to compare against this.
	 * @return True if their matrices are the same, False otherwise
	 */
	bool isEqual(Orientation* other);

	/////////////////////////////////
	// Boundaries
	/////////////////////////////////

	/**
	 * @param new_boundary The new boundary to add to the current list of boundaries.
	 */
	void addBoundary(int new_boundary);

	/**
	 * @brief Empties the current list of boundaries then fills it with new ones.
	 *
	 * @param new_boundaries The vector of new boundaries.
	 */
	void setBoundary(std::vector<int> new_boundaries);

	/**
	 * @return The vector of boundaries.
	 */
	std::vector<int> getBoundary();

	/////////////////////////////////
	// Flips
	/////////////////////////////////

	/**
	 * @param flip The value to set the flipNum to.
	 */
	void setFlipNum(size_t flip);

	/**
	 * @return The flipNum.
	 */
	size_t getFlipNum();

	/** @brief Updates angle_violated status */
	void setAngleViolation(bool angleStatus);

	/** @return angle_violated status */
	bool angleViolated();

	Eigen::Vector3d TB;
	Eigen::Vector3d eaB;

private:

	/**
	 * The set of node indices that this orientation belongs.
	 * An orientation can take place in the boundary of multiple regions
	 */
	std::vector<int> boundary; //connections

	/**
	 * For 2 helices, there are 3 tetrahedrons, hence there are 8 ways to realize a Cayley point.
	 * flipNum determines which orientation of a realization to compute
	 */
	size_t flipNum;

	/** Cartesian coordinates of first three vertices of ActiveConstraintGraph in helixB */
	double fromB[3][3];

	/** Transformed Cartesian coordinates of first three vertices of ActiveConstraintGraph */
	double toB[3][3];

	/** True if the angle between two helices is out of the allowed range, False otherwise */
	bool angle_violated;

};

#endif /* ORIENTATION_H_ */
