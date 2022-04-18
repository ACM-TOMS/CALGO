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
 * CayleyPoint.h
 *
 *  Created on: Feb 22, 2010
 *      Author: James Pence
 */

/*
 Are we using c++11? If we are, we can use delegating constructors:
 PointMultiD::PointMultiD(size_t size, double* values) : PointMultiD::PointMultiD() {...}
 PointMultiD::PointMultiD(double* values) : PointMultiD::PointMultiD(6) {}
 etc.
 Otherwise, code could really be cut down with a private init function. I just
 went ahead and added one.
 */

#ifndef CAYLEYPOINT_H_
#define CAYLEYPOINT_H_

#include <vector>
#include <list>
#include <cstdlib>
#include <cmath>
#include <algorithm>

// #include "Utils.h"
#include "Orientation.h"
#include "PointSet.h"

/**
 * This class represents a multi-dimension point in the Cayley parameter space, and stores the
 * corresponding Cartesian space orientations of the molecular composite.
 */
class CayleyPoint {
public:

	/////////////////////////////////
	// Constructors/Destructors
	/////////////////////////////////

	/** @brief Default constructor */
	CayleyPoint();

	/** @brief Like default constructor, but initilizes the point values */
	CayleyPoint(std::vector<double> values);

	/** @brief Copy constructor */
	CayleyPoint(const CayleyPoint &p);

	/**
	 * @brief Constructor that defines a point with one orientation originating
	 * from another configuration space.
	 * The data values are found from the orientation orient.
	 *
	 * @param orient The orientation that is duplicated and duplicated version is saved to this orient list.
	 * @param paramlines Set of parameter pairs each represented by <helA atom index , helB atom index>
	 */
	CayleyPoint(Orientation* orient,
			std::vector<std::pair<int, int> > paramlines, PointSet* helA,
			PointSet* helB);

	/** @brief Deletes all associated orientations. */
	virtual ~CayleyPoint();

	/////////////////////////////////
	// Other
	/////////////////////////////////

	/**
	 * @param setting The Realizable value for the instance to take.
	 * Should be False if volume is negative, True otherwise.
	 */
	void setRealizable(bool setting);

	/**
	 * @return The Realizable value of the instance. Should be False if volume is
	 * negative, True otherwise.
	 */
	bool isRealizable();

	/**
	 * @return The dimensions of the point, AKA the size of the data vector.
	 */
	size_t dim();

	/**
	 * @param out An array that gets filled with the values of the first 6
	 * dimensions of the instance. If the dimensions of the point are < 6, then
	 * the remaining values are filled with 0.
	 */
	void getPoint(double out[6]);

	/////////////////////////////////
	// Orientation methods
	/////////////////////////////////

	/**
	 * @param orient If it's not NULL and has not already been added, then it
	 * gets added to the instances vector of Orientations. Otherwise it gets
	 * deleted.
	 */
	void addOrientation(Orientation* orient);

	/**
	 * @return The vector of added orientations.
	 */
	std::vector<Orientation*> getOrientations();

	/**
	 * @brief It replaces the instance's vector of Orientations with the new one.
	 *
	 * @param oris The replacement orientations.
	 */
	void setOrientations(std::vector<Orientation*> oris);

	/**
	 * @param ori The Orientation to remove from the point..
	 * @return True if it was found (and removed), False otherwise.
	 */
	bool removeOrientation(Orientation* ori);

	/**
	 * @brief Deletes all of the orientations curretly associated with the instance.
	 * @details [long description]
	 */
	void trim_PointMultiD();

	/**
	 * @return True if Orientations have been added (and have not yet all been
	 * removed), False otherwise.
	 */
	bool hasOrientation();

	/**
	 * @return The boundaries of all of the instance's Orientations combined in
	 * one vector.
	 */
	std::list<int> getBoundaries();

	/**
	 * @brief Prints the data as well as realizable and badAngleN to cout.
	 * @details [long description]
	 */
	void printData();

	/////////////////////////////////
	// Operators
	/////////////////////////////////

	/**
	 * @brief Access a value in the data vector.
	 *
	 * @param index The index of the entry of interest
	 * @return If index < size of vector, the value there is returned. Otherwise,
	 * it returns 0.
	 */
	double operator[](size_t index);

	/////////////////////////////////
	// Numbers
	/////////////////////////////////

	/**
	 * There are set methods for things like badAngleN, collidN
	 * Since the orientations with bad angle or collisions are not saved in the orientation list.
	 * These numbers cannot be calculated in this class.
	 */
	void incrementBadAngleN();
	void incrementCollidN();

	void setBadAngleN(int n);
	void setCollidN(int n);

	int getBadAngleN();
	int getCollidN();

	/////////////////////////////////
	// Public variables
	/////////////////////////////////

	int axis; //for jacobian sampling
	int zIndex;

private:
	/**
	 * @brief Initializes member variables. Called by all constructors.
	 */
	void init_();

	/** value of the Cayley parameters (non-edge lengths) */
	std::vector<double> data;

	/** False if volume is negative, True otherwise. */
	bool realizable;

	/**
	 * set of Cartesian space Orientations of the molecular composite that are
	 * computed by realizing the active constraint graph with given length of edge/non-edges.
	 */
	std::vector<Orientation*> orient;

	/** Number of orientations where the angle between two helices is out of the allowed range */
	int badAngleN;

	/** number of orientations with collisions */
	int collidN;

};

#endif /* CAYLEYPOINT_H_ */
