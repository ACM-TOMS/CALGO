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

#ifndef POINTSET_H_
#define POINTSET_H_

#ifdef WIN32
#include <windows.h>
#endif

#include <vector>
#include <list>

#include "Point.h"
#include "Utils.h"
#include "Stree.h"
#include "PredefinedInteractions.h"

/*
 * This class allows a MolecularUnit to be easily passed and manipulated.
 * It is primarily a sequence of with additional methods to
 * manipulate and get information about the list.
 */
class PointSet {
public:
	/////////////////////////////////
	// Constructors/Destructors
	/////////////////////////////////
	/** @brief Default constructor */
	PointSet();
	/** @brief Constructor with atom list initilization  */
	PointSet(std::vector<Point*> atoms);
	/** @brief Destructor */
	virtual ~PointSet();

	/////////////////////////////////
	// Points of the PointSet
	/////////////////////////////////

	/**
	 * @param a The pointer of an atom to add to the helix.
	 */
	void addAtom(Point *a);

	/**
	 * @return The list of atoms in the helix
	 */
	std::vector<Point*> getAtoms();

	/**
	 * @return True if helix contains no atoms, False otherwise.
	 */
	bool empty();

	/**
	 * @see addAtom
	 * @return The number of atoms added to the the helix.
	 */
	size_t size();

	/**
	 * @param a Pointer to the atom of interest
	 * @return True if helix contains a, False otherwise.
	 */
	bool contains(Point *a);

	/**
	 * @return The index'th atom added to the helix, or NULL if OOB
	 */
	Point* getAtomAt(size_t index);

	/**
	 * @param label The name of the atom of interest
	 * @return The pointer of the atom in the helix with name == label, or NULL
	 */
	Point* getAtomByLabel(std::string label);

	/**
	 * @see getAtomAt
	 * @param a The pointer of the atom of interest
	 * @return The index of atom a in the private list of atoms in the helix, if not in helix -0
	 */
	int getIndexOf(Point* a);

	/////////////////////////////////
	// Atom transformation
	/////////////////////////////////

	/**
	 * @brief This method applies a (supplied) transformation to the atoms of a helix a
	 * specified number of times and outputs a vector of those transformed atoms.
	 *
	 * Transformation is achieved by using fromB, toB and applied on helixB. HelixA is fixed.
	 *
	 * @return The vector of new atoms. Copies of all the original atoms with new locations.
	 * @param from For the transformation.
	 * @param to For the transformation.
	 * @param repeat The number of times to repeat the transformation
	 */
	std::vector<Point*> getXFAtoms(double from[][3], double to[][3],
			size_t repeat = 0);

	/**
	 * @brief Similar to getFXAtoms, but you supply the precalculated transformation
	 * matrix
	 *
	 * @see getXFAtoms
	 * @param trans_mat For the transformation.
	 */
	std::vector<Point*> getXFAtoms(const std::vector<double> &trans_mat,
			size_t repeat = 0);

	/**
	 * @brief Similar to getFXAtoms, but is done to only one atom once
	 * This method applies a transformation to a single atom of a helix and
	 * outputs a new transformed atom;
	 * The transformation applied needs to be supplied.
	 *
	 * @see getXFAtoms
	 * @return The new atom. Copy of the original atoms with a new location.
	 */
	Point* getXFAtom(size_t index, const std::vector<double> &trans_mat,
			size_t repeat = 0);

	/////////////////////////////////
	// Dumbbell methods
	/////////////////////////////////

	/**
	 * @return A vector of pairs that represent atoms in this helix
	 * The distance between atoms of each pair is in the range between RootNodeCreation::min and RootNodeCreation::max
	 * These will be used as candidates for the dumbbells
	 */
	std::vector<std::pair<int, int> > getDumbbells();

	/////////////////////////////////
	// Limits
	/////////////////////////////////

	/**
	 * @brief Call this to have the MolecularUnit determine its max and min xyz values
	 * by considering all atoms it already contains.
	 */
	void calcLimits();

	/**
	 * @brief Overwrites the given arrays with the max and min xyz values of the
	 * the helix atoms.
	 * This function alters the arrays given to it to reflect the maximum and minimum values for the xyz of the helix atoms.
	 * Used for the scale purpose of the display of nodeselection window
	 */
	void getLimits(double max_xyz[3], double min_xyz[3]);

	/////////////////////////////////
	// Other public methods
	/////////////////////////////////

	/**
	 * @brief Changes the location of all the atoms by subtracting off the
	 * average location. Does not change limits.
	 */
	void centerAtoms();

	/////////////////////////////////
	// Sphere Tree methods
	/////////////////////////////////

	/**
	 * @brief Builds a sphere tree using the atoms in the helix.
	 */
	void buildStree();

	/**
	 * @see buildStree
	 * @return The sphere tree based on the atoms in the helix (should call
	 * buildStree first).
	 */
	Stree* getStree();

	/** @brief Creates a molecular unit consists of only contact/interface atoms. */
	void simplify(PredefinedInteractions &distfinder, bool first_label);

	/////////////////////////////////
	// Print
	/////////////////////////////////

	/**
	 * @brief The ostream print method
	 * This operator allows the helix to be sent to the console or another filestream.
	 */
	friend std::ostream& operator<<(std::ostream& os, PointSet& h);

private:
	/*
	 * The collection of points in the PointSet.
	 */
	std::vector<Point*> points;

	/*
	 * The corners of an axis aligned volume containing all the point centers.
	 */
	double hi[3], lo[3];

public:
	Stree* stree;
	std::map<std::string, Point*> index;
};

#endif /*POINTSET_H*/
