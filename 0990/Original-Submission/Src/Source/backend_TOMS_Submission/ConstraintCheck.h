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
 * ConstraintCheck.h
 *
 *  Created on: Oct 20, 2014
 *      Author: Aysegul Ozkan
 */

#ifndef CONSTRAINTCHECK_H_
#define CONSTRAINTCHECK_H_

#include "ActiveConstraintGraph.h"
#include "Orientation.h"

#include <map>
#include <list>
#include <vector>
#include <utility>

/*
 * This class get delegated the duty of determining if the helices intersect.
 * It is designed to check whether any non-active constraints become active or whether they are violated.
 *
 * User have the option to define set of constraints of interest. In that case, new
 * constraint activation check will be done only for this set. For an input Orientation,
 * ConstraintCheck first computes the Cartesian realization for the entire molecular
 * composite then passes it to the subroutines to do user specified constraint check
 * such as steric constraints, tethering constraints or angle constraints.
 */
class ConstraintCheck {
public:

	/** @brief Constructor with ActiveConstraintGraph and PredefinedInteractions initialization */
	ConstraintCheck(ActiveConstraintGraph* cgk, PredefinedInteractions *df);

	/** @brief Destructor */
	virtual ~ConstraintCheck();

	/** @return True if the atoms a and b intersects more than the allowed threshold, False otherwise */
	bool stericViolated(Point* a, Point* b);

	/**
	 * @brief Checks if the atoms a and b are interacting
	 * i.e. the distance between a and b is within the bounding range
	 * If so adds to contactMap
	 *
	 * @param index The index of the atom in the molecular unit
	 *
	 * todo find a better name ifContactAddToList ?
	 */
	void updatePossibleContactList(Point* a, Point* b, int indexA, int indexB);

	/**
	 * @brief Computed the transformation between helix A and helix B
	 * Assumes initially helixA and helixB were aligned paralel
	 * Hence instead of computing the angle between fixed helixA and oriented helixB,
	 * computing the angle between initial helixB and oriented helixB, would do the same job.
	 *
	 * @param fb, tb Determines transformation matrix from original helixB to the oriented helixB
	 *
	 * @return True if the angle between two helices is out of the allowed range, False otherwise
	 */
	bool angleViolationCheck(double fb[3][3], double tb[3][3]);

	/*
	 * @brief This function checks if steric constraint is violated for atom pairs between helixA and helixB
	 * Updates the angle flag of the orientation
	 * Populates the candidate contact list.
	 *
	 * @return True if helix A and helix B overlaps, False otherwise
	 */
	bool stericsViolated(Orientation* ort);

	/*
	 * @return The list of Atom pairs of possible contacts from the most recent
	 * helix collision test for corresponding flip/orientation.
	 * The contacts are sorted by distance between them.
	 */
	std::list<std::pair<int, int> > getContacts(int flipno);

	bool outerGrid; //for Marias case: to check if there is a pair close enough or not

	static std::vector<double> nei_matrix;	// transform matrix of the neibour

private:
	int flip;

	//** The contact list for each flip where contacts are sorted by proximity */
	std::map<double, std::pair<int, int> > contactMap[8];

	ActiveConstraintGraph* currentGraph;
	PointSet *helA, *helB;
	PredefinedInteractions *df;
};

#endif /* CONSTRAINTCHECK_H_ */
