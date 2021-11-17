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
#ifndef PREDEFINEDINTERACTIONS_H_
#define PREDEFINEDINTERACTIONS_H_

#include <string>
#include <vector>
#include <list>
#include <map>
#include <sstream>
#include <iostream>

#include "Point.h"

#include <limits>

/**
 * This class is to hold a certain distance for the pair of atoms in order to define contact/bonding distance,
 * rather than just sum of the radii of the atoms
 */
class PredefinedInteractions {
public:
	/////////////////////////////////
	// Constructors/Destructors
	/////////////////////////////////

	/** @brief Default constructor */
	PredefinedInteractions();

	/**
	 * @brief Populates dist1 where each instance holds labelA, labelB and bonding distance
	 *
	 * @param values Set of rows where each row holds the labels of the two atoms and their bonding distance
	 * @param labelAcol The column index of the label of first atom
	 * @param labelBcol The column index of the label of second atom
	 * @param radiuscol The column index of the bonding distance between two atoms
	 */
	PredefinedInteractions(std::vector<std::vector<std::string> > values,
			size_t labelAcol, size_t labelBcol, size_t radiuscol);

	/**
	 * @brief Populates dist1 where each instance holds labelA, labelB and lower bonding distance
	 * 				and dist2 where each instance holds labelA, labelB and upper bonding distance
	 *
	 * @param values Set of rows where each row holds the labels of the two atoms and their lower and upper bonding distance
	 * @param labelAcol The column index of the label of first atom
	 * @param labelBcol The column index of the label of second atom
	 * @param radiusmincol The column index of the lower bonding distance between two atoms
	 * @param radiusmaxcol The column index of the upper bonding distance between two atoms
	 */
	PredefinedInteractions(std::vector<std::vector<std::string> > values,
			size_t labelAcol, size_t labelBcol, size_t radiusmincol,
			size_t radiusmaxcol);

	/**
	 * @brief Populates dist1 where each instance holds labelA, labelB and bonding distance
	 =	 *
	 * @param values The data that is in matrix format. Labels are stored in the first row and first column.
	 * The bonding distance between atom with ith label and the atom with jth label is stored in i-j th cell of the matrix.
	 */
	PredefinedInteractions(std::vector<std::vector<std::string> > values);

	/** @brief Destructor */
	virtual ~PredefinedInteractions();

	/////////////////////////////////
	// Other public methods
	/////////////////////////////////

	/** @return Sum of radii of atoms atomA and atomB */
	double sumOfRadii(Point * atomA, Point * atomB);

	/** @return Lower bound distance for collision between atoms atomA and atomB */
	double collisionLowerBound(Point * atomA, Point * atomB);

	/**
	 * @return Lower bound distance for bonding between atoms atomA and atomB
	 * if dist2Defined then dist1 is the lower bound distance
	 * if not, then dist1 is the exact value. Hence bondLower() equation needed to be computed.
	 *
	 * If Settings::General::candidate_interactions set to true
	 * Interactions between arbitrary atoms (atoms that are not part of candidate_interactions) does not lead to a bond i.e. not create a child node
	 * Output infinity forces the input pair not to be a contact.
	 */
	double bondingLowerBound(Point * atomA, Point * atomB);

	/**
	 * @return Upper bound distance for bonding between atoms atomA and atomB
	 * if dist2Defined then dist2 is the upper bound distance
	 * if not, then dist1 is the exact value. Hence bondUpper() equation needed to be computed.
	 *
	 * If Settings::General::candidate_interactions set to true
	 * Interactions between arbitrary atoms (atoms that are not part of candidate_interactions) does not lead to a bond i.e. not create a child node
	 * Output -infinity forces the input pair not to be a contact.
	 */
	double bondingUpperBound(Point * atomA, Point * atomB);

	/**
	 * @param atomA, atomB are pointers to the atoms
	 * @return The min(or exact) distance between two atoms which is specified in
	 * distData window.
	 * If the pair of two atoms are not in the dist1 list, it will return 0.
	 */
	double getDist1(Point * atomA, Point * atomB);

	/**
	 * @param atomA, atomB are pointers to the atoms
	 * @return The upper bound distance between two atoms which is specified in
	 * distData window.
	 * If the pair of two atoms are not in the dist2 list or if dist2 is not defined, it will return 0.
	 * But it will not be an obstacle for checking a contact, dist1 can be used to compute upper bound in that case.
	 */
	double getDist2(Point * atomA, Point * atomB);

	/** @brief Assigns other to this PredefinedInteractions */
	void assign(PredefinedInteractions *other);
	bool empty();

	/////////////////////////////////
	// Iterator methods
	/////////////////////////////////
	typedef std::map<std::pair<std::string, std::string>, double>::iterator dist_iterator;
	dist_iterator dist1begin() {
		return dist1.begin();
	}
	dist_iterator dist1end() {
		return dist1.end();
	}
	dist_iterator dist2begin() {
		return dist2.begin();
	}
	dist_iterator dist2end() {
		return dist2.end();
	}

	/////////////////////////////////
	// Helper methods
	/////////////////////////////////
private:
	/** @return The modified distance with user defined threshold cooperating lambda/delta constants for collision check */
	double collisionLower(double distance);

	/** @return The modified distance with user defined threshold cooperating lambda/delta constants for lower distance of bonding */
	double bondLower(double distance);

	/** @return The modified distance with user defined threshold cooperating lambda/delta constants for upper distance of bonding */
	double bondUpper(double distance);

	/** @return True if val takes place in maap, False otherwise */
	bool find(std::map<std::pair<std::string, std::string>, double> maap,
			std::pair<std::string, std::string> val);

private:

	/**
	 * True if dist1 is populated, False otherwise
	 * Settings::General::candidate_interactions is true iff PredefinedInteractions::distTableDefined is true
	 */
	bool distTableDefined;

	/** True if dist2 is populated, False otherwise */
	bool dist2Defined;

	/**
	 * Set of instances where each instance holds labelA, labelB and lower(or exact)bonding distance
	 * If dist2Defined then dist1 holds the lower value, otherwise dist1 holds the exact value.
	 */
	std::map<std::pair<std::string, std::string>, double> dist1;

	/** Set of instances where each instance holds labelA, labelB and upper bonding distance */
	std::map<std::pair<std::string, std::string>, double> dist2;
};

#endif /* PREDEFINEDINTERACTIONS_H_ */
