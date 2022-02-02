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
#ifndef ATOM_H_
#define ATOM_H_

#include <iostream>
#include <vector>
#include <string>

/**
 * This class holds the location and radius of each atom which goes into the
 * helix.
 */

class Point {
public:
	/////////////////////////////////
	// Constructors/Destructors
	/////////////////////////////////

	/** @brief Default constructor */
	Point();

	/** @brief Copy constructor */
	Point(Point* copy);

	/** @brief Constructor with location and radius initialization */
	Point(double location[], double radius);

	/** @brief Constructor with location, radius, name, and atomID initialization */
	Point(double location[], double radius, std::string name,
			std::string atomID);

	/** @brief Destructor, no code */
	virtual ~Point();

	/////////////////////////////////
	// Setters
	/////////////////////////////////

	/**
	 * @brief A setter for the atom's 3D coordinates
	 *
	 * @param x The x coordinate to give the atom.
	 * @param y The y coordinate to give the atom.
	 * @param z The z coordinate to give the atom.
	 */
	void setLocation(double x, double y, double z);

	/**
	 * @param radius The radius to give the atom
	 */
	void setRadius(double radius);

	/**
	 * @param name An identifying name for the atom.
	 */
	void setName(std::string name);

	/**
	 * @param atomID An identifying number (string) for the atom
	 */
	void setPointID(std::string atomID);

	/////////////////////////////////
	// Getters
	/////////////////////////////////

	/**
	 * @brief A getter for the atom's 3D coordinates of the location
	 *
	 * @return A pointer to a double array of length 3 (x,y,z coords)
	 */
	double* getLocation();

	/**
	 * @return The atom's radius
	 */
	double getRadius();

	/**
	 * @return The atom's identifying name
	 */
	std::string getName();

	/**
	 * @return The atom's identifying number (string)
	 */
	std::string getPointID();

	/////////////////////////////////
	// Print
	/////////////////////////////////

	/**
	 * @brief The ostream print method
	 */
	friend std::ostream & operator<<(std::ostream & os, Point & a);

	/////////////////////////////////
	// Other public methods
	/////////////////////////////////

	/**
	 * @brief Saves a line of "data" that was associated with the atom. Useful
	 * if the file came with extra info you want to hold. X,y,z col refer to the
	 * placement of the x,y,z data (linedata[xcol], etc.).
	 */
	void setLine(std::vector<std::string> linedata, size_t xcol, size_t ycol,
			size_t zcol);

	/**
	 * @brief Convert the "data" that was associated with the atom to the pdb_format style
	 * @return A prettily formatted string representation of the line.
	 */
	std::string getLine();

	/**
	 * @see setLine
	 * @return The class variable row, filled with setLine
	 */
	std::vector<std::string> getRow();

private:

	/** Atom's 3D coordinates of the location */
	double loc[3];

	/** radius of the atom */
	double radius;

	std::string name;

	/** The atom's identifying number (string) */
	std::string pointID;

	/** A line of "data" that was associated with the atom. */
	std::vector<std::string> row;

	/** Refer to the placement of the x,y,z data (row[xcol], etc.). */
	size_t xCol, yCol, zCol;
};

#endif /*ATOM_H_*/
