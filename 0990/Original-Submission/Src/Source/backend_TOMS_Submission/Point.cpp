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
#include "Point.h"

#include <stdlib.h>
#include <string>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

using namespace std;

Point::Point() {
	for (int i = 0; i < 3; i++)
		this->loc[i] = 0.0f;
	this->radius = 0.0f;
}

Point::Point(Point* copy) {
	for (int i = 0; i < 3; i++) {
		this->loc[i] = copy->loc[i];
	}
	this->radius = copy->radius;

	this->name = copy->name;
	this->pointID = copy->pointID;
	this->row = copy->row;
	this->xCol = copy->xCol;
	this->yCol = copy->yCol;
	this->zCol = copy->zCol;
}

Point::Point(double location[], double radius) {
	for (int i = 0; i < 3; i++)
		this->loc[i] = location[i];
	this->radius = radius;
	this->name = "";
	this->pointID = "";
}

Point::Point(double location[], double radius, string name, string pointID) {
	for (int i = 0; i < 3; i++)
		this->loc[i] = location[i];
	this->radius = radius;
	this->name = name;
	this->pointID = pointID;
}

Point::~Point() {

}

void Point::setLocation(double x, double y, double z) {
	this->loc[0] = x;
	this->loc[1] = y;
	this->loc[2] = z;
}

void Point::setRadius(double radius) {
	this->radius = radius;
}

void Point::setName(string name) {
	this->name = name;
}

void Point::setPointID(string pointID) {
	this->pointID = pointID;
}

double Point::getRadius() {
	return this->radius;
}

string Point::getName() {
	return this->name;
}

vector<string> Point::getRow() {
	return this->row;
}

string Point::getPointID() {
	return this->pointID;
}

void Point::setLine(vector<string> linedata, size_t xcol, size_t ycol,
		size_t zcol) {
	this->row.assign(linedata.begin(), linedata.end());
	this->xCol = xcol;
	this->yCol = ycol;
	this->zCol = zcol;
}

string Point::getLine() {
	stringstream ss;

	bool pdb_format = false;

	if (!pdb_format) {
		for (size_t i = 0; i < this->row.size(); i++) {
			if (i == this->xCol)
				ss << this->loc[0] << "\t";
			else if (i == this->yCol)
				ss << this->loc[1] << "\t";
			else if (i == this->zCol)
				ss << this->loc[2] << "\t";
			else
				ss << this->row[i] << "\t";
		}

	} else {
		ss << setw(6) << left << this->row[0];
		ss << setw(5) << right << this->row[1];
		ss << "  ";
		ss << setw(3) << left << this->row[2];
		ss << setw(4) << right << this->row[3];
		ss << setw(2) << right << this->row[4];
		ss << setw(4) << right << this->row[5];
		ss << setw(4) << "";
		ss << setw(8) << setiosflags(ios::fixed) << setprecision(3) << right
				<< this->loc[0];
		ss << setw(8) << setiosflags(ios::fixed) << setprecision(3) << right
				<< this->loc[1];
		ss << setw(8) << setiosflags(ios::fixed) << setprecision(3) << right
				<< this->loc[2];
		ss << setw(6) << right << this->row[9];
		ss << setw(6) << right << this->row[10];
		ss << setw(6) << right << this->row[11];

	}

	ss.flush();
	return ss.str();
}

double* Point::getLocation() {
	return this->loc;
}

ostream & operator<<(ostream & os, Point & a) {
	os << "x= " << a.loc[0];
	os << "\ty= " << a.loc[1];
	os << "\tz= " << a.loc[2];
	os << "\name= " << a.name;
	os << "\trad= " << a.radius << endl;
	return os;
}
