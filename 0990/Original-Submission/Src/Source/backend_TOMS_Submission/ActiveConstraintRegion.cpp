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
 * ActiveConstraintRegion.cpp
 *
 *  Created on: 2014
 *      Author: Aysegul Ozkan
 */

#include "ActiveConstraintRegion.h"

#include "Point.h"
#include "Orientation.h"
#include "Settings.h"

#include <cmath>
#include <math.h>
#include <algorithm>
#include <exception>
#include <string>

using namespace std;

#define PI 3.14159265

ActiveConstraintRegion::ActiveConstraintRegion() {
}

ActiveConstraintRegion::ActiveConstraintRegion(
		const ActiveConstraintRegion &acr) {

	//Deep copy everything
	vector<CayleyPoint*> witspc = acr.getWitness();
	for (vector<CayleyPoint*>::iterator iter = witspc.begin();
			iter != witspc.end(); iter++) {
		CayleyPoint *tmp = new CayleyPoint(**iter);
		witspace.push_back(tmp);
	}

	vector<CayleyPoint*> spc = acr.getJustSpace();
	for (vector<CayleyPoint*>::iterator iter = spc.begin(); iter != spc.end();
			iter++) {
		CayleyPoint *tmp = new CayleyPoint(**iter);
		space.push_back(tmp);
	}

}

ActiveConstraintRegion::~ActiveConstraintRegion() {

	trim();
}

void ActiveConstraintRegion::trim() {
	for (vector<CayleyPoint*>::iterator iter = this->space.begin();
			iter != this->space.end(); iter++) {
		try {
			(*iter)->trim_PointMultiD();
			delete (*iter);
		} catch (exception & e) {
			cout << e.what() << endl;
		}
	}
	this->space.clear();

	for (vector<CayleyPoint*>::iterator iter = this->witspace.begin();
			iter != this->witspace.end(); iter++) {
		try {
			(*iter)->trim_PointMultiD(); //this line in unnecessary, does the same job with the following line.
			delete (*iter);
		} catch (exception & e) {
			cout << e.what() << endl;
		}
	}
	this->witspace.clear();
}

void ActiveConstraintRegion::setSpaceVolume(double min[6], double max[6]) {
	if (this->volume.second[0] == -1) { // it isn't already set.
		for (size_t x = 0; x < 6; x++) {
			this->volume.first[x] = min[x];
			this->volume.second[x] = max[x];
		}
	}
}

void ActiveConstraintRegion::getSpaceVolume(double min[6], double max[6]) {
	for (size_t x = 0; x < 6; x++) {
		min[x] = this->volume.first[x];
		max[x] = this->volume.second[x];
	}
}

vector<CayleyPoint*> ActiveConstraintRegion::convertSpace(
		ActiveConstraintRegion* other, vector<pair<int, int> > paramlines,
		PointSet *helA, PointSet *helB) {
	vector<CayleyPoint*> output;
	vector<CayleyPoint*> input = other->getSpace();
	for (size_t i = 0; i < input.size(); i++) {
		if (input[i]->hasOrientation()) {
			vector<Orientation*> ors = input[i]->getOrientations();
			for (size_t j = 0; j < ors.size(); j++) {
				CayleyPoint* pointt = new CayleyPoint(ors[j], paramlines, helA,
						helB);
				output.push_back(pointt);
			}
		}
	}

	return output;
}

vector<CayleyPoint*> ActiveConstraintRegion::getSpace() //todo ??make it efficient, keep another variable composite of both spaces, create it once and call pointer of it many.
{
	vector<CayleyPoint*> output;
	output.assign(this->witspace.begin(), this->witspace.end());
	output.insert(output.end(), this->space.begin(), this->space.end());
	return output;
}

void ActiveConstraintRegion::setWitSpace(vector<CayleyPoint*> input) {
	this->witspace.assign(input.begin(), input.end());
}

void ActiveConstraintRegion::setJustSpace(vector<CayleyPoint*> input) {
	this->space.assign(input.begin(), input.end());
}

vector<CayleyPoint*> ActiveConstraintRegion::getJustSpace() const {
	return this->space;
}

vector<CayleyPoint*> ActiveConstraintRegion::getWitness() const {
	return this->witspace;
}

int ActiveConstraintRegion::getSpaceSize() {
	return this->space.size();
}

int ActiveConstraintRegion::getCombinedSize() {
	return (this->space.size() + this->witspace.size());
}

int ActiveConstraintRegion::getWitnessSize() {
	return this->witspace.size();
}

void ActiveConstraintRegion::insertSpace(CayleyPoint* point) {
	this->space.push_back(point);
}

void ActiveConstraintRegion::insertWitness(CayleyPoint* point) {
	this->witspace.push_back(point);
}

