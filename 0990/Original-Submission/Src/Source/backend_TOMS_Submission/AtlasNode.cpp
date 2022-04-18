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
 * roadNode.cpp
 *
 *  Created on: Feb 2, 2009
 *      Author: James Pence
 */

#include "AtlasNode.h"

#include "Utils.h"

#include <algorithm>
#include <cmath>
using namespace std;



AtlasNode::AtlasNode() {

	this->visited = false;

	for (int i = 0; i < 3; i++) {
		this->loc[i] = 0;
		this->vel[i] = 0;
		this->force[i] = 0;
	}
	this->velMag = 15.0;
	this->forceMag = 1.0;

	this->numID = -1;
	this->numPass = 1;
	this->dim = 0;

	this->complete = false;
	this->noGoodOrientation = true;

	this->constraintGraph = new ActiveConstraintGraph(); // NULL;
	this->region = new ActiveConstraintRegion();

	this->dimWritten = false;
}

AtlasNode::AtlasNode(int ID, bool complete, bool empty, int numDim,
	double* location, vector<int> connection) {


	this->visited = false;

	for (int i = 0; i < 3; i++) {
		this->loc[i] = location[i];
		this->vel[i] = 0;
		this->force[i] = 0;
	}
	this->velMag = 15.0;
	this->forceMag = 1.0;

	this->numID = ID;
	this->numPass = 1;
	this->dim = numDim;

	this->complete = complete;

	if (numDim == 0 && !Settings::Sampling::short_range_sampling) // if it is dimension 0, there is no parameter to sample except 6d. // and virus case. || Settings::General::virus_mer
		this->complete = true;

	this->noGoodOrientation = empty; //it would be better to initialize it to false, by this way the node will be displayed during sample and if there is no good orientation in it after search THEN it will be set to TRUE then it will disappear.

	this->connection.assign(connection.begin(), connection.end());
	this->constraintGraph = NULL;
	this->region = new ActiveConstraintRegion(); // ??????? NULL;

	this->dimWritten = false;
}

AtlasNode::~AtlasNode() {
	this->connection.clear();
	delete this->constraintGraph;
	delete this->region;
}

void AtlasNode::setLocation(double x, double y, double z) {
	this->loc[0] = x;
	this->loc[1] = y;
	this->loc[2] = z;
}

double* AtlasNode::getLocation() {
	return this->loc;
}

double* AtlasNode::getVelocity() {
	return this->vel;
}

double* AtlasNode::getForce() {
	return this->force;
}

double AtlasNode::getForceMag() {
	return this->forceMag;
}

vector<int> AtlasNode::getConnection() {
	return this->connection;
}

void AtlasNode::applyForce() {
	this->forceMag = Utils::mag(this->force);

	for (int i = 0; i < 3; i++) {
		if (i == 2) {
			this->vel[i] += this->force[i] * .01;
		} else {
			this->vel[i] += this->force[i] * .04;
		}
	}
	for (int i = 0; i < 3; i++) {
		this->vel[i] *= .4;
	}

}

void AtlasNode::applyVelocity() {
	for (int i = 0; i < 3; i++) {
		if (fabs(this->vel[i]) > 30) {
			if (this->vel[i] < 0) {
				this->vel[i] = -30;
			} else {
				this->vel[i] = 30;
			}
		}
		this->loc[i] += this->vel[i];
	}
	this->velMag = Utils::mag(this->vel);
}

void AtlasNode::addConnection(int other) {
	for (size_t i = 0; i < this->connection.size(); i++) {
		if (this->connection[i] == other) {
			return;
		}
	}

	this->connection.push_back(other);
}

bool AtlasNode::removeConnection(int other) {
	vector<int>::iterator iter;
	iter = find(this->connection.begin(), this->connection.end(), other);
	if (iter == this->connection.end()) {
		return false;
	} else {
		this->connection.erase(iter);
		return true;
	}
}

bool AtlasNode::isConnectedTo(int other) {
	for (size_t i = 0; i < this->connection.size(); i++) {
		if (this->connection[i] == other) {
			return true;
		}
	}
	return false;
}

void AtlasNode::setComplete(bool setting) {
	this->complete = setting;
}

bool AtlasNode::hasAnyGoodOrientation() {
	return !noGoodOrientation;
}

void AtlasNode::setFoundGoodOrientation(bool setting) {
	noGoodOrientation = !setting;
}

bool AtlasNode::isComplete() {
	return this->complete;
}

int AtlasNode::getID() {
	return this->numID;
}

void AtlasNode::setID(int id) {
	this->numID = id;
}

int AtlasNode::getDim() {
	return this->dim;
}

int AtlasNode::getParamDim() {
	if (this->constraintGraph != NULL) {
		return this->constraintGraph->getParamDim();
	}
	return this->dim;
}

ActiveConstraintRegion* AtlasNode::getACR() {
	return this->region;
}

void AtlasNode::setACR(ActiveConstraintRegion* newRegion) {
	if (this->region != NULL) {
		delete this->region;
	}
	this->region = newRegion;
}

ActiveConstraintGraph* AtlasNode::getCG() {
	return this->constraintGraph;
}

void AtlasNode::setCG(ActiveConstraintGraph* acg) {
	if (this->constraintGraph != NULL) {
		delete this->constraintGraph;
	}
	this->constraintGraph = acg;
	this->dim = acg->getDim();
}

void AtlasNode::trimNode() {
	delete this->region;
	this->region = NULL;

	delete this->constraintGraph;
	this->constraintGraph = NULL;
}
