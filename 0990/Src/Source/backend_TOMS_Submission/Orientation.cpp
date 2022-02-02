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

#include "Orientation.h"
#include "Utils.h"
#include <vector>
#include <cassert>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/Dense"
#include "Eigen/SVD"

using namespace std;
using namespace Eigen;

using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::MatrixXd;
using Eigen::VectorXd;

Orientation::Orientation(Orientation* copy) {
	double fromB[3][3];
	double toB[3][3];
	copy->getFromTo(fromB, toB);

	this->setFromTo(fromB, toB);

	this->boundary = copy->getBoundary();
	this->angle_violated = copy->angle_violated;

	this->flipNum = copy->flipNum;

	this->TB = copy->TB;
	this->eaB = copy->eaB;
}

Orientation::Orientation(double *fromB[3], double *toB[3]) {
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			this->fromB[i][j] = fromB[i][j];
			this->toB[i][j] = toB[i][j];
		}

	this->angle_violated = false;
}

Orientation::Orientation(double fromB[][3], double toB[][3]) {
	this->setFromTo(fromB, toB);
	this->angle_violated = false;
}

Orientation::~Orientation() {
}

bool Orientation::isEqual(Orientation* other) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (this->fromB[i][j] != other->fromB[i][j]
					|| this->toB[i][j] != other->toB[i][j]) {
				return false;
			}
		}
	}

	return true;
}

void Orientation::setFromTo(double fromB[][3], double toB[][3]) {

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			this->fromB[i][j] = fromB[i][j];
			this->toB[i][j] = toB[i][j];
		}

}

void Orientation::getFromTo(double fromB[][3], double toB[][3]) {

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			fromB[i][j] = this->fromB[i][j];
			toB[i][j] = this->toB[i][j];
		}
}

void Orientation::addBoundary(int new_boundary) {
	if (new_boundary >= 0)  // filters out all negative boundaries
		this->boundary.push_back(new_boundary);
}

void Orientation::setBoundary(vector<int> new_boundaries) {
	this->boundary.clear();
	for (int i = 0; i < new_boundaries.size(); i++)
		addBoundary(new_boundaries[i]);
}

vector<int> Orientation::getBoundary() {
	return this->boundary;
}

void Orientation::setFlipNum(size_t flip) {
	this->flipNum = flip;
}

size_t Orientation::getFlipNum() {
	return this->flipNum;
}

void Orientation::setAngleViolation(bool angleStatus) {
	this->angle_violated = angleStatus;
}

bool Orientation::angleViolated() {
	return this->angle_violated;
}

