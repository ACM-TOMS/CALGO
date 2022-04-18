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
 * ConstraintCheck.cpp
 *
 *  Created on: Oct 20, 2014
 *      Author: Aysegul Ozkan
 */

#include "ConstraintCheck.h"

#include "Settings.h"
#include "Point.h"
#include "PointSet.h"
#include "CartesianRealizer.h"
#include "Settings.h"

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/Dense"
#include "Eigen/SVD"
using namespace Eigen;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <fstream>
using namespace std;

#define PI 3.14159265

vector<double> ConstraintCheck::nei_matrix;

ConstraintCheck::ConstraintCheck(ActiveConstraintGraph* cgk,
		PredefinedInteractions *df) {

	this->outerGrid = false;
	this->currentGraph = cgk;

	this->helA = this->currentGraph->getMolecularUnitA();
	this->helB = this->currentGraph->getMolecularUnitB();
	this->df = df;

}

ConstraintCheck::~ConstraintCheck() {
}

bool ConstraintCheck::stericViolated(Point* a, Point* b) {

	double* aLoc = a->getLocation();
	double* bLoc = b->getLocation();
	double sumOfRadius = df->sumOfRadii(a, b);
	double originDist = sqrt(
			((aLoc[0] - bLoc[0]) * (aLoc[0] - bLoc[0]))
					+ ((aLoc[1] - bLoc[1]) * (aLoc[1] - bLoc[1]))
					+ ((aLoc[2] - bLoc[2]) * (aLoc[2] - bLoc[2])));

	double collisionLowerBound = df->collisionLowerBound(a, b);

	if (originDist < 0.9 + sumOfRadius)
		this->outerGrid = false;

	return originDist - collisionLowerBound < 0;

}

void ConstraintCheck::updatePossibleContactList(Point* a, Point* b, int indexA,
		int indexB) {
	double* aLoc = a->getLocation();
	double* bLoc = b->getLocation();
	double originDist = sqrt(
			((aLoc[0] - bLoc[0]) * (aLoc[0] - bLoc[0]))
					+ ((aLoc[1] - bLoc[1]) * (aLoc[1] - bLoc[1]))
					+ ((aLoc[2] - bLoc[2]) * (aLoc[2] - bLoc[2])));

	double bondingLowerBound = df->bondingLowerBound(a, b);
	double bondingUpperBound = df->bondingUpperBound(a, b);
	double middle = (bondingLowerBound + bondingUpperBound) / 2;

	if (bondingLowerBound <= originDist && originDist <= bondingUpperBound)
		(this->contactMap[this->flip])[abs(originDist - middle)] = make_pair(
				indexA, indexB);
}

bool ConstraintCheck::angleViolationCheck(double fb[3][3], double tb[3][3]) {
	if (Settings::Constraint::angleLow >= Settings::Constraint::angleHigh)
		return false;

	Vector3d p1(fb[0][0], fb[0][1], fb[0][2]);
	Vector3d p2(fb[1][0], fb[1][1], fb[1][2]);
	Vector3d p3(fb[2][0], fb[2][1], fb[2][2]);

	Vector3d P1(tb[0][0], tb[0][1], tb[0][2]);
	Vector3d P2(tb[1][0], tb[1][1], tb[1][2]);
	Vector3d P3(tb[2][0], tb[2][1], tb[2][2]);

	Vector3d v1 = p2 - p1;
	Vector3d v2 = p3 - p1;
	Vector3d v3 = v1.cross(v2);

	Vector3d V1 = P2 - P1;
	Vector3d V2 = P3 - P1;
	Vector3d V3 = V1.cross(V2);

	Matrix3d m, M, R;
	m << v1(0), v2(0), v3(0), v1(1), v2(1), v3(1), v1(2), v2(2), v3(2);
	M << V1(0), V2(0), V3(0), V1(1), V2(1), V3(1), V1(2), V2(2), V3(2);

	R = M * m.inverse(); //rotation matrix
	Vector3d t = P1 - R * p1;  //translation

	//compute mean
	vector<Point*> xfHelB = helB->getAtoms();
	Vector3d mean(0, 0, 0);
	for (size_t iter = 0; iter < xfHelB.size(); iter++) {
		double *l = xfHelB[iter]->getLocation();
		Vector3d p(l[0], l[1], l[2]);
		mean += R * p;
	}
	mean = mean / xfHelB.size();

	Vector3d TB = t + mean;

	Vector3d eaB = Utils::RotMatToEuler(R);
	eaB[0] = eaB[0] * 180 / PI;
	eaB[2] = eaB[2] * 180 / PI;
	double cos_eaB1 = eaB[1];

	double angle = acos(eaB[1]) * 180 / PI;
	if (angle < Settings::Constraint::angleLow
			|| angle > Settings::Constraint::angleHigh)
		return true;


	return false;

}

bool ConstraintCheck::stericsViolated(Orientation* ort) {

	if (ort == NULL)
		return true;

	this->outerGrid = true;

	this->flip = ort->getFlipNum();

	map<double, pair<int, int> > oldContacts;// (this->contactMap[this->flip].begin(),this->contactMap[this->flip].end());

	this->contactMap[this->flip].clear();

	bool out = false;

	double fb[3][3], tb[3][3];
	ort->getFromTo(fb, tb);

	bool angle_violation_status = angleViolationCheck(fb, tb);
	ort->setAngleViolation(angle_violation_status);

	vector<Point*> first = this->helA->getAtoms();
	vector<Point*> second = this->helB->getXFAtoms(fb, tb);

	int indexA = 0;

	for (map<double, pair<int, int> >::iterator iter = oldContacts.begin();
			iter != oldContacts.end() && !out; iter++) {
		int ha = (*iter).second.first;
		int hb = (*iter).second.second;

		if (!this->currentGraph->isConnected(ha, hb)) {
			out = (out || stericViolated(first[ha], second[hb]));
			updatePossibleContactList(first[ha], second[hb], ha, hb);
		}
	}

	if (!out) {
		this->contactMap[this->flip].clear();

		for (vector<Point*>::iterator ita = first.begin();
				ita != first.end() && !out; ita++) {
			int indexB = 0;

			for (vector<Point*>::iterator itb = second.begin();
					itb != second.end(); itb++) {
				if (!this->currentGraph->isConnected(indexA, indexB)) {
					out = (stericViolated(*ita, *itb) || out);
					updatePossibleContactList(*ita, *itb, indexA, indexB);
				}
				indexB++;
			}
			indexA++;
		}
	}

	// because getXFAtoms() makes new atoms they need to be deleted
	for (vector<Point*>::iterator ita = second.begin(); ita != second.end();
			ita++) {
		delete *ita;
	}

	if (Settings::Constraint::wholeCollision) { // to consider all atoms not just interface atoms. if candidate_interactions is provided, then molecular unit consisted of only interface atoms
		vector<double> trans_b = Utils::getTransMatrix(fb, tb); //Utils::getSingleTransformFromTwo(fa,ta,fb,tb);
		bool is_nei_collide = is_collide(this->helA->getStree(),
				this->helA->getStree(), trans_b, this->nei_matrix);
		return out || is_nei_collide; // || is_collide(a->getStree(),b->getStree(),fa,ta,fb,tb);  //FOR VIRUS CASE!  //FIX HERE,, I COMMENDED SINCE I DELETE FA, TA
	} else
		return out;

}

list<pair<int, int> > ConstraintCheck::getContacts(int flipno) {
	map<double, pair<int, int> >::iterator iter;
	list<pair<int, int> > output;
	for (iter = this->contactMap[flipno].begin();
			iter != this->contactMap[flipno].end(); iter++) {
		output.push_back((*iter).second);
	}
	return output;
}
