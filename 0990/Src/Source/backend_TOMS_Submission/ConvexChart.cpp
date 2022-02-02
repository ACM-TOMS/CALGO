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
 * Algorithm.cpp
 *
 *  Created on: 2014
 *      Author: Aysegul Ozkan
 */

#include "ConvexChart.h"

#include "Orientation.h"
#include "ActiveConstraintGraph.h"
#include "PointSet.h"
#include "Settings.h"

#include <math.h>
#include <queue>
#include <iomanip>
#include <iostream>
#include <string>
#include <cstdlib>
using namespace std;

#define PI 3.14159265

ConvexChart::ConvexChart(ActiveConstraintGraph* cgk, bool dense,
		CayleyParameterization* cayleyParam, PredefinedInteractions *df) {

	for (size_t i = 0; i < 6; i++) {
		this->minOfMin[i] = 0;
		this->maxOfMax[i] = 0;
	}

	this->currentGraph = cgk;

	this->helA = this->currentGraph->getMolecularUnitA();
	this->helB = this->currentGraph->getMolecularUnitB();
	this->df = df;

	verticesA = this->currentGraph->getVerticesA();
	verticesB = this->currentGraph->getVerticesB();

	contacts = this->currentGraph->getContacts();

	dir = -1;
	bin = 1;

	this->dense = dense;
	for (int i = 0; i < 6; i++)
		this->point[i] = 0;

	for (int i = 0; i < NO; i++)
		for (int j = 0; j < NO; j++) {
			this->param_lengthUpper[i][j] = -1;
			this->param_lengthLower[i][j] = -1;
			this->edge_length[i][j] = -1;
			this->param_length[i][j] = -1;
		}

	this->stepSize = this->currentGraph->getStepSize();
	if (dense) {
		this->stepSize = this->stepSize / 2;
		this->currentGraph->setStepSize(this->stepSize);
	}

	this->parameters = cayleyParam->getParameters();

	this->updateList = cayleyParam->getUpdateList();
	this->boundaryComputationWay = cayleyParam->getBoundaryComputationWay();
	this->partial3tree = cayleyParam->is_partial3tree();
	this->tetrahedra = cayleyParam->getTetras();

	computeFixEdgeDistances();
	stepSize_of_stepGridContact = this->stepSize;
}

ConvexChart::~ConvexChart() {

}

void ConvexChart::computeFixEdgeDistances() {

	for (size_t i = 0; i < this->verticesA.size(); i++) { //inside helix a
		for (size_t j = i + 1; j < this->verticesA.size(); j++) {
			this->edge_length[i][j] = Utils::dist(
					helA->getAtomAt(verticesA[i])->getLocation(),
					helA->getAtomAt(verticesA[j])->getLocation()); //The location retrieved from an atom is the actual pointer don't delete it
			setMinMaxParam(i, j);
		}
	}

	for (size_t i = 0; i < this->verticesB.size(); i++) { //inside helix b
		for (size_t j = i + 1; j < this->verticesB.size(); j++) {
			this->edge_length[i + NO / 2][j + NO / 2] = Utils::dist(
					helB->getAtomAt(verticesB[i])->getLocation(),
					helB->getAtomAt(verticesB[j])->getLocation()); //The location retrieved from an atom is the actual pointer don't delete it
			setMinMaxParam(i + NO / 2, j + NO / 2);
		}
	}

	if (!Settings::Sampling::short_range_sampling) {
		for (size_t i = 0; i < this->contacts.size(); i++) //between the contacts in helixA and helixB
				{
			int v1 = this->contacts[i].first;
			int v2 = this->contacts[i].second;
			this->edge_length[v2][v1] = df->sumOfRadii(getAtom(v1),
					getAtom(v2));
			setMinMaxParam(v2, v1);
		}
	}

}

bool ConvexChart::initializeChart(bool contin,
		ActiveConstraintRegion * region) {
	this->dense = dense;
	if (dense) {
		this->stepSize = this->stepSize / 2;
		this->currentGraph->setStepSize(this->stepSize);
	}

	//if you change the order of the parameters, min-max range can change ??
	for (int i = this->parameters.size() - 1; i >= 0; i--) { //to start from reverse if you put contacts at the end
		int v1 = this->parameters[i].first;
		int v2 = this->parameters[i].second;

		int isContact = find_pair_in_vector(this->contacts,
				this->parameters[i]);
		if (isContact >= 0 && Settings::Sampling::short_range_sampling) {
			this->param_lengthUpper[v1][v2] = df->bondingUpperBound(getAtom(v1),
					getAtom(v2));
			this->param_lengthLower[v1][v2] = df->bondingLowerBound(getAtom(v1),
					getAtom(v2));
		} else {
			computeRange(v1, v2, i);
		}

		this->maxOfMax[i] = this->param_lengthUpper[v1][v2];
		this->minOfMin[i] = this->param_lengthLower[v1][v2];

		if (region->getSpaceSize() == 0 || contin) {
			if (region->getWitness().empty()) //for initial dumbbell  or for a dumbbell which in fact has a witness point but not saved to its file (since we do save the witness only one leaf node)
				this->witnessPoint[i] = this->param_lengthLower[v1][v2];

			if (contin) {
				if (region->getCombinedSize() > 0) //space + witness not empty  //getSpace().empty()
						{
					this->point[i] = (*(region->getSpace().back()))[i]; //CONTINUE FROM LAST POINT BUT KEEP THE ORIGINAL MIDPOINT AT THE SAME TIME.
					if (region->getWitnessSize() > 0)
						this->witnessPoint[i] =
								(*(region->getWitness().front()))[i]; //start from the middle(witness) point
				} else
					this->point[i] = this->witnessPoint[i]; //there should set some value to witnessPoint even if the witness has bad angle and deleted later.
			}
		} else
			this->witnessPoint[i] = this->param_lengthLower[v1][v2]; //probably for the dense case

		if (this->witnessPoint[i] < 0 || std::isnan(this->witnessPoint[i])) { //quick fix and for debug purpose
			this->witnessPoint[i] = this->param_lengthLower[v1][v2];
			cout << "ERRORRRRRRRRs" << endl;
		}

		if (this->witnessPoint[i] < this->param_lengthLower[v1][v2]
				&& this->witnessPoint[i]
						> this->param_lengthLower[v1][v2] - 0.1) {
			this->witnessPoint[i] = this->param_lengthLower[v1][v2];
		} else if (this->witnessPoint[i] > this->param_lengthUpper[v1][v2]
				&& this->witnessPoint[i]
						< this->param_lengthUpper[v1][v2] + 0.1) {
			this->witnessPoint[i] = this->param_lengthUpper[v1][v2];
		}

		if (this->witnessPoint[i] < this->param_lengthLower[v1][v2]
				|| this->witnessPoint[i] > this->param_lengthUpper[v1][v2]) { //quick fix and for debug purpose
			cout << "ERRORRRRRRRRs  witnessPoint is out of boundary "
					<< this->witnessPoint[i] << " "
					<< this->witnessPoint[i] - this->param_lengthLower[v1][v2]
					<< " "
					<< this->param_lengthUpper[v1][v2] - this->witnessPoint[i]
					<< endl; //one reason is; when a parameter becomes a contact, its value is set to be sum of radius. since it is forced to be exact contact, the range for the new parameters may not fit the value of witness orientation.
			this->witnessPoint[i] = this->param_lengthLower[v1][v2];
			this->point[i] = this->witnessPoint[i];
		}

		if (!contin)
			this->point[i] = this->witnessPoint[i];

		this->param_length[v1][v2] = this->point[i];
		this->param_length[v2][v1] = this->point[i];

	} //point can be <= 6 dimensional. If it is 2 dimensional(when there is 2 parameter) then point[2,3,4,5] is null

	for (int i = 0; i < NO; i++)
		for (int j = i; j < NO; j++) {
			if (this->param_lengthLower[i][j] > this->param_lengthUpper[i][j]) {
				cout << "invalid description along parameter " << i << "-" << j
						<< " min " << this->param_lengthLower[i][j] << " max "
						<< this->param_lengthUpper[i][j] << endl;
			}
		}

	this->gridDone = false;

	if (!contin) {

		for (size_t a = this->parameters.size(); a < 6; a++) {
			this->minOfMin[a] = 0;
			this->maxOfMax[a] = 0;
		}
		region->setSpaceVolume(this->minOfMin, this->maxOfMax);
	}

	//FOR JACOBIAN SAMPLING
	//set middle point
	for (int m = 0; m < 6; m++)
		for (size_t p = 0; p < this->parameters.size(); p++)
			this->mjpoint[m][p] = this->witnessPoint[p];

	return true;
}

void ConvexChart::setMinMaxParam(int v1, int v2) {
	this->edge_length[v2][v1] = this->edge_length[v1][v2];
	this->param_lengthUpper[v1][v2] = this->edge_length[v1][v2];
	this->param_lengthUpper[v2][v1] = this->edge_length[v1][v2];
	this->param_lengthLower[v1][v2] = this->edge_length[v1][v2];
	this->param_lengthLower[v2][v1] = this->edge_length[v1][v2];
	this->param_length[v1][v2] = this->edge_length[v1][v2];
	this->param_length[v2][v1] = this->edge_length[v1][v2];
}

void ConvexChart::setRangeByTriangleInequality(int v1, int v2) {

	int paramv1v2No = find_pair_in_vector(this->parameters, make_pair(v1, v2));

	if (this->edge_length[v1][v2] == -1
			|| std::isnan(this->edge_length[v1][v2])) {
		double max = 100;
		for (int i = 0; i < NO; i++) {

			int paramv1iNo = find_pair_in_vector(this->parameters,
					make_pair(v1, i));  // pair = <v, i> is not ordered
			int paramv2iNo = find_pair_in_vector(this->parameters,
					make_pair(v2, i));
			if ((paramv1iNo >= 0 && paramv1iNo < paramv1v2No)
					|| (paramv2iNo >= 0 && paramv2iNo < paramv1v2No)) // if v1-i is a parameter that is not set yet, then skip. order of the parameter is important.
				continue; // todo instead of skipping why dont you use just the min or max value of it to decrease the range of this param
						  //       no, if you use it, it will lead you to wrong range, since the parameter value will be set/fixed to another value in next steps !

			double v1i = this->param_length[v1][i];
			if (v1i == -1)
				v1i = this->param_lengthUpper[v1][i];
			double iv2 = this->param_length[v2][i];
			if (iv2 == -1)
				iv2 = this->param_lengthUpper[v2][i];
			if (v1i != -1 && iv2 != -1 && max > (v1i + iv2)) {
				max = v1i + iv2;
			}
		}

		if (max != 100) {
			this->param_lengthUpper[v1][v2] = max;
			this->param_lengthUpper[v2][v1] = max;
		}

		double min = 0;
		for (int i = 0; i < NO; i++) {
			double localmin = 100;

			int paramv1iNo = find_pair_in_vector(this->parameters,
					make_pair(v1, i));
			int paramv2iNo = find_pair_in_vector(this->parameters,
					make_pair(v2, i));
			if ((paramv1iNo >= 0 && paramv1iNo < paramv1v2No)
					|| (paramv2iNo >= 0 && paramv2iNo < paramv1v2No))
				continue;

			double maxv1i = this->param_length[v1][i];
			double minv1i = this->param_length[v1][i];
			if (maxv1i == -1)
				maxv1i = this->param_lengthUpper[v1][i];
			if (minv1i == -1)
				minv1i = this->param_lengthLower[v1][i];

			double maxiv2 = this->param_length[v2][i];
			double miniv2 = this->param_length[v2][i];
			if (maxiv2 == -1)
				maxiv2 = this->param_lengthUpper[v2][i];
			if (miniv2 == -1)
				miniv2 = this->param_lengthLower[v2][i];

//----------------
			if (minv1i > maxiv2 && minv1i != -1 && maxiv2 != -1)
				localmin = minv1i - maxiv2; // finds the smallest difference of edges not to restrict feasible intervals
			else if (miniv2 > maxv1i && miniv2 != -1 && maxv1i != -1)
				localmin = miniv2 - maxv1i;
			else
				localmin = 0; //otherwise min and max ranges are intermingled, hence at the intersected places the difference of the edges is 0
//----------------

			if (min < localmin && localmin != 100)
				min = localmin;
		}

		Point * atomA = getAtom(v1);
		Point * atomB = getAtom(v2);
		double collisionLowerBound = df->collisionLowerBound(atomA, atomB);

//		if(Settings::Collision::stericConstraint)  	// uncomment to allow collison on parameters
		if (min < collisionLowerBound)
			min = collisionLowerBound;

		if (min != 0) {
			this->param_lengthLower[v1][v2] = min;
			this->param_lengthLower[v2][v1] = min;
		}

		// this case may happen if maxrange is less than sumofradii or by some mysterious reason  // todo print to see it
		if (this->param_lengthLower[v1][v2] > this->param_lengthUpper[v1][v2]) {
			this->param_lengthLower[v1][v2] = this->param_lengthUpper[v1][v2];
			this->param_lengthLower[v2][v1] = this->param_lengthUpper[v1][v2];
		}
	} else
		setMinMaxParam(v1, v2);

}

int ConvexChart::find_pair_in_vector(vector<pair<int, int> > ivector,
		pair<int, int> ipair) {
	int index = -1;
	for (vector<pair<int, int> >::iterator iter = ivector.begin();
			iter != ivector.end(); iter++) {
		index++;
		if (((*iter).first == ipair.first && (*iter).second == ipair.second)
				|| ((*iter).first == ipair.second
						&& (*iter).second == ipair.first))
			return index;
	}
	return -1;
}

vector<vector<int> > ConvexChart::getTetras() {
	return this->tetrahedra;
}

Point * ConvexChart::getAtom(int vertex) {
	if (vertex < NO / 2) {
		vertex = this->currentGraph->getA(vertex);
		return this->helA->getAtomAt(vertex);
	} else {
		vertex = this->currentGraph->getB(vertex - NO / 2);
		return this->helB->getAtomAt(vertex);
	}
}

void ConvexChart::setWitnessPoint(CayleyPoint* wp) {

	wp->getPoint(this->witnessPoint);

}

double ConvexChart::dynamicStepSize(int v1, int v2, int p) {
	double mid = (this->param_lengthUpper[v1][v2]
			+ this->param_lengthLower[v1][v2]) / 2;
	double ratio =
			(mid - this->point[p])
					/ (this->param_lengthUpper[v1][v2]
							- this->param_lengthLower[v1][v2]);

	double linearStep = this->stepSize;
	if (Settings::Sampling::dynamicStepSizeWithin == 1)
		linearStep = this->stepSize + 1.8 * this->stepSize * ratio; // 4. / 3.
	else if (Settings::Sampling::dynamicStepSizeWithin == 2)
		linearStep = this->stepSize - 1.8 * this->stepSize * ratio; // 4. / 3.

	if (this->param_lengthUpper[v1][v2] == this->param_lengthLower[v1][v2])
		linearStep = this->stepSize;

	return linearStep;
}

void ConvexChart::stepGrid() {

	bool dens = false;
	if (dense) {
		for (size_t i = 1; i < this->parameters.size(); i++) {
			double x = (this->point[i] - this->witnessPoint[i])
					/ (2 * this->stepSize);
			int y = ceil(x);
			x = x - y;
			if (abs(x) > 0.0001 && abs(x) < 0.9999) //  abs(x - int(x))>0.0001  //if point is in the middle
				dens = true;
		}
	}

	double temp, temp2;
	bool carry = true;
	int p = 0;
	int v1, v2;
	for (size_t i = 0; i < this->parameters.size(); i++) {
		if (carry) {
			v1 = this->parameters[i].first;
			v2 = this->parameters[i].second;

			double tempStepSize = this->stepSize;
			this->stepSize = dynamicStepSize(v1, v2, p);

			double tmpoint = this->witnessPoint[p];
			if (this->witnessPoint[p] < this->param_lengthLower[v1][v2]) //min and max connection is changing since i compute the range continuously
				tmpoint = this->param_lengthLower[v1][v2];
			else if (this->witnessPoint[p] > this->param_lengthUpper[v1][v2])
				tmpoint = this->param_lengthUpper[v1][v2];

			if (this->point[p] >= tmpoint)
				this->point[p] += this->stepSize;

			else
				this->point[p] -= this->stepSize;

			if (dense && p == 0 && !dens) //  ! if at least one of the other parameters are in the middle of step interval
					{
				temp = tmpoint + this->stepSize;
				temp2 = tmpoint - this->stepSize;
				if (this->point[0] > temp) // != // to start from offset point
					this->point[p] += this->stepSize;
				else if (this->point[0] < temp2) // !=
					this->point[p] -= this->stepSize;
			}

			temp = this->param_lengthUpper[v1][v2] + this->stepSize;
			if (this->point[p] > this->param_lengthUpper[v1][v2]
					&& this->point[p] < temp) {
				this->point[p] = this->param_lengthUpper[v1][v2];
			} else if (this->point[p] > this->param_lengthUpper[v1][v2]) {
				this->point[p] = tmpoint - this->stepSize;
				if (this->point[p] > this->param_lengthUpper[v1][v2]) //maxconnection is changing since i compute the range continuously
					this->point[p] = this->param_lengthUpper[v1][v2];

			}
			carry = false;
			this->param_length[v1][v2] = this->point[p];
			this->param_length[v2][v1] = this->point[p];

			temp = this->param_lengthLower[v1][v2] - this->stepSize;
			if (this->point[p] < this->param_lengthLower[v1][v2]
					&& this->point[p] > temp) {
				this->point[p] = this->param_lengthLower[v1][v2];
				this->param_length[v1][v2] = this->point[p];
				this->param_length[v2][v1] = this->point[p];
				carry = false;
			} else if (this->point[p] < this->param_lengthLower[v1][v2]) {
				carry = true;
				this->point[p] = tmpoint;

				this->param_length[v1][v2] = this->point[p];
				this->param_length[v2][v1] = this->point[p];
				p++;
			}

			this->stepSize = tempStepSize;

		}

	}

	for (int x = p; x >= 0; x--) {
		//should update min max range for parameters from 0th till p'th  not just p'th parameter because if p'th parameters value is changed then that means every previous parameters value is also changed by setting it to mid value.

		if (x < this->parameters.size()) {
			for (int i = 0; i < updateList[x].size(); i++) {
				int paramNo = updateList[x][i];
				if (boundaryComputationWay[paramNo] != -2) // to check if it is not contact
						{
					v1 = this->parameters[paramNo].first;
					v2 = this->parameters[paramNo].second;

					computeRange(v1, v2, paramNo);
				}
			}
		}
	}

	for (int x = p - 1; x >= 0; x--) { //should set it to midpoint after min max range is recomputed
		if (boundaryComputationWay[x] != -2)  // to check if it is not contact
				{
			v1 = this->parameters[x].first;
			v2 = this->parameters[x].second;

			this->point[x] = this->witnessPoint[x]; // to make it always start from real mid point not tmpoint one. since tmpoint is rounded version of the moment when witnessPoint hits the boundary.
			if (this->point[x] < this->param_lengthLower[v1][v2]) {
				this->point[x] = this->param_lengthLower[v1][v2];
				this->param_length[v1][v2] = this->point[x];
				this->param_length[v2][v1] = this->point[x];
			} else if (this->point[x] > this->param_lengthUpper[v1][v2]) {
				this->point[x] = this->param_lengthUpper[v1][v2];
				this->param_length[v1][v2] = this->point[x];
				this->param_length[v2][v1] = this->point[x];
			}
		}
	}
	this->gridDone = carry;

}

void ConvexChart::updateBoundariesDependent(int x) {

	for (int i = 0; i < updateList[x].size(); i++) {
		int paramNo = updateList[x][i];
		if (boundaryComputationWay[paramNo] != -2) // to check if it is not contact
				{
			int v1 = this->parameters[paramNo].first;
			int v2 = this->parameters[paramNo].second;

			computeRange(v1, v2, paramNo);
		}
	}

}

void ConvexChart::updateBoundariesSpecific(int paramNo) {
	int v1 = this->parameters[paramNo].first;
	int v2 = this->parameters[paramNo].second;

	double curr = param_length[v1][v2];

	computeRange(v1, v2, paramNo);

	param_length[v1][v2] = curr; //since setRangeByTetrahedralInequality changes params value
	param_length[v2][v1] = param_length[v1][v2];
}

void ConvexChart::computeRange(int v1, int v2, int paramNo) {
	if (boundaryComputationWay[paramNo] == -1)
		setRangeByTriangleInequality(v1, v2);
	else
		setRangeByTetrahedralInequality(v1, v2,
				boundaryComputationWay[paramNo]); //this method not only sets min, max but also param ! check if that is what u want !!
}

void ConvexChart::updateBoundaries() {

	for (int i = this->parameters.size() - 1; i >= 0; i--) { //to start from reverse if you put contacts at the end
		int v1 = this->parameters[i].first;
		int v2 = this->parameters[i].second;
		int isContact = find_pair_in_vector(this->contacts,
				this->parameters[i]);

		if (isContact >= 0 && Settings::Sampling::short_range_sampling) {
			this->param_lengthUpper[v1][v2] = df->bondingUpperBound(getAtom(v1),
					getAtom(v2));
			this->param_lengthLower[v1][v2] = df->bondingLowerBound(getAtom(v1),
					getAtom(v2));
		} else {
			double curr = param_length[v1][v2];

			computeRange(v1, v2, i);

			param_length[v1][v2] = curr; //since setRangeByTetrahedralInequality changes params value
			param_length[v2][v1] = param_length[v1][v2];
		}
	}
}

void ConvexChart::setInitialPoint(vector<double> p) {

	for (size_t i = 0; i < this->parameters.size(); i++) {
		int v1 = this->parameters[i].first;
		int v2 = this->parameters[i].second;
		this->point[i] = p[i];
		this->param_length[v1][v2] = this->point[i];
		this->param_length[v2][v1] = this->point[i];
	}
	bin = 1;
}

void ConvexChart::setDir() {
	dir = -1;
}

bool ConvexChart::stepAround() {
	stepSize_of_stepGridContact = this->stepSize;

	if (dir == 2 * this->parameters.size() - 1)
		return true; //done
	else {
		dir++;

		int v1 = this->parameters[dir / 2].first;
		int v2 = this->parameters[dir / 2].second;

		if (dir % 2 == 0)
			this->point[dir / 2] += this->stepSize; //forward direction
		else
			this->point[dir / 2] -= this->stepSize; //backward direction
		// allow this->point to be less than edge_lengthLower since the purpose of this method is to search for a colliding config to start binary search in between.

		if (this->point[dir / 2] > this->param_lengthUpper[v1][v2]) {
			stepSize_of_stepGridContact = this->stepSize
					- (this->point[dir / 2] - this->param_lengthUpper[v1][v2]);
			this->point[dir / 2] = this->param_lengthUpper[v1][v2];
		}

		this->param_length[v1][v2] = this->point[dir / 2];
		this->param_length[v2][v1] = this->point[dir / 2];

		return false;
	}
}

void ConvexChart::stepGridBinary(bool valid) {

	int v1 = this->parameters[dir / 2].first;
	int v2 = this->parameters[dir / 2].second;

	int times = -1; //walk through colliding point.
	if (!valid)   // collision
		times = 1;  //walk through valid point

	if (dir % 2 == 0)
		this->point[dir / 2] -= times * stepSize_of_stepGridContact
				/ pow(2., bin);   //forward direction
	else
		this->point[dir / 2] += times * stepSize_of_stepGridContact
				/ pow(2., bin);   //backward direction

	bin++; //to half the step size one more time during binary search

	if (this->point[dir / 2] > this->param_lengthUpper[v1][v2])
		this->point[dir / 2] = this->param_lengthUpper[v1][v2];

	if (this->point[dir / 2] < this->param_lengthLower[v1][v2]) ///
		this->point[dir / 2] = this->param_lengthLower[v1][v2];

	this->param_length[v1][v2] = this->point[dir / 2];
	this->param_length[v2][v1] = this->point[dir / 2];

}

bool ConvexChart::doneGrid() {
	return this->gridDone;
}

vector<double> ConvexChart::getPoint() {
	vector<double> out;
	for (size_t i = 0; i < this->parameters.size(); i++)
		out.push_back(point[i]);
	return out;
}

vector<double> ConvexChart::getParamPoint() {
	vector<double> out;
	for (size_t p = 0; p < this->parameters.size(); p++) {
		int v1 = this->parameters[p].first;
		int v2 = this->parameters[p].second;
		out.push_back(this->param_length[v1][v2]);
	}
	return out;
}

void ConvexChart::setRangeByTetrahedralInequality(int v1, int v2, int tetraNo) {
	setRangeByTriangleInequality(v1, v2);

	int inds[4];
	int no = 0;
	for (int m = 0; m < 4; m++)
		if (this->tetrahedra[tetraNo][m] != v1
				&& this->tetrahedra[tetraNo][m] != v2) {
			inds[no] = this->tetrahedra[tetraNo][m];
			no++;
		}

	inds[2] = v1;
	inds[3] = v2;

	// set the range of the parameter( inds[2]-inds[3] ) according to tetrahedra equality

	if (param_length[inds[0]][inds[1]] == -1
			|| param_length[inds[0]][inds[2]] == -1
			|| param_length[inds[1]][inds[2]] == -1
			|| param_length[inds[0]][inds[3]] == -1
			|| param_length[inds[1]][inds[3]] == -1)
		return; //in order to compute tetrahedral range all those edges have to have values. If the value of the edge will be computed in realization class by known positions of vertices, they may not have a value yet. In this case, we can not compute the range.

	double posi[4][3];
	posi[0][0] = 0;
	posi[0][1] = 0;
	posi[0][2] = 0;
	posi[1][0] = param_length[inds[0]][inds[1]];
	posi[1][1] = 0;
	posi[1][2] = 0;

	double ac = param_length[inds[0]][inds[2]], ab =
			param_length[inds[0]][inds[1]], bc = param_length[inds[1]][inds[2]];
	double cosCAB = (ac * ac + ab * ab - bc * bc) / (2 * ac * ab);
	double sinCAB = sqrt(abs(1 - cosCAB * cosCAB));
	posi[2][0] = param_length[inds[0]][inds[2]] * cosCAB;
	posi[2][1] = param_length[inds[0]][inds[2]] * sinCAB;
	posi[2][2] = 0;

	ac = param_length[inds[0]][inds[3]];
	ab = param_length[inds[0]][inds[1]];
	bc = param_length[inds[1]][inds[3]];
	cosCAB = (ac * ac + ab * ab - bc * bc) / (2 * ac * ab);
	sinCAB = sqrt(abs(1 - cosCAB * cosCAB));
	posi[3][0] = param_length[inds[0]][inds[3]] * cosCAB;
	posi[3][1] = param_length[inds[0]][inds[3]] * sinCAB;
	posi[3][2] = 0;

	double nposi[2];
	nposi[0] = posi[3][0];
	nposi[1] = -1 * posi[3][1];
	double minrange = sqrt(
			(posi[2][0] - posi[3][0]) * (posi[2][0] - posi[3][0])
					+ (posi[2][1] - posi[3][1]) * (posi[2][1] - posi[3][1]));
	double maxrange = sqrt(
			(posi[2][0] - nposi[0]) * (posi[2][0] - nposi[0])
					+ (posi[2][1] - nposi[1]) * (posi[2][1] - nposi[1]));

	if (minrange > maxrange) {
		double temp = minrange;
		minrange = maxrange;
		maxrange = temp;
	}
	if (minrange < param_lengthLower[inds[2]][inds[3]])
		minrange = param_lengthLower[inds[2]][inds[3]];
	if (maxrange > param_lengthUpper[inds[2]][inds[3]])
		maxrange = param_lengthUpper[inds[2]][inds[3]];

	if (minrange > maxrange) //This case may happen if triangular boundary is tighter than tetrahedral boundary.
		minrange = maxrange;

	param_lengthLower[inds[2]][inds[3]] = minrange;
	param_lengthLower[inds[3]][inds[2]] = minrange;
	param_lengthUpper[inds[2]][inds[3]] = maxrange;
	param_lengthUpper[inds[3]][inds[2]] = maxrange;
	param_length[inds[2]][inds[3]] = minrange;
	param_length[inds[3]][inds[2]] = minrange;

}
