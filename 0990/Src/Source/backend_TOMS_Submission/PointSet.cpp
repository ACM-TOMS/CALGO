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
#include "PointSet.h"

#include <algorithm>

#include "Settings.h"

using namespace std;

PointSet::PointSet() {
}

PointSet::PointSet(vector<Point*> atoms) {
	this->points = atoms;
}

PointSet::~PointSet() {
	for (size_t x = 0; x < this->points.size(); x++) {
		delete (this->points[x]);
	}
}

void PointSet::addAtom(Point *a) {
	this->points.push_back(a);
}

bool PointSet::contains(Point *a) {
	return (find(this->points.begin(), this->points.end(), a) != this->points.end());
}

bool PointSet::empty() {
	return this->points.empty();
}

size_t PointSet::size() {
	return this->points.size();
}

Point* PointSet::getAtomAt(size_t index) {
	if (index >= 0 && index < points.size()) {
		return this->points.at(index);
	}
	cerr << "getAtomAt: NULL " << index << endl;
	return NULL;
}

Point* PointSet::getAtomByLabel(string label) {
	for (vector<Point*>::iterator ait = points.begin(); ait != points.end();
			ait++) {
		if ((*ait)->getName() == label)
			return *ait;
	}
	cerr << "can't find the atom with label  " << label << endl;
	return NULL;
}

int PointSet::getIndexOf(Point* a) {
	int output = -1;
	for (size_t i = 0; (i < this->points.size()) && (output < 0); i++) {
		if (this->points.at(i) == a) {
			output = (int) i;
		}
	}
	return output;
}

vector<Point*> PointSet::getAtoms() {

	return this->points;
}

vector<Point*> PointSet::getXFAtoms(double from[][3], double to[][3],
		size_t repeat) {
	vector<double> trans_mat = Utils::getTransMatrix(from, to);
	return getXFAtoms(trans_mat, repeat);
}

vector<Point*> PointSet::getXFAtoms(const vector<double> &trans_mat,
		size_t repeat) {
	vector<Point*> output;
	for (size_t iter = 0; iter < this->points.size(); iter++) {
		Point *current = getXFAtom(iter, trans_mat, repeat);
		output.push_back(current);
	}
	return output;
}

Point* PointSet::getXFAtom(size_t index, const vector<double> &trans_mat,
		size_t repeat) {
	Point * current = new Point(this->points[index]);
	double *l = current->getLocation(); // l is the actual location so don't delete
	Vector newL, preL;
	newL = Vector(l).trans(trans_mat);
	for (size_t mer = 0; mer < repeat; mer++) {
		preL = newL;
		newL = preL.trans(trans_mat);
	}
	current->setLocation(newL[0], newL[1], newL[2]);
	return current;

}

void PointSet::centerAtoms() {
	size_t iter;
	double sumx = 0, sumy = 0, sumz = 0;
	for (iter = 0; iter < this->points.size(); iter++) {
		double *l = (this->points[iter])->getLocation();
		sumx += l[0];
		sumy += l[1];
		sumz += l[2];
	}
	sumx = sumx / this->points.size();
	sumy = sumy / this->points.size();
	sumz = sumz / this->points.size();

	for (iter = 0; iter < this->points.size(); iter++) {
		double *l = (this->points[iter])->getLocation();
		(this->points[iter])->setLocation(l[0] - sumx, l[1] - sumy, l[2] - sumz);
	}
}

vector<pair<int, int> > PointSet::getDumbbells() {
	vector<pair<int, int> > output;
	vector<Point*>::iterator iter1, iter2;

	double minAtom = 0, maxAtom = 0;
	int mainAxis = 2; // assumes z axis is main axis of helix

	bool middleDumbbells = false;
	if (Settings::RootNodeCreation::participatingAtomIndex_low
			< Settings::RootNodeCreation::participatingAtomIndex_high)
		middleDumbbells = true;

	if (middleDumbbells) {
		minAtom =
				this->points[Settings::RootNodeCreation::participatingAtomIndex_low]->getLocation()[mainAxis];
		maxAtom =
				this->points[Settings::RootNodeCreation::participatingAtomIndex_high]->getLocation()[mainAxis];
		cout << "minAtom " << minAtom << " maxAtom " << maxAtom << endl;
	}

	int i = 0;
	for (iter1 = this->points.begin(); iter1 != this->points.end(); iter1++) {
		int j = i;

		for (iter2 = iter1; iter2 != this->points.end(); iter2++) {

			Point *a, *b;
			a = *iter1;
			b = *iter2;

			if (!middleDumbbells
					|| (a->getLocation()[mainAxis] >= minAtom
							&& a->getLocation()[mainAxis] <= maxAtom
							&& b->getLocation()[mainAxis] >= minAtom
							&& b->getLocation()[mainAxis] <= maxAtom)) { // Specifcally looking for  pairs near the middle of the helices

				double curDist = Utils::dist(a->getLocation(),
						b->getLocation());

				if (curDist
						< Settings::RootNodeCreation::initial4DContactSeparation_high
						&& curDist
								> Settings::RootNodeCreation::initial4DContactSeparation_low) {
					pair<int, int> result;
					result.first = i;
					result.second = j;
					output.push_back(result);
				}
			}

			j++;
		}
		i++;
	}
	return output;
}

void PointSet::getLimits(double max_xyz[3], double min_xyz[3]) {
	for (int i = 0; i < 3; i++) {
		max_xyz[i] = this->hi[i];
		min_xyz[i] = this->lo[i];
	}
}

void PointSet::calcLimits() {
	for (size_t i = 0; i < this->points.size(); i++) {
		double *loc;
		loc = this->points[i]->getLocation();
		for (size_t j = 0; j < 3; j++) {
			if (i == 0) {
				this->hi[j] = loc[j];
				this->lo[j] = loc[j];
			}
			if (this->hi[j] < loc[j]) {
				this->hi[j] = loc[j];
			}

			if (this->lo[j] > loc[j]) {
				this->lo[j] = loc[j];
			}
		}
	}
}

void PointSet::buildStree() {
	this->stree = new Stree(this->points);
}

void PointSet::simplify(PredefinedInteractions &distfinder, bool first_label) {
	// build an index for fast query

	for (vector<Point*>::iterator iter = points.begin(); iter != points.end();
			iter++) {
		index[(*iter)->getName()] = (*iter);
	}

	//clear the old atom array
	this->points.clear();

	//for each atom pair in the dist table, find it and
	//insert it into the atom array
	map<string, Point*>::iterator idx_iter;
	for (PredefinedInteractions::dist_iterator iter = distfinder.dist1begin();
			iter != distfinder.dist1end(); iter++) {
		string label;
		if (first_label)
			label = iter->first.first;
		else
			label = iter->first.second;

		idx_iter = index.find(label);
		if (idx_iter == index.end()) {
			cerr << "can't find atom with label " << label << endl;
			continue;
		} else {

			// check if the atom is alread in the atoms vector
			// void duplication
			bool atom_exist = false;
			for (int i = 0; i < points.size(); i++) {
				if (points[i]->getName() == label) {
					atom_exist = true;
					break;
				}
			}

			if (!atom_exist) {
				points.push_back(idx_iter->second);
				cout << idx_iter->second->getLocation() << endl;
			}
		}
	}

	//re-calculate the bound
	this->calcLimits();
}

Stree* PointSet::getStree() {
	return this->stree;
}

ostream & operator<<(ostream & os, PointSet & h) {
	size_t i;
	for (i = 0; i < h.points.size(); i++) {
		os << "[" << i << "]" << *(h.points[i]);
	}
	return os;
}

