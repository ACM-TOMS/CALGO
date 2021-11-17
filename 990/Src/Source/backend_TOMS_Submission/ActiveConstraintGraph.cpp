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
 * ActiveConstraintGraph.cpp
 *
 *  Created on: 2008-2014
 *      Author: Aysegul Ozkan
 */

#include "ActiveConstraintGraph.h"

#include "Point.h"
#include "Orientation.h"
#include "CayleyPoint.h"

#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std;

#define PI 3.14159265

ActiveConstraintGraph::ActiveConstraintGraph() {
}

ActiveConstraintGraph::ActiveConstraintGraph(ActiveConstraintGraph* ptrToCpy) {

	this->helA = ptrToCpy->getMolecularUnitA();
	this->helB = ptrToCpy->getMolecularUnitB();
	this->stepSize = ptrToCpy->stepSize;

	for (vector<pair<int, int> >::iterator iter =
			ptrToCpy->participants.begin();
			iter != ptrToCpy->participants.end(); iter++)
			{
		pair<int, int> temp;
		temp.first = (*iter).first;
		temp.second = (*iter).second;
		this->participants.push_back(temp);

		//first appeared contact is first added to the set
		if (find(this->unordered_setA.begin(), this->unordered_setA.end(),
				temp.first) == this->unordered_setA.end())
			this->unordered_setA.push_back((*iter).first);
		if (find(this->unordered_setB.begin(), this->unordered_setB.end(),
				temp.second) == this->unordered_setB.end())
			this->unordered_setB.push_back((*iter).second);
	}

	completeTo3by3Graph(); //this method will be executed again when a new contact is added later on. since verticesA/B includes vertex for the parameters as well, those vertices will be overwritten by future contact atoms.

}

ActiveConstraintGraph::ActiveConstraintGraph(
		vector<pair<int, int> > participants, PointSet* helA, PointSet* helB) {
	this->stepSize = Settings::Sampling::stepSize;

	this->helA = helA;
	this->helB = helB;

	this->parameters.clear();

	if (participants.empty())
		return;

	for (vector<pair<int, int> >::iterator iter = participants.begin();
			iter != participants.end(); iter++) {
		this->participants.push_back(*iter);

		//first appeared contact is first added to the set
		if (find(this->unordered_setA.begin(), this->unordered_setA.end(),
				(*iter).first) == this->unordered_setA.end())
			this->unordered_setA.push_back((*iter).first);
		if (find(this->unordered_setB.begin(), this->unordered_setB.end(),
				(*iter).second) == this->unordered_setB.end())
			this->unordered_setB.push_back((*iter).second);
	}

	if (!Settings::Statistics::createPseudoAtlas) //during createPseudoAtlas, contactGraph is created for every config. It is too inefficient and takes a lot of time to create parameters everytime. I am assuming, I do import RoadMap.txt file which is the list of graph with parameters from previous easal atlas generations.
	{
		completeTo3by3Graph();
	}

	convertParticipantsToContacts();

}

ActiveConstraintGraph::ActiveConstraintGraph(
		vector<pair<int, int> > participants, PointSet* helA, PointSet* helB,
		vector<pair<int, int> > PL) {
	this->stepSize = Settings::Sampling::stepSize;

	this->helA = helA;
	this->helB = helB;

	this->parameters.clear();

	if (participants.empty())
		return;

	for (vector<pair<int, int> >::iterator iter = participants.begin();
			iter != participants.end(); iter++) {
		this->participants.push_back(*iter);

		//first appeared contact is first added to the set  (for consensus)
		if (find(this->unordered_setA.begin(), this->unordered_setA.end(),
				(*iter).first) == this->unordered_setA.end())
			this->unordered_setA.push_back((*iter).first);
		if (find(this->unordered_setB.begin(), this->unordered_setB.end(),
				(*iter).second) == this->unordered_setB.end())
			this->unordered_setB.push_back((*iter).second);
	}

	this->verticesA.assign(this->unordered_setA.begin(),
			this->unordered_setA.end());
	this->verticesB.assign(this->unordered_setB.begin(),
			this->unordered_setB.end());

	// the vertex that 'first' appear in PL should be added to verticesA/B 'first'.
	for (vector<pair<int, int> >::iterator it = PL.begin(); it != PL.end();
			it++) {
		if (find(this->verticesA.begin(), this->verticesA.end(), (*it).first)
				== this->verticesA.end())
			this->verticesA.push_back((*it).first);

		if (find(this->verticesB.begin(), this->verticesB.end(), (*it).second)
				== this->verticesB.end())
			this->verticesB.push_back((*it).second);
	}

	this->paramlines = PL;
	setParamsFromParamlines();

	convertParticipantsToContacts();
}

ActiveConstraintGraph::~ActiveConstraintGraph() {

	this->participants.clear();
	this->paramlines.clear();

}

void ActiveConstraintGraph::convertParticipantsToContacts() {
	this->contacts.clear();
	for (vector<pair<int, int> >::iterator iter = this->participants.begin();
			iter != this->participants.end(); iter++) {
		pair<int, int> temp;
		temp.first = find_in_setA((*iter).first);
		temp.second = find_in_setB((*iter).second) + NO / 2;
		this->contacts.push_back(temp);
	}
}

void ActiveConstraintGraph::setParamsFromParamlines() {
	this->parameters.clear();
	for (vector<pair<int, int> >::iterator iter = this->paramlines.begin();
			iter != this->paramlines.end(); iter++) {
		int pa = find_in_verticesA((*iter).first);
		int pb = find_in_verticesB((*iter).second) + NO / 2;
		this->parameters.push_back(make_pair(pa, pb));
	}
}

void ActiveConstraintGraph::completeTo3by3Graph() {
	bool closestAtom = true;

	this->verticesA.assign(this->unordered_setA.begin(),
			this->unordered_setA.end());
	this->verticesB.assign(this->unordered_setB.begin(),
			this->unordered_setB.end());

	if (this->unordered_setA.size() <= 2)  //find missing vertices
			{
		vector<int> toAdd;
		if (closestAtom)
			toAdd = completeTo3by3Graph_byClosestAtom(helA, this->verticesA,
					this->unordered_setA.size());
		else
			//perpendicular angle
			toAdd = completeTo3by3Graph_that_makes_perpendicularAngle(helA,
					this->verticesA, this->unordered_setA.size());

		for (int i = 0; i < toAdd.size(); i++)
			this->verticesA.push_back(toAdd[i]);
	}

	if (this->unordered_setB.size() <= 2) //find missing vertices
			{
		vector<int> toAdd;
		if (closestAtom)
			toAdd = completeTo3by3Graph_byClosestAtom(helB, this->verticesB,
					this->unordered_setB.size());
		else
			//perpendicular angle
			toAdd = completeTo3by3Graph_that_makes_perpendicularAngle(helB,
					this->verticesB, this->unordered_setB.size());

		for (int i = 0; i < toAdd.size(); i++)
			this->verticesB.push_back(toAdd[i]);
	}

}

vector<int> ActiveConstraintGraph::completeTo3by3Graph_byClosestAtom(
		PointSet * hel, vector<int> vertices, int size) {
	vector<int> out_vertices;

	int hel_size = hel->getAtoms().size();
	int vertices_1 = -1;

	if (size == 1) //Find 2nd atom that is closest to first atom
			{
		vertices_1 = vertices[0] + 1; //no need to worry about duplicates since there was only 1 atom.
		if (vertices_1 >= hel_size)
			vertices_1 = vertices[0] - 1;

		out_vertices.push_back(vertices_1);
	} else
		vertices_1 = vertices[1];

	// Find arbitrary 3rd atom.
	int h2;
	int i = 0;
	do {
		h2 = i++;
	} while ((vertices[0] == h2 || vertices_1 == h2) && i < hel_size); //while h2 already exists in current vertices  and there is more atoms to search for

	out_vertices.push_back(h2);

	return out_vertices;
}

//todo actually it would be better to choose equilateral triangle. because in perpendicular case, 2 atoms can be so close to each other.
vector<int> ActiveConstraintGraph::completeTo3by3Graph_that_makes_perpendicularAngle(
		PointSet * hel, vector<int> vertices, int size) {
	vector<int> out_vertices;

	int hel_size = hel->getAtoms().size();

	int vertices_1 = -1;
	if (size == 1) //Find 2nd atom
			{
		double minDist = 100;
		for (int h = 0; h < hel_size; h++) {
			if (h == vertices[0])
				continue;
			double dist = Utils::dist(
					hel->getAtomAt(vertices[0])->getLocation(),
					hel->getAtomAt(h)->getLocation());
			if (dist <= minDist) {
				vertices_1 = h;
				minDist = dist;
			}
		}
		if (vertices_1 == -1)
			cout
					<< "ERROR in completeTo3by3Graph_that_makes_perpendicularAngle "
					<< endl;

		out_vertices.push_back(vertices_1);
	} else
		vertices_1 = vertices[1];

	// Find 3rd atom.
	double closestAngle = 90;
	int closestVertex = -1;
	for (int h = 0; h < hel_size; h++) {
		if (vertices[0] == h || vertices_1 == h) // if h already exists in current vertices
			continue;
		double dist_01 = Utils::dist(hel->getAtomAt(vertices[0])->getLocation(),
				hel->getAtomAt(vertices[1])->getLocation());
		double dist_0h = Utils::dist(hel->getAtomAt(vertices[0])->getLocation(),
				hel->getAtomAt(h)->getLocation());
		double dist_1h = Utils::dist(hel->getAtomAt(vertices[1])->getLocation(),
				hel->getAtomAt(h)->getLocation());
		double angle = acos(
				(dist_0h * dist_0h + dist_1h * dist_1h - dist_01 * dist_01)
						/ (2 * dist_0h * dist_1h)) * 180.0 / PI;
		double absangle = abs(90 - angle);
		if (closestAngle >= absangle) {
			closestAngle = absangle;
			closestVertex = h;
		}
	}
	out_vertices.push_back(closestVertex);

	return out_vertices;
}

int ActiveConstraintGraph::getA(int no) {
	if (this->verticesA.size() > no)
		return this->verticesA[no];
	else {
		cout << "getA ERROR no:" << no << " this->verticesA.size() "
				<< this->verticesA.size() << endl;
		return -1;
	}
}

int ActiveConstraintGraph::getB(int no) {
	if (this->verticesB.size() > no)
		return this->verticesB[no];
	else {
		cout << "getB ERROR" << endl << no << " this->verticesB.size() "
				<< this->verticesB.size() << endl;
		return -1;
	}
}

int ActiveConstraintGraph::find_in_setA(int atom) {

	int i = 0;
	for (list<int>::iterator iter = this->unordered_setA.begin();
			iter != this->unordered_setA.end(); iter++) {
		if (*iter == atom)
			return i;
		i++;
	}
	return -1;
}

int ActiveConstraintGraph::find_in_setB(int atom) {

	int i = 0;
	for (list<int>::iterator iter = this->unordered_setB.begin();
			iter != this->unordered_setB.end(); iter++) {
		if (*iter == atom)
			return i;
		i++;
	}
	return -1;
}

int ActiveConstraintGraph::verticesA_size() {
	return this->verticesA.size();
}

int ActiveConstraintGraph::verticesB_size() {
	return this->verticesB.size();
}

int ActiveConstraintGraph::find_in_verticesA(int atom) {

	int i = 0;
	for (vector<int>::iterator iter = this->verticesA.begin();
			iter != this->verticesA.end(); iter++) {
		if (*iter == atom)
			return i;
		i++;
	}
	return -1;
}

int ActiveConstraintGraph::find_in_verticesB(int atom) {

	int i = 0;
	for (vector<int>::iterator iter = this->verticesB.begin();
			iter != this->verticesB.end(); iter++) {
		if (*iter == atom)
			return i;
		i++;
	}
	return -1;
}

bool ActiveConstraintGraph::mayHaveThirdAtom() {

	if (this->participants.size() < 6
			|| (this->unordered_setA.size() > 2
					&& this->unordered_setB.size() > 2))
		return true;
	return false;
}

bool ActiveConstraintGraph::isDependent() {

	for (int i = 0; i < this->participants.size(); i++) {
		int count = 1;
		for (int j = i + 1; j < this->participants.size(); j++)
			if (this->participants[i].first == this->participants[j].first)
				count++;
		if (count >= 4)
			return true;
	}

	for (int i = 0; i < this->participants.size(); i++) {
		int count = 1;
		for (int j = i + 1; j < this->participants.size(); j++)
			if (this->participants[i].second == this->participants[j].second)
				count++;
		if (count >= 4)
			return true;
	}

	return false;
}

int ActiveConstraintGraph::constraintSize() {

	return this->participants.size();
}

size_t ActiveConstraintGraph::getDim() {
	int dim = 6 - this->participants.size();
	return dim;
}

size_t ActiveConstraintGraph::getParamDim() {

	return paramlines.size();
}

void ActiveConstraintGraph::setParamLines(vector<pair<int, int> > parameters) {
	vector<pair<int, int> > parlines;
	for (size_t i = 0; i < parameters.size(); i++) {
		int first = getA(parameters[i].first);
		int second = getB(parameters[i].second - NO / 2);
		parlines.push_back(make_pair(first, second));
	}
	this->paramlines = parlines;
}

vector<pair<int, int> > ActiveConstraintGraph::getParamLines() {
	return this->paramlines;
}

void ActiveConstraintGraph::addContact(pair<int, int> contact) {

	bool contactExisted = false;
	for (vector<pair<int, int> >::iterator iter = this->participants.begin();
			iter != this->participants.end(); iter++)
		if (((*iter).first == contact.first && (*iter).second == contact.second)) //that contact exists before
			contactExisted = true;

	if (contactExisted)
		return;

	this->participants.push_back(contact);

	//first appeared contact is first added to the set
	if (find(this->unordered_setA.begin(), this->unordered_setA.end(),
			contact.first) == this->unordered_setA.end()) //if cannot find, then push
		this->unordered_setA.push_back(contact.first);
	if (find(this->unordered_setB.begin(), this->unordered_setB.end(),
			contact.second) == this->unordered_setB.end())
		this->unordered_setB.push_back(contact.second);

	completeTo3by3Graph();

	convertParticipantsToContacts();
}

// returns true if this contactGraph is parent of input graph.
bool ActiveConstraintGraph::IsParentOf(ActiveConstraintGraph* other) { //if other is child of this, then it should have more contacts
	bool output = true;
	output = output
			&& (this->participants.size() <= other->participants.size());
	output = output && (this->getParamDim() >= other->getParamDim()); //to make sure 5d is parent of 4d even they have both 2 dumbbells

	if (output) { //checks if all contacts of this included in other.
		for (vector<pair<int, int> >::iterator iter =
				this->participants.begin();
				output && iter != this->participants.end(); iter++) {
			bool found = false;
			for (vector<pair<int, int> >::iterator iter2 =
					other->participants.begin();
					!found && iter2 != other->participants.end(); iter2++) {
				if ((*iter).first == (*iter2).first
						&& (*iter).second == (*iter2).second) {
					found = true;
				}
			}
			if (!found)
				output = false; // this CG has a different contact pair than otherCG.
		}
	}

	output = output && (this->helA == other->getMolecularUnitA());
	output = output && (this->helB == other->getMolecularUnitB());
	return output;
}

bool ActiveConstraintGraph::isConnected(int x, int y) {
	bool output = false;
	for (vector<pair<int, int> >::iterator it = this->participants.begin();
			it != this->participants.end(); it++) {
		if (it->first == x && it->second == y) {
			output = true;
		}
	}
	return output;
}

vector<pair<int, int> > ActiveConstraintGraph::getParticipants() {
	return this->participants;
}

vector<pair<int, int> > ActiveConstraintGraph::getContacts() {
	return this->contacts;
}

vector<int> ActiveConstraintGraph::getVerticesA() {
	return this->verticesA;
}

void ActiveConstraintGraph::setVerticesA(vector<int> va) {
	this->verticesA.clear();
	this->verticesA.assign(va.begin(), va.end());
}

void ActiveConstraintGraph::setVerticesB(vector<int> vb) {
	this->verticesB.clear();
	this->verticesB.assign(vb.begin(), vb.end());
}

vector<int> ActiveConstraintGraph::getVerticesB() {
	return this->verticesB;
}

list<int> ActiveConstraintGraph::get_unordered_setA() {
	return this->unordered_setA;
}

list<int> ActiveConstraintGraph::get_unordered_setB() {
	return this->unordered_setB;
}

PointSet* ActiveConstraintGraph::getMolecularUnitA() {
	return this->helA;
}

PointSet* ActiveConstraintGraph::getMolecularUnitB() {
	return this->helB;
}

ostream& operator<<(ostream& os, ActiveConstraintGraph &acg) {

	os << "participants" << endl;
	for (vector<pair<int, int> >::iterator iter = acg.participants.begin();
			iter != acg.participants.end(); iter++)
		os << (*iter).first << "-" << (*iter).second << "   ";
	os << endl;

	os << "unordered_setA" << endl;
	for (list<int>::iterator iter = acg.unordered_setA.begin();
			iter != acg.unordered_setA.end(); iter++)
		os << (*iter) << "   ";
	os << endl;

	os << "unordered_setB" << endl;
	for (list<int>::iterator iter = acg.unordered_setB.begin();
			iter != acg.unordered_setB.end(); iter++)
		os << (*iter) << "   ";
	os << endl;

	os << "verticesA" << endl;
	for (vector<int>::iterator iter = acg.verticesA.begin();
			iter != acg.verticesA.end(); iter++)
		os << (*iter) << "   ";
	os << endl;

	os << "verticesB" << endl;
	for (vector<int>::iterator iter = acg.verticesB.begin();
			iter != acg.verticesB.end(); iter++)
		os << (*iter) << "   ";
	os << endl;

	//////////////////////////

	os << "contacts.size " << acg.contacts.size() << " ";
	for (vector<pair<int, int> >::iterator iter = acg.contacts.begin();
			iter != acg.contacts.end(); iter++)
		os << "c " << iter->first << " " << iter->second << " ";
	os << endl;

	os << "paramlines.size " << acg.paramlines.size() << " ";
	for (vector<pair<int, int> >::iterator iter = acg.paramlines.begin();
			iter != acg.paramlines.end(); iter++)
		os << "p " << iter->first << " " << iter->second << " ";
	os << endl;

	os << "parameters.size " << acg.parameters.size() << " ";
	for (vector<pair<int, int> >::iterator iter = acg.parameters.begin();
			iter != acg.parameters.end(); iter++)
		os << "p " << iter->first << " " << iter->second << " ";
	os << endl;

	os << endl;

	return os;
}

double ActiveConstraintGraph::getStepSize() {
	return this->stepSize;
}

void ActiveConstraintGraph::setStepSize(double siz) {
	this->stepSize = siz;
}

int ActiveConstraintGraph::getCNoinA() {

	return this->unordered_setA.size();

}

int ActiveConstraintGraph::getCNoinB() {

	return this->unordered_setB.size();
}

vector<pair<int, int> > ActiveConstraintGraph::getParameters() {
	return this->parameters;
}

void ActiveConstraintGraph::setParameters(vector<pair<int, int> > parameters) {
	this->parameters = parameters;
}

vector<ActiveConstraintGraph*> ActiveConstraintGraph::getAllAncestors() {
	vector<ActiveConstraintGraph*> ancestors;
	int power_of_two = 1;
	for (int i = 0; i < participants.size(); i++) {
		power_of_two *= 2;
	}

	for (int i = 0; i < power_of_two; i++) {

		vector<pair<int, int> > parts;
		//decode
		int dupi = i;
		for (int j = 0; j < this->participants.size(); j++) {
			if (dupi & 0x01) {
				parts.push_back(participants[j]);
			}
			dupi >>= 1;
		}
		if (parts.size()
				< 6 - Settings::RootNodeCreation::dimension_of_rootNodes)
			continue;
		sort(parts.begin(), parts.end());
		ActiveConstraintGraph* cgnew = new ActiveConstraintGraph(parts,
				this->helA, this->helB);
		// guaranteed >= 0 since it returns unsigned (size_t)
		if (cgnew->getDim()
				<= Settings::RootNodeCreation::dimension_of_rootNodes) { // && cgnew->getDim() >= 0
			ancestors.push_back(cgnew);
		} else
			delete cgnew;
	}
	return ancestors;
}

bool CGcomp(ActiveConstraintGraph* lhs, ActiveConstraintGraph* rhs) {
	//debug
	cout << *lhs << endl;
	cout << *rhs << endl;

	// first compare dim
	if (lhs->getDim() != rhs->getDim())
		return lhs->getDim() > rhs->getDim();
	// same dim, compare each pair
	else {
		for (int i = 0; i < lhs->participants.size(); i++) {
			if (lhs->participants[i] != rhs->participants[i])
				return lhs->participants[i] < rhs->participants[i];
		}
	}
	return false;
}

