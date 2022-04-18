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
 * CayleyParameterization.cpp
 *
 *  Created on: 2014
 *      Author: Aysegul Ozkan
 */

#include "CayleyParameterization.h"

#include "Orientation.h"
#include "Settings.h"
#include "ActiveConstraintRegion.h"

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <queue>
#include <iomanip>      // std::setw#include <time.h>
#include <iostream>     // std::cout#include <algorithm>    // std::next_permutation, std::sortusing namespace std;

#define PI 3.14159265

/*
 * Set of complete 3 trees
 * Each represented by 6 edges where an edge is represented by <helix A "vertex" index, helix B "vertex" index>
 * Note that the edges between vertices of same helix are also part of the complete 3-tree.
 *
 * For the sake of simplicity, 4*3, 5*3 and 6*3 complete 3 trees are not given in the list since their symmetric ones are given.
 */
static int complete3trees[11][6][2] = { { { 0, 0 }, { 0, 1 }, { 0, 2 },
		{ 1, 1 }, { 1, 2 }, { 2, 2 } }, { { 0, 0 }, { 0, 1 }, { 0, 2 },
		{ 1, 1 }, { 1, 2 }, { 2, 0 } }, { { 0, 0 }, { 0, 1 }, { 1, 0 },
		{ 1, 1 }, { 2, 1 }, { 2, 2 } }, // 3*3
		{ { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 1 }, { 1, 2 }, { 2, 3 } }, {
				{ 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 2 }, { 1, 3 }, { 2, 3 } }, {
				{ 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 2 }, { 1, 3 }, { 2, 2 } }, {
				{ 0, 0 }, { 0, 1 }, { 1, 1 }, { 1, 2 }, { 1, 3 }, { 2, 3 } }, // 3*4
		{ { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 2 }, { 1, 3 }, { 2, 4 } }, {
				{ 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 2 }, { 2, 3 }, { 2, 4 } }, {
				{ 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 3 }, { 1, 4 }, { 2, 4 } }, // 3*5
		{ { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 3 }, { 1, 4 }, { 2, 5 } } // 3*6
};

CayleyParameterization::CayleyParameterization(ActiveConstraintGraph* cgk,
		bool basic) {
	this->currentGraph = cgk;

	this->helA = this->currentGraph->getMolecularUnitA();
	this->helB = this->currentGraph->getMolecularUnitB();

	verticesA = this->currentGraph->getVerticesA();
	verticesB = this->currentGraph->getVerticesB();

	this->contacts = this->currentGraph->getContacts();

	cout << this->currentGraph->getParameters().size() << endl;

	if (this->currentGraph->getParameters().size() > 0) //if currentGraph is created before (or maybe use_existing_parametrization for symmetric node)
		this->parameters = this->currentGraph->getParameters();
	else {
		if (Settings::AtlasBuilding::parameterMinDeviation
				&& this->contacts.size() == 1) {
			this->parameters = parameterMinDeviation();
			currentGraph->setVerticesA(this->verticesA);
			currentGraph->setVerticesB(this->verticesB);
		} else
			this->parameters = defineParameters();

		if (Settings::Sampling::short_range_sampling) {
			//add contacts
			for (vector<pair<int, int> >::iterator iter =
					this->contacts.begin(); iter != this->contacts.end();
					iter++)
				this->parameters.push_back(*iter);
		}
	}

	if (!basic) {
		built3tree();
		orderParameters();
	}

	//update currentGraph's parameters to this new ordered parameters
	this->currentGraph->setParamLines(this->parameters);
	this->currentGraph->setParameters(this->parameters);

	cout << *this->currentGraph;

}

CayleyParameterization::~CayleyParameterization() {
}

vector<pair<int, int> > CayleyParameterization::getParameters() {
	return this->parameters;
}

bool CayleyParameterization::is_partial3tree() {
	return this->partial3tree;
}

vector<vector<int> > CayleyParameterization::getTetras() {
	return this->tetrahedra;
}

vector<vector<int> > CayleyParameterization::getUpdateList() {
	return this->updateList;
}

vector<int> CayleyParameterization::getBoundaryComputationWay() {
	return this->boundaryComputationWay;
}

vector<pair<int, int> > CayleyParameterization::defineParameters() {
	vector<pair<int, int> > temp_parameters;

	int dim = this->currentGraph->getDim();
	if (dim == 0)
		//todo: clean params
		return temp_parameters;

	int asize = this->currentGraph->verticesA_size();
	int bsize = this->currentGraph->verticesB_size();

	int indicesA[asize];
	for (int i = 0; i < asize; i++)
		indicesA[i] = i;

	int indicesB[bsize];
	for (int i = 0; i < bsize; i++)
		indicesB[i] = i;

	int begin = 0, end = -1; // the range of the indices where to search in the complete3trees
	bool use_symmetric = false; //if currentGraph is 4*3, 5*3 or 6*3 then use symmetric of the graph given in the list of complete3trees
	if (asize == 3 && bsize == 3) {
		begin = 0;
		end = 2;
	} else if ((asize == 4 && bsize == 3) || (asize == 3 && bsize == 4)) {
		begin = 3;
		end = 6;
		if (asize == 4 && bsize == 3)
			use_symmetric = true;
	} else if ((asize == 5 && bsize == 3) || (asize == 3 && bsize == 5)) {
		begin = 7;
		end = 9;
		if (asize == 5 && bsize == 3)
			use_symmetric = true;
	} else if ((asize == 6 && bsize == 3) || (asize == 3 && bsize == 6)) {
		begin = 10;
		end = 10;
		if (asize == 6 && bsize == 3)
			use_symmetric = true;
	}

	bool matched = false;
	for (int i = begin; i <= end; i++) //check corresponding graph according to number of contact atoms
			{
		matched = true;

		int graph[6][2];
		for (int j = 0; j < 6; j++) {
			if (use_symmetric) {
				graph[j][0] = complete3trees[i][j][1];
				graph[j][1] = complete3trees[i][j][0];
			} //contacts is symmetric of the graph given in the list complete3trees
			else {
				graph[j][0] = complete3trees[i][j][0];
				graph[j][1] = complete3trees[i][j][1];
			}
		}

		// checks if the contacts of currentGraph is part of any isomorphisms of the "graph"
		do {  // permutates indicesA

			do {  // permutates indicesB

				matched = true;
				for (vector<pair<int, int> >::iterator iter =
						this->contacts.begin(); iter != this->contacts.end();
						iter++) {
					bool found_the_contact = false;

					for (size_t j = 0; j < 6; j++)
						if (indicesA[graph[j][0]] == (*iter).first
								&& indicesB[graph[j][1]] + NO / 2
										== (*iter).second) {
							found_the_contact = true;
							break;
						}

					if (!found_the_contact) {
						matched = false;
						break;
					}
				}

			} while (!matched
					&& std::next_permutation(indicesB, indicesB + bsize));
		} while (!matched && std::next_permutation(indicesA, indicesA + asize));

		if (matched) {
			for (int p = 0; p < 6; p++)
				temp_parameters.push_back(
						make_pair(indicesA[graph[p][0]],
								indicesB[graph[p][1]] + NO / 2)); //this also includes contacts too, but in the next step they will be eliminated.

			break;
		}
	}

	if (!matched) {
		cout << "not matched " << endl;
		cout << "contacts ";
		for (vector<pair<int, int> >::iterator iter = this->contacts.begin();
				iter != this->contacts.end(); iter++)
			cout << (*iter).first << "-" << (*iter).second << ", ";
		cout << endl;
	}

// If the contacts are not part of any pre-determined graphs i.e. it is not partial 3-tree, then assign random parameters
	// In non-partial 3-tree case, if Settings::General::runSolver set to true, an external solver such as matlab can be used. Then those 'random' parameters will be used from matlab, etc.
	int all_possible_parameters[7][2] = { { 2, 1 }, { 1, 2 }, { 1, 0 },
			{ 0, 1 }, { 2, 2 }, { 2, 0 }, { 0, 2 } }; //i put 21 12 before 10, 01 to make sure the third atom is used
	if (!matched)
		for (int p = 0; p < 7; p++)
			temp_parameters.push_back(
					make_pair(all_possible_parameters[p][0],
							all_possible_parameters[p][1] + NO / 2));

	//delete contacts
	for (vector<pair<int, int> >::iterator iter = this->contacts.begin();
			iter != this->contacts.end(); iter++) {
		for (size_t j = 0; j < temp_parameters.size(); j++) {
			if (temp_parameters[j].first == (*iter).first
					&& temp_parameters[j].second == (*iter).second) {
				temp_parameters.erase(temp_parameters.begin() + j); // erase the j+1 th element
				break;
			}
		}
	}

	int no_params = 6 - this->contacts.size();
	if (temp_parameters.size() > no_params)
		temp_parameters.erase(temp_parameters.begin() + no_params,
				temp_parameters.end());

	return temp_parameters;

}

// this method finds the parameter set that gives min deviation so you would think to keep in ConvexChart class.
// this method creates verticesA and B, so you may think to keep it in ActiveConstraintGraph however
// its main job is to create parameters thats why we kept it here.
vector<pair<int, int> > CayleyParameterization::parameterMinDeviation() //works for 5d
{
	vector<pair<int, int> > temp_parameters;

	this->verticesA.assign(currentGraph->get_unordered_setA().begin(),
			currentGraph->get_unordered_setA().end());
	this->verticesB.assign(currentGraph->get_unordered_setB().begin(),
			currentGraph->get_unordered_setB().end());

	double volume = 0;
	double stepsize = 0;
	string thisparams = "";

	//assuming 5d case hence there is 1 contact
	int contacta, contactb;
	contacta = this->verticesA[0];
	contactb = this->verticesB[0];

	cout << contacta << "-" << contactb << ": ";

	string edgeset = "0007 0008 0106 0107 0108 0206 0207 0208 "; //all possible non-edges of 3by3 graph. 1 contact hence 3*3-1 = 8 non-edges

	int sizeA = helA->getAtoms().size();
	int sizeB = helB->getAtoms().size();

	// try to find the one with min deviation among all different atom combinations
	// you can assign (size-1 choose 2) options for 2vertices of contact graph. 1 vertex is already known by contact.
	double overallmindeviation = 100;
	for (int i1 = 0; i1 < sizeA; i1++) {
		if (i1 == contacta)
			continue;
		for (int i2 = i1 + 1; i2 < sizeA; i2++) {
			if (i2 == contacta)
				continue;
			for (int j1 = 0; j1 < sizeB; j1++) {
				if (j1 == contactb)
					continue;
				for (int j2 = j1 + 1; j2 < sizeB; j2++) {
					if (j2 == contactb)
						continue;
					int NO = 6;
					int verticesA[NO / 2];
					int verticesB[NO / 2];

					verticesA[0] = contacta;
					verticesB[0] = contactb;
					verticesA[1] = i1;
					verticesA[2] = i2;
					verticesB[1] = j1;
					verticesB[2] = j2;

					double edges[NO / 2][NO / 2];
					double connections[NO][NO];
					double maxconnections[NO][NO];
					double minconnections[NO][NO];

					for (int v1 = 0; v1 < NO; v1++)
						for (int v2 = 0; v2 < NO; v2++) {
							maxconnections[v1][v2] = -1;
							maxconnections[v2][v1] = -1;
							minconnections[v1][v2] = -1;
							minconnections[v2][v1] = -1;
							connections[v1][v2] = -1;
							connections[v2][v1] = -1;
						}
//------------------------------
					//set fixed distances for contacts
					Point * ac = helA->getAtomAt(contacta);
					Point * bc = helB->getAtomAt(contactb);
					double condist = ac->getRadius() + bc->getRadius(); // + Settings::Constraint::contactUpperThreshold;  //use distfinder sumofradii
					maxconnections[0][NO / 2] = condist;
					maxconnections[NO / 2][0] = condist;
					minconnections[0][NO / 2] = condist;
					minconnections[NO / 2][0] = condist;
					connections[0][NO / 2] = condist;
					connections[NO / 2][0] = condist;

//------------------------------
					//set fixed edge lengths inside helix a
					for (size_t v1 = 0; v1 < NO / 2; v1++) {
						for (size_t v2 = v1 + 1; v2 < NO / 2; v2++) {
							connections[v1][v2] =
									Utils::dist(
											helA->getAtomAt(verticesA[v1])->getLocation(),
											helA->getAtomAt(verticesA[v2])->getLocation());	//The location retrieved from an atom is the actual pointer don't delete it
							connections[v2][v1] = connections[v1][v2];
							maxconnections[v1][v2] = connections[v1][v2];
							maxconnections[v2][v1] = connections[v1][v2];
							minconnections[v1][v2] = connections[v1][v2];
							minconnections[v2][v1] = connections[v1][v2];
						}
					}

					//set fixed edge lengths inside helix b
					for (size_t v1 = NO / 2; v1 < NO; v1++) {
						for (size_t v2 = v1 + 1; v2 < NO; v2++) {
							connections[v1][v2] =
									Utils::dist(
											helB->getAtomAt(
													verticesB[v1 - NO / 2])->getLocation(),
											helB->getAtomAt(
													verticesB[v2 - NO / 2])->getLocation());//The location retrieved from an atom is the actual pointer don't delete it
							connections[v2][v1] = connections[v1][v2];
							maxconnections[v1][v2] = connections[v1][v2];
							maxconnections[v2][v1] = connections[v1][v2];
							minconnections[v1][v2] = connections[v1][v2];
							minconnections[v2][v1] = connections[v1][v2];
						}
					}

					//------------------------------
					// set max distances for non-edges
					bool change = true;
					while (change) //if the max distance of an edge is updated, then it may effect the max distance of other edges. So keep update process when there is a change in max distance of any edge.
					{
						change = false;
						for (int v1 = 0; v1 < NO / 2; v1++)
							for (int v2 = NO / 2; v2 < NO; v2++) {
								if (connections[v1][v2] == -1) { //if not contact, then find the max distance for edge v1-v2
									double max = 100;
									for (int i = 0; i < NO; i++)
										if (maxconnections[v1][i] != -1
												&& maxconnections[i][v2] != -1
												&& max
														> (maxconnections[v1][i]
																+ maxconnections[i][v2]))
											max = maxconnections[v1][i]
													+ maxconnections[i][v2];

									if (maxconnections[v1][v2] != max) {
										maxconnections[v1][v2] = max;
										maxconnections[v2][v1] = max;
										change = true;
									}
								}
							}
					}

//------------------------------
					// set min distances for non-edges
					change = true;
					while (change)  //updates min distances
					{
						change = false;
						for (int v1 = 0; v1 < NO / 2; v1++)
							for (int v2 = NO / 2; v2 < NO; v2++) {
								if (connections[v1][v2] == -1
										|| std::isnan(connections[v1][v2])) { // not a contact

									double min = 0;
									for (int i = 0; i < NO; i++) {
										double localmin = 100;

										if (minconnections[v1][i]
												> maxconnections[i][v2]
												&& minconnections[v1][i] != -1
												&& maxconnections[i][v2] != -1)
											localmin = minconnections[v1][i]
													- maxconnections[i][v2];
										else if (minconnections[v2][i]
												> maxconnections[i][v1]
												&& minconnections[v2][i] != -1
												&& maxconnections[i][v1] != -1)
											localmin = minconnections[v2][i]
													- maxconnections[i][v1];
										else
											localmin = 0;

										if (min < localmin && localmin != 100)
											min = localmin;
									}

									Point * a = helA->getAtomAt(verticesA[v1]);
									Point * b = helB->getAtomAt(
											verticesB[v2 - NO / 2]);

									double sumradii = a->getRadius()
											+ b->getRadius();
									double sum =
											Settings::Constraint::collisionLambda
													* sumradii
													+ Settings::Constraint::collisionDelta; //todo use distancefinder collisionLowerBound method

									if (min < sum)
										min = sum;

									if (min != minconnections[v1][v2]) {
										minconnections[v1][v2] = min;
										minconnections[v2][v1] = min;
										change = true;
									}
								}
							}
					}

//------------------------------
					//sets the length of distance interval.   it is positive for non-edges. 0 for contact edges.
					for (int v1 = 0; v1 < NO / 2; v1++)
						for (int v2 = NO / 2; v2 < NO; v2++) {
							edges[v1][v2 - NO / 2] = maxconnections[v1][v2]
									- minconnections[v1][v2];
						}

					double topdown[6][3]; //6 edges of (value, v1, v2).  (v1, v2 is indices of vertices in edges) and ( value is edge[v1][v2] )
					topdown[0][0] = 100;
					topdown[1][0] = 100;
					topdown[2][0] = 100;  //3 mins
					topdown[3][0] = 0;
					topdown[4][0] = 0;
					topdown[5][0] = 0;    //3 maxs

//------------------------------
					//find 3 min and 3 max distance intervals
					for (int v1 = 0; v1 < NO / 2; v1++)
						for (int v2 = 0; v2 < NO / 2; v2++) {
							if (v1 == 0 && v2 == 0) //skip first contact
								continue;

							if (edges[v1][v2] < topdown[0][0]) {
								topdown[2][0] = topdown[1][0];
								topdown[1][0] = topdown[0][0];
								topdown[0][0] = edges[v1][v2];

								topdown[2][1] = topdown[1][1];
								topdown[2][2] = topdown[1][2];
								topdown[1][1] = topdown[0][1];
								topdown[1][2] = topdown[0][2];
								topdown[0][1] = v1;
								topdown[0][2] = v2;
							} else if (edges[v1][v2] < topdown[1][0]) {
								topdown[2][0] = topdown[1][0];
								topdown[1][0] = edges[v1][v2];

								topdown[2][1] = topdown[1][1];
								topdown[2][2] = topdown[1][2];
								topdown[1][1] = v1;
								topdown[1][2] = v2;
							} else if (edges[v1][v2] < topdown[2][0]) {
								topdown[2][0] = edges[v1][v2];
								topdown[2][1] = v1;
								topdown[2][2] = v2;
							}

							if (edges[v1][v2] > topdown[5][0]) {
								topdown[3][0] = topdown[4][0];
								topdown[4][0] = topdown[5][0];
								topdown[5][0] = edges[v1][v2];

								topdown[3][1] = topdown[4][1];
								topdown[3][2] = topdown[4][2];
								topdown[4][1] = topdown[5][1];
								topdown[4][2] = topdown[5][2];
								topdown[5][1] = v1;
								topdown[5][2] = v2;
							} else if (edges[v1][v2] > topdown[4][0]) {
								topdown[3][0] = topdown[4][0];
								topdown[4][0] = edges[v1][v2];

								topdown[3][1] = topdown[4][1];
								topdown[3][2] = topdown[4][2];
								topdown[4][1] = v1;
								topdown[4][2] = v2;
							} else if (edges[v1][v2] > topdown[3][0]) {
								topdown[3][0] = edges[v1][v2];
								topdown[3][1] = v1;
								topdown[3][2] = v2;
							}

						}

//------------------------------
					//3 edges to be removed from top or down. compute the deviations if those 3 edges were removed from the set.
					double sum = 0, squares = 0;
					for (int v1 = 0; v1 < NO / 2; v1++)
						for (int v2 = 0; v2 < NO / 2; v2++) {
							sum += edges[v1][v2];
							squares += edges[v1][v2] * edges[v1][v2];
						}
					int edgeno = NO / 2 * NO / 2 - 3;

					//case 1: remove down 3
					double average1 = (sum
							- (topdown[0][0] + topdown[1][0] + topdown[2][0]))
							/ edgeno;
					double avgsquares1 = (squares
							- (topdown[0][0] * topdown[0][0]
									+ topdown[1][0] * topdown[1][0]
									+ topdown[2][0] * topdown[2][0])) / edgeno;
					double deviation1 = sqrt(avgsquares1 - average1 * average1);

					//case : remove top 3
					double average2 = (sum
							- (topdown[3][0] + topdown[4][0] + topdown[5][0]))
							/ edgeno;
					double avgsquares2 = (squares
							- (topdown[3][0] * topdown[3][0]
									+ topdown[4][0] * topdown[4][0]
									+ topdown[5][0] * topdown[5][0])) / edgeno;
					double deviation2 = sqrt(avgsquares2 - average2 * average2);

					//case 3: remove down 1 top 2
					double average3 = (sum
							- (topdown[0][0] + topdown[4][0] + topdown[5][0]))
							/ edgeno;
					double avgsquares3 = (squares
							- (topdown[0][0] * topdown[0][0]
									+ topdown[4][0] * topdown[4][0]
									+ topdown[5][0] * topdown[5][0])) / edgeno;
					double deviation3 = sqrt(avgsquares3 - average3 * average3);

					//case 4: remove down 2 top 1
					double average4 = (sum
							- (topdown[0][0] + topdown[1][0] + topdown[5][0]))
							/ edgeno;
					double avgsquares4 = (squares
							- (topdown[0][0] * topdown[0][0]
									+ topdown[1][0] * topdown[1][0]
									+ topdown[5][0] * topdown[5][0])) / edgeno;
					double deviation4 = sqrt(avgsquares4 - average4 * average4);

//------------------------------
					// finds the minimum deviation among above 4 cases
					string removeset; //3 non-edges to be removed
					double mindeviation = 100;
					int N = 6;
					//make sure these 3 edges are not from all distinct vertices, to keep partial 3-tree
					//also make sure these 3 edges are not from all same vertex, (the vertex cannot hang out by itself)
					if (mindeviation > deviation1
							&& !(topdown[0][1] != topdown[1][1]
									&& topdown[1][1] != topdown[2][1]
									&& topdown[0][1] != topdown[2][1]
									&& topdown[0][2] != topdown[1][2]
									&& topdown[1][2] != topdown[2][2]
									&& topdown[0][2] != topdown[2][2]))
						if (!((topdown[0][1] == topdown[1][1]
								&& topdown[1][1] == topdown[2][1])
								|| (topdown[0][2] == topdown[1][2]
										&& topdown[1][2] == topdown[2][2]))) {
							mindeviation = deviation1;
							std::stringstream removeset1;
							removeset1 << 0 << topdown[0][1] << 0
									<< topdown[0][2] + N << " " << 0
									<< topdown[1][1] << 0 << topdown[1][2] + N
									<< " " << 0 << topdown[2][1] << 0
									<< topdown[2][2] + N << " ";
							removeset = removeset1.str();
						}
					if (mindeviation > deviation2
							&& !(topdown[3][1] != topdown[4][1]
									&& topdown[4][1] != topdown[5][1]
									&& topdown[3][1] != topdown[5][1]
									&& topdown[3][2] != topdown[4][2]
									&& topdown[4][2] != topdown[5][2]
									&& topdown[3][2] != topdown[5][2]))
						if (!((topdown[3][1] == topdown[4][1]
								&& topdown[4][1] == topdown[5][1])
								|| (topdown[3][2] == topdown[4][2]
										&& topdown[4][2] == topdown[5][2]))) {
							mindeviation = deviation2;
							std::stringstream removeset2;
							removeset2 << 0 << topdown[3][1] << 0
									<< topdown[3][2] + N << " " << 0
									<< topdown[4][1] << 0 << topdown[4][2] + N
									<< " " << 0 << topdown[5][1] << 0
									<< topdown[5][2] + N << " ";
							removeset = removeset2.str();
						}
					if (mindeviation > deviation3
							&& !(topdown[0][1] != topdown[4][1]
									&& topdown[4][1] != topdown[5][1]
									&& topdown[0][1] != topdown[5][1]
									&& topdown[0][2] != topdown[4][2]
									&& topdown[4][2] != topdown[5][2]
									&& topdown[0][2] != topdown[5][2]))
						if (!((topdown[0][1] == topdown[4][1]
								&& topdown[4][1] == topdown[5][1])
								|| (topdown[0][2] == topdown[4][2]
										&& topdown[4][2] == topdown[5][2]))) {
							mindeviation = deviation3;
							std::stringstream removeset3;
							removeset3 << 0 << topdown[0][1] << 0
									<< topdown[0][2] + N << " " << 0
									<< topdown[4][1] << 0 << topdown[4][2] + N
									<< " " << 0 << topdown[5][1] << 0
									<< topdown[5][2] + N << " ";
							removeset = removeset3.str();
						}
					if (mindeviation > deviation4
							&& !(topdown[0][1] != topdown[1][1]
									&& topdown[1][1] != topdown[5][1]
									&& topdown[0][1] != topdown[5][1]
									&& topdown[0][2] != topdown[1][2]
									&& topdown[1][2] != topdown[5][2]
									&& topdown[0][2] != topdown[5][2]))
						if (!((topdown[0][1] == topdown[1][1]
								&& topdown[1][1] == topdown[5][1])
								|| (topdown[0][2] == topdown[1][2]
										&& topdown[1][2] == topdown[5][2]))) {
							mindeviation = deviation4;
							std::stringstream removeset4;
							removeset4 << 0 << topdown[0][1] << 0
									<< topdown[0][2] + N << " " << 0
									<< topdown[1][1] << 0 << topdown[1][2] + N
									<< " " << 0 << topdown[5][1] << 0
									<< topdown[5][2] + N << " ";
							removeset = removeset4.str();
						}

//------------------------------
					//the removeset with min deviation will be removed from main edge set, and that (remaining edge set) will be our parameter set in description class
					string paramset = edgeset;
					for (size_t r = 0; r < removeset.size(); r = r + 5) {
						int pos = paramset.find(removeset.substr(r, 4), 0);
						if (pos >= 0)
							paramset.erase(pos, 5);
					}

					if (overallmindeviation > mindeviation) {
						overallmindeviation = mindeviation;

						thisparams = paramset;

						// to determine verticesA, uses the ordering of params (i.e. the vertex that 'first' appear in params should be added to verticesA 'first'. )
						this->verticesA.clear();
						this->verticesB.clear();
						for (int t = 0; t < paramset.size(); t = t + 5) {
							int v1 = atoi(paramset.substr(t, 2).c_str());
							int v2 = atoi(paramset.substr(t + 2, 2).c_str()) - N
									+ N / 2;

							if (find(this->verticesA.begin(),
									this->verticesA.end(), verticesA[v1])
									== this->verticesA.end())
								this->verticesA.push_back(verticesA[v1]);
							if (find(this->verticesB.begin(),
									this->verticesB.end(), verticesB[v2])
									== this->verticesB.end())
								this->verticesB.push_back(verticesB[v2]);
						}

						//computes the volume of parameter space
						volume = 1;
						for (int t = paramset.size() - 5; t >= 0; t = t - 5) {
							int v1 = atoi(paramset.substr(t, 2).c_str());
							int v2 = atoi(paramset.substr(t + 2, 2).c_str()) - N
									+ N / 2;
							if ((maxconnections[v1][v2] - minconnections[v2][v1])
									!= 0)
								volume = volume
										* (maxconnections[v1][v2]
												- minconnections[v2][v1]);
							//cout <<  v1 << ":" << verticesA[v1]  << "-" << v2+N/2 << ":" << verticesB[v2-N/2] << " min:" << minconnections[v1][v2] << " max:" << maxconnections[v1][v2] << " ";

						}
						stepsize = pow(volume / 1000000, 0.2);// root 1/5 = 0.2 //10 step per dimension hence in 5d 10^5

					}

				}
			}
		}
	}

	double stepFactor = Settings::Sampling::stepSize / 0.3; //assume 0.4 is base, than if we set Settings::Sampling::stepSize to be 0.6, than we will multiply mideviStepSize with this factor 1.5. So that stepsize will change according to user input as well
	stepsize = stepsize * stepFactor;

	currentGraph->setStepSize(stepsize);

	//Settings::Sampling::stepSize = stepsize; //commented on jan30 2014

	for (size_t p = 0; p < thisparams.size(); p = p + 5) {
		int v1 = atoi(thisparams.substr(p, 2).c_str());
		int v2 = atoi(thisparams.substr(p + 2, 2).c_str()); // - 6; //todo: that number 6 to be changed later
		temp_parameters.push_back(make_pair(v1, v2));
	}

	return temp_parameters;
}

//todo maybe you can have a hash table for that.
void CayleyParameterization::built3tree() {
	for (int i = 0; i < NO; i++)
		vertex_position[i] = false;

	for (int i = 0; i < NO; i++)
		for (int j = 0; j < NO; j++) {
			edge[i][j] = false;
			fixedEdge[i][j] = false;
		}

	int pA = this->currentGraph->getCNoinA();
	int pB = this->currentGraph->getCNoinB();

	if (pA < 3)
		pA = 3;
	if (pB < 3)
		pB = 3;

	//set true for distances inside of an helixA
	for (int i = 0; i < pA; i++)
		for (int j = 0; j < pA; j++)
			edge[i][j] = true;

	//set true for distances inside of an helixB
	for (int i = NO / 2; i < NO / 2 + pB; i++)
		for (int j = NO / 2; j < NO / 2 + pB; j++)
			edge[i][j] = true;

	//set true for distances for each contact
	for (size_t i = 0; i < this->contacts.size(); i++) {
		int ind1 = this->contacts[i].first;
		int ind2 = this->contacts[i].second;
		edge[ind1][ind2] = edge[ind2][ind1] = true;
	}

	//set true for distances for each param
	for (size_t i = 0; i < this->parameters.size(); i++) {
		int ind1 = this->parameters[i].first;
		int ind2 = this->parameters[i].second;
		edge[ind1][ind2] = edge[ind2][ind1] = true;
	}

	//builds 3tree by storing tetrahedrons
	bool non_partial3tree = false;
	if (find4Clique()) //base tetrahedron
	{
		//builts new locations onto the base tetrahedron.
		int new_located_vertex = -1;
		do {
			new_located_vertex = identifyRigidVertex();
			if (new_located_vertex != -1)
				identifyFixedEdgesConnectedToVertex(new_located_vertex);
		} while (new_located_vertex != -1);
		//now it should be able to fill all the values of vertices, else it means it needs maple

		for (size_t i = 0; i < 3; i++)
			if (!vertex_position[i] || !vertex_position[i + NO / 2])
				non_partial3tree = true;
	} else
		non_partial3tree = true;

	this->partial3tree = !non_partial3tree;

}

int CayleyParameterization::identifyRigidVertex() {
	vector<int> inds;
	for (int i = 0; i < NO; i++) //finds three indices i,j,k(that their location is known) where this point l is connected to this three indices.
		if (vertex_position[i])
			for (int j = i + 1; j < NO; j++)
				if (vertex_position[j])
					for (int k = j + 1; k < NO; k++)
						if (vertex_position[k]) //these points located before.
							for (int l = 0; l < NO; l++) //find a point with unknown location and then set it
								if (!vertex_position[l] && edge[i][l]
										&& edge[j][l] && edge[k][l]) {
									inds.push_back(l);
									inds.push_back(i);
									inds.push_back(j);
									inds.push_back(k);
									tetrahedra.push_back(inds);
									vertex_position[l] = true;
									return l;
								}
	return -1;
}

void CayleyParameterization::identifyFixedEdgesConnectedToVertex(
		int new_located_vertex) {

	for (int i = 0; i < NO; i++)
		if (vertex_position[i] && !edge[i][new_located_vertex]) {
			edge[i][new_located_vertex] = edge[new_located_vertex][i] = true;
			fixedEdge[i][new_located_vertex] =
					fixedEdge[new_located_vertex][i] = true;
		}
}

bool CayleyParameterization::find4Clique() {

	//find 4 click where all vertices is from helA
	vector<int> inds;
	for (int i = 0; i < NO / 2; i++)
		for (int j = i + 1; j < NO / 2; j++)
			if (edge[i][j])
				for (int k = j + 1; k < NO / 2; k++)
					if (edge[i][k] && edge[j][k])
						for (int l = k + 1; l < NO / 2; l++)
							if (edge[i][l] && edge[j][l] && edge[k][l]) {
								inds.push_back(j);
								inds.push_back(k);
								inds.push_back(l);
								inds.push_back(i);
								tetrahedra.push_back(inds);
								vertex_position[j] = vertex_position[k] =
										vertex_position[l] =
												vertex_position[i] = true;
								return true; //done
							}

	//find 4 click where all vertices is from helB
	for (int i = NO / 2; i < NO; i++)
		for (int j = i + 1; j < NO; j++)
			if (edge[i][j])
				for (int k = j + 1; k < NO; k++)
					if (edge[i][k] && edge[j][k])
						for (int l = k + 1; l < NO; l++)
							if (edge[i][l] && edge[j][l] && edge[k][l]) {
								inds.push_back(j);
								inds.push_back(k);
								inds.push_back(l);
								inds.push_back(i);
								tetrahedra.push_back(inds);
								vertex_position[j] = vertex_position[k] =
										vertex_position[l] =
												vertex_position[i] = true;
								return true; //done
							}

//find 4 click where first vertex is from helA and rest 3 vertices are from helB
	for (int i = 0; i < NO / 2; i++)
		for (int j = NO / 2; j < NO; j++)
			if (edge[i][j])
				for (int k = j + 1; k < NO; k++)
					if (edge[i][k] && edge[j][k])
						for (int l = k + 1; l < NO; l++)
							if (edge[i][l] && edge[j][l] && edge[k][l]) {
								inds.push_back(j);
								inds.push_back(k);
								inds.push_back(l);
								inds.push_back(i);
								tetrahedra.push_back(inds);
								vertex_position[j] = vertex_position[k] =
										vertex_position[l] =
												vertex_position[i] = true;
								return true; //done
							}

	//find 4 click where 2 vertices is from helA and rest 2 vertices are from helB
	for (size_t i = 0; i < NO / 2; i++)
		for (size_t j = i + 1; j < NO / 2; j++)
			if (edge[i][j])
				for (int k = NO / 2; k < NO; k++)
					if (edge[i][k] && edge[j][k])
						for (int l = k + 1; l < NO; l++)
							if (edge[i][l] && edge[j][l] && edge[k][l]) {
								inds.push_back(j);
								inds.push_back(k);
								inds.push_back(l);
								inds.push_back(i);
								tetrahedra.push_back(inds);
								vertex_position[j] = vertex_position[k] =
										vertex_position[l] =
												vertex_position[i] = true;
								return true; //done
							}

	//find 4 click where first 3 vertices is from helA and last vertex are from helB
	for (size_t i = 0; i < NO / 2; i++)
		for (size_t j = i + 1; j < NO / 2; j++)
			if (edge[i][j])
				for (size_t k = j + 1; k < NO / 2; k++)
					if (edge[i][k] && edge[j][k])
						for (int l = NO / 2; l < NO; l++)
							if (edge[i][l] && edge[j][l] && edge[k][l]) {
								inds.push_back(j);
								inds.push_back(k);
								inds.push_back(l);
								inds.push_back(i);
								tetrahedra.push_back(inds);
								vertex_position[j] = vertex_position[k] =
										vertex_position[l] =
												vertex_position[i] = true;
								return true; //done
							}

	return false;

}

void CayleyParameterization::printTetras() {
	cout << "tetrahedras " << endl;
	for (int i = 0; i < this->tetrahedra.size(); i++)
		cout << " " << this->tetrahedra[i][0] << " " << this->tetrahedra[i][1]
				<< " " << this->tetrahedra[i][2] << " "
				<< this->tetrahedra[i][3] << endl;

	cout << "verticesA ";
	for (int i = 0; i < verticesA.size(); i++)
		cout << verticesA[i] << " ";
	cout << endl;

	cout << "verticesB ";
	for (int i = 0; i < verticesB.size(); i++)
		cout << verticesB[i] << " ";
	cout << endl;

	cout << "tetrahedras with real atom indices" << endl;
	for (int i = 0; i < this->tetrahedra.size(); i++) {
		for (int j = 0; j < 4; j++) {
			if (this->tetrahedra[i][j] < 6)
				cout << " " << verticesA[this->tetrahedra[i][j]];
			else
				cout << " " << verticesB[this->tetrahedra[i][j] - 6];
		}
		cout << endl;
	}
}

void CayleyParameterization::orderParameters() {

	if (!this->partial3tree) {
		for (int i = 0; i < this->parameters.size(); i++)
			boundaryComputationWay.push_back(-1); //triangular
		return;
	}

	vector<pair<int, int> > oparameters;
	vector<pair<int, int> > paramss = this->parameters;
	for (int i = this->parameters.size() - 1; i >= 0; i--) { //put contacts at the beginning
		int isContact = Utils::findPair(this->contacts, this->parameters[i]);
		if (isContact >= 0 && Settings::Sampling::short_range_sampling) {
			oparameters.insert(oparameters.begin(), this->parameters[i]);
			paramss.erase(paramss.begin() + i);
		}
	}

	//tetraHedraParams[j][i]=true means: jth tetrahedron includes ith parameter
	bool tetraHedraParams[this->tetrahedra.size()][this->parameters.size()];

	for (int j = 0; j < this->tetrahedra.size(); j++) {
		for (int i = 0; i < this->parameters.size(); i++) {
			int v1 = this->parameters[i].first;
			int v2 = this->parameters[i].second;

			bool found = false;
			for (int m = 0; m < 4 && !found; m++)
				for (int n = 0; n < 4 && !found; n++)
					if (this->tetrahedra[j][m] == v1
							&& this->tetrahedra[j][n] == v2)
						found = true; //the parameter v1-v2 is found at tetrahedra j.

			tetraHedraParams[j][i] = found;
		}
	}

	int degreeOfParams[this->parameters.size()]; //the number of times that parameter exist in different tetrahedra.
	for (int i = 0; i < this->parameters.size(); i++) {
		degreeOfParams[i] = 0;
		for (int j = 0; j < this->tetrahedra.size(); j++) {
			if (tetraHedraParams[j][i])
				degreeOfParams[i]++;
		}
	}

	/*
	 * How to Order Parameters?
	 *
	 * Initially you would think the order of parameters while building 3-tree procedure.
	 *  However that ordering would not be that efficient.
	 *  for example:  consider the graph "0006 0107 ", "0207 0106 0008 0007 "
	 *  assume tetraHedras are created in the order: 0127, 0167, 0678
	 *  so the parameters that take place in tetrahedras are:
	 *  tetrahedra:   parameters
	 *  0127      :   27, 07
	 *  0167      :   07, 16
	 *  0678      :   07, 08
	 *  then lets sort the parameters in the reverse order that take place in tetrahedras
	 *  then you would have the parameters ordered as 08, 16, 07, 27
	 *  if the parameters take place in the same tetrahedra, then there is dependency relation between them.
	 *  which parameter depends on which parameter is defined according to order of parameters.
	 *  hence in this case,  08 depends on 07. 16 depends on 07. 07 depends on 27.
	 *  so whenever we change the value of the parameter 27, we need to update the boundary for parameter 07.
	 *  and whenever we change the value of the parameter 07, we need to update the boundary for parameters 16 and 08.
	 *  this makes dependency graph longer, hence requires more boundary computation
	 *
	 *  another bad parameter ordering is: 16, 07, 08, 27
	 *  remember: computing boundary of the parameter through setRange method uses only one tetrahedra
	 *  and requires all the edges other than input parameter to have a value.
	 *  hence when you change the value of 08, you need to update the boundary of 07.
	 *  however your new boundary of 07, may exceed the boundary defined by 27. hence cause infeasible solution!!!
	 *
	 *
	 *  however if you had sorted the parameters as 08, 16, 27, 07, then 27 would depend on 07
	 *  hence the chance on 27 would not cause a change on 07 and hence would not cause boundary computations on 16 and 08
	 */

	for (int j = 0; j < this->tetrahedra.size(); j++) {    //sets  paramsOrdered
		bool entered = true;
		while (entered) { //to get all parameters in the tetrahedra j.
			entered = false;
			int maxDegree = 0;
			int maxi = 0;
			int maxpos = 0;
			for (int i = 0; i < paramss.size(); i++) {
				int pos = Utils::findPair(this->parameters, paramss[i]);
				if (tetraHedraParams[j][pos]
						&& maxDegree < degreeOfParams[pos]) {
					maxDegree = degreeOfParams[pos];
					maxi = i;
					maxpos = pos;
					entered = true;
				}
			}
			if (entered) {
				oparameters.insert(oparameters.begin(),
						this->parameters[maxpos]);
				paramss.erase(paramss.begin() + maxi);
			}
		}
	}

	// sets  updateList to set which parameters will be updated when one parameter changed
	for (int i = 0; i < oparameters.size(); i++) {
		vector<int> depents;
		int pos = Utils::findPair(this->parameters, oparameters[i]);
		if (pos >= 0) {
			for (int k = i - 1; k >= 0; k--) {
				//updateList should be created in the order which parameter wanted to be updated first (with reverse order of paramsOrdered)
				//the beginning of paramsOrdered will depend on the later parameters.
				int pos2 = Utils::findPair(this->parameters, oparameters[k]);
				bool found = false;
				for (int j = 0; j < this->tetrahedra.size() && !found; j++) {
					if (tetraHedraParams[j][pos] && tetraHedraParams[j][pos2]) //if those two parameters exist in same tetrahedra
							{
						found = true;
						depents.push_back(k);
					}
				}
			}
		}
		updateList.push_back(depents);
	}

	/*
	 * While building 3-trees through tetrahedras, NOT all the edges of the tetrahedras are defined directly.
	 * directly means: through a contact or parameter or fixed distance by the edge being taking place in one helix.
	 * indirect means: the distance of the edge is computed through known 2 vertex location.
	 * If the distance for an edge is computed indirectly, then the range of the parameter that takes place in that tetrahedra
	 * cannot be computed by tetrahedra equalities in ConvexChart class, because the distance for dependent edge is not set yet.
	 * That distance will be computed during realization by CartesianRealizer class.
	 *
	 * example:
	 * suppose you have contacts 06 and 17  and parameters: 07, 18, 27, 28  tetrahedras: 0127, 1278, 0786
	 * so first the vertices 0,1,2,7 are located by the base tetrahedta 0127
	 * then vertex 8 is located that is connected to the face 127
	 * then vertex 6 is located that is connected to the face 078
	 * however face 078 is not complete, edge 08 is not there.
	 * The distance for that edge 08 needs to be computed by the locations of known two vertices 0 and 8
	 *
	 * In order to compute the distance interval for parameter 07 by tetrahedral equality, we need to know edge 08 first.
	 * However at the time of convex chart computations, the cartesian coordinates are not realized yet.
	 * So that distance 08 will not be able to computed. As a consequence, the range for the parameter 07 will not be as tight.
	 */

	// sets boundaryComputationWay
	for (int i = 0; i < oparameters.size(); i++) {
		int pos = Utils::findPair(this->contacts, oparameters[i]);
		if (pos >= 0)
			boundaryComputationWay.push_back(-2);  //sum of radius +- epsilon
		else {
			pos = Utils::findPair(this->parameters, oparameters[i]);
			if (pos >= 0) {
				bool first = true; //to find out first tetrahedra it takes place
				for (int j = 0; j < this->tetrahedra.size() && first; j++) {
					if (tetraHedraParams[j][pos]) {
						if (first) {
							first = false;
							int smallest = i; //pos;
							for (int k = 0; k < this->parameters.size(); k++) {
								if (tetraHedraParams[j][k]) {
									int before = Utils::findPair(oparameters,
											this->parameters[k]);
									if (before < smallest)
										smallest = before;
								}
							}

							bool containsFixedEdge = false;
							for (int m = 0; m < 4; m++)
								for (int n = 0; n < 4; n++)
									if (fixedEdge[this->tetrahedra[j][m]][this->tetrahedra[j][n]])
										containsFixedEdge = true;

							if (!containsFixedEdge && smallest == i) //if it is the last(first) parameter in the first tetra it takes place, then tetrahedra
								boundaryComputationWay.push_back(j); //tetrahedra,  also shows which tetrahedra to use.
							else
								boundaryComputationWay.push_back(-1); //triangular
						}
					}
				}
			}

		}
	}

	this->parameters = oparameters;
}

