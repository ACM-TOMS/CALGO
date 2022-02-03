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
 * CartesianRealizer.cpp
 *
 *  Created on: 2008-2014
 *      Author: Aysegul Ozkan
 */

#include "CartesianRealizer.h"

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

using namespace std;

#define PI 3.14159265

bool CartesianRealizer::partial3tree = true;

vector<vector<int> > CartesianRealizer::tetrahedra;

Vector3d CartesianRealizer::TB;
Vector3d CartesianRealizer::eaB;

bool CartesianRealizer::debug = false;

bool CartesianRealizer::realizable = true;

PointSet* CartesianRealizer::helA;
PointSet* CartesianRealizer::helB;

vector<int> CartesianRealizer::verticesA;
vector<int> CartesianRealizer::verticesB;

ActiveConstraintGraph* CartesianRealizer::graph;
int CartesianRealizer::flipNo = -1;

double * CartesianRealizer::fromA[3];
double * CartesianRealizer::fromB[3];
double * CartesianRealizer::toA[3];
double * CartesianRealizer::toB[3];

double CartesianRealizer::positions[NO][3];     // locations
double CartesianRealizer::edge_length[NO][NO];

Orientation* CartesianRealizer::computeRealization(ActiveConstraintGraph* cgk,
		ConvexChart* descr, int flipno, bool& fail) {

	return computeRealization(cgk, descr, flipno, descr->param_length, fail);
}

Orientation* CartesianRealizer::computeRealization(ActiveConstraintGraph* cgk,
		ConvexChart* descr, int flipno, double paramconnections[][12],
		bool& fail) {
	fail = false;

	flipNo = flipno;

	for (int i = 0; i < NO; i++)
		for (int j = 0; j < NO; j++)
			edge_length[i][j] = paramconnections[i][j];

	graph = cgk;

	helA = graph->getMolecularUnitA();
	helB = graph->getMolecularUnitB();

	verticesA = graph->getVerticesA();
	verticesB = graph->getVerticesB();

	fromA[0] = helA->getAtomAt(verticesA[0])->getLocation();
	fromA[1] = helA->getAtomAt(verticesA[1])->getLocation();
	fromA[2] = helA->getAtomAt(verticesA[2])->getLocation();

	fromB[0] = helB->getAtomAt(verticesB[0])->getLocation();
	fromB[1] = helB->getAtomAt(verticesB[1])->getLocation();
	fromB[2] = helB->getAtomAt(verticesB[2])->getLocation();
	realizable = true;

	for (int i = 0; i < NO; i++) {
		positions[i][0] = -111;
		positions[i][1] = -111;
		positions[i][2] = -111;
	}

	vector<bool> mirrorlist;
	mirrorlist.push_back(flipno % 2 == 1);
	mirrorlist.push_back((flipno / 2) % 2 == 1);
	mirrorlist.push_back((flipno / 4) % 2 == 1);

	for (int i = 0; i < 3; i++) {
		toA[i] = new double[3];
		toB[i] = new double[3];
	}

	partial3tree = descr->partial3tree;
	if (partial3tree) {
		tetrahedra = descr->getTetras();
		size_t mirrorIndex = 1;
		if (setBaseTetra(mirrorlist[0])) //if indices are in the same helix , then no need for mirror
			mirrorIndex = 0;
		if (!realizable) {
			fail = true;
			destructor();
			return NULL;
		}

		//builds new locations onto the base tetrahedron.
		size_t tetraIndex = 1;
		do {
			if (!locateVertex(tetraIndex, mirrorlist[mirrorIndex])) //if locateVertex is true, then indices are in the same helix, then no need for mirror and volume is positive
					{
				if (!realizable) {
					fail = true;
					destructor();
					return NULL;
				}
				computeLengthOfTheEdgesConnectedToVertex(
						tetrahedra[tetraIndex][0]);
				mirrorIndex++;

			}
			tetraIndex++;
		} while (tetraIndex < tetrahedra.size());
		//now it should be able to fill all the values of vertices, else it means it needs maple

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				toA[i][j] = positions[i][j];
				toB[i][j] = positions[i + NO / 2][j];
				if (toA[i][j] == -111 || toB[i][j] == -111
						|| std::isnan(toA[i][j]) || std::isnan(toB[i][j]))
					cout << "eroorrrr toA[i][j]== -111 || toB[i][j]== -111"
							<< endl;
			}
		}

		if (isDistorted()) {
			realizable = false;
			fail = true;
			destructor();
			return NULL;
		}

		// TRANSLATE HELIX B ACCORDING TO HELIX A
		for (int j = 0; j < 3; j++) { //tranlates toB[j]
			double pot[3];
			Utils::matApp(toB[j], toA[0], toA[1], toA[2], fromA[0], fromA[1],
					fromA[2], pot);

			for (int i = 0; i < 3; i++) {
				toB[j][i] = pot[i];
			}
		}

		// FIX HELIX A
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				toA[i][j] = fromA[i][j];
	}

	Orientation * ori_output = getOrienation();
	destructor();
	return ori_output;
}

void CartesianRealizer::destructor() {
	for (int i = 0; i < 3; i++) {
		delete[] toA[i];
		delete[] toB[i];
	}
}

// for better result, additionally we can check the transformation from fromB to toB as well
bool CartesianRealizer::isDistorted() {
	//creates 2 points p1 and p2 where the distance between them is exactly 1
	double *p1 = new double[3];
	p1[0] = 0;
	p1[1] = 0;
	p1[2] = 0;
	double *p2 = new double[3];
	p2[0] = 0;
	p2[1] = 0;
	p2[2] = 1;

	//compute pt1 and pt2, the transformed location of p1 and p2 using the same transformation from fromA to toA
	double pt1[3];
	double pt2[3];
	Utils::matApp(p1, toA[0], toA[1], toA[2], fromA[0], fromA[1], fromA[2],
			pt1);
	Utils::matApp(p2, toA[0], toA[1], toA[2], fromA[0], fromA[1], fromA[2],
			pt2);

	//compute the distance between pt1 and pt2
	double dist = Utils::dist(pt1, pt2);
	delete[] p1;
	delete[] p2;

	//the distance should approximate to 1, otherwise there is a distortion
	if (dist < .9 || dist > 1.1) {
		return true;
	}
	return false;
}

bool CartesianRealizer::locateVertex(int tetrahedron_index, bool mirror) {

	double posi[4][3];

	vector<int> ind = tetrahedra[tetrahedron_index];
	printPositionsAndConnections(ind);

	if (ind[0] < NO / 2 && ind[1] < NO / 2 && ind[2] < NO / 2
			&& ind[3] < NO / 2) //inside of same helix i.e. helix A
					{
		double *from[4];
		for (int i = 0; i < 4; i++) {
			from[i] = helA->getAtomAt(verticesA[ind[i]])->getLocation();
		}
		double pot[3];
		Utils::matApp(from[0], from[1], from[2], from[3], positions[ind[1]],
				positions[ind[2]], positions[ind[3]], pot);

		for (int i = 0; i < 3; i++)
			positions[ind[0]][i] = pot[i];

		return true;
	} else if (ind[0] >= NO / 2 && ind[1] >= NO / 2 && ind[2] >= NO / 2
			&& ind[3] >= NO / 2)  //inside of same helix i.e. helix B
					{
		double *from[4];
		for (int i = 0; i < 4; i++) {
			from[i] =
					helB->getAtomAt(verticesB[ind[i] - NO / 2])->getLocation();
		}
		double pot[3];
		Utils::matApp(from[0], from[1], from[2], from[3], positions[ind[1]],
				positions[ind[2]], positions[ind[3]], pot);

		for (int i = 0; i < 3; i++)
			positions[ind[0]][i] = pot[i];

		return true;
	}

	locateTetrahedron(tetrahedron_index, posi, mirror);
	if (!realizable)
		return false;

	double pot[3];
	Utils::matApp(posi[0], posi[1], posi[2], posi[3], positions[ind[1]],
			positions[ind[2]], positions[ind[3]], pot);

	for (int i = 0; i < 3; i++) {
		positions[ind[0]][i] = pot[i];
		if (debug)
			cout << positions[ind[0]][i] << endl;
	}

	if (debug)
		cout << ind[0] << "-" << ind[1] << "   " << ind[0] << "-" << ind[2]
				<< "   " << ind[0] << "-" << ind[3] << "   " << endl;
	if (debug)
		cout << "fact " << edge_length[ind[0]][ind[1]] << " "
				<< edge_length[ind[0]][ind[2]] << " "
				<< edge_length[ind[0]][ind[3]] << endl;
	if (debug)
		cout << "before " << Utils::dist(posi[0], posi[1]) << " "
				<< Utils::dist(posi[0], posi[2]) << " "
				<< Utils::dist(posi[0], posi[3]) << endl;
	if (debug)
		cout << "later " << Utils::dist(positions[ind[0]], positions[ind[1]])
				<< " " << Utils::dist(positions[ind[0]], positions[ind[2]])
				<< " " << Utils::dist(positions[ind[0]], positions[ind[3]])
				<< endl;

	bool foundNan = false; //if volume is close to 0, even if locateTetrahedron gives correct numbers after matApp, they become nan. matApp doesnot work for volume 0 case. It cannot compute translation matrix for this case!
	for (int i = 0; i < 4; i++)
		if (std::isnan(positions[ind[i]][0]) || std::isnan(positions[ind[i]][1])
				|| std::isnan(positions[ind[i]][2]))
			foundNan = true;

	bool disturbed = false;
	for (int i = 0; i < 4; i++)
		for (int j = i + 1; j < 4; j++)
			if (abs(
					Utils::dist(posi[i], posi[j])
							- Utils::dist(positions[ind[i]], positions[ind[j]]))
					> .6)
				disturbed = true;

	if (disturbed || foundNan) {

		if (debug)
			for (int i = 1; i < 4; i++) {
				cout << ind[i] << ": " << positions[ind[i]][0] << " "
						<< positions[ind[i]][1] << " " << positions[ind[i]][2]
						<< endl;
				for (int j = i + 1; j < 4; j++)
					cout << ind[i] << "-" << ind[j] << ": "
							<< edge_length[ind[i]][ind[j]] << " "
							<< Utils::dist(positions[ind[i]], positions[ind[j]])
							<< endl;
			}

		double * pot = Utils::find4thPoint_inTheSamePlane(positions[ind[1]],
				positions[ind[2]], positions[ind[3]],
				edge_length[ind[0]][ind[1]], edge_length[ind[0]][ind[2]],
				edge_length[ind[0]][ind[3]]);

		for (int i = 0; i < 3; i++)
			positions[ind[0]][i] = pot[i];

		if (debug)
			cout << ind[0] << ": " << positions[ind[0]][0] << " "
					<< positions[ind[0]][1] << " " << positions[ind[0]][2]
					<< endl;

		delete[] pot;

		if (std::isnan(positions[ind[0]][0]) || std::isnan(positions[ind[0]][1])
				|| std::isnan(positions[ind[0]][2])) //if volume is close to 0, even if locateTetrahedron gives correct numbers after matApp, they become nan. matApp doesnot work for volume 0 case. It cannot compute translation matrix for this case!
			realizable = false;

		if (debug)
			cout << ind[0] << "-" << ind[1] << "   " << ind[0] << "-" << ind[2]
					<< "   " << ind[0] << "-" << ind[3] << "   " << endl;
		if (debug)
			cout << "fact " << edge_length[ind[0]][ind[1]] << " "
					<< edge_length[ind[0]][ind[2]] << " "
					<< edge_length[ind[0]][ind[3]] << endl;
		if (debug)
			cout << "later "
					<< Utils::dist(positions[ind[0]], positions[ind[1]]) << " "
					<< Utils::dist(positions[ind[0]], positions[ind[2]]) << " "
					<< Utils::dist(positions[ind[0]], positions[ind[3]])
					<< endl;

		if (abs(
				edge_length[ind[0]][ind[1]]
						- Utils::dist(positions[ind[0]], positions[ind[1]]))
				> .6)
			realizable = false;
		if (abs(
				edge_length[ind[0]][ind[2]]
						- Utils::dist(positions[ind[0]], positions[ind[2]]))
				> .6)
			realizable = false;
		if (abs(
				edge_length[ind[0]][ind[3]]
						- Utils::dist(positions[ind[0]], positions[ind[3]]))
				> .6)
			realizable = false;

	}

	return false;
}

void CartesianRealizer::computeLengthOfTheEdgesConnectedToVertex(int vertex) {

	for (int i = 0; i < NO; i++)
		if (positions[i][0] != -111 && edge_length[i][vertex] == -1) //if ith vertex is located before and if the distance i-vertex is not defined before
				{
			edge_length[i][vertex] = edge_length[vertex][i] = Utils::dist(
					positions[i], positions[vertex]);
			if (debug)
				cout << "edge_length[" << i << "][" << vertex << "] "
						<< edge_length[i][vertex] << endl;
		}
}

void CartesianRealizer::locateTetrahedron(int tetrahedron_index,
		double posi[][3], bool mirror) {
	double position[4][3];
	double edgeLengths[6];
	vector<int> inds = tetrahedra[tetrahedron_index];

	printPositionsAndConnections(inds);

	edgeLengths[0] = edge_length[inds[0]][inds[1]];
	edgeLengths[1] = edge_length[inds[2]][inds[3]];
	edgeLengths[2] = edge_length[inds[0]][inds[2]];
	edgeLengths[3] = edge_length[inds[1]][inds[3]];
	edgeLengths[4] = edge_length[inds[0]][inds[3]];
	edgeLengths[5] = edge_length[inds[1]][inds[2]];

	double vc = Utils::volumeTetra(edgeLengths);
	vc += 0.0001; //if vs= -..e^-...
	if (debug)
		cout << "vc2 " << vc << endl;

	if (vc < 0) {
		realizable = false;
		return;
	}

	bool foundNan = false;
	Utils::lenToTetra(edgeLengths, position[0], position[1], position[2],
			position[3], mirror);
	for (int i = 0; i < 4; i++) {
		posi[i][0] = position[i][0];
		posi[i][1] = position[i][1];
		posi[i][2] = position[i][2];
		if (std::isnan(posi[i][0]) || std::isnan(posi[i][1])
				|| std::isnan(posi[i][2]))
			foundNan = true;
	}

	bool disturbed = false;
	for (int i = 0; i < 4; i++)
		for (int j = i + 1; j < 4; j++)
			if (abs(
					edge_length[inds[i]][inds[j]]
							- Utils::dist(posi[i], posi[j])) > .6) {
				disturbed = true;
			}

	//if you enter here when volume not exactly zero, you may disturb your lengths !
	//so enter here if the lengths are already disturbed.
	if ((vc < 15 && vc >= 0 && disturbed) || foundNan) {
		posi[0][0] = 0;
		posi[0][1] = 0;
		posi[0][2] = 0;
		posi[1][0] = edge_length[inds[0]][inds[1]];
		posi[1][1] = 0;
		posi[1][2] = 0;

		double ac = edge_length[inds[0]][inds[2]], ab =
				edge_length[inds[0]][inds[1]], bc =
				edge_length[inds[1]][inds[2]];
		double cosCAB = (ac * ac + ab * ab - bc * bc) / (2 * ac * ab);
		double sinCAB = sqrt(abs(1 - cosCAB * cosCAB));
		posi[2][0] = edge_length[inds[0]][inds[2]] * cosCAB;
		posi[2][1] = edge_length[inds[0]][inds[2]] * sinCAB;
		posi[2][2] = 0;

		ac = edge_length[inds[0]][inds[3]];
		ab = edge_length[inds[0]][inds[1]];
		bc = edge_length[inds[1]][inds[3]];
		cosCAB = (ac * ac + ab * ab - bc * bc) / (2 * ac * ab);
		sinCAB = sqrt(abs(1 - cosCAB * cosCAB));
		posi[3][0] = edge_length[inds[0]][inds[3]] * cosCAB;
		posi[3][1] = edge_length[inds[0]][inds[3]] * sinCAB;
		posi[3][2] = 0;

		double nposi[2];
		nposi[0] = posi[3][0];
		nposi[1] = -1 * posi[3][1];
		double cd1 = sqrt(
				(posi[2][0] - posi[3][0]) * (posi[2][0] - posi[3][0])
						+ (posi[2][1] - posi[3][1])
								* (posi[2][1] - posi[3][1]));
		double cd2 = sqrt(
				(posi[2][0] - nposi[0]) * (posi[2][0] - nposi[0])
						+ (posi[2][1] - nposi[1]) * (posi[2][1] - nposi[1]));
		double cd = edge_length[inds[2]][inds[3]];

		if (abs(cd - cd1) > abs(cd - cd2))
			posi[3][0] = nposi[0];
		posi[3][1] = nposi[1];

	}
}

bool CartesianRealizer::setBaseTetra(bool mirror) {

	vector<int> inds = tetrahedra.front();
	if (inds[0] < NO / 2 && inds[1] < NO / 2 && inds[2] < NO / 2
			&& inds[3] < NO / 2) //inside of same helix
					{
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 3; j++)
				positions[inds[i]][j] =
						helA->getAtomAt(verticesA[inds[i]])->getLocation()[j];

		}
		return true;
	} else if (inds[0] >= NO / 2 && inds[1] >= NO / 2 && inds[2] >= NO / 2
			&& inds[3] >= NO / 2) //inside of same helix
					{
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 3; j++)
				positions[inds[i]][j] = helB->getAtomAt(
						verticesB[inds[i] - NO / 2])->getLocation()[j];

		}
		return true;
	}

	if (debug)
		cout << "base " << inds[0] << " " << inds[1] << " " << inds[2] << " "
				<< inds[3] << endl;

	double posi[4][3];
	locateTetrahedron(0, posi, mirror);
	for (int i = 0; i < 4; i++) {
		positions[inds[i]][0] = posi[i][0];   //x coordinate
		positions[inds[i]][1] = posi[i][1];   //y coordinate
		positions[inds[i]][2] = posi[i][2];   //z coordinate
	}

	return false;

}

void CartesianRealizer::compute_TB_EAB() {

	Vector3d p1(fromB[0][0], fromB[0][1], fromB[0][2]);
	Vector3d p2(fromB[1][0], fromB[1][1], fromB[1][2]);
	Vector3d p3(fromB[2][0], fromB[2][1], fromB[2][2]);

	Vector3d P1(toB[0][0], toB[0][1], toB[0][2]);
	Vector3d P2(toB[1][0], toB[1][1], toB[1][2]);
	Vector3d P3(toB[2][0], toB[2][1], toB[2][2]);

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

	Vector3d TB_ = t + mean;

	Vector3d eaB_ = Utils::RotMatToEuler(R);
	eaB_[0] = eaB_[0] * 180 / PI;
	eaB_[2] = eaB_[2] * 180 / PI;
	double cos_eaB1 = eaB_[1];

	TB = TB_;
	eaB = eaB_;
}

Orientation* CartesianRealizer::getOrienation() {

	Orientation* output = new Orientation(fromB, toB);
	output->setFlipNum(flipNo);

	return output;
}

void CartesianRealizer::printPositionsAndConnections(vector<int> ind) {

	if (debug)
		cout << "indices " << ind[0] << " " << ind[1] << " " << ind[2] << " "
				<< ind[3] << endl;

	if (debug)
		for (int i = 0; i < 4; i++) {
			cout << ind[i] << ": " << positions[ind[i]][0] << " "
					<< positions[ind[i]][1] << " " << positions[ind[i]][2]
					<< endl;
			for (int j = i + 1; j < 4; j++)
				cout << ind[i] << "-" << ind[j] << ": "
						<< edge_length[ind[i]][ind[j]] << " "
						<< Utils::dist(positions[ind[i]], positions[ind[j]])
						<< endl;
		}

	if (debug)
		cout << endl;
}
