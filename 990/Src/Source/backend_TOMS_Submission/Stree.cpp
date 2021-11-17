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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <math.h>

#include <vector>
#include <set>

#include "Stree.h"
#include <string>

using namespace std;

using Eigen::Vector3f;
using Eigen::Matrix3f;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::MatrixXf;
using Eigen::VectorXf;

int count_test = 0;
set<string> name_set;

bool is_snode_collide(snode* node1, snode* node2, vector<double>& ta,
		vector<double>& tb, double treshold, bool count_mode) {

	if (sqrt(
			(node1->center.trans(ta) - node2->center.trans(tb)).squared_length())
			> node1->radius + node2->radius + treshold) {
		return false;
	} else if (node1->is_leaf() && node2->is_leaf()) {
		count_test++;
		name_set.insert(node1->atom->getName());
		return true;
	} else {
		bool flag = false;
		for (list<snode*>::iterator iter1 = node1->children.begin();
				iter1 != node1->children.end(); iter1++) {
			for (list<snode*>::iterator iter2 = node2->children.begin();
					iter2 != node2->children.end(); iter2++) {

				if (is_snode_collide(*iter1, *iter2, ta, tb, treshold,
						count_mode)) {
					if (count_mode) {
						flag = true;
					} else
						return true;
				}
			}
		}
		return flag;
	}
}

// Only for test
void testMat(double fa[][3], double ta[][3]) {
	vector<double> trans_a = Utils::getTransMatrix(fa, ta);
	double o[3];
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 3; j++) {
			o[j] = rand() % 300 / 50.0;
		}

		double no[3];
		cout << "o " << o[0] << "," << o[1] << "," << o[2] << endl;
		Utils::matApp(o, fa[0], fa[1], fa[2], ta[0], ta[1], ta[2], no);
		cout << "no " << no[0] << "," << no[1] << "," << no[2] << endl;
		Vector v(o);
		v = v.trans(trans_a);
		cout << "v " << v.x() << "," << v.y() << "," << v.z() << endl;
		double err = 0;
		for (int j = 0; j < 3; j++) {
			err += (no[j] - v[j]) * (no[j] - v[j]);
		}
		if (err > 0.001)
			cout << "not pass" << err << endl;
		else
			cout << "pass" << endl;
	}
}

bool is_collide(Stree *t1, Stree *t2, double fa[][3], double ta[][3],
		double fb[][3], double tb[][3], double threshold) {
	vector<double> trans_a = Utils::getTransMatrix(fa, ta);
	vector<double> trans_b = Utils::getTransMatrix(fb, tb);
	count_test = 0;
	bool res = is_snode_collide(&t1->root, &t2->root, trans_a, trans_b,
			threshold, false);
	return res;

}

bool is_collide(Stree *t1, Stree *t2, vector<double> &trans_a,
		vector<double> &trans_b, double threshold) {
	count_test = 0;
	bool res = is_snode_collide(&t1->root, &t2->root, trans_a, trans_b,
			threshold, false);
	return res;

}

int count_collision(Stree *t1, Stree *t2, double fa[][3], double ta[][3],
		double fb[][3], double tb[][3], double threshold) {
	vector<double> trans_a = Utils::getTransMatrix(fa, ta);
	vector<double> trans_b = Utils::getTransMatrix(fb, tb);
	count_test = 0;
	name_set.clear();
	bool res = is_snode_collide(&t1->root, &t2->root, trans_a, trans_b,
			threshold, true);
	return name_set.size();

}

int count_collision(Stree *t1, Stree *t2, vector<double> &trans_a,
		vector<double> &trans_b, double threshold) {
	count_test = 0;
	bool res = is_snode_collide(&t1->root, &t2->root, trans_a, trans_b,
			threshold, true);
	return count_test;

}
