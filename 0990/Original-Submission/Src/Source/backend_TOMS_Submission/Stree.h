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
#ifndef _STRESS_H_
#define _STRESS_H_

#include <cstdlib>
#include <list>
#include <vector>
#include <iostream>
#include "Point.h"
#include "Utils.h"

#include "miniball.h"
#include <Eigen/Core>
#include <Eigen/LU>

using namespace std;
//const double PI = 3.1415926535987f;
#define NUM_BRANCH 4

class snode {
private:
	Vector random_within_sphere(Vector &center, double radius) {
		double x, y, z;
		x = (rand() % 100) / 100.0 * radius + center.x();
		y = (rand() % 100) / 100.0 * radius + center.y();
		z = (rand() % 100) / 100.0 * radius + center.z();
		return Vector(x, y, z);
	}

public:
	Vector center;
	double radius;
	Point* atom;

	std::list<snode*> children;

	snode(Point* atom) :
			center(atom->getLocation()), radius(atom->getRadius()) {
		// self contain, careful when clean up;
		children.push_back(this);
		this->atom = atom;

	}

	snode(std::vector<Point*> &atoms) {
		this->atom = NULL;

		// if the array have only one element
		// create a leaf node
		if (atoms.size() == 1) {
			center = Vector(atoms[0]->getLocation());
			this->radius = atoms[0]->getRadius();
			this->atom = atoms[0];
			children.push_back(this);
			return;
		}

		// create a sphere including all the small spheres
		miniball::Miniball<3> mb;
		miniball::Point<3> p;
		double max_rad = -1;
		for (std::vector<Point*>::iterator iter = atoms.begin();
				iter != atoms.end(); iter++) {
			for (int i = 0; i < 3; i++) {
				p[i] = (*iter)->getLocation()[i];
			}

			// this could be improve to make a tighter sphere
			if ((*iter)->getRadius() > max_rad) {
				max_rad = (*iter)->getRadius();
			}
			mb.check_in(p);
		}
		mb.build();
		this->center = Vector(mb.center()[0], mb.center()[1], mb.center()[2]);
		this->radius = sqrt(mb.squared_radius()) + max_rad;

		//TODO test only
#if 0
		if(this->radius <= 2*max_rad) {
			children.push_back(this);
			return;
		}
#endif
		if (atoms.size() <= NUM_BRANCH) {
			for (int i = 0; i < atoms.size(); i++) {
				children.push_back(new snode(atoms[i]));
			}
			return;
		}

		// K-Mean algorithm to do clustering
		int branchs = NUM_BRANCH;
		std::vector<Vector> v(branchs);
		std::vector<Vector> ov(branchs);
		double err = 100;

		for (int i = 0; i < branchs; i++) {
			v[i] = random_within_sphere(this->center, this->radius);
		}

		int step = 100;
		while (err > 0.1 || step < 0) { // not converge
			step--;
			std::vector<Vector> temp(branchs);	// new centers
			std::vector<int> counts(branchs);
			for (int i = 0; i < branchs; i++) {
				ov[i] = v[i];
				temp[i] = Vector(0, 0, 0);
				counts[i] = 0;
			}

			for (std::vector<Point*>::iterator iter = atoms.begin();
					iter != atoms.end(); iter++) {
				int minidx = 0;
				double mindis = 1000000.f;
				for (int i = 0; i < branchs; i++) {
					double dis =
							(v[i] - Vector((*iter)->getLocation())).squared_length();
					if (dis < mindis) {
						mindis = dis;
						minidx = i;
					}
				}
				temp[minidx] = temp[minidx] + Vector((*iter)->getLocation());
				counts[minidx]++;
			}
			err = 0;
			for (int i = 0; i < branchs; i++) {
				if (counts[i] > 0) {
					v[i] = temp[i] * (1.0 / counts[i]);
					double dis = (v[i] - ov[i]).squared_length();
					if (dis > err)
						err = dis;
				}
			}
		}

		std::vector<std::vector<Point*> > ps(branchs);
		for (std::vector<Point*>::iterator iter = atoms.begin();
				iter != atoms.end(); iter++) {
			int minidx = 0;
			double mindis = 1000000.f;
			for (int i = 0; i < branchs; i++) {
				double dis =
						(v[i] - Vector((*iter)->getLocation())).squared_length();
				if (dis < mindis) {
					mindis = dis;
					minidx = i;
				}
			}
			ps[minidx].push_back(*iter);
		}

		for (int i = 0; i < branchs; i++) {
			if (ps[i].size() == 0)
				continue;
			else
				children.push_back(new snode(ps[i]));

		}
	}

#if 1
	bool is_leaf() {
		return this->atom != NULL;
	}
#endif

	void clear() {
		if (!is_leaf()) {
			for (std::list<snode*>::iterator iter = children.begin();
					iter != children.end(); iter++) {
				(*iter)->clear();
				delete *iter;
			}
		}
		children.clear();
	}

};

class Stree {

public:
	Stree(std::vector<Point*> &pts) :
			root(pts) {
	}

	~Stree() {
		root.clear();
	}
	snode root;
};

bool is_collide(Stree *t1, Stree *t2, double fa[][3], double ta[][3],
		double fb[][3], double tb[][3], double threshold = 0);
bool is_collide(Stree *t1, Stree *t2, std::vector<double> &trans_a,
		std::vector<double> &trans_b, double threshold = 0);
int count_collision(Stree *t1, Stree *t2, double fa[][3], double ta[][3],
		double fb[][3], double tb[][3], double threshold = 0);
int count_collision(Stree *t1, Stree *t2, std::vector<double> &trans_a,
		std::vector<double> &trans_b, double threshold = 0);

#endif
