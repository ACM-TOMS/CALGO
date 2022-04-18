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
#ifndef RATIOCHECKER_H_
#define RATIOCHECKER_H_

#include <iostream>
#include <vector>
#include <string>

// #ifdef __APPLE__
// 	#include <ext/hash_map>
// #else
// 	#include <hash_map>
// #endif
#include <unordered_map>

#include <cassert>

#include "AtlasBuilder.h"

#include "Atlas.h"
#include "CayleySpaceView.h"
#include "SweepView.h"
#include "SaveLoader.h"
#include "Settings.h"

#include "Cull.h"
#include "md5.h"

#include <fstream>

// using namespace __gnu_cxx;

using Eigen::Vector3d;

using namespace std;

class DistinctRealization {
public:

	static vector<double> cal_di(PointSet *hel, PredefinedInteractions *df,
			const DistinctRealization &r) {
		vector<double> res;

		for (PredefinedInteractions::dist_iterator dit = df->dist2begin();
				dit != df->dist2end(); dit++) {
			Vector a = Vector(hel->index.at(dit->first.first)->getLocation());
			Vector b = Vector(hel->index.at(dit->first.second)->getLocation());
			double dis = sqrt((a - b.trans(r.trans_mat)).squared_length());
			res.push_back(dis);
		}
		return res;
	}

	static string cal_flip_string(PointSet *hel, PredefinedInteractions *df,
			const DistinctRealization &r) {
		//vector<int> res;
		stringstream res;
		vector<Point*> atoms;
		vector<Vector> points;
		int count = 0;
		for (PredefinedInteractions::dist_iterator dit = df->dist2begin();
				dit != df->dist2end() && count < 10; dit++) {
			count++;
			Point *atom = hel->index.at(dit->first.first);
			if (find(atoms.begin(), atoms.end(), atom) == atoms.end()) {
				atoms.push_back(atom);
				points.push_back(
						Vector(atom->getLocation()).trans(r.trans_mat));
			}
			atom = hel->index.at(dit->first.second);

			if (find(atoms.begin(), atoms.end(), atom) == atoms.end()) {
				atoms.push_back(atom);
				points.push_back(Vector(atom->getLocation()));
			}
		}

		count = 0;

		for (int i = 0; i < points.size(); i++) {
			for (int j = i + 1; j < points.size(); j++) {
				for (int k = j + 1; k < points.size(); k++) {
					Vector v1 = (points[j] - points[i]).normalize();
					Vector v2 = (points[k] - points[i]).normalize();
					Vector3d vd1, vd2;
					vd1 << v1[0], v1[1], v1[2];
					vd2 << v2[0], v2[1], v2[2];
					for (int l = k + 1; l < points.size(); l++) {
						Vector v3 = (points[l] - points[i]).normalize();
						Vector3d vd3;
						vd3 << v3[0], v3[1], v3[2];
						double resd = vd1.cross(vd2).dot(vd3);
						//cout << resd << endl;
						if (resd > 1e-4) {
							//res.push_back(1);
							res << '1';
						} else if (resd < -1e-4) {
							//res.push_back(0);
							res << '0';
						} else {
							//res.push_back(2);
							res << '2';
						}
					}
				}
			}
		}
		return md5(res.str());
	}
	friend ostream& operator<<(ostream& os, const DistinctRealization& dr);
public:
	vector<double> trans_mat;
	int node, spacenum;
	PointSet *hel;
	vector<double> di;
	//vector<int> flip_string;
	string flip_string;
	int num_collision;
	ActiveConstraintGraph* cgk;
	DistinctRealization(PointSet *hel, PredefinedInteractions *df,
			double fa[][3], double ta[][3], double fb[][3], double tb[][3],
			int node, int spacenum, string* fstring = NULL) {
		trans_mat = Utils::getSingleTransformFromTwo(
				Utils::getTransMatrix(fa, ta), Utils::getTransMatrix(fb, tb));
		this->node = node;
		this->spacenum = spacenum;
		this->hel = hel;
		di = cal_di(hel, df, *this);
		if (fstring == NULL)
			flip_string = cal_flip_string(hel, df, *this);
		else
			flip_string = *fstring;
		//num_collision = count_collision(hel->stree,hel->stree,fa,ta,fb,tb,0.4);
		num_collision = 0;
	}
	DistinctRealization(PointSet *hel, PredefinedInteractions *df,
			const vector<double> &trans_mat, int node, int spacenum,
			string* fstring = NULL) {
		this->trans_mat = trans_mat;
		this->node = node;
		this->spacenum = spacenum;
		this->hel = hel;
		di = cal_di(hel, df, *this);
		if (fstring == NULL)
			flip_string = cal_flip_string(hel, df, *this);
		else
			flip_string = *fstring;
		//num_collision = count_collision(hel->stree,hel->stree,fa,ta,fb,tb,0.4);
		num_collision = 0;
	}

};

class RealizationCompare {
	PointSet *hel;
	vector<Point*> atoms;
	double threshhold;
public:
	RealizationCompare(PointSet* hel, int threshhold = 3) {
		this->hel = hel;
		this->threshhold = threshhold;
		cal_convex_hull();
	}

	void cal_convex_hull() {
		atoms.clear();
		Chull3D chull(hel);
		chull.compute();
		chull.get_atoms(&atoms);
	}

	double difference(const DistinctRealization &r1,
			const DistinctRealization &r2) const {
		double max_dis = 0;
#if 0
		for(map<string, Point*>::iterator iter = hel->index.begin();
				iter != hel->index.end(); iter++) {
			Vector pos = Vector(iter->second->getLocation());
			double dis = sqrt((pos.trans(r1.trans_mat) - pos.trans(r2.trans_mat)).squared_length());
			if(dis > max_dis)
			max_dis = dis;
		}
#endif
#if 0
		for(vector<Point*>::iterator iter = atoms.begin(); iter != atoms.end(); iter++) {
			Vector pos = Vector((*iter)->getLocation());
			double dis = sqrt((pos.trans(r1.trans_mat) - pos.trans(r2.trans_mat)).squared_length());
			if(dis > max_dis)
			max_dis = dis;
		}
#endif
		for (int i = 0; i < r1.di.size(); i++) {
			if (fabs(r1.di[i] - r2.di[i]) > max_dis) {
				max_dis = fabs(r1.di[i] - r2.di[i]);
			}
		}
		return max_dis;
	}

	bool diff_flip(const DistinctRealization &r1,
			const DistinctRealization &r2) const {
		if (r1.flip_string.size() != r2.flip_string.size()) {
			cout << r1.flip_string << endl;
			cout << r2.flip_string << endl;
		}
		assert(r1.flip_string.size() == r2.flip_string.size());

		for (int i = 0; i < r1.flip_string.size(); i++) {
			if (r1.flip_string[i] != r2.flip_string[i])
				return false;
		}
		return true;
	}

	bool operator()(const DistinctRealization &r1,
			const DistinctRealization &r2) const {
		double max_dis = difference(r1, r2);
		return diff_flip(r1, r2) && max_dis < 0.1;

#if 0
		if(max_dis < threshhold)
		return false;
		else
		return r1.trans_mat[0] < r2.trans_mat[0];
#endif

	}
};

class RealizationHash {
	PointSet *hel;
	const DistinctRealization& refer;
	RealizationCompare rc;
public:
	RealizationHash(PointSet *hel, DistinctRealization &re) :
			refer(re), rc(hel) {
		this->hel = hel;
	}

	size_t operator()(const DistinctRealization &r1) const {

		//RealizationCompare* prc = (RealizationCompare*)&rc;
		//double dis = prc->difference(r1,refer);
		//return floor(dis*5);
		int res = 0;
		for (int i = 0; i < r1.flip_string.size(); i++)
			res += (r1.flip_string[i] - '0');
		return res;
//		double len = 0;
//		for(int i = 0; i < r1.di.size(); i++){
//			len += (r1.di[i])*r1.di[i];
//		}
//		return floor(sqrt(len)*3);
	}
};

class VectorCompare {
public:
	bool operator()(const vector<double> &v1, const vector<double>&v2) {
#if 0
		//int size = min(v1.size(),v2.size());
		for(int i = 0; i < v1.size(); i++) {
			double err = (v1[i]-v2[i]);
			if(i%4 == 3 )
			err /= 60;
			if(err < -0.2)
			return true;
			if(err > 0.2)
			return false;
		}
#endif
		return false;
	}

};

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/LU"

using namespace std;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::Quaterniond;

class RatioChecker {
public:
	RatioChecker(SaveLoader* snl, PointSet *hA, PointSet *hB, Atlas *mapview);
	virtual ~RatioChecker();
	void checkRatio();
	PredefinedInteractions* getDistanceFinderFromFile(const char* filename);
private:

	typedef unordered_map<DistinctRealization, int, RealizationHash,
			RealizationCompare> RatioBookkeeper;
	RatioBookkeeper *bookkeeper;

	SaveLoader* snlr;
	PointSet *helxA;
	PointSet *helxB;
	Atlas* mapView;
	vector<vector<pair<int, int> > > readPairs(const char *filename);
	void preComputeCache(DistinctRealization & correct_dr,
			vector<AtlasNode*> nodes, PredefinedInteractions *df);

};

#endif /* RATIOCHECKER_H_ */
