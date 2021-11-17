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
 * RatioChecker.cpp
 *
 *  Created on: Nov 5, 2011
 *      Author: ruijin
 */

#include "RatioChecker.h"

using namespace std;

RatioChecker::RatioChecker(SaveLoader* snl, PointSet *hA, PointSet *hB,
		Atlas *mapview) :
		snlr(snl), helxA(hA), helxB(hB), mapView(mapview) {

}

RatioChecker::~RatioChecker() {
}

void RatioChecker::preComputeCache(DistinctRealization & correct_dr,
		vector<AtlasNode*> nodes, PredefinedInteractions *df) {
	// pre-compute all relizations
	RatioBookkeeper tempbk(300, RealizationHash(helxA, correct_dr),
			RealizationCompare(helxA));
	ofstream odrfile("reals.txt");

	vector<AtlasNode*>::iterator iter;
	for (iter = nodes.begin(); iter != nodes.end(); iter++) {
		if ((*iter)->getDim() == 0) //for specific strata with dim dimension
				{
			tempbk.clear();
			int node = (*iter)->getID();
			cout << node << endl;
			{
				ActiveConstraintGraph *cgK = (*iter)->getCG();
				ActiveConstraintRegion* region = (*iter)->getACR();

				snlr->loadNode(node, region);
				vector<CayleyPoint*> spc = (*iter)->getACR()->getSpace();
				vector<CayleyPoint*>::iterator its;
				int spacenum = 0;
				for (its = spc.begin(); its != spc.end(); its++) {
					spacenum++;
					vector<Orientation*> sol = (*its)->getOrientations();
					for (size_t a = 0; a < sol.size(); a++) {
						double fa[3][3], fb[3][3], ta[3][3], tb[3][3];
						DistinctRealization dr(this->helxA, df, fa, ta, fb, tb,
								node, spacenum);
						RatioBookkeeper::iterator res = tempbk.find(dr);
						if (res == tempbk.end()) {
							tempbk[dr] = 1;
						} else {
							res->second++;
							continue;
						}
					}

				}

				region->trim();
			}
			for (RatioBookkeeper::iterator iter = tempbk.begin();
					iter != tempbk.end(); iter++) {
				odrfile << iter->first << " " << iter->second << endl;
			}
		}

	}

	odrfile.close();
}

vector<vector<pair<int, int> > > RatioChecker::readPairs(const char* filename) {
	// read exclu_pairs from file
	vector<vector<pair<int, int> > > exclu_pairs;
	ifstream efs(filename);
	char buf[255];
	while (!efs.eof()) {
		int v1, v2;
		vector<pair<int, int> > pairs;
		efs.getline(buf, 255);
		bool line_stop = false;
		char *ptr = strtok(buf, " ");
		while (!line_stop && ptr) {
			int res = sscanf(ptr, "%d,%d", &v1, &v2);
			cout << res << endl;
			if (res == 2) {
				cout << v1 << "," << v2 << endl;
				pairs.push_back(make_pair(v1, v2));
			} else {
				line_stop = true;
			}
			ptr = strtok(NULL, " ");
		}
		if (pairs.size() != 0)
			exclu_pairs.push_back(pairs);
	}
	efs.close();
	return exclu_pairs;
}

void RatioChecker::checkRatio() {
	snlr->loadMapView(mapView);
	vector<AtlasNode*> nodes = mapView->getNodes();
	PredefinedInteractions *distfinder = getDistanceFinderFromFile("tt.txt");

	// Hardcode the real one!
	double correct_mat[] = { -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0 };
	vector<double> trans_mat = Utils::getMatrixFromFileFromTo("fromto.txt");
	ofstream corrfs("corr.mat");
	for (int i = 0; i < 12; i++) {
		corrfs << trans_mat[i] << " ";
	}
	corrfs.close();

	DistinctRealization correct_dr(this->helxA, distfinder,
			Utils::getMatrixFromFileFromTo("fromto.txt"), 0, 0);
	RealizationCompare rc(this->helxA, 0.1);
	bookkeeper = new RatioBookkeeper(300, RealizationHash(helxA, correct_dr),
			RealizationCompare(helxA));

	// read exclu_pairs from file
	vector<vector<pair<int, int> > > exclu_pairs = readPairs("pairs.txt");
	// pre-compute all relizations
	preComputeCache(correct_dr, nodes, distfinder);

	for (vector<vector<pair<int, int> > >::iterator pit = exclu_pairs.begin();
			pit != exclu_pairs.end(); pit++) {
		for (int i = 0; i < pit->size(); i++) {
			cout << "dropping constrain " << pit->at(i).first << "---"
					<< pit->at(i).second << endl;
		}
		bookkeeper->clear();
		ifstream idrfile("reals.txt");
		string line;
		int linum = 0;
		while (getline(idrfile, line)) {
			stringstream linestr(line);
			linum++;
			int inode;
			linestr >> inode;
			vector<pair<int, int> > parts =
					(*mapView)[inode]->getCG()->getParticipants();
			bool found = false;
			for (int i = 0; i < pit->size(); i++) {
				if (find(parts.begin(), parts.end(), pit->at(i))
						!= parts.end()) {
					found = true;
					break;
				}
			}

			if (!found) {
				int drcount = 0;
				string fstring;
				vector<double> mat(12);
				for (int i = 0; i < 12; i++) {
					linestr >> mat[i];
				}
				linestr >> fstring;
				linestr >> drcount;
				if (fstring.length() < 25) {
					cout << "[err]" << fstring << endl;
					cout << linum << endl;
					cout << inode << "," << drcount << endl;
				}
				DistinctRealization dr(this->helxA, distfinder, mat, inode, 0,
						&fstring);
				RatioBookkeeper::iterator res = bookkeeper->find(dr);
				if (res == bookkeeper->end()) {
					(*bookkeeper)[dr] = drcount;
				} else {
					res->second += drcount;
					continue;
				}
			} else {
			}
		}

		idrfile.close();
		stringstream sstr;
		sstr << "copies_" << pit - exclu_pairs.begin() << ".txt";
		ofstream copies_file(sstr.str().c_str());
		int all_distinct = 0;
		int all_copies = 0;
		int true_distinct = 0;
		int true_copies = 0;
		int num_collision = 0;
		for (RatioBookkeeper::iterator iter = bookkeeper->begin();
				iter != bookkeeper->end(); iter++) {
			// gather some infomation
			all_copies += iter->second;
			num_collision += iter->first.num_collision;
			// print to file
			copies_file << iter->first.node << "," << iter->first.spacenum
					<< "," << iter->second << ","
					<< rc.difference(iter->first, correct_dr);
			copies_file << ",[";
			for (int i = 0; i < iter->first.di.size(); i++) {
				copies_file << iter->first.di[i] << ",";
			}
			copies_file << "]" << endl;
			copies_file << "[mat] ";
			for (int i = 0; i < 12; i++) {
				copies_file << iter->first.trans_mat[i] << " ";
			}
			copies_file << endl;
			copies_file << "{";

			copies_file << iter->first.flip_string;
			copies_file << "}" << endl;
			copies_file << iter->first.num_collision << endl;
			//debug
			cout << "corr di" << endl;
			for (int i = 0; i < correct_dr.di.size(); i++) {
				cout << correct_dr.di[i] << " ";
			}
			cout << endl;
			// check if this is the correct one
			if (rc.difference(iter->first, correct_dr) < 0.4) {
				cout << "correct one ";
				cout << iter->first.node << "," << iter->first.spacenum << ","
						<< iter->second << ","
						<< rc.difference(iter->first, correct_dr) << ": ";
				for (int i = 0; i < iter->first.trans_mat.size(); i++)
					cout << iter->first.trans_mat.at(i) << ",";

				cout << endl;
				true_copies += iter->second;
				true_distinct++;
			}
			// check if this is satifying
			bool fit = true;
			for (PredefinedInteractions::dist_iterator dis =
					distfinder->dist2begin(); dis != distfinder->dist2end();
					dis++) {
				Point *ata = helxA->index[dis->first.first];
				Point *atb = helxB->index[dis->first.second];
				double max_dist = dis->second;
				double min_dist = distfinder->getDist1(ata, atb);
				//double min_dist = ata->getMinDist(atb);
				double real_dist =
						sqrt(
								(Vector(ata->getLocation())
										- Vector(atb->getLocation()).trans(
												iter->first.trans_mat)).squared_length());
				cout << real_dist << "," << max_dist << "," << min_dist << endl;
				if (real_dist > max_dist || real_dist < min_dist)
					fit = false;

			}
			if (fit) {
				all_distinct++;
				cout << "valid" << endl;
				copies_file << "valid" << endl;
			}
		}

		copies_file.close();
		ofstream ratio_file("../ratio_stat.txt", ios_base::app);
		ratio_file << all_copies << ","; // all zero-dim realizations
		ratio_file << bookkeeper->size() << ","; // all distinct realizations
		ratio_file << all_distinct << ","; // all vaild zero-dim realizations
		ratio_file << true_copies << ","; // all true zero-dim realizations
		ratio_file << true_distinct << ","; // all true distinct zero-dim realizations
		ratio_file << (double) ((num_collision)) / all_distinct << endl;
		ratio_file.close();
		cout << "all zero-dim realizations " << all_copies << endl; // all zero-dim realizations
		cout << "all distinct realizations " << bookkeeper->size() << endl; // all distinct realizations
		cout << "all vaild zero-dim realizations " << all_distinct << endl; // all vaild zero-dim realizations
		cout << "all true zero-dim realizations " << true_copies << endl; // all true zero-dim realizations
		cout << "all true distinct zero-dim realizations " << true_distinct
				<< endl; // all true distinct zero-dim realizations
		cout << "num of everage collisions"
				<< (double) ((num_collision)) / all_distinct << endl;

	}
	delete bookkeeper;
	delete distfinder;
}

PredefinedInteractions* RatioChecker::getDistanceFinderFromFile(
		const char* filename) {
	ifstream ifs(filename, ifstream::in);
	vector<vector<string> > strings;
	while (ifs.good()) {
		vector<string> row;
		string l1, l2, dist;
		ifs >> l1 >> l2 >> dist;
		if (l1 == "")
			continue;
		row.push_back(l1);
		row.push_back(l2);
		row.push_back(dist);
		strings.push_back(row);
	}
	ifs.close();

	return new PredefinedInteractions(strings, 0, 1, 2, true);
}
