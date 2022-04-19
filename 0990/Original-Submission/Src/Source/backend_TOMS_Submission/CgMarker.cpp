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
 * CgMarker.cpp
 *
 *  Created on: Jun 6, 2011
 *      Author: ruijin
 */

#include "CgMarker.h"

#include "PointSet.h"
#include "ActiveConstraintRegion.h"

#include <algorithm>
#include <list>
#include <iterator>
using namespace std;

CgMarker::CgMarker() :
		cgs(CGcomp) {

}

CgMarker::~CgMarker() {
}

void CgMarker::mark(AtlasNode* rnode) {
	// make sure the space is not empty
	ActiveConstraintGraph* cg = rnode->getCG();
	ActiveConstraintRegion* region = rnode->getACR();

	if (region->getSpace().size() == 0
			|| !region->getSpace().at(0)->hasOrientation()) {
		cout << "[mark]empty" << endl;
		return;
	}

	// get all ancestors and mark them with one Orientation from current node
	vector<ActiveConstraintGraph*> ancestors = cg->getAllAncestors();

	Orientation* orr = new Orientation(
			region->getSpace().at(0)->getOrientations().at(0));
	int ref_count = 0;

	for (int i = 0; i < ancestors.size(); i++) {
		pair<CgMarker::iterator, bool> res = cgs.insert(ancestors[i]);

		if (res.second == false) {	// if the ancestor is already in the marker
			delete ancestors[i];
			ancestors[i] = *(res.first);
		}

		if (cgchildren.find(ancestors[i]) == cgchildren.end()) { // initial the children array for newly added cg
			cgchildren[ancestors[i]] = vector<
					pair<ActiveConstraintGraph*, Orientation*> >();
		}
	}

	// The ancestors form a reversed tree, add the edges of this tree into marker
	for (int i = 0; i < ancestors.size(); i++) {
		for (int j = i + 1; j < ancestors.size(); j++) {// for each ancestors

			int parent = -1;
			int child = -1;

			// identify the parent and child
			if (ancestors[i]->getDim() == ancestors[j]->getDim() + 1
					&& ancestors[i]->IsParentOf(ancestors[j])) {
				parent = i;
				child = j;

			} else if (ancestors[j]->getDim() == ancestors[i]->getDim() + 1
					&& ancestors[j]->IsParentOf(ancestors[i])) {
				parent = j;
				child = i;
			}

			// add the edge between parent/child
			if (parent != -1 && child != -1) {
				vector<pair<ActiveConstraintGraph*, Orientation*> >& children =
						cgchildren[ancestors[parent]];
				int count = 0;
				for (vector<Edge>::iterator eit = children.begin();
						eit != children.end(); eit++) {
					if (eit->first == ancestors[child])
						count++;
					if (count >= 20)
						break;
				}
				// limit the number of child
				if (count < 20) {
					children.push_back(make_pair(ancestors[child], orr));
					ref_count++;
				}
			}
		}
	}

	// set the reference count of the created Orientation
	if (ref_count != 0) {
		this->ref_counter[orr] = ref_count;
	} else
		delete orr;
}

void CgMarker::try_release_orientation(Orientation* orr) {
	// decrease the reference counter of the Orientation. Free it if the counter is 0.
	if (this->ref_counter.find(orr) != this->ref_counter.end()) {
		this->ref_counter[orr]--;
		if (this->ref_counter[orr] <= 0) {
			this->ref_counter.erase(orr);
			delete orr;
		}
	}
}

pair<ActiveConstraintGraph*, vector<pair<ActiveConstraintGraph*, Orientation*> > > CgMarker::pop() {
	pair<ActiveConstraintGraph*,
			vector<pair<ActiveConstraintGraph*, Orientation*> > > res;

	// pop a arbitrary node
	res.first = *(cgs.begin());
	cgs.erase(res.first);

	// pop all its children
	res.second = cgchildren[res.first];
	cgchildren.erase(res.first);

	return res;
}

bool CgMarker::empty() {
	return cgs.empty();
}

vector<CgMarker::Edge> CgMarker::getEdge(ActiveConstraintGraph* cg) {
	vector<CgMarker::Edge> res;
	vector<pair<int, int> > parts = cg->getParticipants();

	// make sure order of the participant, the compare function need it
	sort(parts.begin(), parts.end());
	ActiveConstraintGraph* ncg = new ActiveConstraintGraph(parts,
			cg->getMolecularUnitA(), cg->getMolecularUnitB());

	if (cgs.find(ncg) != cgs.end()) {
		res = cgchildren[*(cgs.find(ncg))];
	} else {
		cout << "[mark]can't find cg" << endl;
	}
	delete ncg;
	return res;
}

void CgMarker::clear() {
	cgs.clear();
	cgchildren.clear();
	ref_counter.clear();
}
