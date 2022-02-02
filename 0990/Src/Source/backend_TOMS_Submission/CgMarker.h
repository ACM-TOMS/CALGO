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
#ifndef CGMARKER_H_
#define CGMARKER_H_

#include "ActiveConstraintGraph.h"
#include "Orientation.h"
#include "AtlasNode.h"

#include <set>
#include <utility>
#include <vector>
#include <unordered_map>

/**
 * The class generate simple hash number for both ActiveConstraintGraph and Orientation
 */
class CgHasher {
public:
	size_t operator()(ActiveConstraintGraph* x) const {
		return size_t(x);
	}

	size_t operator()(Orientation* x) const {
		return size_t(x);
	}
};

/**
 * This class marks the ancestors of 0-d nodes.
 */
class CgMarker {
public:

	typedef std::pair<ActiveConstraintGraph*, Orientation*> Edge; // The child node and the boundary Orientation
	typedef std::set<ActiveConstraintGraph*,
			bool (*)(ActiveConstraintGraph*, ActiveConstraintGraph*)>::iterator iterator;

public:
	CgMarker();
	virtual ~CgMarker();

	/**
	 * Mark the ancestor of the given 0-d node
	 *
	 * @param node	the 0-d node
	 */
	void mark(AtlasNode* rnode);

	/**
	 * Check whether the CgMarker is empty.
	 *
	 * @return true is CgMarker is empty
	 */
	bool empty();

	/**
	 * Get a new node from the marker
	 *
	 * @return
	 */
	std::pair<ActiveConstraintGraph*, std::vector<Edge> > pop();

	/**
	 * clear the CgMarker
	 */
	void clear();

	/**
	 * get all the children of a node and the boundary Orientation
	 */
	std::vector<Edge> getEdge(ActiveConstraintGraph* cg);

	/**
	 * try to free the Orientation
	 *
	 * @param orr	the Orientation that need to be free
	 */
	void try_release_orientation(Orientation* orr);
private:

	// set of all contact graphs
	std::set<ActiveConstraintGraph*,
			bool (*)(ActiveConstraintGraph*, ActiveConstraintGraph*)> cgs;

	// the children of each contact graph(node)
	std::unordered_map<ActiveConstraintGraph*, std::vector<Edge>, CgHasher> cgchildren;

	// reference counter of each Orientation
	std::unordered_map<Orientation*, int, CgHasher> ref_counter;
};

#endif /* CGMARKER_H_ */
