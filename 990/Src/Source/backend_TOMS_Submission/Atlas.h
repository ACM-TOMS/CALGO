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
#ifndef ATLAS_H_
#define ATLAS_H_

#define ROADHUDID 3

#include "AtlasNode.h"

#include <vector>
#include <queue>
#include <list>

class vertex {
public:
	int number;
	int npaths;
	vector<int> adjList;
	bool visited;
	int parent;
	vertex(int num, vector<int> con) {
		number = num;
		parent = -1;
		int n = con.size();
		for (int i = 0; i < n; i++) {
			adjList.push_back(con[i]);
		}
	}

};

/**
 * A representation of an assembly configuration space stratification into active constraint regions.
 * Atlas is a directed acyclic graph that represent the relation of active constraint regions.
 */
class Atlas {
public:

	/////////////////////////////////
	// Constructors/Destructors
	/////////////////////////////////

	/** @brief Default constructor */
	Atlas();

	/** @brief Constructor with nodes and rootIndices initialization */
	Atlas(std::vector<AtlasNode*> nodes);

	/** @brief Destructor */
	virtual ~Atlas();

	/////////////////////////////////
	// Getters
	/////////////////////////////////

	/**
	 * @return The nodes of the atlas
	 */
	std::vector<AtlasNode*> getNodes(); //todo if passing takes time maybe try something like const std::vector<AtlasNode*> &all

	/////////////////////////////////
	// Setters
	/////////////////////////////////

	/**
	 * @brief A setter for nodes and rootIndices
	 *
	 * @param nds The nodes of the atlas
	 */
	void setNodes(std::vector<AtlasNode*> nds);

	/////////////////////////////////
	// Other public methods
	/////////////////////////////////

	/**
	 * @param id The node's ID number
	 * @return The atlas node identified by id
	 */
	AtlasNode* operator[](size_t id);

	/**
	 * @brief A getter for the node's identifying number by depth first search.
	 *
	 * @see findNodeAtTheChildsOfParent
	 * @return The ID number of the node labeled by acg.
	 *
	 *todo rename to search? as in paper
	 */
	int getNodeNum(ActiveConstraintGraph* acg);

	/**
	 * @brief A getter for the node's id number by depth first search
	 * The search starts from the node identified by parent to direction of its children.
	 *
	 * @param parent The ID number of the node where the search starts from.
	 * @return The ID number of the node labeled by acg if found, -1 otherwise
	 */
	int findNodeAtTheChildsOfParent(int parent, ActiveConstraintGraph* acg);

	/**
	 * @return The number of nodes in the atlas
	 */
	size_t number_of_nodes();

	/**
	 * @return The list of id numbers of the incomplete nodes where the sampling is not finished yet.
	 */
	std::list<size_t> getUnfinished();

	/**
	 * @param level The dimension of the nodes of interest
	 * @param nodeNum The ID number of the ancestor node where the search starts from.
	 * @return The list of id numbers of the nodes with level DOF that are descendants of the node identified by nodenum.
	 */
	std::list<size_t> getChildren(size_t nodeNum, int level = -1);

	/**
	 * @brief This method is used only for 5d nodes
	 *
	 * @return The ID number of the node that has symmetric constraint graph to the node's constraint graph
	 * todo check if CayleyParametrization already creates symmetric parametrization or not.
	 * also consider having symmetry for all dimensions.
	 */
	int getSymmetricNode(AtlasNode * node);

	/**
	 * @return True if there is an edge/link between nodeA and nodeB, False otherwise.
	 */
	bool isConnected(int nodeA, int nodeB);

	/**
	 * @brief edits outNodeNumber and return success status
	 *
	 * @param acg The ActiveConstraintGraph of the node to be inserted to the atlas
	 * @param outNodeNumber The ID number of the node that is just added to the atlas
	 * @return 1 if successfully added, 0 if not added since existed before
	 */
	int addNode(ActiveConstraintGraph* acg, int& outNodeNumber);

	/**
	 * @brief Put a link between nodes identified by indexA and indexB.
	 */
	void connect(int indexA, int indexB);

	/**
	 * @brief Cleans the nodes and rootIndices. Cleans contact Graphs and Regions of the nodes as well.
	 */
	void cleanAtlas();

	/* 
	 * Method to find the shortest path between regions of the atlas.
	 */
	int findpath(int src, int dst, std::string relativepath);

	/* 
	 * Method to find all the shortest paths between 0D regions of the atlas.
	 */
	void findAllPaths(string relativepath);

	/* 
	 * Method to find the number of paths between regions of the atlas.
	 */
	void findNumpaths(string relativepath);

	/* 
	 * Helper method to build a tree that will be used for finding the shortest path between regions of the atlas.
	 */
	void BuildTree(vector<vertex*>& Graph);

	/* 
	 * Helper method to build a tree that will be used for finding the number of paths between regions of the atlas.
	 */
	void BuildTree(vector<vector<int> > &matrix, vector<vertex*> &Graph);

private:

	/*
	 * The collection of AtlasNodes.
	 */
	std::vector<AtlasNode*> nodes;

	/*
	 * The set of indices of highest dimensional nodes (that take place at the root of atlas)
	 */
	std::vector<size_t> rootIndices;
	// todo do not stick to dimension just keep root nodes. in case low dimensional node is added randomly as a root node.
	// todo create a tree structure (maybe linked list) for atlas.

};

#endif /* ATLAS_H_ */
