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
 * Atlas.cpp
 *
 *  Created on: Feb 22, 2009
 *      Author: Admin
 */

#include "Atlas.h"

#include <time.h>
#include <string>
#include <fstream>

using namespace std;

vector<int> gnpaths;

Atlas::Atlas() {
	this->nodes.reserve(1000);
}

Atlas::Atlas(vector<AtlasNode*> nodes) {
	this->nodes = nodes;
	this->nodes.reserve(1000);

	for (size_t i = 0; i < this->nodes.size(); i++)
		if (this->nodes[i]->getDim()
				== Settings::RootNodeCreation::dimension_of_rootNodes)
			this->rootIndices.push_back(i);

}

Atlas::~Atlas() {
	cleanAtlas();
}

void Atlas::cleanAtlas() {
	for (vector<AtlasNode*>::iterator iter = this->nodes.begin();
			iter != this->nodes.end(); iter++) {

		(*iter)->trimNode();
		delete *iter;
	}

	this->nodes.clear();
	this->rootIndices.clear();
}

vector<AtlasNode*> Atlas::getNodes() {
	return this->nodes; //todo if it is time consuming to pass that vector (if it copies the memory addresses), then pass pointer to the vector
}

void Atlas::setNodes(std::vector<AtlasNode*> nds) {

	cleanAtlas();

	this->nodes = nds;

	for (size_t i = 0; i < this->nodes.size(); i++)
		if (this->nodes[i]->getDim()
				== Settings::RootNodeCreation::dimension_of_rootNodes)
			this->rootIndices.push_back(i);

}

size_t Atlas::number_of_nodes() {
	return this->nodes.size();
}

list<size_t> Atlas::getUnfinished() {
	list<size_t> output;
	for (size_t i = 0; i < this->nodes.size(); i++) {
		if (this->nodes[i]->isComplete() == 0) {
			output.push_back(i);
		}
	}
	return output;
}

list<size_t> Atlas::getChildren(size_t nodeNum, int level) {
	cout << "starting with node " << nodeNum << " looking for level " << level
			<< endl;
	list<size_t> output, temp;
	int curdim = this->nodes[nodeNum]->getDim();
	if (curdim != level)   //to prevent it to return itself
		output.push_back(nodeNum);
	for (int x = 0; x < (curdim - level) && !output.empty(); x++) {
		temp.insert(temp.end(), output.begin(), output.end());
		output.clear();

		for (list<size_t>::iterator iter = temp.begin(); iter != temp.end();
				iter++) {
			vector<int> con = this->nodes[*iter]->getConnection();
			int dim = this->nodes[*iter]->getDim();

			for (vector<int>::iterator it = con.begin(); it != con.end();
					it++) {
				if (this->nodes[*it]->getDim() < dim
						&& this->nodes[*it]->getDim() >= level) {
					output.push_back(*it);
				}
			}

		}
		output.sort();
		output.unique();

		temp.clear();
	}

	output.sort();
	output.unique();
	return output;

}

int Atlas::getSymmetricNode(AtlasNode * node) //works for 5d node
		{
	vector<pair<int, int> > particp = node->getCG()->getParticipants();
	if (particp[0].second != particp[0].first) {
		vector<pair<int, int> > parts;
		parts.push_back(make_pair(particp[0].second, particp[0].first));
		ActiveConstraintGraph* reverseNode = new ActiveConstraintGraph(parts,
				node->getCG()->getMolecularUnitA(),
				node->getCG()->getMolecularUnitB());
		int revnode = getNodeNum(reverseNode);
		delete reverseNode;
		return revnode;
	} else
		return -1;
}


int Atlas::getNodeNum(ActiveConstraintGraph* node) {

	//this code does not work if you created a node with random specific contacts(if it does not have any parents)
	for (size_t i = 0; i < this->rootIndices.size(); i++) {
		int rootIndex = this->rootIndices[i];
		ActiveConstraintGraph * root_cgk = this->nodes[rootIndex]->getCG();
		if (root_cgk->IsParentOf(node)) {
			if (node->getDim() == root_cgk->getDim()) {
				return rootIndex;
			}

			int found = findNodeAtTheChildsOfParent(rootIndex, node);
			if (found != -1) {
				return found;
			}
		}
	}
	return -1;
}

int Atlas::findNodeAtTheChildsOfParent(int parent,
		ActiveConstraintGraph* node) {

	if (node == NULL)
		return -1;
	vector<int> con = this->nodes[parent]->getConnection(); // con = the node numbers connected to this one
	for (vector<int>::iterator it = con.begin(); it != con.end(); it++) {
		ActiveConstraintGraph *cgKk = this->nodes[*it]->getCG(); //cgKk is an actual pointer to the graph of a connected node.
		if (cgKk->getDim() < this->nodes[parent]->getDim()
				&& cgKk->IsParentOf(node)) {
			if (cgKk->getDim() == node->getDim()) //it is good that you do not compare paramdim, because parameters of node is not set yet at the time this method is called :)
				return *it; //found

			int found = findNodeAtTheChildsOfParent(*it, node);
			if (found != -1)
				return found;
		}
	}
	return -1;
}

//returns node number now for efficiency since it takes time to find it and already computed here.
int Atlas::addNode(ActiveConstraintGraph* nodeA, int& output) {

	int added;
	output = this->getNodeNum(nodeA);  // check if it existed before.
	if (output < 0) {
		output = this->nodes.size();
		double loc[3];  // = new double[3];
		if (output == 0) {
			loc[0] = loc[1] = loc[2] = 0.0f;
		} else {
			double *oldLoc = (this->nodes.back())->getLocation();
			loc[2] = (10
					* (Settings::RootNodeCreation::dimension_of_rootNodes
							- nodeA->getDim()));
			loc[1] = oldLoc[1] + (rand() + 1) / ((double) RAND_MAX);
			loc[0] = oldLoc[0] + .1; // to prevent them from begin right on top of each other.
		}
		vector<int> blankConnect;
		AtlasNode *addNode = new AtlasNode(output, false, false,
				nodeA->getDim(), loc, blankConnect); // set empty to be false in order to display it initially
		addNode->setCG(nodeA);

		this->nodes.push_back(addNode);

		if (nodeA->getDim()
				== Settings::RootNodeCreation::dimension_of_rootNodes) // todo check if it has a parent in the atlas, if not put it to the rootindices as well?
			this->rootIndices.push_back(output);

		added = 1; //succesfully added
	} else {

		added = 0; //not added, existed before
	}

	return added;
}

void Atlas::connect(int indexA, int indexB) {
	if (indexA > -1 && indexB > -1 && indexA != indexB) {
		bool newConnection = !(this->isConnected(indexA, indexB));
		if (newConnection) {
			this->nodes[indexA]->addConnection(indexB);
			this->nodes[indexB]->addConnection(indexA);
		}
	}
}

bool Atlas::isConnected(int nodeA, int nodeB) {
	bool a = this->nodes[nodeA]->isConnectedTo(nodeB);
	bool b = this->nodes[nodeB]->isConnectedTo(nodeA);
	if (a != b)
		cerr << "\nasysmetric Roadmap connection\n";
	return a;
}

AtlasNode* Atlas::operator[](size_t id) {
	AtlasNode* output = this->nodes[id];
	return output;
}

void Atlas::BuildTree(vector<vertex*> &Graph) {
	vector<int> adj;
	for (vector<AtlasNode*>::iterator iter = this->nodes.begin();
			iter != this->nodes.end(); iter++) {
		if ((*iter)->getDim() > 1) {
			continue;
		}
		adj.clear();
		vector<int> push = (*iter)->getConnection();
		for (int i = 0; i < push.size(); i++) {
			if (this->nodes[push[i]]->getDim()
					<= Settings::Paths::energyLevelUpperBound) {
				adj.push_back(push[i]);
			}
		}
		vertex *v1 = new vertex((*iter)->getID(), adj);
		Graph.push_back(v1);

	}
}

//Helper method to find the index of a node in the newly built graph
int findindex(vector<vertex*> Graph, int num) {
	for (int i = 0; i < Graph.size(); i++)
		if (Graph[i]->number == num)
			return i;
	return -1;
}

void Atlas::BuildTree(vector<vector<int> > &matrix, vector<vertex*> &Graph) {
	for (int i = 0; i < Graph.size(); i++) {
		vector<int> adj = Graph[i]->adjList;
		for (vector<int>::iterator iter = adj.begin(); iter != adj.end();
				iter++) {
			int index = findindex(Graph, *iter);
			if (index != -1)
				matrix[i][findindex(Graph, *iter)] = 1;
		}
	}

}

int Atlas::findpath(int src, int dst, std::string relativepath) {
	std::queue<int> Q;
	vector<vertex*> Graph;
	Q.push(src);
	std::ofstream ofs;
	clock_t begin = clock();
	std::string pathfile = relativepath + "paths.txt";
	ofs.open(pathfile, std::ofstream::out | std::ofstream::app);

	BuildTree(Graph);

	//ofs<<"Finding the shortest path between "<<src<<" and "<<dst<<endl;
	if (findindex(Graph, src) == -1 || findindex(Graph, dst) == -1) {
		ofs << "The Source and Destination have to be 0D or 1D nodes\n" << endl;
		return 0;
	}

	while (!Q.empty()) {
		int current = Q.front();
		Q.pop();
		vertex *now = Graph[findindex(Graph, current)];
		Graph[findindex(Graph, current)]->visited = true;
		for (int i = 0; i < now->adjList.size(); i++) {
			if (Graph[findindex(Graph, now->adjList[i])]->visited != true) {
				Q.push(now->adjList[i]);
				Graph[findindex(Graph, now->adjList[i])]->visited = true;
				Graph[findindex(Graph, now->adjList[i])]->parent = current;
			}
		}
	}

	int cur = dst;
	bool print = true;
	while (cur != src) {
		if (Graph[findindex(Graph, cur)]->parent == -1) {
			clock_t end = clock();
			double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			return 0;
		}
		if (print == true) {
			ofs << "The path between " << src << " and " << dst << " is:";
			print = false;
		}
		ofs << cur << "--";

		cur = Graph[findindex(Graph, cur)]->parent;
	}
	clock_t end = clock();
	ofs << cur << endl;
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	ofs << "Found the path in " << elapsed_secs << " seconds" << endl;

	ofs.close();
	for (int i = 0; i < Graph.size(); i++)
		delete Graph[i];
	return 1;
}

void Atlas::findAllPaths(std::string relativepath) {
	vector<int> zeroDnodes;
	for (vector<AtlasNode*>::iterator iter = this->nodes.begin();
			iter != this->nodes.end(); iter++) {
		if ((*iter)->getDim() == 0) {
			zeroDnodes.push_back((*iter)->getID());
		}
	}

	for (int i = 0; i < zeroDnodes.size(); i++) {
		for (int j = i + 1; j < zeroDnodes.size(); j++) {
			findpath(zeroDnodes[i], zeroDnodes[j], relativepath);
		}
	}

}

void mul(vector<vector<int> > &matrix1, vector<vector<int> > &matrix2) {
	int num_nodes = matrix1.size();
	vector<vector<int> > square(num_nodes);
	for (int i = 0; i < num_nodes; i++)
		square[i].resize(num_nodes);

	for (int i = 0; i < num_nodes; i++) {
		for (int j = 0; j < num_nodes; j++) {
			for (int k = 0; k < num_nodes; k++) {
				square[i][j] += matrix1[i][k] * matrix2[k][j];
			}
		}
	}

	for (int i = 0; i < num_nodes; i++) {
		for (int j = 0; j < num_nodes; j++) {
			matrix1[i][j] = square[i][j];
		}
	}

}

void Atlas::findNumpaths(std::string relativepath) {
	std::ofstream ofs;
	clock_t begin = clock();
	clock_t end;
	std::string pathfile = relativepath + "path_matrix.tx";
	ofs.open(pathfile, std::ofstream::out | std::ofstream::app);

	if (Settings::Paths::pathLength <= 0) {
		ofs << "Path length must be positive" << endl;
		ofs.close();
		return;
	}

	vector<vertex*> Graph;
	BuildTree(Graph);

	int num_nodes = Graph.size();
	int n = Settings::Paths::pathLength;

	vector<vector<int> > matrix(num_nodes);
	for (int i = 0; i < num_nodes; i++)
		matrix[i].resize(num_nodes);

	vector<vector<int> > identity(num_nodes);
	for (int i = 0; i < num_nodes; i++)
		identity[i].resize(num_nodes);

	for (int i = 0; i < num_nodes; i++) {
		for (int j = 0; j < num_nodes; j++) {
			if (i == j) {
				identity[i][j] = 1;
			} else {
				identity[i][j] = 0;
			}
		}
	}

	BuildTree(matrix, Graph);

	if (Settings::Paths::pathLength == 1) {
		end = clock();
	} else {
		while (n > 0) {
			if (n % 2 == 0) {
				mul(matrix, matrix);
				n /= 2;
			} else {
				mul(identity, matrix);
				n--;
			}
		}
		end = clock();
	}

	ofs << "The path matrix is" << endl;
	for (int i = 0; i < num_nodes; i++) {
		for (int j = 0; j < num_nodes; j++) {
			ofs << matrix[i][j] << " ";
		}
		ofs << endl;
	}

	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	ofs << "Found number of paths in " << elapsed_secs << " seconds" << endl;
	ofs.close();
}
