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
#ifndef ACTIVECONSTRAINTGRAPH_H_
#define ACTIVECONSTRAINTGRAPH_H_

#include "PointSet.h"
#include "Settings.h"

#include <vector>
#include <set>
#include <string>

/*
 * A class used to store the set of active constraints and the corresponding vertices.
 */
class ActiveConstraintGraph {
public:

	/**
	 * maximum number of atoms that can contribute to activeConstraintGraph for 2 rigid body packing case.
	 * For n=2, dof is 6. Hence there can be 6 contact pairs with each owning distinct atoms that results in 12 different atoms.
	 */
	static const int NO = 12;

	/////////////////////////
	// Constructors
	////////////////////////

	/** @brief Default constructor */
	ActiveConstraintGraph();

	/** @brief Constructor with helA, helB, contacts and vertices of constraint graph initialization */
	ActiveConstraintGraph(std::vector<std::pair<int, int> > participants,
			PointSet* helA, PointSet* helB);

	/**
	 * @brief Constructor with helA, helB, contacts and vertices of constraint graph initialization
	 *
	 * This constructor is called after reading a node from file. It sets vertices of constraint graph from existing parameter PL information.
	 * or it can be called from AtlasBuilder to create a symmetric contact graph of pre-existing symmetric node
	 */
	ActiveConstraintGraph(std::vector<std::pair<int, int> > participants,
			PointSet* helA, PointSet* helB,
			std::vector<std::pair<int, int> > PL);

	/**
	 * @brief Constructor with helA, helB, contacts and vertices of constraint graph initialization
	 * This copy constructor is called by only from atlasBuilder class while creating child contactgraph just before a new contact is added.
	 */
	explicit ActiveConstraintGraph(ActiveConstraintGraph* ptrToCpy);

	/** @brief Destructor */
	virtual ~ActiveConstraintGraph();

	/////////////////////////
	// Getter/Setter
	////////////////////////

	/**
	 * @return The number of constraints i.e. number of contacts
	 */
	int constraintSize();

	/**
	 * @return Dimension of active constraint region represented by the contact graph
	 */
	size_t getDim();

	/**
	 * @brief The parameter dimension can be bigger than getDim if we allow sampling on the range for a contact.
	 * i.e. the contact distant is not fixed, it has a range.
	 *
	 * @return Dimension of the sampling space
	 */
	size_t getParamDim();

	/**
	 * @brief Getter of MolecularUnit
	 * @return A pointer to this graph's helixA
	 */
	PointSet* getMolecularUnitA();

	/**
	 * @brief Getter of MolecularUnit
	 * @return A pointer to this graph's helixB
	 */
	PointSet* getMolecularUnitB();

	/**
	 * @return Participants i.e. set of	contacts represented by <helix A "atom" index , helix B atom index>
	 */
	std::vector<std::pair<int, int> > getParticipants();

	/**
	 * @return Contacts i.e. set of	contacts represented by <helix A "vertex" index , helix B vertex index>
	 */
	std::vector<std::pair<int, int> > getContacts();

	/**
	 * @brief Getter for verticesA vector
	 * @return Vector of contact + parameter "atom" indices for the vertices of helixA side of contact graph
	 */
	vector<int> getVerticesA();

	/**
	 * @brief Getter for verticesB vector
	 * @return Vector of contact + parameter "atom" indices for the vertices of helixB side of contact graph
	 */
	vector<int> getVerticesB();

	/** @brief setter for verticesA vector */
	void setVerticesA(vector<int>);

	/** @brief setter for verticesB vector */
	void setVerticesB(vector<int>);

	/** @return The list of contact atom indices for the vertices of helixA side of contact graph */
	list<int> get_unordered_setA();

	/** @return The list of contact atom indices for the vertices of helixB side of contact graph */
	list<int> get_unordered_setB();

	/**
	 * getter for parameters
	 * @see MyDescription
	 */
	vector<pair<int, int> > getParameters();

	/**
	 * setter for parameters
	 * @see MyDescription
	 */
	void setParameters(vector<pair<int, int> > parameters);

	/**
	 * getter/setter of Parameter Lines.
	 * Parameter Line is pair of integer where each integer is the index of the atom in molecular unit.
	 */
	void setParamLines(vector<pair<int, int> > parameters);
	vector<pair<int, int> > getParamLines();

	/**
	 * getter/setter for sampling step size
	 */
	double getStepSize();
	void setStepSize(double siz = Settings::Sampling::stepSize);

	/**
	 * @return The number of different atoms in helix A that are part of contacts of constraint graph
	 */
	int getCNoinA();

	/**
	 * @return The number of different atoms in helix B that are part of contacts of constraint graph
	 */
	int getCNoinB();

	//////////////////////////////////
	// Relation between contact graphs
	//////////////////////////////////

	/**
	 * @return All contacts graphs with less contacts
	 */
	vector<ActiveConstraintGraph*> getAllAncestors();

	/**
	 * @return True if the contact set is a subset of the other's contact set, otherwise False
	 */
	bool IsParentOf(ActiveConstraintGraph* other);

	////////////////////////////
	// Operators
	////////////////////////////

	/** @brief This method proves a process of outputting the contact graph to a stream most likely the console */
	friend ostream& operator<<(ostream& os, ActiveConstraintGraph &cgid);

	//////////////////
	// Misc
	//////////////////

	/**
	 * @brief In order to realize the graph we need at least three atoms from each helix
	 * if it is more than 0 dim, then the parameter SHOULD handle the third atom
	 *
	 * if it is 0 dim, 6 contacts, the participants may NOT contain 3 different atoms in each side.
	 * And since it is 0 dim, there is no room for parameters that could fix this issue.
	 *
	 * @return True if constraint graph currently has (or may have in the future) 3 different vertices from both helix A and B side, False otherwise
	 */
	bool mayHaveThirdAtom();

	/**
	 * @brief 1 atom should be in contact with at most 3 atoms, otherwise fourth is dependent (over-constraint) and prevents realizing due to lack of edges.
	 *
	 * @return True if constraint graph is over constraint, False otherwise
	 */
	bool isDependent();

	/**
	 * @param x y are atoms' indices in helix
	 * @return True if (x,y) is a contact, False otherwise
	 */
	bool isConnected(int x, int y);

	/**
	 * @brief Add a new contact to constraint graph.
	 * Updates participants, contacts, unordered_setA, unordered_setB, verticesA, verticesB
	 *
	 * @param contact is a pair of integer where each integer is the "index of the atom in molecular unit".
	 */
	void addContact(pair<int, int> contact);

	/*
	 * @brief Get "index of the atom in molecular unit" of no'th vertex in verticesA
	 *
	 * these are like library functions that is necessary and do not exist in list class
	 * like [] operator
	 */
	int getA(int no);

	/** @brief Get "index of the atom in molecular unit" of no'th vertex in verticesB */
	int getB(int no);

	/** @return number of contact + parameter "atom"s for the vertices of helix A side of contact graph */
	int verticesA_size();

	/** @return number of contact + parameter "atom"s for the vertices of helix B side of contact graph */
	int verticesB_size();

	////////////////////
	// Friend function
	////////////////////

	/**
	 * @brief Compare two contact graph
	 */
	friend bool CGcomp(ActiveConstraintGraph* lhs, ActiveConstraintGraph* rhs);

	//    bool use_existing_parametrization; //symmetricExisted; //true if its symmetric node is created before. Then we will use same parametrization of preexisting node
	vector<int> independent_directions; //for jacobian sampling, independent rows on jacobian matrix

private:

	/////////////////////////
	// Helper functions
	////////////////////////

	/*
	 * @brief Add atom marker to make sure there are at least 3 atom markers from each molecular composite so that the graph is realizable.
	 * While choosing additional atoms, it has 2 options, choosing the closest atom to each other
	 * or the atoms that leads certain angle.
	 *
	 * @detail There has to exist at least 3 vertex in both side, if the contacts do not provide 3 different atoms,
	 * then we assign additional vertices according to some criteria : either by closest Atom or the one which makes perpendicular angle.
	 * @see completeTo3by3Graph_byClosestAtom, completeTo3by3Graph_that_makes_perpendicularAngle
	 *
	 * an example: for the contacts "06 07", there is only 1 vertex on the left size hence 2 more vertices are needed for left side
	 * and 1 more vertex needed for the right side.
	 *---
	 * additional vertices should be selected to keep symmetry, otherwise sampling density ends up not being well proportional among the nodes which have symmetric contacts !!!
	 * for instance, for atom a, if b is chosen, then for atom b, a should be chosen!
	 *---
	 * verticesA/B does not care about the order of atom indices
	 */
	void completeTo3by3Graph();

	/**
	 * @brief there has to exist at least 3 vertex in both side, this method assigns additional vertices according to closest Atom
	 *
	 * this method assumes the atoms in moleculerUnit are ordered according to location.
	 *
	 * @return vector of parameter "atom"s to be added to the vertices of contact graph
	 */
	vector<int> completeTo3by3Graph_byClosestAtom(PointSet * hel,
			vector<int> vertices, int size);

	/**
	 * @brief there has to exist at least 3 vertex in both side, this method assigns additional vertices according to the one which makes perpendicular angle.
	 *
	 * from one helix side, there are 3 atoms, means 2 lines. Those 2 lines should be perpendicular.
	 *
	 * @return vector of parameter "atom"s to be added to the vertices of contact graph
	 */
	vector<int> completeTo3by3Graph_that_makes_perpendicularAngle(
			PointSet * hel, vector<int> vertices, int size);

	/**
	 * @brief Convert paramLines into parameters format
	 *
	 * uses class member paramlines i.e. set of  <helix A atom index , helix B atom index> pairs
	 * updates class member parameters i.e. set of <verticesA index , verticesB index> pairs
	 *
	 * @see find_in_verticesA, find_in_verticesB
	 */
	void setParamsFromParamlines();

	/**
	 * @brief Convert participants into contacts format
	 * uses class member participants i.e. set of  <helix A atom index , helix B atom index> pairs
	 * updates class member contacts i.e. set of <verticesA index , verticesB index> pairs
	 *
	 * @see find_in_setA, find_in_setB
	 */
	void convertParticipantsToContacts();

	/*
	 * @brief Acts like find() operator
	 *
	 * @param atom Index of the atom in molecular unit
	 * @return vertex index of the participant 'atom' in the unordered_setA
	 * -1 means 'atom' is not part of contacts.
	 */
	int find_in_setA(int atom);

	/*
	 * @param atom Index of the atom in molecular unit
	 * @return vertex index of the participant 'atom' in the unordered_setB
	 * -1 means 'atom' is not part of contacts.
	 */
	int find_in_setB(int atom);

	/*
	 * @param atom Index of the atom in molecular unit
	 * @return vertex index of the 'atom' in the verticesA
	 * -1 means 'atom' is neither contact nor parameter
	 */
	int find_in_verticesA(int atom);

	/*
	 * @param atom Index of the atom in molecular unit
	 * @return vertex index of the 'atom' in the verticesB
	 * -1 means 'atom' is neither contact nor parameter
	 */
	int find_in_verticesB(int atom);

	/////////////////////////
	// Variables
	////////////////////////

	PointSet *helA, *helB;

	/*
	 * Set of Contacts
	 * A Contact is a pair of integer where each integer is the "index of the atom in molecular unit".
	 * i.e. A Contact is represented by <helix A "atom" index , helix B "atom" index>
	 */
	vector<pair<int, int> > participants;
	//todo maybe rename participant to contacts to prevent confusion and then distinguish labeling of index. maybe say contactLines ?

	/*
	 * Set of Contacts
	 * A Contact is a pair of integer where each integer is the "index of the vertex in constraint graph".
	 * i.e. A Contact is represented by <helix A "vertex" index , helix B "vertex" index>
	 *
	 * i th vertex from helix B is represented by i+NO/2 as index number
	 */
	vector<pair<int, int> > contacts;

	/*
	 * Set of atom marker index pairs that represent parameters
	 * A Parameter is a pair of integer where each integer is the "index of the atom in molecular unit".
	 * i.e. A Parameter is represented by <helix A "atom" index , helix B "atom" index>
	 */
	vector<pair<int, int> > paramlines;

	/*
	 * Set of atom marker index pairs that represent parameters
	 * A Parameter is a pair of integer where each integer is the "index of the vertex in constraint graph".
	 * i.e. A Parameter is represented by <helix A "vertex" index , helix B "vertex" index>
	 *
	 * i th vertex from helix B is represented by i+NO/2 as index number
	 */
	vector<pair<int, int> > parameters;

	/*
	 * list of contact atom indices for the vertices of contact graph
	 * contains unique unordered contact indices from helix. Only contacts, not the temporary parameter purpose vertices
	 * long time ago "set" was being used as a data structure, but since ordering of atom indices are not important and even causing complications with vertex indices when new contacts are added, I have changed it to unordered set format.
	 */
	list<int> unordered_setA, unordered_setB;

	/*
	 * vector of contact + parameter "atom" indices for the vertices of contact graph
	 * (atoms ,that are not in the actual contact list but hold for parameters, also exist in this list)
	 * the unique vertex list of the contact graph
	 * if unordered_setA size >= 3, then verticesA is exactly same as unordered_setA, otherwise verticesA has additional parameter vertices to complete 3 vertices.
	 */
	vector<int> verticesA, verticesB;

	double stepSize;

};

/*
 * define a total order of cg, in order to use in set
 */
bool CGcomp(ActiveConstraintGraph* lhs, ActiveConstraintGraph* rhs);

#endif /* ACTIVECONSTRAINTGRAPH_H_ */
