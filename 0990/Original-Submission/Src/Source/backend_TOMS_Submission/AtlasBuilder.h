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
#ifndef ATLASBUILDER_H_
#define ATLASBUILDER_H_

#include "AtlasNode.h"
#include "SaveLoader.h"
#include "CgMarker.h"
#include "ConstraintCheck.h"
#include "ConvexChart.h"

#include <list>
#include <set>

/**
 * The goal of this class is to encapsulate the central algorithm of the MolecularUnit
 * configuration space program.
 * This class populates the ActiveConstraintRegion for each activeConstraintGraph by
 * sampling inside the boundaries of its ConvexChart. It creates and explores only regions
 * that contain at least one Cartesian realization, witness point.
 * If activeConstraintGraph is not partial 3-tree, then the region is populated by ray
 * tracing, i.e. witness points coming from all ancestors.
 */
class AtlasBuilder {
public:

	/**
	 * maximum number of atoms that can contribute to activeConstraintGraph for 2 rigid body packing case.
	 * For n=2, dof is 6. Hence there can be 6 contact pairs with each owning distinct atoms that results in 12 different atoms.
	 */
	static const int NO = 12;

	/**
	 * @brief This method is called to create the child nodes of rnode_prime and next generations.
	 * It adds the connections and adds the witness, but does NOT sample
	 * Descendant nodes will be created by adding contacts from contactList to cgKprime by one by one.
	 * and it will continue recursively
	 *----
	 * Creating descendant nodes will stop when descendant node dimension is 0 i.e. when cgKprime->getK() >= 6 even if contactList.size() is not emptied.
	 * The recursion could continue till overconstrained node, however this will increase the atlas size a lot !
	 *
	 * @param firstPath To show the first leaf node. Witness point will be saved only at the first leaf node. (to prevent memory blow-up)
	 */
	static void createChildNodes(AtlasNode *rnode_prime,
			std::list<std::pair<int, int> > contactList, Orientation* ori,
			Atlas* atlas, SaveLoader* snl, bool firstPath = true);

public:

	/**
	 * @brief This is the constructor method of the AtlasBuilder. Its job is the to some of the very basic setup steps.
	 * The AtlasBuilder holds onto pointers to influence and communicate with the other classes.
	 */
	AtlasBuilder(PointSet* a, PointSet* b, SaveLoader* snl,
			PredefinedInteractions *df, Atlas* atlas);

	/** @brief Destructor */
	virtual ~AtlasBuilder();

	/** @brief Set up the initial contact graphs for the root nodes */
	void setup();

	/**
	 * @brief Main entry of the AtlasBuilder
	 * Starts sampling root nodes one by one.
	 */
	virtual void startAtlasBuilding();

	/**
	 *  @brief This method is the recursive central point of the algorithm.
	 *
	 *  @param rnode				the roadmap node that is being sampled
	 *  @param dense 				is this sampling a refinement
	 *  @param coming_witness_ori 	the orientation at which location the sampling should start
	 *  @param continu 				whether this is a first/new sample or a continued sample.
	 *  @param bret 				is this a Breath First sampling
	 */
	bool sampleTheNode(AtlasNode *rnode, bool dense,
			Orientation* coming_witness_ori, bool continu, bool bret);

	/** @brief add witness point to all its parents if reverseWitness is enabled */
	void addReverseWitness(AtlasNode *rnode, ActiveConstraintGraph *cgK,
			ConvexChart *chart, ActiveConstraintRegion *region,
			Orientation* coming_witness_ori, bool & noGoodOrientation);

	virtual bool MyJacobianSampleRec(AtlasNode *rnode, bool dense,
			Orientation* orrw, bool continu, bool breath, bool bret) {
		return false;
	}
	;
	virtual bool MyJacobianSample(AtlasNode *rnode, bool dense,
			Orientation* orrw, bool continu, bool breath, bool bret) {
		return false;
	}
	;

	PointSet *a, *b;
	PredefinedInteractions *df;

protected:

	/** this value determines if some of the preset comments should be displayed. */
	const bool verbose;

	/** Set of CGs of the root nodes. bool is to show if CG is done or not */
	std::list<std::pair<ActiveConstraintGraph*, bool> > rootGraphs;

	/**
	 * An atlas object that will be poulated by AtlasBuilder.
	 * This object is shared between front-end and back-end of the algorithm.
	 */
	Atlas* atlas;
	SaveLoader* snl;

	CgMarker cgmarker;

protected:
	// helper functions

	/**
	 * @brief Computes Orientations for each flip of constraint graph for the given parametrization/description values.
	 *
	 * @see CartesianRealizer::computeRealization
	 * @return Set of feasible Orientations
	 */
	std::list<Orientation*> findRealizations(ActiveConstraintGraph *cgK,
			ConvexChart* des);

	/**
	 * @param[out] empty	If there are no unfinished rootGraph, empty will be set to true
	 *
	 * @return The next initial contact graph that has not been sampled and set it to be "done"
	 */
	ActiveConstraintGraph* getNextRootGraph(bool &empty);

	/*
	 * @brief Determine the step size proportional to the volume of the sample region
	 */
	void determineStepSizeDynamically(ActiveConstraintGraph *cgK,
			ActiveConstraintRegion * region, bool dense,
			CayleyParameterization * cparam);

	/*
	 * @brief This Method takes the 2 helix's an sets up the initial contact graphs that each holds
	 * two contacts to be used in root nodes
	 * Interaction in the distance table are used as the constraints of the graph
	 */
	void create_initial_contactGraphs_for_virusCase();

	/*
	 * @brief Creates the initial contact graphs of two contacts to be used in root nodes
	 * todo rename dumbell staff to something else
	 */
	void create_initial_4d_contactGraphs_usingDumbbells();

	/*
	 * @brief Creates the initial contact graphs each holds only one contact to be used in root nodes
	 */
	void create_initial_5d_contactGraphs();

	/**
	 * @brief For a given feasible grid point, check neighbor for collision and find the boundary by binary search.
	 *
	 * @see createChildContactGraps_fromTheBoundary
	 */
	void findBoundary(std::list<Orientation*>::iterator ori_on_lattice,
			ConvexChart* desc, ActiveConstraintGraph* cgK, AtlasNode* rnode,
			ActiveConstraintRegion *region, ConstraintCheck* detector,
			bool bret, bool& noGoodOrientation, int& noPoints,
			bool& boundary_ori_found_and_saved);

	/**
	 * @brief Add new contacts to the contact graph cgK and create new child contact graphs.
	 * The contacts are taken from contactList. Each time only 1 new contact is added to the parent rnode's contact graph cgK
	 *
	 * It also populate descendant nodes (without sampling) by createChildNodes method
	 *
	 * @param contactList All contacts that occur at the boundary orientation i.e. orie_on_boundary
	 * @param orie_on_boundary The orientation that is found by findBoundary method
	 *
	 * @see findBoundary
	 * @see createChildNodes
	 *
	 * previous name was addContactsAndCreateChildren
	 */
	void createChildContactGraps_fromTheBoundary(
			std::list<std::pair<int, int> >& contactList,
			std::list<Orientation*>::iterator ori_on_lattice,
			ActiveConstraintGraph* cgK, bool bret,
			Orientation* orie_on_boundary, AtlasNode* rnode,
			bool& noGoodOrientation, int& noPoints);

	/**
	 * @brief Create new child node corresponding to contact graph that is recently created by createChildContactGraps_fromTheBoundary method if they do not exist.
	 * Add it to atlas, and link it to parent node.
	 *
	 * Calls sampleTheNode that will do the sampling and search recursively.
	 *
	 * @see createChildContactGraps_fromTheBoundary
	 * @see sampleTheNode
	 */
	AtlasNode* createChildNode(AtlasNode *rnode,
			ActiveConstraintGraph *cgKprime, Orientation* wit_orr_toSend,
			bool & savedOnce, bool leafWitness, bool & noGoodOrientation,
			int & noPoints, bool bret);

	/**
	 * @brief refine the roadmap by checking missing node and connection in cgmarker
	 * and re-sample corresponding nodes.
	 *
	 * When a 0d node is found by sampling, all its ancestor nodes whould also be found.
	 * But some times some node are missing because of sampling error. So we mark the ancestors
	 * of all 0-d nodes that supposed to be found, and if they are missing, we re-sample its
	 * parent try to find it.
	 */
	void refineMap();

};

#endif /* ATLASBUILDER_H_ */
