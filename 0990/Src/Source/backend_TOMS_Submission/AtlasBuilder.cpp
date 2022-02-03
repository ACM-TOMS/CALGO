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
 * AtlasBuilder.cpp
 *
 *  Created on: 2008
 *      Author: Aysegul Ozkan
 */

#include "AtlasBuilder.h"

#include "PointSet.h"
#include "ActiveConstraintGraph.h"
#include "ActiveConstraintRegion.h"
#include "CayleyParameterization.h"
#include "CartesianRealizer.h"
#include "CayleyPoint.h"
#include "Atlas.h"
#include "Settings.h"
#include "Utils.h"

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/Dense"
#include "Eigen/SVD"
using namespace Eigen;
using Eigen::MatrixXd;

#include <cassert>
#include <vector>
#include <iterator>
#include <fstream>
using namespace std;

#define PI 3.14159265

int num_samples;

AtlasBuilder::AtlasBuilder(PointSet* helA, PointSet* helB, SaveLoader* snl,
		PredefinedInteractions *df, Atlas* atlas) :
		verbose(false) {
	this->snl = snl;
	this->atlas = atlas;

	if (this->verbose)
		cout << *helA << endl;
	if (this->verbose)
		cout << *helB << endl;
	if (this->verbose)
		cout << "AtlasBuilder" << endl;

	this->a = helA;
	this->b = helB;
	this->df = df;
}

AtlasBuilder::~AtlasBuilder() {
	// TODO Auto-generated destructor stub
}

void AtlasBuilder::create_initial_contactGraphs_for_virusCase() {
	if (this->verbose)
		cout << "setup" << endl;

	// for each pair of interaction in the distance table
	for (PredefinedInteractions::dist_iterator dit1 = this->df->dist2begin();
			dit1 != this->df->dist2end(); dit1++) {
		for (PredefinedInteractions::dist_iterator dit2 = dit1;
				dit2 != this->df->dist2end(); dit2++) {
			if (dit1 == dit2)
				continue;

			// get the 4 atoms
			Point* a1 = a->getAtomByLabel(dit1->first.first);
			Point* a2 = a->getAtomByLabel(dit2->first.first);

			Point* b1 = b->getAtomByLabel(dit1->first.second);
			Point* b2 = b->getAtomByLabel(dit2->first.second);

			if (this->verbose) {
				cout << "checking " << a1->getName() << " " << a2->getName()
						<< endl;
				cout << "\t" << b1->getName() << " " << b2->getName() << endl;
			}

			// get the distance of the atom pairs
			double da = Utils::dist(a1->getLocation(), a2->getLocation());
			double db = Utils::dist(b1->getLocation(), b2->getLocation());

			if (this->verbose)
				cout << "da db: " << da << " , " << db << endl;

			// if the distance is outside the acceptable range, ignore it
			if (da < Settings::RootNodeCreation::initial4DContactSeparation_low
					|| da
							> Settings::RootNodeCreation::initial4DContactSeparation_high
					|| db
							< Settings::RootNodeCreation::initial4DContactSeparation_low
					|| db
							> Settings::RootNodeCreation::initial4DContactSeparation_high)
				continue;

			double ar1 = a1->getRadius();
			double ar2 = a2->getRadius();

			double br1 = b1->getRadius();
			double br2 = b2->getRadius();

// create the contact graph of pair of atom and add it into the rootGraph list
			{
				vector<pair<int, int> > parts;
				parts.push_back(
						make_pair(a->getIndexOf(a1), b->getIndexOf(b1)));
				parts.push_back(
						make_pair(a->getIndexOf(a2), b->getIndexOf(b2)));
				ActiveConstraintGraph* initial = new ActiveConstraintGraph(
						parts, this->a, this->b);
				this->rootGraphs.push_back(make_pair(initial, true));

			}

		}
	}

}

void AtlasBuilder::create_initial_4d_contactGraphs_usingDumbbells() {
	if (this->verbose)
		cout << "create_initial_4d_contactGraphs_usingDumbbells" << endl;

	vector<pair<int, int> > vectorA, vectorB; // these are used to store the dumbbells generated within each helix.
	vectorA = a->getDumbbells(); // getting the dumbbell candidates
	vectorB = b->getDumbbells();

	if (this->verbose)
		cout << "got dumbbells setA[" << vectorA.size() << "] setB["
				<< vectorB.size() << "]" << endl;
	vector<pair<int, int> >::iterator iter, it;

	//setup to be bi-incident?
	// ar1 ---da---ar2 (helix a)
	//  |           |
	//  |           |
	// br1 ---db---br2 (helix b)
	//Each of the variables get set and then the comparison is made.

	for (iter = vectorA.begin(); iter != vectorA.end(); iter++)
		for (it = vectorB.begin(); it != vectorB.end(); it++) {

			double z_slide_first = abs(
					(a->getAtomAt((*iter).first))->getLocation()[2]
							- (b->getAtomAt((*it).first))->getLocation()[2]); //i can do this because they are originally aligned to same position in the z axis
			double z_slide_second = abs(
					(a->getAtomAt((*iter).second))->getLocation()[2]
							- (b->getAtomAt((*it).second))->getLocation()[2]); //i can do this because they are originally aligned to same position in the z axis

			double da = Utils::dist(
					(a->getAtomAt((*iter).first))->getLocation(),
					(a->getAtomAt((*iter).second))->getLocation());
			double ar1 = (a->getAtomAt((*iter).first))->getRadius();
			double ar2 = (a->getAtomAt((*iter).second))->getRadius();
			double ar12 = ar1 + ar2;
			//	double ar12 = a->getAtomAt( (*iter).first )->getMinDist( a->getAtomAt( (*iter).second) );  //if( ar12 == -1 )  ar12 = ar1 + ar2;

			double db = Utils::dist((b->getAtomAt((*it).first))->getLocation(),
					(b->getAtomAt((*it).second))->getLocation());
			double br1 = (b->getAtomAt((*it).first))->getRadius();
			double br2 = (b->getAtomAt((*it).second))->getRadius();
			double br12 = br1 + br2;
			//	double br12 = b->getAtomAt( (*it).first )->getMinDist( b->getAtomAt( (*it).second) );

			//are they able to be bi-incident?
			if ((da + ar12 >= db - br12) && (da - ar12 <= db + br12)) {

				if (!Settings::RootNodeCreation::useParticipatingPointZDistance
						|| (z_slide_first
								< Settings::RootNodeCreation::ParticipatingPointZDistance
								&& z_slide_second
										< Settings::RootNodeCreation::ParticipatingPointZDistance)) //( abs((*iter).first - (*it).first)<4  &&   abs((*iter).second - (*it).second)<4  )  //to choose dumbbells proportional to the place in the helix (the other direction is ignored now)
						{
					//first paring
					vector<pair<int, int> > parts; //holds the indices of atoms
					parts.push_back(make_pair((*iter).first, (*it).first));
					parts.push_back(make_pair((*iter).second, (*it).second));
					ActiveConstraintGraph* initial = new ActiveConstraintGraph(
							parts, this->a, this->b);
					this->rootGraphs.push_back(make_pair(initial, true)); //contact_graphs with 2 contacts

				}

				if (Settings::RootNodeCreation::reversePairDumbbells) //reverse pairing
				{
					double z_slide_first_reverse =
							abs(
									(a->getAtomAt((*iter).first))->getLocation()[2]
											- (b->getAtomAt((*it).second))->getLocation()[2]); //i can do this because they are originally aligned to same position in the z axis
					double z_slide_second_reverse =
							abs(
									(a->getAtomAt((*iter).second))->getLocation()[2]
											- (b->getAtomAt((*it).first))->getLocation()[2]); //i can do this because they are originally aligned to same position in the z axis

					if (!Settings::RootNodeCreation::useParticipatingPointZDistance
							|| (z_slide_first_reverse
									< Settings::RootNodeCreation::ParticipatingPointZDistance
									&& z_slide_second_reverse
											< Settings::RootNodeCreation::ParticipatingPointZDistance)) //Settings::RootNodeCreation::closeByDumbbellsAmount=5
							{
						vector<pair<int, int> > parts;
						parts.push_back(make_pair((*iter).first, (*it).second));
						parts.push_back(make_pair((*iter).second, (*it).first));
						ActiveConstraintGraph* initial =
								new ActiveConstraintGraph(parts, this->a,
										this->b);
						this->rootGraphs.push_back(make_pair(initial, true));

					}
				}
			}
		}

}

void AtlasBuilder::create_initial_5d_contactGraphs() {
	if (this->verbose)
		cout << "create_initial_5d_contactGraphs" << endl;

	vector<Point*> helA = a->getAtoms();
	vector<Point*> helB = b->getAtoms();

	if(Settings::Sampling::sampleAllNodes) {
	for (size_t i = 0; i < helA.size(); i++) {
		for (size_t j = 0; j < helB.size(); j++) {

			double z_slide = abs(
					helA[i]->getLocation()[2] - helB[j]->getLocation()[2]); //i can do this because they are originally aligned to same position in the z axis

			if (!Settings::RootNodeCreation::useParticipatingPointZDistance
					|| z_slide
							< Settings::RootNodeCreation::ParticipatingPointZDistance) //( abs(i- j)<4 )  //to choose dumbbells proportional to the place in the helix (the other direction is ignored now)
							{
				vector<pair<int, int> > parts; //holds the indices of atoms
				parts.push_back(make_pair(i, j));
				ActiveConstraintGraph* initial = new ActiveConstraintGraph(
						parts, a, b);
				this->rootGraphs.push_back(make_pair(initial, true)); //contact_graphs with 1 contacts

			}
		}
	}
	} else {
				vector<pair<int, int> > parts; //holds the indices of atoms
				parts.push_back(make_pair(11, 16));
				ActiveConstraintGraph* initial = new ActiveConstraintGraph(parts, a, b);
				this->rootGraphs.push_back(make_pair(initial, true)); //contact_graphs with 1 contacts

	}

}

void AtlasBuilder::setup() {
	if (this->verbose)
		cout << "setup" << endl;

	if (Settings::General::candidate_interactions)
		create_initial_contactGraphs_for_virusCase();
	else if (Settings::RootNodeCreation::dimension_of_rootNodes == 5)
		create_initial_5d_contactGraphs();
	else if (Settings::RootNodeCreation::dimension_of_rootNodes == 4)
		create_initial_4d_contactGraphs_usingDumbbells();

	cout << "this->rootGraphs.size() " << this->rootGraphs.size() << endl;

//--------------------------------

	if (this->verbose)
		cout << "Created " << this->rootGraphs.size() << " contactIDs" << endl;
	for (list<pair<ActiveConstraintGraph*, bool> >::iterator iter =
			this->rootGraphs.begin(); iter != this->rootGraphs.end(); iter++) {
		if (this->verbose) {
			cout << *((*iter).first) << endl;
		}
	}

//--------------------------------

	// reorder rootGraphs to start from the middle through end and then beginning
	// this is done because the middle rootGraphs are more likely be of interest and we want to see them first.

	int rootGraphSize = this->rootGraphs.size();
	list<pair<ActiveConstraintGraph*, bool> >::iterator iterl =
			this->rootGraphs.begin();
	for (int i = 0; i < rootGraphSize / 3; i++) //find the 1/3 rootGraph
		iterl++;

	for (int i = rootGraphSize / 3; i > 0; i--) //from 1/3 through the beginning, push rootGraphs at the end
			{
		this->rootGraphs.push_back(*iterl);
		iterl--;
	}
	for (int i = rootGraphSize / 3; i > 0; i--) //remove first 1/3 rootGraphs
		this->rootGraphs.pop_front();

}

ActiveConstraintGraph* AtlasBuilder::getNextRootGraph(bool &empty) {
	empty = true;
	ActiveConstraintGraph* nextd;
	for (list<pair<ActiveConstraintGraph*, bool> >::iterator it =
			this->rootGraphs.begin(); it != this->rootGraphs.end(); it++) {
		if ((*it).second) //not done
		{
			nextd = (*it).first;
			(*it).second = false; //set it done
			empty = false;  //able to find one more rootGraph
			break;
		}
	}
	return nextd;
}

void AtlasBuilder::startAtlasBuilding() {

	clock_t begin = clock();
	num_samples = 0;
	std::ofstream ofs;
	std::string samplesfile = snl->relativePath + "Samples.txt";
	ofs.open(samplesfile, std::ofstream::out | std::ofstream::app);

	cout << "AtlasBuilder::startAtlasBuilding: this->rootGraphs.size() = "
			<< this->rootGraphs.size() << endl;
	int currentrootGraph = 0;

	bool empty;
	ActiveConstraintGraph* nextrootGraph = getNextRootGraph(empty);

	while (!empty) {
		if (!Settings::AtlasBuilding::stop) {
			int nodenum;
			int success = this->atlas->addNode(nextrootGraph, nodenum);
			if (success == 1) {
				cout
						<< "AtlasBuilder::startAtlasBuilding: We start sampling rootGraph "
						<< currentrootGraph << " out of "
						<< this->rootGraphs.size() << " (node number "
						<< nodenum << ")" << endl;
				AtlasNode * rnode = (*this->atlas)[nodenum];

				int symmNodeNum = this->atlas->getSymmetricNode(rnode);
				if (symmNodeNum != -1) // if there is a node with symmetric contact graph, then we want to exploit its parameters
						{
					AtlasNode * symmNode = (*this->atlas)[symmNodeNum];

					vector<pair<int, int> > particp =
							symmNode->getCG()->getParticipants();
					vector<pair<int, int> > sparticp;
					sparticp.push_back(
							make_pair(particp[0].second, particp[0].first));

					vector<pair<int, int> > plines =
							symmNode->getCG()->getParamLines();
					vector<pair<int, int> > splines;
					for (int i = 0; i < plines.size(); i++)
						splines.push_back(
								make_pair(plines[i].second, plines[i].first));

					ActiveConstraintGraph* reverseNode =
							new ActiveConstraintGraph(sparticp, this->a,
									this->b, splines);
					rnode->setCG(reverseNode);

				}

				// start the sampling
				if (Settings::AtlasBuilding::stop == true) {
					return;
				} else {
					this->sampleTheNode(rnode, false, NULL, false, false);
				}

				if (Settings::AtlasBuilding::stop == true) {
					return;
				}

			}

			ofs << "The total number of samples for node " << nodenum << " is "
					<< num_samples << endl;
			num_samples = 0;
			nextrootGraph = getNextRootGraph(empty);

			currentrootGraph++;
		}
	}

	// try refine the Roadmap Node that is missing from the original sampling
	if (Settings::General::reverseWitness)
		refineMap();

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout
			<< "FINISHED AtlasBuilder::startAtlasBuilding: this->rootGraphs.size() = "
			<< this->rootGraphs.size() << endl;
	cout << "Finding number of paths of length " << Settings::Paths::pathLength
			<< endl;
	ofs << "The sampling was completed in " << elapsed_secs << " seconds"
			<< endl;
	ofs.close();

	std::string datapath = snl->relativePath;
	cout << "Started finding shortest paths" << endl;
	atlas->findAllPaths(snl->relativePath);
	cout << "Finished finding shortest paths" << endl;

	cout << "Started finding the number of paths" << endl;
	atlas->findNumpaths(snl->relativePath);
	cout << "Finished finding the number of paths" << endl;

}

void AtlasBuilder::determineStepSizeDynamically(ActiveConstraintGraph *cgK,
		ActiveConstraintRegion * region, bool dense,
		CayleyParameterization * cparam) {
	if (this->verbose)
		cout << "determineStepSizeDynamically" << endl;

	double testStep = 0.4; //was 0.2 before

	// temporarily disable the steric constraint to allow collision on parameters
	// bool userDefined_stericConstraint = Settings::Constraint::stericConstraint;
	// Settings::Constraint::stericConstraint = false;

	cgK->setStepSize(testStep); //increment this stepsize, if it is causing slowness
	ConvexChart *chart = new ConvexChart(cgK, dense, cparam, this->df); //for volume computation, to make blue points as minimum as possible

	/**
	 * GRID sampling allows some range for contact hence GRID sampling volume per node is proportional with that range of the contact.
	 * i.e. Not all root nodes have same volume by GRID sampling.
	 * If we are not doing short-range sampling and forcing contact to be exact distance,
	 * then in order to have proportional sampling density with GRID sampling,
	 * we should set expectedNumberOfSamples per node to be proportional with that ratio
	 */
	double contact_lengthUpper = df->bondingUpperBound(chart->getAtom(0),
			chart->getAtom(6));
	double contact_lengthLower = df->bondingLowerBound(chart->getAtom(0),
			chart->getAtom(6));
	double range = pow(contact_lengthUpper, 3) - pow(contact_lengthLower, 3);

	double expectedNumberOfSamples = 500.; //10000.; //50000.;
	expectedNumberOfSamples = expectedNumberOfSamples * range;

	/** compute the approximate volume of the region */
	int volume = 0;
	if (chart->initializeChart(false, region) && cparam->is_partial3tree())
		for (; !chart->doneGrid() && !Settings::AtlasBuilding::stop;
				chart->stepGrid())
			volume++;

	delete chart;

	double voltimes = expectedNumberOfSamples / volume;
	double root = 1. / cgK->getDim();
	double pwr = pow(voltimes, root);
	double stepsizeee = testStep / pwr;
	cgK->setStepSize(stepsizeee);

}

void AtlasBuilder::createChildContactGraps_fromTheBoundary(
		list<pair<int, int> >& contactList,
		list<Orientation*>::iterator ori_on_lattice, ActiveConstraintGraph* cgK,
		bool bret, Orientation* orie_on_boundary, AtlasNode* rnode,
		bool& noGoodOrientation, int& noPoints) {

	ActiveConstraintGraph* cgKprime[contactList.size()];
	bool firstPath = true;
	bool savedOnce = false;
	for (int c = 0; c < contactList.size(); c++)  //TRY ALL CONTACT LIST
			{
		cgKprime[c] = new ActiveConstraintGraph(cgK);
		pair<int, int> contact = contactList.front();
		cgKprime[c]->addContact(contact);
		contactList.pop_front();

		if (cgKprime[c]->constraintSize() > 6 || cgKprime[c]->isDependent()
				|| !cgKprime[c]->mayHaveThirdAtom()) {
			delete cgKprime[c];
			cgKprime[c] = NULL;
			contactList.push_back(contact);
			continue;
		}

		//--------------
		Orientation* wit_orr_toSend = new Orientation(orie_on_boundary); // orie_on_boundary->getOrienation() ;
		bool leafWitness = (contactList.size() == 0);
		AtlasNode* child_node = createChildNode(rnode, cgKprime[c],
				wit_orr_toSend, savedOnce, leafWitness, noGoodOrientation,
				noPoints, bret);
		int child_nodeID = child_node->getID();
		orie_on_boundary->addBoundary(child_nodeID);
		(*ori_on_lattice)->addBoundary(child_nodeID);
		//--------------

		createChildNodes(child_node, contactList, wit_orr_toSend, this->atlas,
				this->snl, firstPath);
		firstPath = false; // wit_orr_toSend will be saved through only firstPath

		contactList.push_back(contact); //to let 'for loop' to stop
		delete wit_orr_toSend;
	}

}

void AtlasBuilder::findBoundary(list<Orientation*>::iterator ori_on_lattice,
		ConvexChart* desc, ActiveConstraintGraph* cgK, AtlasNode* rnode,
		ActiveConstraintRegion *region, ConstraintCheck* detector, bool bret,
		bool& noGoodOrientation, int& noPoints,
		bool& boundary_ori_found_and_saved) {
	// save the current grid point so we can come back where we were after binary search
	vector<double> pp = desc->getPoint();

	// get the flip on which we want to run binary search
	int flip = (*ori_on_lattice)->getFlipNum();

	// look up and down for one step in each direction around a valid point.
	// If there is a point with collision, then that means there is a contact
	// boundary in between valid point and colliding point.
	// Then you need to do binary search in between to find out that boundary config.
	while (!desc->stepAround() && !Settings::AtlasBuilding::stop) {
		bool fail;
		num_samples++;
		Orientation* orie_on_boundary = CartesianRealizer::computeRealization(
				cgK, desc, flip, fail); //if needs maple, then it will always return first root not the one with specified flip !

		/// found collision, start binary search
		/// todo ACTUALLY CHECK SHOULD BE DONE AROUND BLUE REGION AS WELL,
		/// IF VOLUME NEGATIVE, SEARCH FOR BOUNDARY IN BETWEEN!!!
		if (Settings::Sampling::binarySearch) {
			if (!fail && detector->stericsViolated(orie_on_boundary)) {
				desc->stepGridBinary(false);  ////walk through valid point

				delete orie_on_boundary;
				num_samples++;
				orie_on_boundary = CartesianRealizer::computeRealization(cgK,
						desc, flip, fail);

				// check for collision and contacts
				// if there is a collision then it should continue binary stepping
				bool binvalid = !detector->stericsViolated(orie_on_boundary);
				list<pair<int, int> > contactList2 = detector->getContacts(
						flip);

				// binary search
				// continue search till you find a valid configuration with new contacts
				// find new contact for new small threshold t // FIXME: what does this line mean?
				int num_bin_iteration = 0;
				while ((contactList2.empty() || !binvalid)
						&& num_bin_iteration < 30) {
					num_bin_iteration++;
					desc->stepGridBinary(binvalid);
					delete orie_on_boundary;
					num_samples++;
					orie_on_boundary = CartesianRealizer::computeRealization(
							cgK, desc, flip, fail);
					binvalid = !detector->stericsViolated(orie_on_boundary);
					contactList2 = detector->getContacts(flip);
				}

				// double check if the valid configuration with new contacts exists.
				// it may not in case it exits the loop because of num_bin_iteration is big
				if (binvalid && contactList2.size() != 0
						&& (Settings::AtlasBuilding::ifBadAngleWitness_createChild
								|| !orie_on_boundary->angleViolated()))      ///
						{
					//									Orientation* orie_on_boundary = sr->getOrienation();
					createChildContactGraps_fromTheBoundary(contactList2,
							ori_on_lattice, cgK, bret, orie_on_boundary, rnode,
							noGoodOrientation, noPoints);
					if (Settings::Saving::saveBoundary) //(since we do not save it to child)   //do not add witness point to parent node!
					{
						vector<double> outt = desc->getPoint(); //to get exact boundary position instead of grid position
						CayleyPoint* p4 = new CayleyPoint(outt); //use p4 INSTEAD OF p4d  to save exact position of point after binary search
						if (!orie_on_boundary->angleViolated()) {
							p4->addOrientation(orie_on_boundary);
							boundary_ori_found_and_saved = true;
						} else {
							delete orie_on_boundary; //orange
//						p4->badAngle = true;
							p4->incrementBadAngleN();
						}
						region->insertSpace(p4); //should it be insertwitness !!! no, witness is the one saved at the child node. we are here saving the boundary config to current node.

						noPoints++;
					} else
						delete orie_on_boundary;

					break; //if found a collision in one direction then stop
				} // if( valid(orie_on_boundary,detector) && contactList2.size() != 0
			} // if( !fail &&  !valid(orie_on_boundary,detector)  )else
		}
		delete orie_on_boundary;

		desc->setInitialPoint(pp);
	} // while( !desc->stepGridContact()  && !Settings::AtlasBuilding::stop)
	desc->setInitialPoint(pp);
	desc->setDir();
}

/*
 * This is an overview of the AtlasBuilder, though the actual AtlasBuilder may have changed a bit.
 *
 * Sample(CG(b1,b2,b3...bm) , k) :: m =12 & k=6
 desc = Get description ( CG(b1,b2,b3...bm) , k )
 :: find good parameters (depends on CG)
 & inequalities (depends on identity of  b1,b2,b3...bm)
 FOR each grid point in (convex)desc
 Real=find cartesian realizations w/out helix(CG(b1,b2,b3...bm),k)
 FOR each r in Real
 IF valid w/helix(r)
 getCG(r) = (CG(b1,b2,b3...bm'),k')
 output r with CG
 IF(CG(b1,b2,b3...bm'),k') != ( CG(b1,b2,b3...bm) , k )
 IF k < 6 and notdone(CG(b1,b2,b3...bm'),k')::CG & k used to identify
 Sample(CG(b1,b2,b3...bm'),k')
 FI
 done(CG(b1,b2,b3...bm'),k')
 FI
 FI
 ROF
 ROF
 done(CG(b1,b2,b3...bm),k)
 */

bool AtlasBuilder::sampleTheNode(AtlasNode *rnode, bool dense,
		Orientation* coming_witness_ori, bool continu, bool bret) {
	int from = rnode->getID(); // from is the node number of the incoming node as stored in the Roadmap object.
	if (this->verbose)
		cout << "AtlasBuilder::MySample: Started node " << from << endl;

	if (dense) //it is now being refined which means it isn't a completed graph anymore until refinement is finished
		rnode->setComplete(false); //later this information will be saved to RoadMap.txt file

	ActiveConstraintRegion *region = rnode->getACR();
	ActiveConstraintGraph *cgK = rnode->getCG();
	ConstraintCheck *detector = new ConstraintCheck(cgK, this->df); //Create a fresh ConstraintCheck detector
	CayleyParameterization *desc = new CayleyParameterization(cgK, false); //Create a fresh description

	if (Settings::Sampling::dynamicStepSizeAmong)
		determineStepSizeDynamically(cgK, region, dense, desc);

	// save node after parameters is set (by creating description object).
	//
	// recent info for the completeness and anygoodorientation will be saved to roadmap.txt
	// later by saveMapView() which is called by some time periods.
	//
	// you could saveNode at the end of this method, but it is better to saveNode before going
	// recursive call not to lose at least any available important information in any crush.
	if (!continu)
		snl->saveNodeHeader(rnode, this->atlas);

	if (!continu)
		this->snl->appendDimension(rnode); //dimension information in the first line of NodeX.txt

	if (continu && !Settings::AtlasBuilding::breadthFirst && !dense) {
		vector<int> con = rnode->getConnection(); // con = the node numbers connected to this one
		size_t dim = rnode->getDim();
		for (vector<int>::iterator it = con.begin(); it != con.end(); it++) {
			AtlasNode *child_node = (*this->atlas)[*it];
			if (child_node->getDim() < dim && !child_node->isComplete()) // if it is a lesser dimension and incomplete
					{
				snl->loadNode(*it, child_node->getACR());
				this->sampleTheNode(child_node, false, NULL, continu, false);
			}
		}
	}

	bool noGoodOrientation = true;
	if (continu)
		noGoodOrientation = !rnode->hasAnyGoodOrientation();

	ConvexChart *chart = new ConvexChart(cgK, dense, desc, this->df);

	// if the coming_witness_ori is not NULL, calculate the parameter value
	if (coming_witness_ori != NULL) {
		CayleyPoint* wp4 = new CayleyPoint(coming_witness_ori,
				cgK->getParamLines(), a, b); //deletion of coming_witness_ori is handled WHERE IT IS CREATED.

		chart->setWitnessPoint(wp4);
		if (!coming_witness_ori->angleViolated()) {
			region->insertWitness(wp4);
			// cgK->insertWitness(wp4); //  add it to the graph // it may be needed as a starting point to sample in the description class
			this->snl->appendWitness(rnode, wp4); // save it to the nodeX.txt  //note: it will 'not' be saved to file 'twice' in next steps again. only space points are saved at the next steps, there is no more witness and witness saving job later.
			noGoodOrientation = false;
		}
	}

	addReverseWitness(rnode, cgK, chart, region, coming_witness_ori,
			noGoodOrientation);

	if (bret) //bret true means this is child node of 'the node that is sampled by breathFirst'. Hence we need to return back.
	{
		if (this->verbose)
			cout << "Returning\n" << endl;
		region->trim();
		delete detector;
		delete desc;
		delete chart;
		return true;
	}

	bret = Settings::AtlasBuilding::breadthFirst;
	int noPoints = 0;
	int loop = 0;

	//continu : Flag of whether this is a first/new sample or a continued sample.
	if (desc->is_partial3tree() && chart->initializeChart(continu, region)) {

		list<Orientation*> real;
		for (; !chart->doneGrid() && !Settings::AtlasBuilding::stop;
				chart->stepGrid()) {
			// find all realization for the grid point
			real = findRealizations(cgK, chart);
			vector<double> out = chart->getPoint();
			CayleyPoint* p4d = new CayleyPoint(out);
			p4d->setRealizable(!real.empty());

			if (!real.empty())
				for (list<Orientation*>::iterator ori_on_lattice = real.begin();
						ori_on_lattice != real.end()
								&& !Settings::AtlasBuilding::stop;
						ori_on_lattice++) {
					// check for steric constraint violation and for any possible future contacts
					if (!detector->stericsViolated(*ori_on_lattice)) {
						/// check for bad angle
						if ((*ori_on_lattice)->angleViolated())
							p4d->incrementBadAngleN();
						else
							noGoodOrientation = false;

						bool boundary_ori_found_and_saved = false;

						//list<pair<int,int> > contactList = detector->getContacts();
						//if( !contactList.empty() ) //it is easy to miss contacts through grid stepping by just 'contactList.empty()' check
						if (Settings::RootNodeCreation::createChildren) //set it false JUST TO CREATE INTERIOR of root nodes
						{
							// find the boundary by binary search
							findBoundary(ori_on_lattice, chart, cgK, rnode,
									region, detector, bret, noGoodOrientation,
									noPoints, boundary_ori_found_and_saved);
						}

						// save both 'ori_on_lattice' and 'orie_on_boundary'
						if (!boundary_ori_found_and_saved) {
							if (!(*ori_on_lattice)->angleViolated()) {
								p4d->addOrientation((*ori_on_lattice));
							} else {
								delete (*ori_on_lattice); //orange
//							p4d->badAngle = true;  //commented since incrementBadAngleN done above
							}
						} else
							delete (*ori_on_lattice);

					} else {
						if ((*ori_on_lattice)->angleViolated())
							p4d->incrementBadAngleN();

						delete (*ori_on_lattice); //collision
						p4d->incrementCollidN();
					}

				} //endif real iterator

			region->insertSpace(p4d);
			noPoints++;
			if (noPoints >= Settings::Saving::savePointsFrequency) {
				this->snl->appendSpacePoints(rnode);
				noPoints = 0;
			}
		}
	}

	if (noPoints != 0) {
		this->snl->appendSpacePoints(rnode);
		noPoints = 0;
	}

	rnode->setFoundGoodOrientation(!noGoodOrientation);

	if (!Settings::AtlasBuilding::stop) {
		rnode->setComplete(true);
		region->trim();
	}

	delete detector;
	delete desc;
	delete chart;

	if (this->verbose)
		cout << "AtlasBuilder::MySample: Finished node " << from << endl;

	return !noGoodOrientation; //should delete this node and the path. // no do not delete, keep it for reference not to sample it again

}

AtlasNode* AtlasBuilder::createChildNode(AtlasNode *rnode,
		ActiveConstraintGraph *child_CG, Orientation* wit_orr_toSend,
		bool & savedOnce, bool leafWitness, bool & noGoodOrientation,
		int & noPoints, bool bret) {
	int from = rnode->getID();

	int child_nodeID;
	int success = this->atlas->addNode(child_CG, child_nodeID);
	AtlasNode* child_node = (*this->atlas)[child_nodeID];

	bool contin = false;
	if (success == 0)   //existed before : may be done or incomplete
			{
		delete child_CG;
		child_CG = NULL;
		child_CG = child_node->getCG(); //you need existing CG to get parameter etc. information
		if (!child_node->isComplete()) //incomplete  //&& !Settings::AtlasBuilding::breadthFirst
			contin = true;

	}

	wit_orr_toSend->addBoundary(from);

//	if(!contin)  //contin=true does not mean they are connected, maybe child_node is created by another node than 'rnode'
	this->atlas->connect(from, child_nodeID);

	if (success == 1 || (contin && !Settings::AtlasBuilding::breadthFirst)) //during 'breathFirst', we create child node but do not sample it. Sampling will be done after parent node is completed. Hence if contin is true, then that means child node is already created then no need to have recursive call.
			{
		if (!savedOnce) //add once
		{
			this->snl->appendSpacePoints(rnode); //save them before going recursion
			noPoints = 0;
			savedOnce = true;
		}

		//MAIN RECURSIVE CALL
		if (this->sampleTheNode(child_node, false, wit_orr_toSend, contin,
				bret)) //if there is a realization with good angle at the childs, then set noGoodOrientation to false
			noGoodOrientation = false;

	} else   //done
	{
		if (Settings::Saving::saveWitnessToFinalChild && leafWitness) //contactList2 being empty means this is the finalChild //do not add witness to every child, it is alot of duplicates
			if (!wit_orr_toSend->angleViolated()) {
				//if(!child_CG->isWitness(from))
				{
					CayleyPoint* wp4 = new CayleyPoint(wit_orr_toSend,
							child_CG->getParamLines(), a, b); //copy of wit_orr_toSend is added to wp4
					this->snl->appendWitness(child_node, wp4);
					wp4->trim_PointMultiD();
					delete wp4;

					child_node->setFoundGoodOrientation(true);
				}
			}

		if (child_node->hasAnyGoodOrientation()) //if the child has any good orientation, then parent should be displayed.
			noGoodOrientation = false;
	}

	return child_node;

}

/*
 * This method creates all subset of permutations NOT combinations. The permutation is necessary in order to connect all parents to all child nodes.
 * 1.tree = a -> ab   2.tree = b -> ba        creating the node ba from the 2. tree is necessary to have a connection between b and ba even though ba=ab is created in the 1. tree
 * But be careful to not save witness to all permutations !!!!!
 * Hence saveWitnessToFinalChild once by only from the first path!
 */
void AtlasBuilder::createChildNodes(AtlasNode *rnode_prime,
		list<pair<int, int> > contactList, Orientation* ori, Atlas* atlas_prime,
		SaveLoader* snlr, bool firstPath) {

	return;

	if (Settings::AtlasBuilding::breadthFirst
			|| !Settings::RootNodeCreation::createChildren)
		return;

	ActiveConstraintGraph *cgKprime = rnode_prime->getCG();
	if (cgKprime->constraintSize() >= 6 || Settings::AtlasBuilding::stop) //base case  // NOT ready to create contactGraph more than 6 contacts yet, you need change all description class for that. Also it causes a lot of empty nodes to be created, which are combinations of over constraints.!
		return;

	int from = rnode_prime->getID();
	if (from == -1)
		cout << "ERRORRR createChildNodes" << endl;

	ActiveConstraintGraph *cgKchilds[contactList.size()];
	for (int c = 0; c < contactList.size(); c++)  //TRY ALL CONTACT LIST
			{

		cgKchilds[c] = new ActiveConstraintGraph(cgKprime);
		pair<int, int> contact = contactList.front();
		cgKchilds[c]->addContact(contact);
		contactList.pop_front();

		if (cgKchilds[c]->isDependent()
				|| (cgKchilds[c]->constraintSize() >= 6
						&& !cgKchilds[c]->mayHaveThirdAtom())) {
			delete cgKchilds[c];
			cgKchilds[c] = NULL;
			contactList.push_back(contact);
			continue;
		}

		int child_nodeID;
		int success = atlas_prime->addNode(cgKchilds[c], child_nodeID);
		AtlasNode* child_node = (*atlas_prime)[child_nodeID];
		if (success == 0) { //exists
			delete cgKchilds[c];
			cgKchilds[c] = NULL;
			cgKchilds[c] = child_node->getCG();
			atlas_prime->connect(from, child_nodeID);
		} else {

			CayleyParameterization *desc = new CayleyParameterization(
					cgKchilds[c]); //set parameters to compute witness point correctly.  // to make it run faster, calls the 2. constructor of description class which is designed to set only parameters.

			atlas_prime->connect(from, child_nodeID);

			if (!Settings::Statistics::createPseudoAtlas) // while creating toy pseudo atlas, saving job is done later for efficiency
			{
				snlr->saveNodeHeader(child_node, atlas_prime); // Saving node is also done after sampling is finished to have more accurate data about completeness and emptiness.
				snlr->appendDimension(child_node);
			}

			delete desc;
			desc = NULL;
		}

		//add witness
		Orientation* orr = new Orientation(ori);
		orr->addBoundary(from);
		if (orr != NULL && !orr->angleViolated()) {
			child_node->setFoundGoodOrientation(true); //even you do not save witness to this node and even there is no other point other than witness, say that there is a good point in this node, because in reality there is a good witness, but you did not save it here and you saved that witness to grandchild, because it had more contacts.

			if (Settings::Saving::saveWitnessToFinalChild && firstPath
					&& (contactList.size() == 0 || child_node->getDim() == 0)) // TO SAVE WITNESS POINT 'ONLY' AT THE first leaf NODE!
					{
				CayleyPoint* wp4 = new CayleyPoint(orr,
						cgKchilds[c]->getParamLines(), snlr->a, snlr->b); //since this->a, b is not static, i used snlr

				if (!Settings::Statistics::createPseudoAtlas) {
					snlr->appendWitness(child_node, wp4); // save it to the file of graph
					wp4->trim_PointMultiD();
					delete wp4;
				} else {
					child_node->getACR()->insertWitness(wp4); // add it to the region,  while creating toy pseudo atlas, saving job is done later for efficiency
				}
			}
		}
		delete orr;
		orr = NULL; //since findPoint creates a new orientation (which is copy of orr) then adds it to wp4, you need to delete this orr!

		createChildNodes(child_node, contactList, ori, atlas_prime, snlr,
				firstPath); //i assumed contactList passed by value, i mean another list is created with same contacts, so the changes to contactList in the recursive call will not affect current contactList
		firstPath = false;
		contactList.push_back(contact);

	}

}

void AtlasBuilder::addReverseWitness(AtlasNode *rnode,
		ActiveConstraintGraph *cgK, ConvexChart *chart,
		ActiveConstraintRegion *region, Orientation* coming_witness_ori,
		bool & noGoodOrientation) {
	if (Settings::General::reverseWitness) {
		vector<CgMarker::Edge> edges = this->cgmarker.getEdge(cgK);
		cout << "[mark]" << *cgK << endl;
		cout << "[mark] getting witness " << edges.size() << endl;
		for (int i = 0; i < edges.size(); i++) {
			Orientation* ornt = edges[i].second; //DO NOT FORGET TO DELETE ORNT
			CayleyPoint* wp4 = new CayleyPoint(ornt, cgK->getParamLines(), a,
					b);
			chart->setWitnessPoint(wp4);
			if (!coming_witness_ori->angleViolated()) {
				region->insertWitness(wp4);
				//				cgK->insertWitness(wp4); //  add it to the graph, it will be needed as a starting point in the description class
				this->snl->appendWitness(rnode, wp4); // save it too the graph
				noGoodOrientation = false;
			}
		}
		if (cgK->getDim() == 0) {
			this->cgmarker.mark(rnode); //cgK
			cout << "[mark] mark 0" << endl;
		}
	}
}

void AtlasBuilder::refineMap() {
	vector<AtlasNode*> nodes = this->atlas->getNodes();
	/*
	 * load all 0-d nodes with their space points and mark their ancestor
	 * the space points are needed to serve as the starting point of the re-sampling
	 * TODO: actually only 1 space point is needed, is there a way to avoid loading the whole node?
	 */
	for (vector<AtlasNode*>::iterator iter = nodes.begin(); iter != nodes.end();
			iter++) {
		if ((*iter)->getDim() == 0) {
			snl->loadNode((*iter)->getID(), (*iter)->getACR());
			// ActiveConstraintGraph* cg = (*iter)->getCG();
			cgmarker.mark((*iter));
		}
	}

	while (!cgmarker.empty()) {
		pair<ActiveConstraintGraph*,
				vector<pair<ActiveConstraintGraph*, Orientation*> > > res =
				cgmarker.pop();
		ActiveConstraintGraph* cg = res.first;
		vector<pair<ActiveConstraintGraph*, Orientation*> > &edges = res.second;

		int nodenum = this->atlas->getNodeNum(cg);
		cout << *(cg) << endl;
		if (nodenum == -1) {
			delete cg;
			continue;
		}

		cout << "refining " << nodenum << endl;

		// check each child, if it is not in the roadmap or the connection doesn't exist, redo sample
		for (int i = 0; i < edges.size(); i++) {
			cout << endl << i << "th child" << endl;
			cout << *(res.second[i].first) << endl;
			if (atlas->getNodeNum(res.second[i].first) == -1
					|| !atlas->isConnected(nodenum,
							atlas->getNodeNum(res.second[i].first))) {

				AtlasNode * rnod = (*this->atlas)[nodenum];
				snl->loadNode(nodenum, rnod->getACR());

				Orientation* orr = res.second[i].second;

				// debug output, print the Orientation
				if (true) {
					cout << endl;
					cout << "o";

					vector<int> boundary = orr->getBoundary();
					int flip = orr->getFlipNum();

					double fb[3][3], tb[3][3];
					orr->getFromTo(fb, tb);
					for (int i = 0; i < 3; i++)
						for (int j = 0; j < 3; j++) {
							cout << " " << fb[i][j];
						}
					for (int i = 0; i < 3; i++)
						for (int j = 0; j < 3; j++) {
							cout << " " << tb[i][j];
						}
					cout << " " << boundary.size();
					for (int i = 0; i < boundary.size(); i++)
						cout << " " << boundary[i];
					cout << " " << flip;
					cout << endl;

				}

				// re-sample the node with the new starting Orientation
				sampleTheNode(rnod, false, edges[i].second, false, false);

				break;
			}
		}
		delete cg;

	}
}

list<Orientation*> AtlasBuilder::findRealizations(ActiveConstraintGraph *cgK,
		ConvexChart* des) {

	list<Orientation*> output;
	num_samples++;

	if (des->partial3tree)
		for (int i = 0; i < 8; i++) //for each flip
				{
			bool fail;
			Orientation *relz = CartesianRealizer::computeRealization(cgK, des,
					i, fail);
			if (!fail) {
				// We do not do angle check here i.e. if( !relz->angleViolated() )
				// Because after binary search around this orientation with bad angle, there may have orientations with good angle
				// Also child nodes can cause small angles so if this realization cause a new contact, then create the child not
				// and keep it as witness at the child node. After that delete this realization from parent(this node), and make it orange.
				output.push_back(relz);

			} else { //volume negative

				delete relz;
				relz = NULL;
				break; //in the next flips it will always get volume negative too
			}
		}

#ifdef USE_MATLAB
	//RUIJIN
	if(!des->partial3tree && Settings::General::runSolver) {
		bool fail;
		num_samples++;
		Orientation *relz = CartesianRealizer::computeRealization(cgK, des, 0, fail); //the first realization created by first constructor, the rest is by second constructor
		if( !fail ) {	//means maple is able to return in specified time
			msol = relz->mapslns;
			helix_base* solver = relz->solver;
			vector<vector<double> > rts = relz->rts;

			if(!Settings::Constraint::checkAngle || relz->angleViolated() )//small anglee
			output.push_back( relz );
			else {
				//delete relz;
				output.push_back( relz );
			}

			for(int i=1; i<msol; i++)
			{
				CartesianRealizer *relz = new CartesianRealizer(des, solver, rts, i);
				if(!Settings::Constraint::checkAngle || relz->angleViolated() ) //small anglee
				output.push_back( relz );
				else {
					//		delete relz;
					output.push_back( relz );
				}
			}

		}
		else
		delete relz;
	}
#endif

	return output;
}
