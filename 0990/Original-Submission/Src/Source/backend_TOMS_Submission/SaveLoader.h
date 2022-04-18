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
#ifndef SAVELOADER_H_
#define SAVELOADER_H_

#include "Atlas.h"

#include <errno.h> //where EBUSY is defined.#include <list>
#include <utility>

#ifdef WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#ifdef __APPLE__
#include <libkern/OSAtomic.h>
#endif

class SaveLoader {
public:

	/**
	 * @brief Constructor with relativePath, helA and helB initialization
	 *
	 * @param relativePath The directory information on where the data files will be saved.
	 * @param a & b are provided so that contact graphs loaded later will have them defined.
	 */
	SaveLoader(string relativePath, PointSet *a, PointSet *b);

	/** @brief Destructor */
	virtual ~SaveLoader();

	/**
	 * @brief This function saves the identifying data of the node and its contact graph to the RoadMap.txt file.
	 * It is not saving the region ACR just basic information.
	 *
	 * Also in every 3 minutes, it saves atlas as fresh RoadMap with most updated information of nodes.
	 * Since some of the header information can change during the process such as
	 * Complete Empty NumberOfConnections and Nodes this node is connected to
	 */
	void saveNodeHeader(AtlasNode* node, Atlas* atlas);

	/**
	 * @brief Writes the following information of the 'node' to the 'outFile'
	 * NodeID Complete nonEmpty ContactSize \"contacts # #\"  ParamDimension \"parameters # #\" StepSize \"Location x y z\" NumberOfConnections \"Nodes this node is connected to # # # ...\"
	 * nonEmpty=0 means there is no realization with good angle, so the node will not be displayed, but still saved not to sample it again.
	 * previous name was writeNodeToRoadMap
	 */
	void writeNodeHeader(AtlasNode* node, ostream &outFile);

	/**
	 * @brief Saves node header of all nodes in the atlas to RoadMap.txt file
	 * Roadmap file is duplicated and previous copies are saved regularly.
	 * So during this saving process, even incomplete nodes are saved too.
	 * Data about completeness and emptiness .. are saved.
	 */
	void saveRoadMap(Atlas* atlas);

	/**
	 * @brief Loads Cayley space points to acr of the node identified by num
	 * @param num ID of the node
	 * @param acr Pointer to ActiveConstraintRegion of the node. Empty initially to be filled.
	 */
	void loadNode(int num, ActiveConstraintRegion* acr);

	/**
	 * @brief Loads a fresh atlas from the RoadMap.txt file
	 * Creates AtlasNodes from scratch by reading the headers of the node from RoadMap.txt file.
	 * ActiveConstraintGraph that labels the AtlasNode is also created since when a new node is
	 * about to be created, its label is checked to prevent duplications of the node and hence re-sampling.
	 *
	 * However ActiveConstraintRegion of the AtlasNode is not loaded to prevent memory blow up.
	 * It will be loaded whenever it is necessary. Such as whenever user click on a node to see its space.
	 *
	 * Instead of creating a new atlas, a pointer is passed to be filled in. So that Atlas is shared by backEnd and frontEnd
	 * in order to ensure that the AtlasBuilder class and viewer classes(GUI) are using the same atlas
	 *
	 * @see loadNextAtlasNode
	 */
	void loadMapView(Atlas* atlas);

	/**
	 * @brief Search the RoadMap.txt for the header of the node identified by nodenum.
	 * If found, creates an AtlasNode with those informations.
	 *
	 * @see loadNextAtlasNode
	 */
	AtlasNode* loadAtlasNode(int nodenum);

	/**
	 * @brief Loads the following header information from the 'inFile' and creates a fresh AtlasNode with its ActiveConstraintGraph.
	 * NodeID Complete nonEmpty ContactSize \"contacts # #\"  ParamDimension \"parameters # #\" StepSize \"Location x y z\" NumberOfConnections \"Nodes this node is connected to # # # ...\"
	 * nonEmpty=0 means there is no realization with good angle, so the node will not be displayed, but still saved not to sample it again.
	 */
	AtlasNode* loadNextAtlasNode(ifstream& inFile);

	/**
	 * @brief Writes all the atoms of hA and hB for a specific orient to the filename.
	 */
	void writeOrientationInPDBformat(PointSet *hA, PointSet *hB,
			Orientation *orien, string filename);

	/**
	 * @brief This function opens the node data file and adds the newly sampled data to it.
	 * i.e. Saves recent CayleyPoints with its orientations except witness points.
	 * Then trims those data from the ActiveConstraintRegion of the node.
	 * todo trims cleans witness points as well. make sure you have saved them before or you do not need them too.
	 */
	void appendSpacePoints(AtlasNode* node);

	/**
	 * @brief This function opens the node data file and adds the witness data to it.
	 * i.e. Saves witness CayleyPoint with its orientations.
	 *
	 * todo maybe you can delete wit after you write it to file in this method, check if it hurts?
	 */
	void appendWitness(AtlasNode* rnode, CayleyPoint* wit);

	/**
	 * @brief This function opens the node data file and writes ParamDim data following the letter d.
	 */
	void appendDimension(AtlasNode* rnode);

	/**
	 * @brief This function opens the node data file and writes the volume information.
	 * It would be useful if just triangular inequalities are used to determine the range of the parameters.
	 * i.e. When chart is not tight
	 */
	void appendVolume(AtlasNode* rnode);

	/**
	 * @brief This function opens the node data file and adds the newly sampled data to it.
	 * i.e. Saves recent CayleyPoints with its orientations  (witness points saved as well).
	 * Then trims those data from the ActiveConstraintRegion of the node.
	 */
	void saveRecentPointsToFile(AtlasNode* rnode);

	/**
	 * @brief For each AtlasNode of atlas, it saves recently sampled points (witness included) to its node.txt files.
	 */
	void saveAtlas(Atlas* atlas);

	/**
	 * @bried Loads CayleyPoints with its orientations from the file and save it to output.
	 *
	 * @see readCayleyPoint
	 * @see readOrientation
	 */
	void loadActiveConstraintRegion(istream &file,
			ActiveConstraintRegion *output);

	/**
	 * @brief Creates a fresh CayleyPoint with the following information read from the file
	 * w/*    \"parameter values x y z w ... \"    \"isRealizable 0/1 \"   zIndex  badAngleN collidN
	 *
	 * It does not loads set of orientations that belongs to CayleyPoint
	 *
	 * @param dim Number of parameters
	 */
	CayleyPoint* readCayleyPoint(istream &file, size_t dim);

	/**
	 * @brief Create a fresh orientation with the following information read from the file
	 * o fbx fby fbz  tbx tby tbz NumberOfBoundaries \"Boundaries this orientation takes place # # # ...\" Flipnumber
	 */
	Orientation* readOrientation(istream &file);

	/**
	 * @brief Write CayleyPoint information to the file
	 * w/*    \"parameter values x y z w ... \"    \"isRealizable 0/1 \"   zIndex  badAngleN collidN
	 * set of orientations
	 *
	 * @see writeOrientation
	 */
	void writeCayleyPoint(CayleyPoint* pnt, ostream &file,
			bool witness = false);

	/**
	 * @brief Write orientation information to the file
	 * o fbx fby fbz  tbx tby tbz NumberOfBoundaries \"Boundaries this orientation takes place # # # ...\" Flipnumber
	 */
	void writeOrientation(Orientation* ornt, ostream &file);

	//todo these methods is not used. check if it can be necessary in future, if not delete it.
	ofstream* getFileFromBuffer(int nodenum);
	std::list<std::pair<ofstream*, int> > openFiles; // <filehandler, nodenumber>
	void readContactGraph0(ifstream &file, ActiveConstraintRegion *output);

	/** The directory information on where the data files will be saved.*/
	string relativePath;

	PointSet *a, *b;
private:


	/** To keep track of time where in each 3 minutes roadmap.txt file will be rewritten. */
	static time_t last_saved_time;

};

#endif /* SAVELOADER_H_ */
