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
#ifndef SETTINGS_H_
#define SETTINGS_H_

#ifdef USE_MATLAB
#include "ruijin/helix_ve.hpp"
#include "ruijin/helix_ev.hpp"
#endif

#include <string>
#include <vector>
#include <array>
// using namespace std;

// inline bool generateDefaultSettings(string filename);

enum Modes {
	Stopped,
	TreeSample,
	ForestSample,
	ForestSampleAS,
	BreadthFirst,
	RefineSampling
};
extern std::string ModeNames[];

namespace Settings {
typedef unsigned int uint;

/**
 * @brief Loads Settings from file and puts them into the global variables
 * in this namespace. If a setting doesn't exit, the program terminates.
 *
 * @return True if file loaded. False otherwise.
 */
bool load(const char *filename);

bool save(const char *filename);

namespace PointSetA {
extern std::string file;
extern std::vector<int> ignored_rows;
extern int x_col;
extern int y_col;
extern int z_col;
extern int radius_col;
extern int label_col;
extern int pointNo_col;
}

namespace PointSetB {
extern std::string file;
extern std::vector<int> ignored_rows;
extern int x_col;
extern int y_col;
extern int z_col;
extern int radius_col;
extern int label_col;
extern int pointNo_col;
}

//todo rename distancedata
namespace DistanceData {
extern std::string file;
extern std::vector<int> ignored_rows;
extern int label1_col;
extern int label2_col;

//todo put a note saying that this bonding value will be used instead of the (r_i + r_j) in the main easal window
extern int radius_col;

//todo if min and max columns are set, then disable easal bonding threshold row in main input window.
extern int radiusMin_col; //min radius
extern int radiusMax_col; // max radius
}

namespace Output {
extern std::string dataDirectory;
}

namespace General {
// previous name was virus_mer
extern bool candidate_interactions; // allows only weak inter-atomic interface interactions.
extern std::string Status;
/**
 * For each of the 3 interfaces (2-fold, 3-fold and 5-fold), we determined the pairs
 * of interacting residues that are conserved in related viruses (10-20 pairs for each interface).
 * Only these were used as the candidate interactions for the assembly.
 */
// Settings::General::candidate_interactions is true iff PredefinedInteractions::distTableDefined is true
// no_gui
extern bool reverseWitness;
}

//todo in the display mention about root node when using following variables
namespace RootNodeCreation {

/** set it false JUST TO sample INTERIOR of root nodes */
extern bool createChildren; //todo put it to user input window or not???

extern int dimension_of_rootNodes;
//extern int    middleDumbbells_low;  // the lowest positioning of the one dumbel atom
//extern int    middleDumbbells_high; // the highest positioning of the one dumbel atom
extern int participatingAtomIndex_low; // the lowest positioning of the one dumbel atom
extern int participatingAtomIndex_high; // the highest positioning of the one dumbel atom
//extern bool   closeByDumbbells;
//extern double closeByDumbbellsAmount;
extern bool useParticipatingPointZDistance;
extern double ParticipatingPointZDistance;
extern bool reversePairDumbbells;
//extern double min;  // the min distance between two atoms of the dumbel inside one helix
//extern double max;  // the max distance between two atoms of the dumbel inside one helix
extern double initial4DContactSeparation_low;
extern double initial4DContactSeparation_high;
}

namespace Sampling {
extern Modes runMode;
extern bool runSample;  //set false  TO CONTINUE and load FROM PREVIOUS RUN
extern double GridXY;
extern double GridZ;

extern double stepSize;
extern bool short_range_sampling;
extern bool dynamicStepSizeAmong;
extern int dynamicStepSizeWithin;
extern bool sampleAllNodes;

extern bool binarySearch;

extern std::array<double, 6> gridSteps;
}

namespace Constraint {
//    extern bool   stericConstraint;
extern bool wholeCollision;

extern double activeLowerLambda;
extern double activeLowerDelta;
extern double activeUpperLambda;
extern double activeUpperDelta;

extern double collisionLambda;
extern double collisionDelta;

extern double angleLow;
extern double angleHigh;
}

namespace AtlasBuilding {
extern bool stop;
extern bool breadthFirst;
extern bool parameterMinDeviation; // todo put the flag to run this method or not into the settings.ini file (not to the GUI). Also remove this variable when submitting to TOMS???

extern bool ifBadAngleWitness_createChild; //todo put it to user input window or not???
}

namespace Saving {
extern int savePointsFrequency;
extern bool saveWitnessToFinalChild;
extern bool saveBoundary;

extern bool saveOutGridOrientation;
}

namespace Statistics {
extern bool run_statistics;

extern std::string folder1Location;
extern std::string folder2Location;
extern std::string folder3Location;
extern std::string gridLocation;
extern bool createPseudoAtlas;
}

namespace Paths {
extern int pathLength;
extern int energyLevelLowerBound;
extern int energyLevelUpperBound;
}

#ifdef USE_MATLAB
vector<helix_base*> solvers;
MatlabEngine *me;
#endif
}

#endif
