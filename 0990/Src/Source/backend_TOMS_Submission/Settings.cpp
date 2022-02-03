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
#include "Settings.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <cstdlib> // MUST come before simpleini
using namespace std;
using namespace Settings;

#include "simpleini/SimpleIni.h"
#include "ThreadShare.h"

std::string ModeNames[] = { "Stopped", "Sample Current Tree", "Cleanup",
		"Auto-Solve", "Breadth First Sampling", "Refine Sampling" };

namespace Settings {
namespace PointSetA {
std::string file = "../files/demo/molecular_unit_A.pdb";
std::vector<int> ignored_rows;
int x_col = 6;
int y_col = 7;
int z_col = 8;
int radius_col = 9;
int label_col = 3;
int pointNo_col = 1;
}

namespace PointSetB {
std::string file = "../files/demo/molecular_unit_B.pdb";
std::vector<int> ignored_rows;
int x_col = 6;
int y_col = 7;
int z_col = 8;
int radius_col = 9;
int label_col = 3;
int pointNo_col = 1;
}

namespace DistanceData {
std::string file = "../files/source_files/union computed desired distances.txt";
std::vector<int> ignored_rows;
int label1_col = 0;
int label2_col = 1;
int radius_col = 2;
int radiusMin_col = -1; //the column number for min radius
int radiusMax_col = -1; //the column number for max radius
}

namespace Output {
std::string dataDirectory = "../data/";
}

namespace General {
std::string Status("We display the status here");
bool candidate_interactions = false; // bool   virus  = false;

// no_gui

bool reverseWitness = false;
}

namespace RootNodeCreation {
bool createChildren = true;
int dimension_of_rootNodes = 5;

int participatingAtomIndex_low = 0;
int participatingAtomIndex_high = 0;
bool useParticipatingPointZDistance = false;
double ParticipatingPointZDistance = 7;
bool reversePairDumbbells = false;
double initial4DContactSeparation_low = 1.8;
double initial4DContactSeparation_high = 7.2;
}

namespace Sampling {
bool runSample = true;
Modes runMode;

double GridXY = 26;
double GridZ = 7; 

double stepSize = 0.2;
bool short_range_sampling = false;
bool dynamicStepSizeAmong = false;
int dynamicStepSizeWithin = 0;
bool sampleAllNodes;
bool binarySearch = false;

bool JacobianSampling = false;
std::array<double, 6> gridSteps = { 1, 1, 1, 10, 0.044658199, 10 };
}

namespace Constraint {

bool wholeCollision = false;

double activeLowerLambda = 0.8;
double activeLowerDelta = 0;
double activeUpperLambda = 1;
double activeUpperDelta = 1;

double collisionLambda = 0.8; 
double collisionDelta = 0;

double angleLow = 0;
double angleHigh = 30;
}

namespace AtlasBuilding {
bool stop = false;
bool breadthFirst = false;
bool parameterMinDeviation = false;

bool ifBadAngleWitness_createChild = false;
}

namespace Saving {
int savePointsFrequency = 10000;
bool saveWitnessToFinalChild = true;
bool saveBoundary = false;

bool saveOutGridOrientation = false;
}

namespace Statistics {
bool createPseudoAtlas = false;
}

namespace Paths {
int pathLength;
int energyLevelLowerBound;
int energyLevelUpperBound;
}

}

namespace ThreadShare {
SaveLoader *save_loader;
AtlasBuilder *atlas_builder;
Atlas *atlas_view;
}

std::string GetValue_PruneComments(CSimpleIniA &si, const char *category,
		const char *key) {
	const char *val = si.GetValue(category, key, NULL);
	if (val == NULL) {
		cout << "Setting \"[" << category << "] " << key << "\" is not present." << endl;
		exit(0);
	}

	std::string ret(val);
	int i;
	for (i = 0; i < ret.size() && ret[i] != ';'; i++)
		;
	ret.erase(i, ret.size() - i);
	return ret;
}

double getDouble(CSimpleIniA &si, const char *category, const char *key) {
	std::string val = GetValue_PruneComments(si, category, key);
	return atof(val.c_str());
}

int getInt(CSimpleIniA &si, const char *category, const char *key) {
	std::string val = GetValue_PruneComments(si, category, key);
	return atoi(val.c_str());
}

bool getBool(CSimpleIniA &si, const char *category, const char *key) {
	std::string val = GetValue_PruneComments(si, category, key);
	if (val.find("true") != val.npos)
		return true;
	return false;
}

std::string getString(CSimpleIniA &si, const char *category, const char *key) {
	std::string val = GetValue_PruneComments(si, category, key);

	size_t quote;
	quote = val.find('"');
	val.erase(0, quote + 1);

	quote = val.rfind('"');
	val.erase(quote, val.size() - quote);

	return val;
}

std::vector<double> getVectorDouble(CSimpleIniA &si, const char *category,
		const char *key) {
	std::string val = GetValue_PruneComments(si, category, key);

	size_t brace;
	brace = val.find('{');
	val.erase(0, brace + 1);

	brace = val.rfind('}');
	val.erase(brace, val.size() - brace);

	std::vector<double> ret;

	//determine if string is empty
	int i;
	for (i = 0; i < val.size() && val[i] != ' '; i++)
		;
	if (i == val.size())
		return ret;

	// read in array
	size_t comma;
	while ((comma = val.find(',')) != val.npos) {
		std::string num = val.substr(0, comma - 1);
		ret.push_back(atof(num.c_str()));
		val.erase(0, comma + 1);
	}
	ret.push_back(atof(val.c_str()));

	return ret;
}

/// gets one value when empty!!!!

std::vector<int> getVectorInt(CSimpleIniA &si, const char *category,
		const char *key) {
	std::string val = GetValue_PruneComments(si, category, key);

	size_t brace;
	brace = val.find('{');
	val.erase(0, brace + 1);

	brace = val.rfind('}');
	val.erase(brace, val.size() - brace);

	std::vector<int> ret;

	// read in array
	size_t comma;
	while ((comma = val.find(',')) != val.npos) {
		std::string num = val.substr(0, comma - 1);
		ret.push_back(atoi(num.c_str()));
		val.erase(0, comma + 1);
	}

	// check if there is any remaining element
	int i;
	for (i = 0; i < val.size() && val[i] != ' '; i++)
		;
	if (i != val.size())
		ret.push_back(atoi(val.c_str()));

	return ret;
}

#include <fstream>
bool Settings::load(const char *filename) {

	CSimpleIniA ini;
	ini.SetUnicode();
	if (ini.LoadFile(filename) < 0) {
		cout << "Failed to open file." << endl;
	}

	if (ini.IsEmpty()) {
		cout << "Empty ini." << endl;
	}

	///////////////////////
	// MolecularUnitA
	///////////////////////
	PointSetA::file = getString(ini, "PointSetA", "file");
	PointSetA::ignored_rows = getVectorInt(ini, "PointSetA", "ignored_rows");
	;
	PointSetA::x_col = getInt(ini, "PointSetA", "x_col");
	PointSetA::y_col = getInt(ini, "PointSetA", "y_col");
	PointSetA::z_col = getInt(ini, "PointSetA", "z_col");
	PointSetA::radius_col = getInt(ini, "PointSetA", "radius_col");
	PointSetA::label_col = getInt(ini, "PointSetA", "label_col");
	PointSetA::pointNo_col = getInt(ini, "PointSetA", "pointNo_col");

	///////////////////////
	// MolecularUnitB
	///////////////////////
	PointSetB::file = getString(ini, "PointSetB", "file");
	PointSetB::ignored_rows = getVectorInt(ini, "PointSetB", "ignored_rows");
	;
	PointSetB::x_col = getInt(ini, "PointSetB", "x_col");
	PointSetB::y_col = getInt(ini, "PointSetB", "y_col");
	PointSetB::z_col = getInt(ini, "PointSetB", "z_col");
	PointSetB::radius_col = getInt(ini, "PointSetB", "radius_col");
	PointSetB::label_col = getInt(ini, "PointSetB", "label_col");
	PointSetB::pointNo_col = getInt(ini, "PointSetB", "pointNo_col");

	///////////////////////
	// DistanceData
	///////////////////////
	DistanceData::file = getString(ini, "DistanceData", "file");
	DistanceData::ignored_rows = getVectorInt(ini, "DistanceData",
			"ignored_rows");
	DistanceData::label1_col = getInt(ini, "DistanceData", "label1_col");
	DistanceData::label2_col = getInt(ini, "DistanceData", "label2_col");
	DistanceData::radius_col = getInt(ini, "DistanceData", "radius_col");
	DistanceData::radiusMin_col = getInt(ini, "DistanceData", "radiusMin_col");
	DistanceData::radiusMax_col = getInt(ini, "DistanceData", "radiusMax_col");

	///////////////////////
	// Output
	///////////////////////
	Output::dataDirectory = getString(ini, "Output", "dataDirectory");

	///////////////////////
	// General
	///////////////////////

	General::candidate_interactions = getBool(ini, "General",
			"candidate_interactions");

	General::reverseWitness = getBool(ini, "General", "reverseWitness");

	///////////////////////
	// RootNodeCreation
	///////////////////////
	RootNodeCreation::createChildren = getBool(ini, "RootNodeCreation",
			"createChildren");
	RootNodeCreation::dimension_of_rootNodes = getInt(ini, "RootNodeCreation",
			"dimension_of_initialContactGraphs");
	RootNodeCreation::useParticipatingPointZDistance = getBool(ini,
			"RootNodeCreation", "useParticipatingPointZDistance");
	RootNodeCreation::ParticipatingPointZDistance = getDouble(ini,
			"RootNodeCreation", "ParticipatingPointZDistance");
	RootNodeCreation::reversePairDumbbells = getBool(ini, "RootNodeCreation",
			"reversePairDumbbells");
	RootNodeCreation::initial4DContactSeparation_low = getDouble(ini,
			"RootNodeCreation", "initial4DcontactSeparationLow");
	RootNodeCreation::initial4DContactSeparation_high = getDouble(ini,
			"RootNodeCreation", "initial4DcontactSeparationHigh");

	///////////////////////
	// Sampling
	///////////////////////
	Sampling::runSample = getBool(ini, "Sampling", "runSample");
	Sampling::GridXY = getDouble(ini, "Sampling", "GridXY");
	Sampling::GridZ = getDouble(ini, "Sampling", "GridZ");

	Sampling::stepSize = getDouble(ini, "Sampling", "stepSize");
	Sampling::short_range_sampling = getBool(ini, "Sampling", "sixDimensions");
	Sampling::dynamicStepSizeAmong = getBool(ini, "Sampling", "dynamicStepSizeAmong");
	Sampling::dynamicStepSizeWithin = getInt(ini, "Sampling", "dynamicStepSizeWithin");
	Sampling::sampleAllNodes = getBool(ini, "Sampling", "sampleAllNodes");
	Sampling::binarySearch = getBool(ini, "Sampling", "binarySearch");

	Sampling::JacobianSampling = getBool(ini, "Sampling", "JacobianSampling");
	std::vector<double> gridSteps_ = getVectorDouble(ini, "Sampling",
			"gridSteps");
	Sampling::gridSteps[0] = gridSteps_[0];
	Sampling::gridSteps[1] = gridSteps_[1];
	Sampling::gridSteps[2] = gridSteps_[2];
	Sampling::gridSteps[3] = gridSteps_[3];
	Sampling::gridSteps[4] = gridSteps_[4];
	Sampling::gridSteps[5] = gridSteps_[5];

	///////////////////////
	// Constraint
	///////////////////////

	Constraint::wholeCollision = getBool(ini, "Constraint", "wholeCollision");

	Constraint::activeLowerLambda = getDouble(ini, "Constraint",
			"activeLowerLambda");
	Constraint::activeLowerDelta = getDouble(ini, "Constraint",
			"activeLowerDelta");

	Constraint::activeUpperLambda = getDouble(ini, "Constraint",
			"activeUpperLambda");
	Constraint::activeUpperDelta = getDouble(ini, "Constraint",
			"activeUpperDelta");

	Constraint::collisionLambda = getDouble(ini, "Constraint",
			"collisionLambda");
	Constraint::collisionDelta = getDouble(ini, "Constraint", "collisionDelta");

	Constraint::angleLow = getDouble(ini, "Constraint", "angleLow");
	Constraint::angleHigh = getDouble(ini, "Constraint", "angleHigh");

	///////////////////////
	// AtlasBuilding
	///////////////////////
	AtlasBuilding::stop = getBool(ini, "AtlasBuilding", "stop");
	AtlasBuilding::breadthFirst = getBool(ini, "AtlasBuilding", "breadthFirst");
	AtlasBuilding::parameterMinDeviation = getBool(ini, "AtlasBuilding",
			"parameterMinDeviation");

	AtlasBuilding::ifBadAngleWitness_createChild = getBool(ini, "AtlasBuilding",
			"ifBadAngleWitness_createChild");

	///////////////////////
	// Saving
	///////////////////////
	Saving::savePointsFrequency = getInt(ini, "Saving", "savePointsFrequency");
	Saving::saveWitnessToFinalChild = getBool(ini, "Saving",
			"saveWitnessToFinalChild");
	Saving::saveBoundary = getBool(ini, "Saving", "saveBoundary");

	///////////////////////
	// Statistics
	///////////////////////
	Statistics::createPseudoAtlas = getBool(ini, "Statistics",
			"createPseudoAtlas");

	///////////////////////
	// Paths
	///////////////////////

	Paths::pathLength = getInt(ini, "Paths", "path_length");
	Paths::energyLevelLowerBound = getInt(ini, "Paths",
			"energy_level_lower_bound");
	Paths::energyLevelUpperBound = getInt(ini, "Paths",
			"energy_level_upper_bound");

	return true;
}

void writeCommentLine(ofstream &outfile, const char *comment) {
	outfile << "; " << comment << endl;
}

void writeNewLines(ofstream &outfile, unsigned int lines) {
	for (int i = 0; i < lines; i++)
		outfile << endl;
}

void writeVarCommentAndEndl(ofstream &outfile, const char *comment) {
	if (strcmp(comment, "") != 0)
		outfile << " ; " << comment;
	outfile << endl;
}

void writeSection(ofstream &outfile, const char *section) {
	outfile << "[" << section << "]" << endl;
}

// string
void writeVar(ofstream &outfile, const char *key, std::string val,
		const char *comment = "") {
	outfile << key << " = \"" << val << "\"";
	writeVarCommentAndEndl(outfile, comment);
}

// bool
void writeVar(ofstream &outfile, const char *key, bool val,
		const char *comment = "") {
	outfile << key << " = " << (val ? "true" : "false");
	writeVarCommentAndEndl(outfile, comment);
}

// int
void writeVar(ofstream &outfile, const char *key, int val, const char *comment =
		"") {
	outfile << key << " = " << val;
	writeVarCommentAndEndl(outfile, comment);
}

// double
void writeVar(ofstream &outfile, const char *key, double val,
		const char *comment = "") {
	outfile << key << " = " << val;
	// outfile << key << " = " << std::showpoint << val;
	// outfile << key << " = " << std::setprecision(1) << val;
	writeVarCommentAndEndl(outfile, comment);
}

// vector<int>
void writeVar(ofstream &outfile, const char *key, std::vector<int> val,
		const char *comment = "") {
	outfile << key << " = {";
	if (val.size() != 0) {
		for (int i = 0; i < val.size() - 1; i++)
			outfile << val[i] << ", ";
		outfile << val[val.size() - 1];
	}
	outfile << "}";
	writeVarCommentAndEndl(outfile, comment);
}

// vector<double>
void writeVar(ofstream &outfile, const char *key, std::vector<double> val,
		const char *comment = "") {
	outfile << key << " = {";
	if (val.size() != 0) {
		for (int i = 0; i < val.size() - 1; i++)
			outfile << val[i] << ", ";
		outfile << val[val.size() - 1];
	}
	outfile << "}";
	writeVarCommentAndEndl(outfile, comment);
}

// vector<double>
void writeVar(ofstream &outfile, const char *key, std::array<double, 6> val,
		const char *comment = "") {
	outfile << key << " = {";
	for (int i = 0; i < 5; i++)
		outfile << val[i] << ", ";
	outfile << val[5];
	outfile << "}";
	writeVarCommentAndEndl(outfile, comment);
}

bool Settings::save(const char *filename) {
	ofstream outfile;

	outfile.open(filename);
	if (!outfile.is_open()) {
		cout << "File \"" << filename << "\" could not be opened." << endl;
		return false;
	}

	writeCommentLine(outfile, "Strings are enclosed in \"\"");
	writeCommentLine(outfile, "Bools are either \"true\" or \"false\"");
	writeCommentLine(outfile,
			"Arrays are enclosed in {} and values are separated by commas");

	writeNewLines(outfile, 2);

	///////////////////////
	// MolecularUnitA
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "PointSetA");
	writeVar(outfile, "file", PointSetA::file);
	writeVar(outfile, "ignored_rows", PointSetA::ignored_rows);
	writeVar(outfile, "x_col", PointSetA::x_col);
	writeVar(outfile, "y_col", PointSetA::y_col);
	writeVar(outfile, "z_col", PointSetA::z_col);
	writeVar(outfile, "radius_col", PointSetA::radius_col);
	writeVar(outfile, "label_col", PointSetA::label_col);
	writeVar(outfile, "pointNo_col", PointSetA::pointNo_col);

	///////////////////////
	// MolecularUnitB
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "PointSetB");
	writeVar(outfile, "file", PointSetB::file);
	writeVar(outfile, "ignored_rows", PointSetB::ignored_rows);
	writeVar(outfile, "x_col", PointSetB::x_col);
	writeVar(outfile, "y_col", PointSetB::y_col);
	writeVar(outfile, "z_col", PointSetB::z_col);
	writeVar(outfile, "radius_col", PointSetB::radius_col);
	writeVar(outfile, "label_col", PointSetB::label_col);
	writeVar(outfile, "pointNo_col", PointSetB::pointNo_col);

	///////////////////////
	// DistanceData
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "DistanceData");
	writeVar(outfile, "file", DistanceData::file);
	writeVar(outfile, "ignored_rows", DistanceData::ignored_rows);
	writeVar(outfile, "label1_col", DistanceData::label1_col);
	writeVar(outfile, "label2_col", DistanceData::label2_col);
	writeVar(outfile, "radius_col", DistanceData::radius_col);
	writeVar(outfile, "radiusMin_col", DistanceData::radiusMin_col);
	writeVar(outfile, "radiusMax_col", DistanceData::radiusMax_col);

	///////////////////////
	// Output
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "Output");
	writeVar(outfile, "dataDirectory", Output::dataDirectory+"/");

	///////////////////////
	// General
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "General");

	writeVar(outfile, "candidate_interactions",
			General::candidate_interactions);
	// General::no_gui
	writeVar(outfile, "reverseWitness", General::reverseWitness);

	///////////////////////
	// RootNodeCreation
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "RootNodeCreation");
	writeVar(outfile, "createChildren", RootNodeCreation::createChildren);
	writeVar(outfile, "dimension_of_initialContactGraphs",
			RootNodeCreation::dimension_of_rootNodes);
	writeVar(outfile, "useParticipatingPointZDistance",
			RootNodeCreation::useParticipatingPointZDistance);
	writeVar(outfile, "ParticipatingPointZDistance",
			RootNodeCreation::ParticipatingPointZDistance);
	writeVar(outfile, "reversePairDumbbells",
			RootNodeCreation::reversePairDumbbells);
	writeVar(outfile, "initial4DcontactSeparationLow",
			RootNodeCreation::initial4DContactSeparation_low);
	writeVar(outfile, "initial4DcontactSeparationHigh",
			RootNodeCreation::initial4DContactSeparation_high);
	///////////////////////
	// Sampling
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "Sampling");
	writeVar(outfile, "runSample", Sampling::runSample);
	writeVar(outfile, "GridXY", Sampling::GridXY);
	writeVar(outfile, "GridZ", Sampling::GridZ);

	writeVar(outfile, "stepSize", Sampling::stepSize);
	writeVar(outfile, "sixDimensions", Sampling::short_range_sampling);
	writeVar(outfile, "dynamicStepSizeAmong", Sampling::dynamicStepSizeAmong);
	writeVar(outfile, "dynamicStepSizeWithin", Sampling::dynamicStepSizeWithin);
	writeVar(outfile, "sampleAllNodes", Sampling::sampleAllNodes);
	writeVar(outfile, "binarySearch", Sampling::binarySearch);

	writeVar(outfile, "JacobianSampling", Sampling::JacobianSampling);
	writeVar(outfile, "gridSteps", Sampling::gridSteps);

	///////////////////////
	// Constraint
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "Constraint");
	writeVar(outfile, "wholeCollision", Constraint::wholeCollision);

	writeVar(outfile, "activeLowerLambda", Constraint::activeLowerLambda);
	writeVar(outfile, "activeLowerDelta", Constraint::activeLowerDelta);

	writeVar(outfile, "activeUpperLambda", Constraint::activeUpperLambda);
	writeVar(outfile, "activeUpperDelta", Constraint::activeUpperDelta);

	writeVar(outfile, "collisionLambda", Constraint::collisionLambda);
	writeVar(outfile, "collisionDelta", Constraint::collisionDelta);

	writeVar(outfile, "angleLow", Constraint::angleLow);
	writeVar(outfile, "angleHigh", Constraint::angleHigh);

	///////////////////////
	// AtlasBuilding
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "AtlasBuilding");
	writeVar(outfile, "stop", AtlasBuilding::stop);
	writeVar(outfile, "breadthFirst", AtlasBuilding::breadthFirst);
	writeVar(outfile, "parameterMinDeviation",
			AtlasBuilding::parameterMinDeviation);

	writeVar(outfile, "ifBadAngleWitness_createChild",
			AtlasBuilding::ifBadAngleWitness_createChild);

	///////////////////////
	// Saving
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "Saving");
	writeVar(outfile, "savePointsFrequency", Saving::savePointsFrequency);
	writeVar(outfile, "saveWitnessToFinalChild",
			Saving::saveWitnessToFinalChild);
	writeVar(outfile, "saveBoundary", Saving::saveBoundary);

	///////////////////////
	// Statistics
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "Statistics");
	writeVar(outfile, "createPseudoAtlas", Statistics::createPseudoAtlas);

	///////////////////////
	// Pahts
	///////////////////////
	writeNewLines(outfile, 2);
	writeSection(outfile, "Paths");
	writeVar(outfile, "path_length", Paths::pathLength);
	writeVar(outfile, "energy_level_upper_bound", Paths::energyLevelUpperBound);
	writeVar(outfile, "energy_level_lower_bound", Paths::energyLevelLowerBound);

	outfile.close();
	return true;
}
