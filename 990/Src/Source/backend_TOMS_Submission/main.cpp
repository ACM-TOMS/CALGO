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
// source
#include "AtlasBuilder.h"
#include "Atlas.h"

#include "Settings.h"
#include "SaveLoader.h"
#include "readIn.h"

// standard
#include <iostream>
#include <map>

using namespace std;

// Only use once settings are loaded
void init_DistanceFinder_from_settings(PredefinedInteractions *df) {
	vector<vector<string> > data = readData(Settings::DistanceData::file);

	if (data.size() == 0)
		return;
	bool matrixInsteadOfColumns = false;
	//Check if the data is in matrix format
	if (data.front().size() == data.size()) {
		for (std::vector<string>::size_type iter = 0; iter != data.size();
				iter++) {
			if (data[iter][1] == data[1][iter])
				matrixInsteadOfColumns = true;
			else {
				matrixInsteadOfColumns = false;
				break;
			}
		}
	}

	PredefinedInteractions *output;
	if (matrixInsteadOfColumns) {
		output = new PredefinedInteractions(data);
	} else if (Settings::DistanceData::radiusMin_col >= 0
			&& Settings::DistanceData::radiusMax_col >= 0) {
		output = new PredefinedInteractions(data,
				Settings::DistanceData::label1_col,
				Settings::DistanceData::label2_col,
				Settings::DistanceData::radiusMin_col,
				Settings::DistanceData::radiusMax_col);
	} else if (Settings::DistanceData::radius_col >= 0) {
		output = new PredefinedInteractions(data,
				Settings::DistanceData::label1_col,
				Settings::DistanceData::label2_col,
				Settings::DistanceData::radius_col);
	} else if (Settings::DistanceData::radiusMin_col >= 0) {
		output = new PredefinedInteractions(data,
				Settings::DistanceData::label1_col,
				Settings::DistanceData::label2_col,
				Settings::DistanceData::radiusMin_col);
	} else {
		output = new PredefinedInteractions();
	}

	df->assign(output);
	delete output;
}

// Only use once settings are loaded and you got your df from init_global_DistanceFinder
PointSet* init_PointSet_A_from_settings(PredefinedInteractions *df) {
	PointSet *muA = readPointSetFromFile(Settings::PointSetA::file,
			Settings::PointSetA::ignored_rows, Settings::PointSetA::x_col,
			Settings::PointSetA::y_col, Settings::PointSetA::z_col,
			Settings::PointSetA::radius_col, Settings::PointSetA::label_col,
			Settings::PointSetA::pointNo_col);
	muA->buildStree();
	if (Settings::General::candidate_interactions)
		muA->simplify(*df, true);

	return muA;
}

// Only use once settings are loaded and you got your df from init_global_DistanceFinder
PointSet* init_PointSet_B_from_settings(PredefinedInteractions *df) {
	PointSet *muB = readPointSetFromFile(Settings::PointSetB::file,
			Settings::PointSetB::ignored_rows, Settings::PointSetB::x_col,
			Settings::PointSetB::y_col, Settings::PointSetB::z_col,
			Settings::PointSetB::radius_col, Settings::PointSetB::label_col,
			Settings::PointSetB::pointNo_col);
	muB->buildStree();
	if (Settings::General::candidate_interactions)
		muB->simplify(*df, true);

	return muB;
}

/**
 * Main entrance of the program.
 * It handles command-line arguments.
 * Starts sampling the atlas by making a call to startAtlasBuilding
 * All the necessary input will be read from a settings file.
 */
int main(int argc, char** argv) {

	// Declare objects to hold Point Sets
	PointSet *psA;
	PointSet *psB;

	// Declare objects to hold Point Sets
	SaveLoader *save_loader;
	AtlasBuilder *atlas_builder;

	Atlas *atlas = new Atlas();
	// end HEAP vars

	// command line arguments
	string settings_file = "./settings.ini";
#ifdef WIN32
	settings_file = "..\\settings.ini";
#endif

	// handle command-line arguments
	for (int i = 1; i < argc; i++) {
		string inp = argv[i];
		string flag, val;

		int eqsign = inp.find("=");
		if (eqsign + 1 == inp.length()) {
			cerr << "No value for flag " << inp.substr(0, eqsign) << "."
					<< endl;
			return 0;
		} else if (eqsign == -1) {
			flag = inp;
			val = "";
		} else {
			flag = inp.substr(0, eqsign);
			val = inp.substr(eqsign + 1, inp.length());
		}

		if (!flag.compare("-settings")) {
			if (!val.compare("")) {
				i++;
				if (i == argc) {
					cerr << "No filename provided for flag " << flag << "."
							<< endl;
					return 0;
				}
				settings_file = argv[i];
			} else {
				cerr << "Invalid value for flag " << flag << "." << endl;
				return 0;
			}
		}
	}

	// load the settings
	cout << "Loading settings from \"" << settings_file << "\"." << endl;
	mkdir(Settings::Output::dataDirectory.c_str(), S_IRWXU | S_IRWXG | S_IRWXO ); //the folder has to exist first
	Settings::load(settings_file.c_str());

	Settings::Sampling::runSample = true;

	// distance data
	PredefinedInteractions df;
	init_DistanceFinder_from_settings(&df);

	// PointSet A and B
	psA = init_PointSet_A_from_settings(&df);
	psB = init_PointSet_B_from_settings(&df);

	// SaveLoader
	save_loader = new SaveLoader(Settings::Output::dataDirectory, psA, psB);

	if (Settings::Constraint::wholeCollision) {
		// reading the neighbor matrix
		ConstraintCheck::nei_matrix = Utils::getMatrixFromFileFromTo("nei.txt");
	} else {
		ConstraintCheck::nei_matrix = Utils::getIdentityMatrix();
	}
	
	std::string s_file = Settings::Output::dataDirectory + "/settings.ini";
	Settings::save(s_file.c_str());

	atlas_builder = new AtlasBuilder(psA, psB, save_loader, &df, atlas);

	atlas_builder->setup();
	cout << "Thread_Main: AtlasBuilder Set up done." << endl;

	cout << "Thread_Main: Calling atlas_builder->startAtlasBuilding()." << endl;
	atlas_builder->startAtlasBuilding();

	cout << "Thread_Main: Calling this->save_loader->saveRoadMap(this->atlas)."
			<< endl;
	save_loader->saveRoadMap(atlas);

	cout << "Thread_Main: Finishes and Exits.." << endl;

	return 0; /* ANSI C requires main to return int. */
}
