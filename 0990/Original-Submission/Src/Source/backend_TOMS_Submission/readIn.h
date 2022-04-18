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
#ifndef READIN_H_
#define READIN_H_

#include "PointSet.h"

#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>

using namespace std;

/**
 * @brief Reads in a file containing a MolecularUnit. Neg col values result in
 * default values
 *
 * @return A pointer to the MU in the file, or NULL if something is wrong.
 */
inline PointSet* readPointSetFromFile(string filename, vector<int> ignored_rows,
		unsigned int x_col, unsigned int y_col, unsigned int z_col,
		int radius_col, int label_col, int atomNo_col) {
	PointSet *ret_mu = new PointSet();
	// open file
	ifstream file;
	file.open(filename.c_str());

	// check that the file was good.
	if (!file.is_open()) {
		cout << "File \"" << filename << "\" does not exist." << endl;
		delete ret_mu;
		return ret_mu;
	}

	// sort the ignored rows vector for simplified read in logic
	sort(ignored_rows.begin(), ignored_rows.end());

	// read file line by line
	for (int row = 0, ignore_idx = 0; !file.eof(); row++) {
		// we've encountered an ignored row
		if (ignore_idx < ignored_rows.size()
				&& row == ignored_rows[ignore_idx]) {
			ignore_idx++;
			continue;
		}

		// read in the line
		string line;
		getline(file, line);

		if (line.empty() || line == " ")
			continue;

		// read in each column
		vector<string> line_data;
		stringstream line_ss(line);
		string s;
		while (line_ss >> s)
			line_data.push_back(s);

		// read in column by column
		string label, atomNo, trash;
		double loc[3], radius;
		for (int col = 0; col < line_data.size(); col++) {
			stringstream col_ss(line_data[col]);

			if (col == atomNo_col)
				col_ss >> atomNo;
			else if (col == label_col)
				col_ss >> label;
			else if (col == x_col)
				col_ss >> loc[0];
			else if (col == y_col)
				col_ss >> loc[1];
			else if (col == z_col)
				col_ss >> loc[2];
			else if (col == radius_col)
				col_ss >> radius;
		}

		//default values
		if (radius_col < 0)
			radius = 1.2;
		if (label_col < 0)
			label = "";
		if (atomNo_col < 0)
			atomNo = "";

		// make and add atom to the molecular unit
		Point *a = new Point(loc, radius, label, atomNo);
		a->setLine(line_data, x_col, y_col, z_col);
		ret_mu->addAtom(a);

	}

	// close file
	file.close();

	// make it calc limits, then return
	ret_mu->calcLimits();
	return ret_mu;
}

/**
 * @brief Read a line of a file into a string stream.
 *
 * @param outputStringStream The stream to read into.
 * @param fileInput The file to read from
 * @return [description]
 */
inline bool readline(stringstream& outputStringStream, ifstream& fileInput) {
	string str;
	do {
		getline(fileInput, str);
	} while (!fileInput.eof() && (str.length() < 2 || isspace(str[0])));
	outputStringStream.clear();
	outputStringStream.str(str);
	outputStringStream.flush();
	if (fileInput.eof()) {
		outputStringStream.clear();
		return false;
	}
	return true;
}

/*
 * The data for the table format is read in with this function
 * creating a vector of rows, each represented as a vector.
 */
inline std::vector<std::vector<std::string> > readData(std::string filename) {
	ifstream file;
	stringstream stringRow;
	std::vector<std::vector<std::string> > output;
	size_t entriesPerRow = 0;
	file.open(filename.c_str());

	if (!file.is_open()) {
		cerr << "Unable to open - \"" << filename << "\"" << endl;
		return output;
	}

	while (readline(stringRow, file)) {
		std::vector<std::string> row;
		while (!stringRow.eof()) {
			std::string element;
			stringRow >> element;
			row.push_back(element);
		}
		output.push_back(row);
	}
	file.close();

	for (size_t y = 0; y < output.size(); y++) {
		if (entriesPerRow < output[y].size()) {
			entriesPerRow = output[y].size();
		}
	}

	for (size_t y = 0; y < output.size(); y++) {
		for (size_t x = output[y].size(); x < entriesPerRow; x++) {
			output[y].push_back("0");
		}
	}

	return output;
}

#endif
