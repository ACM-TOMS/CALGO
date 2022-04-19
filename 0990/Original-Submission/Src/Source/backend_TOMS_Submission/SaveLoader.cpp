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
 * SaveLoader.cpp
 *
 *  Created on: Apr 12, 2009
 *      Author: James Pence
 */

#include "SaveLoader.h"

#include "ActiveConstraintGraph.h"
#include "CayleyPoint.h"
#include "Orientation.h"
#include "PointSet.h"

#include <string>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <limits>

#include <time.h>
using namespace std;

time_t SaveLoader::last_saved_time = time(NULL);

SaveLoader::SaveLoader(string relativePosition, PointSet *a, PointSet *b) {
	if (relativePosition.at(relativePosition.length() - 1) != '/')
		relativePosition.append("/");
	this->relativePath = relativePosition;
#ifdef WIN32
	mkdir(this->relativePath.c_str());
#else
	mkdir(this->relativePath.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#endif
	this->a = a;
	this->b = b;

}

SaveLoader::~SaveLoader() {
}

void SaveLoader::saveNodeHeader(AtlasNode* node, Atlas* atlas) {

	ofstream outFile;
	stringstream name;
	name << this->relativePath << "RoadMap.txt";
	name.flush();
	outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);

	writeNodeHeader(node, outFile);

	time_t current_time = time(NULL); //have static start time variable and local end time, do save if 5 minutes passed.
	int time_elapsed = difftime(current_time, last_saved_time);
	if (time_elapsed > 180) {  //3 minutes
		last_saved_time = current_time;
		saveRoadMap(atlas); //!!!duplicates current file too, so maybe it is inefficient to call it in every connection
	}
}

void SaveLoader::writeNodeHeader(AtlasNode* node, ostream &outFile) {

	ActiveConstraintGraph* cgk = node->getCG();

	outFile << "# " << node->getID() << " " << node->isComplete() << " "
			<< node->hasAnyGoodOrientation() << " ";

	vector<pair<int, int> > contacts = cgk->getParticipants();
	outFile << contacts.size() << " ";
	for (vector<pair<int, int> >::iterator iter = contacts.begin();
			iter != contacts.end(); iter++)
		outFile << "c " << iter->first << " " << iter->second << " ";

	vector<pair<int, int> > par = cgk->getParamLines();
	outFile << par.size() << " ";
	for (vector<pair<int, int> >::iterator iter = par.begin();
			iter != par.end(); iter++)
		outFile << "p " << iter->first << " " << iter->second << " ";

	outFile << cgk->getStepSize() << " ";
	outFile << node->getLocation()[0] << " " << node->getLocation()[1] << " "
			<< node->getLocation()[2] << " ";
	vector<int> con = node->getConnection();
	outFile << con.size();
	for (vector<int>::iterator it = con.begin(); it != con.end(); it++) {
		outFile << " " << *it;
	}
	outFile << "\n";

}

void SaveLoader::saveRecentPointsToFile(AtlasNode* rnode) {

	stringstream name;
	name << this->relativePath << "Node" << rnode->getID() << ".txt";
	name.flush();

	ofstream outFile;
	outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);

	vector<CayleyPoint*> witspc = rnode->getACR()->getWitness();
	for (vector<CayleyPoint*>::iterator iter = witspc.begin();
			iter != witspc.end(); iter++) {
		writeCayleyPoint(*iter, outFile, true);
	}
	vector<CayleyPoint*> spc = rnode->getACR()->getJustSpace();
	for (vector<CayleyPoint*>::iterator iter = spc.begin(); iter != spc.end();
			iter++) {
		writeCayleyPoint(*iter, outFile);
	}

	outFile.flush();
	outFile.close();

	rnode->getACR()->trim();

}

void SaveLoader::appendSpacePoints(AtlasNode* rnode) {

	cout << "."; // a progress indicator on the console.
	cout.flush();
	ofstream outFile;
	stringstream name;
	name << this->relativePath << "Node" << rnode->getID() << ".txt";
	name.flush();

	outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);
	vector<CayleyPoint*> space = rnode->getACR()->getJustSpace();

	for (size_t i = 0; i < space.size(); i++) {
		if (space[i] == NULL)
			continue;
		this->writeCayleyPoint(space[i], outFile);
	}
	outFile.flush();
	outFile.close();
	rnode->getACR()->trim();
}

void SaveLoader::appendWitness(AtlasNode* rnode, CayleyPoint* wit) {

	ofstream outFile;

	stringstream name;
	name << this->relativePath << "Node" << rnode->getID() << ".txt";
	name.flush();
	outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);

	this->writeCayleyPoint(wit, outFile, true);

	outFile.flush();
	outFile.close();

}

void SaveLoader::appendDimension(AtlasNode* rnode) {

	if (rnode->dimWritten)  // 6 February 2016
		return;
	else
		rnode->dimWritten = true;

	ofstream outFile;
	stringstream name;
	name << this->relativePath << "Node" << rnode->getID() << ".txt";
	name.flush();
	outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);
	outFile << "d " << rnode->getParamDim() << "\n"; //endl;
	outFile.flush();
	outFile.close();
}

void SaveLoader::appendVolume(AtlasNode* rnode) {
	double min[6], max[6];
	rnode->getACR()->getSpaceVolume(min, max);
	ofstream outFile;
	stringstream name;
	name << this->relativePath << "Node" << rnode->getID() << ".txt";
	name.flush();
	outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);
	outFile << "v ";
	for (size_t i = 0; i < 6; i++) {
		outFile << min[i] << " ";
	}
	for (size_t i = 0; i < 6; i++) {
		outFile << max[i];
		if (i < 5) {
			outFile << " ";
		}
	}
	outFile << "\n"; //endl;
	outFile.flush();
	outFile.close();
}

ofstream* SaveLoader::getFileFromBuffer(int nodenum) {
	for (list<pair<ofstream*, int> >::iterator iter = openFiles.begin();
			iter != openFiles.end(); iter++) {
		if ((*iter).second == nodenum)
			return (*iter).first;
	}

	stringstream name;
	name << this->relativePath << "Node" << nodenum << ".txt";
	name.flush();

	pair<ofstream*, int> newFile;
	newFile.first = new ofstream((name.str()).c_str(),
			ios_base::app | ios_base::out);
	newFile.second = nodenum;

	openFiles.push_front(newFile);
	if (openFiles.size() > 300) {
		ofstream * outFileLast = openFiles.back().first;
		(*outFileLast).flush();
		(*outFileLast).close();
		delete outFileLast;
		openFiles.pop_back();
	}

	return openFiles.front().first;
}

void SaveLoader::loadNode(int num, ActiveConstraintRegion* output) {

	ifstream inFile;
	stringstream name;
	name << this->relativePath << "Node" << num << ".txt";
	name.flush();

	cout << "SaveLoader::loadNode: Reading node " << num << " from "
			<< name.str() << endl;
	try {
		inFile.open((name.str()).c_str(), ios_base::in);
		if (inFile.is_open()) {
			loadActiveConstraintRegion(inFile, output);
		} else {
			cout << "SaveLoader::loadNode: Failed to open " << name.str()
					<< endl;
		}
		inFile.close();
	} catch (string e) {
		cerr << e;
		inFile.close();
	} catch (exception &e) {
		cerr << e.what();
		inFile.close();
	}
}

void SaveLoader::saveRoadMap(Atlas* atlas) {

	cout << ":";
	cout.flush();
	ofstream outFile;
	stringstream name;
	name << this->relativePath << "RoadMap.txt";
	name.flush();
	for (size_t i = 5; i > 0; i--) {
		stringstream oldername, newername;
		oldername << this->relativePath << "Roadmap" << i << ".txt";
		oldername.flush();
		newername << this->relativePath << "Roadmap" << (i - 1) << ".txt";
		newername.flush();
		if (i > 1) {
			remove(oldername.str().c_str());
			rename(newername.str().c_str(), oldername.str().c_str());
		} else {
			remove(oldername.str().c_str());
			rename(name.str().c_str(), oldername.str().c_str());
		}
	}

	outFile.open(name.str().c_str());
	if (!outFile.is_open())
		cout << "SaveLoader::saveMapView: Failed to create " << name.str()
				<< endl;
	outFile
			<< "/ NodeID Complete nonEmpty ContactSize \"contacts # #\"  ParamDimension \"parameters # #\"  StepSize \"Location x y z\" NumberOfConnections \"Nodes this node is connected to # # # ...\" "
			<< "\n"; //endl;

	size_t size = atlas->number_of_nodes();
	for (size_t iter = 0; iter < size; iter++) //since roadmap is extending during saving process by caller thread,  it is possible that there will be connection nodes, where the node itself is not saved yet.
			{
		AtlasNode* node = (*atlas)[iter];
		writeNodeHeader(node, outFile);
	}

	outFile.flush();
	outFile.close();

}

void SaveLoader::saveAtlas(Atlas* atlas) {

	size_t size = atlas->number_of_nodes();
	for (size_t num = 0; num < size; num++) {

		AtlasNode* node = (*atlas)[num];
		appendDimension(node);

		if (node->getACR()->getCombinedSize() > 0)
			saveRecentPointsToFile(node);
	}
}

void SaveLoader::loadMapView(Atlas* output) {
	bool debug = false;
	if (debug)
		cout << "loadMapView: Begin." << endl;

	if (output != NULL) {
		if (output->number_of_nodes() > 0) {
			output->cleanAtlas();
			cerr << "deleted nodes ";
		}
	} else
		output = new Atlas();

	ifstream inFile;
	stringstream name;
	vector<AtlasNode*> nodes;
	string trash;

	name << this->relativePath << "RoadMap.txt";
	name.flush();
	inFile.open(name.str().c_str());
	inFile >> trash;

	if (debug)
		cout << "loadMapView: Opening file " << name.str() << endl;

	while ((trash[0] == '#' || trash[0] == '/') && !inFile.eof()) {
		if (debug)
			cout << "Trash: " << trash << endl;

		if (trash[0] == '#') {
			AtlasNode* tmp = loadNextAtlasNode(inFile);
			nodes.push_back(tmp);
		} else {
			string trsh;
			getline(inFile, trsh);
		}
		inFile >> trash; //read #
	}
	inFile.close();

	output->setNodes(nodes);

}

AtlasNode* SaveLoader::loadAtlasNode(int nodenum) {

	AtlasNode* tmp;
	ifstream inFile;
	stringstream name;
	string trash;
	int id;
	bool found = false;
	name << this->relativePath << "RoadMap.txt";
	name.flush();
	inFile.open((name.str()).c_str());
	if (inFile.is_open()) {
		inFile >> trash;

		while ((trash[0] == '#' || trash[0] == '/') && !inFile.eof()) {
			if (trash[0] == '#') {
				inFile >> id;
				if (id == nodenum) {
					inFile.unget(); //to make loadAtlasNode able to read the id again
					tmp = loadNextAtlasNode(inFile);
					found = true;
					break;
				} else
					getline(inFile, trash);
			} else
				getline(inFile, trash);

			inFile >> trash; //read #
		}
		inFile.close();
	} else {
		cout << "SaveLoader::loadAtlasNode: Could not open roadmap file "
				<< name.str() << endl;
	}


	if (!found)
		tmp = new AtlasNode();
	return tmp;
}

AtlasNode* SaveLoader::loadNextAtlasNode(ifstream& inFile) {

	string trash;
	int id;
	bool complete;
	bool empty;
	int contactSize;
	vector<pair<int, int> > partic;
	int paramDim;
	vector<pair<int, int> > params;

	double stepsize;
	double* location = new double[3];
	int numberOfConnections;
	vector<int> con;

	inFile >> id >> complete >> empty >> contactSize;

	for (int i = 0; i < contactSize; i++) {
		inFile >> trash; //read 'c'
		pair<int, int> contact;
		inFile >> contact.first >> contact.second;
		partic.push_back(contact);
	}

	inFile >> paramDim;
	for (int i = 0; i < paramDim; i++) {
		inFile >> trash; //read 'p'
		pair<int, int> par;
		inFile >> par.first >> par.second;
		params.push_back(par);
	}

	inFile >> stepsize;
	inFile >> location[0] >> location[1] >> location[2];

	inFile >> numberOfConnections;
	for (int it = 0; it < numberOfConnections; it++) {
		int val;
		inFile >> val;
		con.push_back(val);
	}

	AtlasNode *tmp = new AtlasNode(id, complete, empty, 6 - contactSize,
			location, con);

	ActiveConstraintGraph * output = new ActiveConstraintGraph(partic, this->a,
			this->b, params); //16SEP13

	if (empty)
		tmp->setFoundGoodOrientation(true);
	else
		tmp->setFoundGoodOrientation(false);

	output->setStepSize(stepsize);

	tmp->setCG(output); //necessary since i need contact information to see if a new cgk is sampled before or not (need to go through whole cgks)

	delete[] location;

	return tmp;
}

void SaveLoader::readContactGraph0(ifstream &file,
		ActiveConstraintRegion *output) {
	// get length of file:
	file.seekg(0, file.end);
	int length = file.tellg();
	file.seekg(0, file.beg);

	int bufferSize = length; //8192;
	double iterations = length * 1. / bufferSize;
	cout << "iterations " << iterations << " length " << length << endl;

	char * buffer;
	for (int iter = 0; iter < iterations; iter++) {
		buffer = new char[bufferSize];  //sizeof(buffer)
		file.read(buffer, length);  // read data as a block:
		stringstream block;
		block << buffer;
		loadActiveConstraintRegion(block, output);

		delete[] buffer;

	}
}

void SaveLoader::loadActiveConstraintRegion(istream &file,
		ActiveConstraintRegion *output) {

	time_t Start_t = time(NULL);

	int oris_in_sparseregion = 0, totaloris = 0;
	bool debug = false;
	string element;
	CayleyPoint *currentPoint = NULL;
	CayleyPoint *prevPoint = NULL;
	vector<CayleyPoint*> witnesses;
	int paramDim;
	size_t line = 0;
	try {
		while (!file.eof()) {
			line++;

			if (file.fail()) { //can happen if the previous file reading operation "file >> something" cause error such as mismatch of types
				throw string("previous");
			}

			file >> element;
			if (file.eof())
				break;

			if (file.fail()) {
				throw string("read element");
			}

			switch (element[0]) {

			case 'd': //param dimension
			{
				if (debug)
					cout << "d";
				file >> paramDim;
			}
				break;

			case '*': //Point
			{
				if (debug)
					cout << "*";

				currentPoint = readCayleyPoint(file, paramDim);
				if (currentPoint != NULL) {
					output->insertSpace(currentPoint);
				}

			}
				break;

			case 'w': // witness point
			{
				currentPoint = NULL; //set it null otherwise when there file read error, and it crashes, you are gonna delete previous valid point
				if (debug)
					cout << "w";
				try {
					currentPoint = readCayleyPoint(file, paramDim);
				} catch (string err) {
					file.clear();
					string trash;
					getline(file, trash);
					cerr << "\n[" << line << "]Bad witness " << err << "\n";
					delete currentPoint;
					currentPoint = NULL;
					break;
				}
				if (currentPoint != NULL) {
//							output->insertWitness(currentPoint);
					witnesses.push_back(currentPoint);
				}

			}
				break;

			case 'o': // orientation
			{
				totaloris++;
				if (debug)
					cout << "o";
				Orientation* temp = NULL;
				try {
					temp = readOrientation(file);
				} catch (string err) {
					file.clear();
					cerr << "\n[" << line << "]Bad orientation " << err << "\n";
					file.ignore(numeric_limits<streamsize>::max(), '\n');
					if (temp != NULL) {
						delete temp;
						temp = NULL;
					}
					break;
				}
				if ((temp != NULL) && (currentPoint != NULL)) {
					currentPoint->addOrientation(temp);
				} else {
					delete temp;
					temp = NULL;
					cerr << "\n[" << line
							<< "] orientation deleted because of bad point \n";
				}
			}
				break;
			case '/': {
				if (debug)
					cout << "/";
				file.ignore(numeric_limits<streamsize>::max(), '\n');
			}
				break;

				//for old type node files
			case 'f':
			case 'r':
			case 'b':
			case 'c':
			case 'p':
			case 't':
			case 's':
			case 'v':
			case '-':
				file.ignore(numeric_limits<streamsize>::max(), '\n');
				break;

			default:
				cout << "\nline = \"" << line << " element " << element << "\""
						<< "\n"; //endl;
				cout.flush();
			}
			if (debug)
				cout.flush();
		}
		for (size_t i = 0; i < witnesses.size(); i++) {
			output->insertWitness(witnesses[i]);
		}
	} catch (string err) {
		delete currentPoint;
		cerr << "\n" << line << "file error - " << err << endl;
	}

	cout << "totaloris " << totaloris << endl;

	time_t t1 = time(NULL);
	int time_1 = difftime(t1, Start_t);
	cout << "time_1 " << time_1 << endl;
}

CayleyPoint* SaveLoader::readCayleyPoint(istream &file, size_t dim) {
	CayleyPoint* output = NULL;
	vector<double> pos;
	bool t;
	int badAngleN, collideN;

	double tmp;
	for (size_t i = 0; i < dim; i++) {
		file >> tmp;
		if (file.fail())
			throw string("read pointpos");
		pos.push_back(tmp);
	}

	file >> t;

	int axis = 0; //for jacobian sampling
	file >> axis;

	string temp;
	file >> temp;
	if (temp[0] == 'o') //somehow dimension is 1 but that parameter is not written to the file !!!
			{
		badAngleN = t;
		t = tmp;
		file.unget();
	} else
		badAngleN = atoi(temp.c_str());

	file >> collideN;

	for (size_t i = 0; i < pos.size(); i++) {
		if (std::isnan(pos[i])) {
			throw(string("nan value in Pos"));
		}
	}
	output = new CayleyPoint(pos);

	output->setRealizable(t);
	output->setBadAngleN(badAngleN);
	output->setCollidN(collideN);

	output->zIndex = axis;

	return output;

}

Orientation* SaveLoader::readOrientation(istream &file) {

	double fb[3][3], tb[3][3];
	string strValue;
	bool good = true;

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			file >> strValue;
			if (strValue == "nan") {
				good = false;
			} else {
				fb[i][j] = atof(strValue.c_str());
			}
		}

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			file >> strValue;
			if (strValue == "nan") {
				good = false;
			} else {
				tb[i][j] = atof(strValue.c_str());
			}
		}

	if (file.fail()) {
		throw string("read orient fromto");
	}
	int boundsize, flip;
	vector<int> boundary;
	file >> boundsize;
	int bndry;
	for (int i = 0; i < boundsize; i++) {
		file >> bndry;
		boundary.push_back(bndry);
	}
	file >> flip;

	if (file.fail()) {
		throw string("read orient boundflip");
	}

	Orientation* output = NULL;
	if (good) {
		output = new Orientation(fb, tb);
		output->setBoundary(boundary);
		output->setFlipNum(flip);
	}

	return output;

}

void SaveLoader::writeCayleyPoint(CayleyPoint* pnt, ostream &file,
		bool witness) {

	if (witness) {
		file << "w";
	} else {
		file << "*";
	}
	size_t dim;

	for (size_t i = 0; i < pnt->dim(); i++) {
		file << " " << (*pnt)[i];
	}

	file << " " << pnt->isRealizable();

	file << " " << pnt->zIndex; //for jacobian sampling

	file << " " << pnt->getBadAngleN();
	file << " " << pnt->getCollidN();

	file << "\n"; //endl;
	vector<Orientation*> ornt = pnt->getOrientations();
	for (vector<Orientation*>::iterator iter = ornt.begin(); iter != ornt.end();
			iter++) {
		writeOrientation(*iter, file);

	}
}

void SaveLoader::writeOrientation(Orientation* ornt, ostream &file) {

	int flip = ornt->getFlipNum();

	file << "o";

	vector<int> boundary = ornt->getBoundary();

	double fb[3][3], tb[3][3];
	ornt->getFromTo(fb, tb);

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			file << " " << fb[i][j];
		}

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			file << " " << tb[i][j];
		}
	file << " " << boundary.size();
	for (int i = 0; i < boundary.size(); i++)
		file << " " << boundary[i];
	file << " " << flip;   //the value is incorrect if it witness orientation.
	file << "\n"; // endl;

}

void SaveLoader::writeOrientationInPDBformat(PointSet *hA, PointSet *hB,
		Orientation *orien, string filename) {


	double fb[3][3], tb[3][3];
	orien->getFromTo(fb, tb);
	vector<Point*> atms = hA->getAtoms();
	vector<Point*> btms = hB->getXFAtoms(fb, tb);

	ofstream file;
	file.open(filename.c_str());
	for (size_t i = 0; i < atms.size(); i++) {
		file << atms[i]->getLine() << "\n"; //endl;
		cout << atms[i]->getLine() << "\n"; //endl;
	}
	for (size_t i = 0; i < btms.size(); i++) {
		file << btms[i]->getLine() << "\n"; //endl;
		cout << btms[i]->getLine() << "\n"; //endl;
	}
	file.flush();
	file.close();

	for (size_t i = 0; i < btms.size(); i++) {
		delete btms[i];
	}
}
