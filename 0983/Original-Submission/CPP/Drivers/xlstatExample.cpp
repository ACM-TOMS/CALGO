// This file is self sufficient and should be compiled as a usual cpp project.
//
// Compilation instruction: 
//	g++ -Wextra -pedantic -std=c++11 -o3 xlstatExaple.cpp xlstatExactCoq.cpp -o exactCochran

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "xlstatExactCoq.h"

int loadData(std::string &filename, matTshort &table)
{
	unsigned int szRow = 0;
    std::ifstream  data(filename);

    std::string line;
	vecTshort temp;

	table.clear();

	if (data.is_open())
	{
		while(std::getline(data, line))
		{
			std::stringstream  lineStream(line);
			std::string        cell;
			while(std::getline(lineStream, cell,','))
			{
				temp.push_back(atoi(cell.c_str()));
			}
			if (szRow > 0 && szRow != temp.size())
				return -1;
			table.push_back(temp);
			temp.clear();
		}
	}
	else 
		return -2;
	
	return 0;
}

/*****************************************************************************************************************************************************************
*** \brief  Calling function example
*** \author Arnaud Belletoile
****************************************************************************************************************************************************************/
int main()
{
	int retval = 0;
	matTshort table;
	double pVal = 0.;
	double timeLimit = 30;

	std::string filename;

	filename = "ex1data8x4.csv";
	retval = loadData(filename, table);
	if (retval < 0)
	{
		if (retval == -1)
			std::cout << "Input file does not contain a regular matrix. Program has to stop." << std::endl;
		else 
			std::cout << "Input file not found. Program has to stop." << std::endl;
		return -1;
	}

	std::cout << "*******************************************************************************" << std::endl;
	std::cout << "**** First Example (data source in file " << filename << ")" << std::endl;
	std::cout << "*******************************************************************************" << std::endl << std::endl;
	std::cout << "Entry matrix size is 4x4. This example is the one given in the first section of the publication." << std::endl;
	std::cout << "Expected p-value is : 0.073" << std::endl;

	if (xlstatCoqExact(table, pVal, timeLimit) >= 0)
		std::cout << "Computed p-value = " << pVal << std::endl; 
	else
		std::cout << "Time limit exceeded. The program had to stop." << std::endl;
	
	filename = "ex2data12x7.csv";
	retval = loadData(filename, table);
	if (retval < 0)
	{
		if (retval == -1)
			std::cout << "Input file does not contain a regular matrix. Program has to stop." << std::endl;
		else 
			std::cout << "Input file not found. Program has to stop." << std::endl;
		return -1;
	}

	std::cout << std::endl << std::endl;
	std::cout << "*******************************************************************************" << std::endl;
	std::cout << "**** Second Example (data source in file " << filename << ")" << std::endl;
	std::cout << "*******************************************************************************" << std::endl << std::endl;
	std::cout << "Entry matrix size is 12x7. This is a typical use case." << std::endl;
	std::cout << "Expected p-value is : 0.654" << std::endl;

	if (xlstatCoqExact(table, pVal, timeLimit) >= 0)
		std::cout << "Computed p-value = " << pVal << std::endl; 

	filename = "ex3data100x20.csv";
	retval = loadData(filename, table);
	if (retval < 0)
	{
		if (retval == -1)
			std::cout << "Input file does not contain a regular matrix. Program has to stop." << std::endl;
		else 
			std::cout << "Input file not found. Program has to stop." << std::endl;
		return -1;
	}

	std::cout << std::endl << std::endl;
	std::cout << "*******************************************************************************" << std::endl;
	std::cout << "**** Third Example (data source in file " << filename << ")" << std::endl;
	std::cout << "*******************************************************************************" << std::endl << std::endl;
	std::cout << "Entry matrix size is 100x20. This dataset is really big. The time limit (" << timeLimit << " s.) may be exceeded." << std::endl;

	if (xlstatCoqExact(table, pVal, timeLimit) >= 0)
		std::cout << "Computed p-value = " << pVal << std::endl; 

	return 0;
}
