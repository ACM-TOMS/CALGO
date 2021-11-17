#ifndef DERIVATIVE_SYSTEMS_H
#define DERIVATIVE_SYSTEMS_H

/** \file derivative_systems.hpp */


/**
 \defgroup filecreation File Creation
 
 Methods and functions for creating and parsing Bertini input files.
 */


#include "bertini1/bertini_extensions.hpp"
#include "programConfiguration.hpp"
#include "io/fileops.hpp"



/**
 \brief used in parse_names, this function increments and appends items found to numItems, itemNames, and itemLines.
 
 \ingroup filecreation
 
 \param numItems counter for number of items found
 \param itemNames The found names of the items
 \param itemLines array of integer indices of lines on which the items are found
 \param IN the open file
 \param lineNumber the number of the line on which we are.
 */
void addItems(int *numItems, char ***itemNames, int **itemLines, FILE *IN, int lineNumber);


/**
 \brief Find the lines where items are declared and find names
 
 \ingroup filecreation
 
 \param numItems The number of items found.
 \param itemNames The names of the items found.
 \param itemLines The integer numbers of the lines on which the items are found.
 \param IN the open file in which to search.
 \param name The declaration to look for; e.g., "constants"
 \param num_declarations The number of declarations to find. This is previously established using partitionParse.
 */
void parse_names(int *numItems, char ***itemNames, int **itemLines, FILE *IN, char *name, int num_declarations);


/**
 \brief Get a string containing all the constants in an input file.
 
 \ingroup filecreation
 
 \return String containing the constants as parsed out by this function.
 \param filename The name of the file to parse.
 \param numConstants the number of constants to search for.
 \param consts an array of arrays of strings containing the names of the constants.
 \param lineConstants The number of the lines containing the constants.
 */
std::string just_constants(boost::filesystem::path filename, int numConstants, char **consts, int *lineConstants);


/**
 \brief Write a matrix to an open file as constants, including the declaration and definition.
 
 \ingroup filecreation
 
 \param M The matrix to write
 \param prefix String prefix to which to append numbers to declare the matrix.
 \param OUT an already open file pointer.
 */
void write_matrix_as_constants(mat_mp M, std::string prefix, FILE *OUT);



/**
 \brief Write a vector to an open file as constants, including the declaration and definition.
 
 \ingroup filecreation
 
 \param V The matrix to write
 \param prefix String prefix to which to append numbers to declare the vector.
 \param OUT an already open file pointer.
 */
void write_vector_as_constants(vec_mp V, std::string prefix, FILE *OUT);


#endif
