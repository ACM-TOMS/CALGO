#ifndef _FILEOPS_H
#define _FILEOPS_H


/** 
 \file fileops.hpp 
 
 */

#include <iostream>
#include <fstream>
#include <sys/stat.h>        //  declare the 'stat' structure 
#include <sys/types.h>

#include <mpi.h> // this *cannot* be inside an extern "C"{} wrapper.
 
#include <dirent.h>
#include <sstream>
#include <set>



// #include "bertini1/bertini_headers.hpp"


extern "C" {
#include "partitionParse.h"
}


#include "boost/filesystem.hpp"



/**
\brief Enter a 1-second open-sleep cycle until the desired file is found.  Potentially an infinite loop if the file is never generated.

\param name The name of the file you want to wait on.
*/
void WaitOnGeneratedFile(const std::string & name);


/**
 \brief parse a bertini input file into two separate files, according to the two inputs; also counts the numbers of several types of declarations.
 
 A wrapper around partitionParse.  
 
 If sc_flag is 1, will parse variable_groups and functions *twice* adding bar for those which will correspond to the conjugated...
 
 \return The returned value of partitionParse.
 \param declarations pointer to an array of integers, containing the number of 9 types of declarations.
 \param input_filename The name of the input file to parse.
 \param functions_filename The name of the file into which we should put in the functions.
 \param config_filename The name into which to put the config part of the file
 \param not_sc_flag Flag indicating whether to run in sc mode.
 */
int partition_parse(int **declarations,
					boost::filesystem::path input_filename,
					boost::filesystem::path functions_filename,
					boost::filesystem::path config_filename,
					int not_sc_flag);


/**
 \brief Renames the arr.out deg.out num.out config and preproc_data files to *.bak.
 
 These files restored by restore_bertini_files_dotbak().
 */
void rename_bertini_files_dotbak();


/**
 \brief restores the .bak files which were renamed to *.bak by move_bertini_files_dotbak
 */
void restore_bertini_files_dotbak();

/**
 \brief purges an entire directory of all files except those which start with a period (.)
 
 \param directoryName the name of the directory to empty.
 */
void purge_previous_directory(char *directoryName);



/**
 opens a file to read, throwing runtime_error if cannot.
 
 \return A file pointer, just opened in 'r' mode.
 \param filename the name of the file to open.
 */
FILE *safe_fopen_read(boost::filesystem::path filename);

/**
 opens a file to write, throwing runtime_error if cannot.
 
 \return A file pointer, just opened in 'w' mode.
 \param filename the name of the file to open.
 */
FILE *safe_fopen_write(boost::filesystem::path filename);

/**
 opens a file to append, throwing runtime_error if cannot.
 
 \return A file pointer, just opened in 'a' mode.
 \param filename the name of the file to open.
 */
FILE *safe_fopen_append(boost::filesystem::path filename);




/**
 \brief copies a file character by character.
 
 \param input_file The name of the file to read from
 \param OUTfile The name of the file to copy TO.
 */
void copyfile(boost::filesystem::path input_file, boost::filesystem::path OUTfile);


/**
 \brief copies a file character by character.
 
 \param IN The already open file to read from.
 \param OUT The already open file to write to.
 */
void copyfile(FILE *IN,FILE *OUT);




/**
 \brief quit Bertini_real programs, but not gracefully at all.  This is inspired by Bertini's bexit.
 
 \param errorCode A numeric exit code.
 */
void br_exit(int errorCode);


/**
 \brief deliberately segfault by reading out of range.
 */
void deliberate_segfault();


/**
 \brief Have the user input a value untl it's an integer, return that value.
 
 \ingroup ui
 
 \return The integer the user inputted.
 */
int getInteger();




/**
 \brief Parse a string that has an integer value in string form.
 
 \ingroup ui
 
 
 \param text  The integer value as a string.
 \param results The value to set as a string.
 \return A boolean to indicate whether the parsing was successful or not.
 */
bool parseInteger( std::string const& text, int& results );

/**
 \brief Display a menu option to the user and ask for an integer input within the specified range.
 
 \ingroup ui
 
 \param display_string - The menu as a string.
 \param min_value The minimum value allowed.
 \param max_value The maximum value allowed.
 \return The integer the user specified.
 */
int get_int_choice(std::string display_string,int min_value,int max_value);


/**
 \brief Display a menu option to the user and ask for an integer input within the specified range.
 
 \ingroup ui
 
 \param display_string The menu as a string.
 \param valid_values  A std::set of valid integer values.  all others will be rejected.
 \return The integer the user specified.
 */
int get_int_choice(std::string display_string, const std::set<int> & valid_values);



#endif

