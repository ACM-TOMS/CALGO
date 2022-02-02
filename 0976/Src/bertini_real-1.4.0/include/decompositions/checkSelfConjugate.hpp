#ifndef _CHECK_SELFCONJUGATE_H
#define _CHECK_SELFCONJUGATE_H


/** \file checkSelfConjugate.hpp */


#include "programConfiguration.hpp"
#include "io/fileops.hpp"
#include "bertini1/bertini_extensions.hpp"


/**
 \brief main function for performing self-conjugacy check.
 
 returns 1 if self-conjugate, 0 else.
 
 \param test_point The point to conjugate and test for membership.
 \param program_options The current program state.
 \param input_file the name of the bertini input file to parse out.
 
 \return boolean integer indicating whether the component is self-conjugate.
 
 */
bool checkSelfConjugate(vec_mp test_point,
						BertiniRealConfig & program_options,
						boost::filesystem::path input_file);




/**
 \brief returns the incidence number according to the incidence matrix
 
 \param test_point The point to feed into Bertini.
 \param program_options The current program state.
 \param input_file the name of the bertini input file
 
 \return the incidence number for the component, relative to the incidence_matrix file.
 */
int get_incidence_number(vec_mp test_point, BertiniRealConfig & program_options, boost::filesystem::path input_file);

/**
 \brief write a single point to "member_points"
 
 \param point_to_write the homogeneous point to write
 \return the number 0.  Why?
 */
int write_member_points_singlept(vec_mp point_to_write);

/**
 \brief write a single point, and its complex conjugate, to "member_points"
 
 \param point_to_write The point to write, along with its conjugate.
 \return the number 0.  Why?
 */
int write_member_points_sc(vec_mp point_to_write);

/**
 \brief write the input file to feed bertini to perform membership testing
 
 \param outputFile the name of the file to write
 \param funcInput the name of the func_input file
 \param configInput the name of the config file
 \param tracktype bertini track type.  should be 3.
 */
void membership_test_input_file(boost::filesystem::path outputFile,
                                boost::filesystem::path funcInput,
                                boost::filesystem::path configInput,
                                int  tracktype);

/**
 \brief read the incicence_matrix file.  return the incidence numbers for the member_points
 
 \todo check this function for correctness when the points satisfy multiple components (its on the intersection). The return type probably ought to be a vector of vectors of ints.
 
 \return A std::vector<int> of the indicence numbers for the test points.
 */
std::vector<int> read_incidence_matrix();





#endif
