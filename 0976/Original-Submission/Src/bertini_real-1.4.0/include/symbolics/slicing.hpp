#pragma once

#include "io/fileops.hpp"
#include "nag/witness_set.hpp"
#include "symbolics/derivative_systems.hpp"
/**
 \brief Write a Bertini input file which slices at input linears.
 
 \param input_file The name of the input file to which we are appending functions.
 \param output_file The name of the output file, which will have the sliced system in it.
 \param linears The linear functions to slice with.
 \param num_to_add The number of linears we are adding on.  This should match the number of vec_mp in linears.
 \param W The input witness set.
 */
void create_sliced_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
													vec_mp * linears, int num_to_add,
													const WitnessSet & W);

