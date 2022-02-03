#pragma once

#include "io/fileops.hpp"
#include "nag/witness_set.hpp"
#include "symbolics/derivative_systems.hpp"
/**
 \brief Write a Bertini input file corresponding to the intersection of the input system with a sphere of a given center and radius.
 
 \param input_file The name of the input file to which we are appending functions.
 \param output_file The name of the output file, which will have the sphere system in it.
 \param sphere_radius The radius of the sphere.
 \param sphere_center The center of the sphere
 \param W the input witness set.
 */
void create_sphere_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
													comp_mp sphere_radius,
													vec_mp sphere_center,
													const WitnessSet & W);
													