#ifndef _SAMPLER_H
#define _SAMPLER_H


/** \file sampler.hpp */


#include <boost/timer/timer.hpp>

#include "io/fileops.hpp"

#include "programConfiguration.hpp"
#include "decompositions/checkSelfConjugate.hpp"

#include "nag/solvers/multilintolin.hpp"
#include "decompositions/surface.hpp"
#include "decompositions/curve.hpp"









/**
 \defgroup samplermethods Sampler Methods
	
 Methods common to sampling Bertini_real decompositions.
*/



/**
 \brief light startup processing common to all dimensions, for the sampler.
 
 parse the input file, get the tracker_config_t, set up the PPD, adjust a few sampler settings.
 
 \ingroup samplermethods
 
 \param D the base-class-type Decomposition to refine.
 \param sampler_options The current state of sampler.
 \param solve_options the current state of the solver configuration.
 */
void common_sampler_startup(const Decomposition & D,
							sampler_configuration & sampler_options,
							SolverConfiguration & solve_options);




/**
 \brief get the MPType, name of the directory to sample, and the dimension of the Decomposition.
 
 \todo replace this function which reads Decomposition metadata from a file in the directory.
 
 \return the MPType
 \param Dir_Name the name read in by this function
 \param MPType apparently this function returns this in two ways.
 \param dimension The retrieved dimension.
 */
int get_dir_mptype_dimen(boost::filesystem::path & Dir_Name, int & MPType, int & dimension);




















/**
 \brief Set the linear and point in a witness set to the input L and pts.
 
 \todo remove this function, or make a method of the WitnessSet class.
 
 \param W the witness set to modify.
 \param new_linear the linear to set
 \param new_point the input point 
 */
void set_witness_set_mp(WitnessSet & W, vec_mp new_linear,vec_mp new_point);







/**
 \brief given two points and a projection, estimate a projection value for the point halfway between.
 
 \param result the computed estimated projection value 
 \param left Input one.
 \param right Input two.
 \param pi the linear projection to use.
 */
void estimate_new_projection_value(comp_mp result, vec_mp left, vec_mp right, vec_mp pi);


/**
 \brief given two points and a projection, estimate a projection value for the point halfway between.
 
 \param result the computed estimated projection value
 \param estimated_point The computed estimated point.
 \param left point one
 \param right point two
 \param pi the linear projection to use.
 */
void estimate_new_projection_value(comp_mp result, vec_mp estimated_point, vec_mp left, vec_mp right, vec_mp pi);







/**
 \brief triangulate two ribs, each with at least two entries, by iterating from left to right, and always constructing the more equilateral Triangle of the two candidates.  
 
 the end of the loop simply constructs every Triangle between the two ribs until it reaches the end.  
 
 \param rib1 a rib of integer indices in a VertexSet.
 \param rib2 a rib adjacent to rib1, of integer indices in a VertexSet
 \param V the vertex set into which the ribs index.
 \param real_thresh The threshold of imaginary part, so that a point is thresholded to be real.
 \param current_samples The triangulation being built.
 */
void triangulate_two_ribs_by_angle_optimization(const std::vector< int > & rib1, const std::vector< int > & rib2,
												VertexSet & V, double real_thresh,
												std::vector< Triangle> & current_samples);



/**
 \brief compute square of difference between angle and \f$\pi/3\f$ radians.
 
 \param temp a temporary variable.  comp_mp's are expensive to initialize and clear, so this is passed in for optimization.
 \param length1 the length of one of the adjacent sides
 \param length2 the length of the other adjacent side
 \param dot_prod the dot product of the vectors representing the two adjacent legs of the Triangle.
 \return the square of the difference of the computed angle, and \f$\pi/3\f$.
 */
double compute_square_of_difference_from_sixtydegrees(comp_mp temp, comp_mp length1, comp_mp length2, comp_mp dot_prod);


void ScaleByCycleNum(comp_mp result, comp_mp input, int cycle_num_l, int cycle_num_r);








#endif




