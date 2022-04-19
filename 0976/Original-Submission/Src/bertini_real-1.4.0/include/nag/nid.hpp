#pragma once


#include "nag/witness_set.hpp"
#include "limbo.hpp"
#include "decompositions/checkSelfConjugate.hpp"

/**
 \brief metadata for witness points, for the NumericalIrreducibleDecomposition class.
 
 
 */
class WitnessPointMetadata
{

	int dimension_;
	
	int corank_, typeflag_, multiplicity_, component_number_, deflations_needed_;
	double condition_number_, smallest_nonzero_sing_value_, largest_zero_sing_value_;
	
	
	
	
public:
	
	
	/**
	 \brief get the largest zero singular value
	 \return largest zero singular value
	 */
	inline double largest_zero_sing_value() const
	{
		return largest_zero_sing_value_;
	}
	
	
	
	/**
	 \brief get the smallest nonzero singular value
	 \return smallest nonsingular value
	 */
	inline double smallest_nonzero_sing_value() const
	{
		return smallest_nonzero_sing_value_;
	}
	
	
	
	/**
	 \brief get the condition number
	 \return condition number
	 */
	inline double condition_number() const
	{
		return condition_number_;
	}
	
	

	
	/**
	 \brief get the number of deflations needed
	 \return num_deflations_needed
	 */
	inline int num_deflations_needed() const
	{
		return deflations_needed_;
	}
	
	
	
	/**
	 \brief get the multiplicity
	 \return multiplicity	 */
	inline int multiplicity() const
	{
		return multiplicity_;
	}
	
	
	
	
	/**
	 \brief get the type
	 \return type
	 */
	inline int type() const
	{
		return typeflag_;
	}
	
	/**
	 \brief get the corank
	 \return corank
	 */
	inline int corank() const
	{
		return corank_;
	}
	
	
	/**
	 \brief get the dimension
	 \return dimension
	 */
	inline int dimension() const
	{
		return dimension_;
	}
	
	
	
	/**
	 \brief get the component number
	 \return component number
	 */
	inline int component_number() const
	{
		return component_number_;
	}
	
	
	
	
	/**
	 read the metadata from the witness_data file.  call only at the appropriate point.
	 
	 \param IN a pointer to an open file from which to read.  must be set to the correct point
	 */
	void set_from_file(FILE *IN)
	{
		fscanf(IN,"%lf %d %lf %lf %d %d %d %d",
			   &condition_number_,
			   &corank_,
			   &smallest_nonzero_sing_value_,
			   &largest_zero_sing_value_,
			   &typeflag_, // 10 is nonsingular, 15 is singular
			   &multiplicity_,
			   &component_number_,
			   &deflations_needed_);
	}
	
	
	
	
	/**
	 output to a stream.  only really usable with std::cout.
	 \param os the stream to put this text on.
	 \param s the system randomizer to write.
	 */
	friend std::ostream & operator<<(std::ostream &os, WitnessPointMetadata & s)
	{
		os << "condition_number " << s.condition_number_ << "\n";
		os << "corank " << s.corank_ << "\n";
		os << "smallest_nonzero_sing_value " << s.smallest_nonzero_sing_value_ << "\n";
		os << "largest_zero_sing_value " << s.largest_zero_sing_value_ << "\n";
		os << "typeflag " << s.typeflag_ << "\n";
		os << "multiplicity " << s.multiplicity_ << "\n";
		os << "component_number " << s.component_number_ << "\n";
		os << "deflations_needed " << s.deflations_needed_ << std::endl;
		
		return os;
	}
	
	
	
	
	
	/**
	 \brief constructor, setting the dimension in the process.
	 constructor, setting the dimension in the process.
	 
	 \param new_dim the dimension to set to.
	 */
	WitnessPointMetadata(int new_dim){dimension_ = new_dim;}
	
	
	WitnessPointMetadata(){};
	
	WitnessPointMetadata(const WitnessPointMetadata& other){copy(other);}
	
	WitnessPointMetadata& operator=( const WitnessPointMetadata& other)
	{
		copy(other);
		return *this;
	}
	
	void copy(const WitnessPointMetadata& other)
	{
		dimension_ = other.dimension_;
		corank_ = other.corank_;
		typeflag_ = other.typeflag_;
		multiplicity_ = other.multiplicity_;
		component_number_ = other.component_number_;
		deflations_needed_ = other.deflations_needed_;
		condition_number_ = other.condition_number_;
		smallest_nonzero_sing_value_ = other.smallest_nonzero_sing_value_;
		largest_zero_sing_value_ = other.largest_zero_sing_value_;
	}
	
	
};


/**
 \brief metadata for linears read in from the witness_data file.
 
 metadata for linears read in from the witness_data file.
 */
class WitnessLinearMetadata
{
	int dim_;
	
public:
	
	
	/**
	 get the dimension
	 \return the dimension
	 */
	inline int dimension() const
	{
		return dim_;
	}
	
	
	WitnessLinearMetadata(int new_dim){dim_ = new_dim;}
	
	WitnessLinearMetadata(){};
	WitnessLinearMetadata(const WitnessLinearMetadata& other)
	{
		dim_ = other.dim_;
	}
	
	WitnessLinearMetadata& operator=(const WitnessLinearMetadata &other)
	{
		dim_ = other.dim_;
		return *this;
	}
};









/**
 \brief metadata for patches read in from witness_data
 
 metadata for patches read in from witness_data
 */
class WitnessPatchMetadata
{
	int dim_;
	
public:
	
	/**
	 get the dimension
	 \return the dimension
	 */
	inline int dimension() const
	{
		return dim_;
	}
	
	
	WitnessPatchMetadata(){};
	
	
	/**
	 constructor, setting the dimension in the process
	 \param new_dim the dimension to set.
	 */
	WitnessPatchMetadata(int new_dim){dim_ = new_dim;}

	
	WitnessPatchMetadata(const WitnessPatchMetadata& other)
	{
		dim_ = other.dim_;
	}
	
	WitnessPatchMetadata& operator=(const WitnessPatchMetadata &other)
	{
		dim_ = other.dim_;
		return *this;
	}
	
	
	
};




/**
 \brief a nearly complete class for storing bertini's witness_data file.
 
 This class reads in witness_data, and produces witness_sets based on user's choice.
 */
class NumericalIrreducibleDecomposition : public PatchHolder, public LinearHolder, public PointHolder
{
	
private:
	
	std::vector< WitnessPointMetadata > point_metadata;
	std::vector< WitnessLinearMetadata > linear_metadata;
	std::vector< WitnessPatchMetadata > patch_metadata;
	
	std::vector<int> nonempty_dimensions;
	std::map<int,std::map<int,int> > dimension_component_counter;
	std::map<int,std::map<int,std::vector<int> > > index_tracker;
	
	std::vector< std::vector< int> > homogenization_matrix_;
	int num_variables_;
	
public:
	
	
	/**
	 \brief get the number of variables
	 
	 \return the number of variables
	 */
	inline int num_variables() const
	{
		return num_variables_;
	}
	
	
	void reset()
	{
		reset_points();
		reset_linears();
		reset_patches();
		
		point_metadata.resize(0);
		linear_metadata.resize(0);
		patch_metadata.resize(0);
		
	}
	
	
	/** fills this object with the sets in witness_data.
	 \param T the current state of the tracker configuration
	 */
	void populate(tracker_config_t * T);
	
	
	
	/**
	 \brief outermost method for choosing a witness set to construct.
	 
	 This function uses information stored in BertiniRealConfig to construct a witness set.
	 
	 \param options The current state of the program.  If the user passed in a particular component or dimension, this is how it gets into this method.
	 \return the chosen witness set.  may be empty.
	 */
	WitnessSet choose(BertiniRealConfig & options);
	
	
	/**
	 form a witness set automagically, based on the user's call time options.
	 
	 \param options The current state of the program.  If the user passed in a particular component or dimension, this is how it gets into this method.
	 \return the best possible witness set, based on the user's choices at call time to Bertini_real.
	 */
	WitnessSet best_possible_automatic_set(BertiniRealConfig & options);
	
	/**
	 If there are multiple dimensions and components, then the user needs to choose which he wishes to decompose.  This method is that choice.
	 
	 \param options The current state of the program.  If the user passed in a particular component or dimension, this is how it gets into this method.
	 \return A formed witness set, from the user's choices.
	 */
	WitnessSet choose_set_interactive(BertiniRealConfig & options); // lets the user choose a set, and returns a copy of it.
	
	
	/**
	 retreive from the NID the specific witness set of particular dimension and component number.  If either doesn't exist, then the returned witness set is empty.
	 
	 \return The specific witness set desired.  May be empty if the dimension or component number correspond to something don't exist.
	 \param dim The desired dimension
	 \param comp The desired component number.  Starts at 0.
	 */
	WitnessSet form_specific_witness_set(int dim, int comp)	;
	
	
	
//	friend std::ostream & operator<<(std::ostream &os, witness_data & c)
//	{
//		for (auto iter=dimension_component_counter.begin(); iter!=dimension_component_counter.end(); ++iter) {
//			os << "dimension " << iter->first << std::endl;
//			for (auto jter=iter->second.begin(); jter!=iter->second.end(); ++jter) {
//				os << "\tcomponent " << jter->first << " degree " << jter->second << std::endl;
//			}
//		}
//		
//		
//		for (auto iter=index_tracker.begin(); iter!=index_tracker.end(); ++iter) {
//			os << "dimension " << iter->first << " ";
//			for (auto jter = iter->second.begin(); jter!=iter->second.end(); ++jter) {
//				os << "component " << jter->first << std::endl;
//				for (auto kter = jter->second.begin(); kter!=jter->second.end(); kter++) {
//					os << *kter << " ";
//				}
//				os << std::endl;
//			}
//			os << std::endl;
//		}
//		
//		
//		
//		
//		
//		return os;
//	}
	
	
	
	
	/**
	 print the witness_data to std::cout
	 
	 \todo replace this with a friend operator << ()
	 */
	void print() const
	{
		
		std::cout << "the nonempty dimensions:" << std::endl;
		for (auto iter=nonempty_dimensions.begin(); iter!=nonempty_dimensions.end(); ++iter) {
			std::cout << *iter << " ";
		}
		std::cout << std::endl;
		
		
		
		for (auto iter=dimension_component_counter.begin(); iter!=dimension_component_counter.end(); ++iter) {
			std::cout << "dimension " << iter->first << std::endl;
			for (auto jter=iter->second.begin(); jter!=iter->second.end(); ++jter) {
				std::cout << "\tcomponent " << jter->first << " degree " << jter->second << std::endl;
			}
		}
		
		
		for (auto iter=index_tracker.begin(); iter!=index_tracker.end(); ++iter) {
			std::cout << "dimension " << iter->first << std::endl;
			for (auto jter = iter->second.begin(); jter!=iter->second.end(); ++jter) {
				std::cout << "\tcomponent " << jter->first << ": ";
				for (auto kter = jter->second.begin(); kter!=jter->second.end(); kter++) {
					std::cout << *kter << " ";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
		
	}
	
	
	
	
	
	
	
	
	
	
	
private:
	
	void add_linear_w_meta(vec_mp lin, const WitnessLinearMetadata & meta)
	{
		add_linear(lin);
		linear_metadata.push_back(meta);
	}
	
	void add_patch_w_meta(vec_mp pat, const WitnessPatchMetadata & meta)
	{
		add_patch(pat);
		patch_metadata.push_back(meta);
	}
	
	
	int add_solution(vec_mp pt, const WitnessPointMetadata & meta)
	{
		
		int ind = add_point(pt);
		point_metadata.push_back(meta);
		return ind;
	}
	
	
	
	
	
};
















