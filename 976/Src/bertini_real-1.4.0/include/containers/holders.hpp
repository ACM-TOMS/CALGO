#pragma once



#include "bertini1/bertini_extensions.hpp"
#include "io/fileops.hpp"




/**
 \brief base class for holding a set of vec_mp's.
 
 This class gives a way to commit vec_mp's into a class object.  The major data members are pts_mp and num_pts.  The most important member function is add_point(p)
 */
class PointHolder
{
	

protected:
	
	vec_mp *pts_mp_; ///< an array of vec_mp, which are structs and require manual initialization and clearing.
	size_t num_pts_; ///< the number of stored points.
	
public:
	
	
	/**
	 get the ith point
	 
	 \return a reference to the point at the requested index.
	 \param index the index of the point to retrieve.
	 \throws out_of_range if the requested point does not exist.
	 */
	vec_mp & point(unsigned int index) const
	{
		if (index < num_pts_) {
			return pts_mp_[index];
		}
		else
		{
			throw std::out_of_range("trying to get a point in PointHolder which is out of range");
		}
		
	}
	
	
	/**
	 get the number of points in the point holder
	 
	 \return the current number of points
	 */
	size_t num_points() const
	{
		return num_pts_;
	}
	
	
	/**
	 copy all the points from one PointHolder to this one.
	 
	 \param other an input PointHolder from which to copy all the points.
	 */
	void copy_points(const PointHolder & other)
	{
		
        for (unsigned int ii=0; ii<other.num_pts_; ii++)
			add_point(other.pts_mp_[ii]);
    }
	
	
	/**
	 indicate whether the PointHolder has at zero points in it.
	 
	 \return a boolean, true if num_points==0, false otherwise
	 */
	inline bool has_no_points() const
	{
		if (num_pts_==0) {
			return true;
		}
		else{
			return false;
		}
	}
	
	
	/**
	 indicate whether the PointHolder has at least one point in it.
	 
	 \return a boolean, true if num_points>0, false if num_points==0
	 */
	inline bool has_points() const
	{
		if (num_pts_==0) {
			return false;
		}
		else{
			return true;
		}
	}
	
	
	/**
	 \brief reset to empty container.
	 
	 reset to empty container
	 */
	void reset_points()
	{
		for (unsigned int ii =0; ii<num_pts_; ii++)
			clear_vec_mp(pts_mp_[ii]);
		
		if (num_pts_>0) {
			free(pts_mp_);
		}
		
		num_pts_ = 0;
		pts_mp_ = NULL;
	}
	
	void real_threshold_points(double threshold)
	{
		for (unsigned int ii = 0; ii<num_pts_; ii++)
			real_threshold(pts_mp_[ii], threshold);
	}
	
	/**
	 \brief add a point to the PointHolder
	 
	 \return the index of the new point.
	 \param new_point the point to add.
	 */
	int add_point(vec_mp new_point);
	
	
	/**
	 constructor
	 */
	PointHolder(){
		init();
	}
	
	
	/**
	 destructor
	 */
	~PointHolder(){ // the
		
		clear();
		
	}
	
	/**
	 assignment
	 
	 \param other the other to copy from.
	 */
	PointHolder& operator=( const PointHolder & other) {
		
		reset_points();
		
		copy(other);
		return *this;
	}
	
	/**
	 copy operator.  must be explicitly declared because the underlying c structures use pointers.
	 
	 \param other the other to copy from
	 */
	PointHolder(const PointHolder & other){
		init();
		copy(other);
	}
	
	
	/**
	 copy from one to another.
	 
	 \param other the other to copy from
	 */
	void copy(const PointHolder & other)
	{
		copy_points(other);
	}
	
	
	
private:
	
	
	/**
	 set up the pointers to good value.
	 */
	void init()
	{
		pts_mp_ = NULL; // potential data loss if used improperly.
		num_pts_ = 0;
	}
	
	
	/**
	 purge the held data.
	 */
	void clear(){
		reset_points();
	}
};




/**
 \brief a way to hold vec_mp's as patches.
 
 This class offers a way to hold a set of vec_mp's as patch_mp members in a class with automated initting and clearing.
 */
class PatchHolder
{

protected:
	
	vec_mp *patch_mp_;   ///< an array of patch_mp's
	size_t num_patches_; ///< the number of patches stored in this object.
	
public:
	
	
	/**
	 \brief  query how many patches are being held.
	 query how many patches are being held.
	 
	 \return the number of patches.
	 */
	inline size_t num_patches() const
	{
		return num_patches_;
	}
	
	
	/**
	\brief get a pointer to the patch at index
	 
	 \return a pointer to the i^th patch.
	 \param index The index of the patch you want
	 */
	inline vec_mp & patch(unsigned int index) const
	{
		if (index >= num_patches_) {
			throw std::out_of_range("trying to access an out-of-range patch");
		}
		else
			return patch_mp_[index];
	}
	
	
	/**
	 \brief copy all the patches from another PatchHolder
	 
	 Copy all the stored mp patches from another PatchHolder object, without testing for uniqueness.
	 
	 \param other the PatchHolder from which to copy.
	 */
	void copy_patches(const PatchHolder & other) {
		
        for (unsigned int ii=0; ii<other.num_patches_; ii++)
			add_patch(other.patch_mp_[ii]);
    }
	
	
	
	/**
	 \brief resets the patch holder to empty.
	 
	 Reset the PatchHolder to 0 patches.
	 */
	void reset_patches()
	{
		for (unsigned int ii =0; ii<num_patches_; ii++)
			clear_vec_mp(patch_mp_[ii]);
		
		if (num_patches_>0) {
			free(patch_mp_);
		}
		
		num_patches_ = 0;
		patch_mp_ = NULL;
	}
	
	
	/**
	 \brief add a patch to this object.
	 
	 \return the index of the patch just added.
	 \param new_patch the new patch to add.
	 */
	int add_patch(vec_mp new_patch);
	
	
	/**
	 constructor
	 */
	PatchHolder(){
		init();
	}
	
	
	
	/**
	 destructor
	 */
	~PatchHolder(){ // the destructor
		
		clear();
		
	}
	
	
	/**
	 assignment
	 */
	PatchHolder& operator=( const PatchHolder& other) {
		
		reset_patches();
		
		copy(other);
		return *this;
	}
	
	/**
	 copy operator.  must be explicitly declared because the underlying c structures use pointers.
	*/
	 PatchHolder(const PatchHolder & other){
		init();
		copy(other);
	}
	
	
	/**
	 copy from another
	 
	 \param other the other to copy from.
	 */
	void copy(const PatchHolder & other)
	{
		copy_patches(other);
	}
	
	
	
private:
	
	
	/**
	 initialize
	 */
	void init()
	{
		this->patch_mp_ = NULL;
		this->num_patches_ = 0;
	}
	
	
	/**
	 purge old patches
	 */
	void clear(){
		reset_patches();
	}
};


/**
 \brief class for holding vec_mp's as linears in an automated object.
 
 This class offers automated collection of vec_mp's as L_mp.  This is necessary because the vec_mp type must be initted and cleared manually.
 */
class LinearHolder
{
	
protected:
	
	vec_mp *L_mp_;  ///< a pointer array of vec_mp's as linears.
	size_t num_linears_; ///< the number of linears in this collection.
	
	
	
public:
	
	/**
	 \brief  query how many patches are being held.
	 query how many patches are being held.
	 
	 \return the number of patches.
	 */
	inline size_t num_linears() const
	{
		return num_linears_;
	}
	
	
	/**
	 \brief get a pointer to the patch at index
	 
	 \return a pointer to the i^th patch.
	 \param index The index of the patch you want
	 */
	inline vec_mp & linear(unsigned int index) const
	{
		if (index>= num_linears_) {
			throw std::out_of_range("trying to access out-of-range linear");
		}
		else
			return L_mp_[index];
	}
	
	
	/**
	 \brief copies all linears from another LinearHolder.
	 
	 \param other the LinearHolder from which to copy all the linears.
	 */
	void copy_linears(const LinearHolder & other) {
        for (unsigned int ii=0; ii<other.num_linears_; ii++)
			add_linear(other.L_mp_[ii]);
    }
	
	
	/**
	 reset this object to an empty state.
	 */
	void reset_linears()
	{
		for (unsigned int ii =0; ii<num_linears_; ii++)
			clear_vec_mp(L_mp_[ii]);
		
		if (num_linears_>0) {
			free(L_mp_);
		}
		
		num_linears_ = 0;
		L_mp_ = NULL;
	}
	
	
	/**
	 \brief add a linear to this collection.
	 
	 \return the index of the newly added linear.
	 \param new_linear the vec_mp linear to add.
	 */
	int add_linear(vec_mp new_linear);
	
	
	
	LinearHolder(){
		init();
	}
	
	
	~LinearHolder(){ // the destructor
		
		clear();
		
	}
	
	
	// assignment
	LinearHolder& operator=( const LinearHolder& other) {
		reset_linears();
		copy(other);
		return *this;
	}
	
	
	//copy operator.  must be explicitly declared because the underlying c structures use pointers.
	LinearHolder(const LinearHolder & other){
		init();
		copy(other);
	}
	
	
	
	void copy(const LinearHolder & other)
	{
		copy_linears(other);
	}
	
	
	
private:
	
	
	
	void init()
	{
		this->L_mp_ = NULL;
		this->num_linears_ = 0;
	}
	
	void clear(){
		reset_linears();
	}
};


/** 
 \brief holds the names of variables
 */
class NameHolder
{

	
private:
	
	std::vector< std::string > variable_names; ///< the names.
	
public:
	
	
	/**
	 \brief get how many names there are.
	 
	 \return the number of stored names
	 */
	inline size_t num_var_names() const
	{
		return variable_names.size();
	}
	
	
	/**
	 \brief get the ith variable name
	 \param index The index of the name desired.
	 \return the ith name.
	 */
	inline std::string name(unsigned int index) const
	{
		if (index>=variable_names.size()) {
			throw std::out_of_range("trying to get variable name, requested index is out of bounds");
		}
		else
		{
			return variable_names[index];
		}
		
	}
	
	
	
	/**
	 \brief reset the names to empty.
	 */
	void reset_names()
	{
		variable_names.resize(0);
	}
	
	void copy_names(const NameHolder & nomnom)
	{
		this->variable_names = nomnom.variable_names;
	}
	
	
	/**
	 \brief read variable names from names.out
	 
	 Reads variable names from names.out, which must exist before calling this function.
	 
	 \param num_vars the number of variable names to read.
	 */
	void get_variable_names(int num_vars)
	{
		
		variable_names.resize(num_vars);
		
		std::ifstream fin("names.out");
		
		for (int ii=0; ii<num_vars; ++ii){
			fin >> this->variable_names[ii];
		}
		fin.close();
		
	}
	
	
};









