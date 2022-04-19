/**
 * \file mpi_tree.hpp
 * \author Dhairya Malhotra, dhairya.malhotra@gmail.com
 * \date 12-11-2010
 * \brief This file contains the definition of a base class for a distributed
 * MPI tree.
 */

#include <mpi.h>
#include <vector>
#include <string>

#include <pvfmm_common.hpp>
#include <mortonid.hpp>
#include <tree.hpp>

#ifndef _PVFMM_MPI_TREE_HPP_
#define _PVFMM_MPI_TREE_HPP_

namespace pvfmm{

enum BoundaryType{
  FreeSpace,
  Periodic
};

/**
 * \brief Base class for distributed tree.
 */
template <class TreeNode>
class MPI_Tree: public Tree<TreeNode>{

 public:

  typedef TreeNode Node_t;
  typedef typename Node_t::Real_t Real_t;

  /**
   * \brief Constructor.
   */
  MPI_Tree(MPI_Comm c): Tree<Node_t>() {comm=c;}

  /**
   * \brief Virtual destructor.
   */
  virtual ~MPI_Tree() {}

  /**
   * \brief Initialize the distributed MPI tree.
   */
  virtual void Initialize(typename Node_t::NodeData* data_);

  /**
   * \brief Initialize the distributed MPI tree.
   */
  //void InitData(typename Node_t::NodeData& d1);

  /**
   * \brief Find the prticular node. If subdiv is true then subdivide
   * (non-ghost) nodes to create this node.
   */
  TreeNode* FindNode(MortonId& key, bool subdiv, TreeNode* start=NULL);

  /**
   * \brief Adaptive coarsening of distributed tree.
   */
  virtual void CoarsenTree();

  /**
   * \brief Adaptive refinement of distributed tree.
   */
  virtual void RefineTree();

  /**
   * \brief Redistribute the tree among the processes.
   */
  void RedistNodes(MortonId* loc_min=NULL);

  /**
   * \brief Performs global 2:1 balancing of the tree.
   */
  void Balance21(BoundaryType bndry=FreeSpace);

  /**
   * \brief Performs 2:1 balancing of the tree.
   */
  void Balance21_local(BoundaryType bndry=FreeSpace);

  /**
   * \brief Determines and sets colleagues for each node in the tree.
   * Two nodes are colleagues if they are at the same depth in the tree and
   * share either a face, edge or a vertex.
   */
  void SetColleagues(BoundaryType bndry=FreeSpace, Node_t* node=NULL) ;

  /**
   * \brief Checks if everything is okay with the distributed tree. Only needed
   * for debugging purposes.
   */
  bool CheckTree();

  /**
   * \brief Construct the LET by exchanging ghost octants.
   */
  void ConstructLET(BoundaryType bndry=FreeSpace);

  template <class VTKReal>
  struct VTUData_t{
    typedef VTKReal VTKReal_t;

    // Point data
    std::vector<VTKReal_t> coord;
    std::vector<std::vector<VTKReal_t> > value;
    std::vector<std::string> name;

    // Cell data
    std::vector<int32_t> connect;
    std::vector<int32_t> offset;
    std::vector<uint8_t> types;
  };

  /**
   * \brief Write to a <fname>.vtu file.
   */
  void Write2File(const char* fname, int lod=-1);

  /**
   * \brief Returns a pointer to the comm object.
   */
  const MPI_Comm* Comm() {return &comm;}

 protected:

  /**
   * \brief Returns a vector with the minimum Morton Id of the regions
   * controlled by each processor.
   */
  const std::vector<MortonId>& GetMins();

 private:

  void ConstructLET_Hypercube(BoundaryType bndry=FreeSpace);
  void ConstructLET_Sparse(BoundaryType bndry=FreeSpace);

  MPI_Comm comm;
  std::vector<MortonId> mins;

};

}//end namespace

#include <mpi_tree.txx>

#endif //_PVFMM_MPI_TREE_HPP_
