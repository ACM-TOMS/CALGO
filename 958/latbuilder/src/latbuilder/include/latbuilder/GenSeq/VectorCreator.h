// Copyright (c) 2012 David Munger, Pierre L'Ecuyer, Université de Montréal.
// 
// This file is part of Lattice Builder.
// 
// Lattice Builder is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Lattice Builder is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Lattice Builder.  If not, see <http://www.gnu.org/licenses/>.

#ifndef LATBUILDER__GEN_SEQ__VECTOR_CREATOR_H
#define LATBUILDER__GEN_SEQ__VECTOR_CREATOR_H

#include "latbuilder/GenSeq/Creator.h"
#include "latbuilder/SizeParam.h"
#include "latbuilder/Traversal.h"

#include <vector>

namespace LatBuilder { namespace GenSeq {

namespace detail {
   template <class TRAV>
   class Traversal {
   public:
      template <typename... ARGS>
      Traversal(ARGS&&... args):
         m_trav(std::forward<ARGS>(args)...)
      {}

      TRAV operator()()
      { return m_trav; }

   private:
      TRAV m_trav;
   };

   template <typename RAND>
   class Traversal<LatBuilder::Traversal::Random<RAND>> {
   public:
      Traversal(
            size_t nrand = std::numeric_limits<size_t>::max(),
            RAND rand = RAND()
            ):
         m_nrand(nrand),
         m_rand(std::move(rand))
      {}

      LatBuilder::Traversal::Random<RAND> operator()()
      {
         LatBuilder::Traversal::Random<RAND> trav(m_nrand, m_rand);
         m_rand.jump();
         return trav;
      }

   private:
      size_t m_nrand;
      RAND m_rand;
      LatBuilder::Traversal::Random<RAND> m_trav;
   };
}

/**
 * Creator for vectors of sequences of generator values.
 */
template <typename SEQ>
struct VectorCreator {
   typedef std::vector<SEQ> result_type;

   /**
    * Creates a new sequence object.
    *
    * \param sizeParam        Size parameter for the generator sequences.
    * \param dimension        Dimension of the output vector.
    * \param t                Other arguments to be passed to GenSeq::Creator
    */
   template <LatType L, typename... ARGS>
   static result_type create(
         const SizeParam<L>& sizeParam,
         Dimension dimension, 
         ARGS&&... t
         )
   {
      result_type vec(dimension);
      detail::Traversal<typename SEQ::Traversal> trav(std::forward<ARGS>(t)...);
      for (Dimension j = 0; j < dimension; j++)
         vec[j] = GenSeq::Creator<SEQ>::create(sizeParam, trav());
      return vec;
   }
};


}}

#endif
