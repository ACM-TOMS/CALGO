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

#ifndef LATBUILDER__MERIT_SEQ__INNER_PROD_FAST_H
#define LATBUILDER__MERIT_SEQ__INNER_PROD_FAST_H

#include "latbuilder/MeritSeq/CoordUniformStateCreator.h"
#include "latbuilder/BridgeSeq.h"
#include "latbuilder/BridgeIteratorCached.h"
#include "latbuilder/Storage.h"
#include "latbuilder/CachedSeq.h"
#include "latbuilder/IndexMap.h"
#include "latbuilder/fftw++.h"

#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <memory>
#include <vector>

namespace LatBuilder { namespace MeritSeq {

namespace detail {
   template <class GENSEQ>
   struct IsFastCompatible
   { static constexpr bool value = false; };

   template <Compress COMPRESS>
   struct IsFastCompatible<GenSeq::CyclicGroup<COMPRESS>> 
   { static constexpr bool value = true; };
}

/**
 * FFT-based implementation of the inner product for a circulant sequence of
 * vector with a single vector.
 *
 * Implemented for integer powers of prime bases, as proposed in \cite rCOO06a .
 *
 * Computes the inner product with a second vector for all vectors in the
 * sequence at once.
 */
template <LatType LAT, Compress COMPRESS>
class CoordUniformInnerProdFast {

protected:
   typedef typename fftw<Real>::real_vector FFTRealVector;
   typedef typename fftw<Real>::complex_vector FFTComplexVector;

public:
   typedef Storage<LatType::EMBEDDED, COMPRESS> InternalStorage;
   typedef CoordUniformStateList<LatType::EMBEDDED, COMPRESS> StateList;
   typedef typename Storage<LAT, COMPRESS>::MeritValue MeritValue;

   /**
    * Constructor.
    *
    * \param storage       Storage configuration.
    * \param kernel        Kernel.  Used to create a sequence of
    *                      permuatations of the kernel values evaluated at every
    *                      one-dimensional lattice point.
    */
   template <class K>
   CoordUniformInnerProdFast(
         Storage<LAT, COMPRESS> storage,
         const Kernel::Base<K>& kernel
         ):
      m_storage(std::move(storage)),
      m_internalStorage(asIntenalStorage(this->storage())),
      m_kernelValues(kernel.valuesVector(this->internalStorage())),
      m_levelRanges(cacheLevelRanges()),
      m_circulantFFT(computeCirculantFFT())
   {}

   /**
    * Returns the storage configuration instance.
    */
   const Storage<LAT, COMPRESS>& storage() const
   { return m_storage; }

   /**
    * Returns the internal storage configuration instance.
    */
   const InternalStorage& internalStorage() const
   { return m_internalStorage; }

   /**
    * Returns the vector of kernel values.
    */
   const RealVector& kernelValues() const
   { return m_kernelValues; }

   /**
    * Returns the vector of per-level ranges of indices.
    */
   const std::vector<boost::numeric::ublas::range>& levelRanges() const
   { return m_levelRanges; }

   /**
    * Returns the FFT's of the first column of each circulant submatrix in the
    * horizontal block-circulant matrix.
    */
   const std::vector<FFTComplexVector>& circulantFFT() const
   { return m_circulantFFT; }


private:
   std::vector<boost::numeric::ublas::range> cacheLevelRanges() const
   {
      const auto ranges = m_internalStorage.levelRanges();
      std::vector<boost::numeric::ublas::range> out;
      out.assign(ranges.begin(), ranges.end());
      return out;
   }

   template <class E>
   RealVector computeProdValues(
         const boost::numeric::ublas::vector_expression<E>& ve
         ) const
   {
      // in base 2, on each level, we have 2 circulant half-blocks instead of
      // a single circulant block
      if (not internalStorage().symmetric() and internalStorage().sizeParam().base() == 2)
         throw std::logic_error("not implemented for non-symmetric vectors in base 2");

      const auto& vec = ve();
      using namespace boost::numeric::ublas;

      auto itRange = levelRanges().begin();
      auto itRangePrev = itRange;
      auto itCirculantFFT = circulantFFT().begin();

      RealVector out(vec.size());

      while (itRange != levelRanges().end()) {

         if (itCirculantFFT == circulantFFT().end())
            throw std::logic_error("circulant FFT's have too few levels");

         // select vector range
         vector_range<const RealVector> subvec(vec, *itRange);

         // convert to FFT-compatible vectors
         FFTRealVector rvec(subvec.begin(), subvec.end());

         // compute FFT
         FFTComplexVector cvec = fftw<Real>::fft(rvec);

         // ratio of the number or natural elements to the number of internal
         // elements, multiplied by normalization
         size_t compressionRatio = 1;
         if (internalStorage().symmetric() and (itRange - levelRanges().begin()) >= (internalStorage().sizeParam().base() == 2 ? 2 : 1)) {
            // compressionRatio except if uncompressed level has only one element
            compressionRatio = 2;
         }

         // multiply in Fourier space
         for (size_t i = 0; i < cvec.size(); i++)
            cvec[i] *= compressionRatio * (*itCirculantFFT)[i];

         // inverse transform
         fftw<Real>::ifft(cvec, rvec, true);

         // export to the output vector
         std::copy(rvec.begin(), rvec.end(), &out[itRange->start()]);

         // add contributions from lower levels
         if (itRange != itRangePrev) {
            typedef typename vector_range<RealVector>::size_type size_type;
            vector_range<RealVector> curLevel(out, *itRange);
            vector_range<const RealVector> prevLevel(out, *itRangePrev);
            for (size_type i = 0; i < curLevel.size(); i++)
               curLevel[i] += prevLevel[i % prevLevel.size()];
         }

         itRangePrev = itRange;

         ++itCirculantFFT;
         ++itRange;
      }

      return out;
   }

public:
   /**
    * Sequence of inner product values.
    *
    * \tparam GENSEQ    Type of sequence of generator values.  Must be a binding of
    *                   GenSeq::CyclicGroup.
    */
   template <class GENSEQ>
   class Seq :
      public BridgeSeq<
         Seq<GENSEQ>,                           // self type
         GENSEQ,                                // base type
         MeritValue,                            // value type
         BridgeIteratorCached> {

      static_assert(
            detail::IsFastCompatible<GENSEQ>::value,
            "generator sequence is not compatible with fast product");
   public:

      typedef GENSEQ GenSeq;
      typedef typename Seq::Base Base;
      typedef typename Seq::size_type size_type;

      /**
       * Constructor.
       *
       * \param parent     Parent inner product instance.
       * \param genSeq     Sequence of generator sequences that determines the
       *                   order of the permutations of \c baseVec.
       * \param vec        Second operand in the inner product.
       */
      template <class E>
      Seq(
            const CoordUniformInnerProdFast& parent, 
            GenSeq genSeq,
            const boost::numeric::ublas::vector_expression<E>& vec
            ):
         Seq::BridgeSeq_(std::move(genSeq)),
         m_parent(parent),
         m_values(parent.computeProdValues(vec()))
      {}

      /**
       * Returns the parent inner product of this sequence.
       */
      const CoordUniformInnerProdFast& innerProd() const
      { return m_parent; }

      MeritValue element(const typename Base::const_iterator& it) const
      {
         RealVector mlMerit(m_parent.levelRanges().size());
         for (
               auto itRange = m_parent.levelRanges().begin();
               itRange != m_parent.levelRanges().end();
               ++itRange
               ) {
            boost::numeric::ublas::vector_range<const RealVector> curLevel(m_values, *itRange);
            mlMerit[itRange - m_parent.levelRanges().begin()] = curLevel[(it - it.seq().begin()) % curLevel.size()];
         }

         MeritValue merit = m_parent.storage().createMeritValue(0.0);
         return storeMeritValue(merit, std::move(mlMerit));
      }

   private:
      const CoordUniformInnerProdFast& m_parent;
      RealVector m_values;

      Real& storeMeritValue(Real& dest, const RealVector& src) const
      { return dest = src(src.size() - 1); }

      RealVector& storeMeritValue(RealVector& dest, RealVector src) const
      { return dest = std::move(src); }
   };

   /**
    * Creates a new sequence of inner product values by applying a stride
    * permutation based on \c genSeq to the vector of kernel values, then by
    * computing the inner product with \c vec.
    *
    * \param genSeq     Sequence of generator values.
    * \param vec        Second operand in the inner product.
    */
   template <class GENSEQ, class E>
   Seq<GENSEQ> prodSeq(
         const GENSEQ& genSeq,
         const boost::numeric::ublas::vector_expression<E>& vec
         ) const
   { return Seq<GENSEQ>(*this, genSeq, vec); }


private:

   /**
    * Computes the FFT's of the first column of each circulant submatrix in the
    * horizontal block-circulant matrix.
    */
   std::vector<FFTComplexVector> computeCirculantFFT() const
   {
      const auto ranges = levelRanges();

      std::vector<FFTComplexVector> result(ranges.size());

      for (
            auto itRange = ranges.begin();
            itRange != ranges.end();
            ++itRange
            ) {

         // select level
         boost::numeric::ublas::vector_range<const RealVector> lvec(
               kernelValues(),
               *itRange
               );

         // apply circulant-transpose
         const auto tvec = circulantTranspose(lvec);

         // convert to FFT-compatible vectors
         FFTRealVector rvec(tvec.begin(), tvec.begin() + tvec.size());

         // compute FFT
         result[itRange - ranges.begin()] =
            fftw<Real>::fft(rvec);
      }

      return result;
   }

   InternalStorage asIntenalStorage(const InternalStorage& s)
   { return s; }

   template <LatType L, Compress C>
   InternalStorage asIntenalStorage(const Storage<L, C>& s)
   { return InternalStorage(s.sizeParam().numPoints()); }

private:
   class CirculantTranspose {
   public:
      typedef size_t size_type;
      CirculantTranspose(size_type size): m_size(size) {}
      size_type size() const { return m_size; }
      size_type operator() (size_type i) const { return i == 0 ? 0 : size() - i; }
   private:
      size_type m_size;
   };

   template <class E>
   static
   boost::numeric::ublas::vector_indirect<const E, IndexMap<CirculantTranspose>>
   circulantTranspose(const boost::numeric::ublas::vector_expression<E>& v)
   {
      typedef IndexMap<CirculantTranspose> Map;
      return boost::numeric::ublas::vector_indirect<const E, Map>(
            v(),
            Map(CirculantTranspose(v().size()))
            );
   }


private:
   Storage<LAT, COMPRESS> m_storage;
   InternalStorage m_internalStorage;
   RealVector m_kernelValues;
   std::vector<boost::numeric::ublas::range> m_levelRanges;
   std::vector<FFTComplexVector> m_circulantFFT;
};


}}

#endif
