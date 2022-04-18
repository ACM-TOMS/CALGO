// Copyright (c) 2012, 2013 David Munger, Pierre L'Ecuyer, Université de Montréal.
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

#ifndef LATBUILDER__VECTORIZE_H
#define LATBUILDER__VECTORIZE_H

#include <type_traits>
#include <boost/numeric/ublas/vector.hpp>

namespace LatBuilder {

/**
 * Helpers to vectorize simple operations.
 */
namespace Vectorize {

/**
 * Type traits class that checks if \c T is indexable with [].
 */
template <typename T>
class IsIndexable {
private:
   template <typename> static std::false_type helper(...);
   template <typename Tp> static std::true_type helper(typename std::remove_reference<decltype((*(Tp*)nullptr)[0])>::type*);
   typedef decltype(helper<T>(nullptr)) value_type;
public:
   constexpr static bool value = value_type::value;
};

/**
 * Guesses the result type of the binary operator \c OP given arguments of type
 * \c T1 and \c T2.
 */
template <typename OP, typename T1, typename T2>
struct BinaryOperatorResult {
   typedef typename std::conditional<
      IsIndexable<T1>::value or IsIndexable<T2>::value,
      boost::numeric::ublas::vector<typename OP::result_type>,
      typename OP::result_type
         >::type type;
};

/**
 * Automatic dispatcher for binary operators.
 *
 * \tparam OP        Binary operator.
 * \tparam ISVEC1    \c true if first operand is a vector.
 * \tparam ISVEC2    \c true if second operand is a vector.
 */
template <typename OP, bool ISVEC1, bool ISVEC2>
struct BinaryOperator {
   typedef typename OP::result_type result_type;

   template <typename T1, typename T2>
   static result_type apply(const T1& x1, const T2& x2)
   { return OP::apply(x1, x2); }
};

/**
 * Vector-vector specialization of BinaryOperator.
 */
template <typename OP>
struct BinaryOperator<OP, true, true> {
   typedef boost::numeric::ublas::vector<typename OP::result_type> result_type;

   template <typename E1, typename E2>
   static result_type apply(
         const boost::numeric::ublas::vector_expression<E1>& e1,
         const boost::numeric::ublas::vector_expression<E2>& e2)
   {
      const auto& x1 = e1();
      const auto& x2 = e2();
#ifndef NDEBUG
      if (x1.size() != x2.size())
         throw std::logic_error("Vectorize::BinaryOperator: arguments must have the same size");
#endif
      result_type z(x1.size());
      for (typename result_type::size_type i = 0; i < z.size(); i++)
         z[i] = OP::apply(x1[i], x2[i]);
      return z;
   }
};


/**
 * Vector-scalar specialization of BinaryOperator.
 */
template <typename OP>
struct BinaryOperator<OP, true, false> {
   typedef boost::numeric::ublas::vector<typename OP::result_type> result_type;

   template <typename E1, typename T2>
   static result_type apply(
         const boost::numeric::ublas::vector_expression<E1>& e1,
         const T2& x2)
   {
      const auto& x1 = e1();
      result_type z(x1.size());
      for (typename result_type::size_type i = 0; i < z.size(); i++)
         z[i] = OP::apply(x1[i], x2);
      return z;
   }
};


/**
 * Applies the scalar binary operator \c OP to \c x and \c y.
 *
 * Example definition of \c OP:
 * \code
 * struct Sum {
 *    typedef double result_type;
 *    static constexpr result_type apply(const arg_type& x, const arg_type& y)
 *    { return x + y; }
 * };
 * \endcode
 */
template <typename OP, typename T1, typename T2>
typename BinaryOperatorResult<OP, T1, T2>::type
apply(const T1& x, const T2& y)
{ return BinaryOperator<OP, IsIndexable<T1>::value, IsIndexable<T2>::value>::apply(x, y); }

}}

#endif
