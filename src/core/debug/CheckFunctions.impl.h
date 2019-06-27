//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file CheckFunctions.impl.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================



#include "core/debug/OperatorCheck.h"

#include <algorithm>
#include <type_traits>

#include <type_traits>

/// \cond internal

namespace walberla {
namespace debug {
namespace check_functions_detail {

template< typename T >
inline bool check_nullptr( T * p )
{
   return p == 0;
}

template< typename T >
inline bool check_nullptr( const shared_ptr<T> & p )
{
   return !p;
}

template< typename T >
inline bool check_not_nullptr( T * p )
{
   return p != 0;
}

template< typename T >
inline bool check_not_nullptr( const shared_ptr<T> & p )
{
   return !!p;
}

template< typename T, typename U >
inline bool check_equal( const T & lhs, const U & rhs )
{
   typedef std::integral_constant<bool, std::is_arithmetic<T>::value && std::is_arithmetic<U>::value> truth_type;
   return check_equal( lhs, rhs, truth_type() );
}

template< typename T, typename U >
inline bool check_equal( const T & lhs, const U & rhs, const std::true_type & )
{
   typedef typename math::MathTrait<T,U>::High HighType;
   return static_cast<HighType>(lhs) == static_cast<HighType>(rhs);
}

template< typename T, typename U >
inline bool check_equal( const T & lhs, const U & rhs, const std::false_type & )
{
   return lhs == rhs;
}

template< typename T, typename U >
inline bool check_unequal( const T & lhs, const U & rhs )
{
   typedef std::integral_constant<bool, std::is_arithmetic<T>::value && std::is_arithmetic<U>::value> truth_type;
   return check_unequal( lhs, rhs, truth_type() );
}

template< typename T, typename U >
inline bool check_unequal( const T & lhs, const U & rhs, const std::true_type & )
{
   typedef typename math::MathTrait<T,U>::High HighType;
   return static_cast<HighType>(lhs) != static_cast<HighType>(rhs);
}

template< typename T, typename U >
inline bool check_unequal( const T & lhs, const U & rhs, const std::false_type & )
{
   return lhs != rhs;
}

template< typename T, typename U >
inline bool check_float_equal( const T & lhs, const U & rhs )
{
   static_assert( std::is_floating_point<T>::value,  "First operand type T is not a floating point type!");
   static_assert( std::is_floating_point<U>::value, "Second operand type U is not a floating point type!");

   typedef typename math::MathTrait<T,U>::Low LowType;

   LowType low_lhs = static_cast<LowType>( lhs );
   LowType low_rhs = static_cast<LowType>( rhs );

   return floatIsEqual( low_lhs, low_rhs ) || floatIsEqual( ( low_lhs - low_rhs ) / low_lhs, LowType(0) );
}

template< typename T, typename U >
inline bool check_float_unequal( const T & lhs, const U & rhs )
{
   static_assert( std::is_floating_point<T>::value,  "First operand type T is not a floating point type!");
   static_assert( std::is_floating_point<U>::value, "Second operand type U is not a floating point type!");

   typedef typename math::MathTrait<T,U>::Low LowType;

   LowType low_lhs = static_cast<LowType>( lhs );
   LowType low_rhs = static_cast<LowType>( rhs );

   return !floatIsEqual( low_lhs, low_rhs ) && !floatIsEqual( ( low_lhs - low_rhs ) / low_lhs, LowType(0) );
}

template< typename T, typename U >
inline bool check_float_equal_eps( const T & lhs, const U & rhs,
                                   const typename VectorTrait<typename math::MathTrait<T,U>::LowType>::OutputType epsilon )
{
   static_assert( std::is_floating_point<T>::value,  "First operand type T is not a floating point type!");
   static_assert( std::is_floating_point<U>::value, "Second operand type U is not a floating point type!");

   typedef typename math::MathTrait<T,U>::Low LowType;

   LowType low_lhs = static_cast<LowType>( lhs );
   LowType low_rhs = static_cast<LowType>( rhs );

   return floatIsEqual( low_lhs, low_rhs, epsilon ) || floatIsEqual( ( low_lhs - low_rhs ) / low_lhs, LowType(0), epsilon );
}

template< typename T, typename U >
inline bool check_float_unequal_eps( const T & lhs, const U & rhs,
                                     const typename VectorTrait<typename math::MathTrait<T,U>::LowType>::OutputType epsilon )
{
   static_assert( std::is_floating_point<T>::value,  "First operand type T is not a floating point type!");
   static_assert( std::is_floating_point<U>::value, "Second operand type U is not a floating point type!");

   typedef typename math::MathTrait<T,U>::Low LowType;

   LowType low_lhs = static_cast<LowType>( lhs );
   LowType low_rhs = static_cast<LowType>( rhs );

   return !floatIsEqual( low_lhs, low_rhs, epsilon ) && !floatIsEqual( ( low_lhs - low_rhs ) / low_lhs, LowType(0), epsilon );
}

template< typename T, typename U >
inline bool check_identical( const T & lhs, const U & rhs )
{
   typedef std::integral_constant<bool, std::is_arithmetic<T>::value && std::is_arithmetic<U>::value> truth_type;
   return check_identical( lhs, rhs, truth_type() );
}

template< typename T, typename U >
inline bool check_identical( const T & lhs, const U & rhs, const std::true_type & )
{
   typedef typename math::MathTrait<T,U>::High HighType;
   return isIdentical( static_cast<HighType>(lhs), static_cast<HighType>(rhs) );
}

template< typename T, typename U >
inline bool check_identical( const T & lhs, const U & rhs, const std::false_type & )
{
   return lhs == rhs;
}

template< typename T, typename U >
inline bool check_not_identical( const T & lhs, const U & rhs )
{
   typedef std::integral_constant<bool, std::is_arithmetic<T>::value && std::is_arithmetic<U>::value> truth_type;
   return check_not_identical( lhs, rhs, truth_type() );
}

template< typename T, typename U >
inline bool check_not_identical( const T & lhs, const U & rhs, const std::true_type & )
{
   typedef typename math::MathTrait<T,U>::High HighType;
   return !isIdentical( static_cast<HighType>(lhs), static_cast<HighType>(rhs) );
}

template< typename T, typename U >
inline bool check_not_identical( const T & lhs, const U & rhs, const std::false_type & )
{
   return lhs != rhs;
}

template< typename T, typename U >
inline bool check_less( const T & lhs, const U & rhs )
{
   typedef std::integral_constant<bool, std::is_arithmetic<T>::value && std::is_arithmetic<U>::value> truth_type;
   return check_less( lhs, rhs, truth_type() );
}

template< typename T, typename U >
inline bool check_less( const T & lhs, const U & rhs, const std::true_type & )
{
   typedef typename math::MathTrait<T,U>::High HighType;
   return static_cast<HighType>(lhs) < static_cast<HighType>(rhs);
}

template< typename T, typename U >
inline bool check_less( const T & lhs, const U & rhs, const std::false_type & )
{
   return lhs < rhs;
}

template< typename T, typename U >
inline bool check_greater( const T & lhs, const U & rhs )
{
   typedef std::integral_constant<bool, std::is_arithmetic<T>::value && std::is_arithmetic<U>::value> truth_type;
   return check_greater( lhs, rhs, truth_type() );
}

template< typename T, typename U >
inline bool check_greater( const T & lhs, const U & rhs, const std::true_type & )
{
   typedef typename math::MathTrait<T,U>::High HighType;
   return static_cast<HighType>(lhs) > static_cast<HighType>(rhs);
}

template< typename T, typename U >
inline bool check_greater( const T & lhs, const U & rhs, const std::false_type & )
{
   return lhs > rhs;
}

template< typename T, typename U >
inline bool check_less_equal( const T & lhs, const U & rhs )
{
   typedef std::integral_constant<bool, std::is_arithmetic<T>::value && std::is_arithmetic<U>::value> truth_type;
   return check_less_equal( lhs, rhs, truth_type() );
}

template< typename T, typename U >
inline bool check_less_equal( const T & lhs, const U & rhs, const std::true_type & )
{
   typedef typename math::MathTrait<T,U>::High HighType;
   return static_cast<HighType>(lhs) <= static_cast<HighType>(rhs);
}

template< typename T, typename U >
inline bool check_less_equal( const T & lhs, const U & rhs, const std::false_type & )
{
   return lhs <= rhs;
}

template< typename T, typename U >
inline bool check_greater_equal( const T & lhs, const U & rhs )
{
   typedef std::integral_constant<bool, std::is_arithmetic<T>::value && std::is_arithmetic<U>::value> truth_type;
   return check_greater_equal( lhs, rhs, truth_type() );
}

template< typename T, typename U >
inline bool check_greater_equal( const T & lhs, const U & rhs, const std::true_type & )
{
   typedef typename math::MathTrait<T,U>::High HighType;
   return static_cast<HighType>(lhs) >= static_cast<HighType>(rhs);
}

template< typename T, typename U >
inline bool check_greater_equal( const T & lhs, const U & rhs, const std::false_type & )
{
   return lhs >= rhs;
}



template< typename T >
void check( const char * const expression, const char * const filename, int line, T failFunc )
{
   std::stringstream ss;
   ss << "Assertion failed!\nFile:       " << filename << ":" << line << '\n'
      << "Expression: " << expression << '\n';

   failFunc( ss.str() );
}

template< typename T, typename U >
void check_nullptr( T * p, const char * const ptrExpression, const char * const filename, int line, U failFunc )
{
   std::stringstream ss;
   ss << "Assertion failed!\n"
      << "File:       " << filename << ":" << line << '\n'
      << "Expression: " << ptrExpression << " == 0\n"
      << "Value:      " << ptrExpression << " = ";

   printValue( ss, p ) << '\n';

   failFunc( ss.str() );
}

template< typename T, typename U >
void check_nullptr( const shared_ptr<T> & p, const char * const ptrExpression, const char * const filename, int line, U failFunc )
{
   std::stringstream ss;
   ss << "Assertion failed!\n"
      << "File:       " << filename << ":" << line << '\n'
      << "Expression: " << ptrExpression << " == 0\n"
      << "Value:      " << ptrExpression << " = ";

   printValue( ss, p ) << '\n';

   failFunc( ss.str() );
}

template< typename T >
void check_not_nullptr( const char * const ptrExpression, const char * const filename, int line, T failFunc )
{
   std::stringstream ss;
   ss << "Assertion failed!\n"
      << "File:       " << filename << ":" << line << '\n'
      << "Expression: " << ptrExpression << " != 0" << '\n';
   failFunc( ss.str() );
}

template< typename T, typename U, typename V >
void check_equal( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                  const char * const filename, int line, V failFunc)
{
   printErrorAndExit(lhs, rhs, lhsExpression, rhsExpression, " == ", filename, line, failFunc);
}

template< typename T, typename U, typename V >
void check_unequal( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                    const char * const filename, int line, V failFunc )
{
   printErrorAndExit(lhs, rhs, lhsExpression, rhsExpression, " != ", filename, line, failFunc);
}

template< typename T, typename U, typename V >
void check_float_equal( const T & lhs, const U & rhs,
                        const char * const lhsExpression, const char * const rhsExpression,
                        const char * const filename, int line, V failFunc )
{
   int length = static_cast<int>( std::max( std::strlen( lhsExpression ), std::strlen( rhsExpression ) ) );
   std::stringstream ss;
   ss << "Assertion failed!\n"
      << "File:       " << filename << ":" << line << '\n'
      << "Expression: " << lhsExpression << " == " << rhsExpression << '\n'
      //<< "ULP:        " << distance << '\n'
      << "Values:     "  << std::setw(length) << std::setfill(' ') << lhsExpression << " = ";

   printValue( ss, lhs ) << '\n';

   ss << "            " << std::setw(length) << std::setfill(' ') << rhsExpression << " = ";

   printValue( ss, rhs ) << '\n';

   failFunc( ss.str() );
}

template< typename T, typename U, typename V >
void check_float_unequal( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                          const char * const filename, int line, V failFunc )
{
   int length = static_cast<int>( std::max( std::strlen( lhsExpression ), std::strlen( rhsExpression ) ) );
   std::stringstream ss;
   ss << "Assertion failed!\n"
      << "File:       " << filename << ":" << line << '\n'
      << "Expression: " << lhsExpression << " != " << rhsExpression << '\n'
      //<< "ULP:        " << distance << '\n'
      << "Values:     " << std::setw(length) << std::setfill(' ') << lhsExpression << " = ";

   printValue( ss, lhs ) << '\n';

   ss << "            " << std::setw(length) << std::setfill(' ') << rhsExpression << " = ";

   printValue( ss, rhs ) << '\n';

   failFunc( ss.str() );
}

template< typename T, typename U, typename V >
void check_float_equal_eps( const T & lhs, const U & rhs,
                        const char * const lhsExpression, const char * const rhsExpression,
                        const char * const filename, int line, V failFunc,
                        const typename VectorTrait<typename math::MathTrait<T,U>::LowType>::OutputType epsilon )
{
   int length = static_cast<int>( std::max( std::strlen( lhsExpression ), std::strlen( rhsExpression ) ) );
   std::stringstream ss;
   ss << "Assertion failed!\n"
      << "File:       " << filename << ":" << line << '\n'
      << "Expression: " << lhsExpression << " == " << rhsExpression << '\n'
      << "Epsilon:    " << epsilon << '\n'
      //<< "ULP:        " << distance << '\n'
      << "Values:     " << std::setw(length) << std::setfill(' ') << lhsExpression << " = ";

   printValue( ss, lhs ) << '\n';

   ss << "            " << std::setw(length) << std::setfill(' ') << rhsExpression << " = ";

   printValue( ss, rhs ) << '\n';

   failFunc( ss.str() );
}

template< typename T, typename U, typename V >
void check_float_unequal_eps( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                          const char * const filename, int line, V failFunc,
                          const typename VectorTrait<typename math::MathTrait<T,U>::LowType>::OutputType epsilon )
{
   int length = static_cast<int>( std::max( std::strlen( lhsExpression ), std::strlen( rhsExpression ) ) );
   std::stringstream ss;
   ss << "Assertion failed!\n"
      << "File:       " << filename << ":" << line << '\n'
      << "Expression: " << lhsExpression << " != " << rhsExpression << '\n'
      << "Epsilon:    " << epsilon << '\n'
      //<< "ULP:        " << distance << '\n'
      << "Values:     " << std::setw(length) << std::setfill(' ') << lhsExpression << " = ";

   printValue( ss, lhs ) << '\n';

   ss << "            " << std::setw(length) << std::setfill(' ') << rhsExpression << " = ";

   printValue( ss, rhs ) << '\n';

   failFunc( ss.str() );
}

template< typename T, typename U, typename V >
void check_identical( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                      const char * const filename, int line, V failFunc )
{
   printErrorAndExit(lhs, rhs, lhsExpression, rhsExpression, " == ", filename, line, failFunc);
}

template< typename T, typename U, typename V >
void check_not_identical( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                          const char * const filename, int line, V failFunc )
{
   printErrorAndExit(lhs, rhs, lhsExpression, rhsExpression, " != ", filename, line, failFunc);
}

template< typename T, typename U, typename V >
void check_less( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                 const char * const filename, int line, V failFunc )
{
   printErrorAndExit(lhs, rhs, lhsExpression, rhsExpression, " < ", filename, line, failFunc);
}

template< typename T, typename U, typename V >
void check_greater( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                    const char * const filename, int line, V failFunc )
{
   printErrorAndExit(lhs, rhs, lhsExpression, rhsExpression, " > ", filename, line, failFunc);
}

template< typename T, typename U, typename V >
void check_less_equal( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                       const char * const filename, int line, V failFunc )
{
   printErrorAndExit(lhs, rhs, lhsExpression, rhsExpression, " <= ", filename, line, failFunc);
}

template< typename T, typename U, typename V >
void check_greater_equal( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                          const char * const filename, int line, V failFunc )
{
   printErrorAndExit(lhs, rhs, lhsExpression, rhsExpression, " >= ", filename, line, failFunc);
}



template< typename T, typename U, typename V >
void printErrorAndExit( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                        const char * const opString, const char * const filename, int line, V failFunc )
{
   int length = static_cast<int>( std::max( std::strlen( lhsExpression ), std::strlen( rhsExpression )  ) );
   std::ostringstream ss;
   ss << "Assertion failed!\n"
      << "File:       " << filename << ":" << line << '\n'
      << "Expression: " << lhsExpression << opString << rhsExpression << '\n'
      << "Values:     " << std::setw(length) << std::setfill(' ') << lhsExpression << " = ";

   printValue( ss, lhs ) << '\n';

   ss << "            " << std::setw(length) << std::setfill(' ') << rhsExpression << " = ";

   printValue( ss, rhs ) << '\n';

   failFunc( ss.str() );
}

template< typename T >
std::ostream & printValue( std::ostream & os, const T & value )
{
   return printValue( os, value, has_left_shift<std::ostream&, T&>() );
}

template< typename T >
std::ostream & printValue( std::ostream & os, const T & value, const std::true_type & )
{
   return os << value;
}

template< typename T >
std::ostream & printValue( std::ostream & os, const T & /*value*/, const std::false_type & )
{
   return os << "[N/A: Type can not be streamed to std::ostream]";
}

template< typename T >
std::ostream & printValue( std::ostream & os, const T * value )
{
   return os << value;
}

template< typename T >
std::ostream & printValue( std::ostream & os, T * value )
{
   return os << value;
}

} // namespace check_functions_detail
} // namespace debug
} // namespace walberla

/// \endcond
