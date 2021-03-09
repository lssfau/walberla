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
//! \file SqrtTrait.h
//! \ingroup core
//! \author Klaus Iglberger
//! \brief Header file for the square root trait
//
//======================================================================================================================

#pragma once



namespace walberla {
namespace math {

//======================================================================================================================
//
//  SQUARE ROOT TRAIT
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Base template for the SqrtTrait class.
// \ingroup core
//
// The SqrtTrait template evaluates the return type of the std::sqrt function depending on the
// type of the argument. In case of an integral data type, the return value of the std::sqrt
// function is double, whereas the return type is float for single precision arguments and
// long double for long double precision arguments.
//
// <table border="0" cellspacing="0" cellpadding="1">
//    <tr>
//       <td width="250px"> <b> Template argument T </b> </td>
//       <td width="100px"> <b> Type </b> </td>
//    </tr>
//    <tr>
//       <td>float</td>
//       <td>float</td>
//    </tr>
//    <tr>
//       <td>integral data types and double</td>
//       <td>double</td>
//    </tr>
//    <tr>
//       <td>long double</td>
//       <td>long double</td>
//    </tr>
// </table>
*/
template< typename T >
struct SqrtTrait
{
   using Type = double;  //!< Return type of std::sqrt for integral and double arguments.
};
//**********************************************************************************************************************


//**********************************************************************************************************************
/*! \cond internal */
/*!\brief SqrtTrait<float> specialization.
// \ingroup core
*/
template<>
struct SqrtTrait<float>
{
   using Type = float;  //!< Return type of std::sqrt for float arguments.
};
/*! \endcond */
//**********************************************************************************************************************


//**********************************************************************************************************************
/*! \cond internal */
/*!\brief SqrtTrait<long double> specialization.
// \ingroup core
*/
template<>
struct SqrtTrait<long double>
{
   using Type = long double;  //!< Return type of std::sqrt for long double arguments.
};
/*! \endcond */
//**********************************************************************************************************************



} // namespace math
} // namespace walberla

