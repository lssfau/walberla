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
//! \file Algorithm.h
//! \ingroup core
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/is_base_of.hpp>

namespace walberla {

//=================================================================================================
//
//  POLYMORPHIC COUNT
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Counts the pointer to objects with dynamic type \a D.
 *
 * \param first Iterator to the first pointer of the pointer range.
 * \param last Iterator to the pointer one past the last pointer of the pointer range.
 * \return The number of objects with dynamic type \a D.
 *
 * This function traverses the range \f$ [first,last) \f$ of pointers to objects with static
 * type \a S and counts all polymorphic pointers to objects of dynamic type \a D. Note that
 * in case \a D is not a type derived from \a S, a compile time error is created!
 */
template< typename D    // Dynamic type of the objects
        , typename S >  // Static type of the objects
inline size_t polymorphicCount( S *const * first, S *const * last )
{
   static_assert(boost::is_base_of<S, D>::value && !boost::is_base_of<D, S>::value, "D has to be strictly derived from S");

   size_t count( 0 );
   for( S *const * it=first; it!=last; ++it )
      if( dynamic_cast<D*>( *it ) ) ++count;
   return count;
}
//*************************************************************************************************




//=================================================================================================
//
//  POLYMORPHIC FIND
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Finds the next pointer to an object with dynamic type \a D.
 *
 * \param first Iterator to the first pointer of the pointer range.
 * \param last Iterator to the pointer one past the last pointer of the pointer range.
 * \return The next pointer to an object with dynamic type \a D.
 *
 * This function traverses the range \f$ [first,last) \f$ of pointers to objects with static
 * type \a S until it finds the next polymorphic pointer to an object of dynamic type \a D.
 * Note that in case \a D is not a type derived from \a S, a compile time error is created!
 */
template< typename D    // Dynamic type of the objects
        , typename S >  // Static type of the objects
inline S *const * polymorphicFind( S *const * first, S *const * last )
{
   static_assert(boost::is_base_of<S, D>::value && !boost::is_base_of<D, S>::value, "D has to be strictly derived from S");

   while( first != last && !dynamic_cast<D*>( *first ) ) ++first;
   return first;
}
//*************************************************************************************************

} // namespace
