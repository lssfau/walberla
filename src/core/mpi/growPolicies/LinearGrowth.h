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
//! \file LinearGrowth.h
//! \ingroup core
//! \brief Header file for the LinearGrowth policy classes.
//!
//! Copyright (C) 2009 Klaus Iglberger
//! Taken from "pe Physics Engine"
//
//======================================================================================================================

#pragma once

#include <algorithm>


/// \cond internal

namespace walberla {
namespace mpi {


//======================================================================================================================
//
//  CLASS DEFINITION
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Linear growth policy class.
// \ingroup core
//
// The LinearGrowth policy class implements a linear growth strategy. It can be customized for
// any purpose: the \a Growth template argument specifies the factor of the size growth.
*/
template< size_t Growth >
struct LinearGrowth
{
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t operator()( size_t oldSize, size_t minSize ) const;
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************


//**********************************************************************************************************************
/*! \cond internal */
/*!\brief Specialization of the linear growth policy for 0 growth.
// \ingroup util
*/
template<>
struct LinearGrowth<0>;
/*! \endcond */
//**********************************************************************************************************************


//**********************************************************************************************************************
/*! \cond internal */
/*!\brief Specialization of the linear growth policy for 1 growth.
// \ingroup util
*/
template<>
struct LinearGrowth<1>;
/*! \endcond */
//**********************************************************************************************************************




//======================================================================================================================
//
//  UTILITY FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Returns a new size depending on the given old size and the required minimum size.
//
// \param old The old size.
// \param minimum The required minimum size.
// \return The new size (at least the required minimum size).
*/
template< size_t Growth >
inline size_t LinearGrowth<Growth>::operator()( size_t old, size_t minimum ) const
{
   const size_t needed( std::max<size_t>( old*Growth, minimum ) );
   return ( ( needed )?( 4 * ( (needed-1)/4+1 ) ):( 0 ) );
}
//**********************************************************************************************************************

} // namespace mpi
} // namespace walberla

/// \endcond
