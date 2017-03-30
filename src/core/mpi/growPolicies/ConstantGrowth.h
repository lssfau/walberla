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
//! \file ConstantGrowth.h
//! \ingroup core
//! \brief Header file for the ConstantGrowth policy classes.
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
/*!\brief Constant growth policy class.
// \ingroup core
//
// The ConstantGrowth policy class implements a constant growth strategy. It can be customized
// for any purpose: the \a Growth template argument specifies the constant increase of the given
// size.
*/
template< size_t Growth >
struct ConstantGrowth
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
/*!\brief Specialization of the constant growth policy for 0 growth.
// \ingroup util
*/
template<>
struct ConstantGrowth<0>;
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
inline size_t ConstantGrowth<Growth>::operator()( size_t old, size_t minimum ) const
{
   const size_t needed( std::max<size_t>( old+Growth, minimum ) );
   return ( ( needed )?( 4 * ( (needed-1)/4+1 ) ):( 0 ) );
}
//**********************************************************************************************************************

} // namespace mpi
} // namespace walberla
/// \endcond


