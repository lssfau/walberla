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
//! \file AttachableCast.h
//! \ingroup pe
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <boost/type_traits/is_base_of.hpp>

namespace walberla {
namespace pe {

//=================================================================================================
//
//  ATTACHABLE CAST OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Attachable cast operators */
//@{
template< typename To, typename From > inline To* static_attachable_cast( From* attachable );
template< typename To, typename From > inline To* dynamic_attachable_cast( From* attachable );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Static cast for attachables.
 *
 * \param attachable The attachable to be cast.
 * \return The casted attachable.
 *
 * The static_attachable_cast function is used exactly as the built-in static_cast operator
 * but for attachables.

   \code
   SphereID     sphere     = createSphere ( 1, 0.0, 0.0, 0.0, 1.0, iron );
   AttachableID attachable = attachGravity( sphere, 0.0, 0.0, -9.81 );
   GravityID    gravity    = static_attachable_cast<Gravity>( attachable );
   \endcode
 */
template< typename To, typename From >
inline To* static_attachable_cast( From* attachable )
{
   static_assert(boost::is_base_of<Attachable, From>::value, "From has to be derived from Attachable!");
   static_assert(boost::is_base_of<Attachable, To>::value, "To has to be derived from Attachable!");
   return static_cast<To*>( attachable );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Dynamic cast for attachables.
 *
 * \param attachable The attachable to be cast.
 * \return The casted attachable.
 *
 * The dynamic_attachable_cast function is used exactly as the built-in dynamic_cast operator
 * but for attachables.

   \code
   AttachableID attachable;
   GravityID gravity = dynamic_attachable_cast<Gravity>( attachable );
   \endcode
 */
template< typename To, typename From >
inline To* dynamic_attachable_cast( From* attachable )
{
   static_assert(boost::is_base_of<Attachable, From>::value, "From has to be derived from Attachable!");
   static_assert(boost::is_base_of<Attachable, To>::value, "To has to be derived from Attachable!");
   return dynamic_cast<To*>( attachable );
}
//*************************************************************************************************

} // namespace pe
}
