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
//! \file Attachable.cpp
//! \ingroup pe
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Source file for the Attachable class
//
//======================================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <pe/attachable/Attachable.h>
#include <pe/attachable/AttachableStorage.h>
#include <core/DataTypes.h>

namespace walberla {
namespace pe {

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor of the Attachable class.
 *
 * \param type The type of the rigid body.
 * \param sid The unique system-specific ID of the attachable.
 */
Attachable::Attachable( AttachableType type, id_t sid )
   : type_(type)
   , sid_(sid)
   , bodies_()
{}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the Attachable class.
 */
Attachable::~Attachable()
{
}
//*************************************************************************************************




//=================================================================================================
//
//  MPI FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\fn bool Attachable::isRemote() const
 * \brief Returns whether the attachable is remote or not.
 *
 * \return \a true in case the attachable is remote, \a false if not.
 *
 * This function returns whether the attachable is remote or not. In case the attachable is
 * attached to at least a single local rigid body the function returns \a false. Otherwise
 * it returns \a true.
 */
//*************************************************************************************************




//=================================================================================================
//
//  ATTACHABLE SETUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Detaches the given attachable.
 * \ingroup core
 *
 * \param attachable The detachable to be detached.
 */
WALBERLA_PUBLIC void detach( AttachableID attachable )
{
   ///TDOD REVIEW
//   typedef AttachableStorage AS;

   // WARNING: Using friend relationship to get non-constant reference of attachable storage.
//   AS& attachablestorage( theCollisionSystem()->attachablestorage_ );

   // Removing the attachable from the attachable storage
//   const AS::Iterator pos( attachablestorage.find( attachable ) );
//   attachablestorage.remove( pos );

   delete attachable;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Global output operator for AttachableType.
 *
 * \param os Reference to the output stream.
 * \param type The AttachableType to be put into the stream.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, Attachable::AttachableType type )
{
   switch( type )
   {
      case Attachable::gravityType: os << "gravity force generator"; break;
      case Attachable::springType:  os << "spring force generator";  break;
      default: WALBERLA_ASSERT( false, "Unknown attachable type" ); break;
   }

   return os;
}
//*************************************************************************************************

} // namespace pe
}  // namespace walberla
