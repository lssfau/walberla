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
//! \file Attachable.h
//! \ingroup pe
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the Attachable class
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <core/NonCopyable.h>
#include <pe/Types.h>
#include <core/ptrvector/policies/NoDelete.h>
#include <core/ptrvector/PtrVector.h>
#include <pe/Types.h>

#include "core/UniqueID.h"

namespace walberla {
namespace pe {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Attachable interface class.
 * \ingroup core
 *
 * The Attachable class is the base class for the Attachable concept of the physics engine.
 * Classes deriving from this interface class (which are simply called Attachables) can be
 * attached to rigid bodies and are notified in case the rigid body is destroyed. This
 * concept is for example used for all force generators (see the ForceGenerator class
 * description).
 */
class Attachable : private NonCopyable
{
public:
   //*************************************************************************************************
   //! Type codes of the attachables.
   enum AttachableType {
      gravityType  = 1,  //!< Code for the Gravity force generator.
      springType   = 2   //!< Code for the Spring force generator.
   };
   //*************************************************************************************************
protected:
   //**Type definitions****************************************************************************
   typedef PtrVector<RigidBody,NoDelete>  Bodies;  //!< Vector for attached rigid bodies.
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit Attachable( AttachableType type, id_t sid );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   virtual ~Attachable();
   //@}
   //**********************************************************************************************

public:

   //**Type definitions****************************************************************************
   typedef Bodies::Iterator       Iterator;       //!< Iterator over the currently attached bodies.
   typedef Bodies::ConstIterator  ConstIterator;  //!< ConstIterator over the currently attached bodies.
   //**********************************************************************************************

   //**MPI functions*******************************************************************************
   /*!\name MPI functions */
   //@{
   virtual bool isRemote() const = 0;
   //@}
   //**********************************************************************************************

   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline AttachableType getType()     const;
   inline id_t           getSystemID() const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline Iterator      begin();
   inline ConstIterator begin() const;
   inline Iterator      end  ();
   inline ConstIterator end  () const;
   inline size_t        size () const;
   //@}
   //**********************************************************************************************

protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   const AttachableType type_;  //!< The type code of the attachable.
   id_t sid_;                   //!< The unique system-specific attachable ID.
   Bodies bodies_;              //!< Vector of attached rigid bodies.
   //@}
   //**********************************************************************************************

private:
   //**Attachable setup functions******************************************************************
   /*! \cond internal */
   friend void detach( AttachableID attachable );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  ATTACHABLE SETUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Attachable setup functions */
//@{
void detach( AttachableID attachable );
//@}
//*************************************************************************************************

//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the type of the attachable.
 *
 * \return The type of the attachable.
 */
inline Attachable::AttachableType Attachable::getType() const
{
   return type_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the unique system-specific ID of the attachable.
 *
 * \return The system-specific ID.
 */
inline id_t Attachable::getSystemID() const
{
   return sid_;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns an iterator to the first attached rigid body.
 *
 * \return Iterator to the first attached rigid body.
 */
inline Attachable::Iterator Attachable::begin()
{
   return bodies_.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first attached rigid body.
 *
 * \return Iterator to the first attached rigid body.
 */
inline Attachable::ConstIterator Attachable::begin() const
{
   return bodies_.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last attached rigid body.
 *
 * \return Iterator just past the last attached rigid body.
 */
inline Attachable::Iterator Attachable::end()
{
   return bodies_.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last attached rigid body.
 *
 * \return Iterator just past the last attached rigid body.
 */
inline Attachable::ConstIterator Attachable::end() const
{
   return bodies_.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of attached rigid bodies.
 *
 * \return The number of attached rigid bodies
 */
inline size_t Attachable::size() const
{
   return bodies_.size();
}
//*************************************************************************************************

//=================================================================================================
//
//  ATTACHABLE TYPE UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Attachable type utility functions */
//@{
std::ostream& operator<<( std::ostream& os, Attachable::AttachableType type );
//@}
//*************************************************************************************************

} // namespace pe
}  // namespace walberla
