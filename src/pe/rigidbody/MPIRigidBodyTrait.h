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
//! \file MPIRigidBodyTrait.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockID.h"

#include "Owner.h"

#include <iostream>
#include <set>
#include <vector>

namespace walberla{
namespace pe{

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for a specialization of the RigidBodyTrait class template for MPI parallel solvers.
 *
 * This class adds functionality to track which process owns the rigid body and which processes
 * have remote copies of the body in case it is owned by this process.
 */
class MPIRigidBodyTrait
{
protected:
   //**Type definitions****************************************************************************
   typedef std::vector<Owner>             ShadowOwners;       //!< Vector for remote MPI processes the rigid body is contained in.
   typedef std::set<BlockID>              BlockStates;        //!< Stores the information if neighbor block knows about object.
   //**********************************************************************************************

public:
   //**Type definitions****************************************************************************
   typedef ShadowOwners::iterator        ShadowOwnersIterator;       //!< Iterator over the connected processes.
   typedef ShadowOwners::const_iterator  ConstShadowOwnersIterator;  //!< ConstIterator over the connected processes.
   typedef size_t                        SizeType;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit MPIRigidBodyTrait( );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   virtual ~MPIRigidBodyTrait() = default;
   //@}
   //**********************************************************************************************

public:
   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline Owner& getOwner();
   inline const Owner& getOwner() const;
   inline void  setOwner(const Owner& owner);
   //@}
   //**********************************************************************************************

   //** Functions to manipulate Shadow Owners *****************************************************
   /*!\name shadow owners functions */
   //@{
   inline void                 registerShadowOwner  ( const Owner& owner );
   inline void                 deregisterShadowOwner( const Owner& owner );
   inline bool                 isShadowOwnerRegistered( const Owner& owner ) const;
   inline bool                 hasShadowOwners       () const;
   inline ShadowOwnersIterator      beginShadowOwners();
   inline ConstShadowOwnersIterator beginShadowOwners() const;
   inline ShadowOwnersIterator      endShadowOwners  ();
   inline ConstShadowOwnersIterator endShadowOwners  () const;
   inline SizeType             sizeShadowOwners () const;
   inline void                 clearShadowOwners();
   //@}
   //**********************************************************************************************

   //** functions to store neighboring information ************************************************
   /*!\name cache information functions */
   //@{
   inline void                 setBlockState  ( const BlockID& id );
   inline void                 unsetBlockState( const BlockID& id );
   inline bool                 getBlockState  ( const BlockID& id ) const;
   inline size_t               getBlockStateSize() const;
   inline void                 clearBlockStates();
   //@}
   //**********************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   ShadowOwners shadowOwners_;    //!< Vector of all processes the rigid body intersects with.
   BlockStates  blockStates_;
   Owner        owner_;    //!< Rank of the process owning the rigid body.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the MPIRigidBodyTrait constructor.
 *
 * \param body The body ID of this rigid body.
 */
inline MPIRigidBodyTrait::MPIRigidBodyTrait( )
   : owner_( )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Registering a new remote process the rigid body is contained in.
 *
 * \param process The remote process to be registered with the rigid body.
 * \return void
 *
 * This function registers the given remote process with the rigid body. In case the process is
 * already registered, it is not registered again. This call has linear complexity.
 */
inline void MPIRigidBodyTrait::registerShadowOwner( const Owner& owner )
{
   if( !isShadowOwnerRegistered( owner ) )
      shadowOwners_.push_back( owner );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deregistering a remote process from the rigid body.
 *
 * \param process The remote process to be deregistered from the rigid body.
 * \return void
 *
 * This function deregisters the given remote process from the rigid body. This call has linear
 * complexity.
 */
inline void MPIRigidBodyTrait::deregisterShadowOwner( const Owner& owner )
{
   ShadowOwnersIterator pos = std::find( shadowOwners_.begin(), shadowOwners_.end(), owner );
   if( pos != shadowOwners_.end() )
      shadowOwners_.erase( pos );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks whether the given remote process is registered with the rigid body.
 *
 * \param process The remote process that possibly registered with the rigid body.
 * \return \a true if the given process is registered with the rigid body, \a false if not.
 *
 * This call has linear complexity.
 */
inline bool MPIRigidBodyTrait::isShadowOwnerRegistered( const Owner& owner ) const
{
   return std::find( shadowOwners_.begin(), shadowOwners_.end(), owner ) != shadowOwners_.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether any processes are registered with the rigid body.
 *
 * \return \a true if at least one process is registered with the rigid body, \a false if not.
 */
inline bool MPIRigidBodyTrait::hasShadowOwners() const
{
   return !shadowOwners_.empty();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first remote process the rigid body is contained in.
 *
 * \return Iterator to the first remote process.
 */
inline MPIRigidBodyTrait::ShadowOwnersIterator MPIRigidBodyTrait::beginShadowOwners()
{
   return shadowOwners_.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first remote process the rigid body is contained in.
 *
 * \return Iterator to the first remote process.
 */
inline MPIRigidBodyTrait::ConstShadowOwnersIterator MPIRigidBodyTrait::beginShadowOwners() const
{
   return shadowOwners_.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last process the rigid body is contained in.
 *
 * \return Iterator just past the last remote process.
 */
inline MPIRigidBodyTrait::ShadowOwnersIterator MPIRigidBodyTrait::endShadowOwners()
{
   return shadowOwners_.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last process the rigid body is contained in.
 *
 * \return Iterator just past the last remote process.
 */
inline MPIRigidBodyTrait::ConstShadowOwnersIterator MPIRigidBodyTrait::endShadowOwners() const
{
   return shadowOwners_.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of registered processes.
 *
 * \return The number of registered processes.
 */
inline MPIRigidBodyTrait::SizeType MPIRigidBodyTrait::sizeShadowOwners() const
{
   return shadowOwners_.size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing all registered remote processes from the rigid body.
 *
 * \return void
 *
 * This function clears all references to remote processes the rigid body may be contained in.
 */
inline void MPIRigidBodyTrait::clearShadowOwners()
{
   shadowOwners_.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the rank of the owner process.
 *
 * \return The rank of the owner process.
 *
 * If not yet set by the collision system this function returns -1.
 */
inline Owner& MPIRigidBodyTrait::getOwner()
{
   return owner_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the rank of the owner process.
 *
 * \return The rank of the owner process.
 *
 * If not yet set by the collision system this function returns -1.
 */
inline const Owner& MPIRigidBodyTrait::getOwner() const
{
   return owner_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the rank of the owner process.
 *
 * \return The rank of the owner process.
 *
 * If not yet set by the collision system this function returns -1.
 */
inline void MPIRigidBodyTrait::setOwner(const Owner& owner)
{
   owner_ = owner;
}
//*************************************************************************************************

inline void MPIRigidBodyTrait::setBlockState  ( const BlockID& id )
{
   blockStates_.insert( id );
}

inline void MPIRigidBodyTrait::unsetBlockState( const BlockID& id )
{
   blockStates_.erase( id );
}

inline bool MPIRigidBodyTrait::getBlockState  ( const BlockID& id ) const
{
   return blockStates_.find( id ) != blockStates_.end();
}
inline void MPIRigidBodyTrait::clearBlockStates()
{
   blockStates_.clear();
}

inline size_t MPIRigidBodyTrait::getBlockStateSize() const
{
   return blockStates_.size();
}

}  // namespace pe
}  // namespace walberla
