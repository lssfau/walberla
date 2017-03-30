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
//! \file BodyStorage.h
//! \ingroup pe
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <functional>
#include <map>
#include <vector>
#include <core/NonCopyable.h>
#include <core/debug/Debug.h>
#include <core/ptrvector/policies/PtrDelete.h>
#include <core/ptrvector/PtrVector.h>
#include <pe/rigidbody/RigidBody.h>
#include <pe/Types.h>

#include <boost/function.hpp>

namespace walberla {
namespace pe {




//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Body storage of the rigid body simulation world.
 * \ingroup pe
 *
 * A BodyStorage is a data structure for storing rigid bodies. It supports efficient insertion and
 * deletion operations.
 */
class BodyStorage : private NonCopyable
{
public:
   //**Type definitions****************************************************************************
   //! Container for the bodies contained in the simulation world.
   typedef PtrVector<BodyType, PtrDelete>  Bodies;
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   typedef Bodies::SizeType       SizeType;       //!< Size type of the body storage.
   typedef Bodies::Iterator       Iterator;       //!< Iterator over non-const bodies.
   typedef Bodies::ConstIterator  ConstIterator;  //!< Iterator over constant bodies.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline BodyStorage();
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~BodyStorage();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool          isEmpty () const;
   inline SizeType      size    () const;
   template< typename T >
   inline SizeType      size    () const;
   inline Iterator      begin   ();
   inline ConstIterator begin   () const;
   template< typename T >
   inline typename Bodies::template CastIterator<T> begin();
   template< typename T >
   inline typename Bodies::template ConstCastIterator<T> begin() const;
   inline Iterator      end     ();
   inline ConstIterator end     () const;
   template< typename T >
   inline typename Bodies::template CastIterator<T> end();
   template< typename T >
   inline typename Bodies::template ConstCastIterator<T> end() const;
   inline BodyID        at      ( SizeType index );
   inline ConstBodyID   at      ( SizeType index ) const;
   inline Iterator      find    ( id_t sid );
   inline ConstIterator find    ( id_t sid ) const;
   inline Iterator      find    ( ConstBodyID body );
   inline ConstIterator find    ( ConstBodyID body ) const;
   inline void          validate();
   //@}
   //**********************************************************************************************

   //**Add/Remove functions************************************************************************
   /*!\name Add/Remove functions */
   //@{
   inline void          add     ( BodyID body );
   inline void          remove  ( const id_t sid );
   inline void          remove  ( BodyID body );
   inline ConstIterator remove  ( ConstIterator pos );
   inline Iterator      remove  ( Iterator pos );
   inline void          release ( const id_t sid );
   inline void          release ( BodyID body );
   inline ConstIterator release  ( ConstIterator pos );
   inline Iterator      release  ( Iterator pos );
   inline void          clear   ();
   //@}
   //**********************************************************************************************

   //**Callbacks************************************************************************
   /*!\name Callbacks */
   //@{
   inline void          registerAddCallback     ( const std::string name, const boost::function<void (BodyID)>& func );
   inline void          deregisterAddCallback   ( const std::string name );
   inline void          clearAddCallbacks       ( );

   inline void          registerRemoveCallback     ( const std::string name, const boost::function<void (BodyID)>& func );
   inline void          deregisterRemoveCallback   ( const std::string name );
   inline void          clearRemoveCallbacks       ( );
   //@}
   //**********************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Bodies bodies_;  //!< The rigid bodies contained in the simulation world.
   std::map<id_t, SizeType> bodyIDs_;   //!< The association of system IDs to rigid bodies.

   std::map< std::string, boost::function<void (BodyID)> > addCallbacks_;
   std::map< std::string, boost::function<void (BodyID)> > removeCallbacks_;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The standard constructor.
 */
inline BodyStorage::BodyStorage()
	: bodies_( 1000 )
	, bodyIDs_()
{
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor.
 *
 * The destructor clears all rigid bodies from the storage before destructing it.
 */

inline BodyStorage::~BodyStorage()
{
   WALBERLA_ASSERT_EQUAL(addCallbacks_.size(), 0, "Still add callbacks registered!");
   WALBERLA_ASSERT_EQUAL(removeCallbacks_.size(), 0, "Still remove callbacks registered!");
	clear();
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns \a true if the body storage contains no rigid bodies.
 *
 * \return \a true if the body storage is empty, \a false if it is not.
 */

inline bool BodyStorage::isEmpty() const
{
   return bodies_.isEmpty();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of rigid bodies contained in the body storage.
 *
 * \return The number of rigid bodies.
 */

inline BodyStorage::SizeType BodyStorage::size() const
{
   return bodies_.size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of rigid bodies of type \a T contained in the body storage.
 *
 * \return The number of rigid bodies of type \a T.
 *
 * \b Note: The total number of objects of type \a T is not cached but recalculated each time the
 * function is called. Using the templated version of size() to calculate the total number objects
 * of type \a T is therefore more expensive than using the non-template version of size() to get
 * the total number of pointers in the vector!
 */

template< typename T >  // Cast type
inline BodyStorage::SizeType BodyStorage::size() const
{
   return bodies_.template size<T>();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first contained rigid body.
 *
 * \return Iterator to the first contained rigid body.
 */

inline BodyStorage::Iterator BodyStorage::begin()
{
   return bodies_.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator to the first contained rigid body.
 *
 * \return Constant iterator to the first contained rigid body.
 */

inline BodyStorage::ConstIterator BodyStorage::begin() const
{
   return bodies_.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first contained rigid body.
 *
 * \return Iterator to the first contained rigid body.
 */

template< typename T >  // Cast Type
inline BodyStorage::Bodies::template CastIterator<T> BodyStorage::begin()
{
   return bodies_.template begin<T>();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator to the first contained rigid body.
 *
 * \return Constant iterator to the first contained rigid body.
 */

template< typename T >  // Cast Type
inline BodyStorage::Bodies::template ConstCastIterator<T> BodyStorage::begin() const
{
   return bodies_.template begin<T>();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last contained rigid body.
 *
 * \return Iterator just past the last contained rigid body.
 */

inline BodyStorage::Iterator BodyStorage::end()
{
   return bodies_.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator just past the last contained rigid body.
 *
 * \return Constant iterator just past the last contained rigid body.
 */

inline BodyStorage::ConstIterator BodyStorage::end() const
{
   return bodies_.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last contained rigid body.
 *
 * \return Iterator just past the last contained rigid body.
 */

template< typename T >  // Cast Type
inline BodyStorage::Bodies::template CastIterator<T> BodyStorage::end()
{
   return bodies_.template end<T>();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator just past the last contained rigid body.
 *
 * \return Constant iterator just past the last contained rigid body.
 */

template< typename T >  // Cast Type
inline BodyStorage::Bodies::template ConstCastIterator<T> BodyStorage::end() const
{
   return bodies_.template end<T>();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a handle to the indexed rigid body.
 *
 * \param index Access index. The index has to be in the range \f$[0..size-1]\f$.
 * \return Handle to the accessed rigid body.
 *
 * \b Note: No runtime check is performed to ensure the validity of the access index.
 */

inline BodyID BodyStorage::at( SizeType index )
{
   WALBERLA_ASSERT( index < size(), "Invalid body index" );
   return bodies_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant handle to the indexed rigid body.
 *
 * \param index Access index. The index has to be in the range \f$[0..size-1]\f$.
 * \return Constant handle to the accessed rigid body.
 *
 * \b Note: No runtime check is performed to ensure the validity of the access index.
 */

inline ConstBodyID BodyStorage::at( SizeType index ) const
{
   WALBERLA_ASSERT( index < size(), "Invalid body index" );
   return bodies_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Finding a rigid body with a certain unique system-specific ID.
 *
 * \param sid The unique system-specific ID for the search.
 * \return Iterator to the body with system-specific ID \a sid or an iterator just past the end.
 *
 * This function finds the rigid body with the system-specific ID \a sid. In case the rigid body
 * is found, the function returns an iterator to the body. Otherwise, the function returns an
 * iterator just past the end of the last body contained in the body storage.
 */

inline BodyStorage::Iterator BodyStorage::find( id_t sid )
{
   std::map<id_t, SizeType>::const_iterator pos = bodyIDs_.find( sid );
   if( pos == bodyIDs_.end() )
      return bodies_.end();

   return bodies_.begin() + static_cast<Bodies::Iterator::difference_type>(pos->second);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Finding a rigid body with a certain unique system-specific ID.
 *
 * \param sid The unique system-specific ID for the search.
 * \return Constant iterator to the body with system-specific ID \a sid or a constant iterator just past the end.
 *
 * This function finds the rigid body with the system-specific ID \a sid. In case the rigid body
 * is found, the function returns a constant iterator to the body. Otherwise, the function returns
 * a constant iterator just past the end of the last body contained in the body storage.
 */

inline BodyStorage::ConstIterator BodyStorage::find( id_t sid ) const
{
   std::map<id_t, SizeType>::const_iterator pos = bodyIDs_.find( sid );
   if( pos == bodyIDs_.end() )
      return bodies_.end();

   return bodies_.begin() + static_cast<Bodies::Iterator::difference_type>(pos->second);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Finding a specific rigid body in the body storage.
 *
 * \param body The given rigid body for the search.
 * \return Iterator to the given rigid body or an iterator just past the end.
 *
 * This function finds the rigid body in the body storage. In case the rigid body is found,
 * the function returns an iterator to the body. Otherwise, the function returns an iterator
 * just past the end of the last body contained in the body storage.
 */

inline BodyStorage::Iterator BodyStorage::find( ConstBodyID body )
{
   return find( body->getSystemID() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Finding a specific rigid body in the body storage.
 *
 * \param body The given rigid body for the search.
 * \return Constant iterator to the given rigid body or a constant iterator just past the end.
 *
 * This function finds the rigid body in the body storage. In case the rigid body is found,
 * the function returns a constant iterator to the body. Otherwise, the function returns a
 * constant iterator just past the end of the last body contained in the body storage.
 */

inline BodyStorage::ConstIterator BodyStorage::find( ConstBodyID body ) const
{
   return find( body->getSystemID() );
}
//*************************************************************************************************




//=================================================================================================
//
//  ADD/REMOVE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Adding a rigid body to the body storage.
 *
 * \param body The new rigid body to be added to the body storage.
 * \return void
 *
 * This function adds a rigid body to the body storage. Adding bodies with non-unique system ID or
 * adding the same body multiple times results in undefined behaviour. The time complexity is
 * logarithmic unless reallocation occurs.
 */

inline void BodyStorage::add( BodyID body )
{
   WALBERLA_ASSERT( bodyIDs_.find( body->getSystemID() ) == bodyIDs_.end(), "Body with same system ID already added." );
   bodyIDs_[ body->getSystemID() ] = bodies_.size();
   bodies_.pushBack( body );

   for (auto it = addCallbacks_.begin(); it != addCallbacks_.end(); ++it)
   {
      it->second(body);
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing a rigid body from the body storage.
 *
 * \param sid The system id of the rigid body to be removed.
 * \return Iterator to the body after the erased rigid body.
 *
 * This function removes a body from the body storage. \a sid must be a valid system id.
 * Invalidates all iterators pointing at or past
 * the element to be removed. The time complexity is logarithmic unless reallocation occurs.
 * The last element is swapped to the actual position and the length is reduced by one.
 */

inline void BodyStorage::remove( const id_t sid )
{
   std::map<id_t, SizeType>::iterator it = bodyIDs_.find( sid );
   WALBERLA_ASSERT( it != bodyIDs_.end(), "The body's system ID was not registered." );

   // Move last element to deleted place and update the system ID to index mapping.
   SizeType i = it->second;

   for (auto cb = removeCallbacks_.begin(); cb != removeCallbacks_.end(); ++cb)
   {
      cb->second( bodies_[i] );
   }

   bodyIDs_[ bodies_.back()->getSystemID() ] = i;
   std::swap( bodies_[i], bodies_.back() );
   bodyIDs_.erase( it );
   bodies_.popBack();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing a rigid body from the body storage.
 *
 * \param pos The position of the rigid body to be removed.
 * \return Iterator to the body after the erased rigid body.
 *
 * This function removes a body from the body storage. \a pos must be a valid iterator and the
 * rigid body pointer referred to must be valid. Invalidates all iterators pointing at or past
 * the element to be removed. The time complexity is logarithmic unless reallocation occurs.
 */

inline BodyStorage::ConstIterator BodyStorage::remove( ConstIterator pos )
{
   remove( (*pos)->getSystemID() );
   return pos;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing a rigid body from the body storage.
 *
 * \param pos The position of the rigid body to be removed.
 * \return Iterator to the body after the erased rigid body.
 *
 * This function removes a body from the body storage. \a pos must be a valid iterator and the
 * rigid body pointer referred to must be valid. Invalidates all iterators pointing at or past
 * the element to be removed. The time complexity is logarithmic unless reallocation occurs.
 */

inline BodyStorage::Iterator BodyStorage::remove( Iterator pos )
{
   remove( (*pos)->getSystemID() );
   return pos;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing a rigid body from the body storage.
 *
 * \param body A handle of the rigid body to be removed.
 * \return void
 *
 * This function removes a body from the body storage. \a body must be a valid rigid body poitner
 * and must be registered in the body storage. Invalidates all iterators pointing at or past
 * the element to be removed. The time complexity is logarithmic unless reallocation occurs.
 */

inline void BodyStorage::remove( BodyID body )
{
   remove( body->getSystemID() );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Release a rigid body from the body storage.
 *
 * \param sid The system id of the rigid body to be released.
 * \return Iterator to the body after the erased rigid body.
 *
 * This function releases a body from the body storage. The released rigid body is not destroyed.
 * \a sid must be a valid rigid body system id
 * and must be registered in the body storage. Invalidates all iterators pointing at or past
 * the element to be released. The time complexity is logarithmic unless reallocation occurs.
 * Last element is swapped to the actual position and length is reduced by 1.
 */

inline void BodyStorage::release( const id_t sid )
{
   std::map<id_t, SizeType>::iterator it = bodyIDs_.find( sid );
   WALBERLA_ASSERT( it != bodyIDs_.end(), "The body's system ID was not registered." );

   // Move last element to deleted place and update the system ID to index mapping.
   SizeType i = it->second;

   for (auto cb = removeCallbacks_.begin(); cb != removeCallbacks_.end(); ++cb)
   {
      cb->second( bodies_[i] );
   }

   bodyIDs_[ bodies_.back()->getSystemID() ] = i;
   std::swap( bodies_[i], bodies_.back() );
   bodyIDs_.erase( it );
   bodies_.releaseBack();
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Release a rigid body from the body storage.
 *
 * \param pos The position of the rigid body to be released.
 * \return Iterator to the body after the erased rigid body.
 *
 * This function releases a body from the body storage. The released rigid body is not destroyed.
 * \a body must be a valid rigid body pointer
 * and must be registered in the body storage. Invalidates all iterators pointing at or past
 * the element to be released. The time complexity is logarithmic unless reallocation occurs.
 */

inline BodyStorage::ConstIterator BodyStorage::release( ConstIterator pos )
{
   release( (*pos)->getSystemID() );
   return pos;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Releasing a rigid body from the body storage.
 *
 * \param pos The position of the rigid body to be released.
 * \return Iterator to the body after the erased rigid body.
 *
 * This function releases a body from the body storage. The released rigid body is not destroyed.
 * \a body must be a valid rigid body pointer
 * and must be registered in the body storage. Invalidates all iterators pointing at or past
 * the element to be released. The time complexity is logarithmic unless reallocation occurs.
 */

inline BodyStorage::Iterator BodyStorage::release( Iterator pos )
{
   release( (*pos)->getSystemID() );
   return pos;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Release a rigid body from the body storage.
 *
 * \param body A handle of the rigid body to be released.
 * \return void
 *
 * This function releases a body from the body storage. The released rigid body is not destroyed.
 * \a body must be a valid rigid body pointer
 * and must be registered in the body storage. Invalidates all iterators pointing at or past
 * the element to be released. The time complexity is logarithmic unless reallocation occurs.
 */

inline void BodyStorage::release( BodyID body )
{
   release( body->getSystemID() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing all rigid bodies from the body storage.
 *
 * \return void
 *
 * This function removes all bodies from the body storage. The rigid bodies do not have to be
 * valid anymore that is they can already be deallocated. Invalidates all iterators of this
 * container.
 */

inline void BodyStorage::clear()
{
   for (auto bodyIt = bodies_.begin(); bodyIt != bodies_.end(); ++bodyIt)
   {
      for (auto cb = removeCallbacks_.begin(); cb != removeCallbacks_.end(); ++cb)
      {
         cb->second( *bodyIt );
      }
   }

   bodyIDs_.clear();
   bodies_.clear();
}
//*************************************************************************************************

inline void          BodyStorage::registerAddCallback     ( const std::string name, const boost::function<void (BodyID)>& func )
{
   WALBERLA_ASSERT_EQUAL(addCallbacks_.find(name), addCallbacks_.end(), "Callback '" << name << "' already exists!");
   addCallbacks_.insert( std::make_pair(name, func) );
}

inline void          BodyStorage::deregisterAddCallback   ( const std::string name )
{
   auto res = addCallbacks_.find( name );
   if (res != addCallbacks_.end() )
   {
      addCallbacks_.erase( res );
   }
}

inline void          BodyStorage::clearAddCallbacks       ( )
{
   addCallbacks_.clear();
}

inline void          BodyStorage::registerRemoveCallback     ( const std::string name, const boost::function<void (BodyID)>& func )
{
   WALBERLA_ASSERT_EQUAL(removeCallbacks_.find(name), removeCallbacks_.end(), "Callback '" << name << "' already exists!");
   removeCallbacks_.insert( std::make_pair(name, func) );
}

inline void          BodyStorage::deregisterRemoveCallback   ( const std::string name )
{
   auto res = removeCallbacks_.find( name );
   if (res != removeCallbacks_.end() )
   {
      removeCallbacks_.erase( res );
   }
}

inline void          BodyStorage::clearRemoveCallbacks       ( )
{
   removeCallbacks_.clear();
}


//*************************************************************************************************
/*!\brief Validating the correctness of the body storage data structure.
 *
 * \return void
 *
 * This function validates the data structure in linear time and space. If validation fails
 * assertions are triggered unless the pe is compiled in release mode.
 */

inline void BodyStorage::validate()
{
   std::vector<bool> tmp(bodies_.size());
   std::map<id_t, SizeType>::iterator it = bodyIDs_.begin();
   while( it != bodyIDs_.end() ) {
      WALBERLA_ASSERT(tmp[it->second] == false, "Two system IDs point to the same storage index.");
      tmp[it->second] = true;
      WALBERLA_ASSERT(bodies_[it->second]->getSystemID() == it->first, "The mapping of system ID to storage index is wrong since the system ID does not match with the stored body.");
      ++it;
   }

   WALBERLA_ASSERT( bodyIDs_.size() == bodies_.size(), "The mapping is missing some system IDs." );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Compare if two BodyStorages are equal.
 *
 * Since BodyStorages are uncopyable two BodyStorages are considered equal if their adresses are equal.
 */
inline bool operator==(const BodyStorage& lhs, const BodyStorage& rhs) {return &lhs == &rhs;}

}  // namespace pe
}  // namespace walberla
