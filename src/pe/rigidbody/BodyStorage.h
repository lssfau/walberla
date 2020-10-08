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
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <core/NonCopyable.h>
#include <core/debug/Debug.h>
#include <pe/rigidbody/RigidBody.h>
#include <pe/rigidbody/RigidBodyCastIterator.h>
#include <pe/rigidbody/RigidBodyIterator.h>
#include <pe/Types.h>

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <vector>

namespace walberla {
namespace pe {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Body storage of the rigid body simulation world.
 *
 * A BodyStorage is a data structure for storing rigid bodies. It supports efficient insertion and
 * deletion operations.
 */
class BodyStorage : private NonCopyable
{
public:
   //**Type definitions****************************************************************************
   //! Container for the bodies contained in the simulation world.
   using VectorContainer = std::vector< std::unique_ptr<RigidBody> >;
   using ConstVectorContainer = std::vector< std::unique_ptr<const RigidBody> >;
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   using size_type             = VectorContainer::size_type;           //!< Size type of the body storage.
   using iterator              = RigidBodyIterator;                    //!< Iterator over non-const bodies.
   using const_iterator        = ConstRigidBodyIterator;               //!< Iterator over constant bodies.
   template <typename C> //cast type
   using cast_iterator         = RigidBodyCastIterator<C>;
   template <typename C> //cast type
   using const_cast_iterator   = ConstRigidBodyCastIterator<C>;
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
   inline bool           isEmpty () const;
   inline size_type      size    () const;

   inline iterator       begin   ();
   inline const_iterator begin   () const;
   inline const_iterator cbegin  () const;
   inline iterator       end     ();
   inline const_iterator end     () const;
   inline const_iterator cend    () const;

   template< typename C >
   inline cast_iterator<C> begin();
   template< typename C >
   inline const_cast_iterator<C> begin() const;
   template< typename C >
   inline const_cast_iterator<C> cbegin() const;
   template< typename C >
   inline cast_iterator<C> end();
   template< typename C >
   inline const_cast_iterator<C> end() const;
   template< typename C >
   inline const_cast_iterator<C> cend() const;

   inline RigidBody&         front();
   inline const RigidBody&   front() const;
   inline RigidBody&         back();
   inline const RigidBody&   back() const;

   inline BodyID         at      ( size_type index );
   inline ConstBodyID    at      ( size_type index ) const;
   inline iterator       find    ( id_t sid );
   inline const_iterator find    ( id_t sid ) const;
   inline iterator       find    ( ConstBodyID body );
   inline const_iterator find    ( ConstBodyID body ) const;
   inline void           validate();
   //@}
   //**********************************************************************************************

   //**Add/Remove functions************************************************************************
   /*!\name Add/Remove functions */
   //@{
   [[deprecated("Please wrap the body into a unique_ptr")]]
   inline RigidBody&     add     ( BodyID body );
   inline RigidBody&     add     ( std::unique_ptr<RigidBody>&& body );
   inline iterator       remove  ( const id_t sid );
   inline iterator       remove  ( BodyID body );
   inline const_iterator remove  ( const_iterator pos );
   inline iterator       remove  ( iterator pos );
   inline std::unique_ptr<RigidBody> release ( const id_t sid );
   inline std::unique_ptr<RigidBody> release ( BodyID body );
   inline std::unique_ptr<RigidBody> release ( const_iterator& pos );
   inline std::unique_ptr<RigidBody> release ( iterator& pos );
   inline void           clear   ();
   //@}
   //**********************************************************************************************

   //**Callbacks************************************************************************
   /*!\name Callbacks */
   //@{
   inline void          registerAddCallback     ( const std::string& name, const std::function<void (BodyID)>& func );
   inline void          deregisterAddCallback   ( const std::string& name );
   inline void          clearAddCallbacks       ( );

   inline void          registerRemoveCallback     ( const std::string& name, const std::function<void (BodyID)>& func );
   inline void          deregisterRemoveCallback   ( const std::string& name );
   inline void          clearRemoveCallbacks       ( );
   //@}
   //**********************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VectorContainer           bodies_;    //!< The rigid bodies contained in the simulation world.
   std::map<id_t, size_type> bodyIDs_;   //!< The association of system IDs to rigid bodies.

   std::map< std::string, std::function<void (BodyID)> > addCallbacks_;
   std::map< std::string, std::function<void (BodyID)> > removeCallbacks_;
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
   : bodyIDs_()
{
   bodies_.reserve(1000);
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
   return bodies_.empty();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of rigid bodies contained in the body storage.
 *
 * \return The number of rigid bodies.
 */

inline BodyStorage::size_type BodyStorage::size() const
{
   return bodies_.size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first contained rigid body.
 *
 * \return Iterator to the first contained rigid body.
 */

inline BodyStorage::iterator BodyStorage::begin()
{
   return BodyStorage::iterator(bodies_.begin());
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator to the first contained rigid body.
 *
 * \return Constant iterator to the first contained rigid body.
 */

inline BodyStorage::const_iterator BodyStorage::begin() const
{
   return cbegin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator to the first contained rigid body.
 *
 * \return Constant iterator to the first contained rigid body.
 */

inline BodyStorage::const_iterator BodyStorage::cbegin() const
{
   return BodyStorage::const_iterator(bodies_.cbegin());
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last contained rigid body.
 *
 * \return Iterator just past the last contained rigid body.
 */

inline BodyStorage::iterator BodyStorage::end()
{
   return BodyStorage::iterator(bodies_.end());
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator just past the last contained rigid body.
 *
 * \return Constant iterator just past the last contained rigid body.
 */

inline BodyStorage::const_iterator BodyStorage::end() const
{
   return cend();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator just past the last contained rigid body.
 *
 * \return Constant iterator just past the last contained rigid body.
 */

inline BodyStorage::const_iterator BodyStorage::cend() const
{
   return BodyStorage::const_iterator(bodies_.cend());
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first contained rigid body.
 *
 * \return Iterator to the first contained rigid body.
 */

template< typename C >  // Cast Type
inline BodyStorage::cast_iterator<C> BodyStorage::begin()
{
   return BodyStorage::cast_iterator<C>(bodies_.begin(), bodies_.end());
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator to the first contained rigid body.
 *
 * \return Constant iterator to the first contained rigid body.
 */

template< typename C >  // Cast Type
inline BodyStorage::const_cast_iterator<C> BodyStorage::begin() const
{
   return cbegin<C>();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator to the first contained rigid body.
 *
 * \return Constant iterator to the first contained rigid body.
 */

template< typename C >  // Cast Type
inline BodyStorage::const_cast_iterator<C> BodyStorage::cbegin() const
{
   return BodyStorage::const_cast_iterator<C>(bodies_.begin(), bodies_.end());
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last contained rigid body.
 *
 * \return Iterator just past the last contained rigid body.
 */

template< typename C >  // Cast Type
inline BodyStorage::cast_iterator<C> BodyStorage::end()
{
   return BodyStorage::cast_iterator<C>(bodies_.end(), bodies_.end());
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator just past the last contained rigid body.
 *
 * \return Constant iterator just past the last contained rigid body.
 */

template< typename C >  // Cast Type
inline BodyStorage::const_cast_iterator<C> BodyStorage::end() const
{
   return cend<C>();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a constant iterator just past the last contained rigid body.
 *
 * \return Constant iterator just past the last contained rigid body.
 */

template< typename C >  // Cast Type
inline BodyStorage::const_cast_iterator<C> BodyStorage::cend() const
{
   return BodyStorage::const_cast_iterator<C>(bodies_.end(), bodies_.end());
}
//*************************************************************************************************

inline RigidBody&         BodyStorage::front()
{
   return *bodies_.front();
}
inline const RigidBody&   BodyStorage::front() const
{
   return *bodies_.front();
}
inline RigidBody&         BodyStorage::back()
{
   return *bodies_.back();
}
inline const RigidBody&   BodyStorage::back() const
{
   return *bodies_.back();
}


//*************************************************************************************************
/*!\brief Returns a handle to the indexed rigid body.
 *
 * \param index Access index. The index has to be in the range \f$[0..size-1]\f$.
 * \return Handle to the accessed rigid body.
 *
 * \b Note: No runtime check is performed to ensure the validity of the access index.
 */

inline BodyID BodyStorage::at( size_type index )
{
   WALBERLA_ASSERT( index < size(), "Invalid body index" );
   return bodies_[index].get();
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

inline ConstBodyID BodyStorage::at( size_type index ) const
{
   WALBERLA_ASSERT( index < size(), "Invalid body index" );
   return bodies_[index].get();
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

inline BodyStorage::iterator BodyStorage::find( id_t sid )
{
   std::map<id_t, size_type>::const_iterator pos = bodyIDs_.find( sid );
   if( pos == bodyIDs_.end() )
      return BodyStorage::iterator(bodies_.end());

   return BodyStorage::iterator(bodies_.begin() + static_cast<VectorContainer::iterator::difference_type>(pos->second));
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

inline BodyStorage::const_iterator BodyStorage::find( id_t sid ) const
{
   std::map<id_t, size_type>::const_iterator pos = bodyIDs_.find( sid );
   if( pos == bodyIDs_.end() )
      return BodyStorage::const_iterator(bodies_.end());

   return BodyStorage::const_iterator(bodies_.begin() + static_cast<VectorContainer::iterator::difference_type>(pos->second));
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

inline BodyStorage::iterator BodyStorage::find( ConstBodyID body )
{
   return BodyStorage::iterator(find( body->getSystemID() ));
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

inline BodyStorage::const_iterator BodyStorage::find( ConstBodyID body ) const
{
   return BodyStorage::const_iterator(find( body->getSystemID() ));
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

inline RigidBody& BodyStorage::add( BodyID body )
{
   WALBERLA_ASSERT( bodyIDs_.find( body->getSystemID() ) == bodyIDs_.end(), "Body with same system ID already added." );
   bodyIDs_[ body->getSystemID() ] = bodies_.size();
   bodies_.push_back( std::unique_ptr<RigidBody>(body) );

   for (auto it = addCallbacks_.begin(); it != addCallbacks_.end(); ++it)
   {
      it->second(bodies_.back().get());
   }

   return *bodies_.back();
}
//*************************************************************************************************


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

inline RigidBody& BodyStorage::add( std::unique_ptr<RigidBody>&& body )
{
   WALBERLA_ASSERT( bodyIDs_.find( body->getSystemID() ) == bodyIDs_.end(), "Body with same system ID already added." );
   bodyIDs_[ body->getSystemID() ] = bodies_.size();
   bodies_.push_back( std::move(body) );

   for (auto it = addCallbacks_.begin(); it != addCallbacks_.end(); ++it)
   {
      it->second(bodies_.back().get());
   }

   return *bodies_.back();
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

inline
BodyStorage::iterator BodyStorage::remove( const id_t sid )
{
   std::map<id_t, size_type>::iterator it = bodyIDs_.find( sid );
   WALBERLA_ASSERT( it != bodyIDs_.end(), "The body's system ID was not registered." );

   // Move last element to deleted place and update the system ID to index mapping.
   size_type i = it->second;

   for (auto cb = removeCallbacks_.begin(); cb != removeCallbacks_.end(); ++cb)
   {
      cb->second( bodies_[i].get() );
   }

   bodyIDs_[ bodies_.back()->getSystemID() ] = i;
   std::swap( bodies_[i], bodies_.back() );
   bodyIDs_.erase( it );
   bodies_.pop_back();
   return iterator(bodies_.begin() + int_c(i));
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

inline BodyStorage::const_iterator BodyStorage::remove( const_iterator pos )
{
   return remove( pos->getSystemID() );
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

inline BodyStorage::iterator BodyStorage::remove( iterator pos )
{
   return remove( pos->getSystemID() );
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

inline
BodyStorage::iterator BodyStorage::remove( BodyID body )
{
   return remove( body->getSystemID() );
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

inline std::unique_ptr<RigidBody> BodyStorage::release( const id_t sid )
{
   std::map<id_t, size_type>::iterator it = bodyIDs_.find( sid );
   WALBERLA_ASSERT( it != bodyIDs_.end(), "The body's system ID was not registered." );

   // Move last element to deleted place and update the system ID to index mapping.
   size_type i = it->second;

   std::for_each(removeCallbacks_.begin(), removeCallbacks_.end(), [&](auto& cb){cb.second( bodies_[i].get() );});

   bodyIDs_[ bodies_.back()->getSystemID() ] = i;
   std::swap( bodies_[i], bodies_.back() );
   bodyIDs_.erase( it );
   auto tmp = std::move(bodies_.back());
   bodies_.pop_back();
   return tmp;
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

inline std::unique_ptr<RigidBody> BodyStorage::release( const_iterator& pos )
{
   auto diff = std::distance(bodies_.cbegin(), pos.get());
   auto tmp = release( pos->getSystemID() );
   pos = const_iterator(bodies_.begin() + diff);
   return tmp;
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

inline std::unique_ptr<RigidBody> BodyStorage::release( iterator& pos )
{
   auto diff = std::distance(bodies_.begin(), pos.get());
   auto tmp = release( pos->getSystemID() );
   pos = iterator(bodies_.begin() + diff);
   return tmp;
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

inline std::unique_ptr<RigidBody> BodyStorage::release( BodyID body )
{
   return release( body->getSystemID() );
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
         cb->second( bodyIt->get() );
      }
   }

   bodyIDs_.clear();
   bodies_.clear();
}
//*************************************************************************************************

inline void          BodyStorage::registerAddCallback     ( const std::string& name, const std::function<void (BodyID)>& func )
{
   WALBERLA_ASSERT_EQUAL(addCallbacks_.find(name), addCallbacks_.end(), "Callback '" << name << "' already exists!");
   addCallbacks_.insert( std::make_pair(name, func) );
}

inline void          BodyStorage::deregisterAddCallback   ( const std::string& name )
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

inline void          BodyStorage::registerRemoveCallback     ( const std::string& name, const std::function<void (BodyID)>& func )
{
   WALBERLA_ASSERT_EQUAL(removeCallbacks_.find(name), removeCallbacks_.end(), "Callback '" << name << "' already exists!");
   removeCallbacks_.insert( std::make_pair(name, func) );
}

inline void          BodyStorage::deregisterRemoveCallback   ( const std::string& name )
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
   std::map<id_t, size_type>::iterator it = bodyIDs_.begin();
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
