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
//! \file UniqueID.h
//! \ingroup core
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <core/mpi/MPIManager.h>
#include <core/debug/Debug.h>
#include <core/NonCreateable.h>
#include <core/DataTypes.h>


namespace walberla {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generation of unique IDs in MPI environments.
 * \ingroup core
 *
 * The UniqueID class is responsible for the generation of unique system IDs even during an MPI
 * parallel execution. The unique ID is composed from the number of locally generated objects
 * of type \a T and the rank of the local MPI process.
 */
template< typename T >
class UniqueID : private NonCreateable
{
public:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static inline id_t create();
   static inline id_t createGlobal();
   static inline id_t invalidID();
   //@}
   //**********************************************************************************************
   
   static inline id_t count() { return counter_; }
   static inline id_t globalCount() { return globalCounter_; }

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   static id_t counter_;        //!< System ID counter.
   static id_t globalCounter_;  //!< Global system ID counter.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  DEFINITION AND INITIALIZATION OF THE STATIC MEMBER VARIABLES
//
//=================================================================================================

template< typename T >
id_t UniqueID<T>::counter_( 1 );

template< typename T >
id_t UniqueID<T>::globalCounter_( 1 );




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generates a unique ID.
 *
 * \return The newly generated unique ID.
 */
template< typename T >
inline id_t UniqueID<T>::create()
{
   id_t id( counter_++ );

   if( MPIManager::instance()->bitsNeededToRepresentRank() > 0 ) {
      const size_t shift( sizeof(id_t)*size_t(8) - MPIManager::instance()->bitsNeededToRepresentRank() );
      WALBERLA_ASSERT( shift < sizeof(id_t)*size_t(8), "Invalid shift detected" );

      const int rank ( MPIManager::instance()->worldRank() );
      WALBERLA_ASSERT( rank >= 0, "MPI system is not initialized" );

      if( counter_ >= ( id_t(1) << shift ) ) {
         // Revert overflow of counter
         --counter_;
         throw std::runtime_error( "Unable to create unique ID" );
      }

      id |= ( static_cast<id_t>( rank ) << shift );
   }
   else if( counter_ == size_t(0) ) {
      // Revert overflow of counter
      --counter_;
      throw std::runtime_error( "Unable to create unique ID" );
   }

   return id;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Generates an ID which is the same on all processes but unique on each process.
 *
 * \return The newly generated ID.
 */
template< typename T >
inline id_t UniqueID<T>::createGlobal()
{
   id_t id( globalCounter_++ );

   if( MPIManager::instance()->bitsNeededToRepresentRank() > 0 ) {
      const size_t shift( sizeof(id_t)*size_t(8) - MPIManager::instance()->bitsNeededToRepresentRank() );
      WALBERLA_ASSERT( shift < sizeof(id_t)*size_t(8), "Invalid shift detected" );

      const int size( MPIManager::instance()->numProcesses() );
      WALBERLA_ASSERT( size >= 0, "MPI system is not initialized" );

      if( globalCounter_ >= ( id_t(1) << shift ) ) {
         // Revert overflow of counter
         --globalCounter_;
         throw std::runtime_error( "Unable to create unique ID" );
      }

      id |= ( static_cast<id_t>( size ) << shift );
   }
   else if( globalCounter_ == size_t(0) ) {
      // Revert overflow of counter
      --globalCounter_;
      throw std::runtime_error( "Unable to create unique ID" );
   }

   return id;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns an invalid ID
 *
 * \return invalid ID.
 */
template< typename T >
inline id_t UniqueID<T>::invalidID()
{
   return 0;
}
//*************************************************************************************************

}
