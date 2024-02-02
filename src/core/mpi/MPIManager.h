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
//! \file MPIManager.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Michael Zikeli <michael.zikeli@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/math/Uint.h"
#include "core/mpi/Datatype.h"
#include "core/mpi/MPIOperation.h"
#include "core/mpi/MPIWrapper.h"
#include "core/mpi/Operation.h"
#include "core/singleton/Singleton.h"

#include <map>
#include <typeindex>

namespace walberla {
namespace mpi {


/**
 * Encapsulates MPI Rank/Communicator information
 *
 * \ingroup mpi
 *
 * Every process has two ranks/communicators:
 *   World:  This communicator/rank is valid after calling activateMPI,
 *           usually at the beginning of the program.
 *           This communicator never changes.
 *
 *   Custom: Can be adapted to the block structure.
 *           During the block structure setup, either the Cartesian setup
 *           has to be chosen using createCartesianComm() or the
 *           world communicator has to be used: useWorldComm()
 *
 */
class MPIManager : public singleton::Singleton< MPIManager >
{
public:

   WALBERLA_BEFRIEND_SINGLETON;

   ~MPIManager();

   /**
    * Configures the class, initializes numProcesses, worldRank
    * the rank and comm variables are still invalid, until custom communicator is set up
    * @param abortOnException if true, MPI_Abort is called in case of an uncaught exception
    */
   void initializeMPI( int* argc, char*** argv, bool abortOnException = true );

   void finalizeMPI();

   void resetMPI();

   void abort();

   //** Cartesian Communicator *****************************************************************************************
   /*! \name Cartesian Communicator */
   //@{

   void createCartesianComm( int numberOfProcessors[3], int periodicity[3] );
   void createCartesianComm( const uint_t xProcesses,        const uint_t  yProcesses,        const uint_t zProcesses,
                             const bool   xPeriodic = false, const bool    yPeriodic = false, const bool   zPeriodic = false );

   /// Cartesian coordinates of own rank
   void cartesianCoord( int coordOut[3] ) const;
   /// Cartesian coordinates of given rank
   void cartesianCoord( int rank, int coordOut[3] ) const;

   /// translates Cartesian coordinates to rank
   int cartesianRank( int coords[3] ) const;
   /// translates Cartesian coordinates to rank
   int cartesianRank( const uint_t x, const uint_t y, const uint_t z ) const;

   //@}
   //*******************************************************************************************************************



   //** World Communicator *********************************************************************************************
   /*! \name World Communicator */
   //@{

   void useWorldComm() { WALBERLA_MPI_SECTION() { WALBERLA_ASSERT_EQUAL( rank_, -1 ); rank_ = worldRank_; comm_ = MPI_COMM_WORLD; } }

   //@}
   //*******************************************************************************************************************



   //** Getter Functions  **********************************************************************************************
   /*! \name Getter Function */
   //@{

   int worldRank()    const { WALBERLA_MPI_SECTION() { WALBERLA_ASSERT( isMPIInitialized_ ); } return worldRank_;    }
   int numProcesses() const { WALBERLA_MPI_SECTION() { WALBERLA_ASSERT( isMPIInitialized_ ); } return numProcesses_; }

   int      rank() const { WALBERLA_ASSERT_GREATER_EQUAL(rank_, 0); return rank_; }
   MPI_Comm comm() const { WALBERLA_ASSERT_GREATER_EQUAL(rank_, 0); return comm_; }

   uint_t      bitsNeededToRepresentRank() const { return math::uintMSBPosition( uint_c(numProcesses()) ); }

   bool isMPIInitialized()  const { return isMPIInitialized_; }
   bool hasCartesianSetup() const { return cartesianSetup_;  }
   /// Rank is valid after calling createCartesianComm() or useWorldComm()
   bool rankValid()         const { return rank_ >= 0;       }
   bool hasWorldCommSetup()  const { return rankValid() && !hasCartesianSetup();}

   /// Indicates whether MPI-IO can be used with the current MPI communicator; certain versions of OpenMPI produce
   /// segmentation faults when using MPI-IO with a 3D Cartesian MPI communicator (see waLBerla issue #73)
   bool isCommMPIIOValid() const;

   /// Return the custom MPI_Datatype stored in 'customMPITypes_' and defined by the user and passed to 'commitCustomType'.
   template< typename CType >
   MPI_Datatype getCustomType() const
   {
      WALBERLA_MPI_SECTION()
      {
         return customMPITypes_.at(typeid(CType)).operator MPI_Datatype();
      }
      WALBERLA_NON_MPI_SECTION()
      {
         WALBERLA_ABORT( "This should not be called, if waLBerla is compiled without MPI." );
      }
   }

   /// Return the custom MPI_Op stored in 'customMPIOperation_' and defined by the user and passed to 'commitCustomOperation'.
   template< typename CType >
   MPI_Op getCustomOperation(mpi::Operation op) const
   {
      // FIXME the operation is actually type dependent but implementing this is not straightforward,
      //    compare comment at declaration of 'customMPIOperations_'.
      WALBERLA_MPI_SECTION()
      {
         return customMPIOperations_.at(op).operator MPI_Op();
      }
      WALBERLA_NON_MPI_SECTION()
      {
         WALBERLA_ABORT( "This should not be called, if waLBerla is compiled without MPI." );
      }
   }
   //@}
   //*******************************************************************************************************************

   //** Setter Functions  **********************************************************************************************
   /*! \name Setter Function */
   //@{
   ///! \brief Initializes a custom MPI_Datatype and logs it in the customMPITypes_ map.
   ///! \param argument The argument that is expected by the constructor of mpi::Datatype
   ///     At the point of creation 26.01.2024 this is either MPI_Datatype or const int.
   template < typename CType, class ConstructorArgumentType >
   void commitCustomType( ConstructorArgumentType& argument )
   {
      WALBERLA_MPI_SECTION()
      {
         if (isMPIInitialized_ && !currentlyAborting_)
         {
            static_assert( std::is_same_v<ConstructorArgumentType, const int> || std::is_same_v<ConstructorArgumentType, MPI_Datatype>,
                           "mpi::Datatype has only an constructor for an int value or a MPI_Datatype." );
            [[maybe_unused]] auto worked = std::get< 1 >( customMPITypes_.try_emplace(typeid(CType), argument) );
            WALBERLA_ASSERT(worked, "Emplacement of type " << typeid(CType).name() << " did not work.");
         } else {
            WALBERLA_ABORT( "MPI must be initialized before an new MPI_Datatype can be committed." );
         }
      }
      WALBERLA_NON_MPI_SECTION()
      {
         WALBERLA_ABORT( "This should not be called, if waLBerla is compiled without MPI." );
      }
   }

   ///! \brief Initializes a custom MPI_Op and logs it in the customMPIOperation map
   ///! \param op  A operator, e.g. SUM, MIN.
   ///! \param fct The definition of the MPI_User_function used for this operator.
   template< typename CType >
   void commitCustomOperation( mpi::Operation op, MPI_User_function* fct )
   {
      WALBERLA_MPI_SECTION()
      {
         if (isMPIInitialized_ && !currentlyAborting_)
         {
            [[maybe_unused]] auto worked = std::get< 1 >(customMPIOperations_.try_emplace(op, fct));
            WALBERLA_ASSERT(worked, "Emplacement of operation " << typeid(op).name() << " of type "
                                                                << typeid(CType).name() << " did not work.");
         }
         else { WALBERLA_ABORT("MPI must be initialized before an new MPI_Op can be committed."); }
      }
      WALBERLA_NON_MPI_SECTION()
      {
         WALBERLA_ABORT( "This should not be called, if waLBerla is compiled without MPI." );
      }
   }
   //@}
   //*******************************************************************************************************************

   static std::string getMPIErrorString(int errorCode);
   static std::string getMPICommName(MPI_Comm comm);

private:

   /// Rank in MPI_COMM_WORLD
   int worldRank_{0};

   /// Rank in the custom communicator
   int rank_{-1};

   /// Total number of processes
   int numProcesses_{1};

   /// Use this communicator for all MPI calls
   /// this is in general not equal to MPI_COMM_WORLD
   /// this may change during domain setup, where a custom communicator adapted to the
   /// domain is created
   MPI_Comm comm_;

   /// Indicates whether initializeMPI has been called. If true, MPI_Finalize is called upon destruction
   bool isMPIInitialized_{false};

   /// Indicates whether a Cartesian communicator has been created
   bool cartesianSetup_{false};

   bool currentlyAborting_{false};

   bool finalizeOnDestruction_{false};

   // Singleton
   MPIManager() : comm_(MPI_COMM_NULL) { WALBERLA_NON_MPI_SECTION() { rank_ = 0; } }

/// It is possible to commit own datatypes to MPI, that are not part of the standard. One example would be float16.
/// With these maps, it is possible to track self defined MPI_Datatypes and MPI_Ops, to access them at any time and
/// place in the program, also, they are automatically freed once MPIManager::finalizeMPI is called.
/// To initialize types or operations and add them to the map, the getter functions 'commitCustomType' and
/// 'commitCustomOperation' should be used. This can for example be done e.g. in the specialization of the MPITrait of
/// the newly defined type. For an example see MPIWrapper.cpp
   std::map< std::type_index, walberla::mpi::Datatype > customMPITypes_{};
   std::map< walberla::mpi::Operation, walberla::mpi::MPIOperation > customMPIOperations_{};
   // FIXME this must be type specific as well, but doing so is a bit more complicated.
   //  1. Idea) Combining both maps together e.g. as std::map< typeid(CType),
   //                                                          std::pair< MPI_DataType,
   //                                                                     std::map< Operation,
   //                                                                               MPI_Op > > > customMPITypesWithOps_{};
   //           There the access is quite nasty to look at, but easily doable, the construction however is quite difficult
   //           also one needs to make sure that the type is initialized before the operation.
   //  2. Idea) Leaving it as two maps customMPITypes_ and customMPIOperations,
   //           but storing a pair of typeid and operation as key for the operation map.
   //           This way everything would look nice, but someone needs to implement a comparison operator for this pair.
   //           I personally don't know where to put this comparison operator to, since it should not be part of the manager.
   //  3. Idea) Since this relies on the use of MPITrait<CType> --> { MPI_Datatype, MPI_Op } someone could define a object
   //           to store in the MPIManager there, to keep the MPIManager light and easily understandable.
   //           I'm also not sure if the MPITrait is the right spot for this though.
   //  For more information about the changes done in the code to allow custom defined types and operations,
   //  check out MR !647 ( https://i10git.cs.fau.de/walberla/walberla/-/merge_requests/647 )



}; // class MPIManager





// ROOT SECTION //

#define     WALBERLA_ROOT_SECTION() if( ::walberla::MPIManager::instance()->worldRank() == 0 )
#define WALBERLA_NON_ROOT_SECTION() if( ::walberla::MPIManager::instance()->worldRank() != 0 )



// EXCLUSIVE SECTION //

/**
 * Exclusive section is only executed on the process with the given rank
 * \warning The rank here is with respect to MPIManager::comm() (not w.r.t. MPI_COMM_WORLD)
 *          These ranks are not valid until the block setup has been done.
 *          This means that WALBERLA_ROOT_SECTION() may be different than WALBERLA_EXCLUSIVE_SECTION(0) !!
 *          If an exclusive section is needed before the block setup has been done, the
 *          WALBERLA_EXCLUSIVE_WORLD_SECTION can be used.
 */
#define WALBERLA_EXCLUSIVE_SECTION(exclusive_rank) if( ::walberla::MPIManager::instance()->rank() == (exclusive_rank) )

/**
 * Exclusive section using ranks w.r.t. MPI_COMM_WORLD.
 */
#define WALBERLA_EXCLUSIVE_WORLD_SECTION(exclusive_rank) if( ::walberla::MPIManager::instance()->worldRank() == (exclusive_rank) )





// BARRIERS //

#ifdef WALBERLA_BUILD_WITH_MPI

#define WALBERLA_MPI_WORLD_BARRIER() { MPI_Barrier( MPI_COMM_WORLD ); }
#define WALBERLA_MPI_BARRIER()       { MPI_Barrier( ::walberla::MPIManager::instance()->comm() ); }


#define WALBERLA_CRITICAL_SECTION_START for( int __WALBERLA_r__ = 0; __WALBERLA_r__ <  ::walberla::MPIManager::instance()->numProcesses(); ++__WALBERLA_r__ ) \
                                        { if ( __WALBERLA_r__ == ::walberla::MPIManager::instance()->worldRank() ) {


#define WALBERLA_CRITICAL_SECTION_END } MPI_Barrier( MPI_COMM_WORLD ); }


#else

#define WALBERLA_MPI_WORLD_BARRIER() {}
#define WALBERLA_MPI_BARRIER()       {}

#define WALBERLA_CRITICAL_SECTION_START
#define WALBERLA_CRITICAL_SECTION_END


#endif


} // namespace mpi
} // namespace walberla





//======================================================================================================================
//
//  Export to walberla namespace
//
//======================================================================================================================

namespace walberla {
   using mpi::MPIManager;
}


