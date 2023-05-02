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
//! \file LubricationCorrectionMEM.cpp
//! \ingroup pe_coupling
//! \author Kristina Pickl <kristina.pickl@fau.de>
//! \author Dominik Bartuschat
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/mpi/Reduce.h"
#include "core/timing/RemainingTimeLogger.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SplitPureSweep.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "pe/basic.h"
#include "pe/utility/Distance.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include <functional>


namespace lubrication_correction_mem
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;


//////////////
// TYPEDEFS //
//////////////

// pdf field & flag field
typedef lbm::D3Q19< lbm::collision_model::TRT >  LatticeModel_T;
using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

const uint_t FieldGhostLayers = 1;

// pe body ID field

typedef GhostLayerField< pe::BodyID, 1 >  BodyField_T;

// boundary handling

typedef lbm::FreeSlip< LatticeModel_T, FlagField_T>           FreeSlip_T;
typedef pe_coupling::SimpleBB< LatticeModel_T, FlagField_T > MO_T;

typedef BoundaryHandling< FlagField_T, Stencil_T, FreeSlip_T, MO_T > BoundaryHandling_T;

typedef std::tuple<pe::Sphere, pe::Plane> BodyTypeTuple ;

///////////
// FLAGS //
///////////

const FlagUID    Fluid_Flag( "fluid" );
const FlagUID FreeSlip_Flag( "free slip" );
const FlagUID       MO_Flag( "moving obstacle" );
const FlagUID FormerMO_Flag( "former moving obstacle" );


//////////////////////
// HELPER FUNCTIONS //
//////////////////////


//************************************************************************************************************************************
/*! \brief
 *
 * Includes output of previously set force (lubrication correction + drag force).
 * Compares the force to an analytical calculation within a tolerance of 1%.
 */
//************************************************************************************************************************************

class EvaluateLubricationForce
{
public:
   EvaluateLubricationForce( const shared_ptr< StructuredBlockStorage > & blocks,
                             const BlockDataID & bodyStorageID,
                             const shared_ptr<pe::BodyStorage> & globalBodyStorage,
                             uint_t id1, uint_t id2, pe::Vec3 vel, real_t nu_L, real_t radius,
                             SweepTimeloop* timeloop, bool print, bool sphSphTest, bool sphWallTest )
   : blocks_( blocks ), bodyStorageID_( bodyStorageID ), globalBodyStorage_( globalBodyStorage ),
     id1_( id1 ), id2_( id2 ), vel_( vel ), nu_L_( nu_L ), radius_( radius ),
     timeloop_( timeloop ), print_( print ), sphSphTest_( sphSphTest ), sphWallTest_( sphWallTest ) {}

   void operator()()
   {
      if( sphSphTest_ )
      {
         checkSphSphLubricationForce();
      } else if ( sphWallTest_ )
      {
         checkSphWallLubricationForce();
      }
   }

private:

   void checkSphSphLubricationForce()
   {
      const uint_t timestep (timeloop_->getCurrentTimeStep()+1);

      // variables for output of total force - requires MPI-reduction
      pe::Vec3 forceSphr1(0);
      pe::Vec3 forceSphr2(0);

      // temporary variables
      real_t gap        (0);

      // calculate gap between the two spheres
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto curSphereIt = pe::BodyIterator::begin<pe::Sphere>( *blockIt, bodyStorageID_); curSphereIt != pe::BodyIterator::end<pe::Sphere>(); ++curSphereIt )
         {
            pe::SphereID sphereI = ( curSphereIt.getBodyID() );
            if ( sphereI->getID() == id1_ )
            {
               for( auto blockIt2 = blocks_->begin(); blockIt2 != blocks_->end(); ++blockIt2 )
               {
                  for( auto oSphereIt = pe::BodyIterator::begin<pe::Sphere>( *blockIt2, bodyStorageID_); oSphereIt != pe::BodyIterator::end<pe::Sphere>(); ++oSphereIt )
                  {
                     pe::SphereID sphereJ = ( oSphereIt.getBodyID() );
                     if ( sphereJ->getID() == id2_ )
                     {
                        gap = pe::getSurfaceDistance( sphereI, sphereJ );
                        break;
                     }
                  }
               }
               break;
            }
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::reduceInplace( gap, mpi::MAX );
      }
      WALBERLA_ROOT_SECTION()
      {
         if (gap < real_comparison::Epsilon<real_t>::value )
         {
            // shadow copies have not been synced yet as the spheres are outside the overlap region
            gap = walberla::math::Limits<real_t>::inf();
            WALBERLA_LOG_INFO_ON_ROOT( "Spheres still too far apart to calculate gap!" );
         } else {
            WALBERLA_LOG_INFO_ON_ROOT( "Gap between sphere 1 and 2: " << gap );
         }
      }


      // get force on spheres
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin<pe::Sphere>( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            if( bodyIt->getID() == id1_ )
            {
               forceSphr1 += bodyIt->getForce();
            } else if( bodyIt->getID() == id2_ )
            {
               forceSphr2 += bodyIt->getForce();
            }
         }
      }

      // MPI reduction of pe forces over all processes
      WALBERLA_MPI_SECTION()
      {
         mpi::reduceInplace( forceSphr1, mpi::SUM );
         mpi::reduceInplace( forceSphr2, mpi::SUM );
      }
      WALBERLA_LOG_INFO_ON_ROOT("Total force on sphere " << id1_ << " : " << forceSphr1);
      WALBERLA_LOG_INFO_ON_ROOT("Total force on sphere " << id2_ << " : " << forceSphr2);

      if ( print_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file1;
            std::string filename1("Gap_LubricationForceBody1.txt");
            file1.open( filename1.c_str(), std::ofstream::app );
            file1.setf( std::ios::unitbuf );
            file1.precision(15);
            file1 << gap << " " << forceSphr1[0] << std::endl;
            file1.close();

            std::ofstream file2;
            std::string filename2("Gap_LubricationForceBody2.txt");
            file2.open( filename2.c_str(), std::ofstream::app );
            file2.setf( std::ios::unitbuf );
            file2.precision(15);
            file2 << gap << " " << forceSphr2[0] << std::endl;
            file2.close();
         }
      }


      WALBERLA_ROOT_SECTION()
      {
         if ( timestep == uint_t(1000) )
         {
            // both force x-components should be the same only with inverted signs
            WALBERLA_CHECK_FLOAT_EQUAL ( forceSphr2[0], -forceSphr1[0] );

            // according to the formula from Ding & Aidun 2003
            // F = 3/2 * math::pi * rho_L * nu_L * relative velocity of both spheres * r * r * 1/gap
            // the correct analytically calculated value is 339.2920063998
            // in this geometry setup the relative error is 0.1246489711 %
            real_t analytical = real_c(3.0)/real_c(2.0) * walberla::math::pi * real_c(1.0) * nu_L_ * real_c(2.0) * real_c(vel_[0]) * radius_ * radius_ * real_c(1.0)/gap;
            real_t relErr     = std::fabs( analytical - forceSphr2[0] ) / analytical * real_c(100.0);
            WALBERLA_CHECK_LESS( relErr, real_t(1) );
         }
      }
   }

   void checkSphWallLubricationForce()
   {
      const uint_t timestep (timeloop_->getCurrentTimeStep()+1);

      // variables for output of total force - requires MPI-reduction
      pe::Vec3 forceSphr1(0);

      // temporary variables
      real_t gap        (0);

      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto curSphereIt = pe::BodyIterator::begin<pe::Sphere>( *blockIt, bodyStorageID_); curSphereIt != pe::BodyIterator::end<pe::Sphere>(); ++curSphereIt )
         {
            pe::SphereID sphereI = ( curSphereIt.getBodyID() );
            if ( sphereI->getID() == id1_ )
            {
               for( auto globalBodyIt = globalBodyStorage_->begin(); globalBodyIt != globalBodyStorage_->end(); ++globalBodyIt)
               {
                  if( globalBodyIt->getID() == id2_ )
                  {
                     pe::PlaneID planeJ = static_cast<pe::PlaneID>( globalBodyIt.getBodyID() );
                     gap = pe::getSurfaceDistance(sphereI, planeJ);
                     break;
                  }
               }
               break;
            }
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::reduceInplace( gap, mpi::MAX );
      }
      WALBERLA_ROOT_SECTION()
      {
         if (gap < real_comparison::Epsilon<real_t>::value )
         {
            gap = walberla::math::Limits<real_t>::inf();
            WALBERLA_LOG_INFO_ON_ROOT( "Sphere still too far from wall to calculate gap!" );
         } else {
            WALBERLA_LOG_INFO_ON_ROOT( "Gap between sphere and wall: " << gap );
         }
      }

      // get force on sphere
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin<pe::Sphere>( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            if( bodyIt->getID() == id1_ )
            {
               forceSphr1 += bodyIt->getForce();
            }
         }
      }

      // MPI reduction of pe forces over all processes
      WALBERLA_MPI_SECTION()
      {
         mpi::reduceInplace( forceSphr1, mpi::SUM );
      }
      WALBERLA_LOG_INFO_ON_ROOT("Total force on sphere " << id1_ << " : " << forceSphr1);

      if ( print_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file1;
            std::string filename1("Gap_LubricationForceBody1.txt");
            file1.open( filename1.c_str(), std::ofstream::app );
            file1.setf( std::ios::unitbuf );
            file1.precision(15);
            file1 << gap << " " << forceSphr1[0] << std::endl;
            file1.close();
         }
      }

      WALBERLA_ROOT_SECTION()
      {
         if ( timestep == uint_t(26399) )
         {
            // according to the formula from Ding & Aidun 2003
            // F = 6 * math::pi * rho_L * nu_L * relative velocity of both bodies=relative velocity of the sphere * r * r * 1/gap
            // the correct analytically calculated value is 339.292006996217
            // in this geometry setup the relative error is 0.183515322065561 %
            real_t analytical = real_c(6.0) * walberla::math::pi * real_c(1.0) * nu_L_ * real_c(-vel_[0]) * radius_ * radius_ * real_c(1.0)/gap;
            real_t relErr     = std::fabs( analytical - forceSphr1[0] ) / analytical * real_c(100.0);
            WALBERLA_CHECK_LESS( relErr, real_t(1) );
         }
      }

   }

   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID bodyStorageID_;
   shared_ptr<pe::BodyStorage> globalBodyStorage_;

   uint_t         id1_;
   uint_t         id2_;
   pe::Vec3       vel_;
   real_t         nu_L_;
   real_t         radius_;
   SweepTimeloop* timeloop_;
   bool           print_;
   bool           sphSphTest_;
   bool           sphWallTest_;

};
//************************************************************************************************************************************


class ForceReset
{
   public:
   ForceReset( const shared_ptr< StructuredBlockStorage > & blocks,
               const BlockDataID & bodyStorageID,
               uint_t id1, uint_t id2)
      : blocks_( blocks ), bodyStorageID_( bodyStorageID ),
        id1_( id1 ),  id2_ (id2)
      { }

      void operator()()
      {
         for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
         {
            for( auto bodyIt = pe::BodyIterator::begin<pe::Sphere>( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end<pe::Sphere>(); ++bodyIt )
            {
               if( bodyIt->getID() == id1_ )
               {
                  bodyIt->resetForceAndTorque();
               } else if( bodyIt->getID() == id2_ )
               {
                  bodyIt->resetForceAndTorque();
               }
            }
         }
      }
private:
      shared_ptr< StructuredBlockStorage > blocks_;
      const BlockDataID bodyStorageID_;
      uint_t         id1_;
      uint_t         id2_;
};



/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////


//************************************************************************************************************************************
// class SphSphTestBoundaryHandling is responsible for creating the boundary handling and setting up the geometry at the outer
// domain boundaries for the sphere-sphere lubrication test
//************************************************************************************************************************************
class SphSphTestBoundaryHandling
{
public:

      SphSphTestBoundaryHandling( const BlockDataID & flagField, const BlockDataID & pdfField, const BlockDataID & bodyField ) :
         flagField_( flagField ), pdfField_( pdfField ), bodyField_ (bodyField) {}

      BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

      const BlockDataID flagField_;
      const BlockDataID pdfField_;
      const BlockDataID bodyField_;

}; // class SphSphTestBoundaryHandling



BoundaryHandling_T * SphSphTestBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField = block->getData< FlagField_T >( flagField_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfField_ );
   BodyField_T * bodyField = block->getData< BodyField_T >( bodyField_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "cf boundary handling", flagField, fluid,
                                                           FreeSlip_T( "FreeSlip", FreeSlip_Flag, pdfField, flagField, fluid ),
                                                           MO_T( "MO", MO_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ) );

   const auto freeslip = flagField->getFlag( FreeSlip_Flag );

   CellInterval domainBB = storage->getDomainCellBB();
   domainBB.xMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.xMax() += cell_idx_c( FieldGhostLayers );

   // WEST
//   CellInterval west( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMin(), domainBB.yMax(), domainBB.zMax() );
//   storage->transformGlobalToBlockLocalCellInterval( west, *block );

   // EAST
//   CellInterval east( domainBB.xMax(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
//   storage->transformGlobalToBlockLocalCellInterval( east, *block );

   domainBB.yMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.yMax() += cell_idx_c( FieldGhostLayers );

   // SOUTH
   CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( south, *block );
   handling->forceBoundary( freeslip, south );

   // NORTH
   CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( north, *block );
   handling->forceBoundary( freeslip, north );

   domainBB.zMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.zMax() += cell_idx_c( FieldGhostLayers );

   // BOTTOM
   CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
   storage->transformGlobalToBlockLocalCellInterval( bottom, *block );
   handling->forceBoundary( freeslip, bottom );

   // TOP
   CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( top, *block );
   handling->forceBoundary( freeslip, top );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}
//************************************************************************************************************************************




//************************************************************************************************************************************
// class SphWallTestBoundaryHandling is responsible for creating the boundary handling and setting up the geometry at the outer domain boundaries
// for the sphere-wall lubrication test
//************************************************************************************************************************************
class SphWallTestBoundaryHandling
{
public:

      SphWallTestBoundaryHandling( const BlockDataID & flagField, const BlockDataID & pdfField, const BlockDataID & bodyField ) :
         flagField_( flagField ), pdfField_( pdfField ), bodyField_ (bodyField) {}

      BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

      const BlockDataID flagField_;
      const BlockDataID pdfField_;
      const BlockDataID bodyField_;

}; // class SphWallTestBoundaryHandling



BoundaryHandling_T * SphWallTestBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField = block->getData< FlagField_T >( flagField_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfField_ );
   BodyField_T * bodyField = block->getData< BodyField_T >( bodyField_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "cf boundary handling", flagField, fluid,
                                                           FreeSlip_T( "FreeSlip", FreeSlip_Flag, pdfField, flagField, fluid ),
                                                           MO_T( "MO", MO_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ) );

   const auto freeslip = flagField->getFlag( FreeSlip_Flag );

   CellInterval domainBB = storage->getDomainCellBB();
   domainBB.xMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.xMax() += cell_idx_c( FieldGhostLayers );

   // WEST
   CellInterval west( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMin(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( west, *block );
   handling->forceBoundary( freeslip, west );

   // EAST
   CellInterval east( domainBB.xMax(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( east, *block );
   handling->forceBoundary( freeslip, east );

   domainBB.yMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.yMax() += cell_idx_c( FieldGhostLayers );

   // SOUTH
   CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( south, *block );
   handling->forceBoundary( freeslip, south );

   // NORTH
   CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( north, *block );
   handling->forceBoundary( freeslip, north );

   domainBB.zMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.zMax() += cell_idx_c( FieldGhostLayers );

   // BOTTOM
   CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
   storage->transformGlobalToBlockLocalCellInterval( bottom, *block );
   handling->forceBoundary( freeslip, bottom );

   // TOP
   CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( top, *block );
   handling->forceBoundary( freeslip, top );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}
//************************************************************************************************************************************




//////////
// MAIN //
//////////


//************************************************************************************************************************************
/* \brief Lubrication test for testing the lubrication force between two spheres or between a sphere and a plane.
 *
 * The sphere-sphere lubrication test is executed with the command line argument "--sphSphTest". It has the following domain setup.
 *
 *      ______________________________
 *    P|                              |
 *    E|     ID 1           ID 2      |P
 *    R|      ___            ___      |E
 *    I|     /   \ +u    -u /   \     |R
 *    O|    |  x  |->    <-|  x  |    |I
 *    D|     \___/          \___/     |O
 *    I|                              |D
 *    C|                              |I
 *     |______________________________|C
 *
 *  All other walls are set to freeslip boundaries. The spheres are approaching each other driven by the velocities +/- u.
 *  The domain is periodic in x-direction
 *
 *
 * The sphere-wall lubrication test is executed with the command line argument "--sphWallTest". It has the following domain setup.
 *
 *      ______________________________
 *     |                              |
 *     |       ID 1                   |
 *     |        ___                   |
 *     |       /   \                  |
 *     |    <-|  x  |                 |
 *     |       \___/                  |
 *     |                              |
 *     |                              |
 *     |______________________________|
 *
 * All walls are set to freeslip boundaries. The sphere is approaching the wall driven by the velocity -u.
 *
 * If you want to print out additional gap vs. lubrication + hydrodynamic force plots you have to additionally execute the test with
 * the command line argument "--fileIO".
 *
 */
//************************************************************************************************************************************

int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   auto processes = MPIManager::instance()->numProcesses();
   if( processes != 3 && processes != 6)
   {
      std::cerr << "number of processes must be equal to 3 or 6 !" << std::endl;
      return EXIT_FAILURE;
   }

   bool useFZYX        = false;
   bool useSplitKernel = false;
   bool fullPDFSync    = false;
   bool funcTest       = false;
   bool sphSphTest     = false;
   bool sphWallTest    = false;
   bool fileIO         = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--fzyx" )        == 0 ) { useFZYX        = true; continue; }
      if( std::strcmp( argv[i], "--split" )       == 0 ) { useSplitKernel = true; continue; }
      if( std::strcmp( argv[i], "--full" )        == 0 ) { fullPDFSync    = true; continue; }
      if( std::strcmp( argv[i], "--funcTest" )    == 0 ) { funcTest       = true; continue; }
      if( std::strcmp( argv[i], "--sphSphTest" )  == 0 ) { sphSphTest     = true; continue; }
      if( std::strcmp( argv[i], "--sphWallTest" ) == 0 ) { sphWallTest    = true; continue; }
      if( std::strcmp( argv[i], "--fileIO" )      == 0 ) { fileIO         = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   if( ! ( sphSphTest || sphWallTest || funcTest ) )
   {
      WALBERLA_ABORT("Specify either --sphSphTest, --sphWallTest, or --funcTest !");
   }

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   // parameters not equal for all test cases
   uint_t timesteps (0);     // total time steps of the simulation
   real_t nu_L      (0);     // lattice kinematic viscosity
   real_t dt_SI     (0);     // time step [s]
   bool periodicX   (false); // flag for the periodicity in x-direction

   // specification of the spheres
   real_t   radius   (0);
   pe::Vec3 velocity (0);
   uint_t   id1      (0);
   uint_t   id2      (0);

   // for sphere-sphere test
   if ( sphSphTest )
   {
      timesteps = uint_c(1000);
      nu_L      = real_t(2);
      dt_SI     = real_c(0.3);
      periodicX = true;

      radius    = real_t(6);
      velocity  = pe::Vec3( real_c(0.008), real_t(0), real_t(0) );
      id1       = uint_c(1);
      id2       = uint_c(2);
   }

   // for sphere-wall test ( other possible timestep/kinematic viscosity combinations are included )
   if ( sphWallTest )
   {
      timesteps = uint_c(26399);       // uint_c(13199);       // uint_c(17599);       // uint_c(19800);
      nu_L      = real_t(1)/real_t(4); // real_t(1)/real_t(8); // real_t(1)/real_t(6); // real_t(3)/real_t(16);
      dt_SI     = real_c(0.125);
      periodicX = false;
      radius    = real_t(12);
      velocity  = pe::Vec3 (real_c(-0.001),real_t(0),real_t(0) );
      id1       = uint_c(1);
      id2       = uint_c(56);
   }

   uint_t xBlocks = ( processes == 3 ) ? uint_c(3) : uint_c(6); // number of blocks in x-direction
   uint_t yBlocks = uint_c(1);                                  // number of blocks in y-direction
   uint_t zBlocks = uint_c(1);                                  // number of blocks in z-direction

   // parameters equal for all test cases
   real_t rho_SI = real_c(1000);  // rho [kg/m^3]
   real_t dx_SI  = real_c(1e-3);  // dx [m]
   real_t dx     = real_t(1);     // lattice dx

   uint_t length = uint_c(192);   // length of the domain in x-direction in cells
   uint_t width  = uint_c(128);   // width (and height) of the domain in y- and z-direction in cells

   // only perform one single time-step, set the rest of the parameters to the sphere-sphere test
   // switch to a smaller domain and radius for the funcTest
   if ( funcTest )
   {
      timesteps = uint_c(1);
      length    = uint_c(24);
      width     = uint_c(12);
      nu_L      = real_t(2);
      dt_SI     = real_c(0.3);
      periodicX = true;

      radius    = real_t(2);
      velocity  = pe::Vec3 ( real_c(0.0001),real_t(0),real_t(0) );
      id1       = uint_c(1);
      id2       = uint_c(2);
   }

   uint_t xCells = length / xBlocks; // number of cells in x-direction on each block
   uint_t yCells =  width / yBlocks; // number of cells in y-direction on each block
   uint_t zCells =  width / zBlocks; // number of cells in z-direction on each block

   // Perform missing variable calculations
   real_t nu_SI  = dx_SI * dx_SI / dt_SI * nu_L;                    // kinematic viscosity [m^2/s]
   real_t tau    = real_c(0.5) * ( real_t(6) * nu_L + real_t(1) );
   real_t omega  = real_t(1) / tau;


   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   auto blocks = blockforest::createUniformBlockGrid( xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, dx, true,
                                                      periodicX, false, false );

   /////////////////
   // PE COUPLING //
   /////////////////

   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();

   auto bodyStorageID  = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");
   auto ccdID          = blocks->addBlockData(pe::ccd::createHashGridsDataHandling( globalBodyStorage, bodyStorageID ), "CCD");
   auto fcdID          = blocks->addBlockData(pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");

   // set up collision response, here DEM solver
   pe::cr::DEM cr(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID, nullptr);

   // set up synchronization procedure
   const real_t overlap = real_c( 1.5 ) * dx;
   std::function<void(void)> syncCall = std::bind( pe::syncShadowOwners<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );

   // create the material
   const auto myMat = pe::createMaterial( "myMat", real_c(1.4), real_t(0), real_t(1), real_t(1), real_t(0), real_t(1), real_t(1), real_t(0), real_t(0) );

   // sphere-sphere test
   if ( sphSphTest )
   {
      // add four planes
      pe::createPlane( *globalBodyStorage, walberla::id_t(52), Vector3<real_t>(0, 1,0), Vector3<real_t>(0,0,0), myMat ); // front
      pe::createPlane( *globalBodyStorage, walberla::id_t(53), Vector3<real_t>(0,-1,0), Vector3<real_t>(0,real_c(width),0), myMat ); // back
      pe::createPlane( *globalBodyStorage, walberla::id_t(54), Vector3<real_t>(0,0, 1), Vector3<real_t>(0,0,0), myMat ); // bottom
      pe::createPlane( *globalBodyStorage, walberla::id_t(55), Vector3<real_t>(0,0,-1), Vector3<real_t>(0,0,real_c(width)), myMat ); // top

      // create two approaching spheres
      auto sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, id1, Vector3<real_t> (real_c(50.0), real_c(64.0), real_c(64.0)), radius, myMat );
      if ( sphere != nullptr ) sphere->setLinearVel( velocity);
      sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, id2, Vector3<real_t> (real_c(78.0), real_c(64.0), real_c(64.0)), radius, myMat );
      if ( sphere != nullptr ) sphere->setLinearVel(-velocity);
   }

   // sphere-wall test
   if ( sphWallTest )
   {
      // create six planes
      pe::createPlane( *globalBodyStorage, walberla::id_t(52), Vector3<real_t>(0, 1,0), Vector3<real_t>(0,0,0), myMat ); // front
      pe::createPlane( *globalBodyStorage, walberla::id_t(53), Vector3<real_t>(0,-1,0), Vector3<real_t>(0,real_c(width),0), myMat ); // back
      pe::createPlane( *globalBodyStorage, walberla::id_t(54), Vector3<real_t>(0,0, 1), Vector3<real_t>(0,0,0), myMat ); // bottom
      pe::createPlane( *globalBodyStorage, walberla::id_t(55), Vector3<real_t>(0,0,-1), Vector3<real_t>(0,0,real_c(width)), myMat ); // top
      pe::createPlane( *globalBodyStorage, id2, Vector3<real_t>( 1,0,0), Vector3<real_t>(0,0,0), myMat ); // left
      pe::createPlane( *globalBodyStorage, walberla::id_t(57), Vector3<real_t>(-1,0,0), Vector3<real_t>(real_c(length),0,0), myMat ); // right

      // create one sphere ( other possible coordinates and radii are included)
      // <19.2,64,64> // <25.6,64,64> // <28.8,64,64>  // <38.4,64,64>    //1st: <41.765,64,64> // <27.205,64,64> // <20.88,64,64> (chosen s.th. last gap: 0.05 and initial s=1.1 ) // like Ding: <17,64,64>
      //  6           // 8            // 9             // 12              //1st  // 13.05       // 8.5            // 6.525                                                          // like Ding: 4.25
      auto sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, id1, Vector3<real_t> (real_c(38.4), real_c(64.0), real_c(64.0)), radius, myMat );
      if ( sphere != nullptr ) sphere->setLinearVel( velocity);
   }

   if ( funcTest )
   {
      // add four planes
      pe::createPlane( *globalBodyStorage, walberla::id_t(52), Vector3<real_t>(0, 1,0), Vector3<real_t>(0,0,0), myMat ); // front
      pe::createPlane( *globalBodyStorage, walberla::id_t(53), Vector3<real_t>(0,-1,0), Vector3<real_t>(0,real_c(width),0), myMat ); // back
      pe::createPlane( *globalBodyStorage, walberla::id_t(54), Vector3<real_t>(0,0, 1), Vector3<real_t>(0,0,0), myMat ); // bottom
      pe::createPlane( *globalBodyStorage, walberla::id_t(55), Vector3<real_t>(0,0,-1), Vector3<real_t>(0,0,real_c(width)), myMat ); // top

      // create two approaching spheres
      auto sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, id1, Vector3<real_t> (real_c( 6.0) , real_c(6.0), real_c(6.0)), radius, myMat );
      if ( sphere != nullptr ) sphere->setLinearVel( velocity);
      sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, id2, Vector3<real_t> (real_c(12.0) , real_c(6.0), real_c(6.0)), radius, myMat );
      if ( sphere != nullptr ) sphere->setLinearVel( -velocity);

      // set the same boundary handling as for the sphere-sphereTest
      sphSphTest = true;
   }

   syncCall();

   WALBERLA_LOG_INFO_ON_ROOT_SECTION()
   {
      std::stringstream ss;

      if ( sphSphTest )
      {
         ss << "-------------------------------------------------------\n"
            << "   Parameters for the sphere-sphere lubrication test \n"
            << "-------------------------------------------------------\n";
      }
      if ( sphWallTest )
      {
         ss << "-------------------------------------------------------\n"
            << "   Parameters for the sphere-wall lubrication test \n"
            << "-------------------------------------------------------\n";
      }
      ss << " tau          = " << tau      << "\n"
         << " omega        = " << omega    << "\n"
         << " nu_L         = " << nu_L     << "\n"
         << " nu_SI        = " << nu_SI    << "\n"
         << " rho_SI       = " << rho_SI   << "\n"
         << " dx_SI        = " << dx_SI    << "\n"
         << " dt_SI        = " << dt_SI    << "\n"
         << " radius       = " << radius   << "\n"
         << " velocity     = " << velocity << "\n"
         << "-------------------------------------------------------\n"
         << " length       = " << length   << "\n"
         << " width/height = " << width    << "\n"
         << " xBlocks      = " << xBlocks  << "\n"
         << " yBlocks      = " << yBlocks  << "\n"
         << " zBlocks      = " << zBlocks  << "\n"
         << " xCells       = " << xCells   << "\n"
         << " yCells       = " << yCells   << "\n"
         << " zCells       = " << zCells   << "\n"
         << "-------------------------------------------------------\n";
      WALBERLA_LOG_INFO( ss.str() );
   }

   ////////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // add pdf field

   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );

   BlockDataID pdfFieldID = useFZYX ? lbm::addPdfFieldToStorage( blocks, "pdf field (fzyx)", latticeModel,
                                                                 Vector3< real_t >( real_t(0), real_t(0), real_t(0) ), real_t(1),
                                                                 uint_t(1), field::fzyx ) :
                                      lbm::addPdfFieldToStorage( blocks, "pdf field (fzyx)", latticeModel,
                                                                 Vector3< real_t >( real_t(0), real_t(0), real_t(0) ), real_t(1),
                                                                 uint_t(1), field::fzyx );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::fzyx );

   // add boundary handling & initialize outer domain boundaries (moving walls on the front, back, top, and bottom plane)
   BlockDataID boundaryHandlingID;
   if ( sphSphTest )
   {
      boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
               SphSphTestBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID ), "boundary handling" );
   }
   if ( sphWallTest )
   {
      boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
               SphWallTestBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID ), "boundary handling" );
   }

   // map pe bodies into the LBM simulation
   pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag );

   ///////////////
   // TIME LOOP //
   ///////////////

   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // sweep for updating the pe body mapping into the LBM simulation
   timeloop.add()
      << Sweep( pe_coupling::BodyMapping< LatticeModel_T, BoundaryHandling_T >( blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, MO_Flag, FormerMO_Flag, pe_coupling::selectRegularBodies ), "Body Mapping" );

   // sweep for restoring PDFs in cells previously occupied by pe bodies
   typedef pe_coupling::EquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T > Reconstructor_T;
   Reconstructor_T reconstructor( blocks, boundaryHandlingID, bodyFieldID );
   timeloop.add()
      << Sweep( pe_coupling::PDFReconstruction< LatticeModel_T, BoundaryHandling_T, Reconstructor_T >( blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID,
                                                                                                       reconstructor, FormerMO_Flag, Fluid_Flag ), "PDF Restore" );
   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   std::function< void () > commFunction;
   if( fullPDFSync )
   {
      blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > scheme( blocks );
      scheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) );
      commFunction = scheme;
   }
   else
   {
      blockforest::communication::UniformBufferedScheme< Stencil_T > scheme( blocks );
      scheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) );
      commFunction = scheme;
   }

   // add LBM communication function and boundary handling sweep
   timeloop.add()
      << BeforeFunction( commFunction, "LBM Communication" )
      << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

   // LBM stream collide sweep
   if( useSplitKernel )
      timeloop.add() << Sweep( lbm::SplitPureSweep< LatticeModel_T >( pdfFieldID ), "LBM SPLIT PURE" );
   else
      timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag ) ), "LBM DEFAULT" );

   // add lubrication force correction sweep
   timeloop.addFuncAfterTimeStep( pe_coupling::LubricationCorrection( blocks, globalBodyStorage, bodyStorageID, nu_L ), "Lubrication Force" );

   // perform lubrication evaluation
   timeloop.addFuncAfterTimeStep( EvaluateLubricationForce( blocks, bodyStorageID, globalBodyStorage, id1, id2, velocity,
                                                            nu_L, radius, &timeloop, fileIO, sphSphTest, sphWallTest ), "Evaluate Lubrication Force" );

   // reset the forces and apply a constant velocity
   timeloop.addFuncAfterTimeStep( ForceReset( blocks, bodyStorageID, id1, id2 ), "Reset force on bodies");

   // add pe timesteps
   timeloop.addFuncAfterTimeStep( pe_coupling::TimeStep( blocks, bodyStorageID, cr, syncCall, real_c(1.0), uint_c(1) ), "pe Time Step" );

   // remaining time steps
   timeloop.addFuncAfterTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;
   timeloop.run( timeloopTiming );
   timeloopTiming.logResultOnRoot();

   return 0;
}
} //namespace lubrication_correction_mem

int main( int argc, char **argv ){
   lubrication_correction_mem::main(argc, argv);
}


