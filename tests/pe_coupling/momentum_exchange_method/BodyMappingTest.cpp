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
//! \file BodyMappingTest.cpp
//! \ingroup pe_coupling
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
#include "core/mpi/MPIManager.h"

#include "field/AddToStorage.h"

#include "lbm/boundary/all.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"

#include "pe/basic.h"
#include "pe/utility/DestroyBody.h"
#include "pe/utility/GetBody.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include <pe_coupling/utility/all.h>

#include "vtk/all.h"
#include "field/vtk/all.h"

namespace body_mapping_test
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

// PDF field, flag field & body field
using LatticeModel_T = lbm::D3Q19<lbm::collision_model::SRT>;

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

const uint_t FieldGhostLayers = 1;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
typedef GhostLayerField< pe::BodyID, 1 >  BodyField_T;

// boundary handling
typedef lbm::NoSlip< LatticeModel_T, flag_t > NoSlip_T;
typedef pe_coupling::SimpleBB< LatticeModel_T, FlagField_T >  MO_T;
typedef BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, MO_T > BoundaryHandling_T;

using BodyTypeTuple = std::tuple<pe::Sphere> ;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag ( "fluid" );
const FlagUID MO_Flag ( "moving obstacle" );
const FlagUID FormerMO_Flag ( "former moving obstacle" );
const FlagUID NoSlip_Flag  ( "no slip" );


/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID, const BlockDataID & bodyFieldID ) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), bodyFieldID_ ( bodyFieldID ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID bodyFieldID_;

}; // class MyBoundaryHandling

BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField       = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField        = block->getData< PdfField_T > ( pdfFieldID_ );
   BodyField_T * bodyField       = block->getData< BodyField_T >( bodyFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "fixed obstacle boundary handling", flagField, fluid,
                                                           NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
                                                           MO_T (  "MO_BB",  MO_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ) );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}


class MappingChecker
{
public:
   MappingChecker(const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & bodyStorageID, const shared_ptr<pe::BodyStorage> & globalBodyStorage,
                  const BlockDataID & boundaryHandlingID, const BlockDataID & bodyFieldID, real_t sphereRadius) :
         blocks_( blocks ), bodyStorageID_( bodyStorageID ), globalBodyStorage_( globalBodyStorage ),
         boundaryHandlingID_( boundaryHandlingID ), bodyFieldID_( bodyFieldID ),
         sphereVolume_( math::pi * real_t(4) / real_t(3) * sphereRadius * sphereRadius * sphereRadius )
   { }

   // check the mapping in the inner domain of the block and check mapped volume against real sphere volume
   void operator()(std::string & testIdentifier)
   {
      uint_t cellCounter( uint_t(0) );
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         BoundaryHandling_T * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         FlagField_T *        flagField        = boundaryHandling->getFlagField();

         auto xyzSize = flagField->xyzSize();

         // look for bodies in body storage
         for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt ) {
            for (auto cellIt = xyzSize.begin(); cellIt != xyzSize.end(); ++cellIt) {
               Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter(*blockIt, *cellIt);
               Vector3<real_t> distance = cellCenter - bodyIt->getPosition();
               if (bodyIt->containsPoint(cellCenter)) {
                  WALBERLA_CHECK(boundaryHandling->isBoundary(*cellIt),
                                 testIdentifier << "Invalid mapping in cell " << *cellIt
                                                << " with center at " << cellCenter
                                                << " - expected boundary cell. Distance cell center - body center = "
                                                << distance.length() << ".");
                  ++cellCounter;
               } else {
                  WALBERLA_CHECK(boundaryHandling->isDomain(*cellIt),
                                 testIdentifier << "Invalid mapping in cell " << *cellIt
                                                << " with center at " << cellCenter
                                                << " - expected domain cell. Distance cell center - body center = "
                                                << distance.length() << ".");
               }
            }
         }
         // look for bodies in global body storage
         for( auto bodyIt = globalBodyStorage_->begin(); bodyIt != globalBodyStorage_->end(); ++bodyIt )
         {
            for (auto cellIt = xyzSize.begin(); cellIt != xyzSize.end(); ++cellIt) {
               Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter(*blockIt, *cellIt);
               Vector3<real_t> distance = cellCenter - bodyIt->getPosition();
               if (bodyIt->containsPoint(cellCenter)) {
                  WALBERLA_CHECK(boundaryHandling->isBoundary(*cellIt),
                                 testIdentifier << "Invalid mapping in cell " << *cellIt
                                                << " with center at " << cellCenter
                                                << " - expected boundary cell. Distance cell center - body center = "
                                                << distance.length() << ".");
                  ++cellCounter;
               } else {
                  if( globalBodyStorage_->size() <= uint_t(1) ) {
                     WALBERLA_CHECK(boundaryHandling->isDomain(*cellIt),
                                    testIdentifier << "Invalid mapping in cell " << *cellIt
                                                   << " with center at " << cellCenter
                                                   << " - expected domain cell. Distance cell center - body center = "
                                                   << distance.length() << ".");
                  }
               }
            }
         }
      }
      // mapped volume has to be - approximately - the same as the real volume
      real_t mappedVolume = real_c(cellCounter); // dx=1
      WALBERLA_CHECK(std::fabs( mappedVolume - sphereVolume_ ) / sphereVolume_ <= real_t(0.1),
                     "Mapped volume " << mappedVolume << " does not fit to real sphere volume " << sphereVolume_ << ".");
   }

   // checks only the mapping in the ghost layers
   void checkGhostLayer(std::string & testIdentifier)
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         BoundaryHandling_T * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         FlagField_T *        flagField        = boundaryHandling->getFlagField();

         auto xyzSizeWGl = flagField->xyzSizeWithGhostLayer();

         // look for bodies in body storage
         for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            for( auto cellIt = xyzSizeWGl.begin(); cellIt != xyzSizeWGl.end(); ++cellIt )
            {
               if( flagField->isInInnerPart(*cellIt) ) continue;

               Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter( *blockIt, *cellIt );
               Vector3<real_t> distance = cellCenter - bodyIt->getPosition();
               if( bodyIt->containsPoint(cellCenter))
               {
                  WALBERLA_CHECK( boundaryHandling->isBoundary(*cellIt),
                                  testIdentifier << ": Invalid mapping in ghost layer cell " << *cellIt
                                                 << " with center at " << cellCenter
                                                 << " - expected boundary cell. Distance cell center - body center = "
                                                 << distance.length() << ".");
               }
               else
               {
                  WALBERLA_CHECK( boundaryHandling->isDomain(*cellIt),
                                  testIdentifier << ": Invalid mapping in ghost layer cell " << *cellIt
                                                 << " with center at " << cellCenter
                                                 << " - expected domain cell. Distance cell center - body center = "
                                                 << distance.length() << "." );
               }
            }
         }

         // look for bodies in global body storage
         for( auto bodyIt = globalBodyStorage_->begin(); bodyIt != globalBodyStorage_->end(); ++bodyIt )
         {
            for( auto cellIt = xyzSizeWGl.begin(); cellIt != xyzSizeWGl.end(); ++cellIt )
            {

               if( flagField->isInInnerPart(*cellIt) ) continue;

               Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter( *blockIt, *cellIt );
               Vector3<real_t> distance = cellCenter - bodyIt->getPosition();
               if( bodyIt->containsPoint(cellCenter))
               {
                  WALBERLA_CHECK( boundaryHandling->isBoundary(*cellIt),
                                  testIdentifier << ": Invalid mapping in ghost layer cell " << *cellIt
                                                 << " with center at " << cellCenter
                                                 << " - expected boundary cell. Distance cell center - body center = "
                                                 << distance.length() << ".");
               }
               else
               {
                  if( globalBodyStorage_->size() <= uint_t(1) ) {
                     WALBERLA_CHECK(boundaryHandling->isDomain(*cellIt),
                                    testIdentifier << ": Invalid mapping in ghost layer cell " << *cellIt
                                                   << " with center at " << cellCenter
                                                   << " - expected domain cell. Distance cell center - body center = "
                                                   << distance.length() << ".");
                  }
               }
            }
         }
      }
   }

   // check if the pointer in the body field is the correct one
   void checkBodyField(std::string & testIdentifier)
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {

         BodyField_T * bodyField = blockIt->getData< BodyField_T >( bodyFieldID_ );

         auto xyzSizeWGl = bodyField->xyzSizeWithGhostLayer();

         // look for bodies in body storage
         for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            for( auto cellIt = xyzSizeWGl.begin(); cellIt != xyzSizeWGl.end(); ++cellIt )
            {
               Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter( *blockIt, *cellIt );
               Vector3<real_t> distance = cellCenter - bodyIt->getPosition();
               if( bodyIt->containsPoint(cellCenter) )
               {
                  WALBERLA_CHECK_EQUAL( bodyIt->getSystemID(), (*bodyField)(*cellIt)->getSystemID(),
                                        testIdentifier << "Invalid mapping in cell " << *cellIt
                                                       << " with center at " << cellCenter
                                                       << " - expected body. Distance cell center - body center = "
                                                       << distance.length() << ".");
               }
               else
               {
                  WALBERLA_CHECK_NULLPTR( bodyField->get(*cellIt),
                                          testIdentifier << "Invalid mapping in cell " << *cellIt
                                                         << " with center at " << cellCenter
                                                         << " - expected NULL pointer. Distance cell center - body center = "
                                                         << distance.length() << "." );
               }
            }
         }

         // look for bodies in global body storage
         for( auto bodyIt = globalBodyStorage_->begin(); bodyIt != globalBodyStorage_->end(); ++bodyIt )
         {
            for( auto cellIt = xyzSizeWGl.begin(); cellIt != xyzSizeWGl.end(); ++cellIt )
            {
               Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter( *blockIt, *cellIt );
               Vector3<real_t> distance = cellCenter - bodyIt->getPosition();
               if( bodyIt->containsPoint(cellCenter))
               {
                  WALBERLA_CHECK_EQUAL( bodyIt->getSystemID(), (*bodyField)(*cellIt)->getSystemID(),
                                  testIdentifier << "Invalid mapping in cell " << *cellIt
                                                 << " with center at " << cellCenter
                                                 << " - expected body. Distance cell center - body center = "
                                                 << distance.length() << ".");
               }
               else
               {
                  if( globalBodyStorage_->size() <= uint_t(1) )
                  {
                     WALBERLA_CHECK_NULLPTR( bodyField->get(*cellIt),
                                     testIdentifier << "Invalid mapping in cell " << *cellIt
                                                    << " with center at " << cellCenter
                                                    << " - expected NULL pointer. Distance cell center - body center = "
                                                    << distance.length() << ".");
                  }
               }
            }
         }
      }
   }


private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID bodyStorageID_;
   shared_ptr<pe::BodyStorage> globalBodyStorage_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID bodyFieldID_;
   real_t sphereVolume_;

};

class MappingResetter
{
public:
   MappingResetter(const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & boundaryHandlingID, const BlockDataID & bodyFieldID) :
         blocks_( blocks ), boundaryHandlingID_( boundaryHandlingID ), bodyFieldID_ ( bodyFieldID )
   { }

   void operator()()
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         BoundaryHandling_T * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         FlagField_T *        flagField        = boundaryHandling->getFlagField();
         BodyField_T *        bodyField        = blockIt->getData< BodyField_T >( bodyFieldID_ );

         auto xyzSizeWGl = flagField->xyzSizeWithGhostLayer();

         // reset to domain (fluid)
         boundaryHandling->forceDomain(xyzSizeWGl);

         for( auto cellIt = xyzSizeWGl.begin(); cellIt != xyzSizeWGl.end(); ++cellIt )
         {
            // reset body field
            bodyField->get(*cellIt) = NULL;
         }
      }
   }

private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID bodyFieldID_;
};


/*!\brief Test case for the available body mapping functions
 *
 * There are two types of body mapping functions:
 *  - src/pe_coupling/mapping/BodyMapping.h:
 *      - to map a standard LBM boundary condition to the body, e.g. NoSlip
 *  - src/pe_coupling/momentum_exchange_method/BodyMapping.h:
 *      - to map the coupling boundaries (rc/pe_coupling/momentum_exchange_method/boundary) to the body
 *      - the surface velocity of the mapped body is seen by the fluid
 *      - forces acting on the mapped bodies are calculated and set onto the body
 *
 * The following mapping functions are tested:
 *  - mapAllBodies
 *  - selectRegularBodies
 *  - selectFixedBodies
 *  - selectGlobalBodies
 * The mapping inside a block, at a regular block boarder and at a periodic block boarder is tested.
 * The No-Slip mapping as well as the moving boundary mapping is tested.
 * Regular, global and infinite-mass spheres are tested.
 * After each test, the sphere is destroyed and the mapping is erased from the datastructures.
 *
 * Note that global bodies do not see the periodicity and thus the periodic copy has to be created explicitly.
 *
 */
//////////
// MAIN //
//////////
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   bool writeVTK = false;
   const real_t omega  = real_t(1);
   const real_t dx     = real_t(1);
   const real_t radius = real_t(5);

   ///////////////////////////
   // DATA STRUCTURES SETUP //
   ///////////////////////////

   Vector3<uint_t> blocksPerDirection(uint_t(3), uint_t(1), uint_t(1));
   Vector3<uint_t> cellsPerBlock(uint_t(20), uint_t(20), uint_t(20));
   Vector3<bool> periodicity(true, false, false);

   // create fully periodic domain with refined blocks
   auto blocks = blockforest::createUniformBlockGrid( blocksPerDirection[0], blocksPerDirection[1], blocksPerDirection[2],
                                                      cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],
                                                      dx,
                                                      0, false, false,
                                                      periodicity[0], periodicity[1], periodicity[2],
                                                      false );

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( omega );

   // add PDF field ( uInit = <0.1,0,0>, rhoInit = 1 )
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (zyxf)", latticeModel,
                                                                         Vector3<real_t>(real_t(0)), real_t(1),
                                                                         FieldGhostLayers, field::zyxf );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::zyxf, FieldGhostLayers );

   // add boundary handling
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID ), "boundary handling" );

   // pe body storage
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "Storage");
   auto sphereMaterialID = pe::createMaterial( "sphereMat", real_t(1) , real_t(0.3), real_t(0.2), real_t(0.2), real_t(0.24), real_t(200), real_t(200), real_t(0), real_t(0) );

   // pe coupling
   const real_t overlap = real_t( 1.5 ) * dx;

   // vtk
   auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field" );
   flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );

   // testing utilities
   MappingChecker mappingChecker(blocks, bodyStorageID, globalBodyStorage, boundaryHandlingID, bodyFieldID, radius);
   MappingResetter mappingResetter(blocks, boundaryHandlingID, bodyFieldID);
   pe_coupling::BodyMapping<LatticeModel_T, BoundaryHandling_T> regularBodyMapper(blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, MO_Flag, FormerMO_Flag, pe_coupling::selectRegularBodies );
   pe_coupling::BodyMapping<LatticeModel_T, BoundaryHandling_T> globalBodyMapper(blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, MO_Flag, FormerMO_Flag, pe_coupling::selectGlobalBodies );


   // sphere positions for test scenarios
   Vector3<real_t> positionInsideBlock(real_t(10), real_t(10), real_t(10));
   Vector3<real_t> positionAtBlockBoarder(real_t(19), real_t(10), real_t(10));
   Vector3<real_t> positionAtPeriodicBoarder(real_t(1), real_t(10), real_t(10));

   /////////////////////
   // NO SLIP MAPPING //
   /////////////////////

   //////////////////////////////////
   // TEST SPHERE INSIDE ONE BLOCK //
   //////////////////////////////////

   ////////////////////
   // REGULAR SPHERE //
   ////////////////////
   {
      std::string testIdentifier("Test: regular sphere inside block with no slip mapping ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionInsideBlock, radius, sphereMaterialID, false, true, false);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectRegularBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   /////////////////////
   // GLOBAL SPHERE 1 //
   /////////////////////
   {
      std::string testIdentifier("Test: global sphere inside block with no slip mapping, case 1 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionInsideBlock, radius, sphereMaterialID, true, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectGlobalBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   //////////////////
   // FIXED SPHERE //
   //////////////////
   {
      std::string testIdentifier("Test: fixed sphere inside block with no slip mapping ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionInsideBlock, radius, sphereMaterialID, false, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);

      pe_coupling::mapBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectFixedBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);
   }

   //////////////////////////////////
   // TEST SPHERE AT BLOCK BOARDER //
   //////////////////////////////////

   ////////////////////
   // REGULAR SPHERE //
   ////////////////////
   {
      std::string testIdentifier("Test: regular sphere at block boarder with no slip mapping ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtBlockBoarder, radius, sphereMaterialID, false, true, false);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectRegularBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   /////////////////////
   // GLOBAL SPHERE 1 //
   /////////////////////
   {
      std::string testIdentifier("Test: global sphere at block boarder with no slip mapping, case 1 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtBlockBoarder, radius, sphereMaterialID, true, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectGlobalBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   //////////////////
   // FIXED SPHERE //
   //////////////////
   {
      std::string testIdentifier("Test: fixed sphere at block boarder with no slip mapping ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtBlockBoarder, radius, sphereMaterialID, false, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);

      pe_coupling::mapBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectFixedBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);
   }

   /////////////////////////////////////
   // TEST SPHERE AT PERIODIC BOARDER //
   /////////////////////////////////////

   ////////////////////
   // REGULAR SPHERE //
   ////////////////////
   {
      std::string testIdentifier("Test: regular sphere at periodic boarder with no slip mapping ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtPeriodicBoarder, radius, sphereMaterialID, false, true, false);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectRegularBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   /////////////////////
   // GLOBAL SPHERE 1 //
   /////////////////////
   {
      std::string testIdentifier("Test: global sphere at periodic boarder with no slip mapping, case 1 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtPeriodicBoarder, radius, sphereMaterialID, true, false, true);

      //NOTE: global bodies are not communicated, thus they do not follow periodicity!!!
      //workaround: create the periodic copy explicitly
      Vector3<real_t> positionAtPeriodicBoarderCopy(real_t(1) + real_c(blocksPerDirection[0]) * real_c(cellsPerBlock[0]), real_t(10), real_t(10));
      pe::SphereID sp2 = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                          positionAtPeriodicBoarderCopy, radius, sphereMaterialID, true, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectGlobalBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp2->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   //////////////////
   // FIXED SPHERE //
   //////////////////
   {
      std::string testIdentifier("Test: fixed sphere at periodic boarder with no slip mapping ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtPeriodicBoarder, radius, sphereMaterialID, false, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);

      pe_coupling::mapBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectFixedBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);
   }

   //////////////////////////
   // MOVING BODY MAPPING //
   /////////////////////////

   //////////////////////////////////
   // TEST SPHERE INSIDE ONE BLOCK //
   //////////////////////////////////

   //////////////////////
   // REGULAR SPHERE 1 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere inside block with moving body mapping, case 1 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionInsideBlock, radius, sphereMaterialID, false, true, false);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapMovingBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectRegularBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   //////////////////////
   // REGULAR SPHERE 2 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere inside block with moving body mapping, case 2 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionInsideBlock, radius, sphereMaterialID, false, true, false);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) regularBodyMapper(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   /////////////////////
   // GLOBAL SPHERE 1 //
   /////////////////////
   {
      std::string testIdentifier("Test: global sphere inside block with moving body mapping, case 1 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionInsideBlock, radius, sphereMaterialID, true, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapMovingBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectGlobalBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   //////////////////////
   // GLOBAL SPHERE  2 //
   //////////////////////
   {
      std::string testIdentifier("Test: global sphere inside block with moving body mapping, case 2 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionInsideBlock, radius, sphereMaterialID, true, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) globalBodyMapper(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   //////////////////
   // FIXED SPHERE //
   //////////////////
   {
      std::string testIdentifier("Test: fixed sphere inside block with moving body mapping ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionInsideBlock, radius, sphereMaterialID, false, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);

      pe_coupling::mapMovingBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectFixedBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);
   }

   //////////////////////////////////
   // TEST SPHERE AT BLOCK BOARDER //
   //////////////////////////////////

   //////////////////////
   // REGULAR SPHERE 1 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere at block boarder with moving body mapping, case 1 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtBlockBoarder, radius, sphereMaterialID, false, true, false);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapMovingBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectRegularBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   //////////////////////
   // REGULAR SPHERE 2 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere at block boarder with moving body mapping, case 2 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtBlockBoarder, radius, sphereMaterialID, false, true, false);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) regularBodyMapper(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   /////////////////////
   // GLOBAL SPHERE 1 //
   /////////////////////
   {
      std::string testIdentifier("Test: global sphere at block boarder with moving body mapping, case 1 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtBlockBoarder, radius, sphereMaterialID, true, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapMovingBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectGlobalBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   //////////////////////
   // GLOBAL SPHERE  2 //
   //////////////////////
   {
      std::string testIdentifier("Test: global sphere at block boarder with moving body mapping, case 2 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtBlockBoarder, radius, sphereMaterialID, true, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) globalBodyMapper(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   //////////////////
   // FIXED SPHERE //
   //////////////////
   {
      std::string testIdentifier("Test: fixed sphere at block boarder with moving body mapping ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtBlockBoarder, radius, sphereMaterialID, false, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);

      pe_coupling::mapMovingBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectFixedBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);
   }

   /////////////////////////////////////
   // TEST SPHERE AT PERIODIC BOARDER //
   /////////////////////////////////////

   ///////////////////////
   // REGULAR SPHERE  1 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere at periodic boarder with moving body mapping, case 1 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtPeriodicBoarder, radius, sphereMaterialID, false, true, false);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapMovingBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectRegularBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   //////////////////////
   // REGULAR SPHERE 2 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere at periodic boarder with moving body mapping, case 2 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtPeriodicBoarder, radius, sphereMaterialID, false, true, false);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) regularBodyMapper(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   /////////////////////
   // GLOBAL SPHERE 1 //
   /////////////////////
   {
      std::string testIdentifier("Test: global sphere at periodic boarder with moving body mapping, case 1 ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtPeriodicBoarder, radius, sphereMaterialID, true, false, true);

      //NOTE: global bodies are not communicated, thus they do not follow periodicity!!!
      //workaround: create the periodic copy explicitly
      Vector3<real_t> positionAtPeriodicBoarderCopy(real_t(1) + real_c(blocksPerDirection[0]) * real_c(cellsPerBlock[0]), real_t(10), real_t(10));
      pe::SphereID sp2 = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                          positionAtPeriodicBoarderCopy, radius, sphereMaterialID, true, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);

      pe_coupling::mapMovingBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectGlobalBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp2->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   }

   //////////////////
   // FIXED SPHERE //
   //////////////////
   {
      std::string testIdentifier("Test: fixed sphere at periodic boarder with moving body mapping ");
      WALBERLA_LOG_DEVEL(testIdentifier << " - started");
      pe::SphereID sp = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                                         positionAtPeriodicBoarder, radius, sphereMaterialID, false, false, true);

      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);

      pe_coupling::mapMovingBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectFixedBodies );

      if( writeVTK ) flagFieldVTK->write();

      mappingChecker(testIdentifier);
      mappingChecker.checkGhostLayer(testIdentifier);
      mappingChecker.checkBodyField(testIdentifier);
      mappingResetter();

      WALBERLA_LOG_DEVEL(testIdentifier << " - ended");

      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, sp->getSystemID());
      pe::syncNextNeighbors<BodyTypeTuple>(blocks->getBlockForest(), bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, true);
   }


   return 0;

}

} //namespace body_mapping_test

int main( int argc, char **argv ){
   body_mapping_test::main(argc, argv);
}
