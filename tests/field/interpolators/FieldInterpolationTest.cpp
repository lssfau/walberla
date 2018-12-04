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
//! \file FieldInterpolationTest.cpp
//! \ingroup field
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"
#include "core/math/all.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/interpolators/all.h"

#include <vector>

namespace walberla {

namespace field_interpolation_tests {

const uint_t FieldGhostLayers( 1 );

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

typedef GhostLayerField< real_t, 1>          ScalarField_T;
typedef GhostLayerField< Vector3<real_t>, 1> Vec3Field_T;
typedef GhostLayerField< real_t, 3>          MultiComponentField_T;


const FlagUID Domain_Flag ( "domain" );
const FlagUID Boundary_Flag ( "boundary" );

void initFlagField( FlagField_T * field, IBlock * const /*block*/ )
{
   auto domainFlag = field->getOrRegisterFlag( Domain_Flag );
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( field, field->addFlag( x, y, z, domainFlag ); );
}

void initScalarField( const shared_ptr<StructuredBlockStorage> & blocks, const BlockDataID & scalarFieldID )
{
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto field = blockIt->getData<ScalarField_T>( scalarFieldID );
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(field,
                                 const Vector3<real_t> p = blocks->getBlockLocalCellCenter(*blockIt, Cell(x,y,z)) - Vector3<real_t>(real_t(0.5));
                                 //field->get(x,y,z) = real_t(2);
                                 field->get(x,y,z) = p[0];
      );
   }
}

void initVectorField( const shared_ptr<StructuredBlockStorage> & blocks, const BlockDataID & vectorFieldID )
{
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto field = blockIt->getData<Vec3Field_T>( vectorFieldID );
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(field,
                                 const Vector3<real_t> p = blocks->getBlockLocalCellCenter(*blockIt, Cell(x,y,z)) - Vector3<real_t>(real_t(0.5));
                                 field->get(x,y,z) = Vector3<real_t>(p[0], real_t(2), p[1]);
      );
   }
}

void initMultiComponentField( const shared_ptr<StructuredBlockStorage> & blocks, const BlockDataID & multiComponentFieldID )
{
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto field = blockIt->getData<MultiComponentField_T>( multiComponentFieldID );
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(field,
                                 const Vector3<real_t> p = blocks->getBlockLocalCellCenter(*blockIt, Cell(x,y,z)) - Vector3<real_t>(real_t(0.5));
                                 field->get(x,y,z,0) = p[0];
                                 field->get(x,y,z,1) = real_t(2);
                                 field->get(x,y,z,2) = p[1];
      );
   }
}

void setBoundaryFlags( const shared_ptr<StructuredBlockStorage> & blocks,
                       const BlockDataID & flagFieldID, const BlockDataID & scalarFieldID )
{
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto flagField = blockIt->getData<FlagField_T>( flagFieldID );
      auto valueField = blockIt->getData<ScalarField_T>( scalarFieldID );
      auto domainFlag = flagField->getOrRegisterFlag( Domain_Flag );
      auto boundaryFlag = flagField->getOrRegisterFlag( Boundary_Flag );
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(flagField,
                                                       if( x == 2)
                                                       {
                                                          flagField->removeFlag(x,y,z,domainFlag);
                                                          flagField->addFlag(x,y,z,boundaryFlag);
                                                          valueField->get(x,y,z) = std::numeric_limits<real_t>::max();
                                                       }
      );
   }
}

void testNearestNeighborFieldInterpolator( const shared_ptr<StructuredBlockStorage> & blocks, const BlockDataID & flagFieldID,
                                           const BlockDataID & scalarFieldID, const BlockDataID & vectorFieldID, const BlockDataID & multiComponentFieldID )
{
   // field interpolators
   typedef field::NearestNeighborFieldInterpolator<ScalarField_T, FlagField_T>         ScalarFieldInterpolator_T;
   typedef field::NearestNeighborFieldInterpolator<Vec3Field_T, FlagField_T>           Vec3FieldInterpolator_T;
   typedef field::NearestNeighborFieldInterpolator<MultiComponentField_T, FlagField_T> MultiComponentFieldInterpolator_T;
   BlockDataID scalarFieldInterpolatorID         = field::addFieldInterpolator< ScalarFieldInterpolator_T, FlagField_T >( blocks, scalarFieldID, flagFieldID, Domain_Flag );
   BlockDataID vectorFieldInterpolatorID         = field::addFieldInterpolator< Vec3FieldInterpolator_T, FlagField_T >( blocks, vectorFieldID, flagFieldID, Domain_Flag );
   BlockDataID multiComponentFieldInterpolatorID = field::addFieldInterpolator< MultiComponentFieldInterpolator_T, FlagField_T >( blocks, multiComponentFieldID, flagFieldID, Domain_Flag );

   // check scalar interpolation
   {
      Vector3<real_t> interpolationPoint(real_t(1.9), real_t(2.1), real_t(2.1));
      auto containingBlockID = blocks->getBlock(interpolationPoint);
      if( containingBlockID != nullptr )
      {
         real_t interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<ScalarFieldInterpolator_T>(scalarFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue, real_t(1), "NearestNeighborFieldInterpolator: Scalar interpolation failed");
      }
   }

   // check vector interpolation
   {
      Vector3<real_t> interpolationPoint(real_t(5.4),real_t(2.1),real_t(3.2));
      auto containingBlockID = blocks->getBlock( interpolationPoint );
      if( containingBlockID != nullptr ) {
         Vector3<real_t> interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<Vec3FieldInterpolator_T>(vectorFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[0], real_t(5), "NearestNeighborFieldInterpolator: Vec3[0] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[1], real_t(2), "NearestNeighborFieldInterpolator: Vec3[1] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[2], real_t(2), "NearestNeighborFieldInterpolator: Vec3[2] interpolation failed");
      }
   }

   // check multi component interpolation
   {
      Vector3<real_t> interpolationPoint(real_t(4.4),real_t(2.1),real_t(3.2));
      auto containingBlockID = blocks->getBlock( interpolationPoint );
      if( containingBlockID != nullptr ) {
         std::vector<real_t> interpolatedValue(3, real_t(0));
         auto interPtr = containingBlockID->getData<MultiComponentFieldInterpolator_T>(multiComponentFieldInterpolatorID);
         interPtr->get(interpolationPoint, interpolatedValue.begin());
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[0], real_t(4), "NearestNeighborFieldInterpolator: Multi Component[0] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[1], real_t(2), "NearestNeighborFieldInterpolator: Multi Component[1] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[2], real_t(2), "NearestNeighborFieldInterpolator: Multi Component[2] interpolation failed");
      }
   }
}

void testTrilinearFieldInterpolator( const shared_ptr<StructuredBlockStorage> & blocks, const BlockDataID & flagFieldID,
                                     const BlockDataID & scalarFieldID, const BlockDataID & vectorFieldID, const BlockDataID & multiComponentFieldID )
{
   // field interpolators
   typedef field::TrilinearFieldInterpolator<ScalarField_T, FlagField_T>         ScalarFieldInterpolator_T;
   typedef field::TrilinearFieldInterpolator<Vec3Field_T, FlagField_T>           Vec3FieldInterpolator_T;
   typedef field::TrilinearFieldInterpolator<MultiComponentField_T, FlagField_T> MultiComponentFieldInterpolator_T;
   BlockDataID scalarFieldInterpolatorID         = field::addFieldInterpolator< ScalarFieldInterpolator_T, FlagField_T >( blocks, scalarFieldID, flagFieldID, Domain_Flag );
   BlockDataID vectorFieldInterpolatorID         = field::addFieldInterpolator< Vec3FieldInterpolator_T, FlagField_T >( blocks, vectorFieldID, flagFieldID, Domain_Flag );
   BlockDataID multiComponentFieldInterpolatorID = field::addFieldInterpolator< MultiComponentFieldInterpolator_T, FlagField_T >( blocks, multiComponentFieldID, flagFieldID, Domain_Flag );

   // check scalar interpolation
   {
      Vector3<real_t> interpolationPoint(real_t(1.9), real_t(2.1), real_t(2.1));
      auto containingBlockID = blocks->getBlock(interpolationPoint);
      if( containingBlockID != nullptr )
      {
         real_t interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<ScalarFieldInterpolator_T>(scalarFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue, real_t(1.4), "TrilinearFieldInterpolator: Scalar interpolation failed");
      }
   }

   // check vector interpolation
   {
      Vector3<real_t> interpolationPoint(real_t(5.4),real_t(2.1),real_t(3.2));
      auto containingBlockID = blocks->getBlock( interpolationPoint );
      if( containingBlockID != nullptr ) {
         Vector3<real_t> interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<Vec3FieldInterpolator_T>(vectorFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[0], real_t(4.9), "TrilinearFieldInterpolator: Vec3[0] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[1], real_t(2.0), "TrilinearFieldInterpolator: Vec3[1] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[2], real_t(1.6), "TrilinearFieldInterpolator: Vec3[2] interpolation failed");
      }
   }

   // check multi component interpolation
   {
      Vector3<real_t> interpolationPoint(real_t(4.4),real_t(2.1),real_t(3.2));
      auto containingBlockID = blocks->getBlock( interpolationPoint );
      if( containingBlockID != nullptr ) {
         std::vector<real_t> interpolatedValue(3, real_t(0));
         auto interPtr = containingBlockID->getData<MultiComponentFieldInterpolator_T>(multiComponentFieldInterpolatorID);
         interPtr->get(interpolationPoint, interpolatedValue.begin());
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[0], real_t(3.9), "TrilinearFieldInterpolator: Multi Component[0] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[1], real_t(2.0), "TrilinearFieldInterpolator: Multi Component[1] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[2], real_t(1.6), "TrilinearFieldInterpolator: Multi Component[2] interpolation failed");
      }
   }
}

void testKernelFieldInterpolator( const shared_ptr<StructuredBlockStorage> & blocks, const BlockDataID & flagFieldID,
                                  const BlockDataID & scalarFieldID, const BlockDataID & vectorFieldID, const BlockDataID & multiComponentFieldID )
{
   // field interpolators
   typedef field::KernelFieldInterpolator<ScalarField_T, FlagField_T>         ScalarFieldInterpolator_T;
   typedef field::KernelFieldInterpolator<Vec3Field_T, FlagField_T>           Vec3FieldInterpolator_T;
   typedef field::KernelFieldInterpolator<MultiComponentField_T, FlagField_T> MultiComponentFieldInterpolator_T;
   BlockDataID scalarFieldInterpolatorID         = field::addFieldInterpolator< ScalarFieldInterpolator_T, FlagField_T >( blocks, scalarFieldID, flagFieldID, Domain_Flag );
   BlockDataID vectorFieldInterpolatorID         = field::addFieldInterpolator< Vec3FieldInterpolator_T, FlagField_T >( blocks, vectorFieldID, flagFieldID, Domain_Flag );
   BlockDataID multiComponentFieldInterpolatorID = field::addFieldInterpolator< MultiComponentFieldInterpolator_T, FlagField_T >( blocks, multiComponentFieldID, flagFieldID, Domain_Flag );

   // check scalar interpolation
   {
      Vector3<real_t> interpolationPoint(real_t(1.9), real_t(2.1), real_t(2.1));
      auto containingBlockID = blocks->getBlock(interpolationPoint);
      if( containingBlockID != nullptr )
      {
         real_t interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<ScalarFieldInterpolator_T>(scalarFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue, real_t(1.4), "KernelFieldInterpolator: Scalar interpolation failed");
      }
   }

   // check vector interpolation
   {
      Vector3<real_t> interpolationPoint(real_t(5.4),real_t(2.1),real_t(3.2));
      auto containingBlockID = blocks->getBlock( interpolationPoint );
      if( containingBlockID != nullptr ) {
         Vector3<real_t> interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<Vec3FieldInterpolator_T>(vectorFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[0], real_t(4.9), "KernelFieldInterpolator: Vec3[0] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[1], real_t(2.0), "KernelFieldInterpolator: Vec3[1] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[2], real_t(1.6), "KernelFieldInterpolator: Vec3[2] interpolation failed");
      }
   }

   // check multi component interpolation
   {
      Vector3<real_t> interpolationPoint(real_t(4.4),real_t(2.1),real_t(3.2));
      auto containingBlockID = blocks->getBlock( interpolationPoint );
      if( containingBlockID != nullptr ) {
         std::vector<real_t> interpolatedValue(3, real_t(0));
         auto interPtr = containingBlockID->getData<MultiComponentFieldInterpolator_T>(multiComponentFieldInterpolatorID);
         interPtr->get(interpolationPoint, interpolatedValue.begin());
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[0], real_t(3.9), "KernelFieldInterpolator: Multi Component[0] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[1], real_t(2.0), "KernelFieldInterpolator: Multi Component[1] interpolation failed");
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue[2], real_t(1.6), "KernelFieldInterpolator: Multi Component[2] interpolation failed");
      }
   }
}


void testNearestNeighborFieldInterpolatorAtBoundary( const shared_ptr<StructuredBlockStorage> & blocks,
                                                     const BlockDataID & flagFieldID, const BlockDataID & scalarFieldID ) {
   // field interpolators
   typedef field::NearestNeighborFieldInterpolator<ScalarField_T, FlagField_T> ScalarFieldInterpolator_T;
   BlockDataID scalarFieldInterpolatorID = field::addFieldInterpolator<ScalarFieldInterpolator_T, FlagField_T>(blocks, scalarFieldID, flagFieldID, Domain_Flag);

   // check scalar interpolation close to boundary
   {
      Vector3<real_t> interpolationPoint(real_t(1.9), real_t(2.1), real_t(2.1));
      auto containingBlockID = blocks->getBlock(interpolationPoint);
      if (containingBlockID != nullptr) {
         real_t interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<ScalarFieldInterpolator_T>(scalarFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue, real_t(1),
                                    "NearestNeighborFieldInterpolator: Scalar interpolation near boundary failed");
      }
   }

   // check scalar interpolation inside boundary
   {
      Vector3<real_t> interpolationPoint(real_t(2.7), real_t(2.1), real_t(1.1));
      auto containingBlockID = blocks->getBlock(interpolationPoint);
      if (containingBlockID != nullptr) {
         real_t interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<ScalarFieldInterpolator_T>(scalarFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue, real_t(3),
                                    "NearestNeighborFieldInterpolator: Scalar interpolation inside boundary failed");
      }
   }
}

void testTrilinearFieldInterpolatorAtBoundary( const shared_ptr<StructuredBlockStorage> & blocks,
                                               const BlockDataID & flagFieldID, const BlockDataID & scalarFieldID ) {
   // field interpolators
   typedef field::TrilinearFieldInterpolator<ScalarField_T, FlagField_T> ScalarFieldInterpolator_T;
   BlockDataID scalarFieldInterpolatorID = field::addFieldInterpolator<ScalarFieldInterpolator_T, FlagField_T>(blocks, scalarFieldID, flagFieldID, Domain_Flag);

   // check scalar interpolation close to boundary
   {
      Vector3<real_t> interpolationPoint(real_t(1.9), real_t(2.1), real_t(2.1));
      auto containingBlockID = blocks->getBlock(interpolationPoint);
      if (containingBlockID != nullptr) {
         real_t interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<ScalarFieldInterpolator_T>(scalarFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue, real_t(1),
                                    "TrilinearFieldInterpolator: Scalar interpolation near boundary failed");
      }
   }

   // check scalar interpolation inside boundary
   {
      Vector3<real_t> interpolationPoint(real_t(2.7), real_t(2.1), real_t(1.1));
      auto containingBlockID = blocks->getBlock(interpolationPoint);
      if (containingBlockID != nullptr) {
         real_t interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<ScalarFieldInterpolator_T>(scalarFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         WALBERLA_CHECK_FLOAT_EQUAL(interpolatedValue, real_t(3),
                                    "TrilinearFieldInterpolator: Scalar interpolation inside boundary failed");
      }
   }
}

void testKernelFieldInterpolatorAtBoundary( const shared_ptr<StructuredBlockStorage> & blocks,
                                            const BlockDataID & flagFieldID, const BlockDataID & scalarFieldID ) {
   // field interpolators
   typedef field::KernelFieldInterpolator<ScalarField_T, FlagField_T> ScalarFieldInterpolator_T;
   BlockDataID scalarFieldInterpolatorID = field::addFieldInterpolator<ScalarFieldInterpolator_T, FlagField_T>(blocks, scalarFieldID, flagFieldID, Domain_Flag);

   // check scalar interpolation close to boundary
   {
      Vector3<real_t> interpolationPoint(real_t(1.9), real_t(2.1), real_t(2.1));
      auto containingBlockID = blocks->getBlock(interpolationPoint);
      if (containingBlockID != nullptr) {
         real_t interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<ScalarFieldInterpolator_T>(scalarFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         // kernel interpolation can not extrapolate values from the available ones (see comments in KernelFieldInterpolator.h)
         // it will thus yield a value between the available ones, which are 0 and 1
         WALBERLA_CHECK(interpolatedValue < real_t(1),
                        "KernelFieldInterpolator: Scalar interpolation near boundary failed");
         WALBERLA_CHECK(interpolatedValue > real_t(0),
                        "KernelFieldInterpolator: Scalar interpolation near boundary failed");
      }
   }

   // check scalar interpolation inside boundary
   {
      Vector3<real_t> interpolationPoint(real_t(2.7), real_t(2.1), real_t(1.1));
      auto containingBlockID = blocks->getBlock(interpolationPoint);
      if (containingBlockID != nullptr) {
         real_t interpolatedValue(real_t(0));
         auto interPtr = containingBlockID->getData<ScalarFieldInterpolator_T>(scalarFieldInterpolatorID);
         interPtr->get(interpolationPoint, &interpolatedValue);
         // values has to be between the available ones, i.e. 1 and 3
         WALBERLA_CHECK(interpolatedValue > real_t(1),
                        "KernelFieldInterpolator: Scalar interpolation inside boundary failed");
         WALBERLA_CHECK(interpolatedValue < real_t(3),
                        "KernelFieldInterpolator: Scalar interpolation inside boundary failed");
      }
   }
}

int main(int argc, char **argv) {

   mpi::Environment mpiEnv(argc, argv);
   debug::enterTestMode();

   const uint_t numberOfBlocksInDirection = 2;
   const uint_t numberOfCellsPerBlockInDirection = 4;
   const real_t dx = real_t(1);

   // block storage
   auto blocks = blockforest::createUniformBlockGrid( numberOfBlocksInDirection, numberOfBlocksInDirection, numberOfBlocksInDirection,
                                                      numberOfCellsPerBlockInDirection, numberOfCellsPerBlockInDirection, numberOfCellsPerBlockInDirection,
                                                      dx, 0, false, false,
                                                      false, false, false,
                                                      false );
   // flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers, false, initFlagField );

   // data fields
   BlockDataID scalarFieldID         = field::addToStorage< ScalarField_T >( blocks, "scalar field", real_t(0), field::zyxf, FieldGhostLayers );
   BlockDataID vectorFieldID         = field::addToStorage< Vec3Field_T >( blocks, "vec3 field", Vector3<real_t>(real_t(0)), field::zyxf, FieldGhostLayers );
   BlockDataID multiComponentFieldID = field::addToStorage< MultiComponentField_T >( blocks, "multi component field", real_t(0), field::zyxf, FieldGhostLayers );

   initScalarField(blocks, scalarFieldID);
   initVectorField(blocks, vectorFieldID );
   initMultiComponentField(blocks, multiComponentFieldID );

   // test all interpolators with domain flags everywhere, i.e. without special boundary treatment necessary
   testNearestNeighborFieldInterpolator(blocks, flagFieldID, scalarFieldID, vectorFieldID, multiComponentFieldID);
   testTrilinearFieldInterpolator(blocks, flagFieldID, scalarFieldID, vectorFieldID, multiComponentFieldID);
   testKernelFieldInterpolator(blocks, flagFieldID, scalarFieldID, vectorFieldID, multiComponentFieldID);

   // set some boundary flags in flag field and invalidate the corresponding scalar field values
   setBoundaryFlags(blocks, flagFieldID, scalarFieldID);

   // test all interpolators' behavior close to boundary cells
   testNearestNeighborFieldInterpolatorAtBoundary(blocks, flagFieldID, scalarFieldID);
   testTrilinearFieldInterpolatorAtBoundary(blocks, flagFieldID, scalarFieldID);
   testKernelFieldInterpolatorAtBoundary(blocks, flagFieldID, scalarFieldID);

   return 0;
}

} // namespace field_interpolation_tests
}

int main( int argc, char **argv ){
   walberla::field_interpolation_tests::main(argc, argv);
}