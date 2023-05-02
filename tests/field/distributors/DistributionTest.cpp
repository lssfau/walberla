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
//! \file DistributionTest.cpp
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
#include "field/distributors/all.h"

#include <vector>

namespace distribution_tests {

using namespace walberla;

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

void resetScalarField( const shared_ptr<StructuredBlockStorage> & blocks,
                       const BlockDataID & scalarFieldID )
{
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto sField = blockIt->getData<ScalarField_T>( scalarFieldID );
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(sField,
                                 sField->get(x,y,z) = real_t(0);
      );
   }
}

void resetVectorField( const shared_ptr<StructuredBlockStorage> & blocks,
                       const BlockDataID & vectorFieldID )
{
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto vField = blockIt->getData<Vec3Field_T>( vectorFieldID );
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(vField,
                                                       vField->get(x,y,z) = Vector3<real_t>(real_t(0));
      );
   }
}

void resetMultiCompField( const shared_ptr<StructuredBlockStorage> & blocks,
                          const BlockDataID & multiComponentFieldID )
{
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto mField = blockIt->getData<MultiComponentField_T>( multiComponentFieldID );
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(mField,
                                                       mField->get(x,y,z,0) = real_t(0);
                                                       mField->get(x,y,z,1) = real_t(0);
                                                       mField->get(x,y,z,2) = real_t(0);
      );
   }
}


void setBoundaryFlags( const shared_ptr<StructuredBlockStorage> & blocks,
                       const BlockDataID & flagFieldID )
{
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto flagField = blockIt->getData<FlagField_T>( flagFieldID );
      auto domainFlag = flagField->getOrRegisterFlag( Domain_Flag );
      auto boundaryFlag = flagField->getOrRegisterFlag( Boundary_Flag );
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(flagField,
                                                       if( x == 2)
                                                       {
                                                          flagField->removeFlag(x,y,z,domainFlag);
                                                          flagField->addFlag(x,y,z,boundaryFlag);
                                                       }
      );
   }
}

void getScalarFieldQuantities( const shared_ptr<StructuredBlockStorage> & blocks, const BlockDataID & flagFieldID,
                               const BlockDataID & fieldID,
                               real_t & summarizedValue, real_t & minValue, real_t & maxValue)
{
   real_t sum( real_t(0) );
   real_t min( real_t(0) );
   real_t max( real_t(0) );
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto field = blockIt->getData<ScalarField_T>( fieldID );
      auto flagField = blockIt->getData<FlagField_T>( flagFieldID );
      auto domainFlag = flagField->getFlag( Domain_Flag );

      CellInterval xyzSizeWithGhostLayers = field->xyzSizeWithGhostLayer();
      for( auto cellIt = xyzSizeWithGhostLayers.begin(); cellIt != xyzSizeWithGhostLayers.end(); ++cellIt )
      {
         if( flagField->isFlagSet(*cellIt,domainFlag))
         {
            real_t value = field->get(*cellIt);
            sum += value;
            min = std::min(min, value);
            max = std::max(max, value);
         }
      }
   }
   WALBERLA_MPI_SECTION()
   {
      mpi::allReduceInplace( sum, mpi::SUM );
      mpi::allReduceInplace( min, mpi::MIN );
      mpi::allReduceInplace( max, mpi::MAX );
   }

   summarizedValue = sum;
   minValue = min;
   maxValue = max;
}

void getVectorFieldQuantities( const shared_ptr<StructuredBlockStorage> & blocks, const BlockDataID & flagFieldID,
                               const BlockDataID & fieldID,
                               Vector3<real_t> & summarizedValue, Vector3<real_t> & minValue, Vector3<real_t> & maxValue)
{
   Vector3<real_t> sum( real_t(0) );
   Vector3<real_t> min( real_t(0) );
   Vector3<real_t> max( real_t(0) );
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto field = blockIt->getData<Vec3Field_T>( fieldID );
      auto flagField = blockIt->getData<FlagField_T>( flagFieldID );
      auto domainFlag = flagField->getFlag( Domain_Flag );
      CellInterval xyzSizeWithGhostLayers = field->xyzSizeWithGhostLayer();
      for( auto cellIt = xyzSizeWithGhostLayers.begin(); cellIt != xyzSizeWithGhostLayers.end(); ++cellIt )
      {
         if( flagField->isFlagSet(*cellIt,domainFlag))
         {
            Vector3<real_t> value = field->get(*cellIt);
            sum += value;
            for (size_t i = 0; i < 3; ++i) {
               min[i] = std::min(min[i], value[i]);
               max[i] = std::max(max[i], value[i]);
            }
         }
      }
   }
   WALBERLA_MPI_SECTION()
   {
      for( size_t i = 0; i < 3; ++i) {
         sum[i] = mpi::allReduce( sum[i], mpi::SUM );
         min[i] = mpi::allReduce( min[i], mpi::MIN );
         max[i] = mpi::allReduce( max[i], mpi::MAX );
      }
   }

   summarizedValue = sum;
   minValue = min;
   maxValue = max;
}

void getMultiCompFieldQuantities( const shared_ptr<StructuredBlockStorage> & blocks, const BlockDataID & flagFieldID,
                                  const BlockDataID & fieldID,
                                  std::vector<real_t> & summarizedValue, std::vector<real_t> & minValue, std::vector<real_t> & maxValue )
{
   std::vector<real_t> sum( 3 );
   std::vector<real_t> min( 3 );
   std::vector<real_t> max( 3 );
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto field = blockIt->getData<MultiComponentField_T>( fieldID );
      auto flagField = blockIt->getData<FlagField_T>( flagFieldID );
      auto domainFlag = flagField->getFlag( Domain_Flag );
      CellInterval xyzSizeWithGhostLayers = field->xyzSizeWithGhostLayer();
      for( auto cellIt = xyzSizeWithGhostLayers.begin(); cellIt != xyzSizeWithGhostLayers.end(); ++cellIt )
      {
         if( flagField->isFlagSet(*cellIt,domainFlag))
         {
            real_t value0 = field->get(*cellIt, 0);
            sum[0] += value0;
            min[0] = std::min(min[0], value0);
            max[0] = std::max(max[0], value0);
            real_t value1 = field->get(*cellIt, 1);
            sum[1] += value1;
            min[1] = std::min(min[1], value1);
            max[1] = std::max(max[1], value1);
            real_t value2 = field->get(*cellIt, 2);
            sum[2] += value2;
            min[2] = std::min(min[2], value2);
            max[2] = std::max(max[2], value2);
         }
      }
   }
   WALBERLA_MPI_SECTION()
   {
      mpi::allReduceInplace( sum, mpi::SUM );
      mpi::allReduceInplace( min, mpi::MIN );
      mpi::allReduceInplace( max, mpi::MAX );
   }

   summarizedValue = sum;
   minValue = min;
   maxValue = max;
}

void testNearestNeighborDistributor( const shared_ptr<StructuredBlockStorage> & blocks, const BlockDataID & flagFieldID,
                                     const BlockDataID & scalarFieldID, const BlockDataID & vectorFieldID, const BlockDataID & multiComponentFieldID )
{
   // distributors
   typedef field::NearestNeighborDistributor<ScalarField_T, FlagField_T>         ScalarDistributor_T;
   typedef field::NearestNeighborDistributor<Vec3Field_T, FlagField_T>           Vec3Distributor_T;
   typedef field::NearestNeighborDistributor<MultiComponentField_T, FlagField_T> MultiComponentDistributor_T;
   BlockDataID scalarDistributorID         = field::addDistributor< ScalarDistributor_T, FlagField_T >( blocks, scalarFieldID, flagFieldID, Domain_Flag );
   BlockDataID vectorDistributorID         = field::addDistributor< Vec3Distributor_T, FlagField_T >( blocks, vectorFieldID, flagFieldID, Domain_Flag );
   BlockDataID multiComponentDistributorID = field::addDistributor< MultiComponentDistributor_T, FlagField_T >( blocks, multiComponentFieldID, flagFieldID, Domain_Flag );

   // check scalar distribution
   {
      Vector3<real_t> distributionPoint(real_t(1.9), real_t(2.1), real_t(2.1));
      real_t distributionValue(real_t(100));
      auto containingBlockID = blocks->getBlock(distributionPoint);
      if( containingBlockID != nullptr )
      {
         auto distPtr = containingBlockID->getData<ScalarDistributor_T>(scalarDistributorID);
         distPtr->distribute(distributionPoint, &distributionValue);
      }
      real_t sum(real_t(0)), min(real_t(0)), max(real_t(0));
      getScalarFieldQuantities(blocks, flagFieldID, scalarFieldID, sum, min, max);
      WALBERLA_CHECK_FLOAT_EQUAL(sum, distributionValue, "NearestNeighborDistributor: Sum of scalar distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(min, real_t(0), "NearestNeighborDistributor: Min of scalar distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(max, distributionValue, "NearestNeighborDistributor: Max of scalar distribution failed!" );

      resetScalarField(blocks, scalarFieldID);
   }

   // check vector distribution
   {
      Vector3<real_t> distributionPoint(real_t(5.4),real_t(2.1),real_t(3.2));
      Vector3<real_t> distributionValue(real_t(100), real_t(-10), real_t(1));
      auto containingBlockID = blocks->getBlock(distributionPoint);
      if( containingBlockID != nullptr )
      {
         auto distPtr = containingBlockID->getData<Vec3Distributor_T>(vectorDistributorID);
         distPtr->distribute(distributionPoint, &distributionValue);
      }
      Vector3<real_t> sum(real_t(0)), min(real_t(0)), max(real_t(0));
      getVectorFieldQuantities(blocks, flagFieldID, vectorFieldID, sum, min, max);
      WALBERLA_CHECK_FLOAT_EQUAL(sum[0], distributionValue[0], "NearestNeighborDistributor: Sum of vec[0] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(sum[1], distributionValue[1], "NearestNeighborDistributor: Sum of vec[1] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(sum[2], distributionValue[2], "NearestNeighborDistributor: Sum of vec[2] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(min[0], real_t(0), "NearestNeighborDistributor: Min of vec[0] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(min[1], distributionValue[1], "NearestNeighborDistributor: Min of vec[1] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(min[2], real_t(0), "NearestNeighborDistributor: Min of vec[2] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(max[0], distributionValue[0], "NearestNeighborDistributor: Max of vec[0] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(max[1], real_t(0), "NearestNeighborDistributor: Max of vec[1] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(max[2], distributionValue[2], "NearestNeighborDistributor: Max of vec[2] distribution failed!" );

      resetVectorField(blocks, vectorFieldID);
   }

   // check multi component distribution
   {
      Vector3<real_t> distributionPoint(real_t(4.4),real_t(2.1),real_t(3.2));
      std::vector<real_t> distributionValue(3);
      distributionValue[0] = real_t(100);
      distributionValue[1] = real_t(-10);
      distributionValue[2] = real_t(1);
      auto containingBlockID = blocks->getBlock(distributionPoint);
      if( containingBlockID != nullptr )
      {
         auto distPtr = containingBlockID->getData<MultiComponentDistributor_T>(multiComponentDistributorID);
         distPtr->distribute(distributionPoint, distributionValue.begin());
      }
      std::vector<real_t> sum(3), min(3), max(3);
      getMultiCompFieldQuantities(blocks, flagFieldID, multiComponentFieldID, sum, min, max);
      WALBERLA_CHECK_FLOAT_EQUAL(sum[0], distributionValue[0], "NearestNeighborDistributor: Sum of Multi Component[0] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(sum[1], distributionValue[1], "NearestNeighborDistributor: Sum of Multi Component[1] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(sum[2], distributionValue[2], "NearestNeighborDistributor: Sum of Multi Component[2] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(min[0], real_t(0), "NearestNeighborDistributor: Min of Multi Component[0] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(min[1], distributionValue[1], "NearestNeighborDistributor: Min of Multi Component[1] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(min[2], real_t(0), "NearestNeighborDistributor: Min of Multi Component[2] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(max[0], distributionValue[0], "NearestNeighborDistributor: Max of Multi Component[0] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(max[1], real_t(0), "NearestNeighborDistributor: Max of Multi Component[1] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(max[2], distributionValue[2], "NearestNeighborDistributor: Max of Multi Component[2] distribution failed!" );

      resetMultiCompField(blocks, multiComponentFieldID);
   }
}


void testKernelDistributor( const shared_ptr<StructuredBlockStorage> & blocks, const BlockDataID & flagFieldID,
                            const BlockDataID & scalarFieldID, const BlockDataID & vectorFieldID, const BlockDataID & multiComponentFieldID )
{
   // distributors
   typedef field::KernelDistributor<ScalarField_T, FlagField_T>         ScalarDistributor_T;
   typedef field::KernelDistributor<Vec3Field_T, FlagField_T>           Vec3Distributor_T;
   typedef field::KernelDistributor<MultiComponentField_T, FlagField_T> MultiComponentDistributor_T;
   BlockDataID scalarDistributorID         = field::addDistributor< ScalarDistributor_T, FlagField_T >( blocks, scalarFieldID, flagFieldID, Domain_Flag );
   BlockDataID vectorDistributorID         = field::addDistributor< Vec3Distributor_T, FlagField_T >( blocks, vectorFieldID, flagFieldID, Domain_Flag );
   BlockDataID multiComponentDistributorID = field::addDistributor< MultiComponentDistributor_T, FlagField_T >( blocks, multiComponentFieldID, flagFieldID, Domain_Flag );

   // check scalar distribution
   {
      Vector3<real_t> distributionPoint(real_t(1.9), real_t(2.1), real_t(2.1));
      real_t distributionValue(real_t(100));
      auto containingBlockID = blocks->getBlock(distributionPoint);
      if( containingBlockID != nullptr )
      {
         auto distPtr = containingBlockID->getData<ScalarDistributor_T>(scalarDistributorID);
         distPtr->distribute(distributionPoint, &distributionValue);
      }
      real_t sum(real_t(0)), min(real_t(0)), max(real_t(0));
      getScalarFieldQuantities(blocks, flagFieldID, scalarFieldID, sum, min, max);
      WALBERLA_CHECK_FLOAT_EQUAL(sum, distributionValue, "KernelDistributor: Sum of scalar distribution failed!" );
      WALBERLA_CHECK(min >= real_t(0), "KernelDistributor: Min of scalar distribution failed!" );
      WALBERLA_CHECK(max <= distributionValue, "KernelDistributor: Max of scalar distribution failed!" );

      resetScalarField(blocks, scalarFieldID);
   }

   // check vector distribution
   {
      Vector3<real_t> distributionPoint(real_t(5.4),real_t(2.1),real_t(3.2));
      Vector3<real_t> distributionValue(real_t(100), real_t(-10), real_t(1));
      auto containingBlockID = blocks->getBlock(distributionPoint);
      if( containingBlockID != nullptr )
      {
         auto distPtr = containingBlockID->getData<Vec3Distributor_T>(vectorDistributorID);
         distPtr->distribute(distributionPoint, &distributionValue);
      }
      Vector3<real_t> sum(real_t(0)), min(real_t(0)), max(real_t(0));
      getVectorFieldQuantities(blocks, flagFieldID, vectorFieldID, sum, min, max);
      WALBERLA_CHECK_FLOAT_EQUAL(sum[0], distributionValue[0], "KernelDistributor: Sum of vec[0] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(sum[1], distributionValue[1], "KernelDistributor: Sum of vec[1] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(sum[2], distributionValue[2], "KernelDistributor: Sum of vec[2] distribution failed!" );
      WALBERLA_CHECK(min[0] >= real_t(0), "KernelDistributor: Min of vec[0] distribution failed!" );
      WALBERLA_CHECK(min[1] >= distributionValue[1], "KernelDistributor: Min of vec[1] distribution failed!" );
      WALBERLA_CHECK(min[2] >= real_t(0), "KernelDistributor: Min of vec[2] distribution failed!" );
      WALBERLA_CHECK(max[0] <= distributionValue[0], "KernelDistributor: Max of vec[0] distribution failed!" );
      WALBERLA_CHECK(max[1] <= real_t(0), "KernelDistributor: Max of vec[1] distribution failed!" );
      WALBERLA_CHECK(max[2] <= distributionValue[2], "KernelDistributor: Max of vec[2] distribution failed!" );

      resetVectorField(blocks, vectorFieldID);
   }

   // check multi component distribution
   {
      Vector3<real_t> distributionPoint(real_t(4.4),real_t(2.1),real_t(3.2));
      std::vector<real_t> distributionValue(3);
      distributionValue[0] = real_t(100);
      distributionValue[1] = real_t(-10);
      distributionValue[2] = real_t(1);
      auto containingBlockID = blocks->getBlock(distributionPoint);
      if( containingBlockID != nullptr )
      {
         auto distPtr = containingBlockID->getData<MultiComponentDistributor_T>(multiComponentDistributorID);
         distPtr->distribute(distributionPoint, distributionValue.begin());
      }
      std::vector<real_t> sum(3), min(3), max(3);
      getMultiCompFieldQuantities(blocks, flagFieldID, multiComponentFieldID, sum, min, max);
      WALBERLA_CHECK_FLOAT_EQUAL(sum[0], distributionValue[0], "KernelDistributor: Sum of Multi Component[0] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(sum[1], distributionValue[1], "KernelDistributor: Sum of Multi Component[1] distribution failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(sum[2], distributionValue[2], "KernelDistributor: Sum of Multi Component[2] distribution failed!" );
      WALBERLA_CHECK(min[0] >= real_t(0), "KernelDistributor: Min of Multi Component[0] distribution failed!" );
      WALBERLA_CHECK(min[1] >= distributionValue[1], "KernelDistributor: Min of Multi Component[1] distribution failed!" );
      WALBERLA_CHECK(min[2] >= real_t(0), "KernelDistributor: Min of Multi Component[2] distribution failed!" );
      WALBERLA_CHECK(max[0] <= distributionValue[0], "KernelDistributor: Max of Multi Component[0] distribution failed!" );
      WALBERLA_CHECK(max[1] <= real_t(0), "KernelDistributor: Max of Multi Component[1] distribution failed!" );
      WALBERLA_CHECK(max[2] <= distributionValue[2], "KernelDistributor: Max of Multi Component[2] distribution failed!" );

      resetMultiCompField(blocks, multiComponentFieldID);
   }
}


void testNearestNeighborDistributorAtBoundary( const shared_ptr<StructuredBlockStorage> & blocks,
                                               const BlockDataID & flagFieldID, const BlockDataID & scalarFieldID )
{
   // distributor
   typedef field::NearestNeighborDistributor<ScalarField_T, FlagField_T> ScalarDistributor_T;
   BlockDataID scalarDistributorID = field::addDistributor<ScalarDistributor_T, FlagField_T>(blocks, scalarFieldID, flagFieldID, Domain_Flag);

   // check scalar interpolation close to boundary
   {
      Vector3<real_t> distributionPoint(real_t(1.9), real_t(2.1), real_t(2.1));
      real_t distributionValue(real_t(100));
      auto containingBlockID = blocks->getBlock(distributionPoint);
      if (containingBlockID != nullptr) {
         auto distPtr = containingBlockID->getData<ScalarDistributor_T>(scalarDistributorID);
         distPtr->distribute(distributionPoint, &distributionValue);
      }
      real_t sum(real_t(0)), min(real_t(0)), max(real_t(0));
      getScalarFieldQuantities(blocks, flagFieldID, scalarFieldID, sum, min, max);
      WALBERLA_CHECK_FLOAT_EQUAL(sum, distributionValue, "NearestNeighborDistributor: Sum of scalar distribution near boundary failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(min, real_t(0), "NearestNeighborDistributor: Min of scalar distribution near boundary failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(max, distributionValue, "NearestNeighborDistributor: Max of scalar distribution near boundary failed!" );

      resetScalarField(blocks, scalarFieldID);
   }

   // check scalar interpolation inside boundary
   {
      Vector3<real_t> distributionPoint(real_t(2.7), real_t(2.1), real_t(1.1));
      real_t distributionValue(real_t(100));
      auto containingBlockID = blocks->getBlock(distributionPoint);
      if (containingBlockID != nullptr) {
         auto distPtr = containingBlockID->getData<ScalarDistributor_T>(scalarDistributorID);
         distPtr->distribute(distributionPoint, &distributionValue);
      }
      real_t sum(real_t(0)), min(real_t(0)), max(real_t(0));
      getScalarFieldQuantities(blocks, flagFieldID, scalarFieldID, sum, min, max);
      WALBERLA_CHECK_FLOAT_EQUAL(sum, distributionValue, "NearestNeighborDistributor: Sum of scalar distribution inside boundary failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(min, real_t(0), "NearestNeighborDistributor: Min of scalar distribution inside boundary failed!" );
      WALBERLA_CHECK_FLOAT_EQUAL(max, distributionValue, "NearestNeighborDistributor: Max of scalar distribution inside boundary failed!" );

      resetScalarField(blocks, scalarFieldID);
   }
}

void testKernelDistributorAtBoundary( const shared_ptr<StructuredBlockStorage> & blocks,
                                      const BlockDataID & flagFieldID, const BlockDataID & scalarFieldID )
{
   // distributor
   typedef field::KernelDistributor<ScalarField_T, FlagField_T> ScalarDistributor_T;
   BlockDataID scalarDistributorID = field::addDistributor<ScalarDistributor_T, FlagField_T>(blocks, scalarFieldID, flagFieldID, Domain_Flag);

   // check scalar interpolation close to boundary
   {
      Vector3<real_t> distributionPoint(real_t(1.9), real_t(2.1), real_t(2.1));
      real_t distributionValue(real_t(100));
      auto containingBlockID = blocks->getBlock(distributionPoint);
      if (containingBlockID != nullptr) {
         auto distPtr = containingBlockID->getData<ScalarDistributor_T>(scalarDistributorID);
         distPtr->distribute(distributionPoint, &distributionValue);
      }
      real_t sum(real_t(0)), min(real_t(0)), max(real_t(0));
      getScalarFieldQuantities(blocks, flagFieldID, scalarFieldID, sum, min, max);
      WALBERLA_CHECK_FLOAT_EQUAL(sum, distributionValue, "KernelDistributor: Sum of scalar distribution near boundary failed!" );
      WALBERLA_CHECK(min >= real_t(0), "KernelDistributor: Min of scalar distribution near boundary failed!" );
      WALBERLA_CHECK(max <= distributionValue, "KernelDistributor: Max of scalar distribution near boundary failed!" );

      resetScalarField(blocks, scalarFieldID);
   }

   // check scalar interpolation inside boundary
   {
      Vector3<real_t> distributionPoint(real_t(2.7), real_t(2.1), real_t(1.1));
      real_t distributionValue(real_t(100));
      auto containingBlockID = blocks->getBlock(distributionPoint);
      if (containingBlockID != nullptr) {
         auto distPtr = containingBlockID->getData<ScalarDistributor_T>(scalarDistributorID);
         distPtr->distribute(distributionPoint, &distributionValue);
      }
      real_t sum(real_t(0)), min(real_t(0)), max(real_t(0));
      getScalarFieldQuantities(blocks, flagFieldID, scalarFieldID, sum, min, max);
      WALBERLA_CHECK_FLOAT_EQUAL(sum, distributionValue, "KernelDistributor: Sum of scalar distribution inside boundary failed!" );
      WALBERLA_CHECK(min >= real_t(0), "KernelDistributor: Min of scalar distribution inside boundary failed!" );
      WALBERLA_CHECK(max <= distributionValue, "KernelDistributor: Max of scalar distribution inside boundary failed!" );

      resetScalarField(blocks, scalarFieldID);
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
   BlockDataID scalarFieldID         = field::addToStorage< ScalarField_T >( blocks, "scalar field", real_t(0), field::fzyx, FieldGhostLayers );
   BlockDataID vectorFieldID         = field::addToStorage< Vec3Field_T >( blocks, "vec3 field", Vector3<real_t>(real_t(0)), field::fzyx, FieldGhostLayers );
   BlockDataID multiComponentFieldID = field::addToStorage< MultiComponentField_T >( blocks, "multi component field", real_t(0), field::fzyx, FieldGhostLayers );

   // test all distributors with domain flags everywhere, i.e. without special boundary treatment necessary
   testNearestNeighborDistributor(blocks, flagFieldID, scalarFieldID, vectorFieldID, multiComponentFieldID);
   testKernelDistributor(blocks, flagFieldID, scalarFieldID, vectorFieldID, multiComponentFieldID);

   // set some boundary flags in flag field and invalidate the corresponding scalar field values
   setBoundaryFlags(blocks, flagFieldID );

   // test all distributors' behavior close to boundary cells
   testNearestNeighborDistributorAtBoundary(blocks, flagFieldID, scalarFieldID);
   testKernelDistributorAtBoundary(blocks, flagFieldID, scalarFieldID);

   return 0;
}

} // namespace field_distribution_tests

int main( int argc, char **argv ){
   distribution_tests::main(argc, argv);
}