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
//! \file SuViscoelasticityTest.cpp
//! \ingroup lbm
//! \author Cameron Stewart <cstewart@icp.uni-stuttgart.de>
//! \brief D2Q9 Poiseuille flow simulation with Oldroyd-B viscoelasticity - Checks against analytical formula and
//!        reference data
//
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/BoundaryHandling.h"

#include "core/Environment.h"
#include "core/SharedFunctor.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/communication/PackInfo.h"

#include "lbm/boundary/all.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/lattice_model/D2Q9.h"
#include "lbm/lattice_model/ForceModel.h"
#include "lbm/SuViscoelasticity.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "timeloop/SweepTimeloop.h"

namespace walberla {

//////////////
// TYPEDEFS //
//////////////

typedef GhostLayerField< Vector3<real_t>, 1> ForceField_T;

typedef lbm::collision_model::TRT CollisionModel_T;
typedef lbm::force_model::GuoField< ForceField_T > ForceModel_T;
typedef lbm::D2Q9< CollisionModel_T, false, ForceModel_T >  LatticeModel_T;
typedef LatticeModel_T::Stencil    Stencil_T;

typedef GhostLayerField<Matrix3<real_t>, 1> StressField_T;
typedef GhostLayerField< Vector3<real_t>, 1> VelocityField_T;
typedef lbm::PdfField< LatticeModel_T >  PdfField_T;
typedef walberla::uint8_t flag_t;
typedef FlagField< flag_t > FlagField_T;
const uint_t FieldGhostLayers = 2;

typedef lbm::NoSlip< LatticeModel_T, flag_t> NoSlip_T;
typedef BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T> BoundaryHandling_T;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag( "fluid" );
const FlagUID NoSlip_Flag( "no slip" );

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID ) : flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;

}; // class MyBoundaryHandling


BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField       = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField        = block->getData< PdfField_T > ( pdfFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "moving obstacle boundary handling", flagField, fluid,
                                                           NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ) );

   const auto noSlip = flagField->getFlag(NoSlip_Flag);

   CellInterval domainBB = storage->getDomainCellBB();
   domainBB.xMin() -= cell_idx_c( 1 );
   domainBB.xMax() += cell_idx_c( 1 );

   domainBB.yMin() -= cell_idx_c( 1 );
   domainBB.yMax() += cell_idx_c( 1 );

   // SOUTH
   CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( south, *block );
   handling->forceBoundary( noSlip, south );

   // NORTH
   CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( north, *block );
   handling->forceBoundary( noSlip, north );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}

//////////////////////
// Poiseuille Force //
//////////////////////

template < typename BoundaryHandling_T >
class ConstantForce
{
public:

   ConstantForce( BlockDataID forceFieldId, BlockDataID boundaryHandlingId, real_t force)
         : forceFieldId_( forceFieldId ), boundaryHandlingId_(boundaryHandlingId), force_(force)
   {}

   void operator()( IBlock * block )
   {
      ForceField_T *forceField = block->getData< ForceField_T >(forceFieldId_);
      BoundaryHandling_T *boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingId_ );

      WALBERLA_FOR_ALL_CELLS_XYZ(forceField,
                                 {
                                    Cell cell(x,y,z);
                                    if (boundaryHandling->isDomain(cell)) {
                                       forceField->get(cell)[0] += force_;
                                    }
                                 })
   }

private:

   BlockDataID forceFieldId_, boundaryHandlingId_;
   real_t force_;
};

////////////////////
// Take Test Data //
////////////////////

class TestData
{
public:
   TestData(Timeloop & timeloop, shared_ptr< StructuredBlockForest > blocks, BlockDataID pdfFieldId, BlockDataID stressFieldId, uint_t timesteps, uint_t blockSize, real_t L, real_t H, real_t uExpected)
            : timeloop_(timeloop), blocks_(blocks), pdfFieldId_(pdfFieldId), stressFieldId_(stressFieldId), timesteps_(timesteps), blockSize_(blockSize), t0_(0), t1_(0), L_(L), H_(H), uMax_(0.0), uPrev_(0.0),
            uExpected_(uExpected), uSteady_(0.0) {}

   void operator()() {
      for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {
         PdfField_T *pdf = blockIt.get()->getData<PdfField_T>(pdfFieldId_);

         if (blockIt.get()->getAABB().contains(float_c(L_/2.0), float_c(H_/2.0), 0)) {
            uCurr_ = pdf->getVelocity(int32_c(blockSize_/2), int32_c(blockSize_/2), 0)[0]/uExpected_;
            tCurr_ = timeloop_.getCurrentTimeStep();
            if (tCurr_ == timesteps_ - 1){
               uSteady_ = uCurr_;
            }
            if (maxFlag_ == 0) {
               if (uCurr_ >= uPrev_) {
                  uMax_ = uCurr_;
               } else {
                  t0_ = tCurr_;
                  maxFlag_ = 1;
               }
               uPrev_ = uCurr_;
            }
            else if (maxFlag_ == 1) {
               if ((uCurr_ - 1.0) <= (uMax_ - 1.0) / std::exp(1)) {
                  t1_ = tCurr_ - t0_;
                  maxFlag_ = 2;
               }
            }
         }
      }
   }

   real_t getUSteady(){
      return uSteady_;
   }

   real_t getUMax(){
      return uMax_;
   }

   real_t getT0(){
      return real_c(t0_);
   }

   real_t getT1(){
      return real_c(t1_);
   }

private:
   Timeloop & timeloop_;
   shared_ptr< StructuredBlockForest > blocks_;
   BlockDataID pdfFieldId_, stressFieldId_;
   uint_t timesteps_, blockSize_, t0_, t1_, tCurr_;
   uint_t maxFlag_ = 0;
   real_t L_, H_, uMax_, uPrev_, uCurr_, uExpected_, uSteady_;
};

} // namespace walberla

//////////
// Main //
//////////

int main(int argc, char ** argv ){
   using namespace walberla;

   Environment env( argc, argv );

   // read parameter
   shared_ptr<StructuredBlockForest> blocks = blockforest::createUniformBlockGridFromConfig( env.config() );
   auto parameters = env.config()->getOneBlock( "Parameters" );

   // extract some constants from the parameters
   const real_t eta_s     = parameters.getParameter< real_t > ("eta_s");
   const real_t force     = parameters.getParameter< real_t > ("force");
   const real_t eta_p     = parameters.getParameter< real_t > ("eta_p");
   const real_t lambda_p  = parameters.getParameter< real_t > ("lambda_p");
   const uint_t period    = parameters.getParameter< uint_t > ("period");
   const real_t L         = parameters.getParameter< real_t > ("L");
   const real_t H         = parameters.getParameter< real_t > ("H");
   const uint_t blockSize = parameters.getParameter< uint_t > ("blockSize");
   const bool   shortrun  = parameters.getParameter<   bool > ("shortrun");
   const uint_t timesteps = shortrun ? uint_t(2) : parameters.getParameter< uint_t > ("timesteps");

   // reference data
   const real_t          uExpected       = force*H*H/(real_c(8.0)*(eta_s + eta_p));
   const real_t          uMax            = parameters.getParameter< real_t > ("uMax");
   const real_t          t0              = parameters.getParameter< real_t > ("t0");
   const real_t          t1              = parameters.getParameter< real_t > ("t1");

   // create fields
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", FieldGhostLayers);
   BlockDataID forceFieldId = field::addToStorage<ForceField_T>( blocks, "Force Field", Vector3<real_t>(0.0), field::zyxf, FieldGhostLayers);
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::TRT::constructWithMagicNumber( walberla::lbm::collision_model::omegaFromViscosity(eta_s)), lbm::force_model::GuoField<ForceField_T>( forceFieldId ) );
   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel, Vector3<real_t>(), real_c(1.0), FieldGhostLayers );
   BlockDataID stressId = walberla::field::addToStorage<StressField_T>( blocks, "Stress Field", Matrix3<real_t>(0.0), field::zyxf, FieldGhostLayers);
   BlockDataID stressOldId = walberla::field::addToStorage<StressField_T>( blocks, "Old Stress Field", Matrix3<real_t>(0.0), field::zyxf, FieldGhostLayers);
   BlockDataID velocityId = walberla::field::addToStorage<VelocityField_T> (blocks, "Velocity Field", Vector3<real_t>(0.0), field::zyxf, FieldGhostLayers);

   // add boundary handling
   BlockDataID boundaryHandlingId = blocks->addStructuredBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldId, pdfFieldId ), "boundary handling" );

   // create time loop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< Stencil_T > communication( blocks );
   communication.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldId ) );
   auto testData = make_shared< TestData >(TestData(timeloop, blocks, pdfFieldId, stressId, timesteps, blockSize, L, H, uExpected));

   // structure timeloop
   timeloop.add() << BeforeFunction( communication, "communication" )
                  << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingId ), "boundary handling" );
   timeloop.add() << BeforeFunction( lbm::viscoelastic::Su<LatticeModel_T, BoundaryHandling_T>(blocks, forceFieldId, pdfFieldId, boundaryHandlingId, stressId, stressOldId, velocityId,
                                                                                               lambda_p, eta_p, period, true), "viscoelasticity")
                  << Sweep( ConstantForce<BoundaryHandling_T>(forceFieldId, boundaryHandlingId, force),"Poiseuille Force");
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag ) ),
                            "LB stream & collide" )
                  << AfterFunction(makeSharedFunctor(testData), "test data");

   timeloop.run();


   if(!shortrun)
   {
      // compare to reference data
      real_t errSteady = real_c(fabs(testData->getUSteady() - real_c(1.0))/real_c(1.0));
      real_t errMax = real_c(fabs(testData->getUMax() - uMax)/uMax);
      real_t errt0 = real_c(fabs(testData->getT0() - t0)/t0);
      real_t errt1 = real_c(fabs(testData->getT1() - t1)/t1);

      WALBERLA_LOG_RESULT("Steady State Velocity Error: " << errSteady );
      WALBERLA_LOG_RESULT("Maximum Velocity Error: " << errMax );
      WALBERLA_LOG_RESULT("Time of Maximum Error: " << errt0 );
      WALBERLA_LOG_RESULT("Decay Time Error: " << errt1 );

      // check that errors < 1%
      if (errSteady < 0.01 && errMax < 0.01 && errt0 < 0.01 && errt1 < 0.01){
         WALBERLA_LOG_RESULT("Success" );
         return EXIT_SUCCESS;
      }
      else {
         WALBERLA_LOG_RESULT("Failure" );
         return EXIT_FAILURE;
      }
   } else{
      return EXIT_SUCCESS;
   }

}


