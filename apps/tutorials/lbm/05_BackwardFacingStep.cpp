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
//! \file 05_BackwardFacingStep.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Amin Nabikhani
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "lbm/all.h"
#include "timeloop/all.h"

namespace walberla {

//! [typedefs]
using LatticeModel_T = lbm::D2Q9<lbm::collision_model::SRT>;
using Stencil_T = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;
//! [typedefs]

using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

//**********************************************************************************************************************
/*!
*   Class for determining (and logging) the normalized reattachment length after the backward-facing step
*/
//**********************************************************************************************************************
class ReattachmentLengthFinder
{
public:

   ReattachmentLengthFinder( const shared_ptr< StructuredBlockStorage > & blocks,
                             const BlockDataID & pdfFieldID, const BlockDataID & flagFieldID, const FlagUID & fluid,
                             const std::string & filename, const uint_t checkFrequency, const Vector3< real_t > stepSize) :
      blocks_( blocks ),
      pdfFieldID_( pdfFieldID ), flagFieldID_( flagFieldID ), fluid_( fluid ),
      filename_( filename ), checkFrequency_( checkFrequency ), stepSize_( stepSize ),
      executionCounter_( uint_c(0) )
   {
      // open and close file on root process - purpose: If the file already exists, its content will be erased.
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream fileLocBottom( filename_.c_str());
         fileLocBottom << "Time,[Locations on the BOTTOM Wall that Reattachment Occures (Normalized with Step Height) ]" << std::endl;
         fileLocBottom.close();
      }
   }

   void operator()()
   {

      ++executionCounter_;
      if( ( executionCounter_ - uint_c(1) ) % checkFrequency_ != 0 )
         return;

      // variables for storing the process local reattachment location
      std::vector<real_t> reattachmentLocations;

      // iterate all blocks stored locally on this process
      for( auto block = blocks_->begin(); block != blocks_->end(); ++block )
      {
         // retrieve the pdf field and the flag field from the block
         PdfField_T  *  pdfField = block->getData< PdfField_T >( pdfFieldID_ );
         FlagField_T * flagField = block->getData< FlagField_T >( flagFieldID_ );

         // retrieve the lattice model from the pdf field
         const auto & latticeModel = pdfField->latticeModel();

         // get the flag that marks a cell as being fluid
         auto fluid = flagField->getFlag( fluid_ );

         // iterate all cells of the current block
         for( auto cell = pdfField->beginXYZ(); cell != pdfField->end(); ++cell )
         {
            const cell_idx_t x = cell.x();
            const cell_idx_t y = cell.y();
            const cell_idx_t z = cell.z();

            // If the current cell is marked as being fluid ...
            if( flagField->isFlagSet( x, y, z, fluid ) )
            {
               Cell currentPosition( x, y, z );
               blocks_->transformBlockLocalToGlobalCell( currentPosition, *block );

               // only consider the bottom row
               if(currentPosition[1] == uint_t(0))
               {
                  Vector3< real_t > velocity;
                  Vector3< real_t > vel_left;
                  Vector3< real_t > vel_right;
                  getVelocity( velocity, latticeModel, cell );
                  pdfField->getVelocity( vel_left, x-int_c(1), y, z );
                  pdfField->getVelocity( vel_right, x+int_c(1), y, z );
                  if( (vel_left[0] * vel_right[0]) < real_c(0) && velocity[0] > real_c(0))
                  {
                     real_t xPosition = blocks_->getBlockLocalCellCenter(*block, currentPosition)[0];
                     reattachmentLocations.push_back( (xPosition - stepSize_[0]) / stepSize_[1] );
                  }
               }
            }
         }
      }

      WALBERLA_ROOT_SECTION() {
         std::ofstream fileLocBottom(filename_.c_str(), std::ios_base::app);
         fileLocBottom << executionCounter_ << " ";
         for (auto i = reattachmentLocations.begin(); i != reattachmentLocations.end(); i++) {
            fileLocBottom << *i << " ";
         }
         fileLocBottom << std::endl;
         fileLocBottom.close();
      }
   }

private:

   shared_ptr< StructuredBlockStorage > blocks_;
   BlockDataID  pdfFieldID_;
   BlockDataID flagFieldID_;
   FlagUID fluid_;
   std::string filename_;
   uint_t checkFrequency_;
   Vector3< real_t > stepSize_;

   uint_t executionCounter_;
};


int main( int argc, char ** argv )
{
   walberla::Environment walberlaEnv( argc, argv );
   
   auto blocks = blockforest::createUniformBlockGridFromConfig( walberlaEnv.config() );

   // read parameters
   auto parameters = walberlaEnv.config()->getOneBlock( "Parameters" );
   auto domain     = walberlaEnv.config()->getOneBlock( "DomainSetup" );
   auto boundariesConfig = walberlaEnv.config()->getOneBlock( "Boundaries" );

   //! [Params]
   const real_t Re = parameters.getParameter< real_t >( "Re", real_c( 1000 ) );
   const real_t uLB = parameters.getParameter< real_t >( "uLB", real_c( 0.01 )  );
   //! [Params]
   //! [height]
   const Vector3<real_t> domainsize = domain.getParameter< Vector3<real_t> >( "cellsPerBlock", Vector3<real_t>() );
  //! [height]
   const uint_t timesteps = parameters.getParameter< uint_t > ( "timesteps", uint_c( 100000 ) );   

   const Vector3<real_t> inletVelocity(uLB, real_c(0), real_c(0));
   //!  [Omega]
   const real_t viscosity = (domainsize[1] * uLB) / Re;
   const real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);
   //!  [Omega] 
   WALBERLA_LOG_INFO_ON_ROOT( "Re =  " << Re );
   WALBERLA_LOG_INFO_ON_ROOT( "uLB =  " << uLB );
   WALBERLA_LOG_INFO_ON_ROOT( "timesteps =  " << timesteps );
   WALBERLA_LOG_INFO_ON_ROOT( "viscosity =  " << viscosity );
   WALBERLA_LOG_INFO_ON_ROOT( "omega =  " << omega );
				
   const real_t remainingTimeLoggerFrequency = parameters.getParameter< real_t >( "remainingTimeLoggerFrequency", real_c(3) ); // in seconds

   // create fields
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::SRT( omega ) );
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel, inletVelocity, real_t(1) );
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // create and initialize boundary handling
   const FlagUID fluidFlagUID( "Fluid" );


   using BHFactory = lbm::DefaultBoundaryHandlingFactory< LatticeModel_T, FlagField_T > ;

   BlockDataID boundaryHandlingID = BHFactory::addBoundaryHandlingToStorage( blocks, "boundary handling", flagFieldID, pdfFieldID, fluidFlagUID,
                                                                             inletVelocity,
                                                                             boundariesConfig.getParameter< Vector3<real_t> >( "velocity1", Vector3<real_t>() ),
                                                                             boundariesConfig.getParameter< real_t > ( "pressure0", real_c( 1.0 ) ),
                                                                             boundariesConfig.getParameter< real_t > ( "pressure1", real_c( 1.0 ) ) );
  //! [geomboundary]
   geometry::initBoundaryHandling<BHFactory::BoundaryHandling>( *blocks, boundaryHandlingID, boundariesConfig );
  //! [geomboundary]
   geometry::setNonBoundaryCellsToDomain<BHFactory::BoundaryHandling> ( *blocks, boundaryHandlingID );

   // create time loop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication( blocks );
   communication.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) );

   // add LBM sweep and communication to time loop
   timeloop.add() << BeforeFunction( communication, "communication" )
                  << Sweep( BHFactory::BoundaryHandling::getBlockSweep( boundaryHandlingID ), "boundary handling" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, fluidFlagUID ) ), "LB stream & collide" );

   //add ReattachmentLengthFinder logger

   //!  [Logger]
   std::string loggingFileName = "ReattachmentLengthLogging_Re" + std::to_string(uint_c(Re)) + ".txt";
   uint_t loggingFrequency = parameters.getParameter< uint_t >( "loggingFrequency", uint_c(1) );
   Vector3<real_t> stepSize = boundariesConfig.getOneBlock("Body").getParameter< Vector3<real_t> >("max", Vector3<real_t>() );

   timeloop.addFuncAfterTimeStep( ReattachmentLengthFinder( blocks, pdfFieldID, flagFieldID, fluidFlagUID,
                                                            loggingFileName, loggingFrequency, stepSize ),
                                  "reattachment length finder" );
   //!  [Logger]

   // log remaining time
   timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "remaining time logger" );

   // add VTK output to time loop
   lbm::VTKOutput< LatticeModel_T, FlagField_T >::addToTimeloop( timeloop, blocks, walberlaEnv.config(), pdfFieldID, flagFieldID, fluidFlagUID );

   timeloop.run();

   return EXIT_SUCCESS;
}
}

int main( int argc, char ** argv )
{
   walberla::main(argc, argv);
}
