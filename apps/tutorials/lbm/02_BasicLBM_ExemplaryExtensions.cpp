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
//! \file 02_BasicLBM_ExemplaryExtensions.cpp
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "core/all.h"
#include "blockforest/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "gui/all.h"
#include "lbm/all.h"
#include "timeloop/all.h"



namespace walberla {

using LatticeModel_T = lbm::D2Q9<lbm::collision_model::SRT, false, lbm::force_model::SimpleConstant>;

using Stencil_T = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;

using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;



//**********************************************************************************************************************
/*!
*   Class for determining (and logging) the minimum and maximum velocity as well as the minimum and maximum density.
*/
//**********************************************************************************************************************
class VelDensLogger
{
public:

   VelDensLogger( const std::string & filename, const shared_ptr< StructuredBlockStorage > & blocks,
                  const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const FlagUID & fluid,
                  const uint_t checkFrequency, bool includeImages = false ) :
      executionCounter_( uint_c(0) ), checkFrequency_( checkFrequency ),
      filename_( filename ), includeImages_( includeImages ),
      blocks_( blocks ), pdfFieldId_( pdfFieldId ), flagFieldId_( flagFieldId ), fluid_( fluid )
   {
      // open and close file on root process - purpose: If the file already exists, its content will be erased.
      WALBERLA_ROOT_SECTION()
      {
         std::string txtFile = filename_ + ".txt";
         std::ofstream file( txtFile.c_str() );
         file.close();
      }
   }

   void operator()()
   {
      // Every time this operator is called, the execution counter is increased.
      // The operator is only executed the first time it is called, and every "checkFrequency_"-th time after that.
      ++executionCounter_;
      if( ( executionCounter_ - uint_c(1) ) % checkFrequency_ != 0 )
         return;

      // data type for storing a position together with a floating point value
      using PosValuePair = std::pair<Cell, real_t>;

      // variables for storing the process local minimum/maximum velocity & density
      PosValuePair maxVelocity = std::pair< Cell, real_t >( Cell(), real_c(0) );
      PosValuePair minVelocity = std::pair< Cell, real_t >( Cell(), real_c(0) );
      PosValuePair maxDensity  = std::pair< Cell, real_t >( Cell(), real_c(0) );
      PosValuePair minDensity  = std::pair< Cell, real_t >( Cell(), real_c(0) );

      // iterate all blocks stored locally on this process
      for( auto block = blocks_->begin(); block != blocks_->end(); ++block )
      {
         // retrieve the pdf field and the flag field from the block
         PdfField_T  *  pdfField = block->getData< PdfField_T >( pdfFieldId_ );
         FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId_ );

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
               // ... temporarily save the location of the cell in global cell coordinates ...
               Cell currentPosition( x, y, z );
               blocks_->transformBlockLocalToGlobalCell( currentPosition, *block );

               // ... and get the velocity and density of the cell.
               Vector3< real_t > velocity;
               real_t density = getDensityAndVelocity( velocity, latticeModel, cell );

               // If this is the first cell we check, we initialize the variables that
               // store the process local minimum/maximum velocity & density.
               if( block == blocks_->begin() && cell == pdfField->beginXYZ() )
               {
                  maxVelocity.first  = currentPosition;
                  maxVelocity.second = velocity.length();
                  minVelocity.first  = currentPosition;
                  minVelocity.second = velocity.length();

                  maxDensity.first  = currentPosition;
                  maxDensity.second = density;
                  minDensity.first  = currentPosition;
                  minDensity.second = density;
               }
               // If this is not the first cell, we check if we have to update any of these variables.
               else
               {
                  real_t velMagnitude = velocity.length();
                  if( velMagnitude > maxVelocity.second )
                  {
                     maxVelocity.first  = currentPosition;
                     maxVelocity.second = velMagnitude;
                  }
                  if( velMagnitude < minVelocity.second )
                  {
                     minVelocity.first  = currentPosition;
                     minVelocity.second = velMagnitude;
                  }
                  if( density > maxDensity.second )
                  {
                     maxDensity.first  = currentPosition;
                     maxDensity.second = density;
                  }
                  if( density < minDensity.second )
                  {
                     minDensity.first  = currentPosition;
                     minDensity.second = density;
                  }
               }
            }
         }
      }

      // If no image must be created ...
      if( !includeImages_ )
      {
         // ... the simulation global minimum/maximum velocity & density are calculated by performing MPI reductions.
         mpi::reduceInplace( maxVelocity.second, mpi::MAX );
         mpi::reduceInplace( minVelocity.second, mpi::MIN );
         mpi::reduceInplace( maxDensity.second, mpi::MAX );
         mpi::reduceInplace( minDensity.second, mpi::MIN );

         // Finally, these simulation global values are written to file (only by
         // the root process that holds the results of the MPI reductions).
         WALBERLA_ROOT_SECTION()
         {
            std::string txtFile = filename_ + ".txt";
            std::ofstream file( txtFile.c_str(), std::ios_base::app );
            file << executionCounter_ << " "
                 << maxVelocity.second << " " << minVelocity.second << " "
                 << maxDensity.second << " " << minDensity.second << std::endl;
            file.close();
         }
      }
      // If in addition to the information being written to file an image must be created ...
      else
      {
         // ... all the relevant data (min/max velocity & density WITH additional position information)
         // is gathered on the root process.
         std::vector< cell_idx_t > xPosMaxVelocity = mpi::gather( maxVelocity.first.x() );
         std::vector< cell_idx_t > yPosMaxVelocity = mpi::gather( maxVelocity.first.y() );
         std::vector< real_t >     magMaxVelocity  = mpi::gather( maxVelocity.second );

         std::vector< cell_idx_t > xPosMinVelocity = mpi::gather( minVelocity.first.x() );
         std::vector< cell_idx_t > yPosMinVelocity = mpi::gather( minVelocity.first.y() );
         std::vector< real_t >     magMinVelocity  = mpi::gather( minVelocity.second );

         std::vector< cell_idx_t > xPosMaxDensity = mpi::gather( maxDensity.first.x() );
         std::vector< cell_idx_t > yPosMaxDensity = mpi::gather( maxDensity.first.y() );
         std::vector< real_t >     maxDensityVec  = mpi::gather( maxDensity.second );

         std::vector< cell_idx_t > xPosMinDensity = mpi::gather( minDensity.first.x() );
         std::vector< cell_idx_t > yPosMinDensity = mpi::gather( minDensity.first.y() );
         std::vector< real_t >     minDensityVec  = mpi::gather( minDensity.second );

         WALBERLA_ROOT_SECTION()
         {
            uint_t maxVelocityIndex = uint_c(0);
            uint_t minVelocityIndex = uint_c(0);
            uint_t maxDensityIndex = uint_c(0);
            uint_t minDensityIndex = uint_c(0);

            WALBERLA_ASSERT_EQUAL( xPosMaxVelocity.size(), yPosMaxVelocity.size() );
            WALBERLA_ASSERT_EQUAL( xPosMaxVelocity.size(), magMaxVelocity.size() );
            WALBERLA_ASSERT_EQUAL( xPosMaxVelocity.size(), xPosMinVelocity.size() );
            WALBERLA_ASSERT_EQUAL( xPosMaxVelocity.size(), yPosMinVelocity.size() );
            WALBERLA_ASSERT_EQUAL( xPosMaxVelocity.size(), magMinVelocity.size() );
            WALBERLA_ASSERT_EQUAL( xPosMaxVelocity.size(), xPosMaxDensity.size() );
            WALBERLA_ASSERT_EQUAL( xPosMaxVelocity.size(), yPosMaxDensity.size() );
            WALBERLA_ASSERT_EQUAL( xPosMaxVelocity.size(), maxDensityVec.size() );
            WALBERLA_ASSERT_EQUAL( xPosMaxVelocity.size(), xPosMinDensity.size() );
            WALBERLA_ASSERT_EQUAL( xPosMaxVelocity.size(), yPosMinDensity.size() );
            WALBERLA_ASSERT_EQUAL( xPosMaxVelocity.size(), minDensityVec.size() );

            // The root process determines the global minimum/maximum of the velocity & density.
            for( uint_t i = uint_c(1); i != magMaxVelocity.size(); ++i )
            {
               if( magMaxVelocity[i] > magMaxVelocity[maxVelocityIndex] )
                  maxVelocityIndex = i;
               if( magMinVelocity[i] < magMinVelocity[minVelocityIndex] )
                  minVelocityIndex = i;
               if( maxDensityVec[i] > maxDensityVec[maxDensityIndex] )
                  maxDensityIndex = i;
               if( minDensityVec[i] < minDensityVec[minDensityIndex] )
                  minDensityIndex = i;
            }

            // These simulation global values are written to file (only by
            // the root process that now holds the results).
            std::string txtFile = filename_ + ".txt";
            std::ofstream file( txtFile.c_str(), std::ios_base::app );
            file << executionCounter_ << " "
                 << magMaxVelocity[maxVelocityIndex] << " " << magMinVelocity[minVelocityIndex] << " "
                 << maxDensityVec[maxDensityIndex] << " " << minDensityVec[minDensityIndex] << std::endl;
            file.close();

            // Finally, two images are created: one for displaying the location of the minimum and maximum velocity,
            // and one for displaying the location of the minimum and maximum density.
            geometry::RGBAImage velImage( blocks_->getNumberOfXCells(), blocks_->getNumberOfYCells() );
            cell_idx_t yPosInvertion = cell_idx_c( blocks_->getNumberOfYCells() ) - cell_idx_c(1);
            // max -> red (= green and blue to zero)
            velImage.setElement( xPosMaxVelocity[maxVelocityIndex], yPosInvertion - yPosMaxVelocity[maxVelocityIndex], geometry::RGBAImage::G, real_c(0) );
            velImage.setElement( xPosMaxVelocity[maxVelocityIndex], yPosInvertion - yPosMaxVelocity[maxVelocityIndex], geometry::RGBAImage::B, real_c(0) );
            // min -> green (= red and blue to zero)
            velImage.setElement( xPosMinVelocity[minVelocityIndex], yPosInvertion - yPosMinVelocity[minVelocityIndex], geometry::RGBAImage::R, real_c(0) );
            velImage.setElement( xPosMinVelocity[minVelocityIndex], yPosInvertion - yPosMinVelocity[minVelocityIndex], geometry::RGBAImage::B, real_c(0) );
            std::ostringstream filename;
            filename << filename_ << "_" << executionCounter_ << "_velocity.png";
            velImage.save( filename.str() );

            geometry::RGBAImage denImage( blocks_->getNumberOfXCells(), blocks_->getNumberOfYCells() );
            // max -> red
            denImage.setElement( xPosMaxDensity[maxDensityIndex], yPosInvertion - yPosMaxDensity[maxDensityIndex], geometry::RGBAImage::G, real_c(0) );
            denImage.setElement( xPosMaxDensity[maxDensityIndex], yPosInvertion - yPosMaxDensity[maxDensityIndex], geometry::RGBAImage::B, real_c(0) );
            // min -> green
            denImage.setElement( xPosMinDensity[minDensityIndex], yPosInvertion - yPosMinDensity[minDensityIndex], geometry::RGBAImage::R, real_c(0) );
            denImage.setElement( xPosMinDensity[minDensityIndex], yPosInvertion - yPosMinDensity[minDensityIndex], geometry::RGBAImage::B, real_c(0) );
            filename.str("");
            filename << filename_ << "_" << executionCounter_ << "_density.png";
            denImage.save( filename.str() );
         }
      }
   }

private:

   uint_t executionCounter_;
   uint_t checkFrequency_;

   std::string filename_;
   bool        includeImages_;

   shared_ptr< StructuredBlockStorage > blocks_;
   BlockDataID  pdfFieldId_;
   BlockDataID flagFieldId_;

   FlagUID fluid_;
};






//**********************************************************************************************************************
/*!
*   Class for calculating (and logging) the force exerted by the fluid on an obstacle inside the domain.
*/
//**********************************************************************************************************************
class ForceLogger
{
public:

   ForceLogger( const std::string & filename, const shared_ptr< StructuredBlockStorage > & blocks,
                const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const FlagUID & fluid, const FlagUID & obstacle,
                const uint_t checkFrequency ) :
      executionCounter_( uint_c(0) ), checkFrequency_( checkFrequency ), filename_( filename ), blocks_( blocks ), pdfFieldId_( pdfFieldId ), flagFieldId_( flagFieldId ),
      fluid_( fluid ), obstacle_( obstacle )
   {
      // open and close file on root process - purpose: If the file already exists, its content will be erased.
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file( filename.c_str() );
         file.close();
      }
   }

   void operator()()
   {
      // Every time this operator is called, the execution counter is increased.
      // The operator is only executed the first time it is called, and every "checkFrequency_"-th time after that.
      ++executionCounter_;
      if( ( executionCounter_ - uint_c(1) ) % checkFrequency_ != 0 )
         return;

      // retrieve the global cell bounding box of the entire domain
      const auto & domainCellBB = blocks_->getDomainCellBB();

      // variable for storing the accumulated force acting on the object(s) inside the domain
      Vector3< real_t > force;

      // iterate all blocks stored locally on this process
      for( auto block = blocks_->begin(); block != blocks_->end(); ++block )
      {
         // retrieve the pdf field and the flag field from the block
         PdfField_T  *  pdfField = block->getData< PdfField_T >( pdfFieldId_ );
         FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId_ );

         // get the flag that marks a cell as being fluid & the flag that marks a cell as being an obstacle
         auto fluid    = flagField->getFlag( fluid_ );
         auto obstacle = flagField->getFlag( obstacle_ );

         // domain cell bounding box in block local coordinates
         CellInterval globalDomainInLocalCellCoordinates;
         blocks_->transformGlobalToBlockLocalCellInterval( globalDomainInLocalCellCoordinates, *block, domainCellBB );

         // get the block local cell bounding box
         auto xyzSize = pdfField->xyzSize();
         WALBERLA_ASSERT_EQUAL( xyzSize, flagField->xyzSize() );

         // iterate all cells of the block using the just retrieved block local cell bounding box
         for( auto z = xyzSize.zMin(); z <= xyzSize.zMax(); ++z ) {
            for( auto y = xyzSize.yMin(); y <= xyzSize.yMax(); ++y ) {
               for( auto x = xyzSize.xMin(); x <= xyzSize.xMax(); ++x )
               {
                  // It the current cell is a fluid cell ...
                  if( flagField->isFlagSet( x, y, z, fluid ) )
                  {
                     // ... check all neighboring cells ...
                     for( auto it = Stencil_T::beginNoCenter(); it != Stencil_T::end(); ++it )
                     {
                        auto nx = x + cell_idx_c( it.cx() ); // the block local cell
                        auto ny = y + cell_idx_c( it.cy() ); // index of the
                        auto nz = z + cell_idx_c( it.cz() ); // neighboring cell

                        // ... and if the neighbor is inside the domain (otherwise we might also accumulate the force acting on the domain border!)
                        //     and marked as being an obstacle ...
                        if( globalDomainInLocalCellCoordinates.contains( nx, ny, nz ) && flagField->isFlagSet( nx, ny, nz, obstacle ) )
                        {
                           // ... we can calculate the acting force using the momentum exchange method.
                           const real_t f = real_c(2) * pdfField->get( x, y, z, it.toIdx() );
                           force[0] += real_c( it.cx() ) * f;
                           force[1] += real_c( it.cy() ) * f;
                           force[2] += real_c( it.cz() ) * f;
                        }
                     }
                  }
               }
            }
         }
      }

      // The force is accumulated over all processes. The result is stored only on the root process.
      mpi::reduceInplace( force, mpi::SUM );

      // The root process writes the accumulated force to file.
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file( filename_.c_str(), std::ios_base::app );
         file << executionCounter_ << " " << force[0] << " " << force[1] << " " << force[2] << std::endl;
         file.close();
      }
   }

private:

   uint_t executionCounter_;
   uint_t checkFrequency_;

   std::string filename_;

   shared_ptr< StructuredBlockStorage > blocks_;
   BlockDataID  pdfFieldId_;
   BlockDataID flagFieldId_;

   FlagUID fluid_;
   FlagUID obstacle_;
};






//**********************************************************************************************************************
/*!
*   LBGK/SRT algorithm of the lattice Boltzmann method extended by an additional, simple force term (by Luo).
*
*   Please note: naive, straight-forward, completely unoptimized implementation!
*/
//**********************************************************************************************************************

class SRTStream
{
public:

   SRTStream( const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const FlagUID & fluid ) :
      pdfFieldId_( pdfFieldId ), flagFieldId_( flagFieldId ), fluid_( fluid ) {}

   // This operator can be registered as a "sweep" at the framework and will be executed for each process local block
   void operator()( IBlock * const block )
   {
      // retrieve the pdf field and the flag field from the block
      PdfField_T  *  pdfField = block->getData< PdfField_T >( pdfFieldId_ );
      FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId_ );

      // construct a temporary, second pdf field that will be used for streaming
      // by cloning the pdf field stored at the block
      shared_ptr< PdfField_T > tmpPdfField( pdfField->clone() );

      // get the flag that marks a cell as being fluid
      auto fluid = flagField->getFlag( fluid_ );

      // get the block local cell bounding box
      auto xyzSize = pdfField->xyzSize();
      WALBERLA_ASSERT_EQUAL( xyzSize, flagField->xyzSize() );

      // iterate all cells of the block using the just retrieved block local cell bounding box
      for( auto z = xyzSize.zMin(); z <= xyzSize.zMax(); ++z ) {
         for( auto y = xyzSize.yMin(); y <= xyzSize.yMax(); ++y ) {
            for( auto x = xyzSize.xMin(); x <= xyzSize.xMax(); ++x )
            {
               // It the current cell is a fluid cell, the stream process can be executed.
               if( flagField->isFlagSet( x, y, z, fluid ) )
               {
                  // stream pull:
                  // All relevant pdf values from neighboring cells are pulled into the current cell.
                  // To prevent pdf values from being overwritten, values are fetched from the pdf field
                  // of the block and stored in the temporary, function local pdf field.
                  for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
                     tmpPdfField->get( x, y, z, d.toIdx() ) = pdfField->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );
               }
            }
         }
      }

      // Since the pdf field stored at the block was only read from and streamed values are stored in the function local, temporary pdf field,
      // the data of these two fields must be swapped. After this operation, the post-stream state is stored in the pdf field of the block.
      pdfField->swapDataPointers( tmpPdfField.get() );
   }

private:

   BlockDataID  pdfFieldId_;
   BlockDataID flagFieldId_;

   FlagUID fluid_;
};

class SRTCollideForce
{
public:

   SRTCollideForce( const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const FlagUID & fluid,
             const real_t omega, const Vector3< real_t > & force ) :
      pdfFieldId_( pdfFieldId ), flagFieldId_( flagFieldId ), fluid_( fluid ), omega_( omega ), force_( force ) {}

   // This operator can be registered as a "sweep" at the framework and will be executed for each process local block
   void operator()( IBlock * const block )
   {
      // retrieve the pdf field and the flag field from the block
      PdfField_T  *  pdfField = block->getData< PdfField_T >( pdfFieldId_ );
      FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId_ );

      // get the flag that marks a cell as being fluid
      auto fluid = flagField->getFlag( fluid_ );

      // get the block local cell bounding box
      auto xyzSize = pdfField->xyzSize();
      WALBERLA_ASSERT_EQUAL( xyzSize, flagField->xyzSize() );

      // iterate all cells of the block using the just retrieved block local cell bounding box
      for( auto z = xyzSize.zMin(); z <= xyzSize.zMax(); ++z ) {
         for( auto y = xyzSize.yMin(); y <= xyzSize.yMax(); ++y ) {
            for( auto x = xyzSize.xMin(); x <= xyzSize.xMax(); ++x )
            {
               // It the current cell is a fluid cell, the collide process can be executed.
               if( flagField->isFlagSet( x, y, z, fluid ) )
               {
                  // The collision process requires the velocity and density (required for calculating the
                  // equilibrium distribution) of the current cell.
                  Vector3<real_t> velocity;
                  real_t rho = pdfField->getDensityAndEquilibriumVelocity( velocity, x, y, z );

                  // cell local collision
                  for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
                  {
                     // simple, additional force term by Luo
                     const real_t force_trm = real_c(3.0) * LatticeModel_T::w[ d.toIdx() ] *
                              ( real_c(d.cx()) * force_[0] + real_c(d.cy()) * force_[1] + real_c(d.cz()) * force_[2] );

                     // collision (LBGK/SRT) with additional force term
                     pdfField->get( x, y, z, d.toIdx() ) = ( real_c(1.0) - omega_ ) * pdfField->get( x, y, z, d.toIdx() ) +
                                                                           omega_   * lbm::EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, rho ) +
                                                           force_trm;
                  }
               }
            }
         }
      }
   }

private:

   BlockDataID  pdfFieldId_;
   BlockDataID flagFieldId_;

   FlagUID fluid_;

   real_t omega_;
   Vector3< real_t > force_;
};






//**********************************************************************************************************************
/*!
*   MAIN
*/
//**********************************************************************************************************************
int main( int argc, char ** argv )
{
   walberla::Environment walberlaEnv( argc, argv );

   auto blocks = blockforest::createUniformBlockGridFromConfig( walberlaEnv.config() );

   // read parameters
   auto parameters = walberlaEnv.config()->getOneBlock( "Parameters" );

   const real_t          omega           = parameters.getParameter< real_t >         ( "omega",           real_c( 1.4 ) );
   const Vector3<real_t> initialVelocity = parameters.getParameter< Vector3<real_t> >( "initialVelocity", Vector3<real_t>() );
   const uint_t          timesteps       = parameters.getParameter< uint_t >         ( "timesteps",       uint_c( 10 )  );

   const double remainingTimeLoggerFrequency = parameters.getParameter< double >( "remainingTimeLoggerFrequency", 3.0 ); // in seconds

   // create lattice model

   const Vector3< real_t > force = parameters.getParameter< Vector3<real_t> >( "force", Vector3<real_t>() );
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::SRT( omega ), lbm::force_model::SimpleConstant( force ) );

   // create fields

   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel, initialVelocity, real_t(1) );
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // create and initialize boundary handling
   const FlagUID fluidFlagUID( "Fluid" );

   auto boundariesConfig = walberlaEnv.config()->getOneBlock( "Boundaries" );

   using BHFactory = lbm::DefaultBoundaryHandlingFactory<LatticeModel_T, FlagField_T>;
   BlockDataID boundaryHandlingId = BHFactory::addBoundaryHandlingToStorage( blocks, "boundary handling", flagFieldId, pdfFieldId, fluidFlagUID,
                                                                             boundariesConfig.getParameter< Vector3<real_t> >( "velocity0", Vector3<real_t>() ),
                                                                             boundariesConfig.getParameter< Vector3<real_t> >( "velocity1", Vector3<real_t>() ),
                                                                             boundariesConfig.getParameter< real_t > ( "pressure0", real_c( 1.0 ) ),
                                                                             boundariesConfig.getParameter< real_t > ( "pressure1", real_c( 1.0 ) ) );

   geometry::initBoundaryHandling<BHFactory::BoundaryHandling>( *blocks, boundaryHandlingId, boundariesConfig );
   geometry::setNonBoundaryCellsToDomain<BHFactory::BoundaryHandling> ( *blocks, boundaryHandlingId );

   // create time loop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // add velocity/density logger
   timeloop.addFuncBeforeTimeStep( VelDensLogger( "maxmin", blocks, pdfFieldId, flagFieldId, fluidFlagUID,
                                                  parameters.getParameter< uint_t >( "velDenLoggerFrequency", uint_c(1) ),
                                                  parameters.getParameter< bool >( "velDenIncludeImages", false ) ), "min/max logger" );

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication( blocks );
   communication.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldId ) );

   // add LBM algorithm, force evaluation (must be post-collide/pre-stream!), and communication to time loop
   timeloop.add() << Sweep( SRTCollideForce( pdfFieldId,flagFieldId, fluidFlagUID, omega, force ), "LB collide" );
   timeloop.add() << BeforeFunction( communication, "communication" )
                  << Sweep( BHFactory::BoundaryHandling::getBlockSweep( boundaryHandlingId ), "boundary handling" );
   timeloop.add() << BeforeFunction( ForceLogger( "force.txt", blocks, pdfFieldId, flagFieldId, fluidFlagUID, BHFactory::getNoSlip(),
                                                   parameters.getParameter< uint_t >( "forceLoggerFrequency", uint_c(1) ) ), "force logger" )
                  << Sweep( SRTStream( pdfFieldId,flagFieldId, fluidFlagUID ), "LB stream" );

   // LBM stability check
   timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< PdfField_T, FlagField_T >( walberlaEnv.config(), blocks, pdfFieldId,
                                                                                                             flagFieldId, fluidFlagUID ) ),
                                  "LBM stability check" );

   // log remaining time
   timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "remaining time logger" );

   // add VTK output to time loop
   lbm::VTKOutput< LatticeModel_T, FlagField_T >::addToTimeloop( timeloop, blocks, walberlaEnv.config(), pdfFieldId, flagFieldId, fluidFlagUID );


   // create adaptors, so that the GUI also displays density and velocity
   // adaptors are like fields with the difference that they do not store values
   // but calculate the values based on other fields ( here the PdfField )
   field::addFieldAdaptor<lbm::Adaptor<LatticeModel_T>::Density>       ( blocks, pdfFieldId, "DensityAdaptor" );
   field::addFieldAdaptor<lbm::Adaptor<LatticeModel_T>::VelocityVector>( blocks, pdfFieldId, "VelocityAdaptor" );

   if ( parameters.getParameter<bool>( "useGui", false ) )
   {
      GUI gui ( timeloop, blocks, argc, argv );
      lbm::connectToGui<LatticeModel_T> ( gui );
      gui.run();
   }
   else
      timeloop.run();

   return EXIT_SUCCESS;
}
}

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}