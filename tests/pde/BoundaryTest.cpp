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
//! \file BoundaryTest.cpp
//! \ingroup pde
//! \author Dominik Bartuschat <dominik.bartuschat@fau.de>
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "boundary/BoundaryHandling.h"

#include "core/Abort.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"
#include "field/FlagField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

#include "pde/iterations/CGFixedStencilIteration.h"
#include "pde/iterations/CGIteration.h"
#include "pde/boundary/Dirichlet.h"
#include "pde/boundary/Neumann.h"

#include "stencil/D2Q5.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <cmath>

namespace walberla {



typedef GhostLayerField< real_t, 1 > Field_T;
typedef stencil::D2Q5                Stencil_T;
typedef pde::CGIteration<Stencil_T>::StencilField_T  StencilField_T;

typedef walberla::uint8_t      flag_t;
typedef FlagField < flag_t >   FlagField_T;
typedef pde::Dirichlet< Stencil_T, flag_t >  Dirichlet_T;
typedef pde::Neumann< Stencil_T, flag_t >  Neumann_T;

typedef BoundaryHandling< FlagField_T, Stencil_T, Dirichlet_T, Neumann_T > BoundaryHandling_T;


const FlagUID  Domain_Flag( "domain" );
const FlagUID  Dirichlet_Flag( "dirichlet" );
const FlagUID  Neumann_Flag( "neumann" );



///////////////////////
// BOUNDARY HANDLING //
///////////////////////

//**********************************************************************************************************************
/*!
*   \brief Functor that is used to add a boundary handling object to each block
*
*   The member function "BoundaryHandling_T * operator()( IBlock* const block ) const" is called by the framework in
*   order to add a boundary handling object of type 'BoundaryHandling_T' to each block.
*/
//**********************************************************************************************************************

class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & rhsFieldId, const BlockDataID & stencilFieldId, const BlockDataID & adaptBCStencilFieldId, const BlockDataID & flagFieldId, const StructuredBlockStorage& blocks ) :
      rhsFieldId_( rhsFieldId ), stencilFieldId_ ( stencilFieldId ), adaptBCStencilFieldId_ ( adaptBCStencilFieldId ), flagFieldId_ (flagFieldId), blockStorage_(blocks) {}

   BoundaryHandling_T * operator()( IBlock* const block ) const;

private:

   const BlockDataID rhsFieldId_;
   const BlockDataID stencilFieldId_;
   const BlockDataID adaptBCStencilFieldId_;
   const BlockDataID flagFieldId_;
   const StructuredBlockStorage & blockStorage_;

}; // class MyBoundaryHandling

BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block ) const
{

   Field_T * rhsField = block->getData< Field_T >( rhsFieldId_ );
   StencilField_T * stencilField = block->getData< StencilField_T >( stencilFieldId_ );
   StencilField_T * adaptStencilField = block->getData< StencilField_T >( adaptBCStencilFieldId_ );
   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId_ );

   const flag_t domain = flagField->registerFlag( Domain_Flag ); // register the fluid flag at the flag field

   // A new boundary handling instance that uses the just fetched flag field is created:
   // Additional, internal flags used by the boundary handling will be stored in this flag field.
   return new BoundaryHandling_T( "boundary handling", flagField, domain,
                                  Dirichlet_T( "dirichlet", Dirichlet_Flag, rhsField, stencilField, adaptStencilField, flagField ) ,
                                    Neumann_T( "neumann", Neumann_Flag, rhsField, stencilField, adaptStencilField, flagField, blockStorage_ )
   );
}



void initRHS( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & rhsId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      Field_T * rhs = block->getData< Field_T >( rhsId );
      CellInterval xyz = rhs->xyzSize();
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
      {
         const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );
         rhs->get( *cell ) = real_t(4) * math::pi * math::pi * std::sin( real_t(2) * math::pi * p[0] ) * std::sinh( real_t(2) * math::pi * p[1] );
      }
   }
}



void setRHSConstValue( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & rhsId, const real_t value )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      Field_T * rhs = block->getData< Field_T >( rhsId );
      CellInterval xyz = rhs->xyzSize();
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
      {
         rhs->get( *cell ) = value;
      }
   }
}



void setBoundaryConditionsDirichl( shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & boundaryHandlingId )
{
   CellInterval domainBBInGlobalCellCoordinates = blocks->getDomainCellBB();
   domainBBInGlobalCellCoordinates.expand(Cell(1,1,0));

   // Iterate through all blocks that are allocated on this process and part of the block structure 'blocks'
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {

      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingId );

      CellInterval domainBB( domainBBInGlobalCellCoordinates );
      blocks->transformGlobalToBlockLocalCellInterval( domainBB, *block );

      CellInterval west( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMin(), domainBB.yMax(), domainBB.zMax() );
      CellInterval east( domainBB.xMax(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
      CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );

      // Set north boundary to defined function
      for( auto cell = north.begin(); cell != north.end(); ++cell )
      {

         const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );
         real_t val = std::sin( real_t( 2 ) * math::pi * p[0] ) * std::sinh( real_t( 2 ) * math::pi * p[1] );

         boundaryHandling->forceBoundary( Dirichlet_Flag, cell->x(), cell->y(), cell->z(), pde::Dirichlet< Stencil_T, flag_t >::DirichletBC( val ) );

      }

      // Set all other boundaries to zero
      for( auto & ci : { west, east, south } )
      {
         boundaryHandling->forceBoundary( Dirichlet_Flag, ci, pde::Dirichlet< Stencil_T, flag_t >::DirichletBC( real_t( 0 ) ) );
      }

      // All the remaining cells need to be marked as being fluid. The 'fillWithDomain' call does just that.
      boundaryHandling->fillWithDomain( domainBB );
   }
}



void setBoundaryConditionsMixed( shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & boundaryHandlingId )
{
   CellInterval domainBBInGlobalCellCoordinates = blocks->getDomainCellBB();
   domainBBInGlobalCellCoordinates.expand(Cell(1,1,0));

   // Iterate through all blocks that are allocated on this process and part of the block structure 'blocks'
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {

      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingId );

      CellInterval domainBB( domainBBInGlobalCellCoordinates );
      blocks->transformGlobalToBlockLocalCellInterval( domainBB, *block );

      CellInterval west( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMin(), domainBB.yMax(), domainBB.zMax() );
      CellInterval east( domainBB.xMax(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
      CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );

      // Set north boundary to large value
      for( auto cell = north.begin(); cell != north.end(); ++cell )
      {

         const real_t val = real_t(2);
         boundaryHandling->forceBoundary( Dirichlet_Flag, cell->x(), cell->y(), cell->z(), pde::Dirichlet< Stencil_T, flag_t >::DirichletBC( val ) );

      }

      // Set south boundary to large value
      for( auto cell = south.begin(); cell != south.end(); ++cell )
      {

         const real_t val = real_t(-2);
         boundaryHandling->forceBoundary( Dirichlet_Flag, cell->x(), cell->y(), cell->z(), pde::Dirichlet< Stencil_T, flag_t >::DirichletBC( val ) );

      }

      // Set all other boundaries to homogeneous Neumann BCs
      for( auto & ci : { west, east} )
      {
         boundaryHandling->forceBoundary( Neumann_Flag, ci, pde::Neumann< Stencil_T, flag_t >::NeumannBC( real_t( 0 ) ) );
      }

      // All the remaining cells need to be marked as being fluid. The 'fillWithDomain' call does just that.
      boundaryHandling->fillWithDomain( domainBB );
   }
}



void copyWeightsToStencilField( const shared_ptr< StructuredBlockStorage > & blocks, const std::vector<real_t> & weights, const BlockDataID & stencilId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      StencilField_T * stencil = block->getData< StencilField_T >( stencilId );
      
      WALBERLA_FOR_ALL_CELLS_XYZ(stencil,
         for( auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir )
            stencil->get(x,y,z,dir.toIdx()) = weights[ dir.toIdx() ];
      );
   }
}




int main( int argc, char** argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   const uint_t processes = uint_c( MPIManager::instance()->numProcesses() );
   if( processes != uint_t(1) && processes != uint_t(4) && processes != uint_t(8) )
      WALBERLA_ABORT( "The number of processes must be equal to 1, 4, or 8!" );

   logging::Logging::printHeaderOnStream();
   WALBERLA_ROOT_SECTION() { logging::Logging::instance()->setLogLevel( logging::Logging::PROGRESS ); }

   bool shortrun = false;
   for( int i = 1; i < argc; ++i )
      if( std::strcmp( argv[i], "--shortrun" ) == 0 ) shortrun = true;

   const uint_t xBlocks = ( processes == uint_t(1) ) ? uint_t(1) : ( ( processes == uint_t(4) ) ? uint_t(2) : uint_t(4) );
   const uint_t yBlocks = ( processes == uint_t(1) ) ? uint_t(1) : uint_t(2);
   const uint_t xCells = ( processes == uint_t(1) ) ? uint_t(200) : ( ( processes == uint_t(4) ) ? uint_t(100) : uint_t(50) );
   const uint_t yCells = ( processes == uint_t(1) ) ? uint_t(100) : uint_t(50);
   const real_t xSize = real_t(2);
   const real_t ySize = real_t(1);
   const real_t dx = xSize / real_c( xBlocks * xCells + uint_t(1) );
   const real_t dy = ySize / real_c( yBlocks * yCells + uint_t(1) );
   auto blocks = blockforest::createUniformBlockGrid( math::AABB( real_t(0.5) * dx, real_t(0.5) * dy, real_t(0),
                                                                  xSize - real_t(0.5) * dx, ySize - real_t(0.5) * dy, dx ),
                                                      xBlocks, yBlocks, uint_t(1),
                                                      xCells, yCells, uint_t(1),
                                                      true,
                                                      false, false, false );

   BlockDataID solId = field::addToStorage< Field_T >( blocks, "sol", real_t(0), field::zyxf, uint_t(1) );
   BlockDataID rId = field::addToStorage< Field_T >( blocks, "r", real_t(0), field::zyxf, uint_t(1) );
   BlockDataID dId = field::addToStorage< Field_T >( blocks, "d", real_t(0), field::zyxf, uint_t(1) );
   BlockDataID zId = field::addToStorage< Field_T >( blocks, "z", real_t(0), field::zyxf, uint_t(1) );

   BlockDataID rhsId = field::addToStorage< Field_T >( blocks, "rhs", real_t(0), field::zyxf, uint_t(1) );

   initRHS( blocks, rhsId );

   blockforest::communication::UniformBufferedScheme< Stencil_T > synchronizeD( blocks );
   synchronizeD.addPackInfo( make_shared< field::communication::PackInfo< Field_T > >( dId ) );

   std::vector< real_t > weights( Stencil_T::Size );
   weights[ Stencil_T::idx[ stencil::C ] ] = real_t(2) / ( blocks->dx() * blocks->dx() ) + real_t(2) / ( blocks->dy() * blocks->dy() ) + real_t(4) * math::pi * math::pi;
   weights[ Stencil_T::idx[ stencil::N ] ] = real_t(-1) / ( blocks->dy() * blocks->dy() );
   weights[ Stencil_T::idx[ stencil::S ] ] = real_t(-1) / ( blocks->dy() * blocks->dy() );
   weights[ Stencil_T::idx[ stencil::E ] ] = real_t(-1) / ( blocks->dx() * blocks->dx() );
   weights[ Stencil_T::idx[ stencil::W ] ] = real_t(-1) / ( blocks->dx() * blocks->dx() );


   BlockDataID stencilId = field::addToStorage< StencilField_T >( blocks, "w" );
   BlockDataID stencilBCadaptedId = field::addToStorage< StencilField_T >( blocks, "w" );

   BlockDataID flagId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flagField" );

   copyWeightsToStencilField( blocks, weights, stencilId );
   copyWeightsToStencilField( blocks, weights, stencilBCadaptedId );

   BlockDataID boundaryHandlingId = blocks->addBlockData< BoundaryHandling_T >( MyBoundaryHandling( rhsId, stencilId, stencilBCadaptedId, flagId, *blocks ),
                                                                                "boundary handling" );

   // Test Dirichlet BCs //
   setBoundaryConditionsDirichl( blocks, boundaryHandlingId );

   
   SweepTimeloop timeloop( blocks, uint_t(2) );

   timeloop.add()
            << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingId ), "boundary handling" )
            << AfterFunction(
                     pde::CGIteration< Stencil_T >( blocks->getBlockStorage(), solId, rId, dId, zId, rhsId, stencilBCadaptedId,
                                                    shortrun ? uint_t( 10 ) : uint_t( 10000 ), synchronizeD, real_c( 1e-6 ) ), "CG iteration" );

   timeloop.singleStep();

   // Check values for Dirichlet BCs //
   if( !shortrun )
   {
      Cell cellNearBdry( 75, 2, 0 );
      real_t solNearBdry( real_c(-0.16347) );
      Cell cellNearBdryLrg( 24, 95, 0 );
      real_t solNearBdryLrg( real_c(201.47) );
      Cell cellDomCentr( 100, 50, 0 );
      real_t solDomCentr( real_c(0.37587) );

      for( auto block = blocks->begin(); block != blocks->end(); ++block )
      {
         Field_T * sol = block->getData < Field_T > ( solId );
         if( blocks->getBlockCellBB( *block ).contains( cellNearBdry ) )
         {
            blocks->transformGlobalToBlockLocalCell(cellNearBdry,*block);
            // WALBERLA_LOG_DEVEL( "Solution value at cell " << cellNearBdry << ": " << sol->get(cellNearBdry) );
            WALBERLA_CHECK_LESS( std::fabs( solNearBdry - sol->get( cellNearBdry ) ) / solNearBdry, 0.0001, "Invalid value in cell " << cellNearBdry );
         }
         if( blocks->getBlockCellBB( *block ).contains( cellNearBdryLrg ) )
         {
            blocks->transformGlobalToBlockLocalCell(cellNearBdryLrg,*block);
            // WALBERLA_LOG_DEVEL( "Solution value at cell " << cellNearBdryLrg << ": " << sol->get(cellNearBdryLrg) );
            WALBERLA_CHECK_LESS( std::fabs( solNearBdryLrg - sol->get( cellNearBdryLrg ) ) / solNearBdryLrg, 0.0001, "Invalid value in cell " << cellNearBdryLrg );
         }
         if( blocks->getBlockCellBB( *block ).contains( cellDomCentr ) )
         {
            blocks->transformGlobalToBlockLocalCell(cellDomCentr,*block);
            // WALBERLA_LOG_DEVEL( "Solution value at cell " << cellDomCentr << ": " << sol->get(cellDomCentr) );
            WALBERLA_CHECK_LESS( std::fabs( solDomCentr - sol->get( cellDomCentr ) ) / solDomCentr, 0.0001, "Invalid value in cell " << cellDomCentr );
         }
      }
   }
   else
   {
      Cell cellNearBdry( 75, 2, 0 );
      real_t solNearBdry( real_c(-0.008355) );
      Cell cellNearBdryLrg( 24, 95, 0 );
      real_t solNearBdryLrg( real_c(132.188) );
      Cell cellDomCentr( 100, 50, 0 );
      real_t solDomCentr( 0.017603f );

      for( auto block = blocks->begin(); block != blocks->end(); ++block )
      {
         Field_T * sol = block->getData < Field_T > ( solId );
         if( blocks->getBlockCellBB( *block ).contains( cellNearBdry ) )
         {
             blocks->transformGlobalToBlockLocalCell(cellNearBdry,*block);
             // WALBERLA_LOG_DEVEL( "Solution value at cell " << cellNearBdry << ": " << sol->get(cellNearBdry) );
             WALBERLA_CHECK_LESS( std::fabs( solNearBdry - sol->get( cellNearBdry ) ) / solNearBdry, 0.0001, "Invalid value in cell " << cellNearBdry );
         }
         if( blocks->getBlockCellBB( *block ).contains( cellNearBdryLrg ) )
         {
             blocks->transformGlobalToBlockLocalCell(cellNearBdryLrg,*block);
             // WALBERLA_LOG_DEVEL( "Solution value at cell " << cellNearBdryLrg << ": " << sol->get(cellNearBdryLrg) );
             WALBERLA_CHECK_LESS( std::fabs( solNearBdryLrg - sol->get( cellNearBdryLrg ) ) / solNearBdryLrg, 0.0001, "Invalid value in cell " << cellNearBdryLrg );
         }
         if( blocks->getBlockCellBB( *block ).contains( cellDomCentr ) )
         {
             blocks->transformGlobalToBlockLocalCell(cellDomCentr,*block);
             // WALBERLA_LOG_DEVEL( "Solution value at cell " << cellDomCentr << ": " << sol->get(cellDomCentr) );
             WALBERLA_CHECK_LESS( std::fabs( solDomCentr - sol->get( cellDomCentr ) ) / solDomCentr, 0.0001, "Invalid value in cell " << cellDomCentr );
         }
      }

   }

   if( !shortrun )
   {
      vtk::writeDomainDecomposition( blocks );
      field::createVTKOutput< Field_T >( solId, *blocks, "solution_Dirich" )();
   }


   // Test Mixed BCs and re-setting BCs //
   setBoundaryConditionsMixed( blocks, boundaryHandlingId );
   // set RHS to zero to get 'ramp' as solution
   setRHSConstValue( blocks, rhsId, real_t(0));

   timeloop.singleStep();

   // Check values for mixed  BCs //
   // TODO

   if( !shortrun )
   {
      field::createVTKOutput< Field_T >( solId, *blocks, "solution_Mixed" )();
//      field::createVTKOutput< Field_T >( rhsId, *blocks, "rhs_Mixed" )();
//      field::createVTKOutput< StencilField_T >( stencilId, *blocks, "originalStencils_Mixed" )();
//      field::createVTKOutput< StencilField_T >( stencilBCadaptedId, *blocks, "adaptedStencils_Mixed" )();
   }

   logging::Logging::printFooterOnStream();
   return EXIT_SUCCESS;
}
}

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
