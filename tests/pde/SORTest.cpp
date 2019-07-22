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
//! \file SORTest.cpp
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Abort.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

#include "pde/ResidualNorm.h"
#include "pde/ResidualNormStencilField.h"
#include "pde/iterations/RBGSIteration.h"
#include "pde/sweeps/SORFixedStencil.h"
#include "pde/sweeps/SOR.h"

#include "stencil/D2Q5.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <cmath>

namespace walberla {



typedef GhostLayerField< real_t, 1 > PdeField_T;
using Stencil_T = stencil::D2Q5;
using StencilField_T = pde::SOR<Stencil_T>::StencilField_T;



void initU( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & uId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      if( blocks->atDomainYMaxBorder( *block ) )
      {
         PdeField_T * u = block->getData< PdeField_T >( uId );
         CellInterval xyz = u->xyzSizeWithGhostLayer();
         xyz.yMin() = xyz.yMax();
         for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
         {
            const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );
            u->get( *cell ) = std::sin( real_t(2) * math::pi * p[0] ) * std::sinh( real_t(2) * math::pi * p[1] );
         }
      }
   }
}



void initF( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & fId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      PdeField_T * f = block->getData< PdeField_T >( fId );
      CellInterval xyz = f->xyzSize();
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
      {
         const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );
         f->get( *cell ) = real_t(4) * math::pi * math::pi * std::sin( real_t(2) * math::pi * p[0] ) * std::sinh( real_t(2) * math::pi * p[1] );
      }
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



template <typename Field_T>
void clearField( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & fieldId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      block->getData< Field_T >( fieldId )->set( typename Field_T::value_type() );
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

   BlockDataID uId = field::addToStorage< PdeField_T >( blocks, "u", real_t(0), field::zyxf, uint_t(1) );

   initU( blocks, uId );

   BlockDataID fId = field::addToStorage< PdeField_T >( blocks, "f", real_t(0), field::zyxf, uint_t(1) );

   initF( blocks, fId );
   
   real_t omega = real_t(1.9);

   SweepTimeloop timeloop( blocks, uint_t(1) );

   blockforest::communication::UniformBufferedScheme< Stencil_T > communication( blocks );
   communication.addPackInfo( make_shared< field::communication::PackInfo< PdeField_T > >( uId ) );

   std::vector< real_t > weights( Stencil_T::Size );
   weights[ Stencil_T::idx[ stencil::C ] ] = real_t(2) / ( blocks->dx() * blocks->dx() ) + real_t(2) / ( blocks->dy() * blocks->dy() ) + real_t(4) * math::pi * math::pi;
   weights[ Stencil_T::idx[ stencil::N ] ] = real_t(-1) / ( blocks->dy() * blocks->dy() );
   weights[ Stencil_T::idx[ stencil::S ] ] = real_t(-1) / ( blocks->dy() * blocks->dy() );
   weights[ Stencil_T::idx[ stencil::E ] ] = real_t(-1) / ( blocks->dx() * blocks->dx() );
   weights[ Stencil_T::idx[ stencil::W ] ] = real_t(-1) / ( blocks->dx() * blocks->dx() );

   auto SORFixedSweep = pde::SORFixedStencil< Stencil_T >( blocks, uId, fId, weights, omega );

   timeloop.addFuncBeforeTimeStep( pde::RBGSIteration( blocks->getBlockStorage(), shortrun ? uint_t(10) : uint_t(10000),
                                                       communication,
                                                       SORFixedSweep.getRedSweep(), SORFixedSweep.getBlackSweep(),
                                                       pde::ResidualNorm< Stencil_T >( blocks->getBlockStorage(), uId, fId, weights ),
                                                       real_c(1e-6), uint_t(100) ), "SOR iteration" );

   timeloop.run();
   
   // rerun the test with a stencil field
   
   clearField<PdeField_T>( blocks, uId);
   initU( blocks, uId );
   
   BlockDataID stencilId = field::addToStorage< StencilField_T >( blocks, "w" );
   
   SweepTimeloop timeloop2( blocks, uint_t(1) );
   
   copyWeightsToStencilField( blocks, weights, stencilId );
   
   auto SORSweep = pde::SOR< Stencil_T >( blocks, uId, fId, stencilId, omega );
   
   timeloop2.addFuncBeforeTimeStep( pde::RBGSIteration( blocks->getBlockStorage(), shortrun ? uint_t(10) : uint_t(10000),
                                                        communication,
                                                        SORSweep.getRedSweep(), SORSweep.getBlackSweep(),
                                                        pde::ResidualNormStencilField< Stencil_T >( blocks->getBlockStorage(), uId, fId, stencilId ),
                                                        real_c(1e-6), uint_t(100) ), "SOR iteration" );
   
   timeloop2.run();

   if( !shortrun )
   {
      vtk::writeDomainDecomposition( blocks );
      field::createVTKOutput< PdeField_T >( uId, *blocks, "solution" )();
   }

   logging::Logging::printFooterOnStream();
   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
