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
//! \file DiffusionDirichlet.cpp
//! \ingroup lbm
//! \author Matthias Markl <matthias.markl@fau.de>
//! \brief Calculates temperature profile between two plates with a cosinuidal temperature profile
//!
//! Configuration:
//!  -l  channel length (x and y size)
//!  -w  channel width  (z size)
//!  -o  relaxation parameter -> depending on thermal diffusivity k: omega = 1/(3.k+0.5)
//!  -v  velocity in channel along x direction
//!  -e  E^2 error value
//!  -t  amount of timesteps to simulate
//! Setup:
//!  - periodic in x and z coordinate
//!  - plates in the x-z planes
//!  - no hydrodynamic computations( constant velocity )
//!  - comparison with analytical solution in the last timestep
//
// @see Like, L. et. al.: Boundary Conditions for thermal lattice Boltzmann equation method
//
//======================================================================================================================

#include "lbm/boundary/factories/DefaultDiffusionBoundaryHandling.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/gui/Connection.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/Velocity.h"

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/cell/CellSet.h"
#include "core/cell/CellVector.h"
#include "core/config/Config.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"
#include "core/math/IntegerFactorization.h"
#include "core/mpi/Environment.h"
#include "core/stringToNum.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/communication/PackInfo.h"
#include "field/iterators/FieldPointer.h"
#include "field/vtk/VTKWriter.h"

#include "gui/Gui.h"

#include "geometry/initializer/BoundaryFromDomainBorder.h"

#include "stencil/D3Q19.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <complex>
#include <string>


namespace walberla {

using flag_t = uint8_t;
using vec3_t = Vector3<real_t>;

typedef lbm::D3Q19< lbm::collision_model::SRT, true, lbm::force_model::None, 1 >  LM;

using CommunicationStencil = LM::CommunicationStencil;

using MyPdfField = lbm::PdfField<LM>;
typedef GhostLayerField< real_t, 1         >   ScalarField;
typedef GhostLayerField< Vector3<real_t>, 1>   VectorField;
using MyFlagField = FlagField<flag_t>;

typedef lbm::DefaultDiffusionBoundaryHandlingFactory< LM, MyFlagField >  MyBoundaryHandling;

const FlagUID& getFluidFlag(){ static FlagUID uid( "Fluid" ); return uid; }


class PlugFlow
{
public:
   using cplx_t = std::complex<real_t>;

   PlugFlow( real_t L, real_t H, real_t u, real_t k ) :
      period_( real_t(2)*math::pi/L ),
      lambda_( period_*sqrt( cplx_t(real_t(1), u/k/period_) ) ),
      emH_   ( real_t(1) - exp(-lambda_*H) ),
      epH_   ( real_t(1) - exp(+lambda_*H) ),
      eeH_   ( emH_ - epH_ ){}

   inline real_t operator()( real_t x, real_t y )
   {
      const cplx_t ly  = lambda_*y;
      const cplx_t emy = exp(-ly);
      const cplx_t epy = exp(+ly);
#if 1
      // exact solution
      return std::real( exp( cplx_t(real_t(0),period_*x) ) * ( emH_*epy - epH_*emy ) / eeH_ );
#else
      // integral solution
      return std::imag( exp( cplx_t(real_t(0),period_*x) ) * ( emH_*epy + epH_*emy ) / ( eeH_ * period_ * lambda_ ) );
#endif
   }

private:
   const real_t period_;
   const cplx_t lambda_;
   const cplx_t emH_;
   const cplx_t epH_;
   const cplx_t eeH_;
};


class TestSweep{
public:
   TestSweep( ConstBlockDataID pdfFieldID, real_t omega, real_t velocity, real_t error ) :
      pdfFieldID_(pdfFieldID), k_((real_t(1)/omega-real_c(0.5))/real_t(3)), velocity_(velocity), error_(error){}

   void operator()( IBlock* block );
private:
   ConstBlockDataID pdfFieldID_;

   const real_t k_;
   const real_t velocity_;
   const real_t error_;
};

void TestSweep::operator()( IBlock* block )
{
   auto srcPDF_ = block->getData< MyPdfField >( pdfFieldID_ );

   const uint_t sx = srcPDF_->xSize();
   const uint_t sy = srcPDF_->ySize();
   const uint_t sz = srcPDF_->zSize();

   real_t snumerical( real_t(0) ), sanalytical( real_t(0) );

   PlugFlow pf( real_c(sx), real_c(sy), velocity_, k_ );

   for(uint_t z=uint_t(0u); z<sz; ++z) {
      for(uint_t y=uint_t(0u); y<sy; ++y)
      {
         for(uint_t x=uint_t(0u); x<sx; ++x)
         {
#if 1
            // exact solution
            // integration is not needed, because we only use exponential functions
            const real_t analytical = pf( real_c(x)+real_c(0.5), real_c(y)+real_c(0.5) );
#else
            // integral solution
            const real_t analytical = pf( real_c(x), real_c(y) ) + pf( real_c(x+uint_t(1u)), real_c(y+uint_t(1u)) ) - pf( real_c(x+uint_t(1u)), real_c(y) ) - pf( real_c(x), real_c(y+uint_t(1u)) );
#endif
            const real_t numerical  = srcPDF_->getDensity( cell_idx_c(x), cell_idx_c(y), cell_idx_c(z) );

            snumerical  += real_c( pow(analytical-numerical,2) );
            sanalytical += analytical*analytical;
         }
      }
   }

   const real_t E = real_c( sqrt(snumerical/sanalytical) );
   if( E > error_ )
      WALBERLA_ABORT( "E^2 = " << E << " > " << error_ );

   WALBERLA_LOG_RESULT( "Total Error: " << E << " < " << error_);
}


class CosBoundaryConfiguration : public MyBoundaryHandling::DiffusionDirichlet_T::ScalarConfiguration
{
public:
   CosBoundaryConfiguration( real_t period ) : period_( period ){}

   void val( real_t& _val, cell_idx_t x, cell_idx_t, cell_idx_t ) const override { _val = real_c( cos( period_*( real_c(x) + real_c(0.5) ) ) ); }

private:
   const real_t period_ {};
   const vec3_t vel_ {};
};


int main( int argc, char **argv )
{
   real_t omega  = real_t(  1.2 );
   uint_t length = uint_t( 16u  );
   uint_t width  = uint_t(  1u  );
   real_t velx   = real_t(  0.1 );
   uint_t time   = uint_t( 50u  );
   real_t error  = real_t(  0.1 );

   bool useGui = false;
   bool useVTK = false;

   // --- read arguments --- //
   if( argc > 1 ) {
      std::vector<std::string> args( argv, argv + argc );
      for( uint_t i = 1; i < uint_c(argc); ++i ) {
              if( std::string(argv[i]) == "-o" ) omega  = stringToNum<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-l" ) length = stringToNum<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "-w" ) width  = stringToNum<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "-v" ) velx   = stringToNum<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-t" ) time   = stringToNum<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "-e" ) error  = stringToNum<real_t>( args[++i] );
         else if( std::string(argv[i]) == "--gui" ) useGui = true;
         else if( std::string(argv[i]) == "--vtk" ) useVTK = true;
         else if( argv[i][0] != '-' ){
            WALBERLA_ABORT( "Usage: --option or -option value" );
         } else
            WALBERLA_LOG_WARNING( "Ignore unknown option " << argv[i++] );
      }
   }
   vec3_t velocity( velx, real_t(0), real_t(0) );
      
   auto blockStorage = blockforest::createUniformBlockGrid(
      1,      1,      1,      // blocks/processes in x/y/z direction
      length, length, width,  // cells per block in x/y/z direction
      real_t(1),              // cell size
      true,                   // one block per process
      true,   false,  true,   // periodicity
      false );

   BlockDataID velFieldID   = field::addToStorage<VectorField>( blockStorage, "Velocity field", velocity );
   BlockDataID flagFieldID  = field::addFlagFieldToStorage<MyFlagField>( blockStorage, "Flag field" );

   LM lm = LM( lbm::collision_model::SRT( omega ) );
   BlockDataID srcFieldID   = lbm::addPdfFieldToStorage( blockStorage, "PDF AdDif field", lm, vec3_t(), real_t(0) );

   BlockDataID boundaryHandling = MyBoundaryHandling::addDefaultDiffusionBoundaryHandlingToStorage( blockStorage, "BoundaryHandling", flagFieldID, getFluidFlag(), srcFieldID );

   auto cbc = make_shared<CosBoundaryConfiguration>( real_t(2)*math::pi/real_c(length) );
   geometry::initializer::BoundaryFromDomainBorder<MyBoundaryHandling::BoundaryHandling_T> bfdb( *blockStorage, boundaryHandling );
   bfdb.init( MyBoundaryHandling::getDiffusionDirichletBoundaryUID(), stencil::N, cbc, -1, 1 );
   bfdb.init( MyBoundaryHandling::getDiffusionDirichletBoundaryUID(), stencil::S, cbc, -1, 1 );

   SweepTimeloop timeloop( blockStorage->getBlockStorage(), time );

   blockforest::communication::UniformBufferedScheme< CommunicationStencil > scheme( blockStorage );
   scheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LM > >( srcFieldID ) );
   timeloop.addFuncBeforeTimeStep( scheme, "Communication" );

   timeloop.add() << Sweep( MyBoundaryHandling::BoundaryHandling_T::getBlockSweep(boundaryHandling) );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseAdvectionDiffusionSweep< LM, VectorField, MyFlagField >(
                                                srcFieldID, velFieldID, flagFieldID, getFluidFlag() ) ), "LBM_SRT" );

   if( useGui )
   {
      field::addFieldAdaptor< lbm::Adaptor< LM >::Density                     >( blockStorage, srcFieldID, "E" );
      field::addFieldAdaptor< lbm::Adaptor< LM >::StreamMomentumDensityVector >( blockStorage, srcFieldID, "j" );

      GUI gui( timeloop, blockStorage, argc, argv );
      lbm::connectToGui<LM>( gui );
      gui.run();
   }
   else
   {
      timeloop.run();

      TestSweep ts( srcFieldID, omega, velx, error );
      for( auto block = blockStorage->begin(); block != blockStorage->end(); ++block )
         ts( block.get() );
   }
   
   if( useVTK )
   {
      auto vtkOut = vtk::createVTKOutput_BlockData( *blockStorage, "fields", uint_t(1u), uint_t(0u), false, "vtk_out/DiffusionDirichlet" );
      auto scalarWriter   = make_shared< lbm::DensityVTKWriter<LM>  >( srcFieldID,   "E"   );
      auto velocityWriter = make_shared< field::VTKWriter<VectorField> >( velFieldID, "u"   );
      vtkOut->addCellDataWriter( scalarWriter   );
      vtkOut->addCellDataWriter( velocityWriter );
      vtkOut->write();
   }

   return 0;
}

} // namespace walberla


int main(int argc, char **argv) {
   walberla::debug::enterTestMode();
   walberla::mpi::Environment env( argc, argv );
   return walberla::main( argc, argv );
}
