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
//! \file DiffusionTest.cpp
//! \ingroup lbm
//! \author Matthias Markl <matthias.markl@fau.de>
//! \brief Calculates temporal evolution of a temperature profile driven by a cosinus velocity profile over time
//!
//! Configuration:
//!  -d   physical diffusivity
//!  -dx  space discretization
//!  -dt  time discretization
//!  -dv  delta diffusive value
//!  -v   maximum velocity in dim direction
//!  -dim orientation of block setup, x = 0, y = 1, z = 2
//!  -err maximum mean error
//!  -c   corrected model = 1, normal model = 0
//! Setup:
//!  - periodic in x, y and z coordinate
//!  - no hydrodynamic computations
//!  - comparison with analytical solution in each timestep
//
// @see Chopard et.al.: The lattice Boltzmann advection-diffusion model revisited
//
//======================================================================================================================

#include "lbm/boundary/SimpleDiffusionDirichlet.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/gui/Connection.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/vtk/Density.h"

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/OpenMP.h"
#include "core/waLBerlaBuildInfo.h"
#include "core/cell/CellInterval.h"
#include "core/cell/CellSet.h"
#include "core/cell/CellVector.h"
#include "core/config/Config.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"
#include "core/math/IntegerFactorization.h"
#include "core/mpi/Environment.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/communication/PackInfo.h"
#include "field/iterators/FieldPointer.h"

#include "gui/Gui.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <functional>
#include <string>


namespace walberla {

using flag_t = uint8_t;
using vec3_t = Vector3<real_t>;

typedef GhostLayerField< real_t, 1         >  ScalarField;
typedef GhostLayerField< vec3_t, 1         >  VectorField;
using MyFlagField = FlagField<flag_t>;

typedef lbm::D3Q19< lbm::collision_model::SRT, true, lbm::force_model::Correction<VectorField>, 1 > AdvDiffLatticeModel_Corr;
typedef lbm::D3Q19< lbm::collision_model::SRT, true, lbm::force_model::None,                    1 > AdvDiffLatticeModel_None;

const FlagUID& getFluidFlag(){ static FlagUID uid( "Fluid" ); return uid; }


//********************************************************************************************************************************
// Create Adv Diff PDF Field
template< typename AdvDiffLatticeModel >
BlockDataID createAdvDiffPdfField( shared_ptr<StructuredBlockStorage> blockStorage, const BlockDataID & oldMomDensity, real_t omega );

template< >
BlockDataID createAdvDiffPdfField<AdvDiffLatticeModel_Corr>( shared_ptr<StructuredBlockStorage> blockStorage, const BlockDataID & oldMomDensity, real_t omega ){
   AdvDiffLatticeModel_Corr advDiffLatticeModel = AdvDiffLatticeModel_Corr( lbm::collision_model::SRT( omega ), lbm::force_model::Correction<VectorField>( oldMomDensity ) );
   return lbm::addPdfFieldToStorage( blockStorage, "PDF field", advDiffLatticeModel, vec3_t(), real_t(0) );
}

template< >
BlockDataID createAdvDiffPdfField<AdvDiffLatticeModel_None>( shared_ptr<StructuredBlockStorage> blockStorage, const BlockDataID &, real_t omega ){
   AdvDiffLatticeModel_None advDiffLatticeModel = AdvDiffLatticeModel_None( lbm::collision_model::SRT( omega ) );
   return lbm::addPdfFieldToStorage( blockStorage, "PDF field", advDiffLatticeModel, vec3_t(), real_t(0) );
}
//********************************************************************************************************************************


//********************************************************************************************************************************
// Create Flag Field
MyFlagField* flagFieldCreationFunction( IBlock * const block, StructuredBlockStorage * const storage )
{
   auto field = new MyFlagField(
            storage->getNumberOfXCells( *block ), // number of cells in x direction
            storage->getNumberOfYCells( *block ), // number of cells in y direction
            storage->getNumberOfZCells( *block ), // number of cells in z direction
            1);     // ghost layers

   auto flag = field->registerFlag( getFluidFlag() );
   field->set( flag );

   return field;
}
//********************************************************************************************************************************


//********************************************************************************************************************************
// Helper Functions
void incTimeFunc( uint_t& timestep ){ ++timestep; }

void hydroFunc( IBlock* block, BlockDataID velFieldID, const vec3_t& u, const real_t tperiod, const uint_t& timestep )
{
   block->getData<VectorField>( velFieldID )->set( u*cos( real_c(timestep) * tperiod ) ); 
}

void prepFunc( const real_t u, const real_t dv, const real_t D, const real_t cperiod, const real_t tperiod, const uint_t& timestep, real_t& cosi, real_t& sisi, real_t& sexp )
{
   const real_t argument = u * cperiod / tperiod * real_c( sin( tperiod*real_c(timestep) ) );
   cosi = real_c( cos( argument ) );
   sisi = real_c( sin( argument ) );
   sexp = dv * real_c( exp( - cperiod*cperiod*D*real_c(timestep) ) );
}

template< typename AdvDiffPDFField >
void testFunc( IBlock* block, BlockDataID srcFieldID, const uint_t dim, const real_t v, const real_t cperiod, const real_t& cosi, const real_t& sisi, const real_t& sexp, real_t& E_mean_ )
{
   auto srcField = block->getData<AdvDiffPDFField>( srcFieldID );

   for( auto it = srcField->beginXYZ(); it != srcField->end(); ++it )
   {
      const real_t arg = cperiod * ( real_c( it.cell()[dim] ) + real_c(0.5) );
      const real_t ana = v + sexp * ( real_c( cos(arg) )*cosi + real_c( sin(arg) )*sisi );
      const real_t num = lbm::getDensity( srcField->latticeModel(), it );
      E_mean_ += (ana-num) * (ana-num);
   }
}
//********************************************************************************************************************************

template< typename AdvDiffLatticeModel >
int run( int argc, char **argv )
{
   // typedefs
   using AdvDiffStencil = typename AdvDiffLatticeModel::Stencil;
   using AdvDiffPDFField = lbm::PdfField<AdvDiffLatticeModel>;

#ifdef _OPENMP
   omp_set_num_threads( std::max( omp_get_num_threads()>>1, 4 ) );
#endif

   // --- physical default values --- //
   const real_t size[] = { real_t(1), real_t(1), real_t(1) };
   const real_t v      = real_t(1);

   real_t time = real_t( 0.2   );
   real_t err  = real_t( 0     );
   real_t d    = real_t( 1e-6  );
   real_t dx   = real_t( 0.01  );
   real_t dt   = real_t( 0.001 );
   real_t dv   = real_t( 1.0   );
   real_t u_in = real_t( 1     );
   uint_t dim  = uint_t( 2u    );

   bool useGui = false;
   bool useVTK = false;
   bool quiet  = false;

   // --- read arguments --- //
   if( argc > 1 ) {
      std::vector<std::string> args( argv, argv + argc );
      for( uint_t i = 1; i < uint_c(argc); ++i ) {
              if( std::string(argv[i]) == "-d"      )   d      = string_to_num<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-dim"    )   dim    = string_to_num<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "-dx"     )   dx     = string_to_num<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-dt"     )   dt     = string_to_num<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-dv"     )   dv     = string_to_num<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-v"      )   u_in   = string_to_num<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-t"      )   time   = string_to_num<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-err"    )   err    = string_to_num<real_t>( args[++i] );
         else if( std::string(argv[i]) == "--gui"   )   useGui = true;
         else if( std::string(argv[i]) == "--quiet" )   quiet  = true;
         else if( std::string(argv[i]) == "--vtk"   )   useVTK = true;
         else if( std::string(argv[i]) == "-c"      )   ++i;
         else if( argv[i][0] != '-' ){
            WALBERLA_ABORT( "Usage: -option value" );
         } else
            WALBERLA_LOG_WARNING( "Ignore unknown option " << argv[i++] );
      }
   }

   // --- make dimensionless --- //
   const uint_t timesteps = uint_t( time / dt );

   uint_t cells[]   = { uint_t(1u), uint_t(1u), uint_t(1u) };
   cells[dim] = uint_c( size[dim] / dx );

   if(!quiet) WALBERLA_LOG_RESULT( "Simulate " << cells[dim] << " cells in " << timesteps << " timesteps " );

   const real_t D     = d * dt / dx / dx;
   const real_t tau   = real_t(3) * D + real_c(0.5);
   const real_t omega = real_t(1) / tau;

   vec3_t u;
   u[dim] = u_in;
   u = u / dx * dt;

   if(!quiet) WALBERLA_LOG_RESULT( " -> u   = " << u   );
   if(!quiet) WALBERLA_LOG_RESULT( " -> tau = " << tau );

   const real_t tperiod = real_t(2) * math::M_PI / real_c( timesteps  );
   const real_t cperiod = real_t(2) * math::M_PI / real_c( cells[dim] );

   // --- create blockstorage --- //
   auto blockStorage = blockforest::createUniformBlockGrid(
      1,        1,        1,             // blocks/processes in x/y/z direction
      cells[0], cells[1], cells[2],      // cells per block in x/y/z direction
      dx,                                // cell size
      true,                              // one block per process
      true,     true,    true,           // periodicity
      false );

   // --- create fields --- //
   BlockDataID oldMomDensity = field::addToStorage<VectorField>( blockStorage, "Old Momentum Density", vec3_t(), field::fzyx, uint_t(0u) );

   BlockDataID srcFieldID = createAdvDiffPdfField<AdvDiffLatticeModel>( blockStorage, oldMomDensity, omega );

   BlockDataID velFieldID  = field::addToStorage<VectorField>( blockStorage, "Velocity field", Vector3<real_t>() );
   BlockDataID flagFieldID = blockStorage->template addStructuredBlockData<MyFlagField>( &flagFieldCreationFunction, "Flag field" );

   // --- additional sweep variables --- //
   real_t E_mean_ ( real_t(0)  );
   uint_t timestep( uint_t(0u) );
   real_t cosi(real_t(0)), sisi(real_t(0)), sexp(real_t(0));

   // --- data init --- //
   for( auto block = blockStorage->begin(); block != blockStorage->end(); ++block ){
      auto srcField = block->template getData< AdvDiffPDFField >( srcFieldID    );
      auto oMDField = block->template getData< VectorField     >( oldMomDensity );
      for( auto it = srcField->beginXYZ(); it != srcField->end(); ++it ){
         lbm::setToEquilibrium<AdvDiffLatticeModel>( it, u, v + dv * real_c( cos( ( real_c( it.cell()[dim] ) + real_c(0.5) ) * cperiod ) ) );
         oMDField->get(it.x(),it.y(),it.z()) = u * lbm::getDensity( srcField->latticeModel(), it );
      }
   }

   blockforest::communication::UniformBufferedScheme< AdvDiffStencil > flagScheme( blockStorage );
   flagScheme.addPackInfo( make_shared< field::communication::PackInfo< MyFlagField > >( flagFieldID ) );
   flagScheme.communicate();

   // --- timeloop setup --- //
   SweepTimeloop timeloop( blockStorage->getBlockStorage(), timesteps );

   blockforest::communication::UniformBufferedScheme< AdvDiffStencil > scheme( blockStorage );
   scheme.addPackInfo( make_shared< field::communication::PackInfo<  AdvDiffPDFField > >( srcFieldID ) );
   timeloop.addFuncBeforeTimeStep( scheme, "Communication" );

   using std::ref;

   timeloop.add() << Sweep( std::bind( hydroFunc, std::placeholders::_1, velFieldID, u, tperiod, ref(timestep) ), "Hydro Func" );
   
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseAdvectionDiffusionSweep< AdvDiffLatticeModel, VectorField, MyFlagField >(
                                                srcFieldID, velFieldID, flagFieldID, getFluidFlag() ) ), "LBM_SRT" );
  
   timeloop.add() << BeforeFunction( std::bind( prepFunc, u[dim], dv, D, cperiod, tperiod, ref(timestep), ref(cosi), ref(sisi), ref(sexp) ), "prepare test" )
                  << Sweep         ( std::bind( testFunc<AdvDiffPDFField>, std::placeholders::_1, srcFieldID, dim, v, cperiod, ref(cosi), ref(sisi), ref(sexp), ref(E_mean_) ), "Test Func" ) 
                  << AfterFunction ( std::bind( incTimeFunc, ref(timestep) ), "increment time" );

   // --- run timeloop --- //

   if ( useGui )
   {
      field::addFieldAdaptor< typename lbm::Adaptor<AdvDiffLatticeModel>::Density>                    ( blockStorage, srcFieldID, "E" );
      field::addFieldAdaptor< typename lbm::Adaptor<AdvDiffLatticeModel>::StreamMomentumDensityVector>( blockStorage, srcFieldID, "j" );

      GUI gui ( timeloop, blockStorage, argc, argv );
      lbm::connectToGui<AdvDiffLatticeModel>( gui );

      gui.run();
   }
   else if( !quiet )
   {
      WcTimingPool timeloopTiming;
      timeloop.run( timeloopTiming );
      timeloopTiming.logResultOnRoot();
   }
   else
   {
      timeloop.run( );
   }

   if( useVTK )
   {
      auto vtkOut = vtk::createVTKOutput_BlockData( *blockStorage, "fields", 1u, 0u, false, "vtk_out/DiffusionTest" );
      auto densityWriter  = make_shared<lbm::DensityVTKWriter<AdvDiffLatticeModel>>( srcFieldID, "E" );
      vtkOut->addCellDataWriter( densityWriter );
      vtkOut->write();
   }

   E_mean_ /= real_c(timesteps) * real_c(cells[dim]);
   if( err > real_t(0) && E_mean_ > err ){
      WALBERLA_ABORT( "E_mean = " << E_mean_ << " > " << err );
   } else {
      if(!quiet) WALBERLA_LOG_RESULT( "E^2 = " << E_mean_ << " E^1 = " << sqrt(E_mean_) );
   }

   return 0;
}

} // namespace walberla



int main(int argc, char **argv)
{
   using namespace walberla;

   debug::enterTestMode();
   mpi::Environment env( argc, argv );
   bool corr = true;
   for( int i=0; i<argc; ++i ){
      if( std::string(argv[i]) == "-c" ){
         corr = std::atoi( argv[++i] ) != 0;
         break;
      }
   }

   if( corr )
      return run< AdvDiffLatticeModel_Corr >( argc, argv );
   else
      return run< AdvDiffLatticeModel_None >( argc, argv );
}
