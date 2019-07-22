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
//! \file SimpleDiffusionDirichlet.cpp
//! \ingroup lbm
//! \author Matthias Markl <matthias.markl@fau.de>
//! \brief Calculates temporal evolution of a temperature profile between two plates
//!
//! Configuration:
//!  - CHANNEL_LENGTH  distance between both plates (x)
//!  - CHANNEL_WIDTH   y and z size
//!  - OMEGA           relaxation parameter -> depending on thermal diffusivity k: omega = 1/(3.k+0.5)
//!  - DELTA_SCALAR    temperature difference between plates
//!  - NUM_TIMESTEPS   amount of timesteps to simulate
//! Setup:
//!  - periodic in y and z coordinate
//!  - plates in the x-z planes
//!  - no hydrodynamic computations
//!  - comparison with analytical solution in each timestep
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "boundary/BoundaryHandling.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/SharedFunctor.h"
#include "core/cell/CellInterval.h"
#include "core/cell/CellSet.h"
#include "core/cell/CellVector.h"
#include "core/config/Config.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"
#include "core/math/IntegerFactorization.h"
#include "core/math/Limits.h"
#include "core/mpi/Environment.h"

#include "domain_decomposition/SharedSweep.h"

#include "lbm/boundary/factories/DefaultDiffusionBoundaryHandling.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/refinement/PdfFieldSyncPackInfo.h"
#include "lbm/refinement/TimeStep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "field/iterators/FieldPointer.h"

#include "stencil/D3Q19.h"

#include "timeloop/SweepTimeloop.h"

#include <stdexcept>
#include <array>
#include <functional>
#include <string>

#include "gather/GnuPlotGraphWriter.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/Velocity.h"
#include "vtk/VTKOutput.h"


namespace walberla {

typedef GhostLayerField< real_t, 1 >          ScalarField;
typedef GhostLayerField< Vector3<real_t>, 1 > VectorField;

typedef lbm::D3Q19< lbm::collision_model::SRT, true, lbm::force_model::None, 1 >  LM;

using CommunicationStencil = LM::CommunicationStencil;
using MyPdfField = lbm::PdfField<LM>;

using flag_t = uint8_t;
using MyFlagField = FlagField<flag_t>;

typedef lbm::DefaultDiffusionBoundaryHandlingFactory< LM, MyFlagField > MyBoundaryHandling;

/////////////////////
// BLOCK STRUCTURE //
/////////////////////

static void refinementSelection( SetupBlockForest& forest, const uint_t levels )
{
   const AABB & domain = forest.getDomain();

   const real_t domainxMax = domain.xMax() / real_c( pow( real_t(2), int_c( levels - uint_t(1u) ) ) );

   AABB left( domain.xMin(), domain.yMin(), domain.zMin(),
              domainxMax,    domain.yMax(), domain.zMax() );

   for( auto block = forest.begin(); block != forest.end(); ++block )
      if( block->getAABB().intersects( left ) )
         if( block->getLevel() < ( levels - uint_t(1u) ) )
            block->setMarker( true );

}

static void workloadAndMemoryAssignment( SetupBlockForest& forest ) {

   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      block->setWorkload( numeric_cast< workload_t >( uint_t(1) << block->getLevel() ) );
      block->setMemory( numeric_cast< memory_t >(1) );
   }
}


shared_ptr< StructuredBlockForest > makeStructuredBlockStorage( uint_t length, uint_t width, uint_t refinement )
{
   // initialize SetupBlockForest = determine domain decomposition
    SetupBlockForest sforest;

    uint_t cells[]  = { length, width, width  };
    uint_t blocks[] = { uint_t(1u), uint_t(1u), uint_t(1u) };
    sforest.addRefinementSelectionFunction( std::bind( refinementSelection, std::placeholders::_1, refinement ) );
    sforest.addWorkloadMemorySUIDAssignmentFunction( workloadAndMemoryAssignment );

    sforest.init(
       AABB( real_t(0),        real_t(0),        real_t(0),             // blocks/processes in x/y/z direction
             real_c(cells[0]), real_c(cells[1]), real_c(cells[2]) ),    // cells per block in x/y/z direction
             blocks[0]  , blocks[1]  , blocks[2],                       // one block per process
             false      , true       , true);                           // periodicity

    // calculate process distribution
    const memory_t memoryLimit = math::Limits< memory_t >::inf();

    sforest.balanceLoad( blockforest::StaticLevelwiseCurveBalance(true), uint_c( MPIManager::instance()->numProcesses() ), real_t(0), memoryLimit, true );

    MPIManager::instance()->useWorldComm();

    // create StructuredBlockForest (encapsulates a newly created BlockForest)
    shared_ptr< StructuredBlockForest > sbf =
          make_shared< StructuredBlockForest >( make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, false ),
                   cells[0]/blocks[0], cells[1]/blocks[1], cells[2]/blocks[2] );
    sbf->createCellBoundingBoxes();
    return sbf;
}

/////////////////////
// Boundary        //
/////////////////////

const FlagUID& getFluidFlag()           { static FlagUID uid( "Fluid" );           return uid; }

void setFlags( shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & boundaryHandlingId, uint_t ghostLayers, bool closed, real_t dv )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      auto * boundaryHandling = block->getData< MyBoundaryHandling::BoundaryHandling_T >( boundaryHandlingId );

      const uint_t level = blocks->getLevel(*block);

      CellInterval domainBB = blocks->getDomainCellBB( level );
      blocks->transformGlobalToBlockLocalCellInterval( domainBB, *block );

      //const cell_idx_t width = cell_idx_c( uint_t(1) << level ) - cell_idx_t(1);

      domainBB.xMin() -= cell_idx_c( 1 );
      domainBB.xMax() += cell_idx_c( 1 );
      domainBB.yMin() -= cell_idx_c( ghostLayers );
      domainBB.yMax() += cell_idx_c( ghostLayers );
      domainBB.zMin() -= cell_idx_c( ghostLayers );
      domainBB.zMax() += cell_idx_c( ghostLayers );

      MyBoundaryHandling::SimpleDiffusionDirichlet_T::ScalarConfiguration scl( real_t(1)+dv );
      MyBoundaryHandling::SimpleDiffusionDirichlet_T::ScalarConfiguration scr( real_t(1)    );

      // LEFT
      CellInterval left( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(),
                         domainBB.xMin(), domainBB.yMax(), domainBB.zMax());
      boundaryHandling->forceBoundary( MyBoundaryHandling::getSimpleDiffusionDirichletFlagUID1(), left, scl );

      // RIGHT
      CellInterval right( domainBB.xMax(), domainBB.yMin(), domainBB.zMin(),
                          domainBB.xMax(), domainBB.yMax(), domainBB.zMax());
      boundaryHandling->forceBoundary( MyBoundaryHandling::getSimpleDiffusionDirichletFlagUID2(), right, scr );


      if( closed)
      {
         domainBB.xMin() -= cell_idx_c( ghostLayers );
         domainBB.xMax() += cell_idx_c( ghostLayers );
         domainBB.yMin() -= cell_idx_c( ghostLayers );
         domainBB.yMax() += cell_idx_c( ghostLayers );
         domainBB.zMin() -= cell_idx_c( 1 );
         domainBB.zMax() += cell_idx_c( 1 );

         // TOP
         CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
         boundaryHandling->forceBoundary(MyBoundaryHandling::getFreeDiffusionFlagUID(), top );

         // BOTTOM
         CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
         boundaryHandling->forceBoundary( MyBoundaryHandling::getFreeDiffusionFlagUID(), bottom );
      }

       boundaryHandling->fillWithDomain( ghostLayers );
   }
}


class TestSweep
{
public:
   TestSweep( ConstBlockDataID pdfFieldID, real_t omega, real_t minValue, real_t maxValue, uint_t length, uint_t time, shared_ptr<StructuredBlockForest> & blocks ) :
         pdfFieldID_(pdfFieldID),
         blocks_ (blocks),
         k_((real_t(1)/omega-real_c(0.5))/real_t(3)),
         minValue_(minValue),
         maxValue_(maxValue),
         delta_( maxValue - minValue),
         length_(real_c(length)),
         lengthInv_(real_t(1)/real_c(length)),
         pi_(math::pi),
         piInv_(real_t(1)/math::pi),
         valid_(uint_c(std::ceil(omega*omega*omega*real_t(10)))),
         time_( time ),
         expArray(),
         timestep_( uint_t(0u) ),
         E_max_(maxValue-minValue),
         E_mean_(maxValue-minValue)
{
#ifdef TEST_USES_VTK_OUTPUT
        error_.resize(time_);
#endif
        WALBERLA_UNUSED(maxValue_);
}

   void operator()();
private:
   ConstBlockDataID pdfFieldID_;
   shared_ptr<StructuredBlockForest> & blocks_;
   const real_t k_;
   const real_t minValue_;
   const real_t maxValue_;
   const real_t delta_;
   const real_t length_;
   const real_t lengthInv_;
   const real_t pi_;
   const real_t piInv_;

   const uint_t valid_;
   const uint_t time_;

   std::array<real_t,1000> expArray;

   uint_t timestep_;
   real_t E_max_;
   real_t E_mean_;

#ifdef TEST_USES_VTK_OUTPUT
   std::vector< std::vector<double> > error_;
#endif


   inline real_t diffusion( const real_t x )
   {
      const real_t xiL = x*lengthInv_;

      real_t y = real_t(0);
      real_t f = real_t(1);
      for ( uint_t n = 1; n<uint_t(1000u); ++n ){
         const real_t npi = real_c(n)*pi_;
         f *= -real_t(1);
         y += f/real_c(n) * real_c( sin(npi*xiL) ) * expArray[n-1];
      }
      return delta_*(real_t(2)*y*piInv_ + xiL) + minValue_;
   }


   inline real_t diffusionIntegral( const real_t x )
   {
      const real_t xiL  = x*lengthInv_;

      real_t y = real_t(0);
      real_t f = real_t(1);
      for ( uint_t n = 1; n<uint_t(1000u); ++n ){
         const real_t npi = real_c(n)*pi_;
         f *= -real_t(1);
         y -= f*length_/(real_c(n)*npi) * real_c(cos(npi*xiL)) * expArray[n-uint_t(1u)];
      }
      return delta_*(real_t(2)*y*piInv_ + real_c(0.5)*x*xiL) + x*minValue_;
   }
};

#define ANALYTIC_INTEGRAL

void TestSweep::operator()()
{
   ++timestep_;

   E_mean_ = real_t(0);
   real_t E_max  = real_t(0);
   uint_t blockcount = uint_t(0u);

   //loop block
   for( auto block = blocks_->begin(); block != blocks_->end(); ++block, ++blockcount )
   {
      real_t E_mean = real_t(0);

      Cell localCell, globalCell;
      const uint_t level = blocks_->getLevel(*block);
      const real_t ktiLL = -k_*real_c(timestep_)*lengthInv_*lengthInv_;
      const real_t pow2level = real_c( pow( real_t(2), int_c(level) ) );

      for ( uint_t n = 1; n<1000; ++n ){
           const real_t npi = real_c(n)*pi_;
           expArray[n-1] = real_c( exp(npi*npi*ktiLL) );
      }

      auto srcPDF_ = block->getData< MyPdfField >( pdfFieldID_ );

      const uint_t sx = srcPDF_->xSize();
      const uint_t sy = srcPDF_->ySize();
      const uint_t sz = srcPDF_->zSize();

      real_t numerical, lastNumerical, analytical;// E_max = 0.0;

#ifdef ANALYTIC_INTEGRAL
      localCell = Cell(0,0,0);
      blocks_->transformBlockLocalToGlobalCell( globalCell, *block, localCell );
      real_t lastAnalytical = diffusionIntegral( length_ - real_c(globalCell[0]) / pow2level );
#endif

      for( uint_t x = 0; x < sx; ++x )
      {
         localCell = Cell( cell_idx_c( x ), 0, 0 );
         blocks_->transformBlockLocalToGlobalCell( globalCell, *block, localCell );

#ifdef ANALYTIC_INTEGRAL
         const real_t nextAnalytical = diffusionIntegral( length_ - real_c( globalCell[0] + 1 ) / pow2level );
         analytical = ( lastAnalytical - nextAnalytical ) * pow2level;
         lastAnalytical = nextAnalytical;
#else
         // the single diffusion value does not represent the overall temperature in a cell
         analytical = diffusion( length_ - ( real_c( globalCell[0] ) + 0.5 ) / pow2level );
#endif
         lastNumerical = real_t( 0 );
         for( uint_t z = 0; z<sz; ++z ) {
            for( uint_t y = 0; y<sy; ++y )
            {
               numerical = srcPDF_->getDensity( cell_idx_c( x ), cell_idx_c( y ), cell_idx_c( z ) );

               if( lastNumerical > real_t( 0 ) && std::fabs( numerical - lastNumerical ) > 1e-14 )
                  WALBERLA_ABORT( "Unequal numerical values at " << Cell( x, y, z ) << ", t=" << timestep_ << ": E_num=" << numerical << " != " << lastNumerical );

               const real_t E = real_c( sqrt( ( analytical - numerical )*( analytical - numerical ) ) );

               //if( E > E_max_ && E/E_max_ > 1.05 && E > 0.00001 && timestep_ > valid_ )
               //   WALBERLA_ABORT( "Increasing error at " << Cell(x,y,z) << ", t=" << timestep_ << ": E_num=" << numerical << " != E_ana=" << analytical << " -> err=" << E << " fac=" << E/E_max_ );

               E_max = std::max( E_max, E );
               E_mean += E;
               lastNumerical = numerical;
            }
         }

         //WALBERLA_LOG_RESULT( "Last Analytical: " << analytical << " and Last Numerical: " << lastNumerical );
      }
      //E_max  /=                    pow( real_t(8),real_c(level) );
      E_mean /= real_c( sx*sy*sz ) * real_c( pow( real_t( 8 ), int_c( level ) ) );
      E_mean_ += E_mean;
   }
   E_max_   = E_max;
   E_mean_ /= real_c(blockcount);

#ifdef TEST_USES_VTK_OUTPUT
   error_[timestep_].push_back( real_c(timestep_ ) );
   error_[timestep_].push_back( E_mean_ );
#endif

   //WALBERLA_LOG_RESULT( "Max  Error: " << E_max << " and Mean Error: " << E_mean_ );

   if( timestep_ == time_ ){
#ifdef TEST_USES_VTK_OUTPUT
      const std::string & file_name   = "meanError_r1";
      const std::string & file_ending = "dat";
      gather::GnuPlotGraphWriter gnuWriter(file_name, file_ending);
      gnuWriter.process(error_);
#endif
      WALBERLA_LOG_RESULT( "Max  Error: " << E_max_  << "\n"<< "Mean Error: " << E_mean_ );
   }
}

} // namespace walberla

int main( int argc, char **argv )
{
   using namespace walberla;
   debug::enterTestMode();
   mpi::Environment env( argc, argv );
   logging::Logging::printHeaderOnStream();

   uint_t length = uint_t( 16u  );
   uint_t width  = uint_t( 16u  );
   uint_t time   = uint_t( 10u  );
   real_t dv     = real_t(  1   );
   real_t omega  = real_t(  0.7 );
   bool   closed = false;
   uint_t levels = uint_t(  1u );

   bool useVTK = false;

   if( argc > 1 ) {
      std::vector<std::string> args( argv, argv + argc );
      for( uint_t i = 1; i < uint_c(argc); ++i ) {
              if( std::string(argv[i]) == "-l"    )   length  = string_to_num<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "-w"    )   width   = string_to_num<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "-t"    )   time    = string_to_num<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "-d"    )   dv      = string_to_num<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-o"    )   omega   = string_to_num<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-c"    )   closed  = string_to_num<int>( args[++i] ) != 0;
         else if( std::string(argv[i]) == "-r"    )   levels += string_to_num<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "--vtk" )   useVTK  = true;
         else if( argv[i][0] != '-' ){
            std::cerr << "Usage: -option value" << std::endl; return EXIT_FAILURE;
         } else
            std::cerr << "Ignore unknown option " << argv[i++] << std::endl;
      }
   }

   uint_t ghostLayers = uint_t(4u);
   if(closed || levels == uint_t(1u))
   {
      levels      = uint_t(1u);
      ghostLayers = uint_t(1u);
   }

   auto blockStorage = makeStructuredBlockStorage( length, width, levels );

   LM          lm             = LM( lbm::collision_model::SRT( omega ) );
   BlockDataID advDiffFieldID = lbm::addPdfFieldToStorage( blockStorage, "PDF field", lm, Vector3<real_t>(), real_t(1), ghostLayers );

   BlockDataID flagFieldID    = field::addFlagFieldToStorage<MyFlagField>( blockStorage, "Flag field", ghostLayers );
   BlockDataID velFieldID     = field::addToStorage<VectorField>( blockStorage, "Velocity field", Vector3<real_t>() );

   BlockDataID boundaryHandling = MyBoundaryHandling::addDefaultDiffusionBoundaryHandlingToStorage( blockStorage, "DiffusionboundaryHandling", flagFieldID, getFluidFlag(), advDiffFieldID );

   setFlags( blockStorage, boundaryHandling, ghostLayers, closed, dv );

   SweepTimeloop timeloop( blockStorage->getBlockStorage(), time );
   timeloop.addFuncAfterTimeStep( TestSweep(advDiffFieldID, omega, real_t(1), real_t(1)+dv, length, time, blockStorage ),"check error! ");

   if( levels == uint_t(1u) )
   {
      blockforest::communication::UniformBufferedScheme< CommunicationStencil > scheme( blockStorage );
      scheme.addPackInfo( make_shared< field::communication::PackInfo< MyPdfField > >( advDiffFieldID ) );
      timeloop.addFuncBeforeTimeStep( scheme, "Communication" );

      timeloop.add() << Sweep( MyBoundaryHandling::BoundaryHandling_T::getBlockSweep(boundaryHandling), "Boundary" );
      timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseAdvectionDiffusionSweep< LM, VectorField, MyFlagField >(
                                                   advDiffFieldID, velFieldID, flagFieldID, getFluidFlag() ) ), "LBM_SRT" );
   }
   else
   {
      auto mySweep = lbm::makeCellwiseAdvectionDiffusionSweep< LM, VectorField, MyFlagField >( advDiffFieldID, velFieldID, flagFieldID, getFluidFlag() );

      timeloop.addFuncBeforeTimeStep( makeSharedFunctor( lbm::refinement::makeTimeStep< LM, MyBoundaryHandling::BoundaryHandling_T >(
                                                            blockStorage, mySweep, advDiffFieldID, boundaryHandling  ) ), "LBM refinement time step" );
   }
   
   WcTimingPool timeloopTiming;
   timeloop.run(timeloopTiming );
   timeloopTiming.logResultOnRoot();

   if( useVTK )
   {
      auto vtkOut = vtk::createVTKOutput_BlockData( *blockStorage, "block_data", uint_t(1u), uint_t(0u), false, "vtk_out/SimpleDiffusionTest" );
      auto densityWriter  = make_shared<lbm::DensityVTKWriter<LM>>( advDiffFieldID, "E" );
      vtkOut->addCellDataWriter( densityWriter );
      vtkOut->write();
      vtk::writeDomainDecomposition( blockStorage, "domain_decomposition", "vtk_out/SimpleDiffusionTest" );
      field::createVTKOutput<MyFlagField>( flagFieldID, *blockStorage, "flag_field", uint_t(1u), uint_t(1u), false, "vtk_out/SimpleDiffusionTest" )();
   }

   logging::Logging::printFooterOnStream();

   return 0;
}

