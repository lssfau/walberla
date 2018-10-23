#include "core/Environment.h"
#include "python_coupling/CreateConfig.h"
#include "blockforest/Initialization.h"
#include "lbm/field/PdfField.h"
#include "lbm/field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/AddToStorage.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/vtk/VTKOutput.h"
#include "lbm/PerformanceLogger.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "timeloop/all.h"
#include "core/math/Random.h"
#include "geometry/all.h"
#include "cuda/HostFieldAllocator.h"
#include "cuda/communication/GPUPackInfo.h"
#include "core/timing/TimingPool.h"
#include "core/timing/RemainingTimeLogger.h"
#include "cuda/AddGPUFieldToStorage.h"
#include "cuda/communication/UniformGPUScheme.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "domain_decomposition/SharedSweep.h"

#include "EquivalenceTest_LatticeModel.h"
#include "EquivalenceTest_GPUKernel.h"
#include "EquivalenceTest_GPUPackInfo.h"

using namespace walberla;

using NativeLatticeModel_T = lbm::D3Q19<lbm::collision_model::SRT, false>;
using GeneratedLatticeModel_T = lbm::EquivalenceTest_LatticeModel;

using Stencil_T = GeneratedLatticeModel_T::Stencil;
using CommunicationStencil_T = GeneratedLatticeModel_T::CommunicationStencil;
using NativePdfField_T = lbm::PdfField<NativeLatticeModel_T>;
using GeneratedPdfField_T = lbm::PdfField<GeneratedLatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

using CpuCommScheme_T = blockforest::communication::UniformBufferedScheme<CommunicationStencil_T>;
using GpuCommScheme_T = cuda::communication::UniformGPUScheme<CommunicationStencil_T>;


template<typename PdfField_T>
void initPdfField( const shared_ptr<StructuredBlockForest> &blocks, BlockDataID pdfFieldId )
{
   auto domainBB = blocks->getDomainCellBB();

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto pdfField = blockIt->getData<PdfField_T>( pdfFieldId );
      Cell offset( 0, 0, 0 );
      blocks->transformBlockLocalToGlobalCell( offset, *blockIt );

      WALBERLA_FOR_ALL_CELLS_XYZ( pdfField,
         auto globalX = real_c( offset[0] + x );
         auto globalZ = real_c( offset[2] + z );
         auto xArg = real_c(std::sin(real_c(globalX) / real_t(4) * real_c(domainBB.size(0)) ));
         auto zArg = real_c(std::sin(real_c(globalZ) / real_t(4) * real_c(domainBB.size(2)) ));
         pdfField->setToEquilibrium( x, y, z, Vector3<real_t>( 0.05 * std::sin(xArg), 0, 0.05 * std::cos(zArg)));
      );
   }
}


int main( int argc, char **argv )
{
   mpi::Environment env( argc, argv );

   for( auto cfg = python_coupling::configBegin( argc, argv ); cfg != python_coupling::configEnd(); ++cfg )
   {
      auto config = *cfg;
      auto parameters = config->getOneBlock( "Parameters" );

      auto blocks = blockforest::createUniformBlockGridFromConfig( config );

      const real_t omega = parameters.getParameter<real_t>( "omega", real_c( 1.4 ));
      const uint_t timesteps = parameters.getParameter<uint_t>( "timesteps", uint_c( 50 ));

      // Boundary
      BlockDataID flagFieldId = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );
      const FlagUID fluidFlagUID( "Fluid" );
      geometry::setNonBoundaryCellsToDomain<FlagField_T>( *blocks, flagFieldId, fluidFlagUID );


      // Part 1 : Native walberla
      NativeLatticeModel_T nativeLatticeModel = NativeLatticeModel_T( lbm::collision_model::SRT( omega ));
      BlockDataID pdfFieldNativeId = lbm::addPdfFieldToStorage( blocks, "pdfNative", nativeLatticeModel, field::fzyx );
      initPdfField<NativePdfField_T >( blocks, pdfFieldNativeId );
      CpuCommScheme_T nativeComm( blocks );
      nativeComm.addPackInfo( make_shared< lbm::PdfFieldPackInfo< NativeLatticeModel_T > >( pdfFieldNativeId ) );
      auto nativeSweep = lbm::makeCellwiseSweep< NativeLatticeModel_T , FlagField_T >( pdfFieldNativeId, flagFieldId, fluidFlagUID );

      SweepTimeloop nativeTimeLoop( blocks->getBlockStorage(), timesteps );
      nativeTimeLoop.add() << BeforeFunction( nativeComm, "communication" )
                           << Sweep(makeSharedSweep(nativeSweep), "native stream collide" );
      nativeTimeLoop.run();


      // Part 2: Generated CPU Version
      GeneratedLatticeModel_T generatedLatticeModel = GeneratedLatticeModel_T( omega );
      BlockDataID pdfFieldGeneratedId = lbm::addPdfFieldToStorage( blocks, "pdfGenerated", generatedLatticeModel, field::fzyx );
      initPdfField<GeneratedPdfField_T >( blocks, pdfFieldGeneratedId );
      CpuCommScheme_T cpuComm( blocks );
      cpuComm.addPackInfo( make_shared< lbm::PdfFieldPackInfo< GeneratedLatticeModel_T > >( pdfFieldGeneratedId ) );
      SweepTimeloop cpuTimeLoop( blocks->getBlockStorage(), timesteps );
      cpuTimeLoop.add() << BeforeFunction( cpuComm, "communication" )
                        << Sweep(GeneratedLatticeModel_T::Sweep( pdfFieldGeneratedId ), "generated stream collide on cpu" );
      cpuTimeLoop.run();


      // Part 3: Generated GPU Version
      bool overlapCommunication = parameters.getParameter<bool>( "overlapCommunication", true );
      bool cudaEnabledMPI = parameters.getParameter<bool>( "cudaEnabledMPI", false );

      BlockDataID pdfShadowCPU = lbm::addPdfFieldToStorage( blocks, "cpu shadow field", generatedLatticeModel, field::fzyx );
      initPdfField<GeneratedPdfField_T >( blocks, pdfShadowCPU );

      BlockDataID pdfGpuFieldId = cuda::addGPUFieldToStorage<GeneratedPdfField_T >( blocks, pdfShadowCPU, "pdfs on gpu", true );
      auto defaultKernelStream = overlapCommunication ? cuda::StreamRAII::newStream() : cuda::StreamRAII::defaultStream();
      auto innerKernelStartedEvent = make_shared<cuda::EventRAII>();

      pystencils::EquivalenceTest_GPUKernel cudaLbKernel( pdfGpuFieldId, omega, defaultKernelStream );
      GpuCommScheme_T gpuComm( blocks, innerKernelStartedEvent, cudaEnabledMPI );
      gpuComm.addPackInfo( make_shared<pystencils::EquivalenceTest_GPUPackInfo>( pdfGpuFieldId ));
      auto runCommunication = [&]() { gpuComm(); };

      SweepTimeloop gpuTimeLoop( blocks->getBlockStorage(), timesteps );
      if( !overlapCommunication )
      {
         gpuTimeLoop.add() << BeforeFunction( runCommunication, "gpu communication" )
                           << Sweep( cudaLbKernel, "LB stream & collide gpu" );
      }
      else
      {
         gpuTimeLoop.add() << Sweep( [&]( IBlock *b )
                                  {
                                     cudaEventRecord( *innerKernelStartedEvent, defaultKernelStream );
                                     cudaLbKernel.inner( b );
                                  }, "LBM @ inner" );
         gpuTimeLoop.add() << BeforeFunction( runCommunication, "gpu communication" )
                           << Sweep( [&]( IBlock *b ) { cudaLbKernel.outer( b ); }, "LBM @ outer" );
      }
      gpuTimeLoop.run();
      cuda::fieldCpy<GeneratedPdfField_T, cuda::GPUField<real_t>> (blocks, pdfShadowCPU, pdfGpuFieldId);

      // Compare all three versions
      auto errorCPU = real_t(0);
      auto errorGPU = real_t(0);

      for( auto & block : *blocks )
      {
         auto native = block.getData<NativePdfField_T>( pdfFieldNativeId );
         auto cpu = block.getData<GeneratedPdfField_T >( pdfFieldGeneratedId );
         auto gpu = block.getData<GeneratedPdfField_T>( pdfShadowCPU );

         WALBERLA_FOR_ALL_CELLS_XYZ(native,
            for(cell_idx_t f = 0; f < cell_idx_c(NativeLatticeModel_T::Stencil::Q); ++f )
            {
               errorCPU += std::abs( native->get( x, y, z, f ) - cpu->get( x, y, z, f ));
               errorGPU += std::abs( native->get( x, y, z, f ) - gpu->get( x, y, z, f ));
            }
         )
      }
      mpi::reduceInplace(errorCPU, mpi::SUM);
      mpi::reduceInplace(errorGPU, mpi::SUM);
      auto domainBB = blocks->getDomainCellBB();
      errorCPU /= real_c(domainBB.numCells());
      errorGPU /= real_c(domainBB.numCells());
      WALBERLA_LOG_RESULT_ON_ROOT("CPU Error " << errorCPU );
      WALBERLA_LOG_RESULT_ON_ROOT("GPU Error " << errorGPU );
      WALBERLA_CHECK_FLOAT_EQUAL(errorCPU, real_c(0.0));
      WALBERLA_CHECK_FLOAT_EQUAL(errorGPU, real_c(0.0));
   }

   return 0;
}