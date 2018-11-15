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
#include "cuda/ParallelStreams.h"
#include "core/timing/TimingPool.h"
#include "core/timing/RemainingTimeLogger.h"
#include "cuda/AddGPUFieldToStorage.h"
#include "cuda/communication/UniformGPUScheme.h"
#include "cuda/DeviceSelectMPI.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "domain_decomposition/SharedSweep.h"

#include "UniformGridGPU_LatticeModel.h"
#include "UniformGridGPU_LbKernel.h"
#include "UniformGridGPU_PackInfo.h"
#include "UniformGridGPU_UBB.h"
#include "UniformGridGPU_NoSlip.h"


using namespace walberla;

using LatticeModel_T = lbm::UniformGridGPU_LatticeModel;

using Stencil_T = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;
using CommScheme_T = cuda::communication::UniformGPUScheme<CommunicationStencil_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;



int main( int argc, char **argv )
{
   mpi::Environment env( argc, argv );
   cuda::selectDeviceBasedOnMpiRank();

   for( auto cfg = python_coupling::configBegin( argc, argv ); cfg != python_coupling::configEnd(); ++cfg )
   {
      auto config = *cfg;
      auto blocks = blockforest::createUniformBlockGridFromConfig( config );

      // Reading parameters
      auto parameters = config->getOneBlock( "Parameters" );
      const real_t omega = parameters.getParameter<real_t>( "omega", real_c( 1.4 ));
      const uint_t timesteps = parameters.getParameter<uint_t>( "timesteps", uint_c( 50 ));
      const Vector3<real_t> initialVelocity = parameters.getParameter< Vector3<real_t> >( "initialVelocity", Vector3<real_t>() );

      // Creating fields
      auto latticeModel = LatticeModel_T( omega );
      BlockDataID pdfFieldCpuID = lbm::addPdfFieldToStorage( blocks, "pdfs on CPU", latticeModel, initialVelocity, real_t(1), field::fzyx );
      BlockDataID pdfFieldGpuID = cuda::addGPUFieldToStorage<PdfField_T >( blocks, pdfFieldCpuID, "pdfs on GPU", true );
      BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

      // Boundaries
      const FlagUID fluidFlagUID( "Fluid" );
      auto boundariesConfig = config->getOneBlock( "Boundaries" );
      geometry::initBoundaryHandling<FlagField_T>(*blocks, flagFieldID, boundariesConfig);
      geometry::setNonBoundaryCellsToDomain<FlagField_T>(*blocks, flagFieldID, fluidFlagUID);

      lbm::UniformGridGPU_UBB ubb(blocks, pdfFieldGpuID);
      lbm::UniformGridGPU_NoSlip noSlip(blocks, pdfFieldGpuID);
      //lbm::GeneratedFixedDensity pressure(blocks, pdfFieldGpuID);

      ubb.fillFromFlagField<FlagField_T>( blocks, flagFieldID, FlagUID("UBB"), fluidFlagUID );
      noSlip.fillFromFlagField<FlagField_T>( blocks, flagFieldID, FlagUID("NoSlip"), fluidFlagUID );
      //pressure.fillFromFlagField<FlagField_T>( blocks, flagFieldID, FlagUID("pressure"), fluidFlagUID );



      // Communication setup
      bool overlapCommunication = parameters.getParameter<bool>( "overlapCommunication", true );
      bool cudaEnabledMPI = parameters.getParameter<bool>( "cudaEnabledMPI", false );

      int streamHighPriority = 0;
      int streamLowPriority = 0;
      WALBERLA_CUDA_CHECK( cudaDeviceGetStreamPriorityRange(&streamLowPriority, &streamHighPriority) );

      pystencils::UniformGridGPU_LbKernel lbKernel( pdfFieldGpuID, omega );
      lbKernel.setOuterPriority( streamHighPriority );
      CommScheme_T gpuComm( blocks, cudaEnabledMPI );
      gpuComm.addPackInfo( make_shared<pystencils::UniformGridGPU_PackInfo>( pdfFieldGpuID ));


      auto defaultStream = cuda::StreamRAII::newPriorityStream( streamLowPriority );
      auto innerOuterStreams = cuda::ParallelStreams( streamHighPriority );
      auto boundaryOuterStreams = cuda::ParallelStreams( streamHighPriority );
      auto boundaryInnerStreams = cuda::ParallelStreams( streamHighPriority );

      auto overlapTimeStep = [&]()
      {
         auto innerOuterSection = innerOuterStreams.parallelSection( defaultStream );

         innerOuterSection.run([&]( auto innerStream )
         {
            for( auto &block: *blocks )
            {
               {
                  auto p = boundaryInnerStreams.parallelSection( innerStream );
                  p.run( [&block, &ubb]( cudaStream_t s ) { ubb.inner( &block, s ); } );
                  p.run( [&block, &noSlip]( cudaStream_t s ) { noSlip.inner( &block, s ); } );
               }
               lbKernel.inner( &block, innerStream );
            }
         });

         innerOuterSection.run([&]( auto outerStream )
         {
            gpuComm( outerStream );
            for( auto &block: *blocks )
            {
               {
                  auto p = boundaryOuterStreams.parallelSection( outerStream );
                  p.run( [&block, &ubb]( cudaStream_t s ) { ubb.outer( &block, s ); } );
                  p.run( [&block, &noSlip]( cudaStream_t s ) { noSlip.outer( &block, s ); } );
               }
               lbKernel.outer( &block, outerStream );
            }
         });
      };


      auto boundaryStreams = cuda::ParallelStreams( streamHighPriority );
      auto normalTimeStep = [&]()
      {
         gpuComm();
         for( auto &block: *blocks )
         {
            {
               auto p = boundaryStreams.parallelSection( defaultStream );
               p.run( [&block, &ubb]( cudaStream_t s ) { ubb( &block, s ); } );
               p.run( [&block, &noSlip]( cudaStream_t s ) { noSlip( &block, s ); } );
            }
            lbKernel( &block );
         }
      };

      SweepTimeloop timeLoop( blocks->getBlockStorage(), timesteps );
      std::function<void()> timeStep = overlapCommunication ? std::function<void()>( overlapTimeStep ) :
                                                              std::function<void()>( normalTimeStep );
      timeLoop.add() << BeforeFunction( timeStep  )
                     << Sweep( []( IBlock * ) {}, "time step" );

      // VTK
      uint_t vtkWriteFrequency = parameters.getParameter<uint_t>( "vtkWriteFrequency", 0 );
      if( vtkWriteFrequency > 0 )
      {
         auto vtkOutput = vtk::createVTKOutput_BlockData( *blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                          "simulation_step", false, true, true, false, 0 );
         vtkOutput->addCellDataWriter(
                 make_shared<lbm::VelocityVTKWriter<LatticeModel_T> >( pdfFieldCpuID, "Velocity" ));
         vtkOutput->addCellDataWriter( make_shared<lbm::DensityVTKWriter<LatticeModel_T> >( pdfFieldCpuID, "Density" ));
         vtkOutput->addBeforeFunction(
                 cuda::fieldCpyFunctor<PdfField_T, cuda::GPUField<real_t> >( blocks, pdfFieldCpuID, pdfFieldGpuID ));
         timeLoop.addFuncAfterTimeStep( vtk::writeFiles( vtkOutput ), "VTK Output" );
      }

      auto remainingTimeLoggerFrequency = parameters.getParameter< double >( "remainingTimeLoggerFrequency", 3.0 ); // in seconds
      timeLoop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeLoop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "remaining time logger" );

      lbm::PerformanceLogger<FlagField_T> performanceLogger(blocks, flagFieldID, fluidFlagUID, 500);
      timeLoop.addFuncAfterTimeStep( performanceLogger, "remaining time logger" );

      timeLoop.run();
   }

   return 0;
}