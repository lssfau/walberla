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
//! \file 04_LBComplexGeometry.cpp
//! \author Brendan Waters <brendan.waters@sydney.edu.au>
//! \author Girish Kumatagi <girish.h.kumatagi@fau.de>
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/all.h"

#include "domain_decomposition/all.h"

#include "field/all.h"

#include "geometry/all.h"

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   #include "gpu/AddGPUFieldToStorage.h"
   #include "gpu/DeviceSelectMPI.h"
   #include "gpu/FieldCopy.h"
   #include "gpu/GPUWrapper.h"
   #include "gpu/HostFieldAllocator.h"
   #include "gpu/ParallelStreams.h"
   #include "gpu/communication/UniformGPUScheme.h"
   #include "gpu/timeloop/DeviceSweepTimeloop.h"
   #include "gpu/timing/DeviceTimingPool.h"
#endif

#include "lbm/all.h"

#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/evaluation/PerformanceEvaluation.h"
#include "lbm_generated/communication/UniformGeneratedPdfPackInfo.h"

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   #include "lbm_generated/gpu/UniformGeneratedGPUPdfPackInfo.h"
   #include "lbm_generated/gpu/GPUPdfField.h"
   #include "lbm_generated/gpu/AddToStorage.h"
#endif

#include "mesh_common/DistanceComputations.h"
#include "mesh_common/DistanceFunction.h"
#include "mesh_common/MatrixVectorOperations.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/distance_octree/DistanceOctree.h"
#include "mesh_common/vtk/CommonDataSources.h"
#include "mesh_common/vtk/VTKMeshWriter.h"
#include "mesh/blockforest/BlockForestInitialization.h"
#include "mesh/boundary/BoundaryInfo.h"
#include "mesh/boundary/BoundaryLocation.h"
#include "mesh/boundary/BoundaryLocationFunction.h"
#include "mesh/boundary/BoundarySetup.h"
#include "mesh/boundary/BoundaryUIDFaceDataSource.h"
#include "mesh/boundary/ColorToBoundaryMapper.h"
#include "mesh_common/MeshBodyWallDistanceCallback.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/DictWrapper.h"
#include "python_coupling/PythonCallback.h"

#include "timeloop/all.h"

#include "InfoHeader.h"


namespace walberla
{

constexpr uint_t FieldGhostLayer{1};

using StorageSpecification_T     = lbm::LBComplexGeometryStorageSpecification;
using LBMCommunicationStencil_T  = StorageSpecification_T::CommunicationStencil;
using PdfField_T                 = lbm_generated::PdfField< StorageSpecification_T >;

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   using GPUPdfField_T              = lbm_generated::GPUPdfField< StorageSpecification_T >;
   using LBMPackInfo_T              = lbm_generated::UniformGeneratedGPUPdfPackInfo< GPUPdfField_T >;
   using gpu::communication::UniformGPUScheme;
#else
   using LBMPackInfo_T              = lbm_generated::UniformGeneratedPdfPackInfo< PdfField_T >;
   using blockforest::communication::UniformBufferedScheme;
#endif

using SweepCollection_T = lbm::LBComplexGeometrySweepCollection;

using Mesh_T = mesh::TriangleMesh;

using flag_t      = walberla::uint32_t;
using FlagField_T = FlagField< flag_t >;

template< typename T > inline flag_t flag_c( T t ) { return numeric_cast< flag_t >(t); }

using BoundaryCollection_T = lbm::LBComplexGeometryBoundaryCollection< FlagField_T >;

const FlagUID FluidFlagUID("Fluid");
const FlagUID NoSlipFlagUID("NoSlip");
const FlagUID NoSlipQBBFlagUID("NoSlipQBB");
const FlagUID MeshObjNoSlipQBBFlagUID("ObjNoSlipQBB");
const FlagUID FreeSlipFlagUID("FreeSlip");
const FlagUID UBBFlagUID("UBB");
const FlagUID OutflowFlagUID("Outflow");

const BoundaryUID NoSlipBoundaryUID("NoSlip");

std::unordered_map<std::string, FlagUID> flagMap = {  {"NoSlip", NoSlipFlagUID},
                                                      {"NoSlipQBB", NoSlipQBBFlagUID},
                                                      {"ObjNoSlipQBB", MeshObjNoSlipQBBFlagUID},
                                                      {"FreeSlip", FreeSlipFlagUID},
                                                      {"UBB", UBBFlagUID},
                                                      {"Outflow", OutflowFlagUID} 
                                                      // add more as needed
                                                   };

template< typename MeshType >
void vertexToFaceColor(MeshType& mesh, const typename MeshType::Color& defaultColor)
{
   WALBERLA_CHECK(mesh.has_vertex_colors())
   mesh.request_face_colors();

   for (auto faceIt = mesh.faces_begin(); faceIt != mesh.faces_end(); ++faceIt)
   {
      typename MeshType::Color vertexColor;

      bool useVertexColor = true;

      auto vertexIt = mesh.fv_iter(*faceIt);
      WALBERLA_ASSERT(vertexIt.is_valid())

      vertexColor = mesh.color(*vertexIt);

      ++vertexIt;
      while (vertexIt.is_valid() && useVertexColor)
      {
         if (vertexColor != mesh.color(*vertexIt)) useVertexColor = false;
         ++vertexIt;
      }

      mesh.set_color(*faceIt, useVertexColor ? vertexColor : defaultColor);
   }
}

class BoundarySetup
{
   struct DirectionHandler {
      static uint_t dir2idx(const char& direction) {
            switch (direction) {
               case 'E': case 'W': return 0;
               case 'N': case 'S': return 1;
               case 'T': case 'B': return 2;
               default:
                  WALBERLA_ABORT("Invalid direction provided");
            }
      }
   };

public:

   explicit BoundarySetup( const Config::BlockHandle &config, const Vector3<bool> & periodicity ):
            config_(config), periodicity_(periodicity){}

   void getMappedUID(FlagUID & UID, const std::string & handle) const
   {
      auto it = flagMap.find( handle );
      if (it != flagMap.end()) {
         UID = it->second;
      } else {
         WALBERLA_ABORT("Specified boundary condition: " << handle << " not defined.");
      }
   }

   template< typename FlagField_T >
   void fillFlagFieldFromConfig( const std::shared_ptr<blockforest::StructuredBlockForest>& blocks,
                                 const BlockDataID & flagFieldID ) const
   {

      WALBERLA_LOG_INFO_ON_ROOT("Setting domain boundaries ...")

      Config::Block boundaryBlock( config_.getKey() );

      Config::Blocks borderBlocks;
      config_.getBlocks( "Border", borderBlocks );

      for( auto border = borderBlocks.begin(); border != borderBlocks.end(); ++border )
      {
         const std::string direction  = border->getParameter<std::string>("direction");
         const std::string flagHandle = border->getParameter<std::string>("flag");
         
         std::istringstream dirStream(direction);
         std::string dir;

         while ( std::getline(dirStream, dir, ',') ) {
               // Trim whitespace
               dir.erase(remove_if(dir.begin(), dir.end(), isspace), dir.end());

               if (dir.size() != 1) {
                  WALBERLA_ABORT("Invalid direction provided: " << dir);
               }

               if ( periodicity_[ DirectionHandler::dir2idx(dir[0]) ] ){
                  WALBERLA_LOG_WARNING_ON_ROOT("Specified boundary in direction " << direction 
                                                << " conflicts with periodicity. Skipping boundary condition: "  << flagHandle);
                  continue;
               }

               FlagUID userDefinedFlagUID;
               getMappedUID(userDefinedFlagUID, flagHandle);
               
               addBoundary(boundaryBlock, dir, userDefinedFlagUID.getIdentifier());
         }
      }

      Config::BlockHandle boundariesConfig(&boundaryBlock);

      geometry::initBoundaryHandling< FlagField_T >( *blocks, flagFieldID, boundariesConfig);
   }
   
   // Might need to register flag for object if not already registered for field
   template< typename FlagField_T >
   void registerFlagFieldForMeshObject( const std::shared_ptr<blockforest::StructuredBlockForest>& blocks,
                                        const BlockDataID& flagFieldId, const FlagUID & UID ) const
   {
      for (auto& block : *blocks)
      {
         auto * flagField = block.getData<FlagField_T>(flagFieldId);
         if ( !flagField->flagExists(UID))
         flagField->registerFlag(UID);
      }
   }

   
private:
   static void addBoundary( Config::Block & block, const std::string & dir, const std::string & flag ) 
   {
      auto & border = block.createBlock("Border");
      border.addParameter("direction", dir);
      border.addParameter("walldistance", "-1");
      border.addParameter("flag", flag);
   }

   const Config::BlockHandle config_;
   const Vector3<bool> periodicity_;
};

auto velocityCallback= [](const Cell & /*pos*/, const shared_ptr< StructuredBlockForest > & /*structuredBlockforest*/, IBlock & /*block*/, const real_t & velocity )
{
   const real_t v = velocity;

   const Vector3< real_t > result(real_c(v), real_c(0.0), real_c(0.0));

   return result;
};

auto wallDistanceCallback= [](const Cell& /*fluid*/, const Cell& /*boundary*/, const shared_ptr< StructuredBlockForest >& /*SbF*/, IBlock& /*block*/, const real_t & q )
{
   return q;
};

int main(int argc, char** argv)
{
   Environment env(argc, argv);
   if (!env.config()) { WALBERLA_ABORT("No configuration file specified!"); }

   mpi::MPIManager::instance()->useWorldComm();

   #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      gpu::selectDeviceBasedOnMpiRank();
   #endif

   const std::string input_filename(argv[1]);
   const bool inputIsPython = string_ends_with(input_filename, ".py");

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_WORLD_BARRIER()

      #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         WALBERLA_GPU_CHECK(gpuPeekAtLastError())
      #endif
      
      auto config = *cfg;
      logging::configureLogging(config);
      
      /*------------------------------------------------------------------------------------------------------*
       *-----------------------------------CONFIGURE FROM INPUT PARAMETERS------------------------------------*
       *------------------------------------------------------------------------------------------------------*/

      auto parameters             = config->getOneBlock("Parameters");
      const real_t simulationTime = parameters.getParameter< real_t >("simTime", real_c(60));
      const real_t velocity_LB    = parameters.getParameter<real_t>("uLB", real_c( 0.05 ) );
      const real_t velocity_SI    = parameters.getParameter<real_t>("velocity_SI");
      const real_t viscosity_SI   = parameters.getParameter<real_t>("viscosity_SI");
      const real_t Href_SI        = parameters.getParameter<real_t>("Href_SI");

      // Distance from fluid cell centre to axis aligned wall boundary (i.e. Not for mesh object) 
      const real_t q              = parameters.getParameter<real_t>("domainBoundaryWallDistance", real_c( 0.5 ));
      //------------------------------------------------------------------------------------------------------//

      auto geometryParameters    = config->getOneBlock( "Geometry" );
      const std::string meshFile = geometryParameters.getParameter< std::string >( "meshFile" );
      
      //------------------------------------------------------------------------------------------------------//
      
      auto domainParameters = config->getOneBlock("DomainSetup");
      const Vector3< uint_t > domainScaling = domainParameters.getParameter< Vector3< uint_t > >("domainScaling", Vector3< uint_t >(1));
      const Vector3< bool >   periodicity   = domainParameters.getParameter< Vector3< bool > >("periodic", Vector3< bool >(true));
      const Vector3< uint_t > cellsPerBlock = domainParameters.getParameter< Vector3< uint_t > >("cellsPerBlock");
      const real_t            dx_SI         = domainParameters.getParameter< real_t >("dx_SI", real_t(1));
      
      /*-----------------------------------------------------------------------------------------------------*/
      
      const real_t C_l   = dx_SI;
      const real_t dx_LB = dx_SI/C_l; 

      const real_t dt_SI = (velocity_LB / velocity_SI) * dx_SI;
      const real_t C_t   = dt_SI;
      const real_t dt_LB = dt_SI/C_t;

      const real_t speedOfSound = real_c(real_c(1.0) / std::sqrt( real_c(3.0) ));
      const real_t machNumber = velocity_LB / speedOfSound;

      const real_t lattice_viscosity = viscosity_SI *  C_t/(C_l * C_l);    // Physical units are m^2/s

      const real_t omega = lbm::collision_model::omegaFromViscosity(lattice_viscosity);

      const uint_t timesteps = uint_c( std::ceil( simulationTime/dt_SI ) ) + 1;

      /*-----------------------------------------------------------------------------------------------------*/
      std::ostringstream oss;
      oss <<   "- Physical Parameters:"
            << "\n   + dx (m):               " << dx_SI
            << "\n   + dt (s):               " << dt_SI
            << "\n   + Velocity (m/s):       " << velocity_SI
            << "\n   + Viscosity (m^2/s):    " << viscosity_SI
            << "\n   + Reference Height (m): " << Href_SI;

      oss   << "\n- Scaled LBM Parameters:"
            << "\n   + dx:               " << dx_LB
            << "\n   + dt:               " << dt_LB
            << "\n   + Velocity:         " << velocity_LB
            << "\n   + Viscocity:        " << lattice_viscosity
            << "\n   + Reference Height: " << Href_SI/C_l
            << "\n   + Mach Number:      " << machNumber
            << "\n   + Relaxation Rate:  " << omega;

      oss   << "\n- Time Parameters:"
            << "\n   + No. Course Timesteps:        " << timesteps
            << "\n   + Simulation Time (s) [actual/specified]:  [" << real_t(timesteps)*dt_SI <<"/"<<simulationTime<<"]";

      std::string str = oss.str();

      WALBERLA_LOG_INFO_ON_ROOT( "Physical/LBM Parameter Conversion:\n" << str);
      /*------------------------------------------------------------------------------------------------------
      *-------------------------------CREATE DISTANCE OCTREE FROM MESHFILE-----------------------------------
      *------------------------------------------------------------------------------------------------------*/
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT( meshFile );

      auto mesh = make_shared< Mesh_T >();

      mesh->request_vertex_colors();
      WALBERLA_LOG_INFO_ON_ROOT( "Loading mesh" );
      mesh::readAndBroadcast( meshFile, *mesh);

      vertexToFaceColor( *mesh, Mesh_T::Color(255,255,255) );

      auto bunnyAABB = computeAABB(*mesh);
      
      // Scale to unit size, then to match sim parameters
      const Vector3<real_t> mesh_scaling_factor { 1/real_c(bunnyAABB.ySize()) * Href_SI/C_l};
      mesh::scale(*mesh, mesh_scaling_factor );

      WALBERLA_LOG_INFO_ON_ROOT( "Adding distance info to mesh" );
      auto triDist = make_shared< mesh::TriangleDistance< Mesh_T > >( mesh );
      WALBERLA_LOG_INFO_ON_ROOT( "Building distance octree" );
      auto distanceOctree = make_shared< mesh::DistanceOctree< Mesh_T > >( triDist );
      WALBERLA_LOG_INFO_ON_ROOT( "done. Octree has height " << distanceOctree->height() );

      WALBERLA_ROOT_SECTION(){ distanceOctree->writeVTKOutput("distanceOctree"); }

      auto aabb = computeAABB(*mesh);
      aabb.scale(domainScaling);
      aabb.setCenter(aabb.center() + 0.2 * Vector3< real_t >(aabb.xSize(), 0, 0));

      /*------------------------------------------------------------------------------------------------------
      *-------------------------------------------CREATE BLOCK FOREST-----------------------------------------
      *------------------------------------------------------------------------------------------------------*/
      mesh::ComplexGeometryStructuredBlockforestCreator bfc( aabb, Vector3< real_t >(dx_LB),
                                                             mesh::makeExcludeMeshInterior(distanceOctree, dx_LB));

      bfc.setPeriodicity(periodicity);
      auto blocks = bfc.createStructuredBlockForest( cellsPerBlock );

      /*------------------------------------------------------------------------------------------------------
      *------------------------------------------- DATA FIELDS ----------------------------------------------
      *------------------------------------------------------------------------------------------------------*/

      // CPU Fields
      const StorageSpecification_T StorageSpec = StorageSpecification_T();
      const BlockDataID pdfFieldCpuId      = lbm_generated::addPdfFieldToStorage(blocks, "pdfs", StorageSpec, FieldGhostLayer, field::fzyx);

      const BlockDataID velocityFieldCpuId = field::addToStorage< VelocityField_T >(blocks, "velocity", real_c(0.0), field::fzyx, FieldGhostLayer);
      const BlockDataID densityFieldCpuId  = field::addToStorage< ScalarField_T   >(blocks, "density", real_c(1.0), field::fzyx, FieldGhostLayer);

      const BlockDataID flagFieldId        = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", FieldGhostLayer);

      //------------------------------------------------------------------------------------------------------//
      // Initialise fields before transferring them to GPU (if required)
      //------------------------------------------------------------------------------------------------------//
      
      #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         const BlockDataID pdfFieldGpuId      = lbm_generated::addGPUPdfFieldToStorage< PdfField_T >(blocks, pdfFieldCpuId, StorageSpec, "pdfs_GPU", true);
         const BlockDataID velocityFieldGpuId = gpu::addGPUFieldToStorage< VelocityField_T >(blocks, velocityFieldCpuId,  "velocity_GPU" , true);
         const BlockDataID densityFieldGpuId  = gpu::addGPUFieldToStorage< ScalarField_T   >(blocks, densityFieldCpuId,  "density_GPU", true);
      #endif

      //------------------------------------------------------------------------------------------------------//
      
      const Cell innerOuterSplit = Cell( parameters.getParameter< Vector3<cell_idx_t> >( "innerOuterSplit", Vector3<cell_idx_t>(1, 1, 1))  );
      
      #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         const Vector3< int32_t > gpuBlockSize = parameters.getParameter< Vector3< int32_t > >( "gpuBlockSize", Vector3< int32_t >(32, 1, 1) );

         SweepCollection_T sweepCollection( blocks, densityFieldGpuId, pdfFieldGpuId, velocityFieldGpuId,
                                          gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2], omega, innerOuterSplit);

         int streamHighPriority = 0;
         int streamLowPriority  = 0;

         WALBERLA_GPU_CHECK(gpuDeviceGetStreamPriorityRange(&streamLowPriority, &streamHighPriority))
         
         sweepCollection.setOuterPriority(streamHighPriority);
      #else 
         SweepCollection_T sweepCollection( blocks, densityFieldCpuId, pdfFieldCpuId, velocityFieldCpuId, omega, innerOuterSplit);
      #endif
     

      for (auto& block : *blocks){
         sweepCollection.initialise(&block);
      }

      /*------------------------------------------------------------------------------------------------------
      *------------------------------------------- BOUNDARY HANDLING -----------------------------------------
      *------------------------------------------------------------------------------------------------------*/
      mesh::ColorToBoundaryMapper< Mesh_T > colorToBoundaryMapper(( mesh::BoundaryInfo( NoSlipBoundaryUID ) ));
      
      auto boundaryLocations = colorToBoundaryMapper.addBoundaryInfoToMesh( *mesh );

      mesh::VTKMeshWriter< Mesh_T > meshWriter( mesh, "meshBoundaries", 1 );
      meshWriter.addDataSource( make_shared< mesh::BoundaryUIDFaceDataSource< Mesh_T > >( boundaryLocations ) );
      meshWriter.addDataSource( make_shared< mesh::ColorFaceDataSource< Mesh_T > >() );
      meshWriter.addDataSource( make_shared< mesh::ColorVertexDataSource< Mesh_T > >() );
      meshWriter();

      WALBERLA_LOG_INFO_ON_ROOT( "Setting mesh/domain BC Flags" );

      auto boundariesConfig = config->getOneBlock("Boundaries");

      // Set domain boundaries
      BoundarySetup boundarySetter( boundariesConfig, periodicity );
      boundarySetter.fillFlagFieldFromConfig<FlagField_T>( blocks, flagFieldId);

      // Set mesh object boundaries
      const auto ObjectBC = boundariesConfig.getOneBlock( "ObjectBC" );
      const std::string objectBCflagHandle = ObjectBC.getParameter< std::string >( "flag" );

      FlagUID ObjectBCFlagUID; 
      boundarySetter.getMappedUID(ObjectBCFlagUID, objectBCflagHandle);
      boundarySetter.registerFlagFieldForMeshObject<FlagField_T>( blocks, flagFieldId, ObjectBCFlagUID);
      
      mesh::BoundarySetup boundarySetup( blocks, makeMeshDistanceFunction( distanceOctree ), FieldGhostLayer );
      boundarySetup.setFlag<FlagField_T>(flagFieldId, ObjectBCFlagUID, mesh::BoundarySetup::INSIDE);

      // Set remaining cells to fluid
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, FluidFlagUID);

      WALBERLA_MPI_WORLD_BARRIER();

      // Wall distance for the mesh object
      mesh::MeshBodyWallDistance< Mesh_T > meshWallDistanceCallback( distanceOctree );
      std::function< real_t(const Cell&, const Cell&, const shared_ptr< StructuredBlockForest >&, IBlock&) >
         meshWallDistanceFunctor = meshWallDistanceCallback;

      // Wall distance for the domain boundaries
      std::function< real_t (const Cell& fluid, const Cell& boundary, const shared_ptr< StructuredBlockForest >& SbF, IBlock& block) >
      wallDistanceFunctor = std::bind(wallDistanceCallback, std::placeholders::_1, std::placeholders::_2, 
                                      std::placeholders::_3, std::placeholders::_4, q); 

      // Velocity for the inflow boundaries
      std::function< Vector3< real_t >(const Cell&, const shared_ptr< StructuredBlockForest >&, IBlock&) >
      inflowVelocity = std::bind(velocityCallback, std::placeholders::_1, std::placeholders::_2, 
                                 std::placeholders::_3, velocity_LB);      

      #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         BoundaryCollection_T boundaryCollection( blocks, flagFieldId, pdfFieldGpuId, FluidFlagUID, omega, pdfFieldCpuId, inflowVelocity, wallDistanceFunctor, meshWallDistanceFunctor);
      #else 
         BoundaryCollection_T boundaryCollection( blocks, flagFieldId, pdfFieldCpuId, FluidFlagUID, omega, inflowVelocity, wallDistanceFunctor, meshWallDistanceFunctor  );
      #endif
     
      const std::string vtkFlagStr = parameters.getParameter< std::string >("vtkFlagOutputString");
      auto vtkFlagOutput = vtk::createVTKOutput_BlockData(*blocks, vtkFlagStr, uint_t(1), FieldGhostLayer, false, "vtk_out",
                                                         "flags", false, true, true, false, 0);
      auto flagWriter = make_shared< field::VTKWriter< FlagField_T > >(flagFieldId, "flag");
      vtkFlagOutput->addCellDataWriter(flagWriter);
      vtkFlagOutput->write();

      lbm::BlockForestEvaluation<FlagField_T>( blocks, flagFieldId, FluidFlagUID ).logInfoOnRoot();

      /*------------------------------------------------------------------------------------------------------
      *--------------------------------------------- COMMUNICATION -------------------------------------------
      *------------------------------------------------------------------------------------------------------*/
      #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         const bool cudaEnabledMPI = parameters.getParameter< bool >("cudaEnabledMPI", false);
         UniformGPUScheme< LBMCommunicationStencil_T > communication(blocks, cudaEnabledMPI, false);
         communication.addPackInfo( std::make_shared< LBMPackInfo_T >( pdfFieldGpuId ) );
      #else
         UniformBufferedScheme< LBMCommunicationStencil_T > communication(blocks);
         communication.addPackInfo( std::make_shared< LBMPackInfo_T>( pdfFieldCpuId ) );
      #endif
      /*------------------------------------------------------------------------------------------------------
      *------------------------------------------------- TIMELOOP --------------------------------------------
      *------------------------------------------------------------------------------------------------------*/

      #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         DeviceSynchronizeSweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

         auto defaultStream = gpu::StreamRAII::newPriorityStream(streamLowPriority);
         timeloop.add() << BeforeFunction(communication.getCommunicateFunctor(),                         "LBM Communication")
                        << Sweep(boundaryCollection.getSweep(BoundaryCollection_T::ALL, defaultStream),  "LBM Boundary Conditions");
         timeloop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL, defaultStream),   "LBM StreamCollide");
      #else
         SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

         timeloop.add() << BeforeFunction(communication,                                  "LBM Communication")
                        << Sweep(boundaryCollection.getSweep(BoundaryCollection_T::ALL),  "LBM Boundary Conditions");
         timeloop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL),   "LBM StreamCollide");
      #endif

      //------------------------------------------------------------------------------------------------------//
      
      auto remainingTimeLoggerFrequency = parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(-1.0)); // in seconds
      if (remainingTimeLoggerFrequency > 0){
         timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency), "remaining time logger");
      }

      //------------------------------------------------------------------------------------------------------//
      
      const uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);

      if (vtkWriteFrequency > 0)
      {
         const std::string vtkStr = parameters.getParameter< std::string >("vtkOutputString");

         auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, vtkStr, vtkWriteFrequency, uint_t(0), false, "vtk_out",
                                                         "simulation_step", false, true, true, false, 0);
         
         auto velocityWriter = make_shared< field::VTKWriter< VelocityField_T > >(velocityFieldCpuId, "velocity");
         auto densityWriter  = make_shared< field::VTKWriter< ScalarField_T > >(densityFieldCpuId, "density");

         vtkOutput->addCellDataWriter(velocityWriter);
         vtkOutput->addCellDataWriter(densityWriter);

         vtkOutput->addBeforeFunction([&]() 
         {
            for (auto& block : *blocks)
               sweepCollection.calculateMacroscopicParameters(&block);

            #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
               gpu::fieldCpy< VelocityField_T, gpu::GPUField< real_t > >(blocks, velocityFieldCpuId, velocityFieldGpuId);
               gpu::fieldCpy< ScalarField_T,   gpu::GPUField< real_t > >(blocks, densityFieldCpuId,  densityFieldGpuId);
            #endif  
         });
         timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      }
      
      /*------------------------------------------------------------------------------------------------------
      *------------------------------------------------- RUN SIM --------------------------------------------
      *------------------------------------------------------------------------------------------------------*/

      #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         WALBERLA_GPU_CHECK(gpuPeekAtLastError())

         DeviceSynchronizeTimingPool timingPool;

         WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )
         WALBERLA_GPU_CHECK( gpuPeekAtLastError() )
      #else
         WcTimingPool timingPool;
      #endif

      WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")
      timeloop.run(timingPool);
      
      #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )
      #endif

      WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
      timingPool.unifyRegisteredTimersAcrossProcesses();
      timingPool.logResultOnRoot( timing::REDUCE_TOTAL, true );
   }
   return EXIT_SUCCESS;
}
} // namespace walberla

int main(int argc, char** argv) { walberla::main(argc, argv); }