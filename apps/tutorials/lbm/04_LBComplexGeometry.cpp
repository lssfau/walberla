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
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/all.h"

#include "domain_decomposition/all.h"

#include "field/all.h"

#include "geometry/all.h"

#include "gui/all.h"

#include "lbm/all.h"

#include "mesh_common/DistanceComputations.h"
#include "mesh_common/DistanceFunction.h"
#include "mesh_common/MatrixVectorOperations.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/distance_octree/DistanceOctree.h"
#include "mesh_common/vtk/CommonDataSources.h"
#include "mesh_common/vtk/VTKMeshWriter.h"
#include "mesh/blockforest/BlockExclusion.h"
#include "mesh/blockforest/BlockForestInitialization.h"
#include "mesh/boundary/BoundaryInfo.h"
#include "mesh/boundary/BoundaryLocation.h"
#include "mesh/boundary/BoundaryLocationFunction.h"
#include "mesh/boundary/BoundarySetup.h"
#include "mesh/boundary/BoundaryUIDFaceDataSource.h"
#include "mesh/boundary/ColorToBoundaryMapper.h"

#include "timeloop/all.h"

namespace walberla
{
uint_t numGhostLayers = uint_t(1);

//! [typedefs]
using LatticeModel_T         = lbm::D3Q27< lbm::collision_model::SRT >;
using Stencil_T              = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;
//! [typedefs]

using PdfField_T = lbm::PdfField< LatticeModel_T >;

using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;

//! [vertexToFaceColor]
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
//! [vertexToFaceColor]

int main(int argc, char** argv)
{
   walberla::Environment walberlaEnv(argc, argv);

   mpi::MPIManager::instance()->useWorldComm();

   ///////////////////////
   /// PARAMETER INPUT ///
   ///////////////////////

   // read general simulation parameters
   auto parameters = walberlaEnv.config()->getOneBlock("Parameters");

   real_t omega = parameters.getParameter< real_t >("omega", real_c(1.4));
   const Vector3< real_t > initialVelocity =
      parameters.getParameter< Vector3< real_t > >("initialVelocity", Vector3< real_t >());
   const uint_t timesteps = parameters.getParameter< uint_t >("timesteps", uint_c(10));

   const real_t remainingTimeLoggerFrequency =
      parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(3.0)); // in seconds

   //! [parseDomainParameters]
   // read domain parameters
   auto domainParameters = walberlaEnv.config()->getOneBlock("DomainSetup");

   std::string meshFile = domainParameters.getParameter< std::string >("meshFile");
   //! [parseDomainParameters]

   Vector3< uint_t > domainScaling =
      domainParameters.getParameter< Vector3< uint_t > >("domainScaling", Vector3< uint_t >(1));

   const real_t dx = domainParameters.getParameter< real_t >("dx", real_t(1));
   const Vector3< bool > periodicity =
      domainParameters.getParameter< Vector3< bool > >("periodic", Vector3< bool >(true));
   const Vector3< uint_t > cellsPerBlock = domainParameters.getParameter< Vector3< uint_t > >("cellsPerBlock");

   ////////////////////
   /// PROCESS MESH ///
   ////////////////////

   WALBERLA_LOG_INFO_ON_ROOT("Using mesh from " << meshFile << ".")

   //! [readMesh]
   // read in mesh with vertex colors on a single process and broadcast it
   auto mesh = make_shared< mesh::TriangleMesh >();
   mesh->request_vertex_colors();
   mesh::readAndBroadcast(meshFile, *mesh);
   //! [readMesh]

   // color faces according to vertices
   vertexToFaceColor(*mesh, mesh::TriangleMesh::Color(255, 255, 255));

   //! [triDist]
   // add information to mesh that is required for computing signed distances from a point to a triangle
   auto triDist = make_shared< mesh::TriangleDistance< mesh::TriangleMesh > >(mesh);
   //! [triDist]

   //! [octree]
   // building distance octree
   auto distanceOctree = make_shared< mesh::DistanceOctree< mesh::TriangleMesh > >(triDist);
   //! [octree]

   WALBERLA_LOG_INFO_ON_ROOT("Octree has height " << distanceOctree->height())

   //! [octreeVTK]
   // write distance octree to file
   WALBERLA_ROOT_SECTION()
   {
      distanceOctree->writeVTKOutput("distanceOctree");
   }
   //! [octreeVTK]

   ///////////////////////////
   /// CREATE BLOCK FOREST ///
   ///////////////////////////

   //! [aabb]
   auto aabb = computeAABB(*mesh);
   aabb.scale(domainScaling);
   aabb.setCenter(aabb.center() + 0.2 * Vector3< real_t >(aabb.xSize(), 0, 0));
   //! [aabb]

   //! [bfc]
   // create and configure block forest creator
   mesh::ComplexGeometryStructuredBlockforestCreator bfc(aabb, Vector3< real_t >(dx),
                                                         mesh::makeExcludeMeshInterior(distanceOctree, dx));

   bfc.setPeriodicity(periodicity);
   //! [bfc]


   //! [blockForest]
   // create block forest
   auto blocks = bfc.createStructuredBlockForest(cellsPerBlock);
   //! [blockForest]

   ////////////////////////////////////
   /// CREATE AND INITIALIZE FIELDS ///
   ////////////////////////////////////

   // create fields
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::SRT(omega));
   BlockDataID pdfFieldId =
      lbm::addPdfFieldToStorage(blocks, "pdf field", latticeModel, initialVelocity, real_t(1), numGhostLayers);
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", numGhostLayers);

   /////////////////////////
   /// BOUNDARY HANDLING ///
   /////////////////////////

   //! [DefaultBoundaryHandling]
   // create and initialize boundary handling
   const FlagUID fluidFlagUID("Fluid");

   auto boundariesConfig = walberlaEnv.config()->getOneBlock("Boundaries");

   typedef lbm::DefaultBoundaryHandlingFactory< LatticeModel_T, FlagField_T > BHFactory;

   BlockDataID boundaryHandlingId = BHFactory::addBoundaryHandlingToStorage(
      blocks, "boundary handling", flagFieldId, pdfFieldId, fluidFlagUID,
      boundariesConfig.getParameter< Vector3< real_t > >("velocity0", Vector3< real_t >()),
      boundariesConfig.getParameter< Vector3< real_t > >("velocity1", Vector3< real_t >()),
      boundariesConfig.getParameter< real_t >("pressure0", real_c(1.0)),
      boundariesConfig.getParameter< real_t >("pressure1", real_c(1.0)));

   //! [DefaultBoundaryHandling]

   //! [colorToBoundary]
   // set NoSlip UID to boundaries that we colored
   mesh::ColorToBoundaryMapper< mesh::TriangleMesh > colorToBoundaryMapper(
      (mesh::BoundaryInfo(BHFactory::getNoSlipBoundaryUID())));
   colorToBoundaryMapper.set(mesh::TriangleMesh::Color(255, 255, 255),
                             mesh::BoundaryInfo(BHFactory::getNoSlipBoundaryUID()));

   // mark boundaries
   auto boundaryLocations = colorToBoundaryMapper.addBoundaryInfoToMesh(*mesh);
   //! [colorToBoundary]

   //! [VTKMesh]
   // write mesh info to file
   mesh::VTKMeshWriter< mesh::TriangleMesh > meshWriter(mesh, "meshBoundaries", 1);
   meshWriter.addDataSource(make_shared< mesh::BoundaryUIDFaceDataSource< mesh::TriangleMesh > >(boundaryLocations));
   meshWriter.addDataSource(make_shared< mesh::ColorFaceDataSource< mesh::TriangleMesh > >());
   meshWriter.addDataSource(make_shared< mesh::ColorVertexDataSource< mesh::TriangleMesh > >());
   meshWriter();
   //! [VTKMesh]

   //! [boundarySetup]
   // voxelize mesh
   mesh::BoundarySetup boundarySetup(blocks, makeMeshDistanceFunction(distanceOctree), numGhostLayers);

   // set fluid cells
   boundarySetup.setDomainCells< BHFactory::BoundaryHandling >(boundaryHandlingId, mesh::BoundarySetup::OUTSIDE);

   // set up inflow/outflow/wall boundaries from DefaultBoundaryHandlingFactory
   geometry::initBoundaryHandling< BHFactory::BoundaryHandling >(*blocks, boundaryHandlingId, boundariesConfig);
   // set up obstacle boundaries from file
   boundarySetup.setBoundaries< BHFactory::BoundaryHandling >(
      boundaryHandlingId, makeBoundaryLocationFunction(distanceOctree, boundaryLocations), mesh::BoundarySetup::INSIDE);
   //! [boundarySetup]

   //////////////////////////////////
   /// SET UP SWEEPS AND TIMELOOP ///
   //////////////////////////////////

   // create time loop
   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication(blocks);
   communication.addPackInfo(make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >(pdfFieldId));

   // add LBM sweep and communication to time loop
   timeloop.add() << BeforeFunction(communication, "communication")
                  << Sweep(BHFactory::BoundaryHandling::getBlockSweep(boundaryHandlingId), "boundary handling");
   timeloop.add() << Sweep(
      makeSharedSweep(lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >(pdfFieldId, flagFieldId, fluidFlagUID)),
      "LB stream & collide");

   // LBM stability check
   timeloop.addFuncAfterTimeStep(makeSharedFunctor(field::makeStabilityChecker< PdfField_T, FlagField_T >(
                                    walberlaEnv.config(), blocks, pdfFieldId, flagFieldId, fluidFlagUID)),
                                 "LBM stability check");

   // log remaining time
   timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
                                 "remaining time logger");

   //////////////////
   /// VTK OUTPUT ///
   //////////////////

   // add VTK output to time loop
   auto VTKParams = walberlaEnv.config()->getBlock("VTK");
   uint_t vtkWriteFrequency = VTKParams.getBlock("fluid_field").getParameter("writeFrequency", uint_t(0));
   auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "fluid_field", vtkWriteFrequency, uint_t(0), false,
                                                   "vtk_out", "simulation_step", false, true, true, false, 0);

   field::FlagFieldCellFilter< FlagField_T > fluidFilter(flagFieldId);
   fluidFilter.addFlag(fluidFlagUID);
   vtkOutput->addCellInclusionFilter(fluidFilter);

   auto velocityWriter = make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >(pdfFieldId, "Velocity");
   auto densityWriter  = make_shared< lbm::DensityVTKWriter< LatticeModel_T, float > >(pdfFieldId, "Density");
   vtkOutput->addCellDataWriter(velocityWriter);
   vtkOutput->addCellDataWriter(densityWriter);

   timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");

   // create adaptors, so that the GUI also displays density and velocity
   // adaptors are like fields with the difference that they do not store values
   // but calculate the values based on other fields ( here the PdfField )
   field::addFieldAdaptor< lbm::Adaptor< LatticeModel_T >::Density >(blocks, pdfFieldId, "DensityAdaptor");
   field::addFieldAdaptor< lbm::Adaptor< LatticeModel_T >::VelocityVector >(blocks, pdfFieldId, "VelocityAdaptor");

   //////////////////////
   /// RUN SIMULATION ///
   //////////////////////

   if (parameters.getParameter< bool >("useGui", false))
   {
      GUI gui(timeloop, blocks, argc, argv);
      lbm::connectToGui< LatticeModel_T >(gui);
      gui.run();
   }
   else
      timeloop.run();

   return EXIT_SUCCESS;
}
} // namespace walberla

int main(int argc, char** argv) { walberla::main(argc, argv); }
