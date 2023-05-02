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
//! \file SquirmerTest.cpp
//! \ingroup pe_coupling
//! \author Christian Burkard <buch@icp.uni-stuttgart.de>
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/all.h"
#include "core/SharedFunctor.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/mpi/MPIManager.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/ForceModel.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "stencil/D3Q7.h"

#include "timeloop/SweepTimeloop.h"

#include "pe/rigidbody/SquirmerFactory.h"
#include "pe/communication/rigidbody/Squirmer.h"
#include "pe/basic.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"
#include "lbm/vtk/Velocity.h"
#include "vtk/VTKOutput.h"

#include <functional>

namespace squirmer
{


///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

// PDF field, flag field & body field
typedef lbm::D3Q19<lbm::collision_model::TRT> LatticeModel_T;

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
typedef GhostLayerField<pe::BodyID, 1> PeBodyField_T;

const uint_t FieldGhostLayers = 1;

// boundary handling
typedef pe_coupling::SimpleBB<LatticeModel_T, FlagField_T> MO_BB_T;

typedef BoundaryHandling<FlagField_T, Stencil_T, MO_BB_T> BoundaryHandling_T;

using BodyTypeTuple = std::tuple<pe::Squirmer>;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag("fluid");
const FlagUID MO_BB_Flag("moving obstacle BB");
const FlagUID FormerMO_BB_Flag("former moving obstacle BB");

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

class MyBoundaryHandling {
public:

   MyBoundaryHandling(const BlockDataID &flagFieldID, const BlockDataID &pdfFieldID, const BlockDataID &bodyFieldID) :
         flagFieldID_(flagFieldID), pdfFieldID_(pdfFieldID), bodyFieldID_(bodyFieldID) {}

   BoundaryHandling_T *operator()(IBlock *const block, const StructuredBlockStorage *const storage) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID bodyFieldID_;

}; // class MyBoundaryHandling

////////////////////////////////
// VTK OUTPUT HELPER FUNCTION //
////////////////////////////////

BoundaryHandling_T *
MyBoundaryHandling::operator()(IBlock *const block, const StructuredBlockStorage *const storage) const {
   WALBERLA_ASSERT_NOT_NULLPTR(block);
   WALBERLA_ASSERT_NOT_NULLPTR(storage);

   FlagField_T *flagField = block->getData<FlagField_T>(flagFieldID_);
   PdfField_T *pdfField = block->getData<PdfField_T>(pdfFieldID_);
   PeBodyField_T *bodyField = block->getData<PeBodyField_T>(bodyFieldID_);

   const auto fluid = flagField->flagExists(Fluid_Flag) ? flagField->getFlag(Fluid_Flag) : flagField->registerFlag(
         Fluid_Flag);

   BoundaryHandling_T *handling = new BoundaryHandling_T( "fixed obstacle boundary handling", flagField, fluid,
                                                          MO_BB_T("MO_BB", MO_BB_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ) );

   handling->fillWithDomain(FieldGhostLayers);

   return handling;
}

shared_ptr<vtk::VTKOutput> createFluidFieldVTKWriter(shared_ptr<StructuredBlockForest> &blocks,
                                                     const BlockDataID &pdfFieldId, const BlockDataID &flagFieldId,
                                                     std::string fname) {
   auto pdfFieldVTKWriter = vtk::createVTKOutput_BlockData(blocks, fname, uint_c(1), uint_c(0));

   field::FlagFieldCellFilter<FlagField_T> fluidFilter(flagFieldId);
   fluidFilter.addFlag(Fluid_Flag);
   pdfFieldVTKWriter->addCellInclusionFilter(fluidFilter);

   auto velocityWriter = make_shared<lbm::VelocityVTKWriter<LatticeModel_T, float> >(pdfFieldId, "VelocityFromPDF");
   pdfFieldVTKWriter->addCellDataWriter(velocityWriter);

   return pdfFieldVTKWriter;
}

class ResetPosition {
public:

   explicit ResetPosition(const shared_ptr<StructuredBlockStorage> &blockStorage,
                          const BlockDataID &bodyStorageID,
                          real_t dt)
         : blockStorage_(blockStorage), bodyStorageID_(bodyStorageID), dt_(dt) {}

   void operator()() {
      for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt) {
         for (auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_);
              bodyIt != pe::BodyIterator::end(); ++bodyIt) {
            bodyIt->setPosition(bodyIt->getPosition() - bodyIt->getLinearVel() * dt_);
         }
      }
   }

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID &bodyStorageID_;

   real_t dt_;

}; // class ResetPosition


Vector3<real_t> analytical_velocity(Vector3<real_t> r, real_t b1, real_t b2, Vector3<real_t> e, real_t radius) {
   auto r_len = r.length();
   auto rs = r.getNormalized();
   auto frac = radius / r_len;
   auto frac2 = frac * frac;
   auto frac3 = frac2 * frac;
   auto frac4 = frac3 * frac;

   return b1 * frac3 * ((e * rs) * rs - e / real_c(3.0)) +
          (frac4 - frac2) * b2 * 0.5 * (3 * (e * rs) * (e * rs) - 1) * rs + frac4 * b2 * (e * rs) * ((e * rs) * rs - e);
}


int main(int argc, char **argv) {
   debug::enterTestMode();

   mpi::Environment env(argc, argv);

   auto processes = MPIManager::instance()->numProcesses();

   if (processes != 1 && processes != 2 && processes != 4 && processes != 8) {
      std::cerr << "Number of processes must be equal to either 1, 2, 4, or 8!" << std::endl;
      return EXIT_FAILURE;
   }

   bool shortrun = false;
   bool writevtk = false;
   for (int i = 1; i < argc; ++i) {
      if (std::strcmp(argv[i], "--shortrun") == 0) {
         shortrun = true;
         continue;
      }
      if (std::strcmp(argv[i], "--vtk") == 0) {
         writevtk = true;
         continue;
      }
   }

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   real_t R = 3;
   uint_t L = 24;
   uint_t T = 1000;
   if (shortrun) {
      L = 10;
      T = 10;
   }
   const real_t dx = real_c(1.0);
   const real_t dt = real_c(1.0);
   const real_t omega = lbm::collision_model::omegaFromViscosity(real_c(1.0 / 6.0));
   const real_t v0 = real_c(0.01);
   const real_t beta = real_c(5.0);

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////
   const uint_t XBlocks = (processes >= 2) ? uint_t(2) : uint_t(1);
   const uint_t YBlocks = (processes >= 4) ? uint_t(2) : uint_t(1);
   const uint_t ZBlocks = (processes == 8) ? uint_t(2) : uint_t(1);
   const uint_t XCells = L / XBlocks;
   const uint_t YCells = L / YBlocks;
   const uint_t ZCells = L / ZBlocks;

   // create fully periodic domain
   auto blocks = blockforest::createUniformBlockGrid(XBlocks, YBlocks, ZBlocks, XCells, YCells, ZCells, dx, true,
                                                     true, true, true);

   ////////
   // PE //
   ////////

   auto globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");
   auto ccdID = blocks->addBlockData(pe::ccd::createHashGridsDataHandling(globalBodyStorage, bodyStorageID), "CCD");
   auto fcdID = blocks->addBlockData(
         pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");

   pe::cr::HCSITS cr(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID);

   /////////////////
   // PE COUPLING //
   /////////////////

   // connect to pe
   const real_t overlap = real_c(1.5) * dx;

   if (R > real_c(L) * real_c(0.5) - overlap) {
      std::cerr << "Periodic sphere is too large and would lead to invalid PE state!" << std::endl;
      return EXIT_FAILURE;
   }

   std::function<void(void)> syncCall = std::bind(pe::syncShadowOwners<BodyTypeTuple>,
                                                      std::ref(blocks->getBlockForest()), bodyStorageID,
                                                      static_cast<WcTimingTree *>(nullptr), overlap, false);

   const auto myMat = pe::createMaterial("myMat", real_c(1), real_t(0), real_t(1), real_t(1), real_t(0), real_t(1),
                                         real_t(1), real_t(0), real_t(0));

   // create the squirmer in the middle of the domain
   const Vector3<real_t> position(real_c(L) * real_c(0.5), real_c(L) * real_c(0.5), real_c(L) * real_c(0.5));
   const Vector3<real_t> up(0.0, 0.0, 1.0);
   const Vector3<real_t> orientation(1.0, 0.0, 0.0);
   const auto w = std::acos(up * orientation);
   math::Quaternion<real_t> q(up % orientation, w);
   auto squirmer = pe::createSquirmer(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, position, R, v0,
                                      beta, myMat);
   if (squirmer)
      squirmer->setOrientation(q);

   syncCall();

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::TRT::constructWithMagicNumber(omega));

   // add PDF field ( uInit = <0,0,0>, rhoInit = 1 )
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage<LatticeModel_T>(blocks, "pdf field (fzyx)", latticeModel,
                                                                      FieldGhostLayers);

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>(blocks, "flag field");

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<PeBodyField_T>(blocks, "body field", nullptr, field::fzyx);

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData<BoundaryHandling_T>(
         MyBoundaryHandling(flagFieldID, pdfFieldID, bodyFieldID), "boundary handling");

   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   SweepTimeloop timeloop(blocks->getBlockStorage(), T);

   // sweep for updating the pe body mapping into the LBM simulation
   timeloop.add()
         << Sweep(pe_coupling::BodyMapping<LatticeModel_T, BoundaryHandling_T, pe_coupling::NaNDestroyer<LatticeModel_T>, true>(blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID,
                                                               MO_BB_Flag, FormerMO_BB_Flag, pe_coupling::selectRegularBodies),
                  "Body Mapping");

   // sweep for restoring PDFs in cells previously occupied by pe bodies
   typedef pe_coupling::EquilibriumReconstructor<LatticeModel_T, BoundaryHandling_T> Reconstructor_T;
   Reconstructor_T reconstructor(blocks, boundaryHandlingID, bodyFieldID);
   timeloop.add()
         << Sweep(pe_coupling::PDFReconstruction<LatticeModel_T, BoundaryHandling_T, Reconstructor_T, true>(blocks,
                                                                                                      pdfFieldID,
                                                                                                      boundaryHandlingID,
                                                                                                      bodyStorageID,
                                                                                                      globalBodyStorage,
                                                                                                      bodyFieldID,
                                                                                                      reconstructor,
                                                                                                      FormerMO_BB_Flag,
                                                                                                      Fluid_Flag),
                  "PDF Restore");

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   std::function<void()> commFunction;

   blockforest::communication::UniformBufferedScheme<stencil::D3Q27> scheme(blocks);
   scheme.addPackInfo(make_shared<field::communication::PackInfo<PdfField_T> >(pdfFieldID));
   commFunction = scheme;

   // uses standard bounce back boundary conditions
   pe_coupling::mapMovingBodies<BoundaryHandling_T>(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID,
                                                    MO_BB_Flag, pe_coupling::selectRegularBodies);

   auto sweep = lbm::makeCellwiseSweep<LatticeModel_T, FlagField_T>(pdfFieldID, flagFieldID, Fluid_Flag);

   // collision sweep
   timeloop.add() << Sweep(lbm::makeCollideSweep(sweep), "cell-wise LB sweep (collide)");

   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
   timeloop.add() << BeforeFunction(commFunction, "LBM Communication")
                  << Sweep(BoundaryHandling_T::getBlockSweep(boundaryHandlingID), "Boundary Handling");

   // streaming
   timeloop.add() << Sweep(lbm::makeStreamSweep(sweep), "cell-wise LB sweep (stream)");

   // add pe timesteps
   timeloop.addFuncAfterTimeStep(pe_coupling::TimeStep(blocks, bodyStorageID, cr, syncCall, dt), "pe Time Step");
   timeloop.addFuncAfterTimeStep(ResetPosition(blocks, bodyStorageID, dt), "Reset pe positions");

   timeloop.addFuncAfterTimeStep(RemainingTimeLogger(timeloop.getNrOfTimeSteps()), "Remaining Time Logger");

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WALBERLA_LOG_INFO_ON_ROOT("Starting test...");
   WcTimingPool timeloopTiming;
   timeloop.run(timeloopTiming);
   timeloopTiming.logResultOnRoot();

   if (writevtk) {
      WALBERLA_LOG_INFO_ON_ROOT("Writing out VTK file...");
      auto pdfFieldVTKWriter = createFluidFieldVTKWriter(blocks, pdfFieldID, flagFieldID, "fluid_field");
      auto wf = vtk::writeFiles(pdfFieldVTKWriter);
      wf();
      auto flagFieldVTK = vtk::createVTKOutput_BlockData(blocks, "flag_field", 1, 0);
      flagFieldVTK->addCellDataWriter(make_shared<field::VTKWriter<FlagField_T> >(flagFieldID, "FlagField"));
      vtk::writeFiles(flagFieldVTK)();
   }

   if (squirmer) // this object only exists on the node that contains position
   {
      WALBERLA_ASSERT_FLOAT_EQUAL(squirmer->getSquirmerVelocity(), v0);
      WALBERLA_ASSERT_FLOAT_EQUAL(squirmer->getSquirmerBeta(), beta);
      WALBERLA_ASSERT_FLOAT_EQUAL(squirmer->getRadius(), R);
      WALBERLA_ASSERT_FLOAT_EQUAL(squirmer->getPosition(), position);
      WALBERLA_ASSERT_FLOAT_EQUAL(squirmer->getQuaternion(), q);
   }

   WALBERLA_LOG_INFO_ON_ROOT("Checking result...");
   auto b1 = real_c(1.5) * v0;
   auto b2 = beta * b1;
   auto e = q.rotate(up).getNormalized();
   auto radius = R;
   const auto& squirmer_pos = position;

   real_t abs_tolerance = real_c(0.0026);
   real_t rel_tolerance = real_c(0.10);

   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
      IBlock &block = *blockIt;
      auto src = block.getData<PdfField_T>(pdfFieldID);

      auto flagField = block.getData<FlagField_T>(flagFieldID);
      auto flag = flagField->getFlag(Fluid_Flag);

      WALBERLA_FOR_ALL_CELLS_XYZ(src, {
         Vector3<real_t> pos;
         blocks->getBlockLocalCellCenter(block, Cell(x, y, z), pos[0], pos[1], pos[2]);

         if (flagField->isFlagSet(x, y, z, flag)) //check for fluid/non-boundary
         {
            auto r = pos - squirmer_pos;

            auto v_ana = analytical_velocity(r, b1, b2, e, radius);

            // Superposition with the nearest neighbours (faces)
            for (auto s = stencil::D3Q7::beginNoCenter(); s != stencil::D3Q7::end(); ++s) {
               Vector3<real_t> imidpoint(squirmer_pos);
               imidpoint[0] += real_c(s.cx() * int_c(L));
               imidpoint[1] += real_c(s.cy() * int_c(L));
               imidpoint[2] += real_c(s.cz() * int_c(L));
               v_ana += analytical_velocity(pos - imidpoint, b1, b2, e, radius);
            }

            if (!shortrun && r.length() > 2 *
                                          radius) // Exclude the volume around the squirmer, as the analytical solution is a far field solution only
            {
               auto v_data = src->getVelocity(x, y, z);
               auto diff = v_data - v_ana;

               for (uint_t d = 0; d < 3; ++d) {
                  if (std::abs(diff[d] / v_ana[d]) > rel_tolerance && std::abs(diff[d]) > abs_tolerance) {
                     WALBERLA_LOG_DEVEL(
                           "Difference too large in " << Cell(x, y, z) << " (" << r.length() << "). Expected " << v_ana
                                                      << ", got " << v_data << ", relative "
                                                      << std::abs(diff[d] / v_ana[d]) << ", absolute "
                                                      << std::abs(diff[d]));
                  }
               }
            }

            if (writevtk) {
               // store the value into the PDF field so we can write it to VTK later
               src->setToEquilibrium(x, y, z, v_ana);
            }
         }
      });
   }

   if (writevtk) {
      WALBERLA_LOG_INFO_ON_ROOT("Writing out reference VTK file...");
      auto pdfFieldVTKWriter = createFluidFieldVTKWriter(blocks, pdfFieldID, flagFieldID, "fluid_field_ref");
      auto wf = vtk::writeFiles(pdfFieldVTKWriter);
      wf();
   }

   WALBERLA_LOG_INFO_ON_ROOT("Completed test.");

   return 0;
}

} // namespace squirmer

int main( int argc, char **argv ){
   squirmer::main(argc, argv);
}
