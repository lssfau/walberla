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
//! \file CombinedResolvedUnresolved.cpp
//! \author Matthias König <matthias.k.koenig@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/grid_generator/SCIterator.h"
#include "core/logging/all.h"
#include "core/math/all.h"
#include "core/timing/RemainingTimeLogger.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/StabilityChecker.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/all.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"
#include "lbm/vtk/all.h"

#include "pe/Types.h"
#include "pe/basic.h"
#include "pe/cr/ICR.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "pe_coupling/discrete_particle_methods/all.h"
#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"

#include <functional>

namespace combined_resolved_unresolved
{
///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

//////////////
// TYPEDEFS //
//////////////

// PDF field, flag field & body field

using BodyField_T = GhostLayerField<pe::BodyID, 1>;
using TensorField_T = GhostLayerField<Matrix3<real_t>, 1>;
using Vec3Field_T = GhostLayerField<Vector3<real_t>, 1>;
using ScalarField_T = GhostLayerField<real_t, 1>;
using ForceModel_T = lbm::force_model::GuoField< Vec3Field_T >;

using LatticeModel_T = lbm::D3Q19<lbm::collision_model::SRTField<ScalarField_T>, false, ForceModel_T>;
using Stencil_T  = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField< LatticeModel_T >;

using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;

const uint_t FieldGhostLayers = 1;

// boundary handling
using NoSlip_T = lbm::NoSlip<LatticeModel_T, flag_t>;

using MO_T = pe_coupling::CurvedLinear<LatticeModel_T, FlagField_T>;

using BoundaryHandling_T = BoundaryHandling<FlagField_T, Stencil_T, NoSlip_T, MO_T>;

using BodyTypeTuple = std::tuple<pe::Sphere, pe::Plane>;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag("fluid");
const FlagUID NoSlip_Flag("no slip");
const FlagUID MO_Flag("moving obstacle");
const FlagUID FormerMO_Flag("former moving obstacle");

// interpolation (for PDFs, velocity and/or solid volume fraction)
enum Interpolation { INearestNeighbor, IKernel };

// force distribution
enum Distribution { DNearestNeighbor, DKernel };

// drag correlation
enum DragCorrelation { ErgunWenYu, Tang, Stokes, Felice, Tenneti, NoDrag };

// lift correlation
enum LiftCorrelation { NoLift, Saffman };

// added mass correlation
enum AddedMassCorrelation { NoAM, Finn };

// effective viscosity
enum EffectiveViscosity { None, Rescaled, Eilers };

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////
class MyBoundaryHandling
{
 public:
   MyBoundaryHandling(const BlockDataID& flagFieldID, const BlockDataID& pdfFieldID, const BlockDataID& bodyFieldID)
      : flagFieldID_(flagFieldID), pdfFieldID_(pdfFieldID), bodyFieldID_(bodyFieldID)
   {}

   BoundaryHandling_T* operator()(IBlock* const block, const StructuredBlockStorage* const storage) const;

 private:
   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID bodyFieldID_;

}; // class MyBoundaryHandling

BoundaryHandling_T* MyBoundaryHandling::operator()(IBlock* const block,
                                                   const StructuredBlockStorage* const storage) const
{
   WALBERLA_ASSERT_NOT_NULLPTR(block);
   WALBERLA_ASSERT_NOT_NULLPTR(storage);

   FlagField_T* flagField = block->getData< FlagField_T >(flagFieldID_);
   PdfField_T* pdfField   = block->getData< PdfField_T >(pdfFieldID_);
   BodyField_T* bodyField = block->getData< BodyField_T >(bodyFieldID_);

   const auto fluid =
      flagField->flagExists(Fluid_Flag) ? flagField->getFlag(Fluid_Flag) : flagField->registerFlag(Fluid_Flag);

   BoundaryHandling_T* handling = new BoundaryHandling_T(
      "moving obstacle boundary handling", flagField, fluid, NoSlip_T("NoSlip", NoSlip_Flag, pdfField),
      MO_T("MO", MO_Flag, pdfField, flagField, bodyField, fluid, *storage, *block));

   // boundary conditions (no-slip) are set by mapping the planes into the domain

   handling->fillWithDomain(FieldGhostLayers);

   return handling;
}

// VTK Output
shared_ptr< vtk::VTKOutput > createFluidFieldVTKWriter(shared_ptr< StructuredBlockForest >& blocks,
                                                       const BlockDataID& pdfFieldID, const BlockDataID& flagFieldID,
                                                       uint_t writeFrequency, const std::string& vtkBaseFolder)
{
   auto pdfFieldVTKWriter =
      vtk::createVTKOutput_BlockData(blocks, "fluid_field", writeFrequency, uint_t(1), false, vtkBaseFolder);

   blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync(blocks);
   pdfGhostLayerSync.addPackInfo(make_shared< field::communication::PackInfo< PdfField_T > >(pdfFieldID));
   pdfFieldVTKWriter->addBeforeFunction(pdfGhostLayerSync);

   field::FlagFieldCellFilter< FlagField_T > fluidFilter(flagFieldID);
   fluidFilter.addFlag(Fluid_Flag);
   pdfFieldVTKWriter->addCellInclusionFilter(fluidFilter);

   auto velocityWriter =
      make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >(pdfFieldID, "Velocity (Lattice)");
   auto densityWriter = make_shared< lbm::DensityVTKWriter< LatticeModel_T, float > >(pdfFieldID, "Density (Lattice)");
   pdfFieldVTKWriter->addCellDataWriter(velocityWriter);
   pdfFieldVTKWriter->addCellDataWriter(densityWriter);

   return pdfFieldVTKWriter;
}

void createSphereLattice(StructuredBlockForest& forest, pe::BodyStorage& globalBodyStorage,
                         const BlockDataID& bodyStorageID, const AABB& generationDomain, real_t diameter,
                         real_t solidVolumeFraction, const pe::MaterialID& material, real_t initialZVelocity)
{
   real_t sphereVolume        = math::pi * diameter * diameter * diameter / real_t(6);
   real_t numSpheresDesired   = solidVolumeFraction * generationDomain.volume() / sphereVolume;
   uint_t spheresPerDirection = uint_c(std::cbrt(numSpheresDesired));

   real_t spacing = generationDomain.xSize() / real_c(spheresPerDirection);

   WALBERLA_ASSERT(spacing >= diameter);

   Vector3< real_t > generationOrigin(generationDomain.xMin() + spacing * real_t(0.5),
                                      generationDomain.yMin() + spacing * real_t(0.5),
                                      generationDomain.zMin() + spacing * real_t(0.5));

   uint_t numSpheres(0);

   for (auto it = grid_generator::SCIterator(generationDomain, generationOrigin, spacing);
        it != grid_generator::SCIterator(); ++it)
   {
      pe::SphereID sp = pe::createSphere(globalBodyStorage, forest.getBlockStorage(), bodyStorageID, 0, *it,
                                         diameter * real_t(0.5), material);

      if (sp != nullptr)
      {
         sp->setLinearVel(Vector3< real_t >(real_t(0), real_t(0), initialZVelocity));
         ++numSpheres;
      }
   }

   WALBERLA_MPI_SECTION() { mpi::allReduceInplace(numSpheres, mpi::SUM); }

   WALBERLA_LOG_INFO_ON_ROOT("Created " << numSpheres << " spheres of diameter " << diameter << " in "
                                        << generationDomain << " with a center-to-center spacing of " << spacing
                                        << " ( " << real_t(100) * spacing / diameter << "% of diameter )");
}

pe::MaterialID createSphereMaterial(const std::string& name, real_t diameter, real_t densityRatio)
{
   const real_t restitutionCoeff = real_t(0.88);
   const real_t frictionCoeff    = real_t(0.25);

   real_t sphereVolume        = diameter * diameter * diameter * math::pi / real_t(6);
   const real_t particleMass  = densityRatio * sphereVolume;
   const real_t Mij           = particleMass * particleMass / (real_t(2) * particleMass);
   const real_t lnDryResCoeff = std::log(restitutionCoeff);
   const real_t collisionTime = real_t(10.0);
   const real_t stiffnessCoeff =
      math::pi * math::pi * Mij /
      (collisionTime * collisionTime *
       (real_t(1) - lnDryResCoeff * lnDryResCoeff / (math::pi * math::pi + lnDryResCoeff * lnDryResCoeff)));
   const real_t dampingCoeff =
      -real_t(2) * std::sqrt(Mij * stiffnessCoeff) *
      (std::log(restitutionCoeff) /
       std::sqrt(math::pi * math::pi + (std::log(restitutionCoeff) * std::log(restitutionCoeff))));

   return pe::createMaterial(name, densityRatio, restitutionCoeff, frictionCoeff, frictionCoeff, real_t(0), real_t(200),
                             stiffnessCoeff, dampingCoeff, dampingCoeff);
}

bool selectDPMBodies(pe::BodyID bodyID)
{
   pe::SphereID sphere = static_cast< pe::SphereID >(bodyID);
   real_t radius       = sphere->getRadius();
   return radius <= real_t(0.5);
}

bool selectMEMBodies(pe::BodyID bodyID) { return !selectDPMBodies(bodyID); }

void resetSphereVelocities(const shared_ptr< StructuredBlockForest >& blocks, const BlockDataID& bodyStorageID,
                           Vector3< real_t > vel)
{
   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      for (auto bodyIt = pe::BodyIterator::begin< pe::Sphere >(*blockIt, bodyStorageID);
           bodyIt != pe::BodyIterator::end< pe::Sphere >(); ++bodyIt)
      {
         bodyIt->setAngularVel(real_t(0), real_t(0), real_t(0));
         bodyIt->setLinearVel(vel);
      }
   }
}

class DummySweep
{
 public:
   DummySweep() = default;

   void operator()(IBlock* const /*block*/) {}
};

void emptyFunction() {}

//*******************************************************************************************************************
/*!\brief Simualtion of a strongly heterogeneous sized particulate flow system using combined resolved and unresolved
 * methods.
 *
 * For the coupling of resolved particles the Momentum Exchange Method (MEM) is used, whereas for the
 * unresolved particles the Discrete Particle Method by Rettinger, Ruede - "A Coupled Lattice Boltzmann Method and
 * Discrete Element Method for Discrete Particle Simulations of Particulate Flows" is used.
 *
 * For the default setup of the showcase a 100 × 100 × 100mm box is filled with a viscous fluid.
 * The domain is resolved by a 32 × 32 × 32 lattice. One large MEM particle of diameter Dp_MEM = 35mm
 * - which equates to a resolution of 12 lattice cells - is dropped from a gap height of 37.5mm. There
 * is a layer 4096 DPM spheres of diameter Dp_DPM = 1mm from height 12.5mm to 25mm which
 * generates an average solid volume fraction of svf = 0.017 inside the layer. DPM particles thus have
 * a diameter of Dp_DPM = 0.32 cells in lattice units. For the evaluation of forces and fluid properties
 * by the DPM part of the algorithm 10 interaction subcycles are used.
 * This setup can be customized in the customization section.
 *
 * The algorithm, as well as the setup and the outcome are described in detail in
 * Koenig - "Combining fully resolved and unresolved coupling methods for strongly heterogeneous sized particulate
 * flow simulations" to be published
 *
 */
//*******************************************************************************************************************

int main(int argc, char** argv)
{
   debug::enterTestMode();

   mpi::Environment env(argc, argv);

   ///////////////////
   // Customization //
   ///////////////////

   // simulation control
   uint_t vtkIOFreq           = 0;
   std::string baseFolder_vtk = "vtk_out_Combined";

   // physical setup
   uint_t fluidType = 1;

   uint_t timesteps            = uint_t(5000);
   uint_t interactionSubcycles = uint_t(10);

   Interpolation interpol             = Interpolation::IKernel;
   Distribution dist                  = Distribution::DKernel;
   DragCorrelation dragCorr           = DragCorrelation::Tenneti;
   LiftCorrelation liftCorr           = LiftCorrelation::Saffman;
   AddedMassCorrelation addedMassCorr = AddedMassCorrelation::Finn;
   EffectiveViscosity effVisc         = EffectiveViscosity::Eilers;

   // numerical parameters
   uint_t numberOfCellsInHorizontalDirection = uint_t(32);

   for (int i = 1; i < argc; ++i)
   {
      if (std::strcmp(argv[i], "--vtkIOFreq") == 0)
      {
         vtkIOFreq = uint_c(std::atof(argv[++i]));
         continue;
      }
      if (std::strcmp(argv[i], "--fluidType") == 0)
      {
         fluidType = uint_c(std::atof(argv[++i]));
         continue;
      }
      if (std::strcmp(argv[i], "--resolution") == 0)
      {
         numberOfCellsInHorizontalDirection = uint_c(std::atof(argv[++i]));
         continue;
      }
      if (std::strcmp(argv[i], "--baseFolder_vtk") == 0)
      {
         baseFolder_vtk = argv[++i];
         continue;
      }
      if (std::strcmp(argv[i], "--timesteps") == 0)
      {
         timesteps = uint_c(std::atof(argv[++i]));
         continue;
      }
      if (std::strcmp(argv[i], "--interactionSubcycles") == 0)
      {
         interactionSubcycles = uint_c(std::atof(argv[++i]));
         continue;
      }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   //////////////////////////////////////
   // SIMULATION PROPERTIES in SI units//
   //////////////////////////////////////

   real_t densityFluid_SI;
   real_t dynamicViscosityFluid_SI;
   real_t expectedSettlingVelocity_SI;
   switch (fluidType)
   {
   case 1:
      // Re_p around 1.5
      densityFluid_SI             = real_t(970);
      dynamicViscosityFluid_SI    = real_t(373e-3);
      expectedSettlingVelocity_SI = real_t(0.035986);
      break;
   case 2:
      // Re_p around 4.1
      densityFluid_SI             = real_t(965);
      dynamicViscosityFluid_SI    = real_t(212e-3);
      expectedSettlingVelocity_SI = real_t(0.05718);
      break;
   case 3:
      // Re_p around 11.6
      densityFluid_SI             = real_t(962);
      dynamicViscosityFluid_SI    = real_t(113e-3);
      expectedSettlingVelocity_SI = real_t(0.087269);
      break;
   case 4:
      // Re_p around 31.9
      densityFluid_SI             = real_t(960);
      dynamicViscosityFluid_SI    = real_t(58e-3);
      expectedSettlingVelocity_SI = real_t(0.12224);
      break;
   default:
      WALBERLA_ABORT("Only four different fluids are supported! Choose type between 1 and 4.");
   }
   const real_t kinematicViscosityFluid_SI = dynamicViscosityFluid_SI / densityFluid_SI;

   const real_t diameter_MEM_SI      = real_t(35e-3);
   const real_t densitySphere_MEM_SI = real_t(1120);
   const real_t densityRatio_MEM     = densitySphere_MEM_SI / densityFluid_SI;

   const real_t diameter_DPM_SI      = real_t(1e-3);
   const real_t densitySphere_DPM_SI = real_t(1120);
   const real_t densityRatio_DPM     = densitySphere_DPM_SI / densityFluid_SI;

   const real_t gravitationalAcceleration_SI = real_t(9.81);
   Vector3< real_t > domainSize_SI(real_t(100e-3), real_t(100e-3), real_t(100e-3));
   // shift starting gap a bit upwards to match the reported (plotted) values
   const real_t startingGapSize_SI = real_t(20e-3) + real_t(0.5) * diameter_MEM_SI;

   WALBERLA_LOG_INFO_ON_ROOT("Setup (in SI units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize_SI);
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter_MEM = " << diameter_MEM_SI << ", density = " << densitySphere_MEM_SI
                                                          << ", starting gap size = " << startingGapSize_SI);
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: density = " << densityFluid_SI << ", dyn. visc = " << dynamicViscosityFluid_SI
                                                    << ", kin. visc = " << kinematicViscosityFluid_SI);
   WALBERLA_LOG_INFO_ON_ROOT(" - expected settling velocity = "
                             << expectedSettlingVelocity_SI << " --> Re_p = "
                             << expectedSettlingVelocity_SI * diameter_MEM_SI / kinematicViscosityFluid_SI);

   //////////////////////////
   // NUMERICAL PARAMETERS //
   //////////////////////////

   const real_t dx_SI = domainSize_SI[0] / real_c(numberOfCellsInHorizontalDirection);
   const Vector3< uint_t > domainSize(uint_c(floor(domainSize_SI[0] / dx_SI + real_t(0.5))),
                                      uint_c(floor(domainSize_SI[1] / dx_SI + real_t(0.5))),
                                      uint_c(floor(domainSize_SI[2] / dx_SI + real_t(0.5))));
   const real_t diameter_MEM     = diameter_MEM_SI / dx_SI;
   const real_t sphereVolume_MEM = math::pi / real_t(6) * diameter_MEM * diameter_MEM * diameter_MEM;

   const real_t diameter_DPM = diameter_DPM_SI / dx_SI;

   const real_t sphereVolume_DPM = math::pi / real_t(6) * diameter_DPM * diameter_DPM * diameter_DPM;

   const real_t expectedSettlingVelocity = real_t(0.01);
   const real_t dt_SI                    = expectedSettlingVelocity / expectedSettlingVelocity_SI * dx_SI;

   const real_t viscosity      = kinematicViscosityFluid_SI * dt_SI / (dx_SI * dx_SI);
   const real_t relaxationTime = real_t(1) / lbm::collision_model::omegaFromViscosity(viscosity);

   const real_t gravitationalAcceleration = gravitationalAcceleration_SI * dt_SI * dt_SI / dx_SI;

   const real_t dx                                     = real_t(1);
   const real_t dt                                     = real_t(1);
   const real_t dtInteractionSubCycle                  = real_t(dt) / real_t(interactionSubcycles);
   const real_t dtBodyVelocityTimeDerivativeEvaluation = dtInteractionSubCycle;

   const Vector3< real_t > initialFluidVelocity(real_t(0));

   const uint_t numPeSubCycles = uint_t(1);

   WALBERLA_LOG_INFO_ON_ROOT(" - dx_SI = " << dx_SI << ", dt_SI = " << dt_SI);
   WALBERLA_LOG_INFO_ON_ROOT("Setup (in simulation, i.e. lattice, units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter_MEM = " << diameter_MEM);
   WALBERLA_LOG_INFO_ON_ROOT(" - relaxation time (tau) = " << relaxationTime << ", kin. visc = " << viscosity);
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational acceleration = " << gravitationalAcceleration);
   WALBERLA_LOG_INFO_ON_ROOT(" - expected settling velocity = " << expectedSettlingVelocity << " --> Re_p = "
                                                                << expectedSettlingVelocity * diameter_MEM / viscosity);

   if (vtkIOFreq > 0)
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files to folder \"" << baseFolder_vtk << "\" with frequency "
                                                                    << vtkIOFreq);
   }

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   Vector3< uint_t > numberOfBlocksPerDirection(uint_t(2), uint_t(2), uint_t(2));
   Vector3< uint_t > cellsPerBlockPerDirection(domainSize[0] / numberOfBlocksPerDirection[0],
                                               domainSize[1] / numberOfBlocksPerDirection[1],
                                               domainSize[2] / numberOfBlocksPerDirection[2]);
   for (uint_t i = 0; i < 3; ++i)
   {
      WALBERLA_CHECK_EQUAL(cellsPerBlockPerDirection[i] * numberOfBlocksPerDirection[i], domainSize[i],
                           "Unmatching domain decomposition in direction " << i << "!");
   }

   auto blocks = blockforest::createUniformBlockGrid(numberOfBlocksPerDirection[0], numberOfBlocksPerDirection[1],
                                                     numberOfBlocksPerDirection[2], cellsPerBlockPerDirection[0],
                                                     cellsPerBlockPerDirection[1], cellsPerBlockPerDirection[2], dx, 0,
                                                     false, false, false, false, false, // periodicity
                                                     false);

   WALBERLA_LOG_INFO_ON_ROOT("Domain decomposition:");
   WALBERLA_LOG_INFO_ON_ROOT(" - blocks per direction = " << numberOfBlocksPerDirection);
   WALBERLA_LOG_INFO_ON_ROOT(" - cells per block = " << cellsPerBlockPerDirection);

   // write domain decomposition to file
   if (vtkIOFreq > 0) { vtk::writeDomainDecomposition(blocks, "initial_domain_decomposition", baseFolder_vtk); }

   /////////////////
   // PE COUPLING //
   /////////////////

   // set up pe functionality
   shared_ptr< pe::BodyStorage > globalBodyStorage = make_shared< pe::BodyStorage >();
   pe::SetBodyTypeIDs< BodyTypeTuple >::execute();

   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling< BodyTypeTuple >(), "pe Body Storage");
   auto ccdID = blocks->addBlockData(pe::ccd::createHashGridsDataHandling(globalBodyStorage, bodyStorageID), "CCD");
   auto fcdID = blocks->addBlockData(
      pe::fcd::createGenericFCDDataHandling< BodyTypeTuple, pe::fcd::AnalyticCollideFunctor >(), "FCD");

   // set up collision response, here DEM solver
   pe::cr::DEM cr(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID, nullptr);

   // set up synchronization procedure
   const real_t overlap = real_t(1.5) * dx;
   std::function< void(void) > syncCall =
      std::bind(pe::syncShadowOwners< BodyTypeTuple >, std::ref(blocks->getBlockForest()), bodyStorageID,
                static_cast< WcTimingTree* >(nullptr), overlap, false);

   // create pe bodies

   // bounding planes (global)
   const auto planeMaterial = pe::createMaterial("planeMaterial", 1000_r, 0.8_r, 0.1_r, 0.05_r, 0.2_r, 80_r, 100_r, 10_r, 11_r);
   pe::createPlane(*globalBodyStorage, 0, Vector3< real_t >(1, 0, 0), Vector3< real_t >(0, 0, 0), planeMaterial);
   pe::createPlane(*globalBodyStorage, 0, Vector3< real_t >(-1, 0, 0), Vector3< real_t >(real_c(domainSize[0]), 0, 0),
                   planeMaterial);
   pe::createPlane(*globalBodyStorage, 0, Vector3< real_t >(0, 1, 0), Vector3< real_t >(0, 0, 0), planeMaterial);
   pe::createPlane(*globalBodyStorage, 0, Vector3< real_t >(0, -1, 0), Vector3< real_t >(0, real_c(domainSize[1]), 0),
                   planeMaterial);
   pe::createPlane(*globalBodyStorage, 0, Vector3< real_t >(0, 0, 1), Vector3< real_t >(0, 0, 0), planeMaterial);
   pe::createPlane(*globalBodyStorage, 0, Vector3< real_t >(0, 0, -1), Vector3< real_t >(0, 0, real_c(domainSize[2])),
                   planeMaterial);

   // add the MEM sphere
   const auto peMaterial_MEM = createSphereMaterial("peMaterial_MEM", diameter_MEM, densityRatio_MEM);
   Vector3< real_t > initialPosition(real_t(0.5) * real_c(domainSize[0]), real_t(0.5) * real_c(domainSize[1]),
                                     startingGapSize_SI / dx_SI + real_t(0.5) * diameter_MEM);
   pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, initialPosition,
                    real_t(0.5) * diameter_MEM, peMaterial_MEM);

   // add DPM spheres
   const auto peMaterial_DPM = createSphereMaterial("peMaterial_DPM", diameter_DPM, densityRatio_DPM);
   createSphereLattice(*blocks, *globalBodyStorage, bodyStorageID,
                       AABB(real_t(0), real_t(0), real_t(1.0 * real_t(domainSize[2]) / 8.0_r), real_t(domainSize[0]),
                            real_t(domainSize[1]), real_t(2.0 * real_t(domainSize[2]) / 8.0_r)),
                       diameter_DPM, real_t(0.15), peMaterial_DPM, real_t(0));

   syncCall();

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create force field
   BlockDataID forceFieldID = field::addToStorage< Vec3Field_T >(blocks, "force field", Vector3< real_t >(real_t(0)),
                                                                 field::zyxf, FieldGhostLayers);

   // create omega field
   BlockDataID omegaFieldID =
      field::addToStorage< ScalarField_T >(blocks, "omega field", real_t(0), field::zyxf, FieldGhostLayers);

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T(omegaFieldID, ForceModel_T(forceFieldID));

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >(
      blocks, "pdf field (zyxf)", latticeModel, Vector3< real_t >(real_t(0)), real_t(1), uint_t(1), field::zyxf);
   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

   // add body field
   BlockDataID bodyFieldID = field::addToStorage< BodyField_T >(blocks, "body field", nullptr, field::zyxf);

   // add boundary handling
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
      MyBoundaryHandling(flagFieldID, pdfFieldID, bodyFieldID), "boundary handling");

   // field to store fluid velolcity
   BlockDataID velocityFieldID =
      field::addToStorage< Vec3Field_T >(blocks, "velocity field", initialFluidVelocity, field::zyxf, FieldGhostLayers);

   BlockDataID oldVelocityFieldID = field::addToStorage< Vec3Field_T >(
      blocks, "old velocity field", initialFluidVelocity, field::zyxf, FieldGhostLayers);
   BlockDataID swappedOldVelocityFieldID = field::addToStorage< Vec3Field_T >(
      blocks, "swapped old velocity field", initialFluidVelocity, field::zyxf, FieldGhostLayers);

   // field to store curl of fluid velocity
   BlockDataID velocityCurlFieldID = field::addToStorage< Vec3Field_T >(
      blocks, "velocity curl field", Vector3< real_t >(real_c(0)), field::zyxf, FieldGhostLayers);

   // field to store velocity gradient
   BlockDataID velocityGradientFieldID = field::addToStorage< TensorField_T >(
      blocks, "velocity gradient field", Matrix3< real_t >(real_c(0)), field::zyxf, FieldGhostLayers);

   // field to store time derivative of fluid velocity
   BlockDataID timeDerivativeVelocityFieldID = field::addToStorage< Vec3Field_T >(
      blocks, "time derivative velocity field", Vector3< real_t >(real_c(0)), field::zyxf, FieldGhostLayers);

   // create solid volume fraction field
   BlockDataID svfFieldID =
      field::addToStorage< ScalarField_T >(blocks, "svf field", real_t(0), field::zyxf, FieldGhostLayers);

   // create pressure field
   BlockDataID pressureFieldID =
      field::addToStorage< ScalarField_T >(blocks, "pressure field", real_t(0), field::zyxf, FieldGhostLayers);

   // field to store pressure gradient
   BlockDataID pressureGradientFieldID = field::addToStorage< Vec3Field_T >(
      blocks, "pressure gradient field", Vector3< real_t >(real_c(0)), field::zyxf, FieldGhostLayers);

   BlockDataID dragForceFieldID = field::addToStorage< Vec3Field_T >(
      blocks, "drag force field", Vector3< real_t >(real_t(0)), field::zyxf, FieldGhostLayers);

   BlockDataID amForceFieldID = field::addToStorage< Vec3Field_T >(
      blocks, "am force field", Vector3< real_t >(real_t(0)), field::zyxf, FieldGhostLayers);
   BlockDataID liftForceFieldID = field::addToStorage< Vec3Field_T >(
      blocks, "lift force field", Vector3< real_t >(real_t(0)), field::zyxf, FieldGhostLayers);

   // map planes into the LBM simulation -> act as no-slip boundaries
   pe_coupling::mapBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage,
                                                NoSlip_Flag, pe_coupling::selectGlobalBodies);

   // map pe bodies into the LBM simulation
   pe_coupling::mapMovingBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage,
                                                      bodyFieldID, MO_Flag, selectMEMBodies);

   // communication for synchronizing the interaction force field
   pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication< Vec3Field_T > dragForceComm(
      blocks, dragForceFieldID);

   pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication< Vec3Field_T > amForceComm(
      blocks, amForceFieldID);
   pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication< Vec3Field_T > liftForceComm(
      blocks, liftForceFieldID);

   blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > pdfScheme(blocks);
   pdfScheme.addPackInfo(make_shared< field::communication::PackInfo< PdfField_T > >(pdfFieldID));

   shared_ptr< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > > velocityCommunicationScheme =
      make_shared< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(blocks);
   velocityCommunicationScheme->addPackInfo(
      make_shared< field::communication::PackInfo< Vec3Field_T > >(velocityFieldID));

   shared_ptr< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > > oldVelocityCommunicationScheme =
      make_shared< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(blocks);
   oldVelocityCommunicationScheme->addPackInfo(
      make_shared< field::communication::PackInfo< Vec3Field_T > >(oldVelocityFieldID));

   // setup of the communication for synchronizing the solid volume fraction field between neighboring blocks which
   // takes into account the svf values inside the ghost layers
   shared_ptr< pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication< ScalarField_T > >
      svfCommunicationScheme =
         make_shared< pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication< ScalarField_T > >(
            blocks, svfFieldID);

   // setup of the communication for synchronizing the pressure field between neighboring blocks
   shared_ptr< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > > pressureCommunicationScheme =
      make_shared< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(blocks);
   pressureCommunicationScheme->addPackInfo(
      make_shared< field::communication::PackInfo< ScalarField_T > >(pressureFieldID));

   // setup of the communication for synchronizing the pressure gradient field between neighboring blocks
   shared_ptr< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >
      pressureGradientCommunicationScheme =
         make_shared< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(blocks);
   pressureGradientCommunicationScheme->addPackInfo(
      make_shared< field::communication::PackInfo< Vec3Field_T > >(pressureGradientFieldID));

   // setup of the communication for synchronizing the velocity curl field between neighboring blocks
   shared_ptr< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > > velocityCurlCommunicationScheme =
      make_shared< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(blocks);
   velocityCurlCommunicationScheme->addPackInfo(
      make_shared< field::communication::PackInfo< Vec3Field_T > >(velocityCurlFieldID));

   // communication for synchronizing the velocity gradient values
   shared_ptr< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >
      velocityGradientCommunicationScheme =
         make_shared< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(blocks);
   velocityGradientCommunicationScheme->addPackInfo(
      make_shared< field::communication::PackInfo< TensorField_T > >(velocityGradientFieldID));

   shared_ptr< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >
      dragForceFieldToForceFieldAdder =
         make_shared< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
            blocks, forceFieldID, dragForceFieldID, uint_t(1));
   shared_ptr< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >
      amForceFieldToForceFieldAdder =
         make_shared< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
            blocks, forceFieldID, amForceFieldID, uint_t(1));
   shared_ptr< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >
      liftForceFieldToForceFieldAdder =
         make_shared< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
            blocks, forceFieldID, liftForceFieldID, uint_t(1));

   /////////////////////////////////
   //                             //
   //    CORRELATION FUNCTIONS    //
   //                             //
   /////////////////////////////////

   // drag correlation function
   std::function< Vector3< real_t >(const Vector3< real_t >&, const Vector3< real_t >&, real_t, real_t, real_t,
                                    real_t) >
      dragCorrelationFunction;
   if (dragCorr == DragCorrelation::ErgunWenYu)
   { dragCorrelationFunction = pe_coupling::discrete_particle_methods::dragForceErgunWenYu; }
   else if (dragCorr == DragCorrelation::Tang)
   {
      dragCorrelationFunction = pe_coupling::discrete_particle_methods::dragForceTang;
   }
   else if (dragCorr == DragCorrelation::Stokes)
   {
      dragCorrelationFunction = pe_coupling::discrete_particle_methods::dragForceStokes;
   }
   else if (dragCorr == DragCorrelation::Felice)
   {
      dragCorrelationFunction = pe_coupling::discrete_particle_methods::dragForceFelice;
   }
   else if (dragCorr == DragCorrelation::Tenneti)
   {
      dragCorrelationFunction = pe_coupling::discrete_particle_methods::dragForceTenneti;
   }
   else if (dragCorr == DragCorrelation::NoDrag)
   {
      dragCorrelationFunction = pe_coupling::discrete_particle_methods::noDragForce;
   }
   else
   {
      WALBERLA_ABORT("Drag correlation not yet implemented!");
   }

   // lift correlation function
   std::function< Vector3< real_t >(const Vector3< real_t >&, const Vector3< real_t >&, const Vector3< real_t >&,
                                    real_t, real_t, real_t) >
      liftCorrelationFunction;
   if (liftCorr == LiftCorrelation::NoLift)
   { liftCorrelationFunction = pe_coupling::discrete_particle_methods::noLiftForce; }
   else if (liftCorr == LiftCorrelation::Saffman)
   {
      liftCorrelationFunction = pe_coupling::discrete_particle_methods::liftForceSaffman;
   }
   else
   {
      WALBERLA_ABORT("Lift correlation not yet implemented!");
   }

   // added mass correlation function
   std::function< Vector3< real_t >(const Vector3< real_t >&, const Vector3< real_t >&, real_t, real_t) >
      addedMassCorrelationFunction;
   if (addedMassCorr == AddedMassCorrelation::NoAM)
   { addedMassCorrelationFunction = pe_coupling::discrete_particle_methods::noAddedMassForce; }
   else if (addedMassCorr == AddedMassCorrelation::Finn)
   {
      addedMassCorrelationFunction = pe_coupling::discrete_particle_methods::addedMassForceFinn;
   }
   else
   {
      WALBERLA_ABORT("Added mass correlation not yet implemented!");
   }

   ////////////////
   // EVALUATION //
   ////////////////

   // evaluator for bodies' velocity time derivative
   shared_ptr< pe_coupling::discrete_particle_methods::BodyVelocityTimeDerivativeEvaluator >
      bodyVelocityTimeDerivativeEvaluator =
         make_shared< pe_coupling::discrete_particle_methods::BodyVelocityTimeDerivativeEvaluator >(
            blocks, bodyStorageID, dtBodyVelocityTimeDerivativeEvaluation);
   (*bodyVelocityTimeDerivativeEvaluator)();

   // function used to evaluate the interaction force between fluid and particles
   std::function< void(void) > dragAndPressureForceEvaluationFunction;
   if (interpol == Interpolation::INearestNeighbor)
   {
      if (dist == Distribution::DNearestNeighbor)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T >(
            blocks, dragForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, velocityFieldID, svfFieldID,
            pressureGradientFieldID, dragCorrelationFunction, viscosity, selectDPMBodies);
         dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if (dist == Distribution::DKernel)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::KernelDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T >(
            blocks, dragForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, velocityFieldID, svfFieldID,
            pressureGradientFieldID, dragCorrelationFunction, viscosity, selectDPMBodies);
         dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }
   else if (interpol == Interpolation::IKernel)
   {
      if (dist == Distribution::DNearestNeighbor)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T >(
            blocks, dragForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, velocityFieldID, svfFieldID,
            pressureGradientFieldID, dragCorrelationFunction, viscosity, selectDPMBodies);
         dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if (dist == Distribution::DKernel)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::KernelDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T >(
            blocks, dragForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, velocityFieldID, svfFieldID,
            pressureGradientFieldID, dragCorrelationFunction, viscosity, selectDPMBodies);
         dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }

   // function to evaluate the lift force contribution
   std::function< void(void) > liftForceEvaluationFunction;
   if (interpol == Interpolation::INearestNeighbor)
   {
      if (dist == Distribution::DNearestNeighbor)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::LiftForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr =
            make_shared< IFE_T >(blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, velocityFieldID,
                                 velocityCurlFieldID, liftCorrelationFunction, viscosity, selectDPMBodies);
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if (dist == Distribution::DKernel)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::LiftForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::KernelDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr =
            make_shared< IFE_T >(blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, velocityFieldID,
                                 velocityCurlFieldID, liftCorrelationFunction, viscosity, selectDPMBodies);
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }
   else if (interpol == Interpolation::IKernel)
   {
      if (dist == Distribution::DNearestNeighbor)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::LiftForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr =
            make_shared< IFE_T >(blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, velocityFieldID,
                                 velocityCurlFieldID, liftCorrelationFunction, viscosity, selectDPMBodies);
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if (dist == Distribution::DKernel)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::LiftForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::KernelDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr =
            make_shared< IFE_T >(blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, velocityFieldID,
                                 velocityCurlFieldID, liftCorrelationFunction, viscosity, selectDPMBodies);
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }

   // function to evaluate the added mass contribution
   std::function< void(void) > addedMassEvaluationFunction;
   if (interpol == Interpolation::INearestNeighbor)
   {
      if (dist == Distribution::DNearestNeighbor)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::AddedMassForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T >(
            blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, timeDerivativeVelocityFieldID,
            addedMassCorrelationFunction, bodyVelocityTimeDerivativeEvaluator, selectDPMBodies);
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if (dist == Distribution::DKernel)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::AddedMassForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::KernelDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T >(
            blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, timeDerivativeVelocityFieldID,
            addedMassCorrelationFunction, bodyVelocityTimeDerivativeEvaluator, selectDPMBodies);
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }
   else if (interpol == Interpolation::IKernel)
   {
      if (dist == Distribution::DNearestNeighbor)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::AddedMassForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T >(
            blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, timeDerivativeVelocityFieldID,
            addedMassCorrelationFunction, bodyVelocityTimeDerivativeEvaluator, selectDPMBodies);
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if (dist == Distribution::DKernel)
      {
         using IFE_T = pe_coupling::discrete_particle_methods::AddedMassForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::KernelDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T >(
            blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag, timeDerivativeVelocityFieldID,
            addedMassCorrelationFunction, bodyVelocityTimeDerivativeEvaluator, selectDPMBodies);
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }

   // set up effective viscosity calculation
   std::function< real_t(real_t, real_t) > effectiveViscosityFunction;
   if (effVisc == EffectiveViscosity::None)
   { effectiveViscosityFunction = pe_coupling::discrete_particle_methods::calculateUnchangedEffectiveViscosity; }
   else if (effVisc == EffectiveViscosity::Rescaled)
   {
      effectiveViscosityFunction = pe_coupling::discrete_particle_methods::calculateRescaledEffectiveViscosity;
   }
   else if (effVisc == EffectiveViscosity::Eilers)
   {
      effectiveViscosityFunction = pe_coupling::discrete_particle_methods::calculateEilersEffectiveViscosity;
   }
   else
   {
      WALBERLA_ABORT("Effective viscosity not yet implemented!");
   }
   shared_ptr< pe_coupling::discrete_particle_methods::EffectiveViscosityFieldEvaluator > effectiveViscosityEvaluator =
      walberla::make_shared< pe_coupling::discrete_particle_methods::EffectiveViscosityFieldEvaluator >(
         omegaFieldID, svfFieldID, viscosity, effectiveViscosityFunction);

   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

   shared_ptr< pe_coupling::BodiesForceTorqueContainer > bodiesFTContainer1 =
      make_shared< pe_coupling::BodiesForceTorqueContainer >(blocks, bodyStorageID, selectMEMBodies);
   std::function< void(void) > storeForceTorqueInCont1 =
      std::bind(&pe_coupling::BodiesForceTorqueContainer::store, bodiesFTContainer1);
   shared_ptr< pe_coupling::BodiesForceTorqueContainer > bodiesFTContainer2 =
      make_shared< pe_coupling::BodiesForceTorqueContainer >(blocks, bodyStorageID, selectMEMBodies);
   std::function< void(void) > setForceTorqueOnBodiesFromCont2 =
      std::bind(&pe_coupling::BodiesForceTorqueContainer::setOnBodies, bodiesFTContainer2);

   shared_ptr< pe_coupling::BodiesForceTorqueContainer > bodiesFTContainer3 =
      make_shared< pe_coupling::BodiesForceTorqueContainer >(blocks, bodyStorageID, selectMEMBodies);
   std::function< void(void) > storeForceTorqueInCont3 =
      std::bind(&pe_coupling::BodiesForceTorqueContainer::store, bodiesFTContainer3);
   std::function< void(void) > setForceTorqueOnBodiesFromCont3 =
      std::bind(&pe_coupling::BodiesForceTorqueContainer::setOnBodies, bodiesFTContainer3);

   shared_ptr< pe_coupling::ForceTorqueOnBodiesScaler > forceScaler =
      make_shared< pe_coupling::ForceTorqueOnBodiesScaler >(blocks, bodyStorageID, real_t(1), selectMEMBodies);

   std::function< void(void) > setForceScalingFactorToHalf =
      std::bind(&pe_coupling::ForceTorqueOnBodiesScaler::resetScalingFactor, forceScaler, real_t(0.5));

   // initial svf eval
   if (dist == Distribution::DNearestNeighbor)
   {
      pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator< FlagField_T,
                                                                                 field::NearestNeighborDistributor >
         svfEvaluator(blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag, selectDPMBodies);
      for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
         svfEvaluator(&(*blockIt));
   }
   else
   {
      pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator< FlagField_T, field::KernelDistributor >
         svfEvaluator(blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag, selectDPMBodies);
      for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
         svfEvaluator(&(*blockIt));
   }

   (*svfCommunicationScheme)();

   pe_coupling::discrete_particle_methods::ForceFieldResetter forceFieldResetter(forceFieldID);
   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
      forceFieldResetter(&(*blockIt));

   // Averaging the force/torque over two time steps is said to damp oscillations of the interaction force/torque.
   // See Ladd - " Numerical simulations of particulate suspensions via a discretized Boltzmann equation. Part 1.
   // Theoretical foundation", 1994, p. 302

   timeloop.add() << Sweep(DummySweep(), "Dummy Sweep ")
                  // set force/torque from previous time step (in container2)
                  << AfterFunction(setForceTorqueOnBodiesFromCont2, "Force setting")
                  // average the force/torque by scaling it with factor 1/2 (except in first timestep, there it is 1,
                  // which it is initially)
                  << AfterFunction(SharedFunctor< pe_coupling::ForceTorqueOnBodiesScaler >(forceScaler),
                                   "Force averaging")
                  << AfterFunction(setForceScalingFactorToHalf, "Force scaling adjustment")
                  // swap containers
                  << AfterFunction(pe_coupling::BodyContainerSwapper(bodiesFTContainer1, bodiesFTContainer2),
                                   "Swap FT container")
                  << AfterFunction(storeForceTorqueInCont3, "Force Storing")
                  << AfterFunction(pe_coupling::ForceTorqueOnBodiesResetter(blocks, bodyStorageID, selectMEMBodies),
                                   "Reset MEM Forces");

   if (addedMassCorr != AddedMassCorrelation::NoAM)
   {
      timeloop.add() << Sweep(pe_coupling::discrete_particle_methods::FieldDataSwapper< Vec3Field_T >(
                                 velocityFieldID, swappedOldVelocityFieldID),
                              "Velocity Field Swap");
   }

   if (addedMassCorr != AddedMassCorrelation::NoAM)
   {
      timeloop.add()
         << Sweep(
               pe_coupling::discrete_particle_methods::GNSVelocityFieldEvaluator< LatticeModel_T, BoundaryHandling_T >(
                  oldVelocityFieldID, pdfFieldID, svfFieldID, boundaryHandlingID),
               "Old Velocity Field Evaluation")
         << AfterFunction(SharedFunctor< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(
                             oldVelocityCommunicationScheme),
                          "Old Velocity Field Communication");
      timeloop.add()
         << Sweep(pe_coupling::discrete_particle_methods::VelocityGradientFieldEvaluator< LatticeModel_T,
                                                                                          BoundaryHandling_T >(
                     velocityGradientFieldID, oldVelocityFieldID, boundaryHandlingID),
                  "Velocity Gradient Field Evaluation")
         << AfterFunction(SharedFunctor< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(
                             velocityGradientCommunicationScheme),
                          "Velocity Gradient Field Communication");
      timeloop.add() << Sweep(
         pe_coupling::discrete_particle_methods::VelocityTotalTimeDerivativeFieldEvaluator(
            timeDerivativeVelocityFieldID, oldVelocityFieldID, swappedOldVelocityFieldID, velocityGradientFieldID, dt),
         "Velocity Time Derivative Field Evaluation");
   }

   // subcycling loop begin
   for (uint_t subcycle = 1; subcycle <= interactionSubcycles; ++subcycle)
   {
      timeloop.add() << Sweep(DummySweep(), "Dummy Sweep ")
                     << AfterFunction(setForceTorqueOnBodiesFromCont3, "Force setting");

      timeloop.add() << Sweep(pe_coupling::discrete_particle_methods::ForceFieldResetter(forceFieldID),
                              "Force Field Reset");
      timeloop.add()
         << Sweep(DummySweep(), "Force Field Add")
         << AfterFunction(
               SharedFunctor< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
                  dragForceFieldToForceFieldAdder),
               "Drag Force Field To Force Field Adder")
         << AfterFunction(
               SharedFunctor< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
                  amForceFieldToForceFieldAdder),
               "AM Force Field To Force Field Adder")
         << AfterFunction(
               SharedFunctor< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
                  liftForceFieldToForceFieldAdder),
               "Lift Force Field To Force Field Adder");

      timeloop.add()
         << Sweep(
               pe_coupling::discrete_particle_methods::GNSVelocityFieldEvaluator< LatticeModel_T, BoundaryHandling_T >(
                  velocityFieldID, pdfFieldID, svfFieldID, boundaryHandlingID),
               "Velocity Field Evaluation")
         << AfterFunction(SharedFunctor< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(
                             velocityCommunicationScheme),
                          "Velocity Field Communication");

      if (liftCorr != LiftCorrelation::NoLift)
      {
         timeloop.add()
            << Sweep(pe_coupling::discrete_particle_methods::VelocityCurlFieldEvaluator< LatticeModel_T,
                                                                                         BoundaryHandling_T >(
                        velocityCurlFieldID, velocityFieldID, boundaryHandlingID),
                     "Velocity Curl Field Evaluation")
            << AfterFunction(SharedFunctor< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(
                                velocityCurlCommunicationScheme),
                             "Velocity Curl Field Communication");
      }

      timeloop.add()
         << Sweep(
               pe_coupling::discrete_particle_methods::GNSPressureFieldEvaluator< LatticeModel_T, BoundaryHandling_T >(
                  pressureFieldID, pdfFieldID, svfFieldID, boundaryHandlingID),
               "Pressure Field Evaluation")
         << AfterFunction(SharedFunctor< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(
                             pressureCommunicationScheme),
                          "Pressure Field Communication");

      timeloop.add()
         << Sweep(pe_coupling::discrete_particle_methods::PressureGradientFieldEvaluator< LatticeModel_T,
                                                                                          BoundaryHandling_T >(
                     pressureGradientFieldID, pressureFieldID, boundaryHandlingID),
                  "Pressure Gradient Field Evaluation")
         << AfterFunction(SharedFunctor< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(
                             pressureGradientCommunicationScheme),
                          "Pressure Gradient Field Communication");

      if (liftCorr != LiftCorrelation::NoLift)
      {
         // evaluate Flift
         timeloop.add() << Sweep(pe_coupling::discrete_particle_methods::ForceFieldResetter(liftForceFieldID),
                                 "Lift Force Field Reset")
                        << AfterFunction(liftForceEvaluationFunction, "Lift Force Evaluation")
                        << AfterFunction(liftForceComm, "Lift Force Field Communication");
      }

      if (addedMassCorr != AddedMassCorrelation::NoAM)
      {
         // evaluate Fam
         timeloop.add()
            << Sweep(pe_coupling::discrete_particle_methods::ForceFieldResetter(amForceFieldID), "AM Force Field Reset")
            << AfterFunction(addedMassEvaluationFunction, "Added Mass Force Evaluation")
            << AfterFunction(amForceComm, "Force Field Communication")
            << AfterFunction(
                  SharedFunctor< pe_coupling::discrete_particle_methods::BodyVelocityTimeDerivativeEvaluator >(
                     bodyVelocityTimeDerivativeEvaluator),
                  "Body Velocity Time Derivative Evaluation");
      }

      if (liftCorr != LiftCorrelation::NoLift || addedMassCorr != AddedMassCorrelation::NoAM)
      {
         timeloop.add()
            << Sweep(pe_coupling::discrete_particle_methods::ForceFieldResetter(forceFieldID), "Force Field Reset")
            << AfterFunction(
                  SharedFunctor<
                     pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
                     liftForceFieldToForceFieldAdder),
                  "Lift Force Field To Force Field Adder")
            << AfterFunction(
                  SharedFunctor<
                     pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
                     amForceFieldToForceFieldAdder),
                  "AM Force Field To Force Field Adder")
            << AfterFunction(
                  SharedFunctor<
                     pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
                     dragForceFieldToForceFieldAdder),
                  "Drag Force Field To Force Field Adder");

         timeloop.add()
            << Sweep(pe_coupling::discrete_particle_methods::GNSVelocityFieldEvaluator< LatticeModel_T,
                                                                                        BoundaryHandling_T >(
                        velocityFieldID, pdfFieldID, svfFieldID, boundaryHandlingID),
                     "Velocity Field Evaluation")
            << AfterFunction(SharedFunctor< blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > >(
                                velocityCommunicationScheme),
                             "Velocity Field Communication");
      }

      // evaluate Fdrag
      timeloop.add() << Sweep(pe_coupling::discrete_particle_methods::ForceFieldResetter(dragForceFieldID),
                              "Drag Force Field Reset")
                     << AfterFunction(dragAndPressureForceEvaluationFunction,
                                      "Fluid-Particle Interaction Force Evaluation")
                     << AfterFunction(dragForceComm, "Drag Force Field Communication");

      Vector3< real_t > gravitationalForce_MEM(0, 0, -gravitationalAcceleration * densityRatio_MEM * sphereVolume_MEM);
      Vector3< real_t > gravitationalForce_DPM(0, 0, -gravitationalAcceleration * densityRatio_DPM * sphereVolume_DPM);

      Vector3< real_t > buoyancyForce_MEM(0, 0, gravitationalAcceleration * real_t(1) * sphereVolume_MEM);
      Vector3< real_t > buoyancyForce_DPM(0, 0, gravitationalAcceleration * real_t(1) * sphereVolume_DPM);

      // ext forces on bodies
      timeloop.add()
         << Sweep(DummySweep(), "Dummy Sweep ")
         << AfterFunction(
               pe_coupling::ForceOnBodiesAdder(blocks, bodyStorageID, gravitationalForce_DPM, selectDPMBodies),
               "DPM Gravitational Force Add")

         << AfterFunction(
               pe_coupling::ForceOnBodiesAdder(blocks, bodyStorageID, gravitationalForce_MEM, selectMEMBodies),
               "MEM Gravitational Force Add")
         << AfterFunction(pe_coupling::ForceOnBodiesAdder(blocks, bodyStorageID, buoyancyForce_DPM, selectDPMBodies),
                          "DPM Buoyancy Force (due to gravity) Add")

         << AfterFunction(pe_coupling::ForceOnBodiesAdder(blocks, bodyStorageID, buoyancyForce_MEM, selectMEMBodies),
                          "MEM Buoyancy Force (due to gravity) Add")

         << AfterFunction(pe_coupling::TimeStep(blocks, bodyStorageID, cr, syncCall, dtInteractionSubCycle,
                                                numPeSubCycles, emptyFunction),
                          "Pe Time Step");

      // sweep for updating the pe body mapping into the LBM simulation
      timeloop.add() << Sweep(pe_coupling::BodyMapping< LatticeModel_T, BoundaryHandling_T >(
                                 blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID,
                                 MO_Flag, FormerMO_Flag, selectMEMBodies),
                              "Body Mapping");

      // sweep for restoring PDFs in cells previously occupied by pe bodies
      pe_coupling::SphereNormalExtrapolationDirectionFinder extrapolationFinder(blocks, bodyFieldID);
      using Reconstructor_T = pe_coupling::ExtrapolationReconstructor<LatticeModel_T, BoundaryHandling_T, pe_coupling::SphereNormalExtrapolationDirectionFinder>;
      Reconstructor_T reconstructor(blocks, boundaryHandlingID, bodyFieldID, extrapolationFinder, true);

      timeloop.add() << Sweep(pe_coupling::PDFReconstruction< LatticeModel_T, BoundaryHandling_T, Reconstructor_T >(
                                 blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID,
                                 reconstructor, FormerMO_Flag, Fluid_Flag, selectMEMBodies),
                              "PDF Restore");

      if (dist == Distribution::DNearestNeighbor)
      {
         timeloop.add()
            << Sweep(pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator<
                        FlagField_T, field::NearestNeighborDistributor >(blocks, svfFieldID, bodyStorageID, flagFieldID,
                                                                         Fluid_Flag, selectDPMBodies),
                     "Solid Volume Fraction Field Evaluation")
            << AfterFunction(
                  SharedFunctor<
                     pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication< ScalarField_T > >(
                     svfCommunicationScheme),
                  "Solid Volume Fraction Field Communication");
      }
      else
      {
         timeloop.add()
            << Sweep(
                  pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator< FlagField_T,
                                                                                             field::KernelDistributor >(
                     blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag, selectDPMBodies),
                  "Solid Volume Fraction Field Evaluation")
            << AfterFunction(
                  SharedFunctor<
                     pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication< ScalarField_T > >(
                     svfCommunicationScheme),
                  "Solid Volume Fraction Field Communication");
      }
   }

   timeloop.add() << Sweep(pe_coupling::discrete_particle_methods::ForceFieldResetter(forceFieldID),
                           "Force Field Reset");
   timeloop.add()
      << Sweep(DummySweep(), "Force Field Add")
      << AfterFunction(
            SharedFunctor< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
               dragForceFieldToForceFieldAdder),
            "Drag Force Field To Force Field Adder")
      << AfterFunction(
            SharedFunctor< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
               amForceFieldToForceFieldAdder),
            "AM Force Field To Force Field Adder")
      << AfterFunction(
            SharedFunctor< pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder >(
               liftForceFieldToForceFieldAdder),
            "Lift Force Field To Force Field Adder");

   timeloop.add() << Sweep(makeSharedSweep< pe_coupling::discrete_particle_methods::EffectiveViscosityFieldEvaluator >(
                              effectiveViscosityEvaluator),
                           "Effective Viscosity Evaluation");

   bodiesFTContainer2->store();

   // execute GNS-LBM sweep, boundary handling, PDF communication
   auto sweep = pe_coupling::discrete_particle_methods::makeGNSSweep< LatticeModel_T, FlagField_T >(
      pdfFieldID, svfFieldID, flagFieldID, Fluid_Flag);

   timeloop.add() << Sweep(lbm::makeCollideSweep(sweep), "GNS-LBM sweep (collide)");

   timeloop.add() << BeforeFunction(pdfScheme, "LBM Communication")
                  << Sweep(BoundaryHandling_T::getBlockSweep(boundaryHandlingID), "Boundary Handling");

   timeloop.add() << Sweep(lbm::makeStreamSweep(sweep), "GNS-LBM sweep (stream)");

   // store force/torque from hydrodynamic interactions in container1
   timeloop.addFuncAfterTimeStep(storeForceTorqueInCont1, "Force Storing");

   if (vtkIOFreq != uint_t(0))
   {
      shared_ptr< vtk::VTKOutput > pdfFieldVTKWriter =
         createFluidFieldVTKWriter(blocks, pdfFieldID, flagFieldID, vtkIOFreq, baseFolder_vtk);
      timeloop.addFuncAfterTimeStep(vtk::writeFiles(pdfFieldVTKWriter), "VTK (fluid field data)");

      auto bodyVtkOutput = make_shared< pe::SphereVtkOutput >(bodyStorageID, blocks->getBlockStorage());
      auto bodyVTK       = vtk::createVTKOutput_PointData(bodyVtkOutput, "bodies", vtkIOFreq, baseFolder_vtk);
      timeloop.addFuncAfterTimeStep(vtk::writeFiles(bodyVTK), "VTK (sphere data)");

      timeloop.addFuncAfterTimeStep(field::createVTKOutput< Vec3Field_T, float >(forceFieldID, *blocks, "force_field",
                                                                                 vtkIOFreq, uint_t(0), false,
                                                                                 baseFolder_vtk),
                                    "VTK (force field)");
      timeloop.addFuncAfterTimeStep(field::createVTKOutput< ScalarField_T, float >(
                                       svfFieldID, *blocks, "svf_field", vtkIOFreq, uint_t(0), false, baseFolder_vtk),
                                    "VTK (svf field)");

      timeloop.addFuncAfterTimeStep(field::createVTKOutput< ScalarField_T, float >(omegaFieldID, *blocks, "omega_field",
                                                                                   vtkIOFreq, uint_t(0), false,
                                                                                   baseFolder_vtk),
                                    "VTK (omega field)");
   }

   timeloop.addFuncAfterTimeStep(RemainingTimeLogger(timeloop.getNrOfTimeSteps()), "Remaining Time Logger");

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   // time loop
   for (uint_t i = 0; i < timesteps; ++i)
   {
      // perform a single simulation step
      timeloop.singleStep(timeloopTiming);
   }

   timeloopTiming.logResultOnRoot();

   return EXIT_SUCCESS;
}

} // namespace combined_resolved_unresolved

int main(int argc, char** argv) { combined_resolved_unresolved::main(argc, argv); }
