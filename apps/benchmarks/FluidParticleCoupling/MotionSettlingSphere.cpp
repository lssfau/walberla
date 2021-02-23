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
//! \file MotionSettlingSphere.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "boundary/all.h"

#include "blockforest/all.h"

#include "core/all.h"

#include "domain_decomposition/all.h"

#include "field/AddToStorage.h"
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

#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/CurvedLinear.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/ExtrapolationDirectionFinder.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
#include "lbm_mesapd_coupling/utility/AddForceOnParticlesKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/AddHydrodynamicInteractionKernel.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/OmegaBulkAdaption.h"
#include "lbm_mesapd_coupling/utility/InspectionProbe.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"

#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/domain/BlockForestDataHandling.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/kernel/ExplicitEuler.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/mpi/SyncGhostOwners.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/notifications/HydrodynamicForceTorqueNotification.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/BlockCellDataWriter.h"
#include "vtk/Initialization.h"
#include "vtk/VTKOutput.h"

#include <functional>
#include <memory>

#ifdef WALBERLA_BUILD_WITH_CODEGEN
#include "GeneratedLBM.h"
#endif

namespace motion_settling_sphere{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

//////////////
// TYPEDEFS //
//////////////

using namespace walberla;
using walberla::uint_t;

#ifdef WALBERLA_BUILD_WITH_CODEGEN
using LatticeModel_T = lbm::GeneratedLBM;
#else
using LatticeModel_T = lbm::D3Q19< lbm::collision_model::D3Q19MRT>;
#endif

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

using ScalarField_T = GhostLayerField< real_t, 1>;

const uint_t FieldGhostLayers = 1;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag   ( "fluid" );

const FlagUID UBB_Flag     ( "velocity bounce back" );
const FlagUID Outlet_Flag  ( "outlet" );

const FlagUID MEM_BB_Flag   ( "moving obstacle BB" );
const FlagUID MEM_CLI_Flag  ( "moving obstacle CLI" );
const FlagUID FormerMEM_Flag( "former moving obstacle" );


//////////////////////////////
// Coupling algorithm enums //
//////////////////////////////
enum MEMVariant { BB, CLI, MR };

MEMVariant to_MEMVariant( const std::string& s )
{
   if( s == "BB"  ) return MEMVariant::BB;
   if( s == "CLI" ) return MEMVariant::CLI;
   throw std::runtime_error("invalid conversion from text to MEMVariant");
}

std::string MEMVariant_to_string ( const MEMVariant& m )
{
   if( m == MEMVariant::BB  ) return "BB";
   if( m == MEMVariant::CLI ) return "CLI";
   throw std::runtime_error("invalid conversion from MEMVariant to string");
}


/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////
template <typename ParticleAccessor_T>
class MyBoundaryHandling
{
public:

   using UBB_T = lbm::SimpleUBB< LatticeModel_T, flag_t >;
   using Outlet_T = lbm::SimplePressure< LatticeModel_T, flag_t >;
   using MO_SBB_T = lbm_mesapd_coupling::SimpleBB< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
   using MO_CLI_T = lbm_mesapd_coupling::CurvedLinear< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
   using Type = BoundaryHandling< FlagField_T, Stencil_T, UBB_T, Outlet_T, MO_SBB_T, MO_CLI_T >;


   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID,
                       const BlockDataID & particleFieldID, const shared_ptr<ParticleAccessor_T>& ac,
                       Vector3<real_t> uInfty) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), particleFieldID_( particleFieldID ),
      ac_( ac ), velocity_( uInfty )
   {}

   Type * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      WALBERLA_ASSERT_NOT_NULLPTR( storage );

      FlagField_T * flagField = block->getData< FlagField_T >( flagFieldID_ );
      PdfField_T *  pdfField  = block->getData< PdfField_T > ( pdfFieldID_ );
      auto * particleField = block->getData< lbm_mesapd_coupling::ParticleField_T > ( particleFieldID_ );

      const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

      Type * handling = new Type( "moving obstacle boundary handling", flagField, fluid,
                                  UBB_T( "UBB", UBB_Flag, pdfField, velocity_),
                                  Outlet_T( "Outlet", Outlet_Flag, pdfField, real_t(1) ),
                                  MO_SBB_T (  "MEM_BB",  MEM_BB_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block ),
                                  MO_CLI_T( "MEM_CLI", MEM_CLI_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block )  );

      const auto ubb = flagField->getFlag( UBB_Flag );
      const auto outlet = flagField->getFlag( Outlet_Flag );

      CellInterval domainBB = storage->getDomainCellBB();

      domainBB.xMin() -= cell_idx_c( FieldGhostLayers ); // -1
      domainBB.xMax() += cell_idx_c( FieldGhostLayers ); // cellsX

      domainBB.yMin() -= cell_idx_c( FieldGhostLayers ); // 0
      domainBB.yMax() += cell_idx_c( FieldGhostLayers ); // cellsY+1

      domainBB.zMin() -= cell_idx_c( FieldGhostLayers ); // 0
      domainBB.zMax() += cell_idx_c( FieldGhostLayers ); // cellsZ+1

      // BOTTOM
      // -1........-1........-1
      // cellsX+1..cellsY+1..-1
      CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
      storage->transformGlobalToBlockLocalCellInterval( bottom, *block );
      handling->forceBoundary( ubb, bottom );

      // TOP
      // -1........-1........cellsZ+1
      // cellsX+1..cellsY+1..cellsZ+1
      CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      storage->transformGlobalToBlockLocalCellInterval( top, *block );
      handling->forceBoundary( outlet, top );

      handling->fillWithDomain( FieldGhostLayers );

      return handling;
   }

private:
   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID particleFieldID_;

   shared_ptr<ParticleAccessor_T> ac_;
   const Vector3<real_t> velocity_;

}; // class MyBoundaryHandling

/*
 * Functionality for continuous logging of sphere properties during initial (resting) simulation
 */
template< typename ParticleAccessor_T>
class RestingSphereForceEvaluator
{
public:
   RestingSphereForceEvaluator( const shared_ptr< ParticleAccessor_T > & ac, walberla::id_t sphereUid,
                                uint_t averageFrequency, const std::string & basefolder ) :
         ac_( ac ), sphereUid_( sphereUid ), averageFrequency_( averageFrequency )
   {
      WALBERLA_ROOT_SECTION() {
         // initially write header in file
         std::ofstream file;
         filename_ = basefolder;
         filename_ += "/log_init.txt";
         file.open(filename_.c_str());
         file << "# f_z_current f_z_average f_x f_y\n";
         file.close();
      }
   }

   // evaluate the acting drag force on a fixed sphere
   void operator()(const uint_t timestep)
   {

      // average the force over averageFrequency consecutive timesteps since there are small fluctuations in the force
      auto currentForce = getForce();
      currentAverage_ += currentForce[2];

      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( filename_.c_str(), std::ofstream::app );
         file.precision(8);
         file << timestep << "\t" << currentForce[2] << "\t" << dragForceNew_
              << "\t" << currentForce[0] << "\t" << currentForce[1] << std::endl;
         file.close();
      }

      if( timestep % averageFrequency_ != 0) return;

      dragForceOld_ = dragForceNew_;
      dragForceNew_ = currentAverage_ / real_c( averageFrequency_ );
      currentAverage_ = real_t(0);

   }

   // Return the relative temporal change in the drag force
   real_t getForceDiff() const
   {
      return std::fabs( ( dragForceNew_ - dragForceOld_ ) / dragForceNew_ );
   }

   real_t getDragForce() const
   {
      return dragForceNew_;
   }

private:

   // Obtain the drag force acting on the sphere
   Vector3<real_t> getForce()
   {
      Vector3<real_t> force(real_t( 0 ));
      size_t idx = ac_->uidToIdx(sphereUid_);
      if( idx != ac_->getInvalidIdx())
      {
         if(!mesa_pd::data::particle_flags::isSet( ac_->getFlags(idx), mesa_pd::data::particle_flags::GHOST))
         {
            force = ac_->getHydrodynamicForce(idx);
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( force[0], mpi::SUM );
         mpi::allReduceInplace( force[1], mpi::SUM );
         mpi::allReduceInplace( force[2], mpi::SUM );
      }
      return force;
   }

   shared_ptr< ParticleAccessor_T > ac_;
   const walberla::id_t sphereUid_;

   real_t currentAverage_ = real_t(0);
   uint_t averageFrequency_;
   std::string filename_;
   real_t dragForceOld_ = real_t(0);
   real_t dragForceNew_ = real_t(0);

};

/*
 * Functionality for continuous logging of sphere properties during moving simulation
 */
template< typename ParticleAccessor_T>
class MovingSpherePropertyEvaluator
{
public:
   MovingSpherePropertyEvaluator( const shared_ptr< ParticleAccessor_T > & ac, walberla::id_t sphereUid, const std::string & basefolder,
                                  const Vector3<real_t> & u_infty, const real_t Galileo, const real_t GalileoSim,
                                  const real_t gravity, const real_t viscosity, const real_t diameter,
                                  const real_t densityRatio ) :
         ac_( ac ), sphereUid_( sphereUid ),
         u_infty_( u_infty ), gravity_( gravity ), viscosity_( viscosity ),
      diameter_( diameter ), densityRatio_( densityRatio )
   {

      WALBERLA_ROOT_SECTION()
      {
         std::ofstream fileSetup;
         std::string filenameSetup = basefolder;
         filenameSetup += "/setup.txt";
         fileSetup.open( filenameSetup.c_str() );

         fileSetup << "# setup parameters: \n";
         fileSetup << "processes = " << MPIManager::instance()->numProcesses() << "\n";
         fileSetup << "Galileo number (targeted) = " << Galileo << "\n";
         fileSetup << "Galileo number (simulated) = " << GalileoSim << "\n";
         fileSetup << "gravity = " << gravity << "\n";
         fileSetup << "viscosity = " << viscosity << "\n";
         fileSetup << "diameter = " << diameter << "\n";
         fileSetup << "uInfty = " << u_infty_[2] << "\n";
         real_t u_ref = std::sqrt( std::fabs(densityRatio - real_t(1)) * gravity * diameter );
         fileSetup << "u_{ref} = " << u_ref << "\n";
         fileSetup.close();
      }

      filenameRes_ = basefolder;
      filenameRes_ += "/log_results.txt";
      filenameAdd_ = basefolder;
      filenameAdd_ += "/log_additional.txt";

      WALBERLA_ROOT_SECTION()
      {
         // initially write header in file
         std::ofstream fileRes;
         fileRes.open( filenameRes_.c_str() );
         // raw data
         fileRes << "#\t x_p_x\t x_p_y\t x_p_z\t u_p_x\t u_p_y\t u_p_z\t omega_p_x\t omega_p_y\t omega_p_z\t u_{pr}_x\t u_{pr}_y\t u_{pr}_z\t ";
         // processed data, 14 - 18
         fileRes << "Re\t u_{pV}\t u_{pH}\t omega_{pV}\t omega_{pH}\t alpha\n";
         fileRes.close();

         std::ofstream fileAdd;
         fileAdd.open( filenameAdd_.c_str() );
         fileAdd << "# forceX forceY forceZ torqueX torqueY torqueZ\n";
         fileAdd.close();
      }
   }

   // evaluate and write the sphere properties
   void operator()(const uint_t timestep)
   {
      Vector3<real_t> transVel( real_t(0) );
      Vector3<real_t> angularVel( real_t(0) );
      Vector3<real_t> pos( real_t(0) );

      Vector3<real_t> force( real_t(0) );
      Vector3<real_t> torque( real_t(0) );

      size_t idx = ac_->uidToIdx(sphereUid_);
      if( idx != ac_->getInvalidIdx())
      {
         if(!mesa_pd::data::particle_flags::isSet( ac_->getFlags(idx), mesa_pd::data::particle_flags::GHOST))
         {
            pos = ac_->getPosition(idx);
            transVel = ac_->getLinearVelocity(idx);
            angularVel = ac_->getAngularVelocity(idx);
            force = ac_->getHydrodynamicForce(idx);
            torque = ac_->getHydrodynamicTorque(idx);
         }
      }

      std::vector<real_t> particlePos(3);
      particlePos[0]=pos[0]; particlePos[1]=pos[1]; particlePos[2]=pos[2];

      std::vector<real_t> particleTransVel(3);
      particleTransVel[0]=transVel[0]; particleTransVel[1]=transVel[1]; particleTransVel[2]=transVel[2];

      std::vector<real_t> particleAngularVel(3);
      particleAngularVel[0]=angularVel[0]; particleAngularVel[1]=angularVel[1]; particleAngularVel[2]=angularVel[2];

      std::vector<real_t> particleForce(3);
      particleForce[0]=force[0]; particleForce[1]=force[1]; particleForce[2]=force[2];

      std::vector<real_t> particleTorque(3);
      particleTorque[0]=torque[0]; particleTorque[1]=torque[1]; particleTorque[2]=torque[2];

      // reduce to root
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( particlePos, mpi::SUM );
         mpi::reduceInplace( particleTransVel, mpi::SUM );
         mpi::reduceInplace( particleAngularVel, mpi::SUM );
         mpi::reduceInplace( particleForce, mpi::SUM );
         mpi::reduceInplace( particleTorque, mpi::SUM );
      }

      particleHeight_ = particlePos[2];

      //only root process has all the data
      WALBERLA_ROOT_SECTION()
      {

         Vector3<real_t> u_p_r( particleTransVel[0] - u_infty_[0], particleTransVel[1] - u_infty_[1], particleTransVel[2] - u_infty_[2]);
         real_t u_p_H = std::sqrt( u_p_r[0] * u_p_r[0] + u_p_r[1] * u_p_r[1]);
         real_t u_p_V = u_p_r[2];

         real_t omega_p_H = std::sqrt( particleAngularVel[0] * particleAngularVel[0] + particleAngularVel[1] * particleAngularVel[1] );
         real_t omega_p_V = particleAngularVel[2];

         real_t u_ref = std::sqrt( std::fabs(densityRatio_ - real_t(1)) * gravity_ * diameter_ );
         real_t Reynolds = u_p_r.length() * diameter_ / viscosity_;

         // results
         std::ofstream fileRes;
         fileRes.open( filenameRes_.c_str(), std::ofstream::app );
         fileRes.precision(8);
         fileRes << timestep
                 << "\t" << particlePos[0] << "\t" << particlePos[1] << "\t" << particlePos[2]
                 << "\t" << particleTransVel[0] << "\t" << particleTransVel[1] << "\t" << particleTransVel[2]
                 << "\t" << particleAngularVel[0] << "\t" << particleAngularVel[1] << "\t" << particleAngularVel[2]
                 << "\t" << u_p_r[0] << "\t" << u_p_r[1] << "\t" << u_p_r[2]
                 << "\t" << Reynolds << "\t" << u_p_V/u_ref << "\t" << u_p_H/u_ref
                 << "\t" << omega_p_V*(diameter_/u_ref) << "\t" << omega_p_H*(diameter_/u_ref)
                 << "\t" << std::atan( u_p_H/std::fabs(u_p_V) )
                 << std::endl;
         fileRes.close();

         // additional
         std::ofstream fileAdd;
         fileAdd.open( filenameAdd_.c_str(), std::ofstream::app );
         fileAdd.precision(8);
         fileAdd << timestep
                 << "\t" << particleForce[0] << "\t" << particleForce[1] << "\t" << particleForce[2]
                 << "\t" << particleTorque[0] << "\t" << particleTorque[1] << "\t" << particleTorque[2]
                 << std::endl;
         fileAdd.close();
      }
   }

   real_t getParticleHeight() const
   {
      return particleHeight_;
   }

private:

   shared_ptr< ParticleAccessor_T > ac_;
   const walberla::id_t sphereUid_;

   std::string filenameRes_;
   std::string filenameAdd_;

   Vector3<real_t> u_infty_;
   real_t gravity_;
   real_t viscosity_;
   real_t diameter_;
   real_t densityRatio_;

   real_t particleHeight_;

};

/*
 * Result evaluation for the written VTK files
 *
 * input: vel (fluid velocity), u_p (particle velocity), u_infty (inflow velocity)
 *
 * needed data:
 * relative flow velocity: vel_r = vel - u_p
 * relative particle velocity: u_p_r = u_p - u_infty
 *
 * magnitude of relative particle velocity in horizontal plane: u_p_H = sqrt( (u_p_r)_x^2 + (u_p_r)_y^2 )
 * unit vector of particle motion in horizontal plane: e_p_H = ( (u_p_r)_x, (u_p_r)_y, 0) / u_p_H
 * unit vector perpendicular to e_p_H and e_z:  e_p_Hz_perp = ( -(u_p_r)_y, (u_p_r)_x, 0) / u_p_H
 * unit vector of particle motion: e_p_parallel = u_p_r / ||u_p_r||
 * unit vector perpendicular to e_p_parallel and e_p_Hz_perp: e_p_perp = e_p_Hz_perp x e_p_parallel
 * projected fluid velocities: vel_r_parallel = vel_r * (-e_p_parallel)
 *                             vel_r_perp     = vel_r * e_p_perp
 *                             vel_r_Hz_perp  = vel_r * e_p_Hz_perp
 */
template< typename ParticleAccessor_T>
class VTKInfoLogger
{
public:

   VTKInfoLogger( SweepTimeloop* timeloop, const shared_ptr< ParticleAccessor_T > & ac, walberla::id_t sphereUid,
         const std::string & baseFolder, const Vector3<real_t>& u_infty ) :
         timeloop_( timeloop ), ac_( ac ), sphereUid_( sphereUid ), baseFolder_( baseFolder ), u_infty_( u_infty )
   { }

   void operator()()
   {
      Vector3<real_t> u_p( real_t(0) );

      size_t idx = ac_->uidToIdx(sphereUid_);
      if( idx != ac_->getInvalidIdx())
      {
         if(!mesa_pd::data::particle_flags::isSet( ac_->getFlags(idx), mesa_pd::data::particle_flags::GHOST))
         {
            u_p = ac_->getLinearVelocity(idx);
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( u_p[0], mpi::SUM );
         mpi::allReduceInplace( u_p[1], mpi::SUM );
         mpi::allReduceInplace( u_p[2], mpi::SUM );

      }

      Vector3<real_t> u_p_r = u_p - u_infty_;
      real_t u_p_r_sqr_mag = u_p_r.sqrLength();

      real_t u_p_H = std::sqrt( u_p_r[0] * u_p_r[0] + u_p_r[1] * u_p_r[1]);

      Vector3<real_t> e_p_H (u_p_r[0], u_p_r[1], real_t(0));
      e_p_H /= u_p_H;
      Vector3<real_t> e_p_Hz_perp (-u_p_r[1], u_p_r[0], real_t(0));
      e_p_Hz_perp /= u_p_H;

      Vector3<real_t> e_p_parallel = u_p_r / u_p_r.length();
      Vector3<real_t> e_p_perp     = e_p_Hz_perp%e_p_parallel;

      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;

         std::string filename(baseFolder_);
         filename += "/log_vtk/log_vtk_";
         filename += std::to_string( timeloop_->getCurrentTimeStep() );
         filename += ".txt";

         file.open( filename.c_str() );
         file.precision(8);

         file << "u_p = "           << u_p << "\n";
         file << "u_p_r_2 = "       << u_p_r_sqr_mag << "\n\n";

         file << "e_{pH} = "        << e_p_H        << "\n";
         file << "e_{pparallel} = " << e_p_parallel << "\n";
         file << "e_{pperp} =  "    << e_p_perp     << "\n";
         file << "e_{pHzperp} = "   << e_p_Hz_perp  << "\n";

         file.close();
      }
   }

private:

   SweepTimeloop* timeloop_;

   shared_ptr< ParticleAccessor_T > ac_;
   const walberla::id_t sphereUid_;

   std::string baseFolder_;
   Vector3<real_t> u_infty_;

};


/*
 * This extensive, physical test case simulates a single, heavy sphere in ambient fluid flow.
 * It is based on the benchmark proposed in
 * Uhlman, Dusek - "The motion of a single heavy sphere in ambient fluid: A benchmark for interface-resolved
 *                  particulate flow simulations with significant relative velocities" IJMF (2014),
 *                  doi: 10.1016/j.ijmultiphaseflow.2013.10.010
 * Results for LBM done with waLBerla are published in
 * Rettinger, Ruede - "A comparative study of fluid-particle coupling methods for fully resolved
 *                     lattice Boltzmann simulations" CaF (2017),
 *                     doi: 10.1016/j.compfluid.2017.05.033
 * The setup and the benchmark itself are explained in detail in these two papers.
 * Several logging files are written to evaluate the benchmark quantities.
 *
 * This case is the same as can be found in apps/benchmarks/MotionSingleHeavySphere
 * but using the lbm_mesapd_coupling, instead of the pe_coupling.
 * Compared to this previous work, significantly different outcomes can be observed when using the CLI boundary condition.
 * This can be traced back to the distinct way of force averaging.
 * The force averaging applied in the previous study (two LBM steps, then average) introduces a error in the Galilean invariance.
 * Surprisingly, this error can result in better agreement to the reference solutions.
 *
 */
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   bool noViscosityIterations  = false;
   bool vtkIOInit = false;
   bool longDom   = false;
   bool broadDom   = false;
   bool useOmegaBulkAdaption = false;
   real_t bulkViscRateFactor = real_t(1);
   uint_t averageForceTorqueOverTwoTimeStepsMode = 1; // 0: no, 1: running, 2: 2LBM
   bool conserveMomentum = false;
   bool useDiffusiveScaling = true;
   bool useVelocityVerlet   = true;
   bool initFluidVelocity = true;
   bool vtkOutputGhostlayers = false;

   MEMVariant memVariant = MEMVariant::BB;
   std::string reconstructorType = "Grad"; // Eq, EAN, Ext, Grad

   real_t diameter = real_t(18);
   uint_t XBlocks = uint_t(4);
   uint_t YBlocks = uint_t(4);

   real_t viscosity = real_t(0.01); // used in diffusive scaling
   real_t uIn = real_t(0.02); // used for acoustic scaling

   /*
    * 0: custom
    * 1: Ga = 144
    * 2: Ga = 178.46
    * 3: Ga = 190
    * 4: Ga = 250
    */
   uint_t simulationCase = 1;
   real_t Galileo = real_t(1);
   real_t Re_target = real_t(1);
   real_t timestepsNonDim = real_t(1);


   real_t vtkWriteSpacingNonDim = real_t(3); // write every 3 non-dim timesteps
   uint_t vtkNumOutputFiles = 10; // write only 10 vtk files: the 10 final ones

   std::string folderEnding = "";

   ////////////////////////////
   // COMMAND LINE ARGUMENTS //
   ////////////////////////////

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--case"   ) == 0 ) {simulationCase = uint_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--Ga"  ) == 0 ){Galileo = real_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--Re"  ) == 0 ){Re_target = real_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--timestepsNonDim"  ) == 0 ){timestepsNonDim = real_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--diameter"  ) == 0 ){diameter = real_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--noViscosityIterations" ) == 0 ){noViscosityIterations = true; continue;}
      if( std::strcmp( argv[i], "--viscosity"  ) == 0 ){viscosity = real_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--uIn"  ) == 0 ){uIn = real_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--vtkIOInit" ) == 0 ){vtkIOInit = true; continue;}
      if( std::strcmp( argv[i], "--longDom"   ) == 0 ){longDom = true; continue;}
      if( std::strcmp( argv[i], "--broadDom"   ) == 0 ){broadDom = true; continue;}
      if( std::strcmp( argv[i], "--variant"    ) == 0 ){memVariant = to_MEMVariant( argv[++i] ); continue;}
      if( std::strcmp( argv[i], "--useOmegaBulkAdaption"   ) == 0 ){useOmegaBulkAdaption = true; continue;}
      if( std::strcmp( argv[i], "--bvrf"   ) == 0 ){bulkViscRateFactor = real_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--forceTorqueAverageMode"   ) == 0 ){averageForceTorqueOverTwoTimeStepsMode = uint_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--conserveMomentum"   ) == 0 ){conserveMomentum = true; continue;}
      if( std::strcmp( argv[i], "--useDiffusiveScaling"   ) == 0 ){useDiffusiveScaling = true; continue;}
      if( std::strcmp( argv[i], "--XBlocks"   ) == 0 ){XBlocks = uint_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--YBlocks"   ) == 0 ){YBlocks = uint_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--folderEnding"   ) == 0 ){folderEnding = argv[++i]; continue;}
      if( std::strcmp( argv[i], "--reconstructorType" ) == 0 ) {reconstructorType = argv[++i]; continue;}
      if( std::strcmp( argv[i], "--useEuler" )   == 0 ) {useVelocityVerlet = false; continue;}
      if( std::strcmp( argv[i], "--noFluidVelocityInit" )   == 0 ) {initFluidVelocity = false; continue;}
      if( std::strcmp( argv[i], "--vtkOutputGhostlayers" )   == 0 ) {vtkOutputGhostlayers = true; continue;}
      if( std::strcmp( argv[i], "--vtkNumOutputFiles"   ) == 0 ){vtkNumOutputFiles = uint_c( std::atof( argv[++i] ) ); continue;}
      if( std::strcmp( argv[i], "--vtkWriteSpacingNonDim"  ) == 0 ){vtkWriteSpacingNonDim = real_c( std::atof( argv[++i] ) ); continue;}
      WALBERLA_ABORT("command line argument unknown: " << argv[i] );
   }

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   const real_t radius = real_t(0.5) * diameter;
   uint_t xlength = uint_c( diameter * real_t(5.34) );
   uint_t ylength = xlength;
   uint_t zlength = uint_c( diameter * real_t(16) );

   if(broadDom)
   {
      xlength *= uint_t(2);
      ylength *= uint_t(2);
   }

   if(longDom)
   {
      zlength *= uint_t(3);
   }


   // estimate Reynolds number (i.e. inflow velocity) based on Galileo number
   // values are taken from the original simulation of Uhlmann, Dusek
   switch( simulationCase )
   {
   case 0:
      WALBERLA_LOG_INFO_ON_ROOT("Using custom simulation case -> you have to provide Ga, Re, and time steps via command line arguments!");
      break;
   case 1:
      Galileo = real_t(144);
      Re_target = real_t(185.08);
      timestepsNonDim = real_t(100);
      break;
   case 2:
      Galileo = real_t(178.46);
      Re_target = real_t(243.01);
      timestepsNonDim = real_t(250);
      break;
   case 3:
      Galileo = real_t(190);
      Re_target = real_t(262.71);
      timestepsNonDim = real_t(250);
      break;
   case 4:
      Galileo = real_t(250);
      Re_target = real_t(365.10);
      timestepsNonDim = real_t(510);
      break;
   default:
      WALBERLA_ABORT("Simulation case is not supported!");
   }

   if( useDiffusiveScaling)
   {
      // estimate fitting inflow velocity (diffusive scaling, viscosity is fixed)
      // attention: might lead to large velocities for small diameters
      uIn = Re_target * viscosity / diameter;
   } else {
      // estimate fitting viscosity (acoustic scaling: velocity is fixed)
      // note that this is different to the approach taking in the paper, where an diffusive scaling with visc = 0.01 is taken
      // attention: might lead to small viscosities for small diameters
      viscosity = uIn * diameter / Re_target;
   }


   real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);
   real_t omegaBulk = lbm_mesapd_coupling::omegaBulkFromOmega(omega, bulkViscRateFactor);

   Vector3<real_t> uInfty = Vector3<real_t>( real_t(0), real_t(0), uIn );

   const real_t densityRatio = real_t(1.5);

   const uint_t averageFrequency = uint_c( ( ( uint_c(Galileo) >= 200) ? real_t(500) : real_t(2) ) * diameter / uIn ); // for initial simulation
   const real_t convergenceLimit = real_t(1e-4);
   const real_t convergenceLimitGalileo = real_t(1e-4);
   const real_t dx = real_t(1);
   const real_t magicNumberTRT = lbm::collision_model::TRT::threeSixteenth;

   const uint_t timestepsInit = uint_c( ( ( uint_c(Galileo) >= 200) ? real_t(3000) : real_t(100) ) * diameter / uIn ); // maximum number of time steps for the initial simulation
   const uint_t writeFrequencyInit = uint_t(1000); // vtk write frequency init

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////
   const int numProcs = MPIManager::instance()->numProcesses();
   uint_t ZBlocks = uint_c(numProcs) / ( XBlocks * YBlocks );
   const uint_t XCells = xlength / XBlocks;
   const uint_t YCells = ylength / YBlocks;
   const uint_t ZCells = zlength / ZBlocks;

   if( (xlength != XCells * XBlocks) || (ylength != YCells * YBlocks) || (zlength != ZCells * ZBlocks) )
   {
      WALBERLA_ABORT("Domain partitioning does not fit to total domain size!");
   }

   WALBERLA_LOG_INFO_ON_ROOT("Motion settling sphere simulation of case " << simulationCase << " -> Ga = " << Galileo);
   WALBERLA_LOG_INFO_ON_ROOT("Diameter: " << diameter);
   WALBERLA_LOG_INFO_ON_ROOT("Domain: " << xlength << " x " << ylength << " x " << zlength);
   WALBERLA_LOG_INFO_ON_ROOT("Processes: " << XBlocks << " x " << YBlocks << " x " << ZBlocks);
   WALBERLA_LOG_INFO_ON_ROOT("Subdomains: " << XCells << " x " << YCells << " x " << ZCells);
   WALBERLA_LOG_INFO_ON_ROOT("Tau: " << real_t(1)/omega);
   WALBERLA_LOG_INFO_ON_ROOT("Viscosity: " << viscosity);
   WALBERLA_LOG_INFO_ON_ROOT("uIn: " << uIn);
   WALBERLA_LOG_INFO_ON_ROOT("Re_infty: " << uIn * diameter / viscosity);

   auto blocks = blockforest::createUniformBlockGrid( XBlocks, YBlocks, ZBlocks, XCells, YCells, ZCells, dx, true,
                                                      true, true, false );

   //////////////////
   // RPD COUPLING //
   //////////////////

   auto rpdDomain = std::make_shared<mesa_pd::domain::BlockForestDomain>(blocks->getBlockForestPointer());

   //init data structures
   auto ps = walberla::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = walberla::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor = walberla::make_shared<ParticleAccessor_T >(ps, ss);



   real_t xParticle = real_t(0);
   real_t yParticle = real_t(0);
   real_t zParticle = real_t(0);

   // root determines particle position, then broadcasts it
   WALBERLA_ROOT_SECTION()
   {
      if( simulationCase == 1 )
      {
         xParticle = real_c( xlength ) * real_t(0.5);
         yParticle = real_c( ylength ) * real_t(0.5);
      }
      else if( simulationCase == 4 )
      {
         // add random perturbance for chaotic regime
         walberla::math::seedRandomGenerator( std::mt19937::result_type(std::time(nullptr)) );
         xParticle = real_c( xlength ) * real_t(0.5) + walberla::math::realRandom( real_t(-0.5), real_t(0.5) );
         yParticle = real_c( ylength ) * real_t(0.5) + walberla::math::realRandom( real_t(-0.5), real_t(0.5) );

      }
      else
      {
         //add small perturbance to sphere position to break stability due to numerical symmetry
         real_t perturbance = real_t(0.35);

         xParticle = real_c( xlength ) * real_t(0.5) + perturbance;
         yParticle = real_c( ylength ) * real_t(0.5);
      }

      zParticle = (longDom) ? ( diameter * real_t(16) + real_c( xlength ) ) : real_c( xlength );
   }

   // broadcast to other ranks
   WALBERLA_MPI_SECTION()
   {
      mpi::broadcastObject( xParticle );
      mpi::broadcastObject( yParticle );
      mpi::broadcastObject( zParticle );
   }

   // add sphere
   Vector3<real_t> position( xParticle, yParticle, zParticle );
   auto sphereShape = ss->create<mesa_pd::data::Sphere>( radius );
   ss->shapes[sphereShape]->updateMassAndInertia(densityRatio);

   walberla::id_t sphereUid = 0;
   if (rpdDomain->isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), position ))
   {
      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(position);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);
      sphereUid = p.getUid();
   }
   mpi::allReduceInplace(sphereUid, mpi::SUM);

   // set up synchronization procedure
   const real_t overlap = real_t( 1.5 ) * dx;
   std::function<void(void)> syncCall;
   if( XBlocks <= uint_t(4) )
   {
      WALBERLA_LOG_INFO_ON_ROOT("Using next neighbor sync!")
      syncCall = [&ps,&rpdDomain,overlap](){
         mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
         syncNextNeighborFunc(*ps, *rpdDomain, overlap);
      };
      syncCall();
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT("Using ghost owner sync!")
      syncCall = [&ps,&rpdDomain,overlap](){
         mesa_pd::mpi::SyncGhostOwners syncGhostOwnersFunc;
         syncGhostOwnersFunc(*ps, *rpdDomain, overlap);
      };
      for(uint_t i = 0; i < uint_c(XBlocks/2); ++i) syncCall();
   }


   WALBERLA_LOG_INFO_ON_ROOT("Initial sphere position: " << position);

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // add omega bulk field
   BlockDataID omegaBulkFieldID = field::addToStorage<ScalarField_T>( blocks, "omega bulk field", omegaBulk, field::fzyx);


   // create the lattice model
   real_t lambda_e = lbm::collision_model::TRT::lambda_e( omega );
   real_t lambda_d = lbm::collision_model::TRT::lambda_d( omega, magicNumberTRT );
#ifdef WALBERLA_BUILD_WITH_CODEGEN
   WALBERLA_LOG_INFO_ON_ROOT("Using generated TRT-like lattice model!");
   LatticeModel_T latticeModel = LatticeModel_T(omegaBulkFieldID, lambda_d, lambda_e);
#else
   WALBERLA_LOG_INFO_ON_ROOT("Using waLBerla built-in TRT lattice model and ignoring omega bulk!");
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega, magicNumberTRT ) );
#endif

   // add PDF field
   // initial velocity in domain = inflow velocity
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field", latticeModel, initFluidVelocity ? uInfty : Vector3<real_t>(0), real_t(1), uint_t(1), field::fzyx );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // add particle field
   BlockDataID particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>( blocks, "particle field", accessor->getInvalidUid(), field::fzyx, FieldGhostLayers );


   // add boundary handling & initialize outer domain boundaries
   using BoundaryHandling_T = MyBoundaryHandling<ParticleAccessor_T>::Type;
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling<ParticleAccessor_T>( flagFieldID, pdfFieldID, particleFieldID, accessor, uInfty), "boundary handling" );

   // kernels
   mesa_pd::mpi::ReduceProperty reduceProperty;

   lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;
   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
   lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel averageHydrodynamicForceTorque;
   lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel initializeHydrodynamicForceTorqueForAveragingKernel;

   lbm_mesapd_coupling::RegularParticlesSelector sphereSelector;

   lbm_mesapd_coupling::MovingParticleMappingKernel<BoundaryHandling_T> movingParticleMappingKernel(blocks, boundaryHandlingID, particleFieldID);

   // mapping of sphere required by MEM variants
   // sets the correct flags
   if( memVariant == MEMVariant::CLI )
   {
      ps->forEachParticle(false, sphereSelector, *accessor, movingParticleMappingKernel, *accessor, MEM_CLI_Flag);
   }else {
      ps->forEachParticle(false, sphereSelector, *accessor, movingParticleMappingKernel, *accessor, MEM_BB_Flag);
   }

   // base folder to store all logs and vtk output
   std::string basefolder ("MSHS_");

   basefolder += std::to_string( uint_c( Galileo ) );
   basefolder += "_";
   basefolder += std::to_string( uint_c( diameter ) );
   basefolder += "_MEM_";
   basefolder += MEMVariant_to_string( memVariant );
   basefolder += "_bvrf";
   basefolder += std::to_string(int(bulkViscRateFactor));
   if( useOmegaBulkAdaption ) basefolder += "Adapt";

   basefolder += "_" + reconstructorType;

   if( longDom ) basefolder += "_longDom";
   if( broadDom ) basefolder += "_broadDom";
   if( conserveMomentum ) basefolder += "_conserveMomentum";
   if( !folderEnding.empty() ) basefolder += "_" + folderEnding;

   WALBERLA_LOG_INFO_ON_ROOT("Basefolder for simulation results: " << basefolder);

   // create base directory if it does not yet exist
   filesystem::path tpath( basefolder );
   if( !filesystem::exists( tpath ) )
      filesystem::create_directory( tpath );

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks

   blockforest::communication::UniformBufferedScheme< Stencil_T > optimizedPDFCommunicationScheme( blocks );
   optimizedPDFCommunicationScheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) ); // optimized sync

   blockforest::communication::UniformBufferedScheme< Stencil_T > fullPDFCommunicationScheme( blocks );
   fullPDFCommunicationScheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) ); // full sync

   // omega bulk adaption
   using OmegaBulkAdapter_T = lbm_mesapd_coupling::OmegaBulkAdapter<ParticleAccessor_T, decltype(sphereSelector)>;
   real_t defaultOmegaBulk = lbm_mesapd_coupling::omegaBulkFromOmega(omega, real_t(1));
   real_t adaptionLayerSize = real_t(2);
   shared_ptr<OmegaBulkAdapter_T> omegaBulkAdapter = make_shared<OmegaBulkAdapter_T>(blocks, omegaBulkFieldID, accessor, defaultOmegaBulk, omegaBulk, adaptionLayerSize, sphereSelector);
   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
      (*omegaBulkAdapter)(blockIt.get());
   }


   //////////////////////////////////////
   // TIME LOOP FOR INITIAL SIMULATION //
   //////////////////////////////////////

   // Initialization simulation: fixed sphere and simulate until convergence (of drag force)

   // create the timeloop
   SweepTimeloop timeloopInit( blocks->getBlockStorage(), timestepsInit );

   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
   auto bhSweep = BoundaryHandling_T::getBlockSweep( boundaryHandlingID );
   timeloopInit.add() << BeforeFunction( optimizedPDFCommunicationScheme, "LBM Communication" )
                  << Sweep(bhSweep, "Boundary Handling" );

   // stream + collide LBM step
#ifdef WALBERLA_BUILD_WITH_CODEGEN
   auto lbmSweep = LatticeModel_T::Sweep( pdfFieldID );
   timeloopInit.add() << Sweep( lbmSweep, "LB sweep" );
#else
   auto lbmSweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );
   timeloopInit.add() << Sweep( makeSharedSweep( lbmSweep ), "cell-wise LB sweep" );
#endif

   if( vtkIOInit )
   {
      auto pdfFieldVTKInit = vtk::createVTKOutput_BlockData( blocks, "fluid_field_init", writeFrequencyInit, 0, false, basefolder );

      field::FlagFieldCellFilter< FlagField_T > fluidFilterInit( flagFieldID );
      fluidFilterInit.addFlag( Fluid_Flag );
      pdfFieldVTKInit->addCellInclusionFilter( fluidFilterInit );

      pdfFieldVTKInit->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
      pdfFieldVTKInit->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

      timeloopInit.addFuncAfterTimeStep( vtk::writeFiles( pdfFieldVTKInit ), "VTK (fluid field data)" );
   }

   timeloopInit.addFuncAfterTimeStep( RemainingTimeLogger( timeloopInit.getNrOfTimeSteps(), real_t(30) ), "Remaining Time Logger" );

   ////////////////////////////////
   // EXECUTE INITIAL SIMULATION //
   ////////////////////////////////

   real_t gravity = real_t(1);
   real_t GalileoSim = real_t(1);
   real_t ReynoldsSim = real_t(1);
   real_t u_ref = real_t(1);

   WALBERLA_LOG_INFO_ON_ROOT("Starting initialization phase (sphere is kept fixed).");
   WALBERLA_LOG_INFO_ON_ROOT("Iterating, and adapting the viscosity, until the targeted Galileo number is set.");

   RestingSphereForceEvaluator<ParticleAccessor_T> forceEval( accessor, sphereUid, averageFrequency, basefolder );

   while (true) {
      WcTimingPool timeloopInitTiming;

      WALBERLA_LOG_INFO_ON_ROOT("(Re-)Starting initial simulation.");
      for( uint_t i = 1; i < timestepsInit; ++i ){

         ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );

         timeloopInit.singleStep( timeloopInitTiming );

         reduceProperty.operator()<mesa_pd::HydrodynamicForceTorqueNotification>(*ps);
         forceEval(timeloopInit.getCurrentTimeStep()+1);

         // check if the relative change in the average drag force is below the specified convergence criterion
         if (forceEval.getForceDiff() < convergenceLimit && i > 2 * std::max(averageFrequency, zlength) )
         {
            // conditions to break:
            // - force diff sufficiently small
            // - AND more time steps than twice the domain size in z direction to ensure new information due
            //   to possible viscosity change has traveled twice through domain,
            //   or twice the average frequency to have converged averages

            WALBERLA_LOG_INFO_ON_ROOT("Drag force converged at time " << timeloopInit.getCurrentTimeStep());
            break;
         }
         // if simulation gets unstable
         if( std::isnan(forceEval.getDragForce()) )
         {
            WALBERLA_ABORT("NAN value detected in drag force during initial simulation, exiting....");
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT("Initial simulation has ended.")

      //evaluate the gravitational force necessary to keep the sphere at a approximately fixed position
      gravity = forceEval.getDragForce() / ( (densityRatio - real_t(1) ) * diameter * diameter * diameter * math::pi / real_t(6) );
      GalileoSim = std::sqrt( ( densityRatio - real_t(1) ) * gravity * diameter * diameter * diameter ) / viscosity;
      ReynoldsSim = uIn * diameter / viscosity;
      u_ref = std::sqrt( std::fabs(densityRatio - real_t(1)) * gravity * diameter );

      WALBERLA_LOG_INFO_ON_ROOT("Acting gravity = " << gravity );
      WALBERLA_LOG_INFO_ON_ROOT("Simulated Galileo number = " << GalileoSim );
      WALBERLA_LOG_INFO_ON_ROOT("Targeted Galileo number = " << Galileo );
      WALBERLA_LOG_INFO_ON_ROOT("Reynolds number infty = " << ReynoldsSim );

      if ( noViscosityIterations )
      {
         timeloopInitTiming.logResultOnRoot();
         WALBERLA_LOG_INFO_ON_ROOT("Terminate iterations since viscosity should not be changed, flag \"--noViscosityIterations\"");
         break;
      }

      // if simulated Galileo number is close enough to targeted Galileo number, stop the initial simulation
      if( std::abs( GalileoSim - Galileo ) / Galileo < convergenceLimitGalileo )
      {
         timeloopInitTiming.logResultOnRoot();
         WALBERLA_LOG_INFO_ON_ROOT("Iterations converged, simulated Galileo number is close enough to targeted one");
         break;
      }

      // else update the simulation parameter accordingly and resume simulation
      real_t diff = GalileoSim/Galileo;
      viscosity = diff * viscosity;
      omega = lbm::collision_model::omegaFromViscosity( viscosity );

      lambda_e = lbm::collision_model::TRT::lambda_e( omega );
      lambda_d = lbm::collision_model::TRT::lambda_d( omega, magicNumberTRT );

      // also adapt omega bulk
      omegaBulk = lbm_mesapd_coupling::omegaBulkFromOmega(omega, bulkViscRateFactor);
      defaultOmegaBulk = lbm_mesapd_coupling::omegaBulkFromOmega(omega, real_t(1));
      omegaBulkAdapter->setDefaultOmegaBulk(defaultOmegaBulk);
      omegaBulkAdapter->setAdaptedOmegaBulk(omegaBulk);

      // iterate all blocks with an iterator 'block' and change the collision model
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         // get the field data out of the block
         auto pdf = blockIt->getData< PdfField_T > ( pdfFieldID );
#ifdef WALBERLA_BUILD_WITH_CODEGEN
         pdf->latticeModel().omega_magic_ = lambda_d;
         pdf->latticeModel().omega_visc_ = lambda_e;
#else
         pdf->latticeModel().collisionModel().resetWithMagicNumber( omega, magicNumberTRT );
#endif
         (*omegaBulkAdapter)(blockIt.get());
      }


      WALBERLA_LOG_INFO_ON_ROOT("==> Adapting viscosity:");
      WALBERLA_LOG_INFO_ON_ROOT("New viscosity = " << viscosity );
      WALBERLA_LOG_INFO_ON_ROOT("New omega = " << omega );
      WALBERLA_LOG_INFO_ON_ROOT("New omega bulk = " << omegaBulk );

   }

   if( averageForceTorqueOverTwoTimeStepsMode == 1 )
   {
      // maintain a good initial guess for averaging
      ps->forEachParticle(false, mesa_pd::kernel::SelectLocal(), *accessor, initializeHydrodynamicForceTorqueForAveragingKernel, *accessor );
   }
   ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );

   ///////////////
   // TIME LOOP //
   ///////////////

   // actual simulation: freely moving sphere with acting gravity

   // calculate the number of timesteps

   real_t dtSim = (averageForceTorqueOverTwoTimeStepsMode == 2) ? real_t(2) : real_t(1);
   const real_t t_ref = ( diameter / u_ref );
   const uint_t timesteps = uint_c( timestepsNonDim * t_ref / dtSim );

   // set vtk write frequency accordingly

   const uint_t writeFrequency = uint_c( vtkWriteSpacingNonDim * t_ref ); // vtk write frequency
   const uint_t initWriteCallsToSkip = uint_c(std::max( 0, int( timesteps - ( vtkNumOutputFiles - uint_t(1) ) * writeFrequency - uint_t(1) ) ) ); // write calls to be skipped

   WALBERLA_LOG_INFO_ON_ROOT("Starting simulation timeloop with  " << timesteps << " timesteps!");
   WALBERLA_LOG_INFO_ON_ROOT("Sphere is allowed to move freely under action of gravity");

   SweepTimeloop timeloop( blocks->getBlockStorage(), uint_c(dtSim) * timesteps );


   timeloop.addFuncBeforeTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps(), real_t(30) ), "Remaining Time Logger" );

   // vtk output
   auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, "fluid_field", writeFrequency, (vtkOutputGhostlayers) ? FieldGhostLayers : 0, false, basefolder );
   pdfFieldVTK->setInitialWriteCallsToSkip( initWriteCallsToSkip );

   pdfFieldVTK->addBeforeFunction( fullPDFCommunicationScheme );

   // function to output plane infos for vtk output
   pdfFieldVTK->addBeforeFunction(VTKInfoLogger<ParticleAccessor_T>( &timeloop, accessor, sphereUid, basefolder, uInfty ));

   // create folder for log_vtk files to not pollute the basefolder
   filesystem::path tpath2( basefolder+"/log_vtk" );
   if( !filesystem::exists( tpath2 ) )
      filesystem::create_directory( tpath2 );


   field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
   fluidFilter.addFlag( Fluid_Flag );
   pdfFieldVTK->addCellInclusionFilter( fluidFilter );

   pdfFieldVTK->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
   pdfFieldVTK->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

   timeloop.addFuncBeforeTimeStep( vtk::writeFiles( pdfFieldVTK ), "VTK (fluid field data)" );


   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
   timeloop.add() << BeforeFunction( optimizedPDFCommunicationScheme, "LBM Communication" )
                  << Sweep(bhSweep, "Boundary Handling" );

   // stream + collide LBM step
#ifdef WALBERLA_BUILD_WITH_CODEGEN
   timeloop.add() << Sweep( lbmSweep, "LB sweep" );
#else
   timeloop.add() << Sweep( makeSharedSweep( lbmSweep ), "cell-wise LB sweep" );
#endif


   SweepTimeloop timeloopAfterParticles( blocks->getBlockStorage(), timesteps );

   // sweep for updating the pe body mapping into the LBM simulation
   if( memVariant == MEMVariant::CLI )
      timeloopAfterParticles.add() << Sweep( lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T, BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, MEM_CLI_Flag, FormerMEM_Flag, sphereSelector, conserveMomentum), "Particle Mapping" );
   else
      timeloopAfterParticles.add() << Sweep( lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T, BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, MEM_BB_Flag, FormerMEM_Flag, sphereSelector, conserveMomentum), "Particle Mapping" );

   // sweep for restoring PDFs in cells previously occupied by particles
   if( reconstructorType == "EAN")
   {

      auto sphereNormalExtrapolationDirectionFinder = make_shared<lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder>(blocks);
      auto equilibriumAndNonEquilibriumSphereNormalReconstructor = lbm_mesapd_coupling::makeEquilibriumAndNonEquilibriumReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder, uint_t(3), true);
      auto reconstructionManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMEM_Flag, Fluid_Flag, equilibriumAndNonEquilibriumSphereNormalReconstructor, conserveMomentum);

      timeloopAfterParticles.add() << BeforeFunction( fullPDFCommunicationScheme, "PDF Communication" )
                                   << Sweep( makeSharedSweep(reconstructionManager), "PDF Restore" );

   } else if( reconstructorType == "Ext" )
   {
      auto sphereNormalExtrapolationDirectionFinder = make_shared<lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder>(blocks);
      auto extrapolationSphereNormalReconstructor = lbm_mesapd_coupling::makeExtrapolationReconstructor<BoundaryHandling_T,lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder,true>(blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder, uint_t(3), true);

      timeloopAfterParticles.add() << BeforeFunction( fullPDFCommunicationScheme, "PDF Communication" )
                                   << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMEM_Flag, Fluid_Flag, extrapolationSphereNormalReconstructor, conserveMomentum) ), "PDF Restore" );
   } else if( reconstructorType == "Grad")
   {
      auto gradReconstructor = lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, omega, false, true, true);

      timeloopAfterParticles.add() << BeforeFunction( fullPDFCommunicationScheme, "PDF Communication" )
                                   << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMEM_Flag, Fluid_Flag, gradReconstructor, conserveMomentum) ), "PDF Restore" );
   } else if( reconstructorType == "Eq")
   {
      timeloopAfterParticles.add() << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMEM_Flag, Fluid_Flag, conserveMomentum) ), "PDF Restore" );
   } else {
      WALBERLA_ABORT("Unknown reconstructor type " << reconstructorType);
   }

   // update bulk omega in all cells to adapt to changed particle position
   if( useOmegaBulkAdaption )
   {
      timeloopAfterParticles.add() << Sweep( makeSharedSweep(omegaBulkAdapter), "Omega Bulk Adapter");
   }

   Vector3<real_t> extForcesOnSphere( real_t(0), real_t(0), - gravity * ( densityRatio - real_t(1) ) * diameter * diameter * diameter * math::pi / real_t(6));
   lbm_mesapd_coupling::AddForceOnParticlesKernel addGravitationalForce(extForcesOnSphere);


   mesa_pd::kernel::VelocityVerletPreForceUpdate  vvIntegratorPreForce(dtSim);
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(dtSim);
   mesa_pd::kernel::ExplicitEuler explicitEulerIntegrator(dtSim);

   MovingSpherePropertyEvaluator<ParticleAccessor_T> movingSpherePropertyEvaluator( accessor, sphereUid, basefolder, uInfty, Galileo, GalileoSim, gravity, viscosity, diameter, densityRatio );


   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   const bool useOpenMP = false;

   // time loop
   for (uint_t i = 0; i < timesteps; ++i )
   {
      // perform a single simulation step -> this contains LBM and setting of the hydrodynamic interactions
      timeloop.singleStep( timeloopTiming );

      if( averageForceTorqueOverTwoTimeStepsMode == 2)
      {
         // store current force and torque in fOld and tOld
         ps->forEachParticle(useOpenMP, sphereSelector, *accessor, initializeHydrodynamicForceTorqueForAveragingKernel, *accessor );

         // reset f and t
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );

         timeloop.singleStep( timeloopTiming );


         // f = (f+fOld)/2
         ps->forEachParticle(useOpenMP, sphereSelector, *accessor, averageHydrodynamicForceTorque, *accessor );

         reduceProperty.operator()<mesa_pd::HydrodynamicForceTorqueNotification>(*ps);
      }
      else
      {
         reduceProperty.operator()<mesa_pd::HydrodynamicForceTorqueNotification>(*ps);

         if( averageForceTorqueOverTwoTimeStepsMode == 1 )
         {
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, averageHydrodynamicForceTorque, *accessor );
         }
      }


      timeloopTiming["RPD"].start();

      if( useVelocityVerlet)
      {
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPreForce, *accessor);
         syncCall();
      }

      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addHydrodynamicInteraction, *accessor );
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addGravitationalForce, *accessor );

      reduceProperty.operator()<mesa_pd::ForceTorqueNotification>(*ps);

      if( useVelocityVerlet ) {
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPostForce, *accessor);
      }
      else
      {
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, explicitEulerIntegrator, *accessor);
      }
      syncCall();

      timeloopTiming["RPD"].end();

      // evaluation
      timeloopTiming["Logging"].start();
      movingSpherePropertyEvaluator((i+1)*uint_c(dtSim));
      timeloopTiming["Logging"].end();

      // reset after logging here
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );

      timeloopAfterParticles.singleStep( timeloopTiming );

      if(movingSpherePropertyEvaluator.getParticleHeight() > real_c(zlength) - real_t(3) * diameter)
      {
         pdfFieldVTK->forceWrite(timeloop.getCurrentTimeStep());
         WALBERLA_LOG_WARNING_ON_ROOT("Sphere is too close to outflow, stopping simulation");
         break;
      }

   }

   timeloopTiming.logResultOnRoot();

   return EXIT_SUCCESS;
}

} // namespace motion_settling_sphere

int main( int argc, char **argv ){
   motion_settling_sphere::main(argc, argv);
}
