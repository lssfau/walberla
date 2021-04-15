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
//! \file MotionSingleHeavySphere.cpp
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

#include "pe/basic.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"
#include "pe/cr/ICR.h"
#include "pe/Types.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/partially_saturated_cells_method/all.h"
#include "pe_coupling/utility/all.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/BlockCellDataWriter.h"
#include "vtk/Initialization.h"
#include "vtk/VTKOutput.h"

#include <functional>
#include <memory>


namespace motion_single_heavy_sphere{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

//////////////
// TYPEDEFS //
//////////////

// PDF field, flag field & body field
using LatticeModel_T = lbm::D3Q19<lbm::collision_model::TRT, false>;

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
using BodyField_T = GhostLayerField<pe::BodyID, 1>;

using BodyAndVolumeFraction_T = std::pair<pe::BodyID, real_t>;
using BodyAndVolumeFractionField_T = GhostLayerField<std::vector<BodyAndVolumeFraction_T>, 1>;

const uint_t FieldGhostLayers = 1;

// boundary handling
using UBB_T = lbm::SimpleUBB<LatticeModel_T, flag_t>;
using Outlet_T = lbm::SimplePressure<LatticeModel_T, flag_t>;

using MEM_BB_T = pe_coupling::SimpleBB<LatticeModel_T, FlagField_T>;
using MEM_CLI_T = pe_coupling::CurvedLinear<LatticeModel_T, FlagField_T>;
using MEM_MR_T = pe_coupling::CurvedQuadratic<LatticeModel_T, FlagField_T>;

using BoundaryHandling_T = BoundaryHandling<FlagField_T, Stencil_T, UBB_T, Outlet_T, MEM_BB_T, MEM_CLI_T, MEM_MR_T>;

using BodyTypeTuple = std::tuple<pe::Sphere>;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag   ( "fluid" );

const FlagUID UBB_Flag     ( "velocity bounce back" );
const FlagUID Outlet_Flag  ( "outlet" );

const FlagUID MEM_BB_Flag   ( "moving obstacle BB" );
const FlagUID MEM_CLI_Flag  ( "moving obstacle CLI" );
const FlagUID MEM_MR_Flag   ( "moving obstacle MR" );
const FlagUID FormerMEM_Flag( "former moving obstacle" );


//////////////////////////////
// Coupling algorithm enums //
//////////////////////////////
enum MEMVariant { BB, CLI, MR };

MEMVariant to_MEMVariant( const std::string& s )
{
   if( s == "BB"  ) return MEMVariant::BB;
   if( s == "CLI" ) return MEMVariant::CLI;
   if( s == "MR"  ) return MEMVariant::MR;
   throw std::runtime_error("invalid conversion from text to MEMVariant");
}

std::string MEMVariant_to_string ( const MEMVariant& m )
{
   if( m == MEMVariant::BB  ) return "BB";
   if( m == MEMVariant::CLI ) return "CLI";
   if( m == MEMVariant::MR  ) return "MR";
   throw std::runtime_error("invalid conversion from MEMVariant to string");
}

enum PSMVariant { SC1W1, SC2W1, SC3W1, SC1W2, SC2W2, SC3W2 };

PSMVariant to_PSMVariant( const std::string& s )
{
   if( s == "SC1W1"  ) return PSMVariant::SC1W1;
   if( s == "SC2W1"  ) return PSMVariant::SC2W1;
   if( s == "SC3W1"  ) return PSMVariant::SC3W1;
   if( s == "SC1W2"  ) return PSMVariant::SC1W2;
   if( s == "SC2W2"  ) return PSMVariant::SC2W2;
   if( s == "SC3W2"  ) return PSMVariant::SC3W2;
   throw std::runtime_error("invalid conversion from text to PSMVariant");
}

std::string PSMVariant_to_string ( const PSMVariant& m )
{
   if( m == PSMVariant::SC1W1 ) return "SC1W1";
   if( m == PSMVariant::SC2W1 ) return "SC2W1";
   if( m == PSMVariant::SC3W1 ) return "SC3W1";
   if( m == PSMVariant::SC1W2 ) return "SC1W2";
   if( m == PSMVariant::SC2W2 ) return "SC2W2";
   if( m == PSMVariant::SC3W2 ) return "SC3W2";
   throw std::runtime_error("invalid conversion from PSMVariant to string");
}

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////
class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID, const BlockDataID & bodyFieldID, const BlockDataID & pdfFieldPreColID, const Vector3<real_t> & uInfty) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), bodyFieldID_( bodyFieldID ), pdfFieldPreColID_( pdfFieldPreColID ), velocity_( uInfty )
   {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:
   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID bodyFieldID_;
   const BlockDataID pdfFieldPreColID_;
   const Vector3<real_t> & velocity_;

}; // class MyBoundaryHandling


BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField  = block->getData< PdfField_T > ( pdfFieldID_ );
   BodyField_T * bodyField = block->getData< BodyField_T >( bodyFieldID_ );
   PdfField_T *  pdfFieldPreCol  = block->getData< PdfField_T > ( pdfFieldPreColID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "moving obstacle boundary handling", flagField, fluid,
                                                           UBB_T( "UBB", UBB_Flag, pdfField, velocity_),
                                                           Outlet_T( "Outlet", Outlet_Flag, pdfField, real_t(1) ),
                                                           MEM_BB_T (  "MEM_BB",  MEM_BB_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ),
                                                           MEM_CLI_T( "MEM_CLI", MEM_CLI_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ),
                                                           MEM_MR_T (  "MEM_MR",  MEM_MR_Flag, pdfField, flagField, bodyField, fluid, *storage, *block, pdfFieldPreCol ) );

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

/*
 * Functionality for continuous logging of sphere properties during initial (resting) simulation
 */
class RestingSphereForceEvaluator
{
public:
   RestingSphereForceEvaluator( SweepTimeloop* timeloop, const shared_ptr< StructuredBlockStorage > & blocks,
                                const ConstBlockDataID & bodyStorageID, uint_t averageFrequency, const std::string & basefolder ) :
      timeloop_( timeloop ), blocks_( blocks ), bodyStorageID_( bodyStorageID ), averageFrequency_( averageFrequency )
   {
      // initially write header in file
      std::ofstream file;
      filename_ = basefolder;
      filename_ += "/log_init.txt";
      file.open( filename_.c_str() );
      file << "# fDrag_current fDrag_average\n";
      file.close();
   }

   // evaluate the acting drag force on a fixed sphere
   void operator()()
   {
      const uint_t timestep (timeloop_->getCurrentTimeStep()+1);

      // average the force over averageFrequency consecutive timesteps since there are small fluctuations in the force
      real_t currentForce = calcForce();
      currentAverage_ += currentForce;

      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( filename_.c_str(), std::ofstream::app );
         file.precision(8);
         file << timestep << "\t" << currentForce << "\t" << dragForceNew_ << std::endl;
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

   real_t getForce() const
   {
      return dragForceNew_;
   }

private:

   // Obtain the drag force acting on the sphere by summing up all the process local parts of fZ
   real_t calcForce()
   {
      real_t forceZ = real_t( 0 );
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin<pe::Sphere>(*blockIt, bodyStorageID_ ); bodyIt != pe::BodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            forceZ += bodyIt->getForce()[2];
         }
      }
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( forceZ, mpi::SUM );
      }
      return forceZ;
   }

   SweepTimeloop* timeloop_;

   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID bodyStorageID_;

   real_t currentAverage_ = real_t(0);
   uint_t averageFrequency_;
   std::string filename_;
   real_t dragForceOld_ = real_t(0);
   real_t dragForceNew_ = real_t(0);

};

/*
 * Functionality for continuous logging of sphere properties during moving simulation
 */
class MovingSpherePropertyEvaluator
{
public:
   MovingSpherePropertyEvaluator( SweepTimeloop* timeloop, const shared_ptr< StructuredBlockStorage > & blocks,
                                  const ConstBlockDataID & bodyStorageID, const std::string & basefolder,
                                  const Vector3<real_t> & u_infty, const real_t Galileo, const real_t GalileoSim,
                                  const real_t gravity, const real_t viscosity, const real_t diameter,
                                  const real_t densityRatio, const uint_t numLBMSubCycles ) :
      timeloop_( timeloop ), blocks_( blocks ), bodyStorageID_( bodyStorageID ),
      u_infty_( u_infty ), gravity_( gravity ), viscosity_( viscosity ),
      diameter_( diameter ), densityRatio_( densityRatio ), numLBMSubCycles_( numLBMSubCycles )
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
         fileSetup << "numLBMSubCycles = " << numLBMSubCycles << "\n";
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
   void operator()()
   {
      const uint_t timestep ( timeloop_->getCurrentTimeStep() * numLBMSubCycles_ + 1 );

      Vector3<real_t> transVel( real_t(0) );
      Vector3<real_t> angularVel( real_t(0) );
      Vector3<real_t> pos( real_t(0) );

      Vector3<real_t> force( real_t(0) );
      Vector3<real_t> torque( real_t(0) );

      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin<pe::Sphere>(*blockIt, bodyStorageID_ ); bodyIt != pe::LocalBodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            pos= bodyIt->getPosition();
            transVel = bodyIt->getLinearVel();
            angularVel = bodyIt->getAngularVel();
         }
         for( auto bodyIt = pe::BodyIterator::begin<pe::Sphere>(*blockIt, bodyStorageID_ ); bodyIt != pe::BodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            force += bodyIt->getForce();
            torque += bodyIt->getTorque();
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
         mpi::reduceInplace( particlePos, mpi::SUM );
         mpi::reduceInplace( particleTransVel, mpi::SUM );
         mpi::reduceInplace( particleAngularVel, mpi::SUM );
         mpi::reduceInplace( particleForce, mpi::SUM );
         mpi::reduceInplace( particleTorque, mpi::SUM );
      }

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

private:
   SweepTimeloop* timeloop_;

   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID bodyStorageID_;

   std::string filenameRes_;
   std::string filenameAdd_;

   Vector3<real_t> u_infty_;
   real_t gravity_;
   real_t viscosity_;
   real_t diameter_;
   real_t densityRatio_;
   uint_t numLBMSubCycles_;

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
class VTKInfoLogger
{
public:

   VTKInfoLogger( SweepTimeloop* timeloop, const shared_ptr< StructuredBlockStorage > & blocks,
                  const ConstBlockDataID & bodyStorageID, const std::string & baseFolder,
                  const Vector3<real_t>& u_infty ) :
      timeloop_( timeloop ), blocks_( blocks ), bodyStorageID_( bodyStorageID ), baseFolder_( baseFolder ), u_infty_( u_infty )
   { }

   void operator()()
   {
      Vector3<real_t> u_p( real_t(0) );
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin<pe::Sphere>(*blockIt, bodyStorageID_ ); bodyIt != pe::LocalBodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            u_p = bodyIt->getLinearVel();
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

   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID bodyStorageID_;

   std::string baseFolder_;
   Vector3<real_t> u_infty_;

};

class PDFFieldCopy
{
public:
   PDFFieldCopy( const BlockDataID & pdfFieldID, const BlockDataID & pdfFieldPreColID )
      : pdfFieldID_( pdfFieldID ), pdfFieldPreColID_( pdfFieldPreColID )
   {}
   void operator()( IBlock * block )
   {
      auto pdfField = block->getData< PdfField_T>( pdfFieldID_ );
      auto pdfFieldPreCol = block->getData< PdfField_T>( pdfFieldPreColID_ );
      std::copy( pdfField->data(), pdfField->data() + pdfField->allocSize(), pdfFieldPreCol->data() );
      pdfFieldPreCol->resetLatticeModel( pdfField->latticeModel() );

   }
private:
   BlockDataID pdfFieldID_;
   BlockDataID pdfFieldPreColID_;
};

/*
 * This extensive, physical test case simulates a single, heavy sphere in ambient fluid flow.
 * It is based on the benchmark proposed in
 * Uhlmann, Dusek - "The motion of a single heavy sphere in ambient fluid: A benchmark for interface-resolved
 *                  particulate flow simulations with significant relative velocities" IJMF (2014),
 *                  doi: 10.1016/j.ijmultiphaseflow.2013.10.010
 * Results for LBM done with waLBerla are published in
 * Rettinger, Ruede - "A comparative study of fluid-particle coupling methods for fully resolved
 *                     lattice Boltzmann simulations" CaF (2017),
 *                     doi: 10.1016/j.compfluid.2017.05.033
 * The setup and the benchmark itself are explained in detail in these two papers.
 * Several logging files are written to evaluate the benchmark quantities.
 *
 * It is used to compare different (LBM-specific) fluid-particle coupling approaches:
 *  - momentum exchange method, with different boundary conditions:
 *    - BB (simple bounce-back)
 *    - CLI (central linear interpolated)
 *    - MR (multi reflection)
 *  - partially saturated cells method (i.e. Noble-Torczynski method), with different algorithmic variances:
 *    - solid collision operator variant: 1, 2, or 3
 *    - weighting factor variant: 1, or 2
 *
 * The results obtained by this benchmark should be compared to the reference spectral method results from Uhlmann & Dusek.
 * Furthermore, comparisons to the CFD-IBM (classical Navier-Stokes solver with immersed boundary method)
 * results from Uhlmann & Dusek are available in their paper.
 * New coupling algorithms or algorithmic changes to the exiting approaches should also be cross-checked with the
 * values in Rettinger & Ruede.
 *
 * The test case is split into two parts, since a certain, specified Galileo number should be simulated
 * but can not be chosen a-priori:
 * simulation is started with a fixed viscosity, and an inflow velocity estimated from the reference Reynolds number
 * initial simulation:
 *    run simulation with fixed sphere until force on sphere is converged
 *    chose gravity to balance this interaction force
 *    calculate Galileo number based on this gravitational acceleration
 *    if simulated Galileo number is close enough to targeted Galileo number, terminate, else adapt viscosity and repeat
 *
 * moving simulation:
 *    sphere is allowed to move freely with acting gravity, given by initial simulation
 *
 */
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   bool noViscosityIterations  = false;
   bool vtkIOInit = false;
   bool longDom   = false;

   bool useMEM = false;
   MEMVariant memVariant = MEMVariant::BB;
   bool usePSM = false;
   PSMVariant psmVariant = PSMVariant::SC3W2;

   real_t Galileo = real_t(144);
   real_t diameter = real_t(18);

   ////////////////////////////
   // COMMAND LINE ARGUMENTS //
   ////////////////////////////

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--Galileo"   ) == 0 )
      {
         Galileo   = real_c( std::atof( argv[++i] ) ); // targeted Galileo number
         continue;
      }
      if( std::strcmp( argv[i], "--diameter"  ) == 0 )
      {
         diameter  = real_c( std::atof( argv[++i] ) ); // diameter of the spheres, in cells per diameter
         continue;
      }
      if( std::strcmp( argv[i], "--noViscosityIterations" ) == 0 )
      {
         noViscosityIterations   = true; // keep viscosity unchanged in initial simulation
         continue;
      }
      if( std::strcmp( argv[i], "--vtkIOInit" ) == 0 )
      {
         vtkIOInit = true; // write vtk IO for initial simulation
         continue;
      }
      if( std::strcmp( argv[i], "--longDom"   ) == 0 )
      {
         longDom   = true; // use a domain that is extended by the original domain length in positive and negative streamwise direction
         continue;
      }
      if( std::strcmp( argv[i], "--useMEM"    ) == 0 ) // use momentum exchange coupling and the given MEM variant (BB, CLI, or MR)
      {
         useMEM = true;
         memVariant = to_MEMVariant( argv[++i] );
         continue;
      }
      if( std::strcmp( argv[i], "--usePSM"    ) == 0 ) // use partially saturated cells method and the given PSM variant (SC1W1, SC1W2, SC2W1, SC2W2, SC3W1, or SC3W2)
      {
         usePSM = true;
         psmVariant = to_PSMVariant( argv[++i] );
         continue;
      }
      WALBERLA_ABORT("command line argument unknown: " << argv[i] );
   }
   if( !useMEM && !usePSM ) WALBERLA_ABORT("No coupling method chosen via command line!\nChose either \"--useMEM MEMVARIANT\" or \"--usePSM PSMVARIANT\"!");

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   const real_t radius = real_t(0.5) * diameter;
   const uint_t xlength = uint_c( diameter * real_t(5.34) );
   const uint_t ylength = xlength;
   uint_t zlength = uint_c( diameter * real_t(16) );

   if (longDom)
   {
      zlength *= uint_t(3);
   }

   real_t viscosity = real_t(0.01);
   real_t Re_target = real_t(1);
   real_t timestepsNonDim = real_t(1);

   // estimate Reynolds number (i.e. inflow velocity) based on Galileo number
   // values are taken from the original simulation of Uhlmann, Dusek
   switch( int(Galileo) )
   {
   case 144:
      Re_target = real_t(185.08);
      timestepsNonDim = real_t(100);
      break;
   case 178:
      Re_target = real_t(243.01);
      timestepsNonDim = real_t(250);
      break;
   case 190:
      Re_target = real_t(262.71);
      timestepsNonDim = real_t(250);
      break;
   case 250:
      Re_target = real_t(365.10);
      timestepsNonDim = real_t(510);
      break;
   default:
      WALBERLA_ABORT("Galileo number is different from the usual ones (144, 178, 190, or 250). No estimate of the Reynolds number available. Add this case manually!");
   }

   // estimate fitting inflow velocity (diffusive scaling, viscosity is fixed)
   real_t uIn = Re_target * viscosity / diameter;

   real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);

   Vector3<real_t> uInfty = Vector3<real_t>( real_t(0), real_t(0), uIn );

   const real_t densityRatio = real_t(1.5);

   const uint_t averageFrequency = uint_c( ( ( uint_c(Galileo) >= 200) ? real_t(500) : real_t(2) ) * diameter / uIn ); // for initial simulation
   const real_t convergenceLimit = real_t(1e-4);
   const real_t convergenceLimitGalileo = real_t(1e-4);
   const real_t dx = real_t(1);
   const real_t magicNumberTRT = lbm::collision_model::TRT::threeSixteenth;

   const uint_t numLBMSubCycles = ( useMEM ) ? 2 : 1;
   const uint_t numPeSubCycles  = uint_c(numLBMSubCycles); // dtPE = dtLBM

   const uint_t timestepsInit = uint_c( ( ( uint_c(Galileo) >= 200) ? real_t(3000) : real_t(100) ) * diameter / uIn ); // maximum number of time steps for the initial simulation
   const uint_t writeFrequencyInit = uint_t(1000); // vtk write frequency init

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////
   const int numProcs = MPIManager::instance()->numProcesses();


   uint_t XBlocks;
   uint_t YBlocks;
   uint_t ZBlocks;
   if( numProcs > 192 )
   {
      XBlocks = uint_t(6);
      YBlocks = uint_t(6);
      ZBlocks = uint_c(numProcs) / ( XBlocks * YBlocks );
      WALBERLA_CHECK(numProcs % 36 == 0, "An integer multiple of 36 MPI ranks has to be used due to horizontal periodicity and domain partitioning requirements!");
   } else {
      XBlocks = uint_t(4);
      YBlocks = uint_t(4);
      ZBlocks = uint_c(MPIManager::instance()->numProcesses()) / ( XBlocks * YBlocks );
      WALBERLA_CHECK(numProcs % 16 == 0, "An integer multiple of 16 MPI ranks has to be used due to horizontal periodicity and domain partitioning requirements!");
   }
   const uint_t XCells = xlength / XBlocks;
   const uint_t YCells = ylength / YBlocks;
   const uint_t ZCells = zlength / ZBlocks;

   if( (xlength != XCells * XBlocks) || (ylength != YCells * YBlocks) || (zlength != ZCells * ZBlocks) )
   {
      WALBERLA_ABORT("Domain decomposition does not fit to total domain size!");
   }

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

   /////////////////
   // PE COUPLING //
   /////////////////

   // set up pe functionality
   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   auto bodyStorageID  = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");
   auto ccdID          = blocks->addBlockData(pe::ccd::createHashGridsDataHandling( globalBodyStorage, bodyStorageID ), "CCD");
   auto fcdID          = blocks->addBlockData(pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");

   // set up collision response, here DEM solver
   // in this test case, it is only used for the time integration
   pe::cr::DEM cr(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID, nullptr);

   // set up synchronization procedure
   const real_t overlap = real_t( 1.5 ) * dx;
   std::function<void(void)> syncCall;
   if( XBlocks <= uint_t(4) )
      syncCall = std::bind( pe::syncNextNeighbors<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );
   else
      syncCall = std::bind( pe::syncShadowOwners<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );


   real_t xParticle = real_t(0);
   real_t yParticle = real_t(0);
   real_t zParticle = real_t(0);

   // root determines particle position, then broadcasts it
   WALBERLA_ROOT_SECTION()
   {
      if( int(Galileo) == 144 )
      {
         xParticle = real_c( xlength ) * real_t(0.5);
         yParticle = real_c( ylength ) * real_t(0.5);
      }
      else if( int(Galileo) == 250 )
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
   const auto sphereMaterial = pe::createMaterial( "mySphereMat", densityRatio , real_t(0.5), real_t(0.1), real_t(0.1), real_t(0.24), real_t(200), real_t(200), real_t(0), real_t(0) );
   Vector3<real_t> position( xParticle, yParticle, zParticle );
   pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, position, radius, sphereMaterial );
   syncCall();

   WALBERLA_LOG_INFO_ON_ROOT("Initial sphere position: " << position);

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega, magicNumberTRT ) );

   // add PDF field
   // initial velocity in domain = inflow velocity
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (zyxf)", latticeModel, uInfty, real_t(1), uint_t(1), field::zyxf );

   // add PDF field (needed to store pre collision values for MEM_MR scheme)
   BlockDataID pdfFieldPreColID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "nqOdd field (zyxf)", latticeModel, uInfty, real_t(1), uint_t(1), field::zyxf );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::zyxf );

   // add body and volume fraction field
   BlockDataID bodyAndVolumeFractionFieldID = field::addToStorage< BodyAndVolumeFractionField_T >( blocks, "body and volume fraction field",
                                                                                                   std::vector<BodyAndVolumeFraction_T>(), field::zyxf, 0 );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
                                       MyBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID, pdfFieldPreColID, uInfty ), "boundary handling" );

   if( usePSM ){
      // mapping of sphere required by PSM methods

      pe_coupling::BodyAndVolumeFractionMapping bodyMapping( blocks, globalBodyStorage, bodyStorageID, bodyAndVolumeFractionFieldID );
      bodyMapping();

      //initialization of the PDFs inside the particles, important for PSM M3
      if( psmVariant == PSMVariant::SC1W1 || psmVariant == PSMVariant::SC2W1 || psmVariant == PSMVariant::SC3W1 )
         pe_coupling::initializeDomainForPSM< LatticeModel_T, 1 >( *blocks, pdfFieldID, bodyAndVolumeFractionFieldID );
      else
         pe_coupling::initializeDomainForPSM< LatticeModel_T, 2 >( *blocks, pdfFieldID, bodyAndVolumeFractionFieldID );

   }
   else
   {
      // mapping of sphere required by MEM variants
      // sets the correct flags

      if( memVariant == MEMVariant::CLI )
      {
         pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MEM_CLI_Flag );
      }else if ( memVariant == MEMVariant::MR ){
         pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MEM_MR_Flag );
      }else{
         pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MEM_BB_Flag );
      }
   }


   // base folder to store all logs and vtk output
   std::string basefolder ("MSHS_");

   basefolder += std::to_string( uint_c( Galileo ) );
   basefolder += "_";
   basefolder += std::to_string( uint_c( diameter ) );
   if ( usePSM )
   {
      basefolder += "_PSM_";
      basefolder += PSMVariant_to_string( psmVariant );
   }
   else
   {
      basefolder += "_MEM_";
      basefolder += MEMVariant_to_string( memVariant );
   }

   if( longDom )
      basefolder += "_longDom";

   WALBERLA_LOG_INFO_ON_ROOT("Basefolder for simulation results: " << basefolder);

   // create base directory if it does not yet exist
   filesystem::path tpath( basefolder );
   if( !filesystem::exists( tpath ) )
      filesystem::create_directory( tpath );

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   std::function< void () > commFunction;

   blockforest::communication::UniformBufferedScheme< Stencil_T > scheme( blocks );
   scheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) );
   commFunction = scheme;


   //////////////////////////////////////
   // TIME LOOP FOR INITIAL SIMULATION //
   //////////////////////////////////////

   // Initialization simulation: fixed sphere and simulate until convergence (of drag force)

   // create the timeloop
   SweepTimeloop timeloopInit( blocks->getBlockStorage(), timestepsInit );

   if( usePSM ){

      // add LBM communication function and boundary handling sweep
      timeloopInit.add()
            << BeforeFunction( commFunction, "LBM Communication" )
            << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

      // LBM stream collide sweep
      if( psmVariant == PSMVariant::SC1W1 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,1,1>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloopInit.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      } else if( psmVariant == PSMVariant::SC2W1 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,2,1>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloopInit.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      } else if( psmVariant == PSMVariant::SC3W1 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,3,1>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloopInit.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      } else if( psmVariant == PSMVariant::SC1W2 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,1,2>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloopInit.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      } else if( psmVariant == PSMVariant::SC2W2 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,2,2>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloopInit.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      } else if( psmVariant == PSMVariant::SC3W2 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,3,2>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloopInit.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      }

   } else {

      if( memVariant == MEMVariant::MR )
         timeloopInit.add() << Sweep( PDFFieldCopy( pdfFieldID, pdfFieldPreColID ), "pdf field copy" );

      auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );

      // collision sweep
      timeloopInit.add() << Sweep( lbm::makeCollideSweep( sweep ), "cell-wise LB sweep (collide)" );

      // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
      timeloopInit.add() << BeforeFunction( commFunction, "LBM Communication" )
                         << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

      // streaming & force evaluation
      timeloopInit.add() << Sweep( lbm::makeStreamSweep( sweep ), "cell-wise LB sweep (stream)" );

   }

   // evaluate the drag force acting on the sphere...
   shared_ptr< RestingSphereForceEvaluator > forceEval = walberla::make_shared< RestingSphereForceEvaluator >( &timeloopInit, blocks, bodyStorageID, averageFrequency, basefolder );
   timeloopInit.addFuncAfterTimeStep( SharedFunctor< RestingSphereForceEvaluator >( forceEval ), "Evaluating drag force" );

   // ...then reset the force
   timeloopInit.addFuncAfterTimeStep( pe_coupling::ForceTorqueOnBodiesResetter( blocks, bodyStorageID ), "Resetting force on body");

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

   timeloopInit.addFuncAfterTimeStep( RemainingTimeLogger( timeloopInit.getNrOfTimeSteps(), real_t(120) ), "Remaining Time Logger" );

   ////////////////////////////////
   // EXECUTE INITIAL SIMULATION //
   ////////////////////////////////

   real_t gravity = real_t(1);
   real_t GalileoSim = real_t(1);
   real_t ReynoldsSim = real_t(1);
   real_t u_ref = real_t(1);

   WALBERLA_LOG_INFO_ON_ROOT("Starting initialization phase (sphere is kept fixed).");
   WALBERLA_LOG_INFO_ON_ROOT("Iterating, and adapting the viscosity, until the targeted Galileo number is set.");

   while (true) {
      WcTimingPool timeloopInitTiming;

      WALBERLA_LOG_INFO_ON_ROOT("(Re-)Starting initial simulation.");
      for( uint_t i = 1; i < timestepsInit; ++i ){
         timeloopInit.singleStep( timeloopInitTiming );
         // check if the relative change in the average drag force is below the specified convergence criterion
         if (forceEval->getForceDiff() < convergenceLimit && i > 2 * std::max(averageFrequency, zlength) )
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
         if( std::isnan(forceEval->getForce()) )
         {
            WALBERLA_ABORT("NAN value detected in drag force during initial simulation, exiting....");
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT("Initial simulation has ended.")

            //evaluate the gravitational force necessary to keep the sphere at a approximately fixed position
            gravity = forceEval->getForce() / ( (densityRatio - real_t(1) ) * diameter * diameter * diameter * math::pi / real_t(6) );
      GalileoSim = std::sqrt( ( densityRatio - real_t(1) ) * gravity * diameter * diameter * diameter ) / viscosity;
      ReynoldsSim = uIn * diameter / viscosity;
      u_ref = std::sqrt( std::fabs(densityRatio - real_t(1)) * gravity * diameter );

      WALBERLA_LOG_INFO_ON_ROOT("Acting gravity (= interaction force) = " << gravity );
      WALBERLA_LOG_INFO_ON_ROOT("Simulated Galileo number = " << GalileoSim );
      WALBERLA_LOG_INFO_ON_ROOT("Targeted Galileo number = " << Galileo );
      WALBERLA_LOG_INFO_ON_ROOT("Reynolds number infty = " << ReynoldsSim );

      if ( noViscosityIterations )
      {
         timeloopInitTiming.logResultOnRoot();
         WALBERLA_LOG_INFO_ON_ROOT("Terminate iterations since viscosity should not be change, flag \"--noViscosityIterations\"");
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
      real_t newOmega = lbm::collision_model::omegaFromViscosity( viscosity );

      // iterate all blocks with an iterator 'block' and change the collision model
      for( auto block = blocks->begin(); block != blocks->end(); ++block )
      {
         // get the field data out of the block
         auto pdf = block->getData< PdfField_T > ( pdfFieldID );
         pdf->latticeModel().collisionModel().resetWithMagicNumber( newOmega, magicNumberTRT );
      }

      WALBERLA_LOG_INFO_ON_ROOT("==> Adapting viscosity:");
      WALBERLA_LOG_INFO_ON_ROOT("New viscosity = " << viscosity );
      WALBERLA_LOG_INFO_ON_ROOT("New omega = " << newOmega );

   }



   ///////////////
   // TIME LOOP //
   ///////////////

   // actual simulation: freely moving sphere with acting gravity

   // calculate the number of timesteps
   const real_t t_ref = ( diameter / u_ref );
   const uint_t timestepsLBM = uint_c( timestepsNonDim * t_ref );
   const uint_t timesteps = timestepsLBM / numLBMSubCycles + 1; // total number of time steps for the whole simulation

   // set vtk write frequency accordingly
   const real_t dtWriteNonDim = real_t(3); // write every 3 non-dim timesteps
   const uint_t nVTK = ( int(Galileo) != 178 ) ? 2 : 10; // write only 10 vtk files: the 10 final ones

   const uint_t writeFrequency = uint_c( dtWriteNonDim * t_ref ) / numLBMSubCycles; // vtk write frequency
   const uint_t initWriteCallsToSkip = uint_c(std::max( 0, int( timesteps - ( nVTK - uint_t(1) ) * writeFrequency - uint_t(1) ) ) ); // write calls to be skipped

   WALBERLA_LOG_INFO_ON_ROOT("Starting simulation timeloop with  " << timesteps << " timesteps!");
   WALBERLA_LOG_INFO_ON_ROOT("Sphere is allowed to move freely under action of gravity");

   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   if( usePSM ){

      timeloop.addFuncBeforeTimeStep( pe_coupling::BodyAndVolumeFractionMapping( blocks, globalBodyStorage, bodyStorageID, bodyAndVolumeFractionFieldID ), "volume fraction mapping" );

      for( uint_t lbmcycle = 0; lbmcycle < numLBMSubCycles; ++lbmcycle ){

         // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
         timeloop.add()
               << BeforeFunction( commFunction, "LBM Communication" )
               << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

         // LBM stream collide sweep
         if( psmVariant == PSMVariant::SC1W1 )
         {
            auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,1,1>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
            timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
         } else if( psmVariant == PSMVariant::SC2W1 )
         {
            auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,2,1>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
            timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
         } else if( psmVariant == PSMVariant::SC3W1 )
         {
            auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,3,1>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
            timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
         } else if( psmVariant == PSMVariant::SC1W2 )
         {
            auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,1,2>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
            timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
         } else if( psmVariant == PSMVariant::SC2W2 )
         {
            auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,2,2>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
            timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
         } else if( psmVariant == PSMVariant::SC3W2 )
         {
            auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,3,2>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
            timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
         }
      }

      if( numLBMSubCycles != uint_t(1) )
         timeloop.addFuncAfterTimeStep( pe_coupling::ForceTorqueOnBodiesScaler( blocks, bodyStorageID, real_t(1) / real_c(numLBMSubCycles) ), "Force averaging for several LBM steps" );

   }else{

      // sweep for updating the pe body mapping into the LBM simulation
      if( memVariant == MEMVariant::CLI )
         timeloop.add() << Sweep( pe_coupling::BodyMapping< LatticeModel_T, BoundaryHandling_T >( blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, MEM_CLI_Flag, FormerMEM_Flag ), "Body Mapping" );
      else if ( memVariant == MEMVariant::MR )
         timeloop.add() << Sweep( pe_coupling::BodyMapping< LatticeModel_T, BoundaryHandling_T >( blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, MEM_MR_Flag, FormerMEM_Flag ), "Body Mapping" );
      else
         timeloop.add() << Sweep( pe_coupling::BodyMapping< LatticeModel_T, BoundaryHandling_T >( blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, MEM_BB_Flag, FormerMEM_Flag ), "Body Mapping" );


      // reconstruct missing PDFs
      using ExtrapolationFinder_T = pe_coupling::SphereNormalExtrapolationDirectionFinder;
      ExtrapolationFinder_T extrapolationFinder( blocks, bodyFieldID );
      using Reconstructor_T = pe_coupling::ExtrapolationReconstructor<LatticeModel_T, BoundaryHandling_T, ExtrapolationFinder_T>;
      Reconstructor_T reconstructor( blocks, boundaryHandlingID, bodyFieldID, extrapolationFinder, true );
      timeloop.add() << Sweep( pe_coupling::PDFReconstruction< LatticeModel_T, BoundaryHandling_T, Reconstructor_T >
                               ( blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, reconstructor, FormerMEM_Flag, Fluid_Flag ), "PDF Restore" );

      for( uint_t lbmcycle = 0; lbmcycle < numLBMSubCycles; ++lbmcycle ){

         if( memVariant == MEMVariant::MR )
            timeloop.add() << Sweep( PDFFieldCopy( pdfFieldID, pdfFieldPreColID ), "pdf field copy" );

         auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );

         // collision sweep
         timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "cell-wise LB sweep (collide)" );

         // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
         timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                        << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

         // streaming
         timeloop.add() << Sweep( lbm::makeStreamSweep( sweep ), "cell-wise LB sweep (stream)" );
      }

      if( numLBMSubCycles != uint_t(1) )
         timeloop.addFuncAfterTimeStep( pe_coupling::ForceTorqueOnBodiesScaler( blocks, bodyStorageID, real_t(1) / real_c(numLBMSubCycles) ), "Force averaging for several LBM steps" );

   }

   // add gravity
   Vector3<real_t> extForcesOnSphere( real_t(0), real_t(0), - gravity * ( densityRatio - real_t(1) ) * diameter * diameter * diameter * math::pi / real_t(6));
   timeloop.addFuncAfterTimeStep( pe_coupling::ForceOnBodiesAdder( blocks, bodyStorageID, extForcesOnSphere ), "Add external forces (gravity and buoyancy)" );

   // evaluate the sphere properties
   timeloop.addFuncAfterTimeStep( MovingSpherePropertyEvaluator( &timeloop, blocks, bodyStorageID, basefolder, uInfty, Galileo, GalileoSim, gravity, viscosity, diameter, densityRatio, numLBMSubCycles ), "Evaluate sphere");

   // advance pe rigid body simulation
   timeloop.addFuncAfterTimeStep( pe_coupling::TimeStep( blocks, bodyStorageID, cr, syncCall, real_c(numLBMSubCycles), numPeSubCycles ), "pe Time Step" );

   // vtk output
   auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, "fluid_field", writeFrequency, FieldGhostLayers, false, basefolder );
   pdfFieldVTK->setInitialWriteCallsToSkip( initWriteCallsToSkip );

   blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync( blocks );
   pdfGhostLayerSync.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) );
   pdfFieldVTK->addBeforeFunction( pdfGhostLayerSync );

   // function to output plane infos for vtk output
   pdfFieldVTK->addBeforeFunction(VTKInfoLogger( &timeloop, blocks, bodyStorageID, basefolder, uInfty ));

   // create folder for log_vtk files to not pollute the basefolder
   filesystem::path tpath2( basefolder+"/log_vtk" );
   if( !filesystem::exists( tpath2 ) )
      filesystem::create_directory( tpath2 );


   field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
   fluidFilter.addFlag( Fluid_Flag );
   pdfFieldVTK->addCellInclusionFilter( fluidFilter );

   pdfFieldVTK->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
   pdfFieldVTK->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

   timeloop.addFuncAfterTimeStep( vtk::writeFiles( pdfFieldVTK ), "VTK (fluid field data)" );


   timeloop.addFuncAfterTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps(), real_t(120) ), "Remaining Time Logger" );

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;
   timeloop.run( timeloopTiming );
   timeloopTiming.logResultOnRoot();

   return EXIT_SUCCESS;
}

} // namespace motion_single_heavy_sphere

int main( int argc, char **argv ){
   motion_single_heavy_sphere::main(argc, argv);
}
