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
//! \author Brendan Waters <brendan.waters@sydney.edu.au>
//
//======================================================================================================================
#include "blockforest/all.h"

#include "core/all.h"

#include "domain_decomposition/all.h"

#include "field/all.h"

#include "geometry/all.h"

#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/field/initializer/all.h"
#include "lbm/vtk/VTKOutput.h"

#include "lbm_generated/communication/UniformGeneratedPdfPackInfo.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/evaluation/PerformanceEvaluation.h"

#include "timeloop/all.h"

// Generated files
#include "InfoHeader.h"


namespace walberla
{
constexpr uint_t FieldGhostLayer{2};

using StorageSpecification_T  = lbm::AdvectionDiffusionStorageSpecification;
using CommunicationStencil_T  = StorageSpecification_T::CommunicationStencil;

using PdfField_T              = lbm_generated::PdfField< StorageSpecification_T >;

using SweepCollection_T       = lbm::AdvectionDiffusionSweepCollection;

using FlagField_T             = FlagField< uint16_t >;
using BoundaryCollection_T    = lbm::AdvectionDiffusionBoundaryCollection< FlagField_T >;

using PackInfo_T              = lbm_generated::UniformGeneratedPdfPackInfo< PdfField_T >;

using VelocityField_T         = walberla::field::GhostLayerField<real_t, CommunicationStencil_T::D>;
using ScalarField_T           = walberla::field::GhostLayerField< real_t, 1 >;

using ExactField_T = walberla::field::GhostLayerField< real_t, 1 >;


const FlagUID fluidFlagUID("Fluid");
const FlagUID neumannFlagUID("Neumann");

class ScalarHelperFunctions
{
public:
   ScalarHelperFunctions( shared_ptr< StructuredBlockForest > blocks, const BlockDataID & exactDataId, const BlockDataID & scalarDataId, const real_t & kappa,const real_t & h0, const uint_t & ic, const Vector3< real_t > & lamda ):
   blocks_(blocks), exactDataId_(exactDataId), scalarDataId_(scalarDataId), kappa_(kappa), h0_(h0), ic_(ic), lam(lamda)
   {}

   void initDirichletBoundary()
   {  
      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         Block& b           = dynamic_cast< Block& >(*block);
         uint_t level       = b.getLevel();

         const real_t & dx = blocks_->dx( level );
         const real_t & dy = blocks_->dy( level );
         const real_t & dz = blocks_->dz( level );
         const Vector3< real_t > cell_face_offset_min { -1.0 - 0.5*dx , -1.0 - 0.5*dy, -1.0 - 0.5*dz};
         const Vector3< real_t > cell_face_offset_max { -1.0 + 0.5*dx , -1.0 + 0.5*dy, -1.0 + 0.5*dz};
         
         ScalarField_T* scalarField     = block->getData< ScalarField_T >(scalarDataId_);
         CellInterval xyz = scalarField->xyzSizeWithGhostLayer();


         for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
         {  
            const Vector3< real_t > p = blocks_->getBlockLocalCellCenter(*block, *cellIt);

            scalarField->get(*cellIt) =  real_t(1.0) + Average( p + cell_face_offset_min, p + cell_face_offset_max, dx, dy, dz );

         }
      }
   }

   void initVelField(  const BlockDataID & uId, const Vector3<real_t> & InitVel, const real_t & scaling = real_t(1) )
   {
      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         VelocityField_T* u = block->getData< VelocityField_T >(uId);

         CellInterval xyz = u->xyzSizeWithGhostLayer();

         for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
         {  
            u->get(*cellIt, 0)     = InitVel[0] * scaling;
            u->get(*cellIt, 1)     = InitVel[1] * scaling;
            u->get(*cellIt, 2)     = InitVel[2] * scaling;
         }
      }
   }

   void computeExact( const real_t & time, const real_t & scaling = real_t(1) )
   {
      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         Block& b           = dynamic_cast< Block& >(*block);
         uint_t level       = b.getLevel();

         ExactField_T * exactField     = block->getData< ExactField_T >(exactDataId_);

         CellInterval xyz = exactField->xyzSize();
         
         const real_t & dx = blocks_->dx( level );
         const real_t & dy = blocks_->dy( level );
         const real_t & dz = blocks_->dz( level );
         const Vector3< real_t > cell_face_offset_min { -1.0 - 0.5*dx - time*lam[0]*scaling, -1.0 - 0.5*dy - time*lam[1]*scaling, -1.0 - 0.5*dz - time*lam[2]*scaling};
         const Vector3< real_t > cell_face_offset_max { -1.0 + 0.5*dx - time*lam[0]*scaling, -1.0 + 0.5*dy - time*lam[1]*scaling, -1.0 + 0.5*dz - time*lam[2]*scaling};

         for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
         {  
            Cell globalCell;
            blocks_->transformBlockLocalToGlobalCell(globalCell, *block, *cellIt);
            
            const Vector3< real_t > p = blocks_->getCellCenter( globalCell, level );

            exactField -> get(*cellIt) = real_t(1.0) + Average( p + cell_face_offset_min, p + cell_face_offset_max, dx, dy, dz, time );
         }
        
      }
   }

   void calculateErrorNorms( const real_t & time, const real_t & scaling = real_t(1) )
   {
      computeExact( time , scaling );

      real_t er0 { 0.0 };
      real_t er1 { 0.0 };
      real_t er2 { 0.0 };
   
      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         ExactField_T * exactField  = block->getData< ExactField_T >(exactDataId_);
         ScalarField_T * cField     = block->getData< ScalarField_T >(scalarDataId_);

         CellInterval xyz = exactField->xyzSize();

         Block& b           = dynamic_cast< Block& >(*block);
         uint_t level       = b.getLevel();

         const real_t & dx = blocks_->dx( level );
         const real_t & dy = blocks_->dy( level );
         const real_t & dz = blocks_->dz( level );

         const real_t volume = dx * dy * dx;
         
         for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
         { 
            const real_t exact   = exactField -> get(*cellIt);
            const real_t compute = cField -> get(*cellIt);

            er0  = std::max(er0, std::abs(compute - exact));
            er1 += std::abs(compute - exact) * volume;
            er2 += std::pow(std::abs(compute - exact), 2)  * volume;
         }
      }

      WALBERLA_MPI_WORLD_BARRIER();
      mpi::allReduceInplace( er0, mpi::MAX );
      mpi::allReduceInplace( er1, mpi::SUM );
      mpi::allReduceInplace( er2, mpi::SUM );

      error_norms.push_back( Vector3<real_t> (er0, er1, er2 ) );
   }

   std::vector< Vector3<real_t> > & getErrorNorms( ) 
   {
      return error_norms;
   }

   void writeText ( const real_t & CFL, const real_t & scaling = 1 )
   {
      const CellInterval DomainBB = blocks_ ->getDomainCellBB( uint_t(0) );
      WALBERLA_ASSERT_EQUAL(DomainBB.xMax(), DomainBB.yMax(),  DomainBB.zMax())

      std::ofstream f("LBM_Advection_ic-" + std::to_string(ic_) + "_" + std::to_string(DomainBB.xMax() + 1) +"x"+ std::to_string(DomainBB.yMax() + 1)+"x"+ std::to_string(DomainBB.zMax() + 1) + ".tex");

      if (!f.is_open())
        WALBERLA_ABORT("Error opening file for writing.")  

      f << "CFL   = " << CFL << "\n";
      f << "Lamda X = " << lam[0] * scaling << "\n";
      f << "Lamda Y = " << lam[1] * scaling << "\n";
      f << "Lamda Z = " << lam[2] * scaling << "\n";

      f << std::string(100, '_') << "\n";
      f << "i       | L1         | L2         | L∞         \n";
      f << std::string(100, '_') << "\n";

      std::vector<real_t> L0_errors{};
      std::vector<real_t> L1_errors{};
      std::vector<real_t> L2_errors{};

      for (const auto & vec : error_norms) {
           L0_errors.push_back(vec[0]);
           L1_errors.push_back(vec[1]);
           L2_errors.push_back(vec[2]);
      }

      WALBERLA_ASSERT_EQUAL(L0_errors.size(), L1_errors.size(),  L2_errors.size())

      for (uint_t i = 0; i <= L0_errors.size() - 1; ++i) 
      {
         const real_t L0 = L0_errors[i];
         const real_t L1 = L1_errors[i];
         const real_t L2 = L2_errors[i];
   
         f << std::setw(7) << i + 1 << " | "
          << std::scientific << std::setprecision(4) << L1 << " | "
          << std::scientific << std::setprecision(4) << L2 << " | "
          << std::scientific << std::setprecision(4) << L0 << "\n";

      }
      
      f << std::string(100, '_') << "\n";

   }

private:

   shared_ptr< StructuredBlockForest > blocks_;
   const BlockDataID exactDataId_;
   const BlockDataID scalarDataId_;
   const real_t kappa_;
   const real_t h0_;
   const uint_t ic_;
   const Vector3< real_t > lam;

   std::vector< Vector3<real_t> > error_norms;
   
   real_t Average( const Vector3< real_t > & pmin, const Vector3< real_t > & pmax, const real_t & dx, const real_t & dy, const real_t & dz, const real_t & time = 0 )
   {
      if ( ic_ == 3)
      {
         return - (-std::cos(M_PI*pmax[0]) + std::cos(M_PI*pmin[0]))*\
                (-std::cos(M_PI*pmax[1]) + std::cos(M_PI*pmin[1]))*\
                (-std::cos(M_PI*pmax[2]) + std::cos(M_PI*pmin[2]))/\
                std::pow(M_PI , 3) /(-dx*dy*dz);
      }
      else if ( ic_ == 4)
      {
         return (  (2*std::pow(std::sin(M_PI*pmax[0]), 3)*std::cos(M_PI*pmax[0])+3*std::cos(M_PI*pmax[0])*std::sin(M_PI*pmax[0])-3*M_PI*pmax[0]-2*std::pow(std::sin(M_PI*pmin[0]), 3)*std::cos(M_PI*pmin[0])-3*std::cos(M_PI*pmin[0])*std::sin(M_PI*pmin[0])+3*M_PI*pmin[0]) \
                 *(2*std::pow(std::sin(M_PI*pmax[1]), 3)*std::cos(M_PI*pmax[1])+3*std::cos(M_PI*pmax[1])*std::sin(M_PI*pmax[1])-3*M_PI*pmax[1]-2*std::pow(std::sin(M_PI*pmin[1]), 3)*std::cos(M_PI*pmin[1])-3*std::cos(M_PI*pmin[1])*std::sin(M_PI*pmin[1])+3*M_PI*pmin[1]) \
                 *(2*std::pow(std::sin(M_PI*pmax[2]), 3)*std::cos(M_PI*pmax[2])+3*std::cos(M_PI*pmax[2])*std::sin(M_PI*pmax[2])-3*M_PI*pmax[2]-2*std::pow(std::sin(M_PI*pmin[2]), 3)*std::cos(M_PI*pmin[2])-3*std::cos(M_PI*pmin[2])*std::sin(M_PI*pmin[2])+3*M_PI*pmin[2])) \
                /std::pow(M_PI , 3) /(-dx*dy*dz) / 512;
      }
      else
      {
         // Sixth order Gaussian rule
         const real_t s  = std::sqrt(3.0) / std::sqrt(5.0);
         const real_t w1 = 5.0 / 9.0;
         const real_t w2 = 8.0 / 9.0;

         const Vector3<real_t> cc= 0.5 * (pmin + pmax);

         real_t z = cc[2] - dz * s / 2;
         const real_t INT1 =  w1 * (w1 * V0(cc[0] - dx * s / 2, cc[1] - dy * s / 2, z, time) + w1 * V0(cc[0] + dx * s / 2, cc[1] - dy * s / 2, z, time) +
                                   w2 * V0(cc[0], cc[1] - dy * s / 2, z, time)) + \
                             w1 * (w1 * V0(cc[0] - dx * s / 2, cc[1] + dy * s / 2, z, time) + w1 * V0(cc[0] + dx * s / 2, cc[1] + dy * s / 2, z, time) +
                                   w2 * V0(cc[0], cc[1] + dy * s / 2, z, time)) + \
                             w2 * (w1 * V0(cc[0] - dx * s / 2, cc[1], z, time) + w1 * V0(cc[0] + dx * s / 2, cc[1], z, time) + w2 * V0(cc[0], cc[1], z, time));

         z = cc[2] + dz * s / 2;
         const real_t INT2 =  w1 * (w1 * V0(cc[0] - dx * s / 2, cc[1] - dy * s / 2, z, time) + w1 * V0(cc[0] + dx * s / 2, cc[1] - dy * s / 2, z, time) +
                                   w2 * V0(cc[0], cc[1] - dy * s / 2, z, time)) + \
                             w1 * (w1 * V0(cc[0] - dx * s / 2, cc[1] + dy * s / 2, z, time) + w1 * V0(cc[0] + dx * s / 2, cc[1] + dy * s / 2, z, time) +
                                   w2 * V0(cc[0], cc[1] + dy * s / 2, z, time)) + \
                             w2 * (w1 * V0(cc[0] - dx * s / 2, cc[1], z, time) + w1 * V0(cc[0] + dx * s / 2, cc[1], z, time) + w2 * V0(cc[0], cc[1], z, time));

         z = cc[2];
         const real_t INT3 =  w1 * (w1 * V0(cc[0] - dx * s / 2, cc[1] - dy * s / 2, z, time) + w1 * V0(cc[0] + dx * s / 2, cc[1] - dy * s / 2, z, time) +
                                   w2 * V0(cc[0], cc[1] - dy * s / 2, z, time)) + \
                             w1 * (w1 * V0(cc[0] - dx * s / 2, cc[1] + dy * s / 2, z, time) + w1 * V0(cc[0] + dx * s / 2, cc[1] + dy * s / 2, z, time) +
                                   w2 * V0(cc[0], cc[1] + dy * s / 2, z, time)) + \
                             w2 * (w1 * V0(cc[0] - dx * s / 2, cc[1], z, time) + w1 * V0(cc[0] + dx * s / 2, cc[1], z, time) + w2 * V0(cc[0], cc[1], z, time));

         return (w1 * INT1 + w1 * INT2 + w2 * INT3) / 8;
      }
   }

   real_t V0(real_t x, real_t y, real_t z, real_t time)
   {
      const real_t per = 2.0;

      // Repeat until x is within the range [-1, 1]
      while (x < -1.0) {
         x += per;
      }
      while (x > 1.0) {
         x -= per;
      }

      // Repeat until y is within the range [-1, 1]
      while (y < -1.0) {
         y += per;
      }
      while (y > 1.0) {
         y -= per;
      }

      // Repeat until z is within the range [-1, 1]
      while (z < -1.0) {
         z += per;
      }
      while (z > 1.0) {
         z -= per;
      }

      // Compute V0 based on the selected case (ic)
      switch (ic_) {
      case 1:
      {
         real_t s = x * x + y * y + z * z;
         s = std::exp(-10 * s);
         return std::max(s, 1e-15);
      }
      case 2:
      {
         if (std::abs(x) <= 0.5 && std::abs(y) <= 0.5 && std::abs(z) <= 0.5) {
            return 1.0;
         } else {
            return 0.0;
         }
      }
      case 3:
      {
         return std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z);
      }
      case 4:
      {
         return std::pow(std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z), 4);
      }
      case 5:
      {
         real_t r = std::sqrt(x * x + y * y + z * z);
         if (r <= 0.5) {
            return 1.0;
         } else {
            return 0.0;
         }
      }
      case 6:
      {
         real_t s = x * x + y * y + z * z;
         real_t sigma_d = std::sqrt(2 * kappa_ * time);
         real_t sigma_0 = 1/(2*std::sqrt(5));
         s = std::pow((sigma_0*sigma_0)/((sigma_0*sigma_0) + (sigma_d* sigma_d)), 1.5) * std::exp(- s / (2 * ((sigma_0*sigma_0) + (sigma_d* sigma_d))));
         return std::max(s, 1e-15);
      }
      case 7:
      {
         real_t erf_x = 0.5 * ( std::erf( (x) / std::sqrt( 4 * kappa_ * time + (h0_*h0_) )));
         real_t erf_y = 0.5 * ( std::erf( (y) / std::sqrt( 4 * kappa_ * time + (h0_*h0_) )));
         real_t erf_z = 0.5 * ( std::erf( (z) / std::sqrt( 4 * kappa_ * time + (h0_*h0_) )));

         return erf_x * erf_y * erf_z;
      }
      default:
      {
         WALBERLA_LOG_WARNING("Incompatible Initial Condition Setting to Zero!!")
         return 0.0;
      }

      }
   }

};

void setupBoundaryFlagField(StructuredBlockForest& sbfs, const BlockDataID flagFieldID, const uint_t gl)
{
   for (auto bIt = sbfs.begin(); bIt != sbfs.end(); ++bIt)
   {
      Block& b           = dynamic_cast< Block& >(*bIt);
      uint_t level       = b.getLevel();
      auto flagField     = b.getData< FlagField_T >(flagFieldID);

      uint8_t neumannFlag   = uint8_c(flagField->getOrRegisterFlag(neumannFlagUID));

      for (auto cIt = flagField->beginWithGhostLayerXYZ( cell_idx_c(gl) ); cIt != flagField->end(); ++cIt)
      {
         Cell localCell = cIt.cell();
         Cell globalCell(localCell);
         sbfs.transformBlockLocalToGlobalCell(globalCell, b);

         if       ( globalCell.y() < 0 ) { flagField->addFlag(localCell, neumannFlag); }
         else if  ( globalCell.y() >= cell_idx_c(sbfs.getNumberOfYCells(level)) ){ flagField->addFlag(localCell, neumannFlag); } 
         else if  ( globalCell.z() < 0 ) { flagField->addFlag(localCell, neumannFlag); } 
         else if  ( globalCell.z() >= cell_idx_c(sbfs.getNumberOfYCells(level))) { flagField->addFlag(localCell, neumannFlag); }       
      }
   }
}

int main(int argc, char** argv)
{
   Environment env(argc, argv);
   if (!env.config()) { WALBERLA_ABORT("No configuration file specified!"); }

   mpi::MPIManager::instance()->useWorldComm();

   Config::BlockHandle simulationParamsBlock = env.config()->getOneBlock( "SimulationParameters" );
   const real_t         physical_velocity    = simulationParamsBlock.getParameter< real_t      >( "velMag",                         real_t( 1.0 )  );
   const Vector3<real_t>  VelocityComponents = simulationParamsBlock.getParameter<Vector3<real_t> > ( "velocityComponents",   Vector3<real_t> (0)  );
   const real_t            uLB               = simulationParamsBlock.getParameter< real_t      >( "uLB",   real_t(0.1)    ); 
   const uint_t         advectionPeriods     = simulationParamsBlock.getParameter< uint_t      >( "advectionPeriods"                               );
   const real_t remainingTimeLoggerFrequency = simulationParamsBlock.getParameter< real_t      >( "remainingTimeLoggerFrequency",   real_t(5.0)    ); // in seconds
   const uint_t          ic                  = simulationParamsBlock.getParameter< uint_t      >( "ic"  );

   Config::BlockHandle scalarConfigBlock = env.config()->getOneBlock( "ScalarTransport" );
   const real_t  kappa                   = scalarConfigBlock.getParameter< real_t  >( "kappa",   real_t( 0.0 )  );

   //------------------------------------------------------------------------------------------------------//
   Config::BlockHandle domainBlock = env.config()->getOneBlock( "DomainSetup" );
   const real_t          dx        = domainBlock.getParameter< real_t      >( "dx"); // in seconds
   
   const Vector3<real_t> PhysicalInitVel = VelocityComponents * physical_velocity;
   const Vector3<real_t> InitVel = PhysicalInitVel * uLB;

   //------------------------------------------------------------------------------------------------------//

   
   const real_t dt = (uLB / physical_velocity) * dx;
   const real_t lattice_viscosity = kappa * (dt / (dx * dx));
   const real_t omega = lbm::collision_model::omegaFromViscosity(lattice_viscosity);

   const real_t CFL = ( InitVel[0] + InitVel[1] + InitVel[2] ) ;
   WALBERLA_LOG_INFO_ON_ROOT("dt/dxxx" << (kappa*dt)/(dx*dx));
   WALBERLA_LOG_INFO_ON_ROOT("Lattice viscosity in LBM " << lattice_viscosity );
   
   if (CFL > 1)
      { WALBERLA_ABORT("CFL criterion violated"); }

   const real_t singleAdvectionTime = (VelocityComponents * 2).length() / (PhysicalInitVel.length());
   const real_t simulationTime = real_c(advectionPeriods) * singleAdvectionTime;    // Domain Interval [-1,1]
   
   const uint_t timesteps   = uint_t(round(simulationTime/dt)) + 1;

   std::ostringstream oss;
   oss <<   "- Physical Parameters:"
         << "\n   + dx:               " << dx
         << "\n   + dt:               " << dt
         << "\n   + CFL:              " << CFL
         << "\n   + omega:            " << omega;

   oss <<"\n- Time Parameters:"
      << "\n   + Course Timestep (s):       " << dt
      << "\n   + Simulation Time (s) [actual/specified]:  [" << real_t(timesteps)*dt <<"/"<<simulationTime<<"]"
      << "\n   + No. Course Timesteps:      " << timesteps;
   std::string str = oss.str();

   WALBERLA_LOG_INFO_ON_ROOT( "Physical/LBM Parameter Conversion:\n" << str);

   //------------------------------------------------------------------------------------------------------//
   WALBERLA_LOG_INFO_ON_ROOT("Creating block forest...")

   // domain creation
   std::shared_ptr<StructuredBlockForest> blocks;
   {
      const  Vector3< uint_t >  reqNumBlocks  = domainBlock.getParameter< Vector3< uint_t > >( "blocks"); 
      const  Vector3< uint_t >  reqCellsPerBlock  = domainBlock.getParameter< Vector3< uint_t > >( "cellsPerBlock"); 
      const  Vector3< uint_t >  reqDomainSize {reqNumBlocks[0]*reqCellsPerBlock[0], reqNumBlocks[1]*reqCellsPerBlock[1], reqNumBlocks[2]*reqCellsPerBlock[2]};

      Vector3< uint_t > numBlocks;
      Vector3< uint_t > cellsPerBlock;
      blockforest::calculateCellDistribution(reqDomainSize,
                                             uint_c(mpi::MPIManager::instance()->numProcesses()),
                                             numBlocks, cellsPerBlock);

      WALBERLA_LOG_INFO_ON_ROOT(numBlocks << " , " << cellsPerBlock)
      const  Vector3< uint_t >  periodicity  = domainBlock.getParameter< Vector3< uint_t > >( "periodic");
      const  Vector3< real_t >  domainSize { real_c(reqDomainSize[0]) * dx, real_c(reqDomainSize[1]) * dx, real_c(reqDomainSize[2]) * dx} ;

      SetupBlockForest sforest;

      sforest.addWorkloadMemorySUIDAssignmentFunction( blockforest::uniformWorkloadAndMemoryAssignment );

      sforest.init( AABB(0_r, 0_r, 0_r, real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2])),
                     numBlocks[0], numBlocks[1], numBlocks[2], periodicity[0], periodicity[1], periodicity[2] );

      // calculate process distribution

      const memory_t memoryLimit = numeric_cast< memory_t >( sforest.getNumberOfBlocks() );

      const blockforest::GlobalLoadBalancing::MetisConfiguration< SetupBlock > metisConfig(
         true, false, std::bind( blockforest::cellWeightedCommunicationCost, std::placeholders::_1, std::placeholders::_2,
                                 cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2] ) );

      sforest.calculateProcessDistribution_Default( uint_c( MPIManager::instance()->numProcesses() ), memoryLimit,
                                                      "hilbert", 10, false, metisConfig );

      if( !MPIManager::instance()->rankValid() )
         MPIManager::instance()->useWorldComm();

      // create StructuredBlockForest (encapsulates a newly created BlockForest)

      WALBERLA_LOG_INFO_ON_ROOT("SetupBlockForest created successfully:\n" << sforest)

      sforest.writeVTKOutput("domain_decomposition");

      auto bf = std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, false );

      blocks = std::make_shared< StructuredBlockForest >( bf, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2] );
      blocks->createCellBoundingBoxes();

   }

   //------------------------------------------------------------------------------------------------------//
   // LBM Fields
   const StorageSpecification_T StorageSpec = StorageSpecification_T();
   const BlockDataID pdfFieldId  = lbm_generated::addPdfFieldToStorage(blocks, "pdfs", StorageSpec, FieldGhostLayer, field::fzyx);
   const BlockDataID uFieldId    = field::addToStorage< VelocityField_T >(blocks, "u" , real_t(0.0), field::fzyx, FieldGhostLayer);
   const BlockDataID cFieldId    = field::addToStorage< ScalarField_T >(blocks, "c", real_t(0.0), field::fzyx, FieldGhostLayer);
   BlockDataID exactFieldId      = field::addToStorage< ExactField_T >(blocks, "exact", real_t(0.0), field::fzyx, uint_t(1));

   const BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", FieldGhostLayer);

   //setupBoundaryFlagField(*blocks, flagFieldId, FieldGhostLayer);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, fluidFlagUID);

   BoundaryCollection_T boundaryCollection(blocks, flagFieldId, pdfFieldId, fluidFlagUID);

   // Set velocity field for pdf setter
   ScalarHelperFunctions setup( blocks, exactFieldId, cFieldId, kappa, real_t(0.02), ic, PhysicalInitVel );
   setup.initDirichletBoundary();
   setup.initVelField(  uFieldId,  InitVel );



   SweepCollection_T sweepCollection(blocks, pdfFieldId, uFieldId, cFieldId, omega);
   for (auto& block : *blocks)
   {
      sweepCollection.initialise(&block, uint_t(1));
   }

   pystencils::SIConverter convert2SI( cFieldId, real_t(-1) );
   //------------------------------------------------------------------------------------------------------//
   
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication(blocks);
   communication.addPackInfo(std::make_shared<lbm_generated::UniformGeneratedPdfPackInfo< PdfField_T >>(pdfFieldId)); 

   //------------------------------------------------------------------------------------------------------//

   SweepTimeloop timeLoop(blocks, timesteps );

   // LBM VTK Output
   const uint_t vtkWriteFrequency = simulationParamsBlock.getParameter< uint_t >("vtkWriteFrequency", 0);
   if (vtkWriteFrequency > 0)
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "FluidVTK_dx-" + std::to_string( dx ), vtkWriteFrequency, uint_t(0), false,
                                                      "vtk_out", "simulation_step", false, true, true, false, 0);

      auto cWriter   = make_shared< field::VTKWriter< ScalarField_T > >(cFieldId, "Concentration Density");
      
      vtkOutput->addCellDataWriter(cWriter);

      vtkOutput->addBeforeFunction([&]() {
         for (auto& block : *blocks){
            sweepCollection.calculateMacroscopicParameters(&block);
            convert2SI(&block);
            }
      });

      timeLoop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }


   timeLoop.add() << BeforeFunction(communication, "communication")
                  << Sweep(boundaryCollection.getSweep(BoundaryCollection_T::ALL), "Boundary Conditions");
   timeLoop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL), "LBM StreamCollide");

   
   timeLoop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeLoop.getNrOfTimeSteps() , remainingTimeLoggerFrequency), "remaining time logger" );

   auto calculte_error_norms = [&]()
   { 
      real_t time = real_t( timeLoop.getCurrentTimeStep() + 1 ) * dt;

      if ( ( timeLoop.getCurrentTimeStep() + 1 ) % uint_t(singleAdvectionTime/dt) == 0 )
      {
         for (auto& block : *blocks){
            sweepCollection.calculateMacroscopicParameters(&block);
         }
         setup.calculateErrorNorms( time );
      }
         
   } ;

   timeLoop.addFuncAfterTimeStep( calculte_error_norms , "Error Norms"); 


   //------------------------------------------------------------------------------------------------------//
   //-----------------------------------------------RUN SIMULATION-----------------------------------------//
   //------------------------------------------------------------------------------------------------------//
   WcTimingPool timingPool;
   WALBERLA_LOG_INFO_ON_ROOT( "Starting timeloop" );
   timeLoop.run();
   WALBERLA_LOG_INFO_ON_ROOT( "Timeloop done" );
   timingPool.unifyRegisteredTimersAcrossProcesses();
   timingPool.logResultOnRoot( timing::REDUCE_TOTAL, true );

   setup.writeText( CFL );

   return EXIT_SUCCESS;

   }
}

int main(int argc, char** argv) { walberla::main(argc, argv); }
