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
//! \file SuViscoelasticityTest.cpp
//! \ingroup lbm
//! \author Cameron Stewart <cstewart@icp.uni-stuttgart.de>
//! \brief D2Q9 Poisuielle flow simulation with viscoelasticity
//
//
//======================================================================================================================


#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/BoundaryHandling.h"

#include "core/Abort.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/SharedFunctor.h"
#include "core/Filesystem.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/communication/PackInfo.h"
#include "field/StabilityChecker.h"

#include "geometry/initializer/BoundaryFromDomainBorder.h"
#include "geometry/InitBoundaryHandling.h"

#include "lbm/lattice_model/D2Q9.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/lattice_model/ForceModel.h"
#include "lbm/SuViscoelasticity.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/vtk/all.h"

#include "timeloop/SweepTimeloop.h"

#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <lbm/boundary/factories/ExtendedBoundaryHandlingFactory.h>

using namespace walberla;

//////////////
// TYPEDEFS //
//////////////

typedef GhostLayerField< Vector3<real_t>, 1> ForceField_T;
typedef lbm::D2Q9< lbm::collision_model::TRT, false, lbm::force_model::GuoField< ForceField_T > >  LatticeModel_T;
typedef LatticeModel_T::CommunicationStencil    CommunicationStencil_T;
typedef GhostLayerField<Matrix3<real_t>, 1> StressField_T;
typedef GhostLayerField< Vector3<real_t>, 1> VelocityField_T;
typedef lbm::PdfField< LatticeModel_T >  PdfField_T;
typedef walberla::uint8_t flag_t;
typedef FlagField< flag_t > FlagField_T;
typedef lbm::ExtendedBoundaryHandlingFactory< LatticeModel_T, FlagField_T > BHFactory;

class WriteCheckpoint {
public:
   WriteCheckpoint(Timeloop &timeloop, BlockStorage & blockStorage, uint_t checkpointFrequency, uint_t timesteps,
                   BlockDataID pdfFieldId, BlockDataID stressId, std::string dir, std::string paramFile) : timeloop_(timeloop), blockStorage_(blockStorage),
                                                                     checkpointFrequency_(checkpointFrequency), timesteps_(timesteps),
                                                                     pdfFieldId_(pdfFieldId), stressId_(stressId), dir_(dir + "_chk"), paramFile_(paramFile) {}
   void operator()(){
      if(checkpointFrequency_ > 0 && timeloop_.getCurrentTimeStep() % checkpointFrequency_== 0 && timeloop_.getCurrentTimeStep() > 0){
         Save();
      }
      if(timeloop_.getCurrentTimeStep() == timesteps_-1 && filesystem::is_directory(dir_)){
         filesystem::remove_all(dir_);
      }
   }

   bool canLoad()
   {
      if(!filesystem::is_directory(dir_)){
         std::cout << "No Checkpoint Found" << std::endl;
         return false;
      }

      std::ifstream f1;
      std::ifstream f2;
      f1.open(paramFile_, std::ifstream::binary|std::ifstream::ate);
      f2.open(dir_ + "/" + paramFile_, std::ifstream::binary|std::ifstream::ate);

      if (f1.fail() || f2.fail()) {
         std::cout << "Failed to check parameter file" << std::endl;
         return false; //file problem
      }

      if (f1.tellg() != f2.tellg()) {
         std::cout << "Ignoring checkpoint for different size config file" << std::endl;
         return false; //size mismatch
      }

      f1.seekg(0, std::ifstream::beg);
      f2.seekg(0, std::ifstream::beg);
      if(!std::equal(std::istreambuf_iterator<char>(f1.rdbuf()), std::istreambuf_iterator<char>(), std::istreambuf_iterator<char>(f2.rdbuf()))){
         std::cout << "Ignoring checkpoint for different config file" << std::endl;
         return false;
      }
      return true;
   }

   void Load()
   {
      std::cout << "Loading Checkpoint" << std::endl;
      field::readFromFile<PdfField_T>(dir_ + "/pdfField.chk", blockStorage_, pdfFieldId_);
      field::readFromFile<StressField_T>(dir_ + "/stressField.chk", blockStorage_, stressId_);
      std::ifstream f;
      f.open(dir_ + "/chkpoint");
      uint_t currentTime;
      f >> currentTime;
      timeloop_.setCurrentTimeStep(currentTime + 1);
   }

   void Save(){

      auto tmpdir = dir_ + "~";
      if(filesystem::is_directory(tmpdir)){
         filesystem::remove_all(tmpdir);
      }

      filesystem::create_directory(tmpdir);

      filesystem::copy(paramFile_, tmpdir + "/" + paramFile_);
      field::writeToFile<PdfField_T>(tmpdir + "/pdfField.chk", blockStorage_, pdfFieldId_);
      field::writeToFile<StressField_T>(tmpdir + "/stressField.chk", blockStorage_, stressId_);
      std::ofstream f;
      f.open(tmpdir + "/chkpoint");
      f << timeloop_.getCurrentTimeStep();
      f.close();

      if(filesystem::is_directory(dir_)) {
         filesystem::remove_all(dir_);
      }
      filesystem::rename(tmpdir, dir_);
   }

private:
   Timeloop & timeloop_;
   BlockStorage & blockStorage_;
   uint_t checkpointFrequency_, timesteps_;
   BlockDataID pdfFieldId_, stressId_;
   std::string dir_, paramFile_;
};

class WriteVelocity
{
public:
   WriteVelocity(Timeloop & timeloop, uint_t writeFrequency, shared_ptr< StructuredBlockForest > blocks, BlockDataID pdfFieldId) : timeloop_(timeloop), writeFrequency_(writeFrequency), blocks_(blocks),
   pdfFieldId_(pdfFieldId){}

   void operator()(){
      if(writeFrequency_ > 0 && timeloop_.getCurrentTimeStep() % writeFrequency_==0 ){

         std::ofstream f;
         f.open("vel.dat", std::ios_base::app);
         for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
         {
            PdfField_T * pdf = blockIt.get()->getData< PdfField_T >(pdfFieldId_);

            if(blockIt.get()->getAABB().contains(15,15,0)){
               std::cout<< pdf->getVelocity(5,5,0).length() << std::endl;
               f << std::setprecision(10) << pdf->getVelocity(5,5,0).length() << ",";
            }

         }
         f.close();
      }
   }

private:
   Timeloop & timeloop_;
   uint_t writeFrequency_;
   shared_ptr< StructuredBlockForest > blocks_;
   BlockDataID pdfFieldId_;
};



template < typename BoundaryHandling_T >
class ConstantForce
{
public:

   ConstantForce( BlockDataID forceFieldId, BlockDataID boundaryHandlingId, real_t force)
         : forceFieldId_( forceFieldId ), boundaryHandlingId_(boundaryHandlingId), force_(force)
   {}

   void operator()( IBlock * block )
   {
      ForceField_T *forceField = block->getData< ForceField_T >(forceFieldId_);
      BoundaryHandling_T *boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingId_ );

      WALBERLA_FOR_ALL_CELLS_XYZ(forceField,
                                 {
                                    Cell cell(x,y,z);
                                    if (boundaryHandling->isDomain(cell)) {
                                       forceField->get(cell)[0] += force_;
                                    }
                                 })
   }

private:

   BlockDataID forceFieldId_, boundaryHandlingId_;
   real_t force_;
};

int main(int argc, char ** argv ){

   walberla::Environment env( argc, argv );

   // create block storage from the config file
   shared_ptr<StructuredBlockForest> blocks = blockforest::createUniformBlockGridFromConfig( env.config() );

   // read parameters
   auto parameters = env.config()->getOneBlock( "Parameters" );

   // extract some constants from the parameters
   const std::string paramFile(argv[1]);
   const real_t          eta_s           = parameters.getParameter< real_t >         ( "eta_s",           real_c( 1.4 ) );
   const real_t          force           = parameters.getParameter< real_t >         ( "force",           real_t( 1 )   );
   const real_t          eta_p           = parameters.getParameter< real_t >         ( "eta_p",           real_c( 1.0 ) );
   const uint_t          period          = parameters.getParameter< uint_t >         ( "period",          uint_c( 1   ) );
   const uint_t          writePeriod     = parameters.getParameter< uint_t >         ( "writePeriod",     uint_c( 10  ) );

   const uint_t checkpointFrequency = parameters.getParameter< uint_t > ("checkpointFrequency", uint_c(5000));
   const real_t remainingTimeLoggerFrequency = parameters.getParameter< real_t >( "remainingTimeLoggerFrequency", real_c(3.0) ); // in seconds

   // take parameters from command line args preferably
   char* end;
   const real_t lambda_p = ((argc > 2)?        atof(argv[2])              : parameters.getParameter<real_t>( "lambda_p", real_c(10)));
   const std::string baseFolder = ((argc > 2)? std::string(argv[2])       : parameters.getParameter<std::string>( "baseFolder", std::string("vtk_out")));
   const uint_t timesteps = ((argc > 3)?       strtoul(argv[3],&end,10)   : parameters.getParameter<uint_t>( "timesteps", uint_c(10) ));

   // create fields
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", uint_c(2));
   BlockDataID forceFieldId = field::addToStorage<ForceField_T>( blocks, "Force Field", Vector3<real_t>(0.0), field::zyxf, uint_c(2));
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::TRT::constructWithMagicNumber( walberla::lbm::collision_model::omegaFromViscosity(eta_s)), lbm::force_model::GuoField<ForceField_T>( forceFieldId ) );
   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel, Vector3<real_t>(), real_c(1), uint_c(2) );
   BlockDataID stressId = walberla::field::addToStorage<StressField_T>( blocks, "Stress Field", Matrix3<real_t>(0.0), field::zyxf, uint_c(2));
   BlockDataID stressOldId = walberla::field::addToStorage<StressField_T>( blocks, "Old Stress Field", Matrix3<real_t>(0.0), field::zyxf, uint_c(2));
   BlockDataID velocityId = walberla::field::addToStorage<VelocityField_T> (blocks, "Velocity Field", Vector3<real_t>(0.0), field::zyxf, uint_c(1));

   // create and initialize boundary handling
   const FlagUID fluidFlagUID( "Fluid" );

   auto boundariesConfig = env.config()->getOneBlock( "Boundaries" );

   BlockDataID boundaryHandlingId = BHFactory::addBoundaryHandlingToStorage( blocks, "boundary handling", flagFieldId, pdfFieldId, fluidFlagUID );
   geometry::initBoundaryHandling<BHFactory::BoundaryHandling>( *blocks, boundaryHandlingId, boundariesConfig );
   geometry::setNonBoundaryCellsToDomain<BHFactory::BoundaryHandling> ( *blocks, boundaryHandlingId );


   // create time loop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication( blocks );
   communication.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldId ) );

   timeloop.add() << BeforeFunction( communication, "communication" )
                  << Sweep( BHFactory::BoundaryHandling::getBlockSweep( boundaryHandlingId ), "boundary handling" );
   timeloop.add() << BeforeFunction( lbm::viscoelastic::Su<LatticeModel_T, BHFactory::BoundaryHandling>(blocks, forceFieldId, pdfFieldId, boundaryHandlingId, stressId, stressOldId, velocityId,
                                                                                               lambda_p, eta_p, period, true), "viscoelasticity")
                  << Sweep( ConstantForce<BHFactory::BoundaryHandling>(forceFieldId, boundaryHandlingId, force),"Poiseuille Force");
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, fluidFlagUID ) ),
                            "LB stream & collide" );

   // LBM stability check
   timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< ForceField_T, FlagField_T >( env.config(), blocks, forceFieldId,
                                                                                                             flagFieldId, fluidFlagUID ) ),
                                  "LBM stability check" );

   // log remaining time
   timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "remaining time logger" );

   // Add checkpointing
   //auto checkpointing = WriteCheckpoint(timeloop, blocks->getBlockStorage(), checkpointFrequency, timesteps, pdfFieldId, stressId, baseFolder, paramFile);
   //timeloop.addFuncAfterTimeStep(checkpointing);
   //if(checkpointing.canLoad()){checkpointing.Load();}

   // write midvelocity
   //auto writeVelocity = WriteVelocity(timeloop, writePeriod, blocks, pdfFieldId);
   //timeloop.addFuncAfterTimeStep(writeVelocity);

   // add VTK outputs
   timeloop.addFuncAfterTimeStep(field::createVTKOutput<StressField_T>(stressId, *blocks, "stress", writePeriod, uint_c(0), false, baseFolder, std::string("step"),
                                                                       false, true, true, 0, Set<SUID>::emptySet(), Set<SUID>::emptySet(), true, timeloop.getCurrentTimeStep()));
   timeloop.addFuncAfterTimeStep(field::createVTKOutput<VelocityField_T>(velocityId, *blocks, "velocity", writePeriod, uint_c(0), false, baseFolder, std::string("step"),
                                                                       false, true, true, 0, Set<SUID>::emptySet(), Set<SUID>::emptySet(), true, timeloop.getCurrentTimeStep()));

   auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", writePeriod, 0, false, baseFolder );
   flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldId, "FlagField" ) );
   timeloop.addFuncAfterTimeStep(vtk::writeFiles( flagFieldVTK ), "Flag field VTK");

   timeloop.run();


   return EXIT_SUCCESS;
}
