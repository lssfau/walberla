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
//! \file Poiseuille.cpp
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Calculates Flow in a cylindrical pipe - checks against analytical formula
//!
//! Configuration:
//! \verbatim
//!    Poiseuille {
//!       geometry         box;    box or pipe
//!       diameter         0.01;   radius of pipe in m, or half of box width
//!       aspectRatio      5;      pipeLength = aspectRatio  diameter
//!       cellsPerDiameter 60;     dx = diameter / cellsPerDiameter
//!       viscosity        1e-6;   dynamic viscosity in m^2/s
//!       omega            0.9;    LBM relaxation parameter
//!       pressureBoundary 0;      if true use pressure boundary, otherwise use gravity driven flow
//!       useGui;                  enables graphical user interface
//!    }
//! \endverbatim
//! "rect2D" Setup:
//!    - periodic in y and z coordinate
//!    - no slip at borders of x coordinate
//!    - only a velocity in z direction called $w(x)$
//!    - velocities in x,y are zero
//!    - stationary -> all time derivatives zero
//!    - analytical formula, derived from third momentum equation:
//!       $w(x) = - \frac{\rho g_z - \frac{\partial p}{\partial z} { 2 \nu } ( R^2 - x^2)$
//!       where $R$ is half of the channel width
//!    - reduced for gravity driven: w(x) = - \frac{\rho g_z } { 2 \nu } ( R^2 - x^2)
//!    - pressure difference driven, use $\Delta P = g_z \cdot \rho \cdot \Delta Z $
//! "Pipe" Setup:
//!    - periodic in z direction, cylindrical geometry with flow direction in z
//!      and noslip at pipe boundary
//!    - analytical formulas similar but different factor (4 instead of 2):
//!      $ w(x) = - \frac{\rho g_z } { 4 \nu } ( R^2 - x^2) $
//
//======================================================================================================================

#include "lbm/all.h"

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/BoundaryHandling.h"

#include "core/Abort.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/IntegerFactorization.h"
#include "core/timing/RemainingTimeLogger.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/adaptors/ComponentExtractionAdaptor.h"
#include "field/communication/PackInfo.h"

#include "gather/CellGatherPackInfo.h"
#include "gather/CurveGatherPackInfo.h"
#include "gather/GnuPlotGraphWriter.h"
#include "gather/MPIGatherScheme.h"

#include "geometry/initializer/BoundaryFromDomainBorder.h"
#include "geometry/InitBoundaryHandling.h"

#include "gui/Gui.h"

#include "timeloop/SweepTimeloop.h"

#include <iostream>
#include <vector>


using namespace walberla;


typedef walberla::uint8_t                   flag_t;

// Compressible SRT with constant force
using lbm::collision_model::SRT;
using lbm::force_model::GuoConstant;
typedef lbm::D3Q19< SRT, true, GuoConstant> LM;

// Fields
typedef lbm::PdfField<LM>                   PdfField;
typedef FlagField<flag_t>                   FField;
typedef GhostLayerField<Vector3<real_t>,1 > VectorField;
typedef GhostLayerField<real_t,1 >          ScalarField;

/*
 * TODO evaluation is currently wrong
 */
/*
class PoiseuilleVelocityDataProcessor : public gather::DataProcessor
{
public:
   PoiseuilleVelocityDataProcessor( real_t acceleration_l, real_t viscosity_l, cell_idx_t cellsPerRadius,
                                    bool pipeOrBox,
                                    const std::string & filename = "" )
      : acceleration_l_ ( acceleration_l ), viscosity_l_ ( viscosity_l ),
        pipeOrBox_(pipeOrBox), cellsPerRadius_( real_c( cellsPerRadius) ),
        lastMeanError_( 0 ), lastMaxError_( 0 )
   {
      if ( filename.size() > 0 )
         graphWriter_ = make_shared<gather::GnuPlotGraphWriter>( filename );
   }

   virtual ~PoiseuilleVelocityDataProcessor() {}

   real_t getLastMeanError() const { return lastMeanError_; }
   real_t getLastMaxError() const { return lastMaxError_; }

private:
   real_t acceleration_l_;
   real_t viscosity_l_;
   bool   pipeOrBox_; // true for pipe
   real_t cellsPerRadius_;

   shared_ptr<gather::GnuPlotGraphWriter> graphWriter_;

   real_t lastMeanError_;
   real_t lastMaxError_;

   /// returns relative error per point
   real_t compareToAnalyticalSolution( const std::vector<std::vector<real_t> > & data )
   {
      real_t geometryFactor = pipeOrBox_ ? real_t(4) : real_t(2);

      size_t maxErrPosition = data.size()+5;
      real_t maxErr = -1;
      real_t meanErr = 0;

      for( size_t i = 0; i < data.size(); ++i )
      {
         real_t x = real_c( i ) - cellsPerRadius_ + real_c(0.5);
         real_t analytical = acceleration_l_ / ( geometryFactor * viscosity_l_ ) * ( cellsPerRadius_ * cellsPerRadius_ - x * x );
         real_t simulated  = data[i][1];
         WALBERLA_LOG_DEVEL(i << "  Simulated " << simulated << " - analytical " << analytical );
         real_t diff =  std::abs ( (analytical - simulated) / analytical );
         if ( diff > maxErr ) {
            maxErr = diff;
            maxErrPosition = i;
         }
         meanErr += diff;
      }
      meanErr /= real_c( data.size() );
      WALBERLA_LOG_RESULT("Difference to analytical solution: "
                           << "( mean, max [pos] ) = ( "
                           << meanErr << " ,  " << maxErr << " [ " << maxErrPosition << " ] )" );

      lastMeanError_ = meanErr;
      lastMaxError_ = maxErr;

      return maxErr;
   }


   inline virtual void process(const std::vector<std::vector<real_t> > & data)
   {
      compareToAnalyticalSolution( data );
      if ( graphWriter_ )
         graphWriter_->process( data );
   }

};
*/


int main( int argc, char** argv )
{
   using walberla::uint_t;

   debug::enterTestMode();

   // ----------------------------- Read Configuration ------------------------------------------------

   walberla::Environment env ( argc, argv );
   shared_ptr<Config> config = env.config();

   if ( ! config ) {
      std::cerr << "Call this executable with a configuration file as argument!" << std::endl;
      return 1;
   }

   Config::BlockHandle appConfig = config->getBlock( "Poiseuille" );
   if ( ! appConfig ) {
      std::cerr << "Poiseuille Config Block is missing!" << std::endl;
      return 1;
   }

   real_t diameter         = appConfig.getParameter<real_t>("diameter");
   uint_t cellsPerDiameter = appConfig.getParameter<uint_t>("cellsPerDiameter");
   real_t dx               = diameter / real_c( cellsPerDiameter );
   real_t aspectRatio      = appConfig.getParameter<real_t>("aspectRatio");
   uint_t nrCellsZ         = uint_c( aspectRatio * real_c(cellsPerDiameter) );
   std::string boundary    = appConfig.getParameter<std::string>( "boundary", "pressureDriven" );
   std::string scenario    = appConfig.getParameter<std::string>("scenario", "pipe");
   uint_t timesteps        = appConfig.getParameter<uint_t>("timesteps", 200 );
   //uint_t writeInterval    = appConfig.getParameter<uint_t>( "writeInterval", 50 );
   real_t omega            = appConfig.getParameter<real_t> ("omega");
   real_t viscosity        = appConfig.getParameter<real_t> ("viscosity");
   uint_t nrOfCells [3]    =  { cellsPerDiameter + 2, cellsPerDiameter + 2, nrCellsZ +2 };
   real_t dt               = ( real_t(2) - omega ) * dx * dx / ( real_t(6) * omega * viscosity );
   //real_t viscosity_l      = viscosity * dt / ( dx * dx );

   //real_t requiredAccuracy = appConfig.getParameter<real_t>("requiredAccuracy", real_t(0) );

   // ----------------------------- Create Block Structure  ---------------------------------------------

   if ( scenario == "rect2D" )
      nrOfCells[1] = 1;

   std::vector< real_t > weighting;
   weighting.push_back( real_c( nrOfCells[0]) );
   weighting.push_back( real_c( nrOfCells[1]) );
   weighting.push_back( real_c( nrOfCells[2]) );
   uint_t nrOfProcesses = uint_c( MPIManager::instance()->numProcesses() );
   std::vector<uint_t> processes = math::getFactors( nrOfProcesses, 3, weighting );

   bool xPeriodic = false;
   bool yPeriodic = ( scenario == "rect2D" );
   bool zPeriodic = ( boundary == "forceDriven" );

   uint_t cellsPerBlock[3];
   for( uint_t i = 0; i < 3; ++i )
   {
      if ( nrOfCells[i] % processes[i] == 0 )
         cellsPerBlock[i] = nrOfCells[i] / processes[i];
      else
         cellsPerBlock[i] = (nrOfCells[i] + processes[i]) / processes[i];
   }


   using blockforest::createUniformBlockGrid;

   shared_ptr<StructuredBlockForest>
   blocks = createUniformBlockGrid( processes[0],      processes[1],     processes[2],
                                    cellsPerBlock[0],  cellsPerBlock[1], cellsPerBlock[2],
                                    dx,
                                    true,                             // one block per process
                                    xPeriodic, yPeriodic, zPeriodic,  // periodicity
                                    false );                          // do NOT keep global information



   //  -----------------------------    Fields  --------------------------------------------------------

   auto latticeModel = LM ( SRT( omega ), GuoConstant( Vector3<real_t>() ) ); // force is set by initializer

   BlockDataID pdfFieldID  = lbm::addPdfFieldToStorage( blocks, "PdfField", latticeModel, Vector3<real_t>(0), real_t(1) );
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FField>( blocks, "Flag Field" );


   //typedef lbm::Adaptor<LM>::Density        DensityAdaptor;
   //typedef lbm::Adaptor<LM>::VelocityVector VelocityAdaptor;

   //BlockDataID densityAdaptorID   = field::addFieldAdaptor<DensityAdaptor>  ( blocks, pdfFieldID, "DensityAdaptor" );
   //BlockDataID velocityAdaptorID  = field::addFieldAdaptor<VelocityAdaptor> ( blocks, pdfFieldID, "VelocityAdaptor" );


   typedef lbm::DefaultBoundaryHandlingFactory< LM, FField > BHFactory;

   const FlagUID fluidFlagUID( "Fluid" );
   BlockDataID boundaryHandlingId = BHFactory::addBoundaryHandlingToStorage( blocks, "boundary handling", flagFieldID, pdfFieldID, fluidFlagUID,
                                                                             Vector3<real_t>(), Vector3<real_t>(), real_t(0), real_t(0) );


   //typedef field::ComponentExtractionAdaptor< VelocityAdaptor, 0, 2 > ZVelExtractor;
   //BlockDataID zVelocityAdaptorID = field::addFieldAdaptor<ZVelExtractor> ( blocks, velocityAdaptorID, "Z Vel" );
   //

   typedef lbm::initializer::Poiseuille< BHFactory::BoundaryHandling, LM > PoiseuilleInitializer;
   PoiseuilleInitializer initializer( *blocks, boundaryHandlingId, pdfFieldID,
                                      BHFactory::getNoSlip(), BHFactory::getVelocity0(), BHFactory::getPressure0(), BHFactory::getPressure1() );

   initializer.init( appConfig );
   geometry::setNonBoundaryCellsToDomain<BHFactory::BoundaryHandling> ( *blocks, boundaryHandlingId );

   WALBERLA_ROOT_SECTION() {
      WALBERLA_LOG_INFO( "Timestep chosen as " << dt << " s   ( Omega = " << omega << " )");
      WALBERLA_LOG_INFO( "Block Distribution (x,y,z) " << processes[0] << ", " << processes[1] << ", " << processes[2] );
      WALBERLA_LOG_INFO( "Cells per Block    (x,y,z) " << cellsPerBlock[0] << ", " << cellsPerBlock[1] << ", " << cellsPerBlock[2] );
      WALBERLA_LOG_INFO( "Geometry: " << scenario );
   }



   //  --------------------------   Communication  ------------------------------------------------------

   blockforest::communication::UniformBufferedScheme<stencil::D3Q19> commScheme( blocks );
   //commScheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< stencil::D3Q19, real_t> >( pdfFieldID ) );
   commScheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField > >( pdfFieldID ) );


   //  ----------------------------- Sweep & Timeloop   -------------------------------------------------

   SweepTimeloop timeloop( blocks, timesteps );
   timeloop.add() << BeforeFunction( commScheme, "Communication")
                  << Sweep( BHFactory::BoundaryHandling::getBlockSweep( boundaryHandlingId ), "boundary handling" );

   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LM, FField >( pdfFieldID, flagFieldID, fluidFlagUID ) ), "LBM" );

   //  ----------------------------- Setup Gather Operation   -------------------------------------------------
   /*
   CellInterval cellDomain = blocks->getDomainCellBB();
   cell_idx_t xMid = (cellDomain.max()[0] - cellDomain.min()[0]) / 2;
   cell_idx_t yMid = (cellDomain.max()[1] - cellDomain.min()[1]) / 2;
   cell_idx_t zMid = (cellDomain.max()[2] - cellDomain.min()[2]) / 2;

   CellInterval crossLine = cellDomain;
   crossLine.min()[1] = yMid;
   crossLine.max()[1] = yMid;
   crossLine.min()[2] = zMid;
   crossLine.max()[2] = zMid;

   // Do not gather the noslip boundary cells
   crossLine.min()[0]++;
   crossLine.max()[0]--;

   CellInterval heightLine = cellDomain;
   heightLine.min()[0] = xMid;
   heightLine.max()[0] = xMid;
   heightLine.min()[1] = yMid;
   heightLine.max()[1] = yMid;

   heightLine.min()[2]++;
   heightLine.max()[2]--;

   using namespace gather;

   auto densityWriter  = make_shared<GnuPlotGraphWriter>( "density" );
   auto shearWriter    = make_shared<GnuPlotGraphWriter>( "shear" );
   auto omegaWriter    = make_shared<GnuPlotGraphWriter>( "omega" );

   auto & firstBlock = *blocks->begin();
   auto acceleration_l = firstBlock.getData<PdfField>( pdfFieldID )->latticeModel().forceModel().force().length();
   auto velocityWriter = make_shared<PoiseuilleVelocityDataProcessor>( acceleration_l, viscosity_l,
                                                                       cell_idx_c(cellsPerDiameter / 2 ),
                                                                       (scenario=="pipe"), "velocity" );


   auto densityPacker = make_shared< CellGatherPackInfo<DensityAdaptor,CellInterval> >(
            blocks, densityAdaptorID, heightLine, densityWriter);

   auto velocityPacker = make_shared< CellGatherPackInfo<ZVelExtractor,CellInterval> >(
            blocks, zVelocityAdaptorID, crossLine, velocityWriter);

   gather::MPIGatherScheme gatherScheme( blocks->getBlockStorage(), 0, writeInterval );

   gatherScheme.addPackInfo( densityPacker  );
   gatherScheme.addPackInfo( velocityPacker );

   timeloop.addFuncAfterTimeStep( gatherScheme, "GatherScheme" );


   //  ----------------------------- Timeloop   -------------------------------------------------
   */

   if ( appConfig.isDefined("useGui") ) {
      GUI gui ( timeloop, blocks, argc, argv );
      lbm::connectToGui<LM>( gui );
      gui.run();
   }
   else
   {
      // Remaining Time Logger
      timeloop.addFuncAfterTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "RemainingTimeLogger" );
      timeloop.run();
   }

   /*
   if ( requiredAccuracy > 0.0 )
      WALBERLA_CHECK_LESS( velocityWriter->getLastMeanError(), requiredAccuracy );
   */
   return 0;
}



