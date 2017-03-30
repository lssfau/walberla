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
//! \file PerformanceEvaluation.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/Hostname.h"
#include "core/Set.h"
#include "core/waLBerlaBuildInfo.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"
#include "core/uid/SUID.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/CellCounter.h"
#include "field/FlagUID.h"

#include <cstdlib>
#include <map>
#include <string>
#include <sstream>



namespace walberla {
namespace lbm {


//**********************************************************************************************************************
/*!
*   \brief Class for evaluating the performance of LBM simulations
*/
//**********************************************************************************************************************
template< typename CellCounter_T, typename FluidCellCounter_T >
class PerformanceEvaluationBase
{
public:

   PerformanceEvaluationBase( const weak_ptr< StructuredBlockStorage > & blocks,
                              const CellCounter_T & cellCounter, const FluidCellCounter_T & fluidCellCounter,
                              const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                              const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );
   
   void refresh();

   void logResultOnRoot( const uint_t timeSteps, const double time ) const
   {
      WALBERLA_LOG_RESULT_ON_ROOT( "Simulation performance:\n" << loggingString( timeSteps, time ) );
   }

   void logInfoOnRoot( const uint_t timeSteps, const double time ) const
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Simulation performance:\n" << loggingString( timeSteps, time ) );
   }

   std::string loggingString( const uint_t timeSteps, const double time ) const;
   
   void getResultsForSQLOnRoot( std::map< std::string, int > &         integerProperties,
                                std::map< std::string, double > &      realProperties,
                                std::map< std::string, std::string > & stringProperties,
                                const uint_t timeSteps, const double time );
   
   static int processes() { return mpi::MPIManager::instance()->numProcesses(); }

   int threads() const { return processes() * threadsPerProcess_; }
   int cores()   const { return ( threadsPerCore_ == 0 ) ? 0 : ( threads() / threadsPerCore_ ); }

   uint64_t allFineCells() const
   {
      uint64_t c( uint64_t(0) );
      for( uint_t i = uint_t(0); i < levels_; ++i )
         c += cells_.numberOfCells(i) * uint64_c( math::uintPow8( levels_ - uint_t(1) - i ) );
      return c;
   }

   double mlups( const uint_t timeSteps, const double time ) const
   {
      double m( 0.0 );
      for( uint_t i = uint_t(0); i < levels_; ++i )
         m += double_c( timeSteps * math::uintPow2(i) ) * double_c( cells_.numberOfCells(i) );
      return m / ( time * 1000000.0 );
   }

   double mlupsPerProcess( const uint_t timeSteps, const double time ) const
   {
      return mlups( timeSteps, time ) / processes();
   }

   double mlupsPerCore( const uint_t timeSteps, const double time ) const
   {
      return ( cores() == 0 ) ? 0.0 : ( mlups( timeSteps, time ) / cores() );
   }

   double vMlups( const uint_t timeSteps, const double time ) const
   {
      double m( 0.0 );
      for( uint_t i = uint_t(0); i < levels_; ++i )
         m += double_c( timeSteps * math::uintPow2( levels_ - uint_t(1) ) ) *
              double_c( uint64_c( math::uintPow8( levels_ - uint_t(1) - i ) ) * cells_.numberOfCells(i) );
      return m / ( time * 1000000.0 );
   }

   double vMlupsPerProcess( const uint_t timeSteps, const double time ) const
   {
      return vMlups( timeSteps, time ) / processes();
   }

   double vMlupsPerCore( const uint_t timeSteps, const double time ) const
   {
      return ( cores() == 0 ) ? 0.0 : ( vMlups( timeSteps, time ) / cores() );
   }

   double mflups( const uint_t timeSteps, const double time ) const
   {
      double m( 0.0 );
      for( uint_t i = uint_t(0); i < levels_; ++i )
         m += double_c( timeSteps * math::uintPow2(i) ) * double_c( fluidCells_.numberOfCells(i) );
      return m / ( time * 1000000.0 );
   }

   double mflupsPerProcess( const uint_t timeSteps, const double time ) const
   {
      return mflups( timeSteps, time ) / processes();
   }

   double mflupsPerCore( const uint_t timeSteps, const double time ) const
   {
      return ( cores() == 0 ) ? 0.0 : ( mflups( timeSteps, time ) / cores() );
   }

   double vMflups( const uint_t timeSteps, const double time ) const
   {
      double m( 0.0 );
      for( uint_t i = uint_t(0); i < levels_; ++i )
         m += double_c( timeSteps * math::uintPow2( levels_ - uint_t(1) ) ) *
              double_c( uint64_c( math::uintPow8( levels_ - uint_t(1) - i ) ) * fluidCells_.numberOfCells(i) );
      return m / ( time * 1000000.0 );
   }

   double vMflupsPerProcess( const uint_t timeSteps, const double time ) const
   {
      return vMflups( timeSteps, time ) / processes();
   }

   double vMflupsPerCore( const uint_t timeSteps, const double time ) const
   {
      return ( cores() == 0 ) ? 0.0 : ( vMflups( timeSteps, time ) / cores() );
   }

   static double timeStepsPerSecond( const uint_t timeSteps, const double time ) { return double_c( timeSteps ) / time; }

   double fineTimeStepsPerSecond( const uint_t timeSteps, const double time ) const
   {
      return double_c( timeSteps * math::uintPow2( levels_ - uint_t(1) ) ) / time;
   }

private:

   int threadsPerProcess_;
   int threadsPerCore_;

   weak_ptr< StructuredBlockStorage > blocks_;
   uint_t levels_;

   CellCounter_T cells_;
   FluidCellCounter_T fluidCells_;

}; // class PerformanceEvaluationBase



//**********************************************************************************************************************
/*!
*   \brief Class for evaluating the performance of LBM simulations using fields
*
*   Assumes that in-between creating an object of this class and calling any of the member functions the number of cells
*   and the number of fluid cells do not change! For simulations with static geometry, this is always the case.
*/
//**********************************************************************************************************************
template< typename FlagField_T >
class PerformanceEvaluation : public PerformanceEvaluationBase< field::CellCounter< FlagField_T >, field::CellCounter< FlagField_T > >
{
public:
   PerformanceEvaluation( const weak_ptr< StructuredBlockStorage > & blocks,
                          const ConstBlockDataID & flagFieldId, const Set< FlagUID > & fluid,
                          const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                          const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
                          : PerformanceEvaluationBase< field::CellCounter< FlagField_T >, field::CellCounter< FlagField_T > >(
                              blocks,
                              field::CellCounter< FlagField_T >( blocks, flagFieldId, Set< FlagUID >::emptySet(), requiredSelectors, incompatibleSelectors ),
                              field::CellCounter< FlagField_T >( blocks, flagFieldId, fluid, requiredSelectors, incompatibleSelectors ),
                              requiredSelectors, incompatibleSelectors )
   {
   }
};


template< typename CellCounter_T, typename FluidCellCounter_T >
PerformanceEvaluationBase< CellCounter_T, FluidCellCounter_T >::PerformanceEvaluationBase(
                                                                   const weak_ptr< StructuredBlockStorage > & blocks,
                                                                   const CellCounter_T & cellCounter, const FluidCellCounter_T & fluidCellCounter,
                                                                   const Set<SUID> & /*requiredSelectors*/, const Set<SUID> & /*incompatibleSelectors*/ )
   : threadsPerProcess_( 1 ), threadsPerCore_( 0 ),
     blocks_( blocks ),
     cells_( cellCounter ),
     fluidCells_( fluidCellCounter )
{
#ifdef _OPENMP
   if( std::getenv( "OMP_NUM_THREADS" ) == NULL )
      WALBERLA_ABORT( "If you are using a version of the program that was compiled with OpenMP you have to "
                      "specify the environment variable \'OMP_NUM_THREADS\' accordingly!" );
   threadsPerProcess_ = std::atoi( std::getenv( "OMP_NUM_THREADS" ) );
#endif

   if( std::getenv( "THREADS_PER_CORE" ) )
      threadsPerCore_ = std::atoi( std::getenv( "THREADS_PER_CORE" ) );

   refresh();
}



template< typename CellCounter_T, typename FluidCellCounter_T >
void PerformanceEvaluationBase< CellCounter_T, FluidCellCounter_T >::refresh()
{
   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'PerformanceEvaluation' for a block storage object that doesn't exist anymore" );
   
   levels_ = blocks->getNumberOfLevels();
   
   cells_();
   fluidCells_();
}



template< typename CellCounter_T, typename FluidCellCounter_T >
std::string PerformanceEvaluationBase< CellCounter_T, FluidCellCounter_T >::loggingString( const uint_t timeSteps, const double time ) const
{
   std::ostringstream oss;

   std::string na( "n/a *)" );

   std::ostringstream threadsPerCoreString;
   threadsPerCoreString << threadsPerCore_;

   std::ostringstream coresString;
   coresString << cores();

   oss <<   "- processes:   " << processes()
      << "\n- threads:     " << threads() << " (threads per process = " << threadsPerProcess_
      << ", threads per core = " << ( ( threadsPerCore_ == 0 ) ? na : threadsPerCoreString.str() ) << ")"
      << "\n- cores:       " << ( ( threadsPerCore_ == 0 ) ? na : coresString.str() )
      << "\n- time steps:  " << timeSteps;

   if( levels_ > uint_t(1) )
   {
      oss << " (on the coarsest grid, " << ( timeSteps * math::uintPow2( levels_ - uint_t(1) ) ) << " on the finest grid)";
   }

   oss << "\n- time:        " << time << " sec"
      << "\n- cells:       " << cells_.numberOfCells();

   if( levels_ > uint_t(1) )
   {
      oss << " (" << allFineCells() << " if everything were fine -> data reduction by factor of "
         << ( real_c( allFineCells() ) / real_c( cells_.numberOfCells() ) ) << ")";
   }

   oss << "\n- fluid cells: " << fluidCells_.numberOfCells() << " ("
      << ( real_c(100) * real_c( fluidCells_.numberOfCells() ) / real_c( cells_.numberOfCells() ) ) << " % of all cells)";

   if( levels_ > uint_t(1) )
   {
      oss << "\n- distribution of cells to different grid levels:";
      for( uint_t i = uint_t(0); i < levels_; ++i )
         oss << "\n   + level " << i <<": " << cells_.numberOfCells(i) << " cells (" << fluidCells_.numberOfCells(i) << " fluid cells = "
         << ( real_c(100) * real_c( fluidCells_.numberOfCells(i) ) / real_c( cells_.numberOfCells(i) ) )
         << " % of all cells on this level)";
   }

   std::ostringstream mlupsPerCoreString;
   mlupsPerCoreString << mlupsPerCore( timeSteps, time );

   std::ostringstream mflupsPerCoreString;
   mflupsPerCoreString << mflupsPerCore( timeSteps, time );

   oss << "\n- performance: " << mlups( timeSteps, time ) << " MLUPS (million lattice cell updates per second)"
      << "\n               " << mlupsPerProcess( timeSteps, time ) << " MLUPS / process"
      << "\n               " << ( ( threadsPerCore_ == 0 ) ? na : mlupsPerCoreString.str() ) << " MLUPS / core"
      << "\n               " << mflups( timeSteps, time ) << " MFLUPS (million fluid lattice cell updates per second)"
      << "\n               " << mflupsPerProcess( timeSteps, time ) << " MFLUPS / process"
      << "\n               " << ( ( threadsPerCore_ == 0 ) ? na : mflupsPerCoreString.str() ) << " MFLUPS / core"
      << "\n               " << timeStepsPerSecond( timeSteps, time ) << " time steps / second";

   if( levels_ > uint_t(1) )
   {
      std::ostringstream vMlupsPerCoreString;
      vMlupsPerCoreString << vMlupsPerCore( timeSteps, time );

      std::ostringstream vMflupsPerCoreString;
      vMflupsPerCoreString << vMflupsPerCore( timeSteps, time );

      oss << "\n- 'virtual' performance (if everything were fine): " << vMlups( timeSteps, time ) << " MLUPS (million lattice cell updates per second)"
         << "\n                                                   " << vMlupsPerProcess( timeSteps, time ) << " MLUPS / process"
         << "\n                                                   " << ( ( threadsPerCore_ == 0 ) ? na : vMlupsPerCoreString.str() ) << " MLUPS / core"
         << "\n                                                   " << vMflups( timeSteps, time ) << " MFLUPS (million fluid lattice cell updates per second)"
         << "\n                                                   " << vMflupsPerProcess( timeSteps, time ) << " MFLUPS / process"
         << "\n                                                   " << ( ( threadsPerCore_ == 0 ) ? na : vMflupsPerCoreString.str() ) << " MFLUPS / core"
         << "\n                                                   " << fineTimeStepsPerSecond( timeSteps, time ) << " fine time steps / second";
   }

   oss << "\n- build / run information:"
      << "\n   + host machine:   " << getHostName()
      << "\n   + build machine:  " << WALBERLA_BUILD_MACHINE
      << "\n   + git SHA1:       " << WALBERLA_GIT_SHA1
      << "\n   + build type:     " << WALBERLA_BUILD_TYPE
      << "\n   + compiler flags: " << WALBERLA_COMPILER_FLAGS;

   if( threadsPerCore_ == 0 )
      oss << "\n\n  *) only available if environment variable 'THREADS_PER_CORE' is set";

   return oss.str();
}



template< typename CellCounter_T, typename FluidCellCounter_T >
void PerformanceEvaluationBase< CellCounter_T, FluidCellCounter_T >::getResultsForSQLOnRoot( std::map< std::string, int > &         integerProperties,
                                                                                             std::map< std::string, double > &      realProperties,
                                                                                             std::map< std::string, std::string > & stringProperties,
                                                                                             const uint_t timeSteps, const double time )
{
   WALBERLA_NON_ROOT_SECTION()
   {
      return;
   }

   integerProperties[ "levels" ]            = int_c( levels_ );
   integerProperties[ "processes" ]         = processes();
   integerProperties[ "threads" ]           = threads();
   integerProperties[ "cores" ]             = cores();
   integerProperties[ "threadsPerProcess" ] = threadsPerProcess_;
   integerProperties[ "threadsPerCore" ]    = threadsPerCore_;

   integerProperties[ "timeSteps" ] = int_c( timeSteps );
   if( levels_ > uint_t(1) )
      integerProperties[ "fineTimeSteps" ] = int_c( timeSteps * math::uintPow2( levels_ - uint_t(1) ) );

   realProperties[ "time" ] = real_c( time );

   realProperties[ "cells" ] = real_c( cells_.numberOfCells() );
   if( levels_ > uint_t(1) )
      realProperties[ "refinementCellsReduction" ] = real_c( allFineCells() ) / real_c( cells_.numberOfCells() );
   realProperties[ "fluidCells" ] = real_c( fluidCells_.numberOfCells() );

   if( levels_ > uint_t(1) )
   {
      for( uint_t i = uint_t(0); i < levels_; ++i )
      {
         std::ostringstream cells_i;
         std::ostringstream fluidCells_i;

         cells_i << "cells_" << i;
         fluidCells_i << "fluidCells_" << i;

         realProperties[ cells_i.str() ] = real_c( cells_.numberOfCells(i) );
         realProperties[ fluidCells_i.str() ] = real_c( fluidCells_.numberOfCells(i) );
      }
   }

   realProperties[ "MLUPS" ]              = double_c( mlups( timeSteps, time ) );
   realProperties[ "MLUPS_process" ]      = double_c( mlupsPerProcess( timeSteps, time ) );
   realProperties[ "MLUPS_core" ]         = double_c( mlupsPerCore( timeSteps, time ) );
   realProperties[ "MFLUPS" ]             = double_c( mflups( timeSteps, time ) );
   realProperties[ "MFLUPS_process" ]     = double_c( mflupsPerProcess( timeSteps, time ) );
   realProperties[ "MFLUPS_core" ]        = double_c( mflupsPerCore( timeSteps, time ) );
   realProperties[ "timeStepsPerSecond" ] = double_c( timeStepsPerSecond( timeSteps, time ) );

   if( levels_ > uint_t(1) )
   {
      realProperties[ "vMLUPS" ]                 = double_c( vMlups( timeSteps, time ) );
      realProperties[ "vMLUPS_process" ]         = double_c( vMlupsPerProcess( timeSteps, time ) );
      realProperties[ "vMLUPS_core" ]            = double_c( vMlupsPerCore( timeSteps, time ) );
      realProperties[ "vMFLUPS" ]                = double_c( vMflups( timeSteps, time ) );
      realProperties[ "vMFLUPS_process" ]        = double_c( vMflupsPerProcess( timeSteps, time ) );
      realProperties[ "vMFLUPS_core" ]           = double_c( vMflupsPerCore( timeSteps, time ) );
      realProperties[ "fineTimeStepsPerSecond" ] = double_c( fineTimeStepsPerSecond( timeSteps, time ) );
   }

   stringProperties[ "hostMachine" ]   = std::string( getHostName() );
   stringProperties[ "buildMachine" ]  = std::string( WALBERLA_BUILD_MACHINE );
   stringProperties[ "gitVersion" ]    = std::string( WALBERLA_GIT_SHA1 );
   stringProperties[ "buildType" ]     = std::string( WALBERLA_BUILD_TYPE );
   stringProperties[ "compilerFlags" ] = std::string( WALBERLA_COMPILER_FLAGS );
}

} // namespace lbm
} // namespace walberla
