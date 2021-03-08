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
//! \file PerformanceMeter.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "PerformanceMeter.h"

#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/mpi/Reduce.h"

#include <iomanip>
#include <sstream>


namespace walberla {
namespace timeloop {


   //*******************************************************************************************************************
   /*! Creates a PerformanceMeter
    *    To actually measure performance, the class has to be connected to a timeloop.
    *    see getBeforeFunction() and getAfterFunction()
    * \param blockStorage block storage is needed to retrieve the FlagField for cell counting
    *******************************************************************************************************************/
   PerformanceMeter::PerformanceMeter( StructuredBlockStorage & blockStorage)
      : blockStorage_ ( blockStorage ), firstTimingStartStopCall_( true )
   {
   }


   //*******************************************************************************************************************
   /*! Returns a pointer to the function that starts all measurements
    * Use this function to connect a PerformanceMeter to a timeloop
    \code
       timeloop.addFuncBeforeTimeStep( perfMeter.getBeforeFunction() );
    \endcode
    *******************************************************************************************************************/
   std::function<void () > PerformanceMeter::getBeforeFunction()
   {
      return std::bind ( &PerformanceMeter::timingStart, this );
   }


   //*******************************************************************************************************************
   /*! Returns a pointer to the function that ends all measurements
    * Use this function to connect a PerformanceMeter to a timeloop
    \code
       timeloop.addFuncAfterTimeStep ( perfMeter.getAfterFunction() );
    \endcode
    *******************************************************************************************************************/
   std::function<void () > PerformanceMeter::getAfterFunction()
   {
      return std::bind ( &PerformanceMeter::timingEnd, this );
   }


   //*******************************************************************************************************************
   /*! Adds a performance measurement, associated with all cells (usually MLUPS)
    * \param name          name of the performance number (f.e. MLUPS)
    * \param countFunction this function is called very countFreq'th  call and has to return the
    *                      number of 'active' cells on the given block
    * \param countFreq     frequency of counting, if zero the counting happens once in the beginning
    * \param scaling       The performance number is multiplied with this scaling factor before printing.
    *                      Example MLUPS , for the "M" a factor of 1e-6 is needed
    *
    * Behaves like addMeasurement() function above, with n= total number of cells
    *
    *******************************************************************************************************************/
   void PerformanceMeter::addMeasurement ( const std::string & name, const CountFunction& countFunction,
                                           uint_t countFreq, real_t scaling  )
   {
      measurements_.emplace_back( countFunction, name, scaling, countFreq );
      uint_t cells = 0;
      for( auto block = blockStorage_.begin(); block != blockStorage_.end(); ++block )
         cells += countFunction( *block );

      measurements_.back().avgCellsPerTimeStep = real_c ( cells );
   }


   //*******************************************************************************************************************
   /*! Adds a performance measurement, associated with all cells (usually MLUPS)
    * \param name      name of the performance number (f.e. MLUPS)
    * \param scaling   The performance number is multiplied with this scaling factor before printing.
    *                  Example MLUPS , for the "M" a factor of 1e-6 is needed
    *
    * Behaves like addMeasurement() function above, with n= total number of cells
    *
    *******************************************************************************************************************/
   void PerformanceMeter::addMeasurement ( const std::string & name, real_t scaling )
   {
      measurements_.emplace_back( CountFunction(), name, scaling, uint_t(0) );

      uint_t cellsOnProcess = 0;
      for( auto block = blockStorage_.begin(); block != blockStorage_.end(); ++block )
      {
        cellsOnProcess += blockStorage_.getNumberOfXCells( *block ) *
                          blockStorage_.getNumberOfYCells( *block ) *
                          blockStorage_.getNumberOfZCells( *block );
      }

      measurements_.back().avgCellsPerTimeStep = real_c ( cellsOnProcess );
      measurements_.back().counts = 1;
   }



   //*******************************************************************************************************************
   /*! Logs all added performance measurements on the root process
    *******************************************************************************************************************/
   void PerformanceMeter::logResultOnRoot( )
   {
      if( measurements_.empty())
         return;

      std::stringstream ss;
      print( ss, 0);
      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_LOG_RESULT( ss.str() );
      }
   }




   //*******************************************************************************************************************
   /*! Prints all added measurements to the given output stream on specified process
    * \param os            Output stream where results are printed
    * \param targetRank    the MPI world rank of the process that should print.
    *                      or a negative value, if all process should print.
    *******************************************************************************************************************/
   void PerformanceMeter::print( std::ostream & os, int targetRank )
   {
      if( measurements_.empty())
         return;

      std::vector<real_t> totalNrCells;
      reduce( totalNrCells, targetRank );
      if ( MPIManager::instance()->worldRank() == targetRank || targetRank < 0 )
      {
         os << "Performance Results: " << std::endl;
         for( uint_t i = 0; i < measurements_.size(); ++i ) {
            real_t res = totalNrCells[i] / real_c( timer_.average() );
            res *= measurements_[i].scaling;
            os << " - ";
            os << std::setw (16) << std::left  << measurements_[i].name ;
            os << std::setw (8)  << std::right << res << std::endl;
         }
      }
   }


   //*******************************************************************************************************************
   /*! Collects the result of the PerformanceMeter from all MPI_Processes
    *
    * \param targetRank  MPI world rank of the process where the PerformanceMeter is reduced to
    *                    on all other processes a null pointer is returned. If targetRank < 0 all processes
    *                    have a valid result
    *
    * \return  a map of measurement-name to measurement value if (worldRank == targetRank ) || targetRank < 0
    *          and a null pointer otherwise
    *******************************************************************************************************************/
   shared_ptr< std::map<std::string, real_t> > PerformanceMeter::getReduced ( int targetRank )
   {
      using ResultMap = std::map<std::string, real_t>;
      shared_ptr < ResultMap > res;

      std::vector<real_t> totalNrCells;
      reduce( totalNrCells, targetRank );

      if ( MPIManager::instance()->worldRank() == targetRank || targetRank < 0 )
      {
         res = make_shared< ResultMap > ();

         for( uint_t i = 0; i < measurements_.size(); ++i )
         {
            real_t value = totalNrCells[i] / real_c( timer_.average() );
            value *= measurements_[i].scaling;
            (*res) [ measurements_[i].name ] = value;
         }
      }

      return res;
   }


   //*******************************************************************************************************************
   /*! Call operator which starts/stops the measurement
    *
    * When the PerformanceMeter is added to a timeloop this operator is called.
    * One measurement period, goes from one operator() call to the next call of operator()
    * -> Add this before all sweeps in the timeloop.
    *
    * If more fine grained control is needed use getBeforeFunction() and getAfterFunction()
    * then only the sweeps between these functions are measured
    *******************************************************************************************************************/
   void PerformanceMeter::operator()()
   {
      timer_.end();

      // When created, the timer starts, but we have one invalid measurement
      // which is discarded
      if( firstTimingStartStopCall_ ) {
         timer_.reset();
         firstTimingStartStopCall_ = false;
      }
      else
         updateCellCounts();

      timer_.start();
   }

   //*******************************************************************************************************************
   /*! Resets Timings. Added measurements are not cleared.
    *******************************************************************************************************************/
   void PerformanceMeter::clear()
   {
      timer_.reset();
      firstTimingStartStopCall_ = true;
   }


   //*******************************************************************************************************************
   /*! This function is wrapped and returned in getBeforeFunction()
    *******************************************************************************************************************/
   void PerformanceMeter::timingStart()
   {
      timer_.start();
   }

   //*******************************************************************************************************************
   /*! This function is wrapped and returned in getAfterFunction()
    *
    * Stops the timing, if dynamic flags are set, the cells are counted and a running average of the
    * cell count is computed.
    *******************************************************************************************************************/
   void PerformanceMeter::timingEnd()
   {
      timer_.end();
      updateCellCounts();
   }

   //*******************************************************************************************************************
   /*! Loops over measurements and updates them if necessary
    *******************************************************************************************************************/
   void PerformanceMeter::updateCellCounts()
   {
      WALBERLA_ASSERT_GREATER( timer_.getCounter(), 0 );

      uint_t ts = timer_.getCounter();

      for( auto measureIt = measurements_.begin(); measureIt != measurements_.end(); ++measureIt )
      {
         if( measureIt->countingFreq > 0 &&  ts % measureIt->countingFreq == 0 )
         {
            uint_t cells = 0;
            for( auto block = blockStorage_.begin(); block != blockStorage_.end(); ++block )
               cells += measureIt->countFunction( *block );

            real_t cellsLastTimeStep = real_c ( cells );

            measureIt->counts += 1;
            real_t rCounts = real_c( measureIt->counts );
            measureIt->avgCellsPerTimeStep = (rCounts - real_t(1.0) ) / rCounts * measureIt->avgCellsPerTimeStep +
                                                        real_t(1.0)   / rCounts * cellsLastTimeStep;
         }
      }
   }



   //*******************************************************************************************************************
   /*! Communicates the total number of cells among processes
    * \param[out] reduced    vector that after returning holds the total number of cells for each Measurement
    *                        on the process with has world rank = targetRank, or on all processes if targetRank < 0
    * \param[in]  targetRank see above
    *******************************************************************************************************************/
   void PerformanceMeter::reduce ( std::vector<real_t> & reduced, int targetRank )
   {
      if ( measurements_.empty() )
         return;

      reduced.clear();
      for( auto cellType = measurements_.begin(); cellType != measurements_.end(); ++cellType )
         reduced.push_back ( cellType->avgCellsPerTimeStep );

      WALBERLA_NON_MPI_SECTION() {
         return;
      }

      if( targetRank >= 0 )
         mpi::reduceInplace( reduced, mpi::SUM, targetRank, MPI_COMM_WORLD );
      else
         mpi::allReduceInplace( reduced, mpi::SUM, MPI_COMM_WORLD );

   }




} // namespace timeloop
} // namespace walberla



