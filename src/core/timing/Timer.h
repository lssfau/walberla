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
//! \file Timer.h
//! \ingroup core
//! \author Klaus Iglberger
//!
//! Copyright (C) 2009 Klaus Iglberger
//! Taken from "pe Physics Engine" with small changes
//
//======================================================================================================================

#pragma once

#include "CpuPolicy.h"
#include "ReduceType.h"
#include "WcPolicy.h"
#include "core/DataTypes.h"

#include "core/mpi/RecvBuffer.h"
#include "core/mpi/Reduce.h"
#include "core/mpi/SendBuffer.h"

#include <iomanip>
#include <iostream>
#include <limits>


namespace walberla {
namespace timing {


//======================================================================================================================
//
//  CLASS DEFINITION
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Progress timer for time and performance measurements.
// \ingroup timing
//
// The Timer class offers timing & benchmarking functionality for all kinds of applications.
// The following example code demonstrates the use of the WcTimer class, which combines the
// Timer class template with the WcPolicy for wall clock time measurements, for a single time
// measurement:

   \code
   // Creating a new wall clock timer immediately starts a new time measurement
   WcTimer timer;

   ...  // program or code fragment to be measured

   // Stopping the time measurement
   timer.end();

   // Evaluation of the measured time
   double time = timer.last();
   \endcode

// The timer class additionally offers the functionality to start several time measurements in
// order to evaluate minimal, maximal or average times. The next example demonstrates a possible
// setup for such a series of time measurements:

   \code
   // Creating a new wall clock timer
   WcTimer timer;

   ...  // Additional setup code

   // Starting 10 wall clock time measurements
   for( unsigned int i=0; i<10; ++i ) {
      timer.start();
      ...  // program or code fragment to be measured
      timer.end();
   }

   // After the measurements, the desired timing results can be calculated, as for instance the
   // average wall clock time
   double average = timer.average();
   \endcode
*/
template< typename TP >  // Timing policy
class Timer
{
   template< typename T,     // Element type of SendBuffer
             typename G,     // Growth policy of SendBuffer
             typename TP2 >  // Element type of vector
   friend mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const Timer<TP2> & t );

   template< typename T,     // Element type  of RecvBuffer
             typename TP2 >  // Element type of vector
   friend mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, Timer<TP2> & t );
public:
   //**Type definitions*************************************************************************************************
   typedef TP  TimingPolicy;  //!< Timing policy of the Timer.
   //*******************************************************************************************************************

   //**Constructor******************************************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline Timer();
   explicit inline Timer( uint_t counter, double min, double max, double total, double sumOfSquares );

   //@}
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   // No explicitly declared destructor.
   //*******************************************************************************************************************

   //**Timing functions*************************************************************************************************
   /*!\name Timing functions */
   //@{
   inline void start();
   inline void end  ();
   inline void reset();
   //@}
   //*******************************************************************************************************************

   //**Get functions****************************************************************************************************
   /*!\name Get functions */
   //@{
   inline uint_t getCounter() const;
   //@}
   //*******************************************************************************************************************

   //**Time evaluation functions****************************************************************************************
   /*!\name Time evaluation functions */
   //@{
   inline double total()        const;
   inline double sumOfSquares() const;
   inline double average()      const;
   inline double variance()     const;
   inline double min()          const;
   inline double max()          const;
   inline double last()         const;
   //@}
   //*******************************************************************************************************************

   //**Utility functions************************************************************************************************
   /*!\name Utility functions for Timers */
   //@{
   inline void merge( const Timer & other );
   //@}
   //*******************************************************************************************************************

private:

   uint_t counter_;      //!< Number of performed time measurements.
   double start_;        //!< Start of the current time measurement.
   double end_;          //!< End of the current time measurement.
   double time_;         //!< The total elapsed time of all measurements.
   double sumOfSquares_; //!< Sum of each (time measurement)^2
   double min_;          //!< The minimal time of all measurements.
   double max_;          //!< The maximal time of all measurements.
   double last_;         //!< The last measured time.
};
//**********************************************************************************************************************

//======================================================================================================================
//
//  CONSTRUCTOR
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Constructor of the Timer class.
//
// The creation of a new timer immediately starts a new time measurement. It is possible to
// either restart the time measurement at a specific point of time or to continue the time
// measurement and to end it via the end() function.
*/
template< typename TP >  // Timing policy
inline Timer<TP>::Timer()
   : counter_     ( 0   )
   , start_       ( 0.0 )
   , end_         ( 0.0 )
   , time_        ( 0.0 )
   , sumOfSquares_( 0.0 )
   , min_         ( std::numeric_limits<double>::max() )
   , max_         ( 0.0 )
   , last_        ( 0.0 )
{
   // Starting the time measurement
   start();
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Constructor of the Timer class.
// Initializes the timer with
// \param counter      number of timing measurements already done
// \param min          the minimum time of the measurements
// \param max          the maximum time of the measurements
// \param total        the total time of all measurements
// \param sumOfSquares each measurement time squared, then summed
*/
template< typename TP>
inline Timer<TP>::Timer( uint_t _counter, double _min, double _max,
                         double _total, double _sumOfSquares )
   : counter_     ( _counter      )
   , start_       ( 0.0           )
   , end_         ( 0.0           )
   , time_        ( _total        )
   , sumOfSquares_( _sumOfSquares )
   , min_         ( _min          )
   , max_         ( _max          )
   , last_        (  0.0          )
{
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  TIMING FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Starting a single time measurement.
//
// \return void
//
// This function starts a single time measurement.
*/
template< typename TP >  // Timing policy
inline void Timer<TP>::start()
{
   // Starting the time measurement and calculating a time stamp
   start_ = TimingPolicy::getTimestamp();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Ending a single time measurement.
//
// \return void
//
// This function ends the currently running time measurement and performs the necessary
// statistical calculations.
*/
template< typename TP >  // Timing policy
inline void Timer<TP>::end()
{
   // Stopping the time measurement and calculating a time stamp
   end_ = TimingPolicy::getTimestamp();

   // Increasing the counter
   ++counter_;

   // Calculating the wall clock and CPU time
   const double diff( end_ - start_ );

   // Average time measurement
   time_ += diff;

   // Update sumOfSquares for variance calculation
   sumOfSquares_ += (diff * diff);

   // Minimum time measurement
   if( diff < min_ ) min_ = diff;

   // Maximum time measurement
   if( diff  > max_ ) max_ = diff;

   // Last time measurement
   last_ = diff;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Resetting the timer.
//
// \return void
//
// This function completely resets the timer and all information on the performed time
// measurements. In order to start a new time measurement, the start() function has to
// be used.
*/
template< typename TP >  // Timing policy
inline void Timer<TP>::reset()
{
   counter_      = 0;
   start_        = 0.0;
   end_          = 0.0;
   time_         = 0.0;
   sumOfSquares_ = 0.0;
   min_          = std::numeric_limits<double>::max();
   max_          = 0.0;
   last_         = 0.0;
   start();
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  GET FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Returns the total number of time measurements performed by this timer.
//
// \return The number of performed time measurements.
*/
template< typename TP >  // Timing policy
inline uint_t Timer<TP>::getCounter() const
{
   return counter_;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  TIME EVALUATION FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Returns the total elapsed time of all performed time measurements.
//
// \return The total elapsed time of all time measurements.
*/
template< typename TP >  // Timing policy
inline double Timer<TP>::total() const
{
   return time_;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the sum of each (time measurement)^2
//
// Can be used to calculate the variance using the formula
// Var(x) = 1/N *  ( E(x^2) - E(x)^2)
//
// \return sum of each (time measurement)^2
*/
template< typename TP >  // Timing policy
inline double Timer<TP>::sumOfSquares() const
{
   return sumOfSquares_;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the average time of all performed time measurements.
//
// \return The average time.
*/
template< typename TP >  // Timing policy
inline double Timer<TP>::average() const
{
   return time_ / real_c (counter_ );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the variance of all performed time measurements.
//
// \return The variance 1/N * ( E(x^2) - E(x)^2 ) where E is the expectation or average value
*/
template< typename TP >  // Timing policy
inline double Timer<TP>::variance() const
{
   double dblCounter           = real_c (counter_ );
   double expectedSumOfSquares = sumOfSquares_ / dblCounter;
   double expectedSumSquared   = (time_ / dblCounter) * (time_/dblCounter);

   return expectedSumOfSquares - expectedSumSquared;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the minimal time of all performed time measurements.
//
// \return The minimal time.
*/
template< typename TP >  // Timing policy
inline double Timer<TP>::min() const
{
   return min_;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the maximal time of all performed time measurements.
//
// \return The maximal time.
*/
template< typename TP >  // Timing policy
inline double Timer<TP>::max() const
{
   return max_;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the last measured time.
//
// \return The last measured time.
*/
template< typename TP >  // Timing policy
inline double Timer<TP>::last() const
{
   return last_;
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Adds results of an other Timer to own results
//
*/
template< typename TP >  // Timing policy
inline void Timer<TP>::merge( const Timer<TP> & other )
{
   counter_      += other.counter_;
   time_         += other.time_;
   sumOfSquares_ += other.sumOfSquares_;
   min_           = std::min( min_, other.min_ );
   max_           = std::max( max_, other.max_ );
}
//**********************************************************************************************************************


//======================================================================================================================
//
//  REDUCTION
//
//======================================================================================================================

//**********************************************************************************************************************
/*! Returns a reduced Timer, holding information from all processes
 *
 * \param timer      Timer which should be reduced
 * \param rt         Specified the method how the reduction is done. See documentation for ReduceType
 * \param targetRank the world rank of the target process. Or negative value for an all-reduction
 *                   operation
 *
 * \return  a nonzero shared pointer to a Timer on the process with rank "targetRank"
 *          and a zero pointer otherwise. For the all-reduction a valid timer is returned on all processes.
 **********************************************************************************************************************/
template< typename TP >
shared_ptr<Timer<TP> > getReduced( Timer<TP>& timer, ReduceType rt, int targetRank )
{
   WALBERLA_NON_MPI_SECTION() {
      return make_shared<Timer<TP> >( timer );
   }

   double val; //value to be reduced
   switch (rt)
   {
   case REDUCE_MIN  :
      val = timer.min();
      break;
   case REDUCE_AVG  :
      val = timer.average();
      break;

   case REDUCE_MAX  :
      val = timer.max();
      break;

   case REDUCE_TOTAL:
      val = timer.total();
      break;

   default:
      WALBERLA_ABORT( "Unknown reduce type" );
      break;
   }

   double min;
   double max;
   double total;
   double sumOfSquares;

   if (targetRank >= 0)
   {
      min          = mpi::reduce(val, mpi::MIN, targetRank);
      max          = mpi::reduce(val, mpi::MAX, targetRank);
      total        = mpi::reduce(val, mpi::SUM, targetRank);
      sumOfSquares = mpi::reduce(val*val, mpi::SUM, targetRank);
   } else
   {
      min          = mpi::allReduce(val, mpi::MIN);
      max          = mpi::allReduce(val, mpi::MAX);
      total        = mpi::allReduce(val, mpi::SUM);
      sumOfSquares = mpi::allReduce(val*val, mpi::SUM);
   }

   //uint_t counter, double min, double max, double total, double sumOfSquares
   if ( targetRank < 0 || targetRank == MPIManager::instance()->worldRank() )
      return make_shared<Timer<TP> >( mpi::MPIManager::instance()->numProcesses(), min, max, total, sumOfSquares  );

   return nullptr;
}


//======================================================================================================================
//
//  OSTREAM OVERLOAD
//
//======================================================================================================================

template< typename TP >  // Timing policy
std::ostream & operator<< ( std::ostream & os, const Timer<TP> & timer )
{
   os << std::fixed << std::setprecision(3) <<
         "average: " << timer.average() <<
         " | min: " << timer.min() <<
         " | max: " << timer.max() <<
         " | variance: " << timer.variance();
   return os;
}

//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename TP >  // Element type of vector
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const Timer<TP> & t )
{
   buf.addDebugMarker( "ti" );
   buf << t.counter_;
   buf << t.start_;
   buf << t.end_;
   buf << t.time_;
   buf << t.sumOfSquares_;
   buf << t.min_;
   buf << t.max_;
   buf << t.last_;
   return buf;
}

template< typename T,    // Element type  of RecvBuffer
          typename TP >  // Element type of vector
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, Timer<TP> & t )
{
   buf.readDebugMarker( "ti" );
   buf >> t.counter_;
   buf >> t.start_;
   buf >> t.end_;
   buf >> t.time_;
   buf >> t.sumOfSquares_;
   buf >> t.min_;
   buf >> t.max_;
   buf >> t.last_;
   return buf;
}

} //namespace timing

typedef timing::Timer<timing::CpuPolicy>  CpuTimer;
typedef timing::Timer<timing::WcPolicy>    WcTimer;

} // namespace walberla
