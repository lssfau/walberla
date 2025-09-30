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
//! \file TimingPool.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "Timer.h"
#include "ReduceType.h"

#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"
#include "core/Abort.h"
#include "core/logging/Logging.h"
#include "core/mpi/SetReduction.h"
#include "core/mpi/MPIManager.h"

#include <iostream>
#include <map>
#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>


namespace walberla {
namespace timing {

template< typename TP > // Timing policy
class ScopeTimer;


/***********************************************************************************************************************
 * \brief Collection of named timers, with communication functionality.
 *
 * \ingroup timing
 *
 **********************************************************************************************************************/
template< typename TP >  // Timing policy
class TimingPool
{
   template< typename T,     // Element type of SendBuffer
             typename G,     // Growth policy of SendBuffer
             typename TP2 >  // Element type of vector
   friend mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const TimingPool<TP2> & tp );

   template< typename T,     // Element type  of RecvBuffer
             typename TP2 >  // Element type of vector
   friend mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, TimingPool<TP2> & tp );
public:

   //**Construction & Destruction***************************************************************************************
   /*! \name Construction & Destruction */
   //@{
   TimingPool() = default;
   //@}
   //*******************************************************************************************************************

   //** Access to Timer Objects ****************************************************************************************
   /*! \name Access to Timer Objects */
   //@{
   using iterator = typename std::map<std::string, Timer<TP>>::iterator;
   using const_iterator = typename std::map<std::string, Timer<TP>>::const_iterator;

         iterator begin()       { return timerMap_.begin(); }
   const_iterator begin() const { return timerMap_.begin(); }
         iterator end()         { return timerMap_.end();   }
   const_iterator end()   const { return timerMap_.end();   }

   inline       Timer<TP> & operator[]   ( const std::string & n )      ;
   inline const Timer<TP> & operator[]   ( const std::string & n ) const;
   inline bool              timerExists  ( const std::string & n ) const;
   inline bool              registerTimer( const std::string & n );
   ScopeTimer<TP>           getScopeTimer( const std::string & n );
   //@}
   //*******************************************************************************************************************

   //** Reduction ******************************************************************************************************
   /*! \name Reduction */
   //@{
   shared_ptr<TimingPool<TP> > getReduced ( ReduceType rt = REDUCE_TOTAL, int targetWorldRank = 0 ) const;
   //@}
   //*******************************************************************************************************************


   //** Utility ********************************************************************************************************
   /*! \name Utility */
   //@{
   void print       ( std::ostream & os ) const;
   void printMatlab ( std::ostream & os ) const;
   void printCSV    ( std::ostream & os, const std::string& separator = ", " ) const;

   void unifyRegisteredTimersAcrossProcesses();
   void logResultOnRoot ( ReduceType rt = REDUCE_TOTAL, const bool unifyRegisteredTimers = false );
   void merge ( const TimingPool<TP> & tpToMerge, bool mergeDuplicates = true );
   void clear ();
   bool empty() const { return timerMap_.empty(); }
   //@}
   //*******************************************************************************************************************

protected:

   //** Reduction Helper ***********************************************************************************************
   /*! \name Reduction Helper */
   //@{

   void mpiReduce( std::vector<double> & min,
                   std::vector<double> & max,
                   std::vector<double> & val,
                   std::vector<double> & valSq,
                   std::vector<uint32_t> & count,
                   int targetRank,
                   shared_ptr<TimingPool<TP> > & out ) const;
   //@}
   //*******************************************************************************************************************



   std::map<std::string, Timer<TP> > timerMap_;
};



//======================================================================================================================
//
//  REDUCTION
//
//======================================================================================================================

//**********************************************************************************************************************
/*! Returns a reduced TimingPool, holding information from all processes
 *
 * This function has to be called on all processes. The TimingPool has to have the same number
 * of timers on each processor, and the timers have to have the same name.
 *
 * \param rt         Specified the method how the reduction is done. See documentation for ReduceType
 * \param targetRank the world rank of the target process. Or negative value for an all-reduction
 *                   operation
 *
 * \return  a nonzero shared pointer to a timing pool on the process with rank "targetRank"
 *          and a zero pointer otherwise.
 *          If targetRank < 0, a valid pointer is returned on all processes (MPI_Allreduce operation)
 **********************************************************************************************************************/
template< typename TP >
shared_ptr<TimingPool<TP> > TimingPool<TP>::getReduced( ReduceType rt, int targetRank ) const
{
   WALBERLA_NON_MPI_SECTION() {
      return make_shared<TimingPool<TP> >( *this );
   }

   if( timerMap_.empty() )
      return make_shared<TimingPool<TP> >();

   std::vector<double> min;
   std::vector<double> max;
   std::vector<double> vals;
   std::vector<double> valsSq;
   std::vector<uint32_t> count;
   vals.reserve   ( timerMap_.size() );
   valsSq.reserve ( timerMap_.size() );

   shared_ptr<TimingPool<TP> > ret;
   switch (rt)
   {
   case REDUCE_MIN  :
      for ( const auto &timer : timerMap_ | std::views::values ) {
         const double val = timer.min();
         vals  .push_back( val     );
         valsSq.push_back( val*val );
         count .push_back( 1 );
      }
      mpiReduce ( vals, vals, vals, valsSq, count, targetRank, ret );
      break;
   case REDUCE_AVG  :
      for ( const auto &timer : timerMap_ | std::views::values ) {
         const double val = timer.average();
         vals  .push_back( val     );
         valsSq.push_back( val*val );
         count .push_back( 1 );
      }
      mpiReduce ( vals, vals, vals, valsSq, count, targetRank, ret  );
      break;

   case REDUCE_MAX  :
      for ( const auto &timer : timerMap_ | std::views::values ) {
         const double val = timer.max();
         vals  .push_back( val     );
         valsSq.push_back( val*val );
         count .push_back( 1 );
      }
      mpiReduce ( vals, vals, vals, valsSq, count, targetRank, ret  );
      break;

   case REDUCE_TOTAL:
      min.reserve ( timerMap_.size() );
      max.reserve ( timerMap_.size() );

      for ( const auto &timer : timerMap_ | std::views::values ) {
         vals  .push_back( timer.total() );
         valsSq.push_back( timer.sumOfSquares() );
         min   .push_back( timer.min() );
         max   .push_back( timer.max() );
         count .push_back( uint32_c(timer.getCounter()) );
      }
      mpiReduce ( min, max, vals, valsSq, count, targetRank, ret );
      break;

   default:
      WALBERLA_ABORT( "Unknown reduce type" )
      break;
   }

   return ret;
}



//**********************************************************************************************************************
/*! Executes an mpi reduction of the given vectors, and stores the result in a new TimingPool
 **********************************************************************************************************************/
template< typename TP >
void TimingPool<TP>::mpiReduce( std::vector<double> & min,
                                 std::vector<double> & max,
                                 std::vector<double> & val,
                                 std::vector<double> & valSq,
                                 std::vector<uint32_t> & count,
                                 int targetRank,
                                 shared_ptr<TimingPool<TP> > & out ) const
{
   WALBERLA_ASSERT_EQUAL( min.size(), timerMap_.size() );
   WALBERLA_ASSERT_EQUAL( max.size(), timerMap_.size() );
   WALBERLA_ASSERT_EQUAL( val.size(), timerMap_.size() );
   WALBERLA_ASSERT_EQUAL( valSq.size(), timerMap_.size() );
   WALBERLA_ASSERT_EQUAL( count.size(), timerMap_.size() );

   // Target vectors where reduced values are stored
   std::vector<double> minRed;
   std::vector<double> maxRed;
   std::vector<double> sumRed;
   std::vector<double> sumSqRed;
   std::vector<uint32_t> countRed;

   const int rank = MPIManager::instance()->worldRank();

   if ( targetRank < 0 || targetRank == rank )
   {
      minRed.  resize( timerMap_.size() );
      maxRed.  resize( timerMap_.size() );
      sumRed.  resize( timerMap_.size() );
      sumSqRed.resize( timerMap_.size() );
      countRed.resize( timerMap_.size() );
   }

   if( targetRank >= 0 )
   {
      void * minTarget = targetRank == rank ? &minRed.front() : nullptr;
      void * maxTarget = targetRank == rank ? &maxRed.front() : nullptr;
      void * sumTarget = targetRank == rank ? &sumRed.front() : nullptr;
      void * sumSqTarget = targetRank == rank ? &sumSqRed.front() : nullptr;
      void * countTarget = targetRank == rank ? &countRed.front() : nullptr;

      MPI_Reduce( &min.front(), minTarget,
                 int_c(min.size()), MPITrait<double>::type(), MPI_MIN, targetRank,MPI_COMM_WORLD );

      MPI_Reduce( &max.front(), maxTarget,
                 int_c(max.size()), MPITrait<double>::type(), MPI_MAX, targetRank,MPI_COMM_WORLD );

      MPI_Reduce( &val.front(), sumTarget,
                 int_c(val.size()), MPITrait<double>::type(), MPI_SUM, targetRank,MPI_COMM_WORLD );

      MPI_Reduce( &valSq.front(), sumSqTarget,
                 int_c(valSq.size()), MPITrait<double>::type(), MPI_SUM, targetRank,MPI_COMM_WORLD );

      MPI_Reduce( &count.front(), countTarget,
                 int_c(count.size()), MPITrait<uint32_t>::type(), MPI_SUM,targetRank,MPI_COMM_WORLD );
   }
   else
   {
      MPI_Allreduce( &min.front(), &minRed.front(),
                    int_c(min.size()), MPITrait<double>::type(), MPI_MIN, MPI_COMM_WORLD );

      MPI_Allreduce( &max.front(), &maxRed.front(),
                    int_c(max.size()), MPITrait<double>::type(), MPI_MAX, MPI_COMM_WORLD );

      MPI_Allreduce( &val.front(), &sumRed.front(),
                    int_c(val.size()), MPITrait<double>::type(), MPI_SUM, MPI_COMM_WORLD );

      MPI_Allreduce( &valSq.front(), &sumSqRed.front(),
                    int_c(valSq.size()), MPITrait<double>::type(), MPI_SUM, MPI_COMM_WORLD );

      MPI_Allreduce( &count.front(), &countRed.front(),
                    int_c(count.size()), MPITrait<uint32_t>::type(), MPI_SUM, MPI_COMM_WORLD );
   }



   // On root the timing map is replaced by reduced timers
   if ( targetRank < 0 || targetRank == MPIManager::instance()->worldRank() )
   {
      out = make_shared<TimingPool<TP> > ();

      uint_t idx = 0;
      for( const auto &name : timerMap_ | std::views::keys )
      {
         (*out)[name] = Timer<TP>( countRed[idx],
                                   minRed[idx],
                                   maxRed[idx],
                                   sumRed[idx],
                                   sumSqRed[idx] );
         idx++;
      }
   }

}



//======================================================================================================================
//
//  UTILITY
//
//======================================================================================================================


//**********************************************************************************************************************
/*! Prints the timing pool to the given ostream. No reduction is done!
 *
 * If the result of all processes should be printed, first a reduction has to be done:
 * \code
     tp.print(cout); // prints only the local time measurements
     shared_ptr< WcTimingPool> reduced = tp.getReduced(REDUCE_TOTAL, 0);
     WALBERLA_ROOT_SECTION () {
        reduced->print(cout);
     }
   \endcode
 *
 **********************************************************************************************************************/
template< typename TP >
void TimingPool<TP>::print( std::ostream & os ) const
{
   using std::endl;
   using std::setw;
   using std::fixed;
   using std::left;
   using std::right;

   int firstColumn = 0;
   double totalTime = 0.0;
   for ( const auto &[name, timer] : timerMap_ )
   {
      if ( int_c(name.size()) > firstColumn)
         firstColumn = int_c( name.size() );

      totalTime += timer.total();
   }

   firstColumn += 3;

   const int OTHER_COLUMNS  = 15;
   const int PERCENT_COLUMN = 9;
   const int NR_OF_COLUMNS  = 6;
   const int LINE_WIDTH     = firstColumn +  PERCENT_COLUMN + NR_OF_COLUMNS * OTHER_COLUMNS + NR_OF_COLUMNS + 2;

   os << setw(firstColumn)    << std::left  << "Timer"      << "|";
   os << setw(PERCENT_COLUMN) << std::right << "% "         << "|";
   os << setw(OTHER_COLUMNS)  << std::right << "Total"      << "|";
   os << setw(OTHER_COLUMNS)  << std::right << "Average"    << "|";
   os << setw(OTHER_COLUMNS)  << std::right << "Count"      << "|";
   os << setw(OTHER_COLUMNS)  << std::right << "Min"        << "|";
   os << setw(OTHER_COLUMNS)  << std::right << "Max"        << "|";
   os << setw(OTHER_COLUMNS)  << std::right << "Variance"   << "|";

   os << endl;

   for(int i=0; i< LINE_WIDTH; ++i)
      os << "-";

   os << endl;

   for ( const auto &[name, timer] : timerMap_ )
   {
      const uint_t count = timer.getCounter();
      const double percentage = ( count == uint_t(0) ) ? 0.0 : ( timer.total() / totalTime * 100.0 );

      os << setw(firstColumn-1) << std::left  << name  << " "       << "|";
      os << setw(PERCENT_COLUMN-2) << std::right << std::fixed << std::setprecision(2) << percentage << "% |";
      os << setw(OTHER_COLUMNS) << std::right << ( ( count == uint_t(0) ) ? 0.0 : timer.total() )    << "|";
      os << setw(OTHER_COLUMNS) << std::right << std::setprecision(3) << ( ( count == uint_t(0) ) ? 0.0 : timer.average() )  << "|";
      os << setw(OTHER_COLUMNS) << std::right << timer.getCounter() << "|";
      os << setw(OTHER_COLUMNS) << std::right << ( ( count == uint_t(0) ) ? 0.0 : timer.min() )      << "|";
      os << setw(OTHER_COLUMNS) << std::right << ( ( count == uint_t(0) ) ? 0.0 : timer.max() )      << "|";
      os << setw(OTHER_COLUMNS) << std::right << ( ( count == uint_t(0) ) ? 0.0 : timer.variance() ) << "|";
      os << endl;
   }
}


//**********************************************************************************************************************
/*! Prints the TimingPool in a format that can be read  by MATLAB for plotting
 *
 **********************************************************************************************************************/
template< typename TP >
void TimingPool<TP>::printMatlab( std::ostream & os ) const
{
   os << "name = { ";
   for ( auto i = timerMap_.begin(); i != timerMap_.end(); ++i )
      os <<  (i == timerMap_.begin() ? "'" : ", '" ) <<  i->first << "'";
   os << "}; " << std::endl;

   os << "average = [ ";
   for ( auto i = timerMap_.begin(); i != timerMap_.end(); ++i )
      os <<  (i == timerMap_.begin() ? "" : "," ) <<  i->second.average();
   os << "]; " << std::endl;

   os << "count = [ ";
   for ( auto i = timerMap_.begin(); i != timerMap_.end(); ++i )
      os <<  (i == timerMap_.begin() ? "" : "," ) <<  i->second.getCounter();
   os << "]; " << std::endl;


   os << "count = [ ";
   for ( auto i = timerMap_.begin(); i != timerMap_.end(); ++i )
      os <<  (i == timerMap_.begin() ? "" : "," ) <<  i->second.getCounter();
   os << "]; " << std::endl;

   os << "min = [ ";
   for ( auto i = timerMap_.begin(); i != timerMap_.end(); ++i )
      os <<  (i == timerMap_.begin() ? "" : "," ) <<  i->second.min();
   os << "]; " << std::endl;

   os << "max = [ ";
   for ( auto i = timerMap_.begin(); i != timerMap_.end(); ++i )
      os <<  (i == timerMap_.begin() ? "" : "," ) <<  i->second.max();
   os << "]; " << std::endl;

   os << "variance = [ ";
   for ( auto i = timerMap_.begin(); i != timerMap_.end(); ++i )
      os <<  (i == timerMap_.begin() ? "" : "," ) <<  i->second.variance();
   os << "]; " << std::endl;

   os << "total = [ ";
   for ( auto i = timerMap_.begin(); i != timerMap_.end(); ++i )
      os <<  (i == timerMap_.begin() ? "" : "," ) <<  i->second.total();
   os << "]; " << std::endl;
}

//**********************************************************************************************************************
/*! Prints the TimingPool in a format that can be read  by MATLAB for plotting
 *
 **********************************************************************************************************************/
template< typename TP >
void TimingPool<TP>::printCSV    ( std::ostream & os, const std::string& separator ) const
{
   os << "#";
   for ( auto i = timerMap_.begin(); i != timerMap_.end(); )
   {
      os << i->first;
      ++i;
      if (i != timerMap_.end())
         os << separator;
   }
   os << std::endl;
   for ( auto i = timerMap_.begin(); i != timerMap_.end(); )
   {
      os << i->second.average() << separator << i->second.variance();
      ++i;
      if (i != timerMap_.end())
         os << separator;
   }
   os << std::endl;
}



//**********************************************************************************************************************
/*! Ensures that on each process the same timers (= timers with identical names) are registered.
 *
 * Has to be called by all processes, not just on root process.
 **********************************************************************************************************************/
template< typename TP >
void TimingPool<TP>::unifyRegisteredTimersAcrossProcesses()
{
   std::vector< std::string > names;
   for( auto timer = begin(); timer != end(); ++timer )
      names.push_back( timer->first );

   auto gatheredNames = mpi::allReduceSet( names, mpi::UNION );
   for(auto & gatheredName : gatheredNames)
      if( !timerExists(gatheredName) )
         registerTimer( gatheredName );
}



//**********************************************************************************************************************
/*! Reduces the timing pool, and logs reduced result on root.
 *
 * Has to be called by all processes, not just on root process.
 **********************************************************************************************************************/
template< typename TP >
void TimingPool<TP>::logResultOnRoot( ReduceType rt, const bool unifyRegisteredTimers )
{
   if( unifyRegisteredTimers )
      unifyRegisteredTimersAcrossProcesses();

   auto reduced = this->getReduced( rt );
   WALBERLA_ROOT_SECTION() {
      WALBERLA_LOG_RESULT( *reduced );
   }
}



//**********************************************************************************************************************
/*! Merges a second TimingPool into this TimingPool
 *
 *\param tpToMerge        the second TimingPool, which is merged
 *\param mergeDuplicates  parameter has an effect when both TimingPool's have Timers with same name
 *                        if true:  the timers with the same name are merged together
 *                        if false: the timer in this TimingPool is not changed, the content
 *                                  of the second timer with the same name is neglected
 *
 **********************************************************************************************************************/
template< typename TP >
void TimingPool<TP>::merge ( const TimingPool<TP> & tpToMerge, bool mergeDuplicates )
{
   if ( mergeDuplicates )
   {
      for( const auto &timerName : timerMap_ | std::views::keys )
      {
         if ( timerExists( timerName ) )
         {
            // Merge the timer to the existing
            const Timer<TP> & t1 =   (*this)[timerName];
            const Timer<TP> & t2 = tpToMerge[timerName];
            (*this)[timerName] = Timer<TP> ( t1.getCounter()   + t2.getCounter(),
                                                std::min( t1.min(), t2.min() ),
                                                std::max( t1.max(), t2.max() ),
                                                t1.total()        + t2.total(),
                                                t1.sumOfSquares() + t2.sumOfSquares() );
         }
      }
   }

   // Insert skips already existing timers
   timerMap_.insert( tpToMerge.timerMap_.begin(), tpToMerge.timerMap_.end() );
}



//**********************************************************************************************************************
/*! Resets all timers
 */
//**********************************************************************************************************************
template< typename TP >
void TimingPool<TP>::clear ()
{
   for( auto &timer : timerMap_ | std::views::values ) {
      timer.reset();
   }
}


//======================================================================================================================
//
//  ACCESS TO TIMER OBJECTS
//
//======================================================================================================================

template< typename TP >
inline const Timer<TP> & TimingPool<TP>::operator[] ( const std::string & n )  const
{
   if(!timerExists(n))
      throw std::out_of_range("No such timer");
   return (timerMap_.find(n))->second;
}

template< typename TP >
inline Timer<TP> & TimingPool<TP>::operator[] ( const std::string & n )
{
   return timerMap_[n];
}


template< typename TP >
inline bool TimingPool<TP>::timerExists( const std::string & n ) const
{
   return timerMap_.find(n) != timerMap_.end();
}

template< typename TP >
inline bool TimingPool<TP>::registerTimer( const std::string & n )
{
   // timer is started in constructor Timer<TP>()
   // - but doesn't matter since a second call to start replaces the value
   //  only when end() is called a new timing is added
   auto ret = timerMap_.insert( std::pair<std::string, Timer<TP> >( n, Timer<TP>() ) );

   return ret.second; // false if element already existed
}


//======================================================================================================================
//
//  SCOPE TIMER
//
//======================================================================================================================


template< typename TP > // Timing policy
class ScopeTimer
{
public:
   ScopeTimer( Timer<TP> & timer ) : timer_( timer ) { timer_.start(); }
   ~ScopeTimer() { timer_.end(); }
private:
   Timer<TP> & timer_;
};

template< typename TP >
inline ScopeTimer<TP> TimingPool<TP>::getScopeTimer( const std::string & n )
{
   return ScopeTimer<TP>( ( *this )[n] );
}


//======================================================================================================================
//
//  OSTREAM OVERLOAD
//
//======================================================================================================================

template< typename TP >  // Timing policy
std::ostream & operator<< ( std::ostream & os, const TimingPool<TP> & tp ) {
   tp.print(os);
   return os;
}

} // namespace timing
} // namespace walberla

//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace timing {

   template< typename T,    // Element type of SendBuffer
             typename G,    // Growth policy of SendBuffer
             typename TP >  // Element type of vector
   mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const timing::TimingPool<TP> & tp )
   {
      buf.addDebugMarker( "tp" );
      buf << tp.timerMap_;
      return buf;
   }

   template< typename T,    // Element type  of RecvBuffer
             typename TP >  // Element type of vector
   mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, timing::TimingPool<TP> & tp )
   {
      buf.readDebugMarker( "tp" );
      buf >> tp.timerMap_;
      return buf;
   }

} //namespace timing
} //namespace walberla


//======================================================================================================================
//
//  EXPORTS
//
//======================================================================================================================

namespace walberla {
   using WcTimingPool = timing::TimingPool<timing::WcPolicy>;
   using CpuTimingPool = timing::TimingPool<timing::CpuPolicy>;
}
