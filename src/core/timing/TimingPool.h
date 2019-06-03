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

#include <iostream>
#include <map>
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
   TimingPool() {}
   //@}
   //*******************************************************************************************************************

   //** Access to Timer Objects ****************************************************************************************
   /*! \name Access to Timer Objects */
   //@{
   typedef typename std::map<std::string, Timer<TP> >::iterator       iterator;
   typedef typename std::map<std::string, Timer<TP> >::const_iterator const_iterator;

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
   typedef timing::TimingPool<timing::WcPolicy>  WcTimingPool;
   typedef timing::TimingPool<timing::CpuPolicy> CpuTimingPool;
}
