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
//! \file TimingNode.h
//! \ingroup core
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "Timer.h"
#include "ReduceType.h"

#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/mpi/SetReduction.h"


#include <algorithm>
#include <iomanip>
#include <map>
#include <string>


namespace walberla {
namespace timing {

/***********************************************************************************************************************
 * \brief Hierarchical structure of timers.
 *
 * Timer class to time different code snippets. Nested timers can be created and will be output
 * in a tree like structure. Also works MPI parallel by using the reduce function but is NOT
 * threadsafe!
 *
 * \ingroup timing
 *
 **********************************************************************************************************************/
template< typename TP >  // Timing policy
struct TimingNode
{
   /// Creates and initialises the timing structure
   TimingNode();
   TimingNode(const TimingNode& tt);
   TimingNode<TP>& operator=(const TimingNode<TP>& tt);

   void swap(TimingNode<TP>& tt);

   /// Implementation of the recursive print function
   void printImpl(std::ostream & os, const std::string& prefix, const double totalTime, const double parentTime, const int firstColumn) const;

   /// Pointer to the parent node.
   TimingNode* last_;

   /// Node timer
   Timer<TP> timer_;
   /// Collection of sub timers
   std::map<std::string, TimingNode> tree_;
};

template< typename TP >  // Timing policy
TimingNode<TP>::TimingNode()
   : last_(nullptr)
{

}

template< typename TP >  // Timing policy
TimingNode<TP>::TimingNode(const TimingNode<TP>& tn)
   : last_(nullptr)
   , timer_(tn.timer_)
   , tree_(tn.tree_)
{
   // reconstruct hierarchy
   for (auto it = tree_.begin(); it != tree_.end(); ++it)
   {
      it->second.last_ = this;
   }
}

template< typename TP >  // Timing policy
TimingNode<TP>& TimingNode<TP>::operator=(const TimingNode<TP>& tn)
{
   TimingNode<TP> tmp (tn);
   tmp.swap(*this);
   return *this;
}

template< typename TP >  // Timing policy
void TimingNode<TP>::swap(TimingNode<TP>& tt)
{
    std::swap(tree_, tt.tree_);
    // reconstruct hierarchy
    for (auto it = tree_.begin(); it != tree_.end(); ++it)
    {
       it->second.last_ = this;
    }
}

/// Finds the specified timer in the timing hierarchy
/// \param name timer name which may include more than one hierarchy separated by "."
/// \code findTimer(tn, "firstLevel.secondLevel.thirdLevel.timerName"); \endcode
/// \relates TimingNode
template< typename TP >  // Timing policy
const Timer<TP>& findTimer( const TimingNode<TP>& tn, const std::string& name)
{
   auto pos = name.find_first_of('.');
   if (pos != std::string::npos)
   {
      WALBERLA_ASSERT_UNEQUAL( tn.tree_.find(name.substr(0, pos)), tn.tree_.end(), "Could not find timer: " << name.substr(0, pos) );
      return findTimer( tn.tree_.at(name.substr(0, pos)), name.substr(pos+1, std::string::npos));
   } else
   {
      WALBERLA_ASSERT_UNEQUAL( tn.tree_.find(name), tn.tree_.end(), "Could not find timer: " << name );
      return tn.tree_.at(name).timer_;
   }
}

/// Checks if the specified timer exists in the timing hierarchy
/// \param name timer name which may include more than one hierarchy separated by "."
/// \code timerExists(tn, "firstLevel.secondLevel.thirdLevel.timerName"); \endcode
/// \relates TimingNode
template< typename TP >  // Timing policy
bool timerExists( const TimingNode<TP>& tn, const std::string& name )
{
   auto pos = name.find_first_of('.');
   if (pos != std::string::npos)
   {
      if( tn.tree_.find(name.substr(0, pos)) != tn.tree_.end() )
      {
         return timerExists( tn.tree_.at(name.substr(0, pos)), name.substr(pos+1, std::string::npos));
      }
      else
      {
         return false;
      }
   }
   else
   {
      return tn.tree_.find(name) != tn.tree_.end();
   }
}

/// Resets the timer in the TimingNode structure and all sub timers
/// \relates TimingNode
template< typename TP >  // Timing policy
void reset( TimingNode<TP>& tn)
{
   tn.timer_.reset();

   for (auto it = tn.tree_.begin(); it != tn.tree_.end(); ++it)
   {
      reset(it->second);
   }
}

/// Utility function to find the length of the longest timer name
/// \relates TimingNode
template< typename TP >  // Timing policy
size_t getLongestTimerNameLength(const TimingNode<TP>& tn)
{
   size_t max = 0;
   for ( auto it = tn.tree_.begin(); it != tn.tree_.end(); ++it )
   {
      max = std::max( max, it->first.length() );
      max = std::max( max + 1, getLongestTimerNameLength(it->second) );
   }
   return max;
}

/// Utility function to find the depth of a TimingNode tree
/// \relates TimingNode
template< typename TP >  // Timing policy
size_t getHierarchyDepth(const TimingNode<TP>& tn)
{
   size_t depth = 1;
   size_t maxDepth = 0;
   for ( auto it = tn.tree_.begin(); it != tn.tree_.end(); ++it )
   {
      maxDepth = std::max( maxDepth, getHierarchyDepth(it->second) );
   }
   return depth + maxDepth;
}

/// Collects all the timing data from different processes
/// \attention Will overwrite the current timer data! Multiple calls will accumulate!
/// \relates TimingNode
template< typename TP >  // Timing policy
void reduceInplace( TimingNode<TP>& tn, ReduceType rt = REDUCE_TOTAL, int targetRank = 0 )
{
   if (mpi::MPIManager::instance()->numProcesses() == 1) return;
   if (tn.tree_.empty()) return;

   std::vector<double> min;
   std::vector<double> max;
   std::vector<double> vals;
   std::vector<double> valsSq;
   std::vector<uint32_t> count;
   vals.reserve   ( tn.tree_.size() );
   valsSq.reserve ( tn.tree_.size() );

   switch (rt)
   {
      case REDUCE_MIN  :
         for ( auto i = tn.tree_.begin(); i != tn.tree_.end(); ++i ) {
            const double val = i->second.timer_.min();
            vals  .push_back( val     );
            valsSq.push_back( val*val );
            count .push_back( 1 );
         }
         min = vals;
         max = vals;
         break;
      case REDUCE_AVG  :
         for ( auto i = tn.tree_.begin(); i != tn.tree_.end(); ++i ) {
            const double val = i->second.timer_.average();
            vals  .push_back( val     );
            valsSq.push_back( val*val );
            count .push_back( 1 );
         }
         min = vals;
         max = vals;
         break;

      case REDUCE_MAX  :
         for ( auto i = tn.tree_.begin(); i != tn.tree_.end(); ++i ) {
            const double val = i->second.timer_.max();
            vals  .push_back( val     );
            valsSq.push_back( val*val );
            count .push_back( 1 );
         }
         min = vals;
         max = vals;
         break;

      case REDUCE_TOTAL:
         min.reserve ( tn.tree_.size() );
         max.reserve ( tn.tree_.size() );

         for ( auto i = tn.tree_.begin(); i != tn.tree_.end(); ++i ) {
            vals  .push_back( i->second.timer_.total() );
            valsSq.push_back( i->second.timer_.sumOfSquares() );
            min   .push_back( i->second.timer_.min() );
            max   .push_back( i->second.timer_.max() );
            count .push_back( uint32_c(i->second.timer_.getCounter()) );
         }
         break;

      default:
         WALBERLA_ABORT( "Unknown reduce type" );
         break;
   }

   WALBERLA_ASSERT_EQUAL(vals.size(), valsSq.size());
   WALBERLA_ASSERT_EQUAL(vals.size(), min.size());
   WALBERLA_ASSERT_EQUAL(vals.size(), max.size());
   WALBERLA_ASSERT_EQUAL(vals.size(), count.size());

   WALBERLA_DEBUG_SECTION()
   {
      //checking if all timing trees contain the same number of elements
      std::vector<uint32_t> lens;
      lens.push_back(uint32_c(vals.size()));
      lens = mpi::allGatherv(lens);
      std::for_each(lens.begin(), lens.end(), [&](const uint32_t& v){WALBERLA_UNUSED(v); WALBERLA_ASSERT_EQUAL( v, vals.size(), "Different number of TimingTree nodes detected! All TimingTrees need to have the same timers for reduction!");});
   }

   // Target vectors where reduced values are stored
   std::vector<double> minRed;
   std::vector<double> maxRed;
   std::vector<double> sumRed;
   std::vector<double> sumSqRed;
   std::vector<uint32_t> countRed;

   const int rank = MPIManager::instance()->worldRank();

   if ( targetRank < 0 || targetRank == rank )
   {
      minRed.  resize( tn.tree_.size() );
      maxRed.  resize( tn.tree_.size() );
      sumRed.  resize( tn.tree_.size() );
      sumSqRed.resize( tn.tree_.size() );
      countRed.resize( tn.tree_.size() );
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

      MPI_Reduce( &vals.front(), sumTarget,
                  int_c(vals.size()), MPITrait<double>::type(), MPI_SUM, targetRank,MPI_COMM_WORLD );

      MPI_Reduce( &valsSq.front(), sumSqTarget,
                  int_c(valsSq.size()), MPITrait<double>::type(), MPI_SUM, targetRank,MPI_COMM_WORLD );

      MPI_Reduce( &count.front(), countTarget,
                  int_c(count.size()), MPITrait<uint32_t>::type(), MPI_SUM, targetRank,MPI_COMM_WORLD );
   }
   else
   {
      MPI_Allreduce( &min.front(), &minRed.front(),
                     int_c(min.size()), MPITrait<double>::type(), MPI_MIN, MPI_COMM_WORLD );

      MPI_Allreduce( &max.front(), &maxRed.front(),
                      int_c(max.size()), MPITrait<double>::type(), MPI_MAX, MPI_COMM_WORLD );

      MPI_Allreduce( &vals.front(), &sumRed.front(),
                     int_c(vals.size()), MPITrait<double>::type(), MPI_SUM, MPI_COMM_WORLD );

      MPI_Allreduce( &valsSq.front(), &sumSqRed.front(),
                     int_c(valsSq.size()), MPITrait<double>::type(), MPI_SUM, MPI_COMM_WORLD );

      MPI_Allreduce( &count.front(), &countRed.front(),
                     int_c(count.size()), MPITrait<uint32_t>::type(), MPI_SUM, MPI_COMM_WORLD );
   }

   // On root the timing map is replaced by reduced timers

   uint_t idx = 0;
   for( auto it = tn.tree_.begin(); it != tn.tree_.end(); ++it )
   {
      if ( targetRank < 0 || targetRank == MPIManager::instance()->worldRank() )
      {
         it->second.timer_ = Timer<TP>( countRed[idx], minRed[idx], maxRed[idx], sumRed[idx], sumSqRed[idx] );
      }
      reduceInplace( it->second, rt, targetRank );
      idx++;
   }
}



/// Makes sure all nodes on all processes have the same children
/// \relates TimingNode
template< typename TP >  // Timing policy
void synchronizeEntries( TimingNode<TP>& tn )
{
   std::vector<std::string> childNames;
   for( auto i = tn.tree_.begin(); i != tn.tree_.end(); ++i )
   {
      childNames.push_back( i->first );
   }
   
   std::vector<std::string> globalChildNames = mpi::allReduceSet( childNames, mpi::UNION );
   
   std::vector<std::string> missingChildNames( globalChildNames.size() - childNames.size() );
   
   std::set_difference( globalChildNames.begin(), globalChildNames.end(), childNames.begin(), childNames.end(), missingChildNames.begin() );

   for( auto it = missingChildNames.begin(); it != missingChildNames.end(); ++it )
   {
      tn.tree_[*it].last_ = &tn; // insert missing child and setup correct pointer to parent node
   }

   // recurse into children
   for( auto it = tn.tree_.begin(); it != tn.tree_.end(); ++it )
   {
      synchronizeEntries( it->second );
   }
}

template< typename TP >  // Timing policy
void TimingNode<TP>::printImpl(std::ostream & os, const std::string& prefix, const double totalTime, const double parentTime, const int firstColumn) const
{
   using std::endl;
   using std::setw;
   using std::fixed;
   using std::left;
   using std::right;

   const int OTHER_COLUMNS  = 15;
   const int PERCENT_COLUMN = 9;
//   const int NR_OF_COLUMNS  = 6;
//   const int LINE_WIDTH     = firstColumn +  PERCENT_COLUMN + NR_OF_COLUMNS * OTHER_COLUMNS + NR_OF_COLUMNS + 2;

   for ( auto i = tree_.begin(); i != tree_.end(); ++i )
   {
      const uint_t count = i->second.timer_.getCounter();
      const double percentageTotal  = ( count == uint_t(0) ) ? 0.0 : ( i->second.timer_.total() / totalTime * 100.0 );
      const double percentageParent = ( count == uint_t(0) ) ? 0.0 : ( i->second.timer_.total() / parentTime * 100.0 );

      os << setw(firstColumn-1)    << std::left  << prefix + i->first  << " "       << "|";

      os << setw(PERCENT_COLUMN-2) << std::right << std::fixed << std::setprecision(2) << percentageTotal << "% |";
      os << setw(PERCENT_COLUMN-2) << std::right << std::fixed << std::setprecision(2) << percentageParent << "% |";
      os << setw(OTHER_COLUMNS)    << std::right << ( ( count == uint_t(0) ) ? 0.0 : i->second.timer_.total() )    << "|";

      os << setw(OTHER_COLUMNS)    << std::right << ( ( count == uint_t(0) ) ? 0.0 : i->second.timer_.average() )  << "|";
      os << setw(OTHER_COLUMNS)    << std::right << i->second.timer_.getCounter() << "|";
      os << setw(OTHER_COLUMNS)    << std::right << ( ( count == uint_t(0) ) ? 0.0 : i->second.timer_.min() )      << "|";
      os << setw(OTHER_COLUMNS)    << std::right << ( ( count == uint_t(0) ) ? 0.0 : i->second.timer_.max() )      << "|";
      os << setw(OTHER_COLUMNS)    << std::right << ( ( count == uint_t(0) ) ? 0.0 : i->second.timer_.variance() ) << "|";
      os << endl;

      i->second.printImpl(os, " " + prefix, totalTime, i->second.timer_.total(), firstColumn);
   }
}

/// \relates TimingNode
template< typename TP >  // Timing policy
std::ostream& operator<<(std::ostream& os, const TimingNode<TP>& tn)
{
   using std::endl;
   using std::setw;
   using std::fixed;
   using std::left;
   using std::right;

   int firstColumn = 5;
   firstColumn = std::max( firstColumn, static_cast<int> (getLongestTimerNameLength(tn)) );
   double totalTime = 0.0;
   for ( auto i = tn.tree_.begin(); i != tn.tree_.end(); ++i )
   {
      totalTime += i->second.timer_.total();
   }

   firstColumn += 3;

   const int OTHER_COLUMNS  = 15;
   const int PERCENT_COLUMN = 9;
   const int NR_OF_COLUMNS  = 6;
   const int LINE_WIDTH     = firstColumn +  2*PERCENT_COLUMN + NR_OF_COLUMNS * OTHER_COLUMNS + NR_OF_COLUMNS + 2;

   os << setw(firstColumn)    << std::left  << "Timer"      << "|";
   os << setw(PERCENT_COLUMN) << std::right << "Abs % "     << "|";
   os << setw(PERCENT_COLUMN) << std::right << "Rel % "     << "|";
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

   tn.printImpl(os, "", totalTime, totalTime, firstColumn);

   return os;
}


namespace internal {
/// adds a sub timer containing the remainder of all other subtimers on the same hierarchy level
/// \relates TimingNode
template< typename TP >
// Timing policy
void addRemainderNodes(timing::TimingNode<TP> &tn) {
   if (tn.tree_.empty()) {
      return;
   }
   double remainder = tn.timer_.total();
   for (auto i = tn.tree_.begin(); i != tn.tree_.end(); ++i) {
      remainder -= i->second.timer_.total();
      addRemainderNodes(i->second);
   }
   if (tn.last_ != nullptr) {
      WALBERLA_ASSERT( tn.tree_.find("Remainder") == tn.tree_.end());
      tn.tree_["Remainder"].timer_ = timing::Timer<TP>(1, 0, 0, remainder, 0);
      tn.tree_["Remainder"].last_ = &tn;
   }
}
} /// namespace internal
}

using WcTimingNode = timing::TimingNode<timing::WcPolicy>;
using CpuTimingNode = timing::TimingNode<timing::CpuPolicy>;

}
