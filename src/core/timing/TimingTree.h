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
//! \file TimingTree.h
//! \ingroup core
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "TimingNode.h"

#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"

#include <algorithm>
#include <iostream>
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
class TimingTree
{
public:
   /// Creates and initialises the timing structure
   TimingTree();
   TimingTree(const TimingTree& tt);
   TimingTree<TP>& operator=(const TimingTree<TP>& tt);

   void swap(TimingTree<TP>& tt);

   /// Starts a timer beneath the current hierarchy level
   void start(const std::string& name);
   /// Stops the last started timer and jumps back one hierarchy level
   void stop(const std::string& name);
   /// Checks if specified timer is currently running.
   bool isTimerRunning(const std::string& name) const;
   /// Resets the the timing hierarchy
   void reset();

   /// Collects all the timing data from different processes
   ///
   /// \param rt type of reduction used
   /// \param targetRank
   /// \param callSynchronize call synchronize() before doing the reduction
   /// \return
   ///
   /// \attention Setting callSynchronize to false can lead to deadlocks during reduction if different ranks
   ///            have different timing nodes! Setting it to true causes an additional communication.
   TimingTree<TP> getReduced( ReduceType rt = REDUCE_TOTAL, int targetRank = 0, bool callSynchronize = true ) const;

   /// Adds entries which only exist on other processes. Has to be collectively called on all processes!
   void synchronize();

   /// Returns the raw tree data structure
   const TimingNode<TP>& getRawData() const;

   const Timer<TP>& operator[](const std::string& name) const;
   inline bool timerExists ( const std::string & n ) const;

   /// Returns the name of the currently running timer
   /// Might be expensive due to value search.
   std::string getCurrentTimerName() const;

   /// Returns a copy of the timing tree containing the remaining time as a subnode
   TimingTree< TP > getCopyWithRemainder() const;

private:
   /// Tree data structure
   TimingNode<TP>  root_;
   /// Pointer to the current hierarchy level.
   TimingNode<TP>* current_;
};

/// \relates TimingTree
template< typename TP >  // Timing policy
std::ostream& operator<<( std::ostream& os, const TimingTree<TP>& tt )
{
   os << tt.getRawData();
   return os;
}

template< typename TP >  // Timing policy
TimingTree<TP>::TimingTree()
      : current_(&root_)
{

}

template< typename TP >  // Timing policy
TimingTree<TP>::TimingTree(const TimingTree<TP>& tt)
      : root_(tt.root_)
      , current_(&root_)
{
   WALBERLA_ASSERT_EQUAL(tt.current_, &tt.root_, "Copying is only allowed for stopped TimingTrees!\nTimer still running: " << getCurrentTimerName() );
}

template< typename TP >  // Timing policy
TimingTree<TP>& TimingTree<TP>::operator=(const TimingTree<TP>& tt)
{
   TimingTree<TP> tmp (tt);
   tmp.swap(*this);
   return *this;
}

template< typename TP >  // Timing policy
void TimingTree<TP>::swap(TimingTree<TP>& tt)
{
   std::swap(current_, tt.current_);
   std::swap(root_, tt.root_);
}

/// \param name timer name. '.' is not allowed in the timer name!
template< typename TP >  // Timing policy
void TimingTree<TP>::start(const std::string& name)
{
   if (name.find_first_of(".") != std::string::npos)
   {
      WALBERLA_LOG_WARNING("'.' not allowed in timer name!");
   }
   auto tmp = current_;
   current_ = &(current_->tree_[name]);
   current_->last_ = tmp; //if node is created by previous call last_ is not set
   current_->timer_.start();
//   WALBERLA_LOG_DEVEL("Timer started: " << name);
}

template< typename TP >  // Timing policy
void TimingTree<TP>::stop(const std::string& name)
{
   if (name.find_first_of(".") != std::string::npos)
   {
      WALBERLA_LOG_WARNING("'.' not allowed in timer name!");
   }
   WALBERLA_ASSERT_NOT_NULLPTR( current_->last_ );
   auto timerIt = current_->last_->tree_.find(name);
   if ((timerIt == current_->last_->tree_.end()) || (&(timerIt->second) != current_))
   {
      WALBERLA_LOG_WARNING("Trying to stop timer which is not running: " << name << "\nCurrently Running: " << getCurrentTimerName() );
      const auto nan = std::numeric_limits<double>::quiet_NaN();
      current_->tree_[name].timer_ = Timer<TP>(1, nan, nan, nan, nan);
   } else
   {
      current_->timer_.end();
      current_ = current_->last_;
//      WALBERLA_LOG_DEVEL("Timer stopped: " << name);
   }
}

template< typename TP >  // Timing policy
bool TimingTree<TP>::isTimerRunning(const std::string& name) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( current_->last_ );
   auto timerIt = current_->last_->tree_.find(name);
   return !((timerIt == current_->last_->tree_.end()) || (&(timerIt->second) != current_));
}

template< typename TP >  // Timing policy
void TimingTree<TP>::reset()
{
   walberla::timing::reset(root_);
}

template< typename TP >  // Timing policy
TimingTree<TP> TimingTree<TP>::getReduced( ReduceType rt, int targetRank, bool callSynchronize ) const
{
   TimingTree<TP> tt(*this);
   if (callSynchronize)
   {
      tt.synchronize();
   }
   reduceInplace( tt.root_, rt, targetRank );
   return tt;
}

template< typename TP >  // Timing policy
void TimingTree<TP>::synchronize()
{
   synchronizeEntries( root_ );
}

template< typename TP >  // Timing policy
const TimingNode<TP>& TimingTree<TP>::getRawData() const
{
   return root_;
}

/// Finds the specified timer in the timing hierarchy
/// \param name timer name which may include more than one hierarchy separated by "."
/// \code tt["firstLevel.secondLevel.thirdLevel.timerName"].total(); \endcode
template< typename TP >  // Timing policy
const Timer<TP>& TimingTree<TP>::operator[](const std::string& name) const
{
   return findTimer(root_, name);
}

/// Checks if the specified timer exists in the timing hierarchy
/// \param name timer name which may include more than one hierarchy separated by "."
/// \code tt.timerExists("firstLevel.secondLevel.thirdLevel.timerName"); \endcode
template< typename TP >  // Timing policy
bool TimingTree<TP>::timerExists(const std::string& name) const
{
   return walberla::timing::timerExists(root_, name);
}

template< typename TP >  // Timing policy
std::string TimingTree<TP>::getCurrentTimerName() const
{
   if (current_ == current_->last_)
   {
      return "No timer running!";
   }
   for (auto it = current_->last_->tree_.begin(); it != current_->last_->tree_.end(); ++it)
   {
      if (&(it->second) == current_)
         return it->first;
   }
   return "No timer found!";
}

template < typename TP > // Timing policy
TimingTree< TP > TimingTree< TP >::getCopyWithRemainder() const
{
   TimingTree< TP > tt( *this );
   timing::internal::addRemainderNodes< TP >( tt.root_ );
   return tt;
}

}

using WcTimingTree = timing::TimingTree<timing::WcPolicy>;
using CpuTimingTree = timing::TimingTree<timing::CpuPolicy>;
}
