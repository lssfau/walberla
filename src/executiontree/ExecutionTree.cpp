//==============================================================================================================================================================
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
//! \file TaskTree.cpp
//! \ingroup executiontree
//! \author Martin Bauer <martin.bauer@fau.de>
//
//==============================================================================================================================================================


#include <sstream>
#include <iostream>
#include "core/logging/Logging.h"
#include "core/OpenMP.h"
#include "ExecutionTree.h"


namespace walberla {
namespace executiontree {

using timeloop::ITimeloop;


// --------------------------- Logging Integration of Loop node  -----------------------------------------------------------------------------------------------


class LoggingStamp : public logging::Logging::CustomStamp
{
public:
   explicit LoggingStamp( const ITimeloop & timeloop ) : timeloop_( timeloop ) {}
   std::string stamp() override
   {
      std::ostringstream oss;
      int indention;

      if( timeloop_.getNrOfTimeSteps() > 0 )
         indention = int_c( std::ceil( std::log10( real_c( timeloop_.getNrOfTimeSteps() ) ) ) );
      else if( timeloop_.getCurrentTimeStep() > 0 )
         indention = int_c( std::ceil( std::log10( real_c( timeloop_.getCurrentTimeStep() ) ) ) );
      else
         indention = 0;

      oss << std::setw( indention )
          << std::setfill(' ') << std::right << timeloop_.getCurrentTimeStep();
      return std::string("[") + oss.str() + std::string("]");
   }
   uint_t maxStampWidth() override
   {
      if( timeloop_.getNrOfTimeSteps() > 0 )
         return uint_c( std::ceil( std::log10( real_c( timeloop_.getNrOfTimeSteps() ) ) ) ) + uint_c(2);
      else if( timeloop_.getCurrentTimeStep() > 0 )
         return uint_c( std::ceil( std::log10( real_c( timeloop_.getCurrentTimeStep() ) ) ) ) + uint_c(2);
      else
         return uint_c(2);
   }
private:
   const ITimeloop & timeloop_;
};

class LoggingStampManager
{
public:
   LoggingStampManager( const shared_ptr< LoggingStamp > & stamp, const bool useCustomStamp )
           : useCustomStamp_( useCustomStamp )
   {
      if( useCustomStamp_ )
         logging::Logging::instance()->addCustomStamp( stamp );
   }
   ~LoggingStampManager()
   {
      if( useCustomStamp_ )
         logging::Logging::instance()->clearCustomStamp();
   }
private:
   const bool useCustomStamp_;
};


// --------------------------- Printing ------------------------------------------------------------------------------------------------------------------------

void printNode( std::ostream &os, const IFunctionNode &node, uint_t indentation )
{
   for ( uint_t i = 0; i < indentation; ++i )
      os << " ";

   os << node.getName() << std::endl;
   for ( auto &c : node.getChildren())
      printNode( os, *c, indentation + 4 );
}

std::ostream &operator<<( std::ostream &os, const IFunctionNode &node )
{
   printNode( os, node, 0 );
   return os;
}

// --------------------------- Node class implementation -------------------------------------------------------------------------------------------------------


EveryNth::EveryNth( const IFunctionNodePtr &node, uint_t interval, bool onFirst, uint_t startValue )
        : wrapped_( node ), interval_( interval ), onFirst_( onFirst ), calls_( startValue ) {}


void EveryNth::operator()()
{
   if ( calls_ == 0 && !onFirst_ ) {
      ++calls_;
      return;
   }

   if (( calls_ % interval_ ) == 0 )
      ( *wrapped_ )();
   ++calls_;
}

std::string EveryNth::getName() const
{
   std::stringstream ss;
   ss << "every " << interval_ << "th step:";
   return ss.str();
}


Sequence::Sequence( std::initializer_list< IFunctionNodePtr > initializerList, const std::string &name, const TimingTreePtr &timingTree, bool parallel )
        : name_( name ), timingTree_( timingTree ), parallel_( parallel )
{
   for ( auto &e : initializerList )
      children_.push_back( e );
}

void Sequence::operator()()
{
#ifdef WALBERLA_BUILD_WITH_OPENMP
   if( parallel_ )
   {
      if ( timingTree_ )
         timingTree_->start( name_ );

      int threads = int_c( children_.size() );
      #pragma omp parallel num_threads( threads )
      {

         ( *children_[ uint_c( omp_get_thread_num() ) ] )();
      }

      if ( timingTree_ )
         timingTree_->stop( name_ );

      return;
   }
#endif
   WALBERLA_UNUSED(parallel_);

   if ( timingTree_ )
      timingTree_->start( name_ );

   for ( auto &el : children_ )
   {
      ( *el )();
   }

   if ( timingTree_ )
      timingTree_->stop( name_ );
}


Loop::Loop( const IFunctionNodePtr &body, uint_t iterations, bool logTimeStep )
        : body_( body ), currentIteration_( 0 ), iterations_( iterations ), stop_( false ), logTimeStep_( logTimeStep ) {}


void Loop::singleStep()
{
   LoggingStampManager raii( make_shared<LoggingStamp>( *this ), logTimeStep_ );
   ( *body_ )();
   ++currentIteration_;
}

void Loop::operator()()
{
   LoggingStampManager raii( make_shared<LoggingStamp>( *this ), logTimeStep_ );

   for ( ; currentIteration_ < iterations_; ++currentIteration_ )
   {
      ( *body_ )();
      if ( stop_ )
      {
         stop_ = false;
         break;
      }
   }
}

void Loop::synchronizedStop( bool stopVar )
{
   stop_ = stopVar;
   mpi::allReduceInplace( stop_, mpi::LOGICAL_OR );
}

std::string Loop::getName() const
{
   std::stringstream ss;
   ss << "Loop [" << iterations_ << "]";
   return ss.str();
}


} // namespace tasktree
} // namespace walberla