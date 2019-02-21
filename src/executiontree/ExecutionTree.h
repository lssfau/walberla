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
//! \file TaskTree.h
//! \ingroup executiontree
//! \author Martin Bauer <martin.bauer@fau.de>
//
//==============================================================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "timeloop/ITimeloop.h"
#include "core/timing/TimingTree.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <deque>
#include <string>
#include <initializer_list>
#include <functional>


namespace walberla {
namespace executiontree {


// -------------------------------------- Forward Declarations ------------------------------------------------------------------------------------------------

class IFunctionNode;
using IFunctionNodePtr = shared_ptr<IFunctionNode>;
using TimingTreePtr = shared_ptr<WcTimingTree>;

class EveryNth;
class Sequence;
class Loop;

template< typename FunctorType > class Functor;
template< typename FunctorType > class SharedFunctor;
template< typename FunctorType > class Sweep;
template< typename FunctorType > class SharedSweep;


// -------------------------------------- Public Interface     ------------------------------------------------------------------------------------------------


/*! Creates a functor node around any callable object. The wrapped functor is copied.
 *
 * \param any callable object. The object is copied - if its state has to be modified later, pass a shared_ptr to a functor instead
 * \param name optional name of the functor node
 * \param timingTree optional timing tree object to time all executions of this functor
 */
template<typename FunctorType>
IFunctionNodePtr functor( FunctorType t, const std::string &name = "", const shared_ptr< WcTimingTree > &timingTree = nullptr );


/*! Combine multiple task nodes into a (named) sequence
 *
 * \param initializerList list of tasks that are executed in the passed order
 * \param name optional sequence name, used for printing and for labeling time measurements
 * \param timingTree optional timing tree object
 */
shared_ptr< Sequence > sequence( std::initializer_list< IFunctionNodePtr > initializerList,
                                 const std::string &name = "",
                                 const TimingTreePtr &timingTree = nullptr );


/*! All subtasks of this region are executed in parallel using OpenMP */
shared_ptr< Sequence > parallelSequence( std::initializer_list< IFunctionNodePtr > initializerList,
                                         const std::string &name = "",
                                         const TimingTreePtr &timingTree = nullptr );



/*! Note that runs its contents only every n'th call
 *
 * \param node task that is only run every n'th call
 * \param name the interval i.e. "n"
 * \param onFirst if false the task is not run at the first call
 * \param startValue initial call counter
 */
shared_ptr< EveryNth > everyNth( const IFunctionNodePtr &node,
                                 uint_t interval,
                                 bool onFirst = false,
                                 uint_t startValue = 0 );

/*! Runs the child node for the given amount of iterations */
shared_ptr< Loop > loop( const IFunctionNodePtr &body, uint_t iterations, bool logTimeStep = true );

std::ostream &operator<<( std::ostream &os, const IFunctionNode &node );


// -------------------------------------- Node Classes --------------------------------------------------------------------------------------------------------


class IFunctionNode
{
public:
   virtual ~IFunctionNode() {}
   virtual void operator()() = 0;
   virtual const std::string getName() const = 0;
   virtual const std::deque< shared_ptr< IFunctionNode > > getChildren() const { return {}; }
};


template<typename FunctorType>
class Functor : public IFunctionNode
{
public:
   Functor(const FunctorType &functor,
           const std::string &name,
           const TimingTreePtr & timingTree );

   const std::string getName() const override { return name_ != "" ? name_ : "Functor"; };
   void operator() () override;

private:
   FunctorType functor_;
   std::string name_;
   shared_ptr< WcTimingTree > timingTree_;
};


class EveryNth : public IFunctionNode
{
public:
   EveryNth( const IFunctionNodePtr &node, uint_t interval, bool onFirst = false, uint_t startValue = 0 );

   void operator()() override;
   const std::string getName() const override;
   const std::deque< shared_ptr< IFunctionNode > > getChildren() const override { return { wrapped_ }; }

private:
   IFunctionNodePtr wrapped_;
   uint_t interval_;
   bool onFirst_;
   uint_t calls_;
};

class Sequence : public IFunctionNode
{
public:
   Sequence( std::initializer_list< IFunctionNodePtr > initializerList, const std::string &name,
             const TimingTreePtr &timingTree = nullptr, bool parallel = false );

   void operator()() override;

   void push_back( const IFunctionNodePtr &fct ) { children_.push_back( fct ); }
   void push_front( const IFunctionNodePtr &fct ) { children_.push_front( fct ); }
   const std::string getName() const override { return name_ != "" ? name_ : "Sequence"; };
   const std::deque< IFunctionNodePtr > getChildren() const override { return children_; };

private:
   std::string name_;
   std::deque< IFunctionNodePtr > children_;
   shared_ptr< WcTimingTree > timingTree_;
   bool parallel_;
};


class Loop : public IFunctionNode, public timeloop::ITimeloop
{
public:
   Loop( const IFunctionNodePtr &body, uint_t iterations, bool logTimeStep = true );

   void operator()() override;
   void run() override { ( *this )(); }
   void singleStep() override;

   void synchronizedStop( bool stopVal ) override;
   void stop() override { stop_ = true; }
   void setBody( const IFunctionNodePtr &body ) { body_ = body; }
   void setCurrentTimeStep( uint_t ts ) override { currentIteration_ = ts; };
   uint_t getCurrentTimeStep() const override { return currentIteration_; }
   uint_t getNrOfTimeSteps() const override { return iterations_; }

   const std::deque< shared_ptr< IFunctionNode > > getChildren() const override { return { body_ }; }
   const std::string getName()  const override;

private:
   IFunctionNodePtr body_;
   uint_t currentIteration_;
   uint_t iterations_;
   bool stop_;
   bool logTimeStep_;
};




} // namespace executiontree
} // namespace walberla


#include "ExecutionTree.impl.h"