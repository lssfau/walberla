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
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//==============================================================================================================================================================

#pragma once

#include "executiontree/ExecutionTree.h"
#include "ParallelStreams.h"

#include <cuda_runtime.h>

#ifdef CUDART_VERSION
#if CUDART_VERSION <= 9020
cudaError_t cudaLaunchHostFunc( cudaStream_t,  void(CUDART_CB* )( void*  userData ), void* ) {
        static bool printedWarning = false;
        if( ! printedWarning ) {
                WALBERLA_LOG_WARNING_ON_ROOT("Timing of CUDA functions only implemented for CUDA versions >= 10.0" );
                printedWarning = true;
        }
        return cudaSuccess;
}
#endif
#endif

namespace walberla {
namespace executiontree {

// -------------------------------------- Forward Declarations ------------------------------------------------------------------------------------------------

using executiontree::IFunctionNode;
using executiontree::IFunctionNodePtr;
using executiontree::TimingTreePtr;

class SequenceCUDA;
class IFunctionNodeCUDA;
template<typename FunctorClass> class FunctorCUDA;
using IFunctionNodeCUDAPtr = shared_ptr<IFunctionNodeCUDA>;


// -------------------------------------- Public Interface     ------------------------------------------------------------------------------------------------

template<typename FunctorType>
IFunctionNodeCUDAPtr functorCUDA( const FunctorType & t, const std::string &name = "", const TimingTreePtr &timingTree = nullptr );


shared_ptr< SequenceCUDA > sequenceCUDA( std::initializer_list< IFunctionNodeCUDAPtr > initializerList,
                                         const std::string &name, cudaStream_t defaultStream = 0, bool parallel = false, int priority = 0,
                                         const TimingTreePtr &timingTree = nullptr );


// -------------------------------------- Node Classes --------------------------------------------------------------------------------------------------------


class IFunctionNodeCUDA : public IFunctionNode
{
public:
   virtual void operator()( cudaStream_t ) = 0;
};

template<typename FunctorClass>
void CUDART_CB functorCUDAStartTimer(void *data)
{
   auto functor = reinterpret_cast<FunctorClass *>( data );
   functor->timingTree_->start( functor->getName() );
}

template<typename FunctorClass>
void CUDART_CB functorCUDAStopTimer(void *data)
{
   auto functor = reinterpret_cast<FunctorClass *>( data );
   functor->timingTree_->stop( functor->getName() );
}

template<typename FunctorType>
class FunctorCUDA : public IFunctionNodeCUDA
{
public:
   FunctorCUDA( const FunctorType &functor,
                const std::string &name,
                const TimingTreePtr &timingTree )
      : functor_( functor ), name_( name ), timingTree_( timingTree ) {}

   void operator() (cudaStream_t stream) override
   {
      if ( timingTree_ )
      {
         WALBERLA_CUDA_CHECK( cudaLaunchHostFunc( stream, functorCUDAStartTimer<FunctorCUDA<FunctorType> >, this ) );
         executiontree::internal::Caller<FunctorType>::call( functor_, stream );
         WALBERLA_CUDA_CHECK( cudaLaunchHostFunc( stream, functorCUDAStopTimer<FunctorCUDA<FunctorType> >, this ) );
      }
      else
         executiontree::internal::Caller<FunctorType>::call( functor_, stream );
   }

   const std::string getName() const override { return name_ != "" ? name_ : "FunctorCUDA"; };
   void operator() () override {  (*this)( 0 );  }

private:
   friend void CUDART_CB functorCUDAStartTimer<FunctorCUDA<FunctorType> >(void *data);
   friend void CUDART_CB functorCUDAStopTimer<FunctorCUDA<FunctorType> >(void *data);

   FunctorType functor_;
   std::string name_;
   shared_ptr< WcTimingTree > timingTree_;
};


class SequenceCUDA : public IFunctionNodeCUDA
{
public:
   SequenceCUDA( std::initializer_list< IFunctionNodeCUDAPtr > initializerList, const std::string &name, cudaStream_t defaultStream,
                 bool parallel = false, int priority=0,
                 const TimingTreePtr &timingTree = nullptr)
      : name_( name ), defaultStream_( defaultStream), timingTree_( timingTree ), parallelStreams_( priority ), parallel_( parallel ), priority_(priority)
   {
      for ( auto &e : initializerList )
         children_.push_back( e );
   }


   void operator() (cudaStream_t stream) override
   {
      if ( timingTree_ ) {
         WALBERLA_CUDA_CHECK( cudaLaunchHostFunc( stream, functorCUDAStartTimer< SequenceCUDA >, this ));
      }

      if( parallel_ )
      {
         auto parallelSection = parallelStreams_.parallelSection( stream );
         for ( auto &el : children_ )
         {
            ( *el )( parallelSection.stream());
            parallelSection.next();
         }
      }
      else
         for ( auto &el : children_ )
            (*el)( stream );

      if ( timingTree_ ) {
         WALBERLA_CUDA_CHECK( cudaLaunchHostFunc( stream, functorCUDAStopTimer< SequenceCUDA >, this ));
      }
   }

   void operator() () override {  (*this)( defaultStream_ );  }
   void push_back( const IFunctionNodeCUDAPtr &fct ) { children_.push_back( fct ); }
   void push_front( const IFunctionNodeCUDAPtr &fct ) { children_.push_front( fct ); }
   const std::string getName() const override { return name_ != "" ? name_ : "ParallelSequenceCUDA"; };
   const std::deque< IFunctionNodePtr > getChildren() const override {
      std::deque< IFunctionNodePtr > result;
      for( auto & c : children_ )
         result.push_back( c );
      return result;
   };

private:
   friend void CUDART_CB functorCUDAStartTimer< SequenceCUDA >( void *data );
   friend void CUDART_CB functorCUDAStopTimer< SequenceCUDA >( void *data );

   std::string name_;
   cudaStream_t defaultStream_;
   std::deque< IFunctionNodeCUDAPtr > children_;
   shared_ptr< WcTimingTree > timingTree_;
   cuda::ParallelStreams parallelStreams_;
   bool parallel_;
   int priority_;
};


template<typename FunctorType>
IFunctionNodeCUDAPtr functorCUDA( const FunctorType & t, const std::string &name, const shared_ptr< WcTimingTree > &timingTree )
{
   return make_shared<FunctorCUDA<FunctorType> >( t, name, timingTree );
}


shared_ptr< SequenceCUDA > sequenceCUDA( std::initializer_list< IFunctionNodeCUDAPtr > initializerList,
                                         const std::string &name, cudaStream_t defaultStream, bool parallel, int priority,
                                         const TimingTreePtr &timingTree )
{
   return make_shared< SequenceCUDA >( initializerList, name, defaultStream, parallel, priority, timingTree );
}


} // namespace executiontree
} // namespace walberla
