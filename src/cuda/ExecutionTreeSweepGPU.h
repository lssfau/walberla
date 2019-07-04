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
//! \file ExecutionTreeSweepGPU.h
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//==============================================================================================================================================================

#pragma once

#include "domain_decomposition/IBlock.h"
#include "executiontree/ExecutionTree.h"
#include "ExecutionTreeGPU.h"

namespace walberla {
namespace executiontree {


template<typename FunctorType>
IFunctionNodeCUDAPtr sweepCUDA( BlockStorage &bs, const FunctorType & t, const std::string &name = "", const TimingTreePtr &timingTree = nullptr );

template<typename FunctorType>
IFunctionNodeCUDAPtr sweepCUDA( const shared_ptr< StructuredBlockStorage > &bs, const FunctorType & t, const std::string &name = "",
                                const TimingTreePtr &tt = nullptr );


template<typename FunctorType>
class SweepCUDA : public IFunctionNodeCUDA
{
public:
   SweepCUDA( BlockStorage &bs,
              const FunctorType &functor,
              const std::string &name,
              const TimingTreePtr &timingTree )
      : blockStorage_( bs ),
        functor_( functor ),
        name_( name ),
        timingTree_( timingTree ) {}

   SweepCUDA( const shared_ptr <StructuredBlockStorage> &bs,
              const FunctorType &functor,
              const std::string &name,
              const TimingTreePtr &timingTree )
      : blockStorage_( bs->getBlockStorage()),
        functor_( functor ),
        name_( name ),
        timingTree_( timingTree ) {}

   void operator() () override {  (*this)( 0 );  }

   void operator()( cudaStream_t stream ) override
   {
      if ( timingTree_ )
      {
         for ( auto &block: blockStorage_ )
         {
            timingTree_->start( name_ );
            executiontree::internal::Caller<FunctorType>::call( functor_, &block, stream );
            timingTree_->stop( name_ );
         }
      }
      else
         for ( auto &block: blockStorage_ )
            executiontree::internal::Caller<FunctorType>::call( functor_, &block, stream );
   }

   const std::string getName() const override { return name_ != "" ? name_ : "Sweep"; };

private:
   BlockStorage &blockStorage_;

   FunctorType functor_;
   std::string name_;
   TimingTreePtr timingTree_;
};

template<typename FunctorType>
IFunctionNodeCUDAPtr sweepCUDA( BlockStorage &bs, FunctorType t, const std::string &name, const shared_ptr< WcTimingTree > &timingTree )
{
   return make_shared<SweepCUDA<FunctorType> >( bs, t, name, timingTree );
}

template<typename FunctorType>
IFunctionNodeCUDAPtr sweepCUDA( const shared_ptr< StructuredBlockStorage > &bs, const FunctorType & t, const std::string &name,
                                const TimingTreePtr &timingTree )
{
   return make_shared<SweepCUDA<FunctorType> >( bs, t, name, timingTree );
}


} // namespace executiontree
} // namespace walberla
