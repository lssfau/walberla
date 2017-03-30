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
//! \file VCycles.h
//! \ingroup pde
//! \author Dominik Bartuschat <dominik.bartuschat@fau.de>
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/uid/SUID.h"

#include "domain_decomposition/IBlock.h"

#include "field/GhostLayerField.h"

#include "pde/sweeps/RBGSFixedStencil.h"
#include "pde/sweeps/RBGS.h"

#include <boost/function.hpp>

#include <vector>


namespace walberla {
namespace pde {

template< typename Stencil_T >
class VCycles
{
public:

   typedef GhostLayerField< real_t, 1 > PdeField_T;
   typedef std::vector< real_t  > Weight_T;
   typedef GhostLayerField< real_t, Stencil_T::Size >  StencilField_T;

   VCycles( shared_ptr< StructuredBlockForest > blocks, const BlockDataID & uFieldId, const BlockDataID & fFieldId,
            const Weight_T weights,
            const uint_t iterations, const uint_t numLvl,
            const uint_t preSmoothingIters, const uint_t postSmoothingIters,
            const uint_t coarseIters, const boost::function< real_t () > & residualNorm,
            const real_t residualNormThreshold = real_t(0), const uint_t residualCheckFrequency = uint_t(1),
            const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
            const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );

   VCycles( shared_ptr< StructuredBlockForest > blocks, const BlockDataID & uFieldId, const BlockDataID & fFieldId,
            const BlockDataID & stencilFieldId,
            const uint_t iterations, const uint_t numLvl,
            const uint_t preSmoothingIters, const uint_t postSmoothingIters,
            const uint_t coarseIters, const boost::function< real_t () > & residualNorm,
            const real_t residualNormThreshold = real_t(0), const uint_t residualCheckFrequency = uint_t(1),
            const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
            const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );

   void operator()();
   void VCycle();

   uint_t iterationsPerformed() const { return iterationsPerformed_; }
   bool   thresholdReached() const { return thresholdReached_; }
   const std::vector<real_t> & convergenceRate() { return convergenceRate_; }
   static Vector3<uint_t> getSizeForLevel( const uint_t level, const shared_ptr< StructuredBlockStorage > & blocks, IBlock * const block );

protected:

   void coarsenStencilFields();

   StructuredBlockForest & blocks_;
   std::vector< Weight_T  > weights_;

   uint_t iterations_;
   uint_t numLvl_;
   uint_t preSmoothingIters_;
   uint_t postSmoothingIters_;
   uint_t coarseIters_;
   real_t residualNormThreshold_;
   uint_t residualCheckFrequency_;

   uint_t iterationsPerformed_;
   bool thresholdReached_;

   boost::function< real_t() > residualNorm_;
   std::vector<real_t> convergenceRate_;

   std::vector<BlockDataID> uId_;
   std::vector<BlockDataID> fId_;
   std::vector<BlockDataID> rId_;
   std::vector<BlockDataID> stencilId_;
   BlockDataID dId_, zId_;

   std::vector< shared_ptr<pde::RBGSFixedStencil< Stencil_T > > > RBGSFixedSweeps_;
   std::vector< shared_ptr<pde::RBGS< Stencil_T > > >             RBGSSweeps_;
   std::vector< boost::function< void() > > RBGSIteration_;
   boost::function< void() > CGIteration_;
   std::vector<boost::function< void(IBlock *) > > computeResidual_, restrict_, zeroize_, prolongateAndCorrect_;

   std::vector< blockforest::communication::UniformBufferedScheme< Stencil_T > > communication_;

   Set< SUID > requiredSelectors_;
   Set< SUID > incompatibleSelectors_;
};



} // namespace pde
} // namespace walberla


#include "VCycles.impl.h"
