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
//! \file Permeability.h
//! \ingroup lbm
//! \author Tobias Schruff <schruff@iww.rwth-aachen.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"
#include "core/math/AABB.h"
#include "core/logging/Logging.h"
#include "core/config/Config.h"
#include "core/mpi/Reduce.h"
#include "core/mpi/Broadcast.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include <string>

namespace walberla {
namespace lbm {
namespace evaluations {


template< typename PdfField_T, typename BoundaryHandling_T >
class Permeability
{
public:

   using PdfField = PdfField_T;
   using BoundaryHandling = BoundaryHandling_T;

   using LatticeModel = typename PdfField_T::LatticeModel;
   using FlagField = typename BoundaryHandling_T::FlagField;
   using flag_t = typename FlagField::flag_t;

   Permeability( real_t viscosity, const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId, const FlagUID & fluid,
                 const shared_ptr<blockforest::StructuredBlockForest> blocks );

   /*!
    *  \brief Initializes the permeability evaluation with the parameters given in the config file.
    *
    *  Configuration block:
    *
    *  \code{.unparsed}
    *  <some_name>
    *  {
    *     sampleVolume  [<0, 0, 0> , <100, 100, 100>]; [optional, default = domain]
    *     flowAxis      2;                             [required, 0-2]
    *     calcFrequency 100;                           [required, unsigned int]
    *     convCrit      1E-20;                         [optional, default = 1E-20]
    *  }
    *  \endcode
    *
    *
    *  \param config The config block handle containing the parameters.
    */
   void init( const Config::BlockHandle & config );
   
   /*!
    *  \brief Initializes the permeability evaluation with the given set of parameters.
    *
    *  \param sampleVolume  The volume to be considered during permeability evaluation.
    *  \param flowAxis      The flow axis (0 = x-axis, 1 = y-axis, 2 = z-axis) along which the permeability is evaluated. 
    *  \param calcFrequency The frequency (in time steps) in which the permeability is evaluated.
    *  \param convCrit      The delta value (du/dt) at which the flow is considered to have reached steady state.
    */
   void init( const AABB & sampleVolume, uint_t flowAxis, uint_t calcFrequency, real_t convCrit = real_t(1.0E-20) );
   
   /*!
    *  \brief Initializes the permeability evaluation with the given set of parameters.
    *
    *  \param sampleVolume  The global cell interval to be considered during permeability evaluation.
    *  \param flowAxis      The flow axis (0 = x-axis, 1 = y-axis, 2 = z-axis) along which the permeability is evaluated. 
    *  \param calcFrequency The frequency (in time steps) in which the permeability is evaluated.
    *  \param convCrit      The delta value (du/dt) at which the flow is considered to have reached steady state.
    */
   void init( const CellInterval & sampleVolume, uint_t flowAxis, uint_t calcFrequency, real_t convCrit = real_t(1.0E-20) );

   real_t convCriterion() const { return convCrit_;              }
   real_t currentDelta()  const { return delta_;                 }
   real_t currentValue()  const { return k_;                     }
   
   bool   hasConverged()  const { return ( delta_ < convCrit_ ); }

	/*!
	 *  \brief Main evaluation procedure.
	 */
   void operator()();

private:

   void initSampleVolume();

   const real_t nu_;
   const BlockDataID & pdfFieldId_;
   const BlockDataID & boundaryHandlingId_;
   const FlagUID     & fluid_;
   const shared_ptr<StructuredBlockForest> blocks_;
   
   uint_t time_;
   uint_t numSampleFluidNodes_;
   uint_t numC0FluidNodes_;
   uint_t numC1FluidNodes_;
   uint_t flowAxis_;
   uint_t interval_; 

   CellInterval sampleVolume_;
   CellInterval c0_;
   CellInterval c1_;
   
   real_t delta_;
   real_t convCrit_;
   real_t lastU_;
   real_t k_;
   
   bool initialized_;
};


} // namespace evaluations
} // namespace lbm
} // namespace walberla

#include "Permeability.impl.h"
