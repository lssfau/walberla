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
//! \file VCycles.h.impl
//! \ingroup pde
//! \author Dominik Bartuschat <dominik.bartuschat@fau.de>
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#include "VCycles.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "pde/sweeps/Multigrid.h"

#include "core/math/Limits.h"
#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"
#include "pde/ResidualNorm.h"
#include "pde/ResidualNormStencilField.h"
#include "pde/iterations/RBGSIteration.h"
#include "pde/iterations/CGFixedStencilIteration.h"
#include "pde/iterations/CGIteration.h"

#include <functional>

namespace walberla {
namespace pde {

template< typename Stencil_T, typename OperatorCoarsening_T, typename Restrict_T, typename ProlongateAndCorrect_T >
VCycles< Stencil_T, OperatorCoarsening_T, Restrict_T, ProlongateAndCorrect_T >::VCycles(
		shared_ptr< StructuredBlockForest > blocks, const BlockDataID & uFieldId, const BlockDataID & fFieldId,
        const Weight_T weights, const uint_t iterations, const uint_t numLvl,
        const uint_t preSmoothingIters, const uint_t postSmoothingIters,
        const uint_t coarseIters, const std::function< real_t () > & residualNorm,
        const real_t residualNormThreshold, const uint_t residualCheckFrequency,
        const Set<SUID> & requiredSelectors,
        const Set<SUID> & incompatibleSelectors ) :
   blocks_( *blocks ), weights_(numLvl,weights), iterations_(iterations), numLvl_(numLvl), preSmoothingIters_(preSmoothingIters),
   postSmoothingIters_(postSmoothingIters), coarseIters_(coarseIters),
   residualNormThreshold_( residualNormThreshold ), residualCheckFrequency_( residualCheckFrequency ),
   iterationsPerformed_( uint_t(0) ), thresholdReached_( false ), residualNorm_( residualNorm ), convergenceRate_(), stencilId_(),
   requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
{

   static_assert(std::is_same<OperatorCoarsening_T, CoarsenStencilFieldsDCA<Stencil_T>>::value, "Use of weight requires DCA, use constructor with stencil field if you want to employ GCA");

   // Set up fields for finest level
   uId_.push_back( uFieldId );
   fId_.push_back( fFieldId );
   rId_.push_back( field::addToStorage< PdeField_T >( blocks, "r_0", real_t(0), field::zyxf, uint_t(1) ) );

   // Check that coarsest grid has more than one cell per dimension
   auto   block = blocks->begin();
   uint_t xLvl  = blocks->getNumberOfXCells( *block );
   uint_t yLvl  = blocks->getNumberOfYCells( *block );
   uint_t zLvl  = blocks->getNumberOfZCells( *block );

   for( uint_t i = 1; i<numLvl; ++i ){

      if( xLvl % 2 != 0 || yLvl % 2 != 0 || (Stencil_T::D == uint_t(3) && zLvl % 2 != 0) )
         WALBERLA_ABORT("Number of multigrid levels is too high, not possible to further refine for CCMG!")

      xLvl = (xLvl+1)/2;
      yLvl = (yLvl+1)/2;
      if( Stencil_T::D == uint_t(3) )
         zLvl = (zLvl+1)/2;

      WALBERLA_LOG_DEVEL_ON_ROOT("Setting up coarser grid (level " << i << ") with size "<< xLvl << "x" << yLvl << "x" << zLvl );

      if( xLvl < 2 || yLvl < 2|| (Stencil_T::D == uint_t(3) && zLvl< 2) )
         WALBERLA_ABORT("Points per dimension on multigrid level " << i << " is lower than 2");

   }

   real_t scalingFactor = real_t(0.25); // scaling by ( 1/h^2 )^lvl

   // Set up fields for coarser levels
   for ( uint_t lvl = 1; lvl < numLvl; ++lvl )
   {
      auto getSize = std::bind(VCycles<Stencil_T, OperatorCoarsening_T, Restrict_T, ProlongateAndCorrect_T>::getSizeForLevel, lvl, std::placeholders::_1, std::placeholders::_2);
      uId_.push_back( field::addToStorage< PdeField_T >( blocks, "u_"+std::to_string(lvl), getSize, real_t(0), field::zyxf, uint_t(1) ) );
      fId_.push_back( field::addToStorage< PdeField_T >( blocks, "f_"+std::to_string(lvl), getSize, real_t(0), field::zyxf, uint_t(1) ) );
      rId_.push_back( field::addToStorage< PdeField_T >( blocks, "r_"+std::to_string(lvl), getSize, real_t(0), field::zyxf, uint_t(1) ) );

      for ( auto & w: weights_[lvl] )
      {
         w *= scalingFactor;
      }
      scalingFactor *= real_t(0.25);
   }

   // Set up fields for CG on coarsest level
   auto getFineSize = std::bind(VCycles<Stencil_T, OperatorCoarsening_T, Restrict_T, ProlongateAndCorrect_T>::getSizeForLevel, numLvl-1, std::placeholders::_1, std::placeholders::_2);
   dId_ = field::addToStorage< PdeField_T >( blocks, "d", getFineSize, real_t(0), field::zyxf, uint_t(1) );
   zId_ = field::addToStorage< PdeField_T >( blocks, "z", getFineSize, real_t(0), field::zyxf, uint_t(1) );

   // Set up communication
   for ( uint_t lvl = 0; lvl < numLvl-1; ++lvl )
   {
      communication_.emplace_back( blocks );
      communication_[lvl].addPackInfo( make_shared< field::communication::PackInfo< PdeField_T > >( uId_[lvl] ) );
   }

   // Set up communication for CG on coarsest level
   communication_.emplace_back( blocks );
   communication_.back().addPackInfo( make_shared< field::communication::PackInfo< PdeField_T > >( dId_ ) );

   // calculate residual norm on the finest level
   residualNorm_ = ResidualNorm<Stencil_T>( blocks->getBlockStorage(), uId_[0], fId_[0], weights_[0],
                                            requiredSelectors, incompatibleSelectors );

   // Set up RBGS iterations
   for ( uint_t lvl = 0; lvl < numLvl-1; ++lvl )
   {
      RBGSFixedSweeps_.push_back( walberla::make_shared<RBGSFixedStencil< Stencil_T > >( blocks, uId_[lvl], fId_[lvl], weights_[lvl] ) );
      RBGSIteration_.push_back( RBGSIteration(blocks->getBlockStorage(), preSmoothingIters_, communication_[lvl],
                                RBGSFixedSweeps_.back()->getRedSweep(), RBGSFixedSweeps_.back()->getBlackSweep(), [](){ return real_t(1.0); }, 0, 1,
                                requiredSelectors, incompatibleSelectors ) );
   }

   // Set up restriction and prolongation
   for (uint_t lvl = 0; lvl < numLvl-1; ++lvl)
   {
      computeResidual_.push_back( ComputeResidualFixedStencil< Stencil_T >( blocks, uId_[lvl], fId_[lvl], weights_[lvl], rId_[lvl] ) );
      restrict_.push_back( Restrict_T(blocks, rId_[lvl], fId_[lvl+1] ) );
      zeroize_.push_back( Zeroize(blocks, uId_[lvl+1]) );
      prolongateAndCorrect_.push_back( ProlongateAndCorrect_T(blocks, uId_[lvl+1], uId_[lvl]) );
   }

   // Set up CG coarse-grid iteration
   CGIteration_ = CGFixedStencilIteration< Stencil_T >( blocks->getBlockStorage(), uId_.back(), rId_.back(), dId_, zId_, fId_.back(), weights_.back(),
                                                        coarseIters_, communication_.back(), real_t(10)*math::Limits<real_t>::epsilon(),
                                                        requiredSelectors, incompatibleSelectors );
}



template< typename Stencil_T, typename OperatorCoarsening_T, typename Restrict_T, typename ProlongateAndCorrect_T >
VCycles< Stencil_T, OperatorCoarsening_T, Restrict_T, ProlongateAndCorrect_T >::VCycles(
		shared_ptr< StructuredBlockForest > blocks, const BlockDataID & uFieldId, const BlockDataID & fFieldId,
        const BlockDataID & stencilFieldId, const OperatorCoarsening_T & operatorCoarsening,
        const uint_t iterations, const uint_t numLvl,
        const uint_t preSmoothingIters, const uint_t postSmoothingIters,
        const uint_t coarseIters, const std::function< real_t () > & residualNorm,
        const real_t residualNormThreshold, const uint_t residualCheckFrequency,
        const Set<SUID> & requiredSelectors,
        const Set<SUID> & incompatibleSelectors ) :
   blocks_( *blocks ), weights_(), iterations_(iterations), numLvl_(numLvl), preSmoothingIters_(preSmoothingIters),
   postSmoothingIters_(postSmoothingIters), coarseIters_(coarseIters),
   residualNormThreshold_( residualNormThreshold ), residualCheckFrequency_( residualCheckFrequency ),
   iterationsPerformed_( uint_t(0) ), thresholdReached_( false ), residualNorm_( residualNorm ), convergenceRate_(),
   requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
{
   // Set up fields for finest level
   uId_.push_back( uFieldId );
   fId_.push_back( fFieldId );
   rId_.push_back( field::addToStorage< PdeField_T >( blocks, "r_0", real_t(0), field::zyxf, uint_t(1) ) );
   stencilId_.push_back( stencilFieldId );

   // Check that coarsest grid has more than one cell per dimension
   auto   block = blocks->begin();
   uint_t xLvl  = blocks->getNumberOfXCells( *block );
   uint_t yLvl  = blocks->getNumberOfYCells( *block );
   uint_t zLvl  = blocks->getNumberOfZCells( *block );

   for( uint_t i = 1; i<numLvl; ++i ){

      if( xLvl % 2 != 0 || yLvl % 2 != 0 || (Stencil_T::D == uint_t(3) && zLvl % 2 != 0) )
         WALBERLA_ABORT("Number of multigrid levels is too high, not possible to further refine for CCMG!")

      xLvl = (xLvl+1)/2;
      yLvl = (yLvl+1)/2;
      if( Stencil_T::D == uint_t(3) )
         zLvl = (zLvl+1)/2;

      WALBERLA_LOG_DEVEL_ON_ROOT("Setting up coarser grid (level " << i << ") with size "<< xLvl << "x" << yLvl << "x" << zLvl );

      if( xLvl < 2 || yLvl < 2|| (Stencil_T::D == uint_t(3) && zLvl< 2) )
         WALBERLA_ABORT("Points per dimension on multigrid level " << i << " is lower than 2");

   }

   // Set up fields for coarser levels
   for ( uint_t lvl = 1; lvl < numLvl; ++lvl )
   {
      auto getSize = std::bind(VCycles<Stencil_T, OperatorCoarsening_T, Restrict_T, ProlongateAndCorrect_T>::getSizeForLevel, lvl, std::placeholders::_1, std::placeholders::_2);
      uId_.push_back( field::addToStorage< PdeField_T >( blocks, "u_"+std::to_string(lvl), getSize, real_t(0), field::zyxf, uint_t(1) ) );
      fId_.push_back( field::addToStorage< PdeField_T >( blocks, "f_"+std::to_string(lvl), getSize, real_t(0), field::zyxf, uint_t(1) ) );
      rId_.push_back( field::addToStorage< PdeField_T >( blocks, "r_"+std::to_string(lvl), getSize, real_t(0), field::zyxf, uint_t(1) ) );
      stencilId_.push_back( field::addToStorage< StencilField_T >( blocks, "w_"+std::to_string(lvl), getSize, real_t(0), field::zyxf, uint_t(1) ) );
   }

   // CoarsenStencilFieldsDCA<Stencil_T>( blocks, stencilId_, numLvl, uint_t(2)) ();  // scaling by ( 1/h^2 )^lvl
   // CoarsenStencilFieldsGCA<Stencil_T>( blocks, stencilId_, numLvl, real_t(2)) ();
   operatorCoarsening(stencilId_);

   // Set up fields for CG on coarsest level
   auto getFineSize = std::bind(VCycles<Stencil_T, OperatorCoarsening_T, Restrict_T, ProlongateAndCorrect_T>::getSizeForLevel, numLvl-1, std::placeholders::_1, std::placeholders::_2);
   dId_ = field::addToStorage< PdeField_T >( blocks, "d", getFineSize, real_t(0), field::zyxf, uint_t(1) );
   zId_ = field::addToStorage< PdeField_T >( blocks, "z", getFineSize, real_t(0), field::zyxf, uint_t(1) );

   // Set up communication
   for ( uint_t lvl = 0; lvl < numLvl-1; ++lvl )
   {
      communication_.emplace_back( blocks );
      communication_[lvl].addPackInfo( make_shared< field::communication::PackInfo< PdeField_T > >( uId_[lvl] ) );
   }

   // Set up communication for CG on coarsest level
   communication_.emplace_back( blocks );
   communication_.back().addPackInfo( make_shared< field::communication::PackInfo< PdeField_T > >( dId_ ) );

   // calculate residual norm on the finest level
   residualNorm_ = ResidualNormStencilField<Stencil_T>( blocks->getBlockStorage(), uId_[0], fId_[0], stencilId_[0],
                                                        requiredSelectors, incompatibleSelectors );

   // Set up RBGS iterations
   for ( uint_t lvl = 0; lvl < numLvl-1; ++lvl )
   {
      RBGSSweeps_.push_back( walberla::make_shared<RBGS< Stencil_T > >( blocks, uId_[lvl], fId_[lvl], stencilId_[lvl] ) );
      RBGSIteration_.push_back( RBGSIteration(blocks->getBlockStorage(), preSmoothingIters_, communication_[lvl],
                                RBGSSweeps_.back()->getRedSweep(), RBGSSweeps_.back()->getBlackSweep(), [](){ return real_t(1.0); }, 0, 1,
                                requiredSelectors, incompatibleSelectors ) );
   }

   // Set up restriction and prolongation
   for (uint_t lvl = 0; lvl < numLvl-1; ++lvl)
   {
      computeResidual_.push_back( ComputeResidual< Stencil_T >( blocks, uId_[lvl], fId_[lvl], stencilId_[lvl], rId_[lvl] ) );
      restrict_.push_back( Restrict_T(blocks, rId_[lvl], fId_[lvl+1] ) );
      zeroize_.push_back( Zeroize(blocks, uId_[lvl+1]) );
      prolongateAndCorrect_.push_back( ProlongateAndCorrect_T(blocks, uId_[lvl+1], uId_[lvl]) );
   }

   // Set up CG coarse-grid iteration
   CGIteration_ = CGIteration< Stencil_T >( blocks->getBlockStorage(), uId_.back(), rId_.back(), dId_, zId_, fId_.back(), stencilId_.back(),
                                            coarseIters_, communication_.back(), real_t(10)*math::Limits<real_t>::epsilon(),
                                            requiredSelectors, incompatibleSelectors );
}



template< typename Stencil_T, typename OperatorCoarsening_T, typename Restrict_T, typename ProlongateAndCorrect_T >
void VCycles< Stencil_T, OperatorCoarsening_T, Restrict_T, ProlongateAndCorrect_T >::operator()()
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Starting VCycle iteration with a maximum number of " << iterations_ << " cycles and " << numLvl_ << " levels" );
   thresholdReached_ = false;

   real_t residualNorm_old = real_t(1);

   for( uint_t i = 0; i < iterations_; ++i )
   {
      // compute L2 norm of residual as termination criterion
      if( residualNormThreshold_ > real_t(0) && residualCheckFrequency_ > uint_t(0) )
      {
         if( (i % residualCheckFrequency_) == uint_t(0) )
         {
            communication_[0]();
            const real_t residualNorm = residualNorm_();

            WALBERLA_CHECK( math::finite(residualNorm), "Non-finite residual norm detected during the Multigrid V-cycle, "
                                                        "the simulation has probably diverged." );

            // calculating convergence rate
            convergenceRate_.push_back( residualNorm / residualNorm_old );
            residualNorm_old = residualNorm;
            if( i > 0 )
            {
               WALBERLA_LOG_PROGRESS_ON_ROOT("** Multigrid convergence rate: "<< convergenceRate_.back());
            }

            WALBERLA_LOG_PROGRESS_ON_ROOT( "Residual norm after " << i << " Multigrid V-cycles: " << residualNorm );
            if( residualNorm < residualNormThreshold_ )
            {
               WALBERLA_LOG_PROGRESS_ON_ROOT( "Aborting Multigrid V (residual norm threshold reached):"
                                              "\n  residual norm threshold: " << residualNormThreshold_ <<
                                              "\n  residual norm:           " << residualNorm );
               thresholdReached_ = true;
               iterationsPerformed_ = i;
               break;
            }
         }
      }

      // perform one V-cycle
      VCycle();
      
      iterationsPerformed_ = i+1;
   }

   WALBERLA_LOG_PROGRESS_ON_ROOT( "Multigrid V-cycle " << (thresholdReached_ ? "finished" : "aborted") << " after " << iterationsPerformed_ << " iterations" );
}



template< typename Stencil_T, typename OperatorCoarsening_T, typename Restrict_T, typename ProlongateAndCorrect_T >
void VCycles< Stencil_T, OperatorCoarsening_T, Restrict_T, ProlongateAndCorrect_T >::VCycle()
{
   // pre-smoothen -- go from fine to coarse
   for (uint_t l = 0; l < numLvl_-1; ++l)
   {
      WALBERLA_LOG_PROGRESS_ON_ROOT("Communicating, smoothing and communicating on level: "<< l);

      // smoothen
      RBGSIteration_[l]();
      communication_[l]();

      WALBERLA_LOG_PROGRESS_ON_ROOT( "Restricting residual from level "<< l << " to rhs on level " << l+1 );
      WALBERLA_LOG_PROGRESS_ON_ROOT( "Zeroing solution on level " << l+1 );

      // compute and restrict residual
      for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
      {
         computeResidual_[l](&*block);
         restrict_[l](&*block);
         zeroize_[l](&*block);
      }
   }

   WALBERLA_LOG_PROGRESS_ON_ROOT("Solving coarsest grid, level "<< numLvl_ - 1 );

   // solve coarsest level
   CGIteration_();

   // correct and post-smoothen -- go from coarse to fine
   for (uint_t ll = 0; ll < numLvl_-1; ++ll)
   {
      uint_t l = numLvl_-2 - ll;

      WALBERLA_LOG_PROGRESS_ON_ROOT( "Prolongating solution from level "<< l+1 << " to " << l );

      // prolongate and correct solution
      for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
      {
         prolongateAndCorrect_[l](&*block);
      }

      WALBERLA_LOG_PROGRESS_ON_ROOT("Communicating, smoothing and communicating on level: "<< l);

      // smoothen
      RBGSIteration_[l]();
      communication_[l]();
   }
}



template< typename Stencil_T, typename OperatorCoarsening_T, typename Restrict_T, typename ProlongateAndCorrect_T >
Vector3<uint_t> VCycles< Stencil_T, OperatorCoarsening_T, Restrict_T, ProlongateAndCorrect_T >::getSizeForLevel( const uint_t level, const shared_ptr< StructuredBlockStorage > & blocks, IBlock * const block )
{
   Vector3<uint_t> cells( blocks->getNumberOfXCells( *block ), blocks->getNumberOfYCells( *block ), blocks->getNumberOfZCells( *block ) );

   if( level == 0 )
      return cells;

   WALBERLA_ASSERT_EQUAL(cells[0] % (uint_t(2) << uint_c(level-1)), 0, "can only coarsen an even number of cells!");
   WALBERLA_ASSERT_EQUAL(cells[1] % (uint_t(2) << uint_c(level-1)), 0, "can only coarsen an even number of cells!");

   cells[0] = (cells[0] >> level);
   cells[1] = (cells[1] >> level);
   if( Stencil_T::D == 3 )
   {
      WALBERLA_ASSERT_EQUAL(cells[2] % (uint_t(2) << uint_c(level-1)), 0, "can only coarsen an even number of cells!");
      cells[2] = (cells[2] >> level);
   }

   return cells;
}



} // namespace pde
} // namespace walberla
