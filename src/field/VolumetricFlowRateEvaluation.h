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
//! \file VolumetricFlowRateEvaluation.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/uid/SUID.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/EvaluationFilter.h"
#include "field/iterators/IteratorMacros.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <string>
#include <type_traits>



namespace walberla {
namespace field {

namespace internal {

using FlowRateSolution_T = std::function<real_t ()>;
using FlowRateVelocitySolution_T = std::function<Vector3<real_t> (const Vector3<real_t> &)>;

const std::string     volumetricFlowRateEvaluationFilename("flowrate.dat");
const real_t          volumetricFlowRateEvaluationNormalization( real_t(1) );
const Vector3<bool>   volumetricFlowRateEvaluationAxis( Vector3<bool>( true, false, false ) );
const Vector3<real_t> volumetricFlowRateEvaluationPoint( Vector3<real_t>( real_c(0.5) ) );

const std::string volumetricFlowRateEvaluationConfigBlock("VolumetricFlowRateEvaluation");

}



//**********************************************************************************************************************
/*!
*   \brief Class for evaluating the volumetric flow rate of a simulation
*
*   \section docVolumetricFlowRateEvaluation Volumetric Flow Rate Evaluation
*
*   Class for evaluating the volumetric flow rate of a simulation. The flow rate can either just be evaluated or it can
*   be compared with an exact solution in case an exact solution exists and the exact solution is provided via a
*   callback function.
*
*   Do not create objects of class VolumetricFlowRateEvaluation directly, better use one of the various
*   'makeVolumetricFlowRateEvaluation' functions below!
*
*   Template parameters:
*   - VelocityField_T: the field storing the simulation velocities (must be a field that stores Vector3 values)
*   - Filter_T: the type of the evaluation filter (see \ref docEvaluationFilter in 'EvaluationFilter.h')
*
*   Parameters for setting up and controlling flow rate evaluation:
*   - blocks: the block storage
*   - fieldId: block data ID of the velocity field
*   - filter: the evaluation filter that indicates which cells are processed
*   - plot frequency: the plotting interval - used for saving the data to file.
*                     If set to '0', no plotting data is created.
*   - log frequency: the logging interval - used for logging the data via the Logging singleton.
*                    If set to '0', no logging is performed.
*   - solution: the solution callback function - must return the exact solution for the flow rate
*   - velocity solution: the velocity solution callback function - must return the exact solution for the velocity when
*                        called with a position inside the domain
*   - filename: the name of the file that stores the data for plotting
*   - normalization factor: an optional factor the simulation values and the solution values are multiplied with
*   - domain normalization: By default, the surface area required for the calculation of the volumetric flow rate
*                           is calculated using the size of the domain as returned by blocks->getDomain(). However,
*                           you can overwrite the size of the domain via the domain normalization parameter: This
*                           parameter specifies a new, different size for the domain.
*   - axis: boolean Vector3 where only one component must be true
*           -> the flow rate through a surface _perpendicular_ to this axis will be calculated
*   - surface point: a point that specifies the position of the surface - this point is always relative to the entire
*                    domain with [<0,0,0>,<1,1,1>], i.e., all coordinates of this point (x, y, and z) must be in [0;1]
*   - required and incompatible selectors
*
*   You do not have to specify an evaluation filter! If you do not specify any filter, _all_ cells are processed and no
*   cell is excluded.
*
*   If you want to use a flag field as evaluation filter, fitting 'makeVolumetricFlowRateEvaluation' functions already
*   exist. These functions need an additional template parameter FlagField_T and you have to provide the block data ID
*   of the flag field together with a set of flag UIDs that specify which cells need to be processed.
*
*   There also exist 'makeVolumetricFlowRateEvaluation' functions that take configuration file data as an additional
*   parameter in order to parse the configuration file for setting up and controlling flow rate evaluation.
*   The configuration file block looks like as follows:
*
*   \code
*   VolumetricFlowRateEvaluation
*   {
*      plotFrequency       [unsigned integer]; // the plot frequency
*      logFrequency        [unsigned integer]; // the log frequency
*      filename            [string]; // the name of the file that stores the data for plotting
*      normalization       [floating point value]; // normalization factor
*      domainNormalization [Vector3: <real,real,real>]; // domain normalization
*      axis                [Vector3: <bool,bool,bool>]; // axis perpendicular to the flow rate evaluation surface
*      point               [Vector3: <real,real,real>]; // surface point
*   }
*   \endcode
*
*   Example:
*
*   \code
*   VolumetricFlowRateEvaluation
*   {
*      plotFrequency 10;
*      logFrequency  1000;
*      filename      FlowRate.txt;
*      axis          <true,false,false>;
*      point         <0.5,0.5,0.5>;
*   }
*   \endcode
*
*   Note that the shared pointer returned by all 'makeVolumetricFlowRateEvaluation' functions can be captured by a
*   SharedFunctor for immediate registration at a time loop (see field::makeSharedFunctor).
*/
//**********************************************************************************************************************

template< typename VelocityField_T,
          typename Filter_T = DefaultEvaluationFilter >
class VolumetricFlowRateEvaluation
{
public:

   using Solution_T = internal::FlowRateSolution_T;
   using VelocitySolution_T = internal::FlowRateVelocitySolution_T;

   VolumetricFlowRateEvaluation( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                                 const Filter_T & filter,
                                 const uint_t plotFrequency, const uint_t logFrequency,
                                 const Solution_T & _solution = Solution_T(),
                                 const VelocitySolution_T & velocitySolution = VelocitySolution_T(),
                                 const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                 const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), fieldId_( fieldId ), filter_( filter ), solution_( _solution ), velocitySolution_( velocitySolution ),
      executionCounter_( uint_t(0) ), plotFrequency_( plotFrequency ), logFrequency_( logFrequency ),
      filename_( internal::volumetricFlowRateEvaluationFilename ),
      normalizationFactor_( internal::volumetricFlowRateEvaluationNormalization ),
      axis_( internal::volumetricFlowRateEvaluationAxis ), surfacePoint_( internal::volumetricFlowRateEvaluationPoint ),
      flowRate_( real_t(0) ), velocitySolutionFlowRate_( real_t(0) ),
      requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {
      auto _blocks = blocks.lock();
      WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to access 'VolumetricFlowRateEvaluation' for a block storage object that doesn't exist anymore" );
      domainNormalization_ = _blocks->getDomain().sizes();
   }

   VolumetricFlowRateEvaluation( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                                 const uint_t plotFrequency, const uint_t logFrequency,
                                 const Solution_T & _solution = Solution_T(),
                                 const VelocitySolution_T & velocitySolution = VelocitySolution_T(),
                                 const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                 const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), fieldId_( fieldId ), filter_( Filter_T() ), solution_( _solution ), velocitySolution_( velocitySolution ),
      executionCounter_( uint_t(0) ), plotFrequency_( plotFrequency ), logFrequency_( logFrequency ),
      filename_( internal::volumetricFlowRateEvaluationFilename ),
      normalizationFactor_( internal::volumetricFlowRateEvaluationNormalization ),
      axis_( internal::volumetricFlowRateEvaluationAxis ), surfacePoint_( internal::volumetricFlowRateEvaluationPoint ),
      flowRate_( real_t(0) ), velocitySolutionFlowRate_( real_t(0) ),
      requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {
      static_assert( (std::is_same< Filter_T, DefaultEvaluationFilter >::value),
                     "This constructor is only available if DefaultEvaluationFilter is set as filter type!" );

      auto _blocks = blocks.lock();
      WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to access 'VolumetricFlowRateEvaluation' for a block storage object that doesn't exist anymore" );
      domainNormalization_ = _blocks->getDomain().sizes();
   }

   void setFilename( const std::string & filename ) { filename_ = filename; }
   void setNormalizationFactor( const real_t f ) { normalizationFactor_ = f; }
   void setDomainNormalization( const Vector3<real_t> & d ) { domainNormalization_ = d; }
   void setAxis( const Vector3< bool > & axis ) { WALBERLA_ASSERT( axis_[0] || axis_[1] || axis_[2] ); axis_ = axis; }
   void setSurfacePoint( const Vector3<real_t> & p ) { surfacePoint_ = p; }

   real_t flowRate() const { return flowRate_; }
   real_t velocitySolutionFlowRate() const { return velocitySolutionFlowRate_; }

   real_t solution() const
   {
      if( !solution_ )
         return real_t(0);

      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'VolumetricFlowRateEvaluation' for a block storage object that doesn't exist anymore" );

      real_t _solution = normalizationFactor_ * solution_();
      const auto & domainAABB = blocks->getDomain();
      if( axis_[0] )
      {
         _solution *= ( domainNormalization_[1] / domainAABB.ySize() ) * ( domainNormalization_[2] / domainAABB.zSize() );
      }
      else if( axis_[1] )
      {
         _solution *= ( domainNormalization_[0] / domainAABB.xSize() ) * ( domainNormalization_[2] / domainAABB.zSize() );
      }
      else
      {
         _solution *= ( domainNormalization_[0] / domainAABB.xSize() ) * ( domainNormalization_[1] / domainAABB.ySize() );
      }
      return _solution;
   }

   void operator()();

private:

   weak_ptr< StructuredBlockStorage > blocks_;
   ConstBlockDataID fieldId_;

   Filter_T filter_;

   Solution_T solution_;
   VelocitySolution_T velocitySolution_;

   uint_t executionCounter_;

   uint_t plotFrequency_;
   uint_t  logFrequency_;

   std::string filename_;

   real_t normalizationFactor_;
   Vector3<real_t> domainNormalization_;

   Vector3< bool > axis_;
   Vector3<real_t> surfacePoint_;

   real_t flowRate_;
   real_t velocitySolutionFlowRate_;

   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;

}; // class VolumetricFlowRateEvaluation



template< typename VelocityField_T, typename Filter_T >
void VolumetricFlowRateEvaluation< VelocityField_T, Filter_T >::operator()()
{
   if( logFrequency_ == uint_t(0) && ( plotFrequency_ == uint_t(0) || filename_.empty() ) )
      return;

   ++executionCounter_;

   const bool plot = ( plotFrequency_ != uint_t(0) && ( executionCounter_ - uint_c(1) ) % plotFrequency_ == uint_t(0) && !filename_.empty() );
   const bool log  = ( logFrequency_  != uint_t(0) && ( executionCounter_ - uint_c(1) ) % logFrequency_  == uint_t(0) );

   if( !log && !plot )
      return;

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'VolumetricFlowRateEvaluation' for a block storage object that doesn't exist anymore" );

   const auto & domainAABB = blocks->getDomain();

   Vector3<real_t> sp = domainAABB.minCorner();
   sp[0] += surfacePoint_[0] * domainAABB.xSize();
   sp[1] += surfacePoint_[1] * domainAABB.ySize();
   sp[2] += surfacePoint_[2] * domainAABB.zSize();

   real_t _flowRate( real_t(0) );
   real_t _velocitySolutionFlowRate( real_t(0) );

   for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
   {
      const VelocityField_T * field = block->template getData< const VelocityField_T >( fieldId_ );

      const auto level = blocks->getLevel( *block );

      const real_t normalizedDx = domainNormalization_[0] * blocks->dx(level) / domainAABB.xSize();
      const real_t normalizedDy = domainNormalization_[1] * blocks->dy(level) / domainAABB.ySize();
      const real_t normalizedDz = domainNormalization_[2] * blocks->dz(level) / domainAABB.zSize();

      filter_( *block );

      if( axis_[0] )
      {
         const auto factor = normalizationFactor_ * normalizedDy * normalizedDz;
         if( !velocitySolution_ )
         {
            WALBERLA_FOR_ALL_CELLS_XYZ_OMP( field, omp parallel for schedule(static) reduction(+:_flowRate),

               if( filter_(x,y,z) )
               {
                  const auto aabb = blocks->getBlockLocalCellAABB( *block, Cell(x,y,z) );
                  const auto center = aabb.center();
                  if( aabb.contains( sp[0], center[1], center[2] ) )
                     _flowRate += field->get(x,y,z)[0] * factor;
               }
            )
         }
         else
         {
            WALBERLA_FOR_ALL_CELLS_XYZ_OMP( field, omp parallel for schedule(static) reduction(+:_flowRate) reduction(+:_velocitySolutionFlowRate),

               if( filter_(x,y,z) )
               {
                  const auto aabb = blocks->getBlockLocalCellAABB( *block, Cell(x,y,z) );
                  const auto center = aabb.center();
                  if( aabb.contains( sp[0], center[1], center[2] ) )
                  {
                     _flowRate += field->get(x,y,z)[0] * factor;
                     _velocitySolutionFlowRate += velocitySolution_(center)[0] * factor;
                  }
               }
            )
         }
      }
      else if( axis_[1] )
      {
         const auto factor = normalizationFactor_ * normalizedDx * normalizedDz;
         if( !velocitySolution_ )
         {
            WALBERLA_FOR_ALL_CELLS_XYZ_OMP( field, omp parallel for schedule(static) reduction(+:_flowRate),

               if( filter_(x,y,z) )
               {
                  const auto aabb = blocks->getBlockLocalCellAABB( *block, Cell(x,y,z) );
                  const auto center = aabb.center();
                  if( aabb.contains( center[0], sp[1], center[2] ) )
                     _flowRate += field->get(x,y,z)[1] * factor;
               }
            )
         }
         else
         {
            WALBERLA_FOR_ALL_CELLS_XYZ_OMP( field, omp parallel for schedule(static) reduction(+:_flowRate) reduction(+:_velocitySolutionFlowRate),

               if( filter_(x,y,z) )
               {
                  const auto aabb = blocks->getBlockLocalCellAABB( *block, Cell(x,y,z) );
                  const auto center = aabb.center();
                  if( aabb.contains( center[0], sp[1], center[2] ) )
                  {
                     _flowRate += field->get(x,y,z)[1] * factor;
                     _velocitySolutionFlowRate += velocitySolution_(center)[1] * factor;
                  }
               }
            )
         }
      }
      else
      {
         const auto factor = normalizationFactor_ * normalizedDx * normalizedDy;
         if( !velocitySolution_ )
         {
            WALBERLA_FOR_ALL_CELLS_XYZ_OMP( field, omp parallel for schedule(static) reduction(+:_flowRate),

               if( filter_(x,y,z) )
               {
                  const auto aabb = blocks->getBlockLocalCellAABB( *block, Cell(x,y,z) );
                  const auto center = aabb.center();
                  if( aabb.contains( center[0], center[1], sp[2] ) )
                     _flowRate += field->get(x,y,z)[2] * factor;
               }
            )
         }
         else
         {
            WALBERLA_FOR_ALL_CELLS_XYZ_OMP( field, omp parallel for schedule(static) reduction(+:_flowRate) reduction(+:_velocitySolutionFlowRate),

               if( filter_(x,y,z) )
               {
                  const auto aabb = blocks->getBlockLocalCellAABB( *block, Cell(x,y,z) );
                  const auto center = aabb.center();
                  if( aabb.contains( center[0], center[1], sp[2] ) )
                  {
                     _flowRate += field->get(x,y,z)[2] * factor;
                     _velocitySolutionFlowRate += velocitySolution_(center)[2] * factor;
                  }
               }
            )
         }
      }
   }

   mpi::reduceInplace( _flowRate, mpi::SUM );
   mpi::reduceInplace( _velocitySolutionFlowRate, mpi::SUM );

   flowRate_ = _flowRate;
   velocitySolutionFlowRate_ = _velocitySolutionFlowRate;
   const real_t _solution = solution();

   WALBERLA_ROOT_SECTION()
   {
      const auto & id = blocks->getBlockDataIdentifier( fieldId_ );

      if( plot && executionCounter_ == uint_t(1) )
      {
         std::ofstream file( filename_.c_str() );
         file << "# volumetric flow rate evaluation of data '" << id <<  "'\n"
              << "# step [1], flow rate (simulation) [2], flow rate (\"discrete\" solution) [3], flow rate (exact solution) [4]"
                 ", rel. error (base = \"discrete\" solution) [5], rel. error (base = exact solution) [6]" << std::endl;
         if( !velocitySolution_ )
            file << "ATTENTION: \"discrete\" solution values not available and thus not valid!" << std::endl;
         if( !solution_ )
            file << "ATTENTION: exact solution values not available and thus not valid!" << std::endl;
         file.close();
      }

      if( log )
      {
         if( velocitySolution_ && solution_ )
         {
            WALBERLA_LOG_INFO( "Evaluation of the volumetric flow rate [data = '" << id << "']:" <<
                               "\n- exact solution:      " << _solution <<
                               "\n- \"discrete\" solution: " <<  velocitySolutionFlowRate_ << " (using reference velocities evaluated at cell centers)" <<
                               "\n- simulation:          " <<  flowRate_ <<  " (using computed velocities evaluated at cell centers)" <<
                               "\n- rel. error (of the simulation):"
                               "\n   + " << ( std::abs( ( flowRate_ - velocitySolutionFlowRate_ ) / velocitySolutionFlowRate_ ) ) << " (base = \"discrete\" solution)" <<
                               "\n   + " << ( std::abs( ( flowRate_ - _solution ) / _solution ) ) << " (base = exact solution)" );
         }
         else if( velocitySolution_ )
         {
            WALBERLA_ASSERT( !solution_ );
            WALBERLA_LOG_INFO( "Evaluation of the volumetric flow rate [data = '" << id << "']:" <<
                               "\n- \"discrete\" solution: " <<  velocitySolutionFlowRate_ << " (using reference velocities evaluated at cell centers)" <<
                               "\n- simulation:          " <<  flowRate_ <<  " (using computed velocities evaluated at cell centers)" <<
                               "\n- rel. error (of the simulation): " << ( std::abs( flowRate_ - velocitySolutionFlowRate_ ) / velocitySolutionFlowRate_ ) );

         }
         else if( solution_ )
         {
            WALBERLA_ASSERT( !velocitySolution_ );
            WALBERLA_LOG_INFO( "Evaluation of the volumetric flow rate [data = '" << id << "']:" <<
                               "\n- exact solution: " << _solution <<
                               "\n- simulation:     " <<  flowRate_ <<
                               "\n- rel. error (of the simulation): " << ( std::abs( flowRate_ - _solution ) / _solution ) );

         }
         else
         {
            WALBERLA_ASSERT( !solution_ && !velocitySolution_ );
            WALBERLA_LOG_INFO( "Evaluation of the volumetric flow rate [data = '" << id << "']: " << flowRate_ );
         }
      }

      if( plot )
      {
         std::ofstream file( filename_.c_str(), std::ofstream::out | std::ofstream::app );
         file << ( executionCounter_ - uint_t(1) ) << " " << flowRate_ << " " << velocitySolutionFlowRate_ << " " << _solution << " "
              << ( isIdentical( velocitySolutionFlowRate_, real_t(0) ) ? real_t(0) : std::abs( ( flowRate_ - velocitySolutionFlowRate_ ) / velocitySolutionFlowRate_ ) ) << " "
              << ( isIdentical( _solution, real_t(0) ) ? real_t(0) : std::abs( ( flowRate_ - _solution ) / _solution ) ) << std::endl;
         file.close();
      }
   }
}



///////////////////////////////////////////////////////////////////////////
// makeVolumetricFlowRateEvaluation functions without configuration file //
///////////////////////////////////////////////////////////////////////////
                                 
template< typename VelocityField_T >
shared_ptr< VolumetricFlowRateEvaluation< VelocityField_T > > makeVolumetricFlowRateEvaluation( const weak_ptr< StructuredBlockStorage > & blocks,
                                                                                                const ConstBlockDataID & velocityFieldId,
                                                                                                const uint_t plotFrequency, const uint_t logFrequency,
                                                                                                const internal::FlowRateSolution_T & solution = internal::FlowRateSolution_T(),
                                                                                                const internal::FlowRateVelocitySolution_T & velocitySolution = internal::FlowRateVelocitySolution_T(),
                                                                                                const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                                                const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using FR_T = VolumetricFlowRateEvaluation<VelocityField_T>;
   return shared_ptr< FR_T >( new FR_T( blocks, velocityFieldId, plotFrequency, logFrequency, solution, velocitySolution, requiredSelectors, incompatibleSelectors ) );
}

template< typename VelocityField_T, typename FlagField_T >
shared_ptr< VolumetricFlowRateEvaluation< VelocityField_T, FlagFieldEvaluationFilter<FlagField_T> > >
makeVolumetricFlowRateEvaluation( const weak_ptr< StructuredBlockStorage > & blocks,
                                  const ConstBlockDataID & velocityFieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                                  const uint_t plotFrequency, const uint_t logFrequency,
                                  const internal::FlowRateSolution_T & solution = internal::FlowRateSolution_T(),
                                  const internal::FlowRateVelocitySolution_T & velocitySolution = internal::FlowRateVelocitySolution_T(),
                                  const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                  const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using FR_T = VolumetricFlowRateEvaluation<VelocityField_T, FlagFieldEvaluationFilter<FlagField_T>>;
   return shared_ptr< FR_T >( new FR_T( blocks, velocityFieldId, FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                        plotFrequency, logFrequency, solution, velocitySolution, requiredSelectors, incompatibleSelectors ) );
}

template< typename VelocityField_T, typename Filter_T >
shared_ptr< VolumetricFlowRateEvaluation< VelocityField_T, Filter_T > >
makeVolumetricFlowRateEvaluation( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & velocityFieldId,
                                  const Filter_T & filter, const uint_t plotFrequency, const uint_t logFrequency,
                                  const internal::FlowRateSolution_T & solution = internal::FlowRateSolution_T(),
                                  const internal::FlowRateVelocitySolution_T & velocitySolution = internal::FlowRateVelocitySolution_T(),
                                  const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                  const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using FR_T = VolumetricFlowRateEvaluation<VelocityField_T, Filter_T>;
   return shared_ptr< FR_T >( new FR_T( blocks, velocityFieldId, filter, plotFrequency, logFrequency, solution, velocitySolution, requiredSelectors, incompatibleSelectors ) );
}



/////////////////////////////////////////////////////////////////////
// makeVolumetricFlowRateEvaluation functions + configuration file //
/////////////////////////////////////////////////////////////////////

namespace internal {

inline void volumetricFlowRateEvaluationConfigParser( const Config::BlockHandle & parentBlockHandle, const std::string & configBlockName,
                                                      uint_t & defaultPlotFrequency, uint_t & defaultLogFrequency,
                                                      std::string & defaultFilename, real_t & defaultNormalizationFactor, Vector3<real_t> & defaultDomainNormalization,
                                                      Vector3<bool> & defaultAxis, Vector3<real_t> & defaultSurfacePoint )
{
   if( parentBlockHandle )
   {
      Config::BlockHandle block = parentBlockHandle.getBlock( configBlockName );
      if( block )
      {
         defaultPlotFrequency = block.getParameter< uint_t >( "plotFrequency", defaultPlotFrequency );
         defaultLogFrequency = block.getParameter< uint_t >( "logFrequency", defaultLogFrequency );
         defaultFilename = block.getParameter< std::string >( "filename", defaultFilename );
         defaultNormalizationFactor = block.getParameter< real_t >( "normalization", defaultNormalizationFactor );
         defaultDomainNormalization = block.getParameter< Vector3<real_t> >( "domainNormalization", defaultDomainNormalization );
         defaultAxis = block.getParameter< Vector3<bool> >( "axis", defaultAxis );
         defaultSurfacePoint = block.getParameter< Vector3<real_t> >( "point", defaultSurfacePoint );
      }
   }
}

inline void volumetricFlowRateEvaluationConfigParser( const shared_ptr< Config > & config, const std::string & configBlockName,
                                                      uint_t & defaultPlotFrequency, uint_t & defaultLogFrequency,
                                                      std::string & defaultFilename, real_t & defaultNormalizationFactor, Vector3<real_t> & defaultDomainNormalization,
                                                      Vector3<bool> & defaultAxis, Vector3<real_t> & defaultSurfacePoint )
{
   if( !!config )
      volumetricFlowRateEvaluationConfigParser( config->getGlobalBlock(), configBlockName, defaultPlotFrequency, defaultLogFrequency,
                                                defaultFilename, defaultNormalizationFactor, defaultDomainNormalization, defaultAxis, defaultSurfacePoint );
}

} // namespace internal

#define WALBERLA_FIELD_MAKE_VOLUMETRIC_FLOW_RATE_EVALUATION_CONFIG_PARSER( config ) \
   uint_t defaultPlotFrequency = uint_t(0); \
   uint_t defaultLogFrequency = uint_t(0); \
   std::string defaultFilename = internal::volumetricFlowRateEvaluationFilename; \
   real_t defaultNormalizationFactor = internal::volumetricFlowRateEvaluationNormalization; \
   auto _blocks = blocks.lock(); \
   WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to execute 'makeVolumetricFlowRateEvaluation' for a block storage object that doesn't exist anymore" ); \
   Vector3<real_t> defaultDomainNormalization = _blocks->getDomain().sizes(); \
   Vector3<bool> defaultAxis = internal::volumetricFlowRateEvaluationAxis; \
   Vector3<real_t> defaultSurfacePoint = internal::volumetricFlowRateEvaluationPoint; \
   internal::volumetricFlowRateEvaluationConfigParser( config, configBlockName, defaultPlotFrequency, defaultLogFrequency, defaultFilename, \
                                                       defaultNormalizationFactor, defaultDomainNormalization, defaultAxis, defaultSurfacePoint );

#define WALBERLA_FIELD_MAKE_VOLUMETRIC_FLOW_RATE_EVALUATION_SET_AND_RETURN() \
   evaluation->setFilename( defaultFilename ); \
   evaluation->setNormalizationFactor( defaultNormalizationFactor ); \
   evaluation->setDomainNormalization( defaultDomainNormalization ); \
   evaluation->setAxis( defaultAxis ); \
   evaluation->setSurfacePoint( defaultSurfacePoint ); \
   return evaluation;

template< typename VelocityField_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< VolumetricFlowRateEvaluation< VelocityField_T > > makeVolumetricFlowRateEvaluation( const Config_T & config,
                                                                                                const weak_ptr< StructuredBlockStorage > & blocks,
                                                                                                const ConstBlockDataID & velocityFieldId,
                                                                                                const internal::FlowRateSolution_T & solution = internal::FlowRateSolution_T(),
                                                                                                const internal::FlowRateVelocitySolution_T & velocitySolution = internal::FlowRateVelocitySolution_T(),
                                                                                                const std::string & configBlockName = internal::volumetricFlowRateEvaluationConfigBlock,
                                                                                                const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                                                const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_VOLUMETRIC_FLOW_RATE_EVALUATION_CONFIG_PARSER( config )
   using FR_T = VolumetricFlowRateEvaluation<VelocityField_T>;
   auto evaluation = shared_ptr< FR_T >( new FR_T( blocks, velocityFieldId, defaultPlotFrequency, defaultLogFrequency,
                                                   solution, velocitySolution, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_VOLUMETRIC_FLOW_RATE_EVALUATION_SET_AND_RETURN()
}

template< typename VelocityField_T, typename FlagField_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< VolumetricFlowRateEvaluation< VelocityField_T, FlagFieldEvaluationFilter<FlagField_T> > >
makeVolumetricFlowRateEvaluation( const Config_T & config,
                                  const weak_ptr< StructuredBlockStorage > & blocks,
                                  const ConstBlockDataID & velocityFieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                                  const internal::FlowRateSolution_T & solution = internal::FlowRateSolution_T(),
                                  const internal::FlowRateVelocitySolution_T & velocitySolution = internal::FlowRateVelocitySolution_T(),
                                  const std::string & configBlockName = internal::volumetricFlowRateEvaluationConfigBlock,
                                  const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                  const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_VOLUMETRIC_FLOW_RATE_EVALUATION_CONFIG_PARSER( config )
   using FR_T = VolumetricFlowRateEvaluation<VelocityField_T, FlagFieldEvaluationFilter<FlagField_T>>;
   auto evaluation = shared_ptr< FR_T >( new FR_T( blocks, velocityFieldId, FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                                   defaultPlotFrequency, defaultLogFrequency, solution, velocitySolution, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_VOLUMETRIC_FLOW_RATE_EVALUATION_SET_AND_RETURN()
}

template< typename VelocityField_T, typename Filter_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< VolumetricFlowRateEvaluation< VelocityField_T, Filter_T > >
makeVolumetricFlowRateEvaluation( const Config_T & config,
                                  const weak_ptr< StructuredBlockStorage > & blocks,
                                  const ConstBlockDataID & velocityFieldId, const Filter_T & filter,
                                  const internal::FlowRateSolution_T & solution = internal::FlowRateSolution_T(),
                                  const internal::FlowRateVelocitySolution_T & velocitySolution = internal::FlowRateVelocitySolution_T(),
                                  const std::string & configBlockName = internal::volumetricFlowRateEvaluationConfigBlock,
                                  const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                  const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_VOLUMETRIC_FLOW_RATE_EVALUATION_CONFIG_PARSER( config )
   using FR_T = VolumetricFlowRateEvaluation<VelocityField_T, Filter_T>;
   auto evaluation = shared_ptr< FR_T >( new FR_T( blocks, velocityFieldId, filter, defaultPlotFrequency, defaultLogFrequency,
                                                   solution, velocitySolution, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_VOLUMETRIC_FLOW_RATE_EVALUATION_SET_AND_RETURN()
}



#undef WALBERLA_FIELD_MAKE_VOLUMETRIC_FLOW_RATE_EVALUATION_CONFIG_PARSER
#undef WALBERLA_FIELD_MAKE_VOLUMETRIC_FLOW_RATE_EVALUATION_SET_AND_RETURN

} // namespace field
} // namespace walberla
