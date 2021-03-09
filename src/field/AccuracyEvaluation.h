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
//! \file AccuracyEvaluation.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/OpenMP.h"
#include "core/logging/Logging.h"
#include "core/debug/CheckFunctions.h"
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

const std::string accuracyEvaluationFilename("accuracy.dat");
const std::string accuracyEvaluationConfigBlock("AccuracyEvaluation");

template< typename T >
inline real_t accuracyEvaluationAbsError( const T & error ) { return std::abs( error ); }

template<>
inline real_t accuracyEvaluationAbsError( const Vector3<real_t> & error ) { return error.length(); }

}



//**********************************************************************************************************************
/*!
*   \brief Class for evaluating the accuracy of a simulation
*
*   \section docAccuracyEvaluation Accuracy Evaluation
*
*   Class for evaluating the accuracy of a simulation by comparing simulation values with values provided by a solution
*   function (can be the analytical solution, if one exists). Comparison is performed by evaluating discrete weighted
*   L1, L2, and Lmax norms of the error (error = difference between simulation and solution values). For each cell that
*   is evaluated, the weight corresponds to the volume fraction of the cell in relation to the entire domain.
*
*   Do not create objects of class AccuracyEvaluation directly, better use one of the various 'makeAccuracyEvaluation'
*   functions below!
*
*   Template parameters:
*   - Field_T: the field storing the simulation values (also works if the field stores data of type Vector3)
*   - SolutionFunction_T: type of the solution function - must return Field_T::value_type and must take one parameter
*                         of type Vector3<real_t> that corresponds to a position inside the simulation domain
*                         (simulation domain = returned by calling getDomain() at the block storage)
*   - Filter_T: the type of the evaluation filter (see \ref docEvaluationFilter in 'EvaluationFilter.h')
*
*   Parameters for setting up and controlling accuracy evaluation:
*   - blocks: the block storage
*   - fieldId: block data ID of the field
*   - solution: the solution callback function - must return the solution when called with a position inside the domain
*   - filter: the evaluation filter that indicates which cells are processed
*   - plot frequency: the plotting interval - used for saving the data to file.
*                     If set to '0', no plotting data is created.
*   - log frequency: the logging interval - used for logging the data via the Logging singleton.
*                    If set to '0', no logging is performed.
*   - filename: the name of the file that stores the data for plotting
*   - normalization factor: an optional factor the simulation values and the solution values are multiplied with
*   - required and incompatible selectors
*
*   You do not have to specify an evaluation filter! If you do not specify any filter, _all_ cells are processed and no
*   cell is excluded.
*
*   If you want to use a flag field as evaluation filter, fitting 'makeAccuracyEvaluation' functions already exist.
*   These functions need an additional template parameter FlagField_T and you have to provide the block data ID of the
*   flag field together with a set of flag UIDs that specify which cells need to be processed.
*
*   There also exist 'makeAccuracyEvaluation' functions that take configuration file data as an additional parameter in
*   order to parse the configuration file for setting up and controlling accuracy evaluation. The configuration file
*   block looks like as follows:
*
*   \code
*   AccuracyEvaluation
*   {
*      plotFrequency [unsigned integer]; // the plot frequency
*      logFrequency  [unsigned integer]; // the log frequency
*      filename      [string]; // the name of the file that stores the data for plotting
*      normalization [floating point value]; // normalization factor
*   }
*   \endcode
*
*   Example:
*
*   \code
*   AccuracyEvaluation
*   {
*      plotFrequency 10;
*      logFrequency  1000;
*      filename      Accuracy.txt;
*   }
*   \endcode
*
*   Note that the shared pointer returned by all 'makeAccuracyEvaluation' functions can be captured by a SharedFunctor
*   for immediate registration at a time loop (see field::makeSharedFunctor).
*/
//**********************************************************************************************************************

template< typename Field_T,
          typename SolutionFunction_T = std::function< typename Field_T::value_type ( const Vector3< real_t > & ) >,
          typename Filter_T = DefaultEvaluationFilter >
class AccuracyEvaluation
{
public:

   AccuracyEvaluation( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                       const SolutionFunction_T & solution, const Filter_T & filter,
                       const uint_t plotFrequency, const uint_t logFrequency,
                       const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                       const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), fieldId_( fieldId ), solution_( solution ), filter_( filter ),
      executionCounter_( uint_t(0) ), plotFrequency_( plotFrequency ), logFrequency_( logFrequency ),
      filename_( internal::accuracyEvaluationFilename ), normalizationFactor_( real_t(1) ),
      requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {}

   AccuracyEvaluation( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                       const SolutionFunction_T & solution,
                       const uint_t plotFrequency, const uint_t logFrequency,
                       const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                       const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), fieldId_( fieldId ), solution_( solution ), filter_( Filter_T() ),
      executionCounter_( uint_t(0) ), plotFrequency_( plotFrequency ), logFrequency_( logFrequency ),
      filename_( internal::accuracyEvaluationFilename ), normalizationFactor_( real_t(1) ),
      L1_( real_t(0) ), L2_( real_t(0) ), Lmax_( real_t(0) ),
      requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {
      static_assert( (std::is_same< Filter_T, DefaultEvaluationFilter >::value),
                     "This constructor is only available if DefaultEvaluationFilter is set as filter type!" );
   }

   void setNormalizationFactor( const real_t f ) { normalizationFactor_ = f; }
   void setFilename( const std::string & filename ) { filename_ = filename; }

   real_t L1() const { return L1_; }
   real_t L2() const { return L2_; }
   real_t Lmax() const { return Lmax_; }

   void operator()();

private:

   weak_ptr< StructuredBlockStorage > blocks_;
   ConstBlockDataID fieldId_;

   SolutionFunction_T solution_;
   Filter_T filter_;

   uint_t executionCounter_;

   uint_t plotFrequency_;
   uint_t  logFrequency_;

   std::string filename_;

   real_t normalizationFactor_;

   real_t L1_;
   real_t L2_;
   real_t Lmax_;

   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;

}; // class AccuracyEvaluation



template< typename Field_T, typename SolutionFunction_T, typename Filter_T >
void AccuracyEvaluation< Field_T, SolutionFunction_T, Filter_T >::operator()()
{
   if( logFrequency_ == uint_t(0) && ( plotFrequency_ == uint_t(0) || filename_.empty() ) )
      return;

   ++executionCounter_;

   const bool plot = ( plotFrequency_ != uint_t(0) && ( executionCounter_ - uint_c(1) ) % plotFrequency_ == uint_t(0) && !filename_.empty() );
   const bool log  = ( logFrequency_  != uint_t(0) && ( executionCounter_ - uint_c(1) ) % logFrequency_  == uint_t(0) );

   if( !log && !plot )
      return;

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'AccuracyEvaluation' for a block storage object that doesn't exist anymore" );

   const auto & domainAABB = blocks->getDomain();

   real_t _L1( real_t(0) );
   real_t _L2( real_t(0) );
   real_t _Lmax( real_t(0) );

   for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
   {
      const Field_T * field = block->template getData< const Field_T >( fieldId_ );

      const auto level = blocks->getLevel( *block );
      const real_t volumeFraction = ( blocks->dx( level ) * blocks->dy( level ) * blocks->dz( level ) ) / domainAABB.volume();

      filter_( *block );

#ifdef _OPENMP

      std::vector< real_t > lmax( numeric_cast<size_t>( omp_get_max_threads() ), real_t(0) );

      #pragma omp parallel
      {
         real_t & threadMax = lmax[ numeric_cast<size_t>( omp_get_thread_num() ) ];
         WALBERLA_FOR_ALL_CELLS_XYZ_OMP( field, omp for schedule(static) reduction(+:_L1) reduction(+:_L2),

            if( filter_(x,y,z) )
            {
               Vector3< real_t > center = blocks->getBlockLocalCellCenter( *block, Cell(x,y,z) );

               const auto error = field->get(x,y,z) - solution_( center );
               const real_t diff = internal::accuracyEvaluationAbsError( error ) * normalizationFactor_;

               _L1 += diff * volumeFraction;
               _L2 += diff * diff * volumeFraction;
               threadMax = std::max( threadMax, diff );
            }
         )
      }

      for( auto v = lmax.begin(); v != lmax.end(); ++v )
         _Lmax = std::max( _Lmax, *v );

#else

      WALBERLA_FOR_ALL_CELLS_XYZ( field,

         if( filter_(x,y,z) )
         {
            Vector3< real_t > center = blocks->getBlockLocalCellCenter( *block, Cell(x,y,z) );

            const auto error = field->get(x,y,z) - solution_( center );
            const real_t diff = internal::accuracyEvaluationAbsError( error ) * normalizationFactor_;
            
            _L1 += diff * volumeFraction;
            _L2 += diff * diff * volumeFraction;
            _Lmax = std::max( _Lmax, diff );
         }
      )

#endif

   }

   mpi::reduceInplace( _L1, mpi::SUM );
   mpi::reduceInplace( _L2, mpi::SUM );
   mpi::reduceInplace( _Lmax, mpi::MAX );
   _L2 = std::sqrt( _L2 );

   L1_ = _L1;
   L2_ = _L2;
   Lmax_ = _Lmax;

   WALBERLA_ROOT_SECTION()
   {
      const auto & id = blocks->getBlockDataIdentifier( fieldId_ );

      if( plot && executionCounter_ == uint_t(1) )
      {
         std::ofstream file( filename_.c_str() );
         file << "# accuracy evaluation of data '" << id <<  "'\n"
              << "# step [1], L1 [2], L2 [3], Lmax [4]" << std::endl;
         file.close();
      }

      if( log )
      {
         WALBERLA_LOG_INFO( "Evaluation of accuracy (weighted norms of the error) [data = '" << id << "']:" <<
                            "\n - L1:   " << L1_ <<
                            "\n - L2:   " << L2_ <<
                            "\n - Lmax: " << Lmax_ );
      }

      if( plot )
      {
         std::ofstream file( filename_.c_str(), std::ofstream::out | std::ofstream::app );
         file << ( executionCounter_ - uint_t(1) ) << " " << L1_ << " " << L2_ << " " << Lmax_ << std::endl;
         file.close();
      }
   }
}



/////////////////////////////////////////////////////////////////
// makeAccuracyEvaluation functions without configuration file //
/////////////////////////////////////////////////////////////////

template< typename Field_T, typename SolutionFunction_T >
shared_ptr< AccuracyEvaluation< Field_T, SolutionFunction_T > > makeAccuracyEvaluation( const weak_ptr< StructuredBlockStorage > & blocks,
                                                                                        const ConstBlockDataID & fieldId, const SolutionFunction_T & solution,
                                                                                        const uint_t plotFrequency, const uint_t logFrequency,
                                                                                        const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                                        const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using AE_T = AccuracyEvaluation<Field_T, SolutionFunction_T>;
   return shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, plotFrequency, logFrequency, requiredSelectors, incompatibleSelectors ) );
}

template< typename Field_T, typename FlagField_T, typename SolutionFunction_T >
shared_ptr< AccuracyEvaluation< Field_T, SolutionFunction_T, FlagFieldEvaluationFilter<FlagField_T> > >
makeAccuracyEvaluation( const weak_ptr< StructuredBlockStorage > & blocks,
                        const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                        const SolutionFunction_T & solution,
                        const uint_t plotFrequency, const uint_t logFrequency,
                        const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                        const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using AE_T = AccuracyEvaluation<Field_T, SolutionFunction_T, FlagFieldEvaluationFilter<FlagField_T>>;
   return shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                        plotFrequency, logFrequency, requiredSelectors, incompatibleSelectors ) );
}

template< typename Field_T, typename Filter_T, typename SolutionFunction_T >
shared_ptr< AccuracyEvaluation< Field_T, SolutionFunction_T, Filter_T > >
makeAccuracyEvaluation( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                        const SolutionFunction_T & solution, const Filter_T & filter,
                        const uint_t plotFrequency, const uint_t logFrequency,
                        const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                        const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using AE_T = AccuracyEvaluation<Field_T, SolutionFunction_T, Filter_T>;
   return shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, filter, plotFrequency, logFrequency, requiredSelectors, incompatibleSelectors ) );
}



///////////////////////////////////////////////////////////
// makeAccuracyEvaluation functions + configuration file //
///////////////////////////////////////////////////////////

namespace internal {

inline void accuracyEvaluationConfigParser( const Config::BlockHandle & parentBlockHandle, const std::string & configBlockName,
                                            uint_t & defaultPlotFrequency, uint_t & defaultLogFrequency,
                                            std::string & defaultFilename, real_t & defaultNormalizationFactor )
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
      }
   }
}

inline void accuracyEvaluationConfigParser( const shared_ptr< Config > & config, const std::string & configBlockName,
                                            uint_t & defaultPlotFrequency, uint_t & defaultLogFrequency,
                                            std::string & defaultFilename, real_t & defaultNormalizationFactor )
{
   if( !!config )
      accuracyEvaluationConfigParser( config->getGlobalBlock(), configBlockName,
                                      defaultPlotFrequency, defaultLogFrequency, defaultFilename, defaultNormalizationFactor );
}

} // namespace internal

#define WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_CONFIG_PARSER( config ) \
   uint_t defaultPlotFrequency = uint_t(0); \
   uint_t defaultLogFrequency = uint_t(0); \
   std::string defaultFilename = internal::accuracyEvaluationFilename; \
   real_t defaultNormalizationFactor = uint_t(1); \
   internal::accuracyEvaluationConfigParser( config, configBlockName, defaultPlotFrequency, defaultLogFrequency, defaultFilename, defaultNormalizationFactor );

#define WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_SET_AND_RETURN() \
   evaluation->setNormalizationFactor( defaultNormalizationFactor ); \
   evaluation->setFilename( defaultFilename ); \
   return evaluation;

template< typename Field_T, typename SolutionFunction_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< AccuracyEvaluation< Field_T, SolutionFunction_T > > makeAccuracyEvaluation( const Config_T & config,
                                                                                        const weak_ptr< StructuredBlockStorage > & blocks,
                                                                                        const ConstBlockDataID & fieldId, const SolutionFunction_T & solution,
                                                                                        const std::string & configBlockName = internal::accuracyEvaluationConfigBlock,
                                                                                        const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                                        const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_CONFIG_PARSER( config )
   using AE_T = AccuracyEvaluation<Field_T, SolutionFunction_T>;
   auto evaluation = shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, defaultPlotFrequency, defaultLogFrequency, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_SET_AND_RETURN()
}

template< typename Field_T, typename FlagField_T, typename SolutionFunction_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< AccuracyEvaluation< Field_T, SolutionFunction_T, FlagFieldEvaluationFilter<FlagField_T> > >
makeAccuracyEvaluation( const Config_T & config,
                        const weak_ptr< StructuredBlockStorage > & blocks,
                        const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                        const SolutionFunction_T & solution,
                        const std::string & configBlockName = internal::accuracyEvaluationConfigBlock,
                        const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                        const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_CONFIG_PARSER( config )
   using AE_T = AccuracyEvaluation<Field_T, SolutionFunction_T, FlagFieldEvaluationFilter<FlagField_T>>;
   auto evaluation = shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                                   defaultPlotFrequency, defaultLogFrequency, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_SET_AND_RETURN()
}

template< typename Field_T, typename Filter_T, typename SolutionFunction_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< AccuracyEvaluation< Field_T, SolutionFunction_T, Filter_T > >
makeAccuracyEvaluation( const Config_T & config,
                        const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                        const SolutionFunction_T & solution, const Filter_T & filter,
                        const std::string & configBlockName = internal::accuracyEvaluationConfigBlock,
                        const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                        const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_CONFIG_PARSER( config )
   using AE_T = AccuracyEvaluation<Field_T, SolutionFunction_T, Filter_T>;
   auto evaluation = shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, filter,
                                                   defaultPlotFrequency, defaultLogFrequency, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_SET_AND_RETURN()
}



#undef WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_CONFIG_PARSER
#undef WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_SET_AND_RETURN

} // namespace field
} // namespace walberla
