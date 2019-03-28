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
//! \file MassEvaluation.h
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
const std::string massEvaluationFilename("mass.dat");
const std::string massEvaluationConfigBlock("MassEvaluation");
}



//**********************************************************************************************************************
/*!
*   \brief Class for evaluating the evolution of the mass of a simulation
*
*   \section docMassEvaluation Mass Evaluation
*
*   Do not create objects of class MassEvaluation directly, better use one of the various 'makeMassEvaluation'
*   functions below!
*
*   Template parameters:
*   - DensityField_T: a scalar field that stores the density. If the density is not explicitly stored in a field but
*                     can be calculated from the data stored in another field, you do not have to create a density
*                     field but you can use a field adaptor.
*   - Filter_T: the type of the evaluation filter (see \ref docEvaluationFilter in 'EvaluationFilter.h')
*
*   Parameters for setting up and controlling mass evaluation:
*   - blocks: the block storage
*   - fieldId: block data ID of the density field (might be a field adaptor)
*   - filter: the evaluation filter that indicates which cells are processed
*   - plot frequency: the plotting interval - used for saving the data to file.
*                     If set to '0', no plotting data is created.
*   - log frequency: the logging interval - used for logging the data via the Logging singleton.
*                    If set to '0', no logging is performed.
*   - filename: the name of the file that stores the data for plotting
*   - domain normalization: In order to evaluate the mass, the density of a cell is multiplied with the volume of the
*                           cell. By default, the volume of a cell corresponds the the simulation space as given by
*                           the domain bounding box (blocks->getDomain()). However, you can overwrite the size of the
*                           domain via the domain normalization parameter. If you do so, the volume of each cell will
*                           correspond to this 'normalized' domain.
*   - required and incompatible selectors
*
*   You do not have to specify an evaluation filter! If you do not specify any filter, _all_ cells are processed and no
*   cell is excluded.
*
*   If you want to use a flag field as evaluation filter, fitting 'makeMassEvaluation' functions already exist. These
*   functions need an additional template parameter FlagField_T and you have to provide the block data ID of the flag
*   field together with a set of flag UIDs that specify which cells need to be processed.
*
*   There also exist 'makeMassEvaluation' functions that take configuration file data as an additional parameter in
*   order to parse the configuration file for setting up and controlling mass evaluation. The configuration file block
*   looks like as follows:
*
*   \code
*   MassEvaluation
*   {
*      plotFrequency [unsigned integer]; // the plot frequency
*      logFrequency  [unsigned integer]; // the log frequency
*      filename      [string]; // the name of the file that stores the data for plotting
*      domain        [Vector3: <real,real,real>]; // domain normalization
*   }
*   \endcode
*
*   Example:
*
*   \code
*   MassEvaluation
*   {
*      plotFrequency 10;
*      logFrequency  1000;
*      filename      Mass.txt;
*      domain        <1,1,1>;
*   }
*   \endcode
*
*   Note that the shared pointer returned by all 'makeMassEvaluation' functions can be captured by a SharedFunctor for
*   immediate registration at a time loop (see field::makeSharedFunctor).
*/
//**********************************************************************************************************************

template< typename DensityField_T, typename Filter_T = DefaultEvaluationFilter, bool Pseudo2D = false >
class MassEvaluation
{
public:

   MassEvaluation( const weak_ptr< StructuredBlockStorage > & blocks,
                   const ConstBlockDataID & fieldId, const Filter_T & filter,
                   const uint_t plotFrequency, const uint_t logFrequency,
                   const std::string & filename = internal::massEvaluationFilename,
                   const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                   const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), filter_( filter ),
      executionCounter_( uint_c(0) ), plotFrequency_( plotFrequency ), logFrequency_( logFrequency ), filename_( filename ),
      fieldId_( fieldId ),
      initialMass_( real_t(0) ), minMass_( real_t(0) ), maxMass_( real_t(0) ),
      requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {
      auto _blocks = blocks.lock();
      WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to access 'MassEvaluation' for a block storage object that doesn't exist anymore" );
      domainNormalization_ = _blocks->getDomain().sizes();
   }

   MassEvaluation( const weak_ptr< StructuredBlockStorage > & blocks,
                   const ConstBlockDataID & fieldId,
                   const uint_t plotFrequency, const uint_t logFrequency,
                   const std::string & filename = internal::massEvaluationFilename,
                   const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                   const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), filter_( Filter_T() ),
      executionCounter_( uint_c(0) ), plotFrequency_( plotFrequency ), logFrequency_( logFrequency ), filename_( filename ),
      fieldId_( fieldId ),
      initialMass_( real_t(0) ), minMass_( real_t(0) ), maxMass_( real_t(0) ),
      requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {
      static_assert( (std::is_same< Filter_T, DefaultEvaluationFilter >::value),
                     "This constructor is only available if DefaultEvaluationFilter is set as filter type!" );

      auto _blocks = blocks.lock();
      WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to access 'MassEvaluation' for a block storage object that doesn't exist anymore" );
      domainNormalization_ = _blocks->getDomain().sizes();
   }
   
   void setDomainNormalization( const Vector3<real_t> & d ) { domainNormalization_ = d; }
   
   void operator()();
   
private:

   weak_ptr< StructuredBlockStorage > blocks_;
   
   Filter_T filter_;

   uint_t executionCounter_;

   uint_t plotFrequency_;
   uint_t  logFrequency_;

   std::string filename_;

   ConstBlockDataID fieldId_;
   
   Vector3<real_t> domainNormalization_;

   real_t initialMass_;
   real_t minMass_;
   real_t maxMass_;

   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;

}; // class MassEvaluation



template< typename DensityField_T, typename Filter_T, bool Pseudo2D >
void MassEvaluation< DensityField_T, Filter_T, Pseudo2D >::operator()()
{
   if( logFrequency_ == uint_t(0) && ( plotFrequency_ == uint_t(0) || filename_.empty() ) )
      return;

   ++executionCounter_;

   const bool plot = ( plotFrequency_ != uint_t(0) && ( executionCounter_ - uint_c(1) ) % plotFrequency_ == uint_t(0) && !filename_.empty() );
   const bool log  = ( logFrequency_  != uint_t(0) && ( executionCounter_ - uint_c(1) ) % logFrequency_  == uint_t(0) );

   if( !log && !plot )
      return;

   real_t mass( real_t(0) );

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'MassEvaluation' for a block storage object that doesn't exist anymore" );

   const auto & domain = blocks->getDomain();
   for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
   {
      const DensityField_T * const density = block->template getData< const DensityField_T  >( fieldId_ );

      const uint_t level = blocks->getLevel( *block );
      const real_t volume = ( blocks->dx( level ) * domainNormalization_[0] / domain.xSize() ) *
                            ( blocks->dy( level ) * domainNormalization_[1] / domain.ySize() ) *
                            ( Pseudo2D ? domainNormalization_[2] : ( blocks->dz( level ) * domainNormalization_[2] / domain.zSize() ) );

      filter_( *block );
      
      WALBERLA_FOR_ALL_CELLS_XYZ_OMP( density, omp parallel for schedule(static) reduction(+:mass),

         if( filter_(x,y,z) )
            mass += density->get(x,y,z) * volume;

      ) // WALBERLA_FOR_ALL_CELLS_XYZ_OMP
   }

   mpi::reduceInplace( mass, mpi::SUM );

   WALBERLA_ROOT_SECTION()
   {
      const auto & id = blocks->getBlockDataIdentifier( fieldId_ );

      if( executionCounter_ == uint_t(1) )
      {
         initialMass_ = mass;
         minMass_ = mass;
         maxMass_ = mass;

         if( plot )
         {
            std::ofstream file( filename_.c_str() );
            file << "# mass evaluation of data '" << id <<  "'\n"
                 << "# step [1], mass [2], min. mass [3], max. mass [4], current deviation from initial mass [5], max. deviation from initial mass [6]" << std::endl;
            file.close();
         }
      }

      minMass_ = std::min( minMass_, mass );
      maxMass_ = std::max( maxMass_, mass );

      const real_t currentDeviation = std::abs( mass - initialMass_ ) / std::abs( initialMass_ );
      const real_t maxDeviation = std::max( std::abs( minMass_ - initialMass_ ), std::abs( maxMass_ - initialMass_ ) ) / std::abs( initialMass_ );

      if( log )
      {
         WALBERLA_LOG_INFO( "Evaluation of mass [data = '" << id << "']:" <<
                            "\n - current mass: " << mass <<
                            "\n -     min mass: " << minMass_ << " (since beginning of evaluation)" <<
                            "\n -     max mass: " << maxMass_ << " (since beginning of evaluation)" <<
                            "\n - initial mass: " << initialMass_ <<
                            "\n - current deviation from initial mass: " << ( real_c(100) * currentDeviation ) << " %" <<
                            "\n -     max deviation from initial mass: " << ( real_c(100) * maxDeviation ) << " %" );
      }

      if( plot )
      {
         std::ofstream file( filename_.c_str(), std::ofstream::out | std::ofstream::app );
         file << ( executionCounter_ - uint_t(1) ) << " " << mass << " " << minMass_ << " " << maxMass_
                                                   << " " << currentDeviation << " " << maxDeviation << std::endl;
         file.close();
      }
   }
}



/////////////////////////////////////////////////////////////
// makeMassEvaluation functions without configuration file //
/////////////////////////////////////////////////////////////

template< typename DensityField_T >
shared_ptr< MassEvaluation< DensityField_T > > makeMassEvaluation( const weak_ptr< StructuredBlockStorage > & blocks,
                                                                   const ConstBlockDataID & fieldId,
                                                                   const uint_t plotFrequency, const uint_t logFrequency,
                                                                   const std::string & filename = internal::massEvaluationFilename,
                                                                   const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                   const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   typedef MassEvaluation< DensityField_T > ME_T;
   return shared_ptr< ME_T >( new ME_T( blocks, fieldId, plotFrequency, logFrequency, filename, requiredSelectors, incompatibleSelectors ) );
}

template< typename DensityField_T, typename FlagField_T >
shared_ptr< MassEvaluation< DensityField_T, FlagFieldEvaluationFilter<FlagField_T> > >
makeMassEvaluation( const weak_ptr< StructuredBlockStorage > & blocks,
                    const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                    const uint_t plotFrequency, const uint_t logFrequency,
                    const std::string & filename = internal::massEvaluationFilename,
                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   typedef MassEvaluation< DensityField_T, FlagFieldEvaluationFilter<FlagField_T> > ME_T;
   return shared_ptr< ME_T >( new ME_T( blocks, fieldId, FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                        plotFrequency, logFrequency, filename, requiredSelectors, incompatibleSelectors ) );
}

template< typename DensityField_T, typename Filter_T >
shared_ptr< MassEvaluation< DensityField_T, Filter_T > > makeMassEvaluation( const weak_ptr< StructuredBlockStorage > & blocks,
                                                                             const ConstBlockDataID & fieldId, const Filter_T & filter,
                                                                             const uint_t plotFrequency, const uint_t logFrequency,
                                                                             const std::string & filename = internal::massEvaluationFilename,
                                                                             const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                             const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   typedef MassEvaluation< DensityField_T, Filter_T > ME_T;
   return shared_ptr< ME_T >( new ME_T( blocks, fieldId, filter, plotFrequency, logFrequency, filename, requiredSelectors, incompatibleSelectors ) );
}



///////////////////////////////////////////////////////
// makeMassEvaluation functions + configuration file //
///////////////////////////////////////////////////////

namespace internal {

inline void massEvaluationConfigParser( const Config::BlockHandle & parentBlockHandle, const std::string & configBlockName,
                                        uint_t & defaultPlotFrequency, uint_t & defaultLogFrequency,
                                        std::string & defaultFilename, Vector3<real_t> & defaultDomainNormalization )
{
   if( parentBlockHandle )
   {
      Config::BlockHandle block = parentBlockHandle.getBlock( configBlockName );
      if( block )
      {
         defaultPlotFrequency = block.getParameter< uint_t >( "plotFrequency", defaultPlotFrequency );
         defaultLogFrequency = block.getParameter< uint_t >( "logFrequency", defaultLogFrequency );
         defaultFilename = block.getParameter< std::string >( "filename", defaultFilename );
         defaultDomainNormalization = block.getParameter< Vector3<real_t> >( "domain", defaultDomainNormalization );
      }
   }
}

inline void massEvaluationConfigParser( const shared_ptr< Config > & config, const std::string & configBlockName,
                                        uint_t & defaultPlotFrequency, uint_t & defaultLogFrequency,
                                        std::string & defaultFilename, Vector3<real_t> & defaultDomainNormalization )
{
   if( !!config )
      massEvaluationConfigParser( config->getGlobalBlock(), configBlockName,
                                  defaultPlotFrequency, defaultLogFrequency, defaultFilename, defaultDomainNormalization );
}

} // namespace internal

#define WALBERLA_FIELD_MAKE_MASS_EVALUATION_CONFIG_PARSER( config ) \
   uint_t defaultPlotFrequency = uint_t(0); \
   uint_t defaultLogFrequency = uint_t(0); \
   std::string defaultFilename = internal::massEvaluationFilename; \
   auto _blocks = blocks.lock(); \
   WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to execute 'makeMassEvaluation' for a block storage object that doesn't exist anymore" ); \
   Vector3<real_t> defaultDomainNormalization( _blocks->getDomain().sizes() ); \
   internal::massEvaluationConfigParser( config, configBlockName, defaultPlotFrequency, defaultLogFrequency, defaultFilename, defaultDomainNormalization );

#define WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN() \
   evaluation->setDomainNormalization( defaultDomainNormalization ); \
   return evaluation;

template< typename DensityField_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< MassEvaluation< DensityField_T > > makeMassEvaluation( const Config_T & config,
                                                                   const weak_ptr< StructuredBlockStorage > & blocks,
                                                                   const ConstBlockDataID & fieldId,
                                                                   const std::string & configBlockName = internal::massEvaluationConfigBlock,
                                                                   const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                   const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_CONFIG_PARSER( config )
   typedef MassEvaluation< DensityField_T > ME_T;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, defaultPlotFrequency, defaultLogFrequency, defaultFilename,
                                                   requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN()
}

template< typename DensityField_T, typename FlagField_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< MassEvaluation< DensityField_T, FlagFieldEvaluationFilter<FlagField_T> > >
makeMassEvaluation( const Config_T & config,
                    const weak_ptr< StructuredBlockStorage > & blocks,
                    const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                    const std::string & configBlockName = internal::massEvaluationConfigBlock,
                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_CONFIG_PARSER( config )
   typedef MassEvaluation< DensityField_T, FlagFieldEvaluationFilter<FlagField_T> > ME_T;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                                   defaultPlotFrequency, defaultLogFrequency, defaultFilename, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN()
}

template< typename DensityField_T, typename Filter_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< MassEvaluation< DensityField_T, Filter_T > > makeMassEvaluation( const Config_T & config,
                                                                             const weak_ptr< StructuredBlockStorage > & blocks,
                                                                             const ConstBlockDataID & fieldId, const Filter_T & filter,
                                                                             const std::string & configBlockName = internal::massEvaluationConfigBlock,
                                                                             const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                             const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_CONFIG_PARSER( config )
   typedef MassEvaluation< DensityField_T, Filter_T > ME_T;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, filter, defaultPlotFrequency, defaultLogFrequency, defaultFilename,
                                                   requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN()
}


#ifndef KEEP_WALBERLA_FIELD_MAKE_MASS_EVALUATION
#undef WALBERLA_FIELD_MAKE_MASS_EVALUATION_CONFIG_PARSER
#undef WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN
#endif

} // namespace field
} // namespace walberla
