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
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#define KEEP_WALBERLA_FIELD_MAKE_MASS_EVALUATION
#include "field/MassEvaluation.h"



namespace walberla {
namespace lbm {



namespace internal {

Vector3<real_t> massEvaluationDomain( const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level )
{
   return Vector3<real_t>( real_c( blocks->getNumberOfXCells(level) ),
                           real_c( blocks->getNumberOfYCells(level) ),
                           real_c( blocks->getNumberOfZCells(level) ) );
}

}



//**********************************************************************************************************************
/*!
*   \section docLBMMassEvaluation Mass Evaluation for the LBM
*
*   Some 'makeMassEvaluation' functions for evaluating the mass evolution of LBM-based simulations.
*   For an extensive documentation see \ref docMassEvaluation in 'field/MassEvaluation.h'.
*
*   Please note that the 'makeMassEvaluation' functions for the LBM pre-configure the domain normalization such that
*   the volume of a cell on level 'level' (0 by default) is equal to 1. This is in accordance with dimensionless
*   lattice units.
*
*   If you do not have a dedicated field for storing the density of your LBM-based simulation, you can use
*
*   \code
*   BlockDataID densityId = field::addFieldAdaptor< lbm::Adaptor<LatticeModel_T>::Density >( blocks, pdfFieldId,
*                                                                                            "density adaptor" );
*   \code
*
*   for creating a density adaptor. The type of this adaptor/field is 'lbm::Adaptor<LatticeModel_T>::Density'.
*   An example for using a SharedFunctor that can be used to register the mass evaluation at a time loop might look
*   like as follows:
*
*   \code
*   makeSharedFunctor( lbm::makeMassEvaluation< lbm::Adaptor<LatticeModel_T>::Density >( config, blocks,
*                                                                                        uint_t(0), densityId ) )
*   \code
*
*   Note that for this example the plot and log frequency can be controlled via the configuration file. In this example,
*   _all_ cells are processed. If not all of your cells are actually fluid cells, you should use a 'makeMassEvaluation'
*   function that only processes cells marked as fluid cells in a flag field:
*
*   \code
*   makeSharedFunctor( lbm::makeMassEvaluation< lbm::Adaptor<LatticeModel_T>::Density, FlagField_T >(
                          config, blocks, uint_t(0), densityId, flagFieldId, fluidFlagUID ) )
*   \code
*/
//**********************************************************************************************************************



/////////////////////////////////////////////////////////////
// makeMassEvaluation functions without configuration file //
/////////////////////////////////////////////////////////////

template< typename DensityField_T >
shared_ptr< walberla::field::MassEvaluation< DensityField_T > > makeMassEvaluation( const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                                                                                    const ConstBlockDataID & fieldId,
                                                                                    const uint_t plotFrequency, const uint_t logFrequency,
                                                                                    const std::string & filename = walberla::field::internal::massEvaluationFilename,
                                                                                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using ME_T = walberla::field::MassEvaluation<DensityField_T>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, plotFrequency, logFrequency, filename, requiredSelectors, incompatibleSelectors ) );
   evaluation->setDomainNormalization( internal::massEvaluationDomain( blocks, level ) );
   return evaluation;
}

template< typename DensityField_T, typename FlagField_T >
shared_ptr< walberla::field::MassEvaluation< DensityField_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T> > >
makeMassEvaluation( const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                    const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & fluid,
                    const uint_t plotFrequency, const uint_t logFrequency,
                    const std::string & filename = walberla::field::internal::massEvaluationFilename,
                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using ME_T = walberla::field::MassEvaluation<DensityField_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, fluid ),
                                                   plotFrequency, logFrequency, filename, requiredSelectors, incompatibleSelectors ) );
   evaluation->setDomainNormalization( internal::massEvaluationDomain( blocks, level ) );
   return evaluation;
}

template< typename DensityField_T, typename Filter_T >
shared_ptr< walberla::field::MassEvaluation< DensityField_T, Filter_T > > makeMassEvaluation( const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                                                                                              const ConstBlockDataID & fieldId, const Filter_T & filter,
                                                                                              const uint_t plotFrequency, const uint_t logFrequency,
                                                                                              const std::string & filename = walberla::field::internal::massEvaluationFilename,
                                                                                              const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                                              const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using ME_T = walberla::field::MassEvaluation<DensityField_T, Filter_T>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, filter, plotFrequency, logFrequency, filename, requiredSelectors, incompatibleSelectors ) );
   evaluation->setDomainNormalization( internal::massEvaluationDomain( blocks, level ) );
   return evaluation;
}

// for pseudo 2D:

template< typename DensityField_T, bool Pseudo2D >
shared_ptr< walberla::field::MassEvaluation< DensityField_T, field::DefaultEvaluationFilter, Pseudo2D > >
makeMassEvaluation( const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                    const ConstBlockDataID & fieldId,
                    const uint_t plotFrequency, const uint_t logFrequency,
                    const std::string & filename = walberla::field::internal::massEvaluationFilename,
                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using ME_T = walberla::field::MassEvaluation<DensityField_T, field::DefaultEvaluationFilter, Pseudo2D>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, plotFrequency, logFrequency, filename, requiredSelectors, incompatibleSelectors ) );
   evaluation->setDomainNormalization( internal::massEvaluationDomain( blocks, level ) );
   return evaluation;
}

template< typename DensityField_T, typename FlagField_T, bool Pseudo2D >
shared_ptr< walberla::field::MassEvaluation< DensityField_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, Pseudo2D > >
makeMassEvaluation( const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                    const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & fluid,
                    const uint_t plotFrequency, const uint_t logFrequency,
                    const std::string & filename = walberla::field::internal::massEvaluationFilename,
                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using ME_T = walberla::field::MassEvaluation<DensityField_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, Pseudo2D>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, fluid ),
                                                   plotFrequency, logFrequency, filename, requiredSelectors, incompatibleSelectors ) );
   evaluation->setDomainNormalization( internal::massEvaluationDomain( blocks, level ) );
   return evaluation;
}

template< typename DensityField_T, typename Filter_T, bool Pseudo2D >
shared_ptr< walberla::field::MassEvaluation< DensityField_T, Filter_T, Pseudo2D > >
makeMassEvaluation( const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                    const ConstBlockDataID & fieldId, const Filter_T & filter,
                    const uint_t plotFrequency, const uint_t logFrequency,
                    const std::string & filename = walberla::field::internal::massEvaluationFilename,
                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using ME_T = walberla::field::MassEvaluation<DensityField_T, Filter_T, Pseudo2D>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, filter, plotFrequency, logFrequency, filename, requiredSelectors, incompatibleSelectors ) );
   evaluation->setDomainNormalization( internal::massEvaluationDomain( blocks, level ) );
   return evaluation;
}



///////////////////////////////////////////////////////
// makeMassEvaluation functions + configuration file //
///////////////////////////////////////////////////////

#define WALBERLA_LBM_MAKE_MASS_EVALUATION_CONFIG_PARSER( config ) \
   uint_t defaultPlotFrequency = uint_t(0); \
   uint_t defaultLogFrequency = uint_t(0); \
   std::string defaultFilename = walberla::field::internal::massEvaluationFilename; \
   Vector3<real_t> defaultDomainNormalization( internal::massEvaluationDomain( blocks, level ) ); \
   walberla::field::internal::massEvaluationConfigParser( config, configBlockName, defaultPlotFrequency, defaultLogFrequency, defaultFilename, defaultDomainNormalization );

template< typename DensityField_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< walberla::field::MassEvaluation< DensityField_T > > makeMassEvaluation( const Config_T & config,
                                                                                    const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                                                                                    const ConstBlockDataID & fieldId,
                                                                                    const std::string & configBlockName = walberla::field::internal::massEvaluationConfigBlock,
                                                                                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_LBM_MAKE_MASS_EVALUATION_CONFIG_PARSER( config )
   using ME_T = walberla::field::MassEvaluation<DensityField_T>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, defaultPlotFrequency, defaultLogFrequency, defaultFilename, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN()
}

template< typename DensityField_T, typename FlagField_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< walberla::field::MassEvaluation< DensityField_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T> > >
makeMassEvaluation( const Config_T & config,
                    const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                    const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & fluid,
                    const std::string & configBlockName = walberla::field::internal::massEvaluationConfigBlock,
                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_LBM_MAKE_MASS_EVALUATION_CONFIG_PARSER( config )
   using ME_T = walberla::field::MassEvaluation<DensityField_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, fluid ),
                                                   defaultPlotFrequency, defaultLogFrequency, defaultFilename, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN()
}

template< typename DensityField_T, typename Filter_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< walberla::field::MassEvaluation< DensityField_T, Filter_T > > makeMassEvaluation( const Config_T & config,
                                                                                              const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                                                                                              const ConstBlockDataID & fieldId, const Filter_T & filter,
                                                                                              const std::string & configBlockName = walberla::field::internal::massEvaluationConfigBlock,
                                                                                              const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                                              const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_LBM_MAKE_MASS_EVALUATION_CONFIG_PARSER( config )
   using ME_T = walberla::field::MassEvaluation<DensityField_T, Filter_T>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, filter, defaultPlotFrequency, defaultLogFrequency, defaultFilename, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN()
}

// for pseudo 2D:

template< typename DensityField_T, bool Pseudo2D, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< walberla::field::MassEvaluation< DensityField_T, field::DefaultEvaluationFilter, Pseudo2D > >
makeMassEvaluation( const Config_T & config,
                    const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                    const ConstBlockDataID & fieldId,
                    const std::string & configBlockName = walberla::field::internal::massEvaluationConfigBlock,
                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_LBM_MAKE_MASS_EVALUATION_CONFIG_PARSER( config )
   using ME_T = walberla::field::MassEvaluation<DensityField_T, field::DefaultEvaluationFilter, Pseudo2D>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, defaultPlotFrequency, defaultLogFrequency, defaultFilename, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN()
}

template< typename DensityField_T, typename FlagField_T, bool Pseudo2D, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< walberla::field::MassEvaluation< DensityField_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, Pseudo2D > >
makeMassEvaluation( const Config_T & config,
                    const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                    const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & fluid,
                    const std::string & configBlockName = walberla::field::internal::massEvaluationConfigBlock,
                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_LBM_MAKE_MASS_EVALUATION_CONFIG_PARSER( config )
   using ME_T = walberla::field::MassEvaluation<DensityField_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, Pseudo2D>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, fluid ),
                                                   defaultPlotFrequency, defaultLogFrequency, defaultFilename, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN()
}

template< typename DensityField_T, typename Filter_T, bool Pseudo2D, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< walberla::field::MassEvaluation< DensityField_T, Filter_T, Pseudo2D > >
makeMassEvaluation( const Config_T & config,
                    const shared_ptr< StructuredBlockStorage > & blocks, const uint_t level,
                    const ConstBlockDataID & fieldId, const Filter_T & filter,
                    const std::string & configBlockName = walberla::field::internal::massEvaluationConfigBlock,
                    const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_LBM_MAKE_MASS_EVALUATION_CONFIG_PARSER( config )
   using ME_T = walberla::field::MassEvaluation<DensityField_T, Filter_T, Pseudo2D>;
   auto evaluation = shared_ptr< ME_T >( new ME_T( blocks, fieldId, filter, defaultPlotFrequency, defaultLogFrequency, defaultFilename, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN()
}



#undef WALBERLA_FIELD_MAKE_MASS_EVALUATION_CONFIG_PARSER
#undef WALBERLA_FIELD_MAKE_MASS_EVALUATION_SET_AND_RETURN
#undef KEEP_WALBERLA_FIELD_MAKE_MASS_EVALUATION

} // namespace lbm
} // namespace walberla
