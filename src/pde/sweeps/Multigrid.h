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
//! \file Multigrid.h
//! \ingroup pde
//! \author Dominik Bartuschat <dominik.bartuschat@fau.de>
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include "field/GhostLayerField.h"



namespace walberla {
namespace pde {



//**********************************************************************************************************************
/*!
 *   \brief Restriction sweep for multigrid
 *
 *   \tparam Stencil_T The stencil used for the discrete operator
 */
//**********************************************************************************************************************
template< typename Stencil_T >
class Restrict
{
public:

   typedef GhostLayerField< real_t, 1 > Field_T;

   //*******************************************************************************************************************
   /* \param blocks the block storage where the fields are stored
    * \param fineFieldId the block data id of the fine field
    * \param coarseFieldId the block data id of the coarse field
    *******************************************************************************************************************/
   Restrict( const shared_ptr< domain_decomposition::StructuredBlockStorage > & blocks,
                const BlockDataID & fineFieldId, const BlockDataID & coarseFieldId ) :
      blocks_( blocks ), fineFieldId_( fineFieldId ), coarseFieldId_( coarseFieldId ) {}

   void operator()( IBlock * const block ) const;

private:

   shared_ptr< domain_decomposition::StructuredBlockStorage > blocks_;
   BlockDataID fineFieldId_;
   BlockDataID coarseFieldId_;
};



//**********************************************************************************************************************
/*!
 *   \brief Prolongation and correction sweep for multigrid
 *
 *   \tparam Stencil_T The stencil used for the discrete operator
 */
//**********************************************************************************************************************
template< typename Stencil_T >
class ProlongateAndCorrect
{
public:

   typedef GhostLayerField< real_t, 1 > Field_T;

   //*******************************************************************************************************************
   /* \param blocks the block storage where the fields are stored
    * \param coarseFieldId the block data id of the coarse field
    * \param fineFieldId the block data id of the fine field
    *******************************************************************************************************************/
   ProlongateAndCorrect( const shared_ptr< domain_decomposition::StructuredBlockStorage > & blocks,
                 const BlockDataID & coarseFieldId, const BlockDataID & fineFieldId ) :
   blocks_( blocks ), fineFieldId_( fineFieldId ), coarseFieldId_( coarseFieldId ) {}

   void operator()( IBlock * const block ) const;

private:

   shared_ptr< domain_decomposition::StructuredBlockStorage > blocks_;
   BlockDataID fineFieldId_;
   BlockDataID coarseFieldId_;
};



//**********************************************************************************************************************
/*!
 *   \brief Residual calculation sweep for multigrid with stencil field
 *
 *   \tparam Stencil_T The stencil used for the discrete operator
 */
//**********************************************************************************************************************
template< typename Stencil_T >
class ComputeResidual
{
public:

   typedef GhostLayerField< real_t, 1 >                Field_T;
   typedef GhostLayerField< real_t, Stencil_T::Size >  StencilField_T;

   //*******************************************************************************************************************
   /* \param blocks the block storage where the fields are stored
    * \param uId the block data id of the solution field
    * \param fId the block data id of the right-hand side field
    * \param stencilId the block data id of the stencil field
    * \param rId the block data id of the residual field
    *******************************************************************************************************************/
   ComputeResidual( const shared_ptr< domain_decomposition::StructuredBlockStorage > & blocks, const BlockDataID & uId,
                    const BlockDataID & fId, const BlockDataID & stencilId, const BlockDataID & rId )
            : blocks_( blocks ), uId_( uId ), fId_( fId ), stencilId_( stencilId ), rId_( rId )
   { }

   void operator()( IBlock * const block ) const;

private:

   shared_ptr< domain_decomposition::StructuredBlockStorage > blocks_;
   BlockDataID uId_;
   BlockDataID fId_;
   BlockDataID stencilId_;
   BlockDataID rId_;
};



//**********************************************************************************************************************
/*!
 *   \brief Residual calculation sweep for multigrid with fixed stencil weights
 *
 *   \tparam Stencil_T The stencil used for the discrete operator
 */
//**********************************************************************************************************************
template< typename Stencil_T >
class ComputeResidualFixedStencil
{
public:

   typedef GhostLayerField< real_t, 1 > Field_T;
   typedef std::vector< real_t  > Weight_T;

   //*******************************************************************************************************************
   /* \param blocks the block storage where the fields are stored
    * \param uId the block data id of the solution field
    * \param fId the block data id of the right-hand side field
    * \param weights vector of stencil weights for the discrete operator
    * \param rId the block data id of the residual field
    *******************************************************************************************************************/
   ComputeResidualFixedStencil( const shared_ptr< domain_decomposition::StructuredBlockStorage > & blocks, const BlockDataID & uId,
                                const BlockDataID & fId, const std::vector <real_t> & weights, const BlockDataID & rId )
            : blocks_( blocks ), uId_( uId ), fId_( fId ), weights_( weights ), rId_( rId )
   { }

   void operator()( IBlock * const block ) const;

private:

   shared_ptr< domain_decomposition::StructuredBlockStorage > blocks_;
   BlockDataID uId_;
   BlockDataID fId_;
   Weight_T weights_;
   BlockDataID rId_;
};



//**********************************************************************************************************************
/*!
 *   \brief Sweep that sets all values in a field to zero
 */
//**********************************************************************************************************************
class Zeroize
{
public:

   typedef GhostLayerField< real_t, 1 > Field_T;

   //*******************************************************************************************************************
   /* \param blocks the block storage where the field is stored
    * \param fieldId the block data id of the field
    *******************************************************************************************************************/
   Zeroize( const shared_ptr< domain_decomposition::StructuredBlockStorage > & blocks, const BlockDataID & fieldId )
            : blocks_( blocks ), fieldId_( fieldId )
   { }

   void operator()( IBlock * const block ) const
   {
      block->getData< Field_T >( fieldId_ )->setWithGhostLayer( real_t(0) );
   }

private:

   shared_ptr< domain_decomposition::StructuredBlockStorage > blocks_;
   BlockDataID fieldId_;
};



//**********************************************************************************************************************
/*!
 *   \brief Direct Coarsening Approach for the stencil field
 *
 *   \tparam Stencil_T The stencil used for the discrete operator
 */
//**********************************************************************************************************************
template< typename Stencil_T >
class CoarsenStencilFieldsDCA
{
public:

   typedef GhostLayerField< real_t, Stencil_T::Size >  StencilField_T;

   //*******************************************************************************************************************
   /* \param blocks the block storage where the fields are stored
    * \param numLvl number of grid levels to use (including the finest level)
    * \param operatorOrder the order of the (continuum) differential operator, e.g. 2 for Laplace
    *******************************************************************************************************************/
   CoarsenStencilFieldsDCA( shared_ptr< StructuredBlockForest > blocks,
                      const uint_t numLvl, const uint_t operatorOrder,
                      const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                      const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
         : blocks_( blocks ), numLvl_(numLvl), operatorOrder_(operatorOrder),
           requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
   { }

   //*******************************************************************************************************************
   /* \param stencilFieldId a vector of the block data ids of the stencil field for all levels (finest first)
    *******************************************************************************************************************/
   void operator()(const std::vector<BlockDataID> & stencilFieldId) const;

private:

   shared_ptr< domain_decomposition::StructuredBlockStorage > blocks_;

   uint_t numLvl_;
   uint_t operatorOrder_;

   Set< SUID > requiredSelectors_;
   Set< SUID > incompatibleSelectors_;

};



//**********************************************************************************************************************
/*!
 *   \brief Galerkin Coarsening Approach for the stencil field
 *
 *   \tparam Stencil_T The stencil used for the discrete operator
 */
//**********************************************************************************************************************
template< typename Stencil_T >
class CoarsenStencilFieldsGCA
{
public:

   typedef GhostLayerField< real_t, Stencil_T::Size >  StencilField_T;

   //*******************************************************************************************************************
   /* \param blocks the block storage where the fields are stored
    * \param numLvl number of grid levels to use (including the finest level)
    * \param overrelaxFact overrelaxation factor, e.g. 2 for the Poisson equation
    *******************************************************************************************************************/
   CoarsenStencilFieldsGCA( shared_ptr< StructuredBlockForest > blocks,
                      const uint_t numLvl, const real_t overrelaxFact = real_t(1),
                      const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                      const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
         : blocks_( blocks ), numLvl_(numLvl), overrelaxFact_(overrelaxFact),
           requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
   { }

   //*******************************************************************************************************************
   /* \param stencilFieldId a vector of the block data ids of the stencil field for all levels (finest first)
    *******************************************************************************************************************/
   void operator()(const std::vector<BlockDataID> & stencilFieldId) const;

private:

   shared_ptr< domain_decomposition::StructuredBlockStorage > blocks_;

   uint_t numLvl_;
   real_t overrelaxFact_;

   Set< SUID > requiredSelectors_;
   Set< SUID > incompatibleSelectors_;

};


} // namespace pde
} // namespace walberla

#include "Multigrid.impl.h"
