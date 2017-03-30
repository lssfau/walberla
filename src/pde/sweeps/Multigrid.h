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



template< typename Stencil_T >
class Restrict
{
public:

   typedef GhostLayerField< real_t, 1 > Field_T;

   Restrict( const shared_ptr< domain_decomposition::StructuredBlockStorage > & blocks,
                const BlockDataID & fineFieldId, const BlockDataID & coarseFieldId ) :
      blocks_( blocks ), fineFieldId_( fineFieldId ), coarseFieldId_( coarseFieldId ) {}

   void operator()( IBlock * const block ) const;

private:

   shared_ptr< domain_decomposition::StructuredBlockStorage > blocks_;
   BlockDataID fineFieldId_;
   BlockDataID coarseFieldId_;
};



template< typename Stencil_T >
class ProlongateAndCorrect
{
public:

   typedef GhostLayerField< real_t, 1 > Field_T;

   ProlongateAndCorrect( const shared_ptr< domain_decomposition::StructuredBlockStorage > & blocks,
                 const BlockDataID & coarseFieldId, const BlockDataID & fineFieldId ) :
   blocks_( blocks ), fineFieldId_( fineFieldId ), coarseFieldId_( coarseFieldId ) {}

   void operator()( IBlock * const block ) const;

private:

   shared_ptr< domain_decomposition::StructuredBlockStorage > blocks_;
   BlockDataID fineFieldId_;
   BlockDataID coarseFieldId_;
};



template< typename Stencil_T >
class ComputeResidual
{
public:

   typedef GhostLayerField< real_t, 1 >                Field_T;
   typedef GhostLayerField< real_t, Stencil_T::Size >  StencilField_T;

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



template< typename Stencil_T >
class ComputeResidualFixedStencil
{
public:

   typedef GhostLayerField< real_t, 1 > Field_T;
   typedef std::vector< real_t  > Weight_T;

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



class Zeroize
{
public:

   typedef GhostLayerField< real_t, 1 > Field_T;

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



} // namespace pde
} // namespace walberla

#include "Multigrid.impl.h"
