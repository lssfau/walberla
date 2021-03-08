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
//! \file PdfFieldInitializer.h
//! \ingroup lbm
//! \author Tobias Schruff <tobias.schruff@gmail.com>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"
#include "core/config/Config.h"
#include "lbm/field/PdfField.h"

#include "ExprSystemInitFunction.h"


namespace walberla {
namespace lbm {
namespace initializer {


template< typename LatticeModel_T >
class PdfFieldInitializer
{
   
public:

   PdfFieldInitializer( const BlockDataID & pdfFieldId, const shared_ptr<StructuredBlockForest> & blocks );

   template< typename InitFunc >
   void initDensity( InitFunc & func ) const;

   template< typename InitFunc >
   void initVelocity( InitFunc & func ) const;

   template< typename InitFunc >
   void initDensityAndVelocity( InitFunc & func ) const;

   void initFromConfig( const Config::BlockHandle & config ) const;


private:

   using PdfField_T = lbm::PdfField<LatticeModel_T>;

   const BlockDataID                       pdfFieldId_;
   const shared_ptr<StructuredBlockForest> blocks_;
};

} // namespace initializer
} // namespace lbm
} // namespace walberla


#include "PdfFieldInitializer.impl.h"

