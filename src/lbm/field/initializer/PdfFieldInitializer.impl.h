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
//! \file PdfFieldInitializer.impl.h
//! \ingroup lbm
//! \author Tobias Schruff <tobias.schruff@gmail.com>
//
//======================================================================================================================


namespace walberla {
namespace lbm {
namespace initializer {


template< typename LatticeModel_T >
PdfFieldInitializer< LatticeModel_T >::PdfFieldInitializer( const BlockDataID & pdfFieldId, const shared_ptr<StructuredBlockForest> & blocks )
   : pdfFieldId_( pdfFieldId ), blocks_( blocks ) { }


template< typename LatticeModel_T >
template< typename InitFunc >
void PdfFieldInitializer< LatticeModel_T >::initDensity( InitFunc & func ) const
{
   for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
   {
      auto pdfField = blockIt->template getData<PdfField_T>( pdfFieldId_ );
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );

      for( auto cellIt = pdfField->beginWithGhostLayerXYZ(); cellIt != pdfField->end(); ++cellIt )
      {
         Cell globalCell( cellIt.cell() );
         blocks_->transformBlockLocalToGlobalCell( globalCell, *blockIt );

         const Vector3<real_t> velocity = pdfField->getVelocity( cellIt.cell() );
         pdfField->setDensityAndVelocity( cellIt.cell(), velocity, func( globalCell ) );
      }
   }
}


template< typename LatticeModel_T >
template< typename InitFunc >
void PdfFieldInitializer< LatticeModel_T >::initVelocity( InitFunc & func ) const
{
   for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
   {
      auto pdfField = blockIt->template getData<PdfField_T>( pdfFieldId_ );
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );

      for( auto cellIt = pdfField->beginWithGhostLayerXYZ(); cellIt != pdfField->end(); ++cellIt )
      {
         Cell globalCell( cellIt.cell() );
         blocks_->transformBlockLocalToGlobalCell( globalCell, *blockIt );
         
         const real_t density = pdfField->getDensity( cellIt.cell() );
         pdfField->setDensityAndVelocity( cellIt.cell(), func( globalCell ), density );
      }
   }
}


template< typename LatticeModel_T >
template< typename InitFunc >
void PdfFieldInitializer< LatticeModel_T >::initDensityAndVelocity( InitFunc & func ) const
{
   for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
   {
      auto pdfField = blockIt->template getData<PdfField_T>( pdfFieldId_ );
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );

      for( auto cellIt = pdfField->beginWithGhostLayerXYZ(); cellIt != pdfField->end(); ++cellIt )
      {
         Cell globalCell( cellIt.cell() );
         blocks_->transformBlockLocalToGlobalCell( globalCell, *blockIt );
         
         const std::vector<real_t> densityAndVelocity = func( globalCell );
         WALBERLA_ASSERT_EQUAL( densityAndVelocity.size(), 4 );

         const real_t           density( densityAndVelocity[0] );
         const Vector3<real_t> velocity( densityAndVelocity[1], densityAndVelocity[2], densityAndVelocity[3] );

         pdfField->setDensityAndVelocity( cellIt.cell(), velocity, density );
      }
   }
}


template< typename LatticeModel_T >
void PdfFieldInitializer< LatticeModel_T >::initFromConfig( const Config::BlockHandle & config ) const
{
   ExprSystemInitFunction initFunction( blocks_ );
   if( !initFunction.parse( config ) )
      return;

   initDensityAndVelocity( initFunction );
}


} // namespace initializer
} // namespace lbm
} // namespace walberla

