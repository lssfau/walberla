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
//! \file PdfFieldDisplayAdaptor.impl.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "PdfFieldDisplayAdaptor.h"

#include "core/debug/Debug.h"

#include "gui/BlockSliceView/CellView.h"

#include "stencil/D2Q9.h"
#include "stencil/Directions.h"


namespace walberla {
namespace lbm {

   template< typename field_t, typename St>
   PdfFieldDisplayAdaptor<field_t, St>::PdfFieldDisplayAdaptor( ConstBlockDataID scalarFieldID )
         : gui::FieldDisplayAdaptor<field_t>( scalarFieldID ),
           displayProperties_ ( NULL )
   {
   }

   template< typename field_t, typename St>
   PdfFieldDisplayAdaptor<field_t, St>::~PdfFieldDisplayAdaptor()
   {
      if ( displayProperties_)
         delete displayProperties_;
   }

   template< typename field_t, typename St>
   void PdfFieldDisplayAdaptor<field_t, St>::addConfigurationItem( QTreeWidgetItem * parentItem )
   {
      if ( ! displayProperties_ )
      {
         QStringList options;
         options << "Text";

         displayProperties_ = new gui::DisplayPropertiesItem( options, parentItem );

         QObject::connect( displayProperties_, SIGNAL( optionChanged() ),
                           this,               SIGNAL( requestRedraw() ) );
      }
   }


   template< typename field_t, typename St>
   void PdfFieldDisplayAdaptor<field_t, St>::configure( const IBlock * block, int sliceDim, cell_idx_t sliceValue,
                                       Vector3<uint_t> & innerSize, cell_idx_t & ghostLayers )
   {
      // Set name of tree item
      WALBERLA_ASSERT_NOT_NULLPTR( displayProperties_ );
      const std::string & name = block->getBlockStorage().getBlockDataIdentifier( blockDataId_ );
      displayProperties_->setText( 0, QString::fromStdString(name) );

      gui::FieldDisplayAdaptor<field_t>::configure( block, sliceDim, sliceValue, innerSize, ghostLayers );
   }

   template< typename field_t, typename Stencil>
   void PdfFieldDisplayAdaptor<field_t, Stencil>::draw( QVector<QVector< gui::CellView*> > & grid, int nrOfGhostLayers )
   {
      using namespace stencil;
      WALBERLA_ASSERT_GREATER_EQUAL( sliceDim_, 0 );
      WALBERLA_ASSERT_LESS_EQUAL( sliceDim_, 2 );
      WALBERLA_ASSERT_NOT_NULLPTR( displayProperties_ ); // call addConfigurationItem first!

      if ( ! displayProperties_->isEnabled() )
         return;

      for( auto it = field_->beginSliceXYZ( sliceInterval_ ) ; it != field_t::staticConstEnd; ++it )
      {
         Cell permutedCell = gui::FieldDisplayAdaptor<field_t>::permuteCoordAccordingToSlice( it.cell(), sliceDim_ );
         gui::CellView * cell = grid[ permutedCell.x() + nrOfGhostLayers ] [ permutedCell.y() + nrOfGhostLayers ];

         for( auto d = stencil::D2Q9::begin(); d != stencil::D2Q9::end(); ++d )
         {
            auto f = Stencil::idx[ map2Dto3D[sliceDim_][*d] ];
            if ( f == stencil::INVALID_DIR )
               cell->setPDF( *d, 0 );
            else
               cell->setPDF( *d, real_c( it.getF( f ) ) );
         }
      }


   }



} // namespace lbm
} // namespace walberla



