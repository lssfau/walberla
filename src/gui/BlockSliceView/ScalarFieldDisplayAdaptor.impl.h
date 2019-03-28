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
//! \file ScalarFieldDisplayAdaptor.impl.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "CellView.h"
#include "ScalarFieldDisplayAdaptor.h"

#include "core/debug/Debug.h"

#include "field/GhostLayerField.h"

#include <cassert>
#include <limits>
#include <type_traits>


namespace walberla {
namespace gui {


template< typename field_t>
ScalarFieldDisplayAdaptor<field_t>::ScalarFieldDisplayAdaptor( ConstBlockDataID scalarFieldID )
      : FieldDisplayAdaptor<field_t>( scalarFieldID ),
        f(0),
        displayProperties_ ( NULL )
{
}

template< typename field_t>
ScalarFieldDisplayAdaptor<field_t>::~ScalarFieldDisplayAdaptor()
{
   if ( displayProperties_)
      delete displayProperties_;
}

template< typename field_t>
void ScalarFieldDisplayAdaptor<field_t>::configure( const IBlock * block, int sliceDim, cell_idx_t sliceValue,
                                    Vector3<uint_t> & innerSize, cell_idx_t & ghostLayers )
{
   // Set name of tree item
   WALBERLA_ASSERT( displayProperties_ );
   const std::string & name = block->getBlockStorage().getBlockDataIdentifier( blockDataId_ );
   displayProperties_->setText(0, QString::fromStdString(name) );

   FieldDisplayAdaptor<field_t>::configure( block, sliceDim, sliceValue, innerSize, ghostLayers );
}


template< typename field_t>
void ScalarFieldDisplayAdaptor<field_t>::addConfigurationItem( QTreeWidgetItem * parentItem )
{
   if ( ! displayProperties_ )
   {
      QStringList options;

      if( field_t::F_SIZE == 1) {
         options << "Numeric" << "Color Map";
      }
      else
      {
         for( uint_t i=0; i < field_t::F_SIZE; ++i ) {
            options << QString("%1 Numeric").arg(i);
            options << QString("%1 Color Map").arg(i);
         }
      }

      displayProperties_ = new DisplayPropertiesItem( options, parentItem );
      displayProperties_->enableGradientSelect();

      QObject::connect( displayProperties_, SIGNAL( optionChanged() ),
                        this,               SIGNAL( requestRedraw() ) );
   }
}

template< typename field_t>
void ScalarFieldDisplayAdaptor<field_t>::draw( QVector<QVector<CellView*> > & grid, int nrOfGhostLayers )
{
   WALBERLA_ASSERT_GREATER_EQUAL( sliceDim_, 0 );
   WALBERLA_ASSERT_LESS_EQUAL( sliceDim_, 2);
   WALBERLA_ASSERT( displayProperties_ ); // call addConfigurationItem first!
   if ( ! displayProperties_->isEnabled() )
      return;

   // displayType is option index in the QStringList passed to the constructor of displayProperties
   int comboSelection = displayProperties_->getComboBoxSelection();
   int displayType = comboSelection % 2;
   f = comboSelection / 2;
   assert ( displayType >=0 && displayType <2 );

   T min = std::numeric_limits<T>::max();
   T max = std::numeric_limits<T>::min();
   if ( displayType == 1 )
      for( auto it = field_->beginSliceXYZ( sliceInterval_ ) ; it != field_->end(); ++it )
      {
         if ( *it >= std::numeric_limits<T>::max() ) continue;

         if (*it < min ) min = *it;
         if (*it > max ) max = *it;
      }


   for( auto it = field_->beginSliceXYZ( sliceInterval_ ) ; it != field_->end(); ++it )
   {
      Cell c = FieldDisplayAdaptor<field_t>::permuteCoordAccordingToSlice( it.cell(), sliceDim_ );
      CellView * cell = grid[ c.x() + nrOfGhostLayers ] [ c.y() + nrOfGhostLayers ];

      if ( displayType == 0 )
         drawScalarFieldNumeric( cell, it );
      else if( displayType == 1 )
         drawScalarFieldColormap( cell, it , min ,max );
   }

}

template< typename field_t>
void ScalarFieldDisplayAdaptor<field_t>::drawScalarFieldNumeric( CellView * cell,
                                                                 const typename field_t::const_iterator & it)
{
   if( std::is_same<T,float>::value || std::is_same<T,double>::value )
      cell->setText(QString("%1").arg(real_c(it.getF(f)), 0, 'g', 6));
   else if ( it.getF(f) < std::numeric_limits<T>::max() )
      cell->setText( QString("%1").arg(it.getF(f))  );
   else
      cell->setText("");
}


template< typename field_t>
void ScalarFieldDisplayAdaptor<field_t>::drawScalarFieldColormap( CellView * cell,
                                                                  const typename field_t::const_iterator & it,
                                                                  T min, T max )
{
   real_t normVal = 0;
   if ( fabs( real_c(max) - real_c(min) ) > 1e-7 )
       normVal = real_c(it.getF(f) - min) / real_c(max - min);

   if ( it.getF(f) < std::numeric_limits<T>::max() )
      cell->setBrush( displayProperties_->getColorFromColormap( normVal ) );
   else
      cell->setBrush( QBrush(Qt::red) );
}


} // namespace gui
} // namespace walberla


