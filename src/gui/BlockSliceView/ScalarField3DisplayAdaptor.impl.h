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
//! \file ScalarField3DisplayAdaptor.impl.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "CellView.h"
#include "ScalarField3DisplayAdaptor.h"

#include "core/debug/Debug.h"

#include "field/GhostLayerField.h"

#include <cassert>
#include <cmath>
#include <type_traits>

namespace walberla {
    namespace gui {


        template< typename field_t>
        ScalarField3DisplayAdaptor<field_t>::ScalarField3DisplayAdaptor( ConstBlockDataID scalarFieldID )
                : FieldDisplayAdaptor< field_t >( scalarFieldID ),
                  displayProperties_ ( NULL )
        {
        }

        template< typename field_t>
        ScalarField3DisplayAdaptor<field_t>::~ScalarField3DisplayAdaptor()
        {
           if ( displayProperties_)
              delete displayProperties_;
        }

        template< typename field_t>
        void ScalarField3DisplayAdaptor<field_t>::configure( const IBlock * block, int sliceDim, cell_idx_t sliceValue,
                                                            Vector3<uint_t> & innerSize, cell_idx_t & ghostLayers )
        {
           // Set name of tree item
           WALBERLA_ASSERT( displayProperties_ );
           const std::string & name = block->getBlockStorage().getBlockDataIdentifier( blockDataId_ );
           displayProperties_->setText(0, QString::fromStdString(name) );

           base_t::configure( block, sliceDim, sliceValue, innerSize, ghostLayers );
        }


        template< typename field_t>
        void ScalarField3DisplayAdaptor<field_t>::addConfigurationItem( QTreeWidgetItem * parentItem )
        {
           if ( ! displayProperties_ )
           {
              QStringList options;
              options << "Numeric" << "Color Map" << "Arrow";

              displayProperties_ = new DisplayPropertiesItem( options, parentItem );
              displayProperties_->enableGradientSelect();

              QObject::connect( displayProperties_, SIGNAL( optionChanged() ),
                                this,               SIGNAL( requestRedraw() ) );
           }
        }

        template< typename field_t>
        void ScalarField3DisplayAdaptor<field_t>::draw( QVector<QVector<CellView*> > & grid, int nrOfGhostLayers )
        {
           WALBERLA_ASSERT_GREATER_EQUAL( sliceDim_, 0 );
           WALBERLA_ASSERT_LESS_EQUAL( sliceDim_, 2);
           WALBERLA_ASSERT( displayProperties_ ); // call addConfigurationItem first!
           if ( ! displayProperties_->isEnabled() )
              return;

           // displayType is option index in the QStringList passed to the constructor of displayProperties
           int displayType = displayProperties_->getComboBoxSelection();
           assert ( displayType >=0 && displayType <3 );

           length_t min = std::numeric_limits<length_t>::max();
           length_t max = std::numeric_limits<length_t>::min();
           if ( displayType == 1 || displayType == 2) // for colormap and numeric, min & max are needed
              for( auto it = field_->beginSliceXYZ( sliceInterval_ ) ; it != field_->end(); ++it )
              {
                 length_t val = std::sqrt( it[0]*it[0] + it[1]*it[1] + it[2]*it[2] );
                 if (val < min ) min = val;
                 if (val > max ) max = val;
              }


           for( auto it = field_->beginSliceXYZ( sliceInterval_ ) ; it != field_->end(); ++it )
           {
              Cell c = base_t::permuteCoordAccordingToSlice( it.cell(), sliceDim_ );
              CellView * cell = grid[ c.x() + nrOfGhostLayers ] [ c.y() + nrOfGhostLayers ];

              if ( displayType == 0 )
                 drawVectorFieldNumeric( cell, it );
              else if( displayType == 1 )
                 drawVectorFieldColormap( cell, it, min ,max );
              else if( displayType == 2 )
                 drawVectorFieldArrow( cell, it, max );
           }

        }

        template< typename field_t>
        void ScalarField3DisplayAdaptor<field_t>::drawVectorFieldNumeric( CellView * cell, const typename field_t::const_iterator & it )
        {
           if( std::is_same<T,float>::value || std::is_same<T,double>::value )
              cell->setText( QString("%1\n%2\n%3").arg( real_c( it[0]), 0,'g',6)
                                                  .arg( real_c( it[1]), 0,'g',6)
                                                  .arg( real_c( it[2]), 0,'g',6) );
           else
              cell->setText( QString("%1\n%2\n%3").arg( it[0])
                                                  .arg( it[1])
                                                  .arg( it[2]) );
        }


        template< typename field_t>
        void ScalarField3DisplayAdaptor<field_t>::drawVectorFieldColormap( CellView * cell, const typename field_t::const_iterator & it,
                                                                          typename ScalarField3DisplayAdaptor<field_t>::length_t min,
                                                                          typename ScalarField3DisplayAdaptor<field_t>::length_t max )
        {
           real_t length = std::sqrt( it[0]*it[0] + it[1]*it[1] + it[2]*it[2] );
           real_t normVal = real_c( length - min ) / real_c ( max-min );
           if ( math::finite( normVal ) )
              cell->setBrush( displayProperties_->getColorFromColormap( normVal ) );
        }

        template< typename field_t>
        void ScalarField3DisplayAdaptor<field_t>::drawVectorFieldArrow(CellView * cell, const typename field_t::const_iterator & it,
                                                                      typename ScalarField3DisplayAdaptor<field_t>::length_t max )
        {
           WALBERLA_ASSERT( sliceDim_ <= 2 && sliceDim_ >=0 );

           uint_t firstDim, secondDim;
           if ( sliceDim_ == 0 )      firstDim = 1, secondDim = 2;
           else if ( sliceDim_ == 1 ) firstDim = 0, secondDim = 2;
           else                       firstDim = 0, secondDim = 1;

           cell->setArrow( real_c( it[ firstDim  ]  ) / max,
                           real_c( it[ secondDim ] ) / max );
        }



    } // namespace gui
} // namespace walberla




