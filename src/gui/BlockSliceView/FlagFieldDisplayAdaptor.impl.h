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
//! \file FlagFieldDisplayAdaptor.impl.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "CellView.h"
#include "FlagFieldDisplayAdaptor.h"

#include "core/debug/Debug.h"

#include <vector>


namespace walberla {
namespace gui {

   template< typename field_t>
   QMap< QString, QColor > FlagFieldDisplayAdaptor<field_t>::createFlagNameToColorMap()
   {
      QMap< QString, QColor > result;
      result["noslip"   ] = QColor(0,0,0,150);
      result["liquid"   ] = QColor(0,130,230,150);
      result["gas"]       = QColor(0,80,130,20);
      result["interface"] = QColor(0,40,230,150);
      result["pressure0"] = QColor(210,60,0,150);
      result["pressure1"] = QColor(0,210,60,150);
      result["pressure2"] = QColor(210,60,0,150);
      result["ubb"]       = QColor(230,0,210,150);
      result["velocity0"] = QColor(130,0,230,150);
      result["velocity1"] = QColor(230,0,210,150);

      return result;
   }
   template< typename field_t>
   QMap< QString, QColor > FlagFieldDisplayAdaptor<field_t>::flagNameToColor = createFlagNameToColorMap();


   template< typename field_t>
   QVector< QColor > FlagFieldDisplayAdaptor<field_t>::createDefaultFlagColors()
   {
      QVector< QColor > result;
      result.append( QColor(175,220,0,150)  );
      result.append( QColor(220,70,70,150)  );
      result.append( QColor(220,50,150)  );
      result.append( QColor(170,170,170,150)  );
      result.append( QColor(0,14,210,150)  );
      result.append( QColor(205,102,29,150)  );
      return result;
   }
   template< typename field_t>
   QVector< QColor > FlagFieldDisplayAdaptor<field_t>::defaultFlagColors = createDefaultFlagColors();

   template< typename field_t>
   int FlagFieldDisplayAdaptor<field_t>::defaultFlagColorCounter = 0;



   template< typename field_t>
   FlagFieldDisplayAdaptor<field_t>::FlagFieldDisplayAdaptor( ConstBlockDataID scalarFieldID )
         : FieldDisplayAdaptor<field_t>( scalarFieldID ),
           displayProperties_ ( NULL ),
           flagField_( NULL ),
           lastBlock_( NULL )
   {
   }

   template< typename field_t>
   FlagFieldDisplayAdaptor<field_t>::~FlagFieldDisplayAdaptor()
   {
      if ( displayProperties_)
         delete displayProperties_;
   }

   template< typename field_t>
   void FlagFieldDisplayAdaptor<field_t>::addConfigurationItem( QTreeWidgetItem * parentItem )
   {
      if ( ! displayProperties_ )
      {
         displayProperties_ = new DisplayPropertiesItem( QStringList(), parentItem, false );
         QObject::connect( displayProperties_, SIGNAL( optionChanged() ),
                           this,               SIGNAL( requestRedraw() ) );
      }
   }


   template< typename field_t>
   void FlagFieldDisplayAdaptor<field_t>::configure( const IBlock * block, int sliceDim, cell_idx_t sliceValue,
                                       Vector3<uint_t> & innerSize, cell_idx_t & ghostLayers )
   {
      // Set name of tree item
      WALBERLA_ASSERT( displayProperties_ );
      const std::string & name = block->getBlockStorage().getBlockDataIdentifier( blockDataId_ );
      displayProperties_->setText( 0, QString::fromStdString(name) );

      FieldDisplayAdaptor<field_t>::configure( block, sliceDim, sliceValue, innerSize, ghostLayers );

      if ( block == lastBlock_ ) {
         return;
      }

      lastBlock_ = block;

      qDeleteAll( flagProperties_ );
      flagProperties_.clear();

      flagField_ = dynamic_cast<const FlagField<T> *> ( field_ );
      if ( flagField_ )
      {
         std::vector<FlagUID> allFlags;
         flagField_->getAllRegisteredFlags( allFlags );

         QStringList options;
         options << "Color" << "Text";
         for( auto i = allFlags.begin(); i != allFlags.end(); ++i )
         {
            flagProperties_.push_back( new DisplayPropertiesItem( options, displayProperties_ ) );
            QObject::connect( flagProperties_.back(), SIGNAL( optionChanged() ),
                              this,                   SIGNAL( requestRedraw() ) );

            flagProperties_.back()->enableColorSelect();
            QString flagName = QString::fromStdString( i->getIdentifier() );
            QColor flagColor;
            if( flagNameToColor.contains(flagName.toLower()) )
            {
               flagColor = flagNameToColor[flagName.toLower()];
            }
            else
            {
               flagColor = defaultFlagColors[ defaultFlagColorCounter ];
               defaultFlagColorCounter = ( defaultFlagColorCounter + 1 ) % defaultFlagColors.size();
            }
            flagProperties_.back()->setText( 0, flagName );
            flagProperties_.back()->setColor( flagColor );
         }
      }

   }

   template< typename field_t>
   void FlagFieldDisplayAdaptor<field_t>::draw( QVector<QVector<CellView*> > & grid, int nrOfGhostLayers )
   {

      using namespace stencil;
      WALBERLA_ASSERT_GREATER_EQUAL( sliceDim_, 0 );
      WALBERLA_ASSERT_LESS_EQUAL( sliceDim_, 2);
      WALBERLA_ASSERT( displayProperties_ ); // call addConfigurationItem first!

      for( auto flagIt = flagProperties_.begin(); flagIt != flagProperties_.end(); ++flagIt )
      {
         if ( !(*flagIt)->isEnabled() )
            continue;

         const FlagUID flagUID( ( *flagIt )->text( 0 ).toAscii().constData() );
         const T flag = flagField_->getFlag( flagUID );
         const int displayType = (*flagIt)->getComboBoxSelection();
         assert ( displayType >=0 && displayType <2 );

         for( auto it = field_->beginSliceXYZ( sliceInterval_ ) ; it != field_->end(); ++it )
         {
            Cell permutedCell = FieldDisplayAdaptor<field_t>::permuteCoordAccordingToSlice( it.cell(), sliceDim_ );
            CellView * cell = grid[ permutedCell.x() + nrOfGhostLayers ] [ permutedCell.y() + nrOfGhostLayers ];

            if ( !isFlagSet(*it,flag) )
               continue;

            if ( displayType == 0 )
               cell->setBrush( (*flagIt)->getSelectedColor() );
            else if ( displayType == 1 )
               cell->addText( ( *flagIt )->text( 0 ).at( 0 ) );
         }
      }
   }





} // namespace gui
} // namespace walberla



