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
//! \file DisplayPropertiesItem.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "DisplayPropertiesItem.h"
#include "qtgradientdialog.h"

#include "core/debug/Debug.h"

#include <QColorDialog>
#include <QDebug>
#include <QInputDialog>
#include <QMenu>
#include <QPainter>
#include <cassert>


namespace walberla {
namespace gui {



//======================================================================================================================
//
//  FieldTreeWidget Implementation
//
//======================================================================================================================


DisplayPropertiesItem::DisplayPropertiesItem( const QStringList & options, QTreeWidgetItem * parItem, bool checkbox)
   : QTreeWidgetItem(parItem),
     gradientSelectEnabled_(false),
     colorSelectEnabled_   (false),
     colorMapGradient_(0,0,100,0),
     colorMap_( 100, 1, QImage::Format_ARGB32 )
{
   WALBERLA_ASSERT( parItem );
   // the parent has to be already added to a treeWidget
   // because the ComboBox has to be set as itemWidget which is a method of the treeWidget
   WALBERLA_ASSERT( parItem->treeWidget() );

   // By default the item is not selected
   if ( checkbox )
      setCheckState( 0, Qt::Unchecked );

   // Specify a default colormap
   colorMapGradient_.setColorAt( 1.0, QColor(0,    0, 255) );
   colorMapGradient_.setColorAt( 0.5, QColor(0,  255,   0) );
   colorMapGradient_.setColorAt(   0, QColor(255,255, 255) );
   updateColormapFromGradient();

   dispColor_ = QColor::fromHsvF( 1.0* qrand() / RAND_MAX, 1.0, 0.9 );
   dispColor_.setAlphaF( 0.6 );

   if ( ! options.empty() )
   {
      // fill combobox
      ComboBoxItem * cmb = new ComboBoxItem( this, 1 );
      for( auto i = options.begin(); i != options.end(); ++i )
         cmb->addItem( *i  );

      connect( cmb, SIGNAL(currentIndexChanged(int)), this, SIGNAL( comboBoxSelectionChanged(int) ) );
      connect( cmb, SIGNAL(currentIndexChanged(int)), this, SIGNAL( optionChanged() ) );


      treeWidget()->setItemWidget( this, 1, cmb );
   }
}



void DisplayPropertiesItem::updateColormapFromGradient()
{
   QPainter painter( &colorMap_ );
   painter.eraseRect( 0, 0, colorMap_.width(), colorMap_.height() );
   painter.fillRect ( 0, 0, colorMap_.width(), colorMap_.height(), QBrush( colorMapGradient_ ) );
}


void DisplayPropertiesItem::showContextMenu( QPoint globalPos )
{
   QMenu   * contextMenu = new QMenu();
   QAction * actSelectColor = NULL;
   QAction * actSelectGradient = NULL;

   if( colorSelectEnabled_ )
      actSelectColor = contextMenu->addAction(tr("Select Color"));
   if( gradientSelectEnabled_ )
      actSelectGradient = contextMenu->addAction(tr("Select ColorMap"));

   QAction * chosenOption = contextMenu->exec(globalPos);

   if( chosenOption == actSelectColor )
      showColorSelect();
   else if ( chosenOption == actSelectGradient )
      showGradientSelect();
}

QColor DisplayPropertiesItem::getColorFromColormap(double val)
{
   WALBERLA_ASSERT_GREATER_EQUAL( val, 0 );
   WALBERLA_ASSERT_LESS_EQUAL( val, 1 );
   QColor result( colorMap_.pixel( (int)(val * 99), 0 ) );
   return result;
}

bool DisplayPropertiesItem::isEnabled() const
{
   return (checkState(0) == Qt::Checked );
}



int DisplayPropertiesItem::getComboBoxSelection()
{
   return dynamic_cast<ComboBoxItem* > ( treeWidget()->itemWidget(this,1) )->currentIndex();
}


void DisplayPropertiesItem::showGradientSelect()
{
   bool ok = false;

   QGradient newColorMap = QtGradientDialog::getGradient(&ok,colorMapGradient_,0,"Select Color Map" );
   if( ok )
   {
      QLinearGradient lastGrad;
      colorMapGradient_ = QLinearGradient(0,0,100,0);
      for(int i=0; i< newColorMap.stops().size(); ++i)
         colorMapGradient_.setColorAt(newColorMap.stops()[i].first, newColorMap.stops()[i].second );

      updateColormapFromGradient();

      emit colormapSelected();
      emit optionChanged();
   }
}

void DisplayPropertiesItem::showColorSelect()
{
   QColor selectedColor = QColorDialog::getColor( dispColor_, NULL, "Select Color", QColorDialog::ShowAlphaChannel );
   if ( selectedColor.isValid() ) {
      setColor( selectedColor );
      emit colorSelected( dispColor_ );
      emit optionChanged();
   }
}

void DisplayPropertiesItem::setColor( QColor c )
{
   dispColor_ = c;
}

} // namespace gui
} // namespace walberla

