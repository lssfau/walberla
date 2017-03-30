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
//! \file DisplayPropertiesItem.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include <QColor>
#include <QComboBox>
#include <QImage>
#include <QLinearGradient>
#include <QTreeWidgetItem>


namespace walberla {
namespace gui {


   //*******************************************************************************************************************
   /*! Represents on line in a QTreeWidget. Provides a combobox, and optional color/gradient selection.
   *
   *  A FieldTreeWidgetItem can be added to a QTreeWidget. In its first column text can be displayed,
   *  in the second column a combobox is shown with user defined entries.
   *  The context menu (right click) of the item provides color and colormap selection ( if enabled )
   *
   *  This class is used to configure the display options of a waLBerla::field
   *   - combobox is used to select display style ( color, arrow ... )
   *   - the color selectors are used to select the colormap for a color coding of the field values, or
   *     the color for a specific flag in a FlagField
   *
   *  \ingroup gui
   *
   */
   //*******************************************************************************************************************
   class DisplayPropertiesItem : public QObject, public QTreeWidgetItem
   {
   Q_OBJECT

   public:
      /// Constructs the item with the given combobox entries
      /// \param parent the parent item, must not be NULL
      DisplayPropertiesItem( const QStringList & comboBoxEntries, QTreeWidgetItem * parent, bool checkbox=true );
      virtual ~DisplayPropertiesItem() {}

      void enableColorSelect   ( bool enable = true ) { colorSelectEnabled_    = enable; }
      void enableGradientSelect( bool enable = true ) { gradientSelectEnabled_ = enable; }

      /// Returns currently selected color
      QColor getSelectedColor()   const   { return dispColor_; }

      /// Returns a color for value between 0 and 1, according to current colormap
      QColor getColorFromColormap ( double value );

      /// Sets the display color
      void setColor( QColor c );

      /// Returns index of selected combobox item (as passed in constructor)
      int    getComboBoxSelection();

      /// True if checkbox before text is selected
      bool   isEnabled() const;

      virtual void showContextMenu( QPoint globalPos );

   signals:
      /// Called whenever a option has changed (combobox, color, colormap)
      void optionChanged();

      /// Called when user selected a new color
      void colorSelected ( const QColor & newColor );

      /// Called when user selected a new colormap, use getColorFromColormap()
      /// to retrieve colormap entries
      void colormapSelected ();

      /// Called when a new combobox item was selected
      void comboBoxSelectionChanged ( int selectedOption );


   protected:

      virtual void showGradientSelect();
      virtual void showColorSelect();


      bool gradientSelectEnabled_;
      bool colorSelectEnabled_;

      void updateColormapFromGradient();

      QColor          dispColor_;
      QLinearGradient colorMapGradient_;
      QImage          colorMap_;

   };






   /****************************************************************************************************************//**
   * Helper class for inserting a combobox into a QListWidget
   ********************************************************************************************************************/
   class ComboBoxItem : public QComboBox
   {
      Q_OBJECT
      public:
          ComboBoxItem(QTreeWidgetItem* it, int col)
          {
               this->item = it;
               this->column = col;

               connect(this, SIGNAL(currentIndexChanged(int)), SLOT(changeItem(int)));
           }

          virtual ~ComboBoxItem() {}

      public slots:
          void changeItem(int index)
          {
            if(index >=0)
               item->setData(this->column, Qt::UserRole, this->itemText(index));
          }

      private:
          QTreeWidgetItem *item;
          int column;
   };


} // namespace gui
} // namespace walberla

