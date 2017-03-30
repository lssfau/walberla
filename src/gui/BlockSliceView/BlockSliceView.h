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
//! \file BlockSliceView.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef Q_MOC_RUN
#include "CellView.h"
#include "Gui.h"
#include "ISliceChangeListener.h"
#include "ui_BlockSliceViewProperties.h"

#include "blockforest/StructuredBlockForest.h"
#endif

#include <QGraphicsScene>
#include <QGraphicsView>
#include <QPixmap>
#include <QStackedWidget>


namespace walberla {
namespace gui {

   class BlockSliceViewProperties;


   //*******************************************************************************************************************
   /*!
   * View to display a slice of a field ( which is stored inside a block)
   *
   * \ingroup gui
   *
   * QGraphicsView based widget that can display various block data as slices.
   */
   //*******************************************************************************************************************
   class BlockSliceView : public QGraphicsView
   {
   Q_OBJECT

   public:

      /*************************************************************************************************************//**
      *  Constructor
      *  @param blockForest   Currently only StructuredBlockForest is supported here. "Structured" since the slice
      *                       views show cell based data structures. BlockForest and not BlockStorage is used because
      *                       we query for global domain information.
      *  @param gui           the GUI object is used as a factory for DisplayAdaptors which encapsulate the details of
      *                       how to draw specific block data
      *  @param propertyStack the stacked widget is used as a container for a widget where display properties can be
      *                       edited. The BlockSliceView makes sure that its own widget is on top, whenever the view
      *                       has the focus
      *  @param blockView3D   pointer to a 3D block view. The BlockSliceView uses the 3D view to display planes to
      *                       indicate the current slice
      *
      *****************************************************************************************************************/
      BlockSliceView( StructuredBlockForest & blockForest,
                      const GUI & gui,
                      QStackedWidget * propertyStack,
                      ISliceChangeListener * blockView3D,
                      QWidget * parent = 0 );

      ~BlockSliceView();


      void setBlock( IBlock * block );

   public slots:
      void redraw();

   protected:
      // Drag & Drop
      virtual void dropEvent( QDropEvent * ev );
      virtual void dragEnterEvent( QDragEnterEvent * ev );
      virtual void dragMoveEvent( QDragMoveEvent * ev );

      // Zoom
      virtual void wheelEvent( QWheelEvent *event );
      virtual void scaleView ( qreal scaleFactor );

   private slots:
      void paintGrid();
      void onSliceValueChange();

   private:
      void reset();

      virtual bool eventFilter( QObject * obj, QEvent * ev );
      virtual void focusInEvent ( QFocusEvent * event );
      virtual void focusOutEvent ( QFocusEvent * event );

      void updateSliceIndicatorIn3DView();

      //Slice Display
      ISliceChangeListener * sliceChangeListener_;
      ISliceChangeListener::SliceID sliceId_;
      bool sliceIdValid_;


      // View Properties
      QStackedWidget * propertyStack_;
      QWidget * propertyWidget_;
      Ui::BlockSliceViewProperties propertyUi_;


      QGraphicsScene * scene;
      QVector<QVector<CellView*> > cells;
      int nrOfGhostLayers_;


      StructuredBlockForest & blockForest_;
      IBlock * block_;

      QGraphicsRectItem * ghostLayerIndicator_;
      std::vector<DisplayAdaptor*> displayItems_;
   };


} // namespace gui
} // namespace walberla



