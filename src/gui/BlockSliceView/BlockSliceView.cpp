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
//! \file BlockSliceView.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "BlockSliceView.h"
#include "DisplayAdaptor.h"
#include "DisplayPropertiesItem.h"
#include "GuiUtil.h"

#include "core/Abort.h"
#include "core/debug/Debug.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q19.h"

#include <QDebug>
#include <QDebug>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QFont>
#include <QGLWidget>
#include <QMimeData>
#include <QResizeEvent>
#include <QWheelEvent>


namespace walberla {
namespace gui {



BlockSliceView::BlockSliceView( StructuredBlockForest & forest,
                                const GUI & gui,
                                QStackedWidget * propertyStack,
                                ISliceChangeListener * sliceChangeListener,
                                QWidget * par )
   : QGraphicsView(par),
     sliceChangeListener_( sliceChangeListener ),
     sliceIdValid_( false ),
     propertyStack_ ( propertyStack),
     blockForest_( forest ),
     block_( NULL ),
     ghostLayerIndicator_ ( NULL )
{
   setRenderHints( QPainter::Antialiasing );
   setViewport( new QGLWidget( QGLFormat(QGL::SampleBuffers))    );
   setViewportUpdateMode( QGraphicsView::FullViewportUpdate );
   scene = new QGraphicsScene( this );
   setScene( scene );

   setDragMode( QGraphicsView::ScrollHandDrag );

   // Setup PropertyUi
   propertyWidget_ = new QWidget( this );
   propertyUi_.setupUi( propertyWidget_ );
   propertyWidget_->setEnabled( false );
   propertyWidget_->show();
   propertyStack->addWidget( propertyWidget_ );
   propertyStack->setCurrentWidget( propertyWidget_ );

   // Connect slice change signals
   connect( propertyUi_.spnSlice, SIGNAL( valueChanged(int) ),
            this,                 SLOT  ( onSliceValueChange()   ) );

   connect( propertyUi_.cmbSlice, SIGNAL( currentIndexChanged(int) ),
            this,                 SLOT  ( paintGrid()            ) );

   propertyUi_.cmbSlice->setCurrentIndex(2);

   // create display items / properties
   propertyWidget_->setVisible( false );
   gui.createAdaptors( displayItems_ );
   for( auto i = displayItems_.begin(); i != displayItems_.end(); ++i )
   {
      connect( (*i), SIGNAL(requestRedraw()), this , SLOT( redraw() ) );
      (*i)->addConfigurationItem( propertyUi_.optionsTree->invisibleRootItem() );
   }

   propertyUi_.optionsTree->installEventFilter(this);

   // For enable/disable of display adaptors
   connect( propertyUi_.optionsTree, SIGNAL(itemClicked ( QTreeWidgetItem * , int  ) ), this, SLOT( redraw() ) );
}

BlockSliceView::~BlockSliceView()
{
   if( sliceIdValid_ && sliceChangeListener_ )
      sliceChangeListener_->removeSlice( sliceId_ );

   // Delete display items
   for( auto i = displayItems_.begin(); i != displayItems_.end(); ++i )
      delete *i;

   // Delete cells
   for(int i=0; i < cells.size(); ++i )
      for(int j=0; j< cells[i].size(); ++j)
         delete cells[i][j];

}


void BlockSliceView::paintGrid()
{
   if ( displayItems_.size() < 1 )
      return;

   // Configure display items

   // Query Display Items for the size and number of ghost layers
   Vector3<uint_t> innerSize ( uint_t(0) );
   nrOfGhostLayers_ = 0;

   cell_idx_t sliceVal = 0;
   int        sliceDim = propertyUi_.cmbSlice->currentIndex ();

   for ( auto i = displayItems_.begin(); i != displayItems_.end(); ++i )
   {
      Vector3<uint_t> curSize;
      cell_idx_t curNrOfGhostLayers;
      (*i)->configure( block_, sliceDim, sliceVal, curSize, curNrOfGhostLayers );

      innerSize[0] = std::max( innerSize[0], curSize[0] );
      innerSize[1] = std::max( innerSize[1], curSize[1] );
      innerSize[2] = std::max( innerSize[2], curSize[2] );
      nrOfGhostLayers_ = std::max( nrOfGhostLayers_, curNrOfGhostLayers );
   }


   // Delete existing cells
   for( int i=0; i < cells.size(); ++i )
      for( int j=0; j< cells[i].size(); ++j )
         delete cells[i][j];

   double bps = CellView::BLOCK_PIXEL_SIZE;

   int xSize = int_c( innerSize[0] ) + 2 * nrOfGhostLayers_;
   int ySize = int_c( innerSize[1] ) + 2 * nrOfGhostLayers_;

   cells.resize( xSize );
   for( int i=0; i< xSize; ++i )
   {
      cells[i].resize(ySize);
      for(int j=0; j< ySize; ++j)
      {
           cells[i][j] = new CellView();
           cells[i][j]->setToolTip( QString("%1, %2").arg(i - nrOfGhostLayers_).arg(j - nrOfGhostLayers_) );
           scene->addItem( cells[i][j] );
           cells[i][j]->setPos(QPointF(i*bps, (ySize-j-1)*bps));
      }
   }


   // Draw red rectangle around inner region
   if ( ghostLayerIndicator_ )
      delete ghostLayerIndicator_;
   ghostLayerIndicator_ = new QGraphicsRectItem( nrOfGhostLayers_ * bps,     nrOfGhostLayers_ * bps,
                                                real_c( innerSize[0] )* bps, real_c( innerSize[1] )* bps );
   ghostLayerIndicator_->setPen( QPen( QBrush(Qt::red), 2)  );
   scene->addItem( ghostLayerIndicator_ );


   // Fit the grid into the view, but leave additional border of ADDITITONAL_SPACE
   const double ADDITIONAL_SPACE = 0.2 * bps;
   fitInView(-ADDITIONAL_SPACE,
             ySize * bps +ADDITIONAL_SPACE,
             xSize  * bps + 2*ADDITIONAL_SPACE,
             ySize  * bps + 2*ADDITIONAL_SPACE,
             Qt::KeepAspectRatio);


   updateSliceIndicatorIn3DView();

   propertyUi_.spnSlice->setMaximum( int_c( innerSize[2] - 1 + uint_c( nrOfGhostLayers_ ) ) );
   propertyUi_.spnSlice->setMinimum( - nrOfGhostLayers_ );
   propertyUi_.spnSlice->setValue ( 0 );


   redraw();
}


void BlockSliceView::updateSliceIndicatorIn3DView()
{
   if ( ! sliceChangeListener_ )
      return;

   // Update slice indicator in 3D view
   if ( sliceIdValid_ && sliceChangeListener_ ) {
      sliceChangeListener_->removeSlice( sliceId_ );
      sliceIdValid_ = false;
   }

   const cell_idx_t sliceVal = propertyUi_.spnSlice->value ();
   const int        sliceDim = propertyUi_.cmbSlice->currentIndex ();

   cell_idx_t range = cell_idx_c( propertyUi_.spnSlice->maximum() +1 - propertyUi_.spnSlice->minimum() - 2*nrOfGhostLayers_ );
   double relSliceVal = (1.0 * sliceVal + 0.5) / range;
   sliceId_ =  sliceChangeListener_->addSlice( block_, sliceDim , relSliceVal );
   sliceIdValid_ = true;

}

void BlockSliceView::onSliceValueChange()
{
   cell_idx_t sliceVal = propertyUi_.spnSlice->value ();
   int        sliceDim = propertyUi_.cmbSlice->currentIndex ();

   for (auto i = displayItems_.begin() ; i != displayItems_.end(); ++i )
   {
      Vector3<uint_t> sizeDummy;
      cell_idx_t glDummy;
      (*i)->configure( block_, sliceDim, sliceVal, sizeDummy, glDummy );
   }
   redraw();
   updateSliceIndicatorIn3DView();
}

void BlockSliceView::redraw()
{
   if( ! block_ )
      return;

   reset();

   for ( auto i = displayItems_.begin(); i != displayItems_.end(); ++i )
      (*i)->draw( cells, nrOfGhostLayers_ );
}

void BlockSliceView::reset()
{
   for(int i=0; i< cells.size(); ++i)
      for(int j=0; j< cells[i].size(); ++j)
         cells[i][j]->reset();
}

void BlockSliceView::setBlock(IBlock * block)
{
   auto blockId = dynamic_cast<blockforest::Block*>( block )->getId();
   uint_t bx,by,bz;
   blockForest_.getRootBlockCoordinates( bx,by,bz, blockId );
   nativeParentWidget()->setWindowTitle( QString("Block (%1,%2,%3)").arg( int_c(bx)).arg(int_c(by)).arg(int_c(bz) ) );

   block_ = block;
   propertyWidget_->setVisible( true );

   paintGrid();
   propertyWidget_->setEnabled( true );
}

void BlockSliceView::focusInEvent ( QFocusEvent *  )
{
   propertyStack_->setCurrentWidget( propertyWidget_ );
   if ( !block_ )
      propertyWidget_->setVisible( false );

   if( sliceIdValid_ && sliceChangeListener_ )
      sliceChangeListener_->setSliceActive( sliceId_, true );
}

void BlockSliceView::focusOutEvent ( QFocusEvent *  )
{
   if( sliceIdValid_ && sliceChangeListener_)
      sliceChangeListener_->setSliceActive( sliceId_, false );
}



bool BlockSliceView::eventFilter( QObject * obj, QEvent * ev )
{
   QTreeWidget * treeWidget = dynamic_cast<QTreeWidget*> ( obj );

   if ( treeWidget && ev->type() == QEvent::ContextMenu )
   {
      QContextMenuEvent * menuEvent = dynamic_cast<QContextMenuEvent*> (ev);
      QTreeWidgetItem * it = treeWidget->itemAt( treeWidget->viewport()->mapFromGlobal( menuEvent->globalPos() ) );


      DisplayPropertiesItem * propItem = dynamic_cast<DisplayPropertiesItem*>(it);
      if(!propItem)
         return false;

      QPoint globalPos = menuEvent->globalPos();
      propItem->showContextMenu(globalPos);
      return true;
   }

   return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Zoom  ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BlockSliceView::scaleView(qreal scaleFactor)
{
    qreal factor = matrix().scale(scaleFactor, scaleFactor).mapRect(QRectF(0, 0, 1, 1)).width();
    if (factor < 0.007 || factor > 1000)
        return;

    scale(scaleFactor, scaleFactor);
}

void BlockSliceView::wheelEvent(QWheelEvent *ev)
{
    scaleView(  std::pow((double)2, ev->delta() / 240.0));
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Drag & Drop //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void BlockSliceView::dragEnterEvent(QDragEnterEvent * ev)
{
    const QMimeData * md = ev->mimeData();
    md->hasFormat(BLOCK_MIMETYPE) ? ev->accept() : ev->ignore();
}

void BlockSliceView::dragMoveEvent(QDragMoveEvent * ev)
{
    QWidget::dragMoveEvent( ev );
}

void BlockSliceView::dropEvent(QDropEvent *ev)
{
    if( ev->source()==this )
        return;

    const QMimeData * md = ev->mimeData();
    WALBERLA_ASSERT( md->hasFormat( BLOCK_MIMETYPE ) );

    void * ptr = getPointerFromMimeData(md,BLOCK_MIMETYPE);
    WALBERLA_ASSERT_NOT_NULLPTR( ptr );

    IBlock * block = static_cast<IBlock*>(ptr);

    ev->acceptProposedAction();

    setBlock(block);
}




} // namespace gui
} // namespace walberla
