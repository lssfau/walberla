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
//! \file BlockTreeView.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "BlockTreeView.h"
#include "GuiUtil.h"

#include "core/debug/Debug.h"

#include <QAction>
#include <QDebug>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QHeaderView>
#include <QMenu>
#include <QMessageBox>


namespace walberla {
namespace gui {



BlockTreeView::BlockTreeView( const GUI & gui, QWidget * par )
   :  QTreeView(par), gui_(gui), block_(0)
{
   setAcceptDrops(true);
   setAlternatingRowColors(true);
}

void BlockTreeView::setBlock( IBlock * b )
{
   block_ = b;
   onDataChange();
}

void BlockTreeView::setPropTreeModel( shared_ptr<PropertyTree>  p )
{
   propTreeModel_ = p;
   setModel(propTreeModel_->getModel());
   header()->resizeSections(QHeaderView::ResizeToContents);
   onDataChange();
}


void BlockTreeView::onDataChange()
{
   if(!block_ || !propTreeModel_ )
      return;

   propTreeModel_->setBlock(block_);
}

void BlockTreeView::dragEnterEvent( QDragEnterEvent * ev )
{
    const QMimeData * md = ev->mimeData();
    md->hasFormat( BLOCK_MIMETYPE ) ? ev->accept() : ev->ignore();
}

void BlockTreeView::dragMoveEvent( QDragMoveEvent * ev )
{
    QWidget::dragMoveEvent(ev);
}

void BlockTreeView::dropEvent( QDropEvent * ev )
{
    if(ev->source()==this )
        return;

    const QMimeData * md = ev->mimeData();
    WALBERLA_ASSERT( md->hasFormat(BLOCK_MIMETYPE) );

    QMenu * contextMenu = new QMenu();

    uint_t chosenTreeModel=0;
    if( gui_.getPropertyTrees().size() == 0)
    {
       // report that no model was registered
       QMessageBox::warning(this,"Cannot Display Block Information",
                            "Cannot Display Block Information - No PropertyTree View registered at GuiManager");
       return;
    }
    else if ( gui_.getPropertyTrees().size() == 1)
    {
       // when there is only one registered model, use it
       chosenTreeModel=0;
    }
    else
    {
       // when more than one model is registered, display a context menu
       for( uint_t i = 0; i < gui_.getPropertyTrees().size(); ++i )
       {
          QAction * a = contextMenu->addAction(QString::fromStdString(  gui_.getPropertyTrees()[i]->getName() ) );
          a->setData( QVariant::fromValue(i) );
       }

       QAction * chosenOption = contextMenu->exec(mapToGlobal(ev->pos()));

       if(!chosenOption)
          return;

       bool ok = false;
       chosenTreeModel = chosenOption->data().toUInt(&ok);
       WALBERLA_ASSERT( ok );
    }

    void * ptr = getPointerFromMimeData(md,BLOCK_MIMETYPE);
    WALBERLA_ASSERT_NOT_NULLPTR( ptr );
    Block * block = static_cast<Block*>(ptr);
    ev->acceptProposedAction();
    setBlock(block);

    // set property tree, create a new copy for this widget
    setPropTreeModel( gui_.getPropertyTrees()[chosenTreeModel]->create() );
}



} //namespace gui
} //namespace walberla
