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
//! \file BlockViewText.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "BlockViewText.h"
#include "GuiUtil.h"

#include <QDrag>
#include <QMouseEvent>

#include "blockforest/BlockForest.h"

namespace walberla {
namespace gui {

BlockViewText::BlockViewText( QWidget * par  )
   : QListWidget( par ), blockForest_( NULL )
{

}

void BlockViewText::onBlockForestChange()
{
   this->setup( blockForest_ );
}

void BlockViewText::setup( BlockForest * blockForest)
{
   using blockforest::Block;
   blockForest_ = blockForest;

   std::vector< shared_ptr<IBlockID>  > allBlocks;
   blockForest_->getAllBlocks( allBlocks );

   this->clear();

   for( auto i = allBlocks.begin(); i != allBlocks.end(); ++i )
   {
      shared_ptr<blockforest::BlockID> bfBlockId = dynamic_pointer_cast< blockforest::BlockID > ( *i );

      if ( blockForest_->blockExistsLocally( **i ) )
      {
         Block * blockPtr = blockForest->getBlock( *bfBlockId );
         auto blockState = blockPtr->getState();
         QString stateStr;
         for( auto stateIt = blockState.begin(); stateIt != blockState.end(); ++stateIt ) {
            stateStr += QString("/") + QString::fromStdString( stateIt->getIdentifier() );
         }

         uint_t bx,by,bz;
         blockForest_->getRootBlockCoordinates( bx,by,bz, *bfBlockId );
         auto listItem = new QListWidgetItem( QString("%1, %2, %3  ").arg(bx).arg(by).arg(bz) + stateStr, this );
         listItem->setData( Qt::UserRole, VPtr<Block>::asQVariant( blockPtr) );
      }
   }

}

void BlockViewText::mousePressEvent(QMouseEvent * me)
{
   QListWidgetItem * curItem = itemAt( me->pos() );
   if ( curItem )
   {
      Block * blockPtr = VPtr<Block>::asPtr(curItem->data( Qt::UserRole ) );
      QMimeData * md = NULL;

      md = createMimeDataFromPointer( blockPtr , BLOCK_MIMETYPE );
      QDrag *drag = new QDrag(this);

      drag->setMimeData(md);
      drag->exec();

   }
   else
      return QListWidget::mousePressEvent( me );

}




} // namespace gui
} // namespace walberla


