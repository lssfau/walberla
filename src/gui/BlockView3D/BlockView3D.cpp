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
//! \file BlockView3D.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "BlockView3D.h"
#include "GuiUtil.h"

#include "core/debug/Debug.h"

#include <QDebug>
#include <QWheelEvent>
#include <qglbuilder.h>
#include <qglcube.h>
#include <qglteapot.h>


namespace walberla {
namespace gui {


BlockView3D::BlockView3D(QWidget * par)
   : QGLView(par),shrinkingFactor_(0.9)
{
   // Enable anti-aliasing
   QGLFormat newFormat = this->format();
   newFormat.setSampleBuffers(true);
   setFormat(newFormat);

   setOption(QGLView::ObjectPicking, true);
   //setOption(QGLView::ShowPicking,true);

   camera()->setEye( QVector3D(25,15, 90));
}

BlockView3D::~BlockView3D()
{
   qDeleteAll(blocks_);
}

void BlockView3D::onBlockForestChange()
{
   //for( auto it = registeredBlocks_.begin(); it != registeredBlocks_.end(); ++it )
   //   deregisterObject( *it );
   //registeredBlocks_.clear();
   //qDeleteAll(blocks_);
   //this->setup( blockForest_ );
}

void BlockView3D::setup( BlockForest * blockForest)
{
   using blockforest::Block;
   blockForest_ = blockForest;

   int objectId=0;

   std::vector< shared_ptr<IBlockID>  > allBlocks;
   blockForest_->getAllBlocks( allBlocks );


   for( auto i = allBlocks.begin(); i != allBlocks.end(); ++i )
   {
      shared_ptr<blockforest::BlockID> bfBlockId = dynamic_pointer_cast< blockforest::BlockID > ( *i );

      BlockDisplayObject * displ = new BlockDisplayObject( blockForest_, bfBlockId, this );
      displ->setShrinkingFactor( shrinkingFactor_ );
      displ->setParent(this);
      blocks_.append(displ);

      if ( blockForest_->blockExistsLocally( **i ) )
      {
         registerObject(objectId,displ);
         registeredBlocks_.append( objectId );
         connect(displ, SIGNAL(hoverChanged()), this,SLOT(updateGL()));
         connect(displ, SIGNAL(hoverChanged()),this,SLOT(update()));
      }
      objectId++;
   }

   updateGL();
}



void BlockView3D::initializeGL(QGLPainter * painter)
{
   painter->setStandardEffect(QGL::LitMaterial);
}

void BlockView3D::paintGL(QGLPainter * painter)
{
   glEnable(GL_BLEND);

   const AABB & bb = blockForest_->getDomain();
   QVector3D bbMid  ( 0.5 * (bb.max(0) + bb.min(0)),
                      0.5 * (bb.max(1) + bb.min(1)),
                      0.5 * (bb.max(2) + bb.min(2)));

   double maxExtend = std::max( std::max( bbMid.x(), bbMid.y() ), bbMid.z() );

   painter->modelViewMatrix().push();
   painter->modelViewMatrix().scale( 15.0 / maxExtend );
   painter->modelViewMatrix().translate(-bbMid);

   foreach(BlockDisplayObject * blkDisp, blocks_)
      blkDisp->draw(painter);

   painter->modelViewMatrix().pop();
}


void BlockView3D::wheelEvent (QWheelEvent * we)
{
   if (we->modifiers() & Qt::ControlModifier )
   {

      shrinkingFactor_ += we->delta() / 5000.0;
      if (shrinkingFactor_ > 1 )
         shrinkingFactor_ = 1;
      else if (shrinkingFactor_ < 0.01 )
         shrinkingFactor_ = 0.01;

      foreach(BlockDisplayObject * blkDisp, blocks_)
         blkDisp->setShrinkingFactor(shrinkingFactor_);


      updateGL();
   }
   else
      QGLView::wheelEvent(we);
}

void BlockView3D::mousePressEvent(QMouseEvent * me)
{
   QObject * pickedObj = objectForPoint(me->pos());
   BlockDisplayObject * bdo = dynamic_cast<BlockDisplayObject*> (pickedObj );

   if(!bdo)
      return QGLView::mousePressEvent(me);

   QMimeData * md = NULL;

   md = createMimeDataFromPointer( blockForest_->getBlock( * bdo->blockId() ) , BLOCK_MIMETYPE );
   QDrag *drag = new QDrag(this);

   drag->setMimeData(md);
   drag->exec();
}




BlockView3D::SliceID BlockView3D::addSlice ( IBlock * block, int sliceDim, double position )
{
   for( auto i = blocks_.begin(); i != blocks_.end(); ++i )
      if( * (*i)->blockId() == block->getId() )
      {
         BlockDisplayObject * bdi = *i;
         SliceID sId =  bdi->addSlice(sliceDim, position );
         sliceToDisplayObject[sId] = (*i);
         updateGL();
         return sId;
      }

   WALBERLA_ASSERT( false ); //block not found
   updateGL();
   return -1;
}

void BlockView3D::removeSlice  ( SliceID id )
{
   sliceToDisplayObject[id]->removeSlice(id);
   updateGL();
}

void BlockView3D::setSliceActive ( SliceID id, bool active )
{
   sliceToDisplayObject[id]->setSliceActive(id, active);
   updateGL();
}





} // namespace gui
} // namespace walberla


