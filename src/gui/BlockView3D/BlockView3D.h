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
//! \file BlockView3D.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Header file for 3D Block View
//
//======================================================================================================================

#pragma once

#ifndef Q_MOC_RUN
#include "BlockDisplayObject.h"
#include "ISliceChangeListener.h"
#include "blockforest/BlockForest.h"
#endif

#include <QMainWindow>
#include <qbox3d.h>
#include <qglview.h>


namespace walberla {
namespace gui {


   /********************************************************************************************************************
    * Shows a 3D view of the blocks contained in a BlockForest
    *
    * \ingroup gui
    *
    * - the user can drag single blocks out of the view
    * - the dragged MIME data carries a blockforest::Block pointer
    * - the distance of the Blocks can be modified by Ctrl+MouseWheel
    *
    * Usage: - call setup() with a BlockForest pointer
    *        - accept drop events, for MIME data type BLOCK_MIMETYPE
    *
    *******************************************************************************************************************/
   class BlockView3D : public QGLView, public ISliceChangeListener
   {
   Q_OBJECT

   public:
      BlockView3D(QWidget * parent =0);
      virtual ~BlockView3D();

      /// Shows the contents of a BlockField
      /// use this as initialization function of the widget
      void setup( BlockForest * blockForest);

      //** Slice Indicators        *************************************************************************************
      /*! \name Slice Indicators  */
      //@{
      typedef BlockDisplayObject::SliceID SliceID;
      virtual SliceID addSlice       ( IBlock * block, int sliceDim, double position );
      virtual void    removeSlice    ( SliceID id );
      virtual void    setSliceActive ( SliceID id, bool active = true );
      //@}
      //****************************************************************************************************************

   public slots:
      void onBlockForestChange();

   protected:
      void paintGL(QGLPainter * painter);
      void initializeGL(QGLPainter * painter);

      // Event handling
      virtual void wheelEvent     (QWheelEvent * we);
      virtual void mousePressEvent(QMouseEvent * me);

   private:
      /// The bounding box of the blocks is scaled down to generated
      /// free space between blocks, which makes picking easier
      /// the factor can be changed dynamically by using Ctrl+MouseWheel
      double shrinkingFactor_;

      BlockForest * blockForest_;

      QList<BlockDisplayObject*> blocks_;
      QList<int> registeredBlocks_;
      QMap<SliceID, BlockDisplayObject*> sliceToDisplayObject;
   };

} // namespace gui
} // namespace walberla


