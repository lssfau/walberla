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
//! \file BlockViewText.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef Q_MOC_RUN
#include "Gui.h"
#endif

#include <QListWidget>


namespace walberla {
namespace blockforest {
   class BlockForest;
}
}

namespace walberla {
namespace gui {


   /********************************************************************************************************************
    * Shows a list view of the blocks contained in a BlockForest
    *
    * Analog to BlockView3D, but uglier :) however with to QGl dependency
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
   class BlockViewText : public QListWidget
   {
   Q_OBJECT

   public:
      BlockViewText( QWidget * parent = 0 );

      /// Shows the contents of a BlockField
      /// use this a initialization function of the widget
      void setup( BlockForest * blockForest);

   public slots:
      void onBlockForestChange();

   protected:
      virtual void mousePressEvent(QMouseEvent * me);

      blockforest::BlockForest * blockForest_;

   };





} // namespace gui
} // namespace walberla



