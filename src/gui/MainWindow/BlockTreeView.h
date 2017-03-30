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
//! \file BlockTreeView.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef Q_MOC_RUN
#include "Gui.h"
#include "PropertyTree.h"
#include "domain_decomposition/IBlock.h"
#endif

#include <QTreeView>


namespace walberla {
namespace gui {


   class BlockTreeView : public QTreeView
   {
      Q_OBJECT
      public:
         BlockTreeView( const GUI & gui, QWidget * parent = 0 );
         virtual ~BlockTreeView() {}

         void setBlock( IBlock * b );
         void setPropTreeModel( shared_ptr<PropertyTree>  propTreeModel );

      public slots:
         void onDataChange();

      protected:
         virtual void dragEnterEvent(QDragEnterEvent * ev);
         virtual void dragMoveEvent(QDragMoveEvent * ev);
         virtual void dropEvent(QDropEvent *ev);

         const GUI & gui_;

         IBlock * block_;
         shared_ptr<PropertyTree>  propTreeModel_;

   };

} // namespace gui
} // namespace walberla

