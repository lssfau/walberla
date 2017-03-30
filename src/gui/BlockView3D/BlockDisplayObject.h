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
//! \file BlockDisplayObject.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef Q_MOC_RUN
#include "blockforest/BlockForest.h"
#include "ISliceChangeListener.h"
#endif

#include "qglabstractscene.h"
#include "qglpainter.h"
#include "qmatrix4x4stack.h"

#include <QEvent>
#include <QObject>


class QGLView;
class QGLSceneNode;

namespace walberla {
namespace gui {

    //******************************************************************************************************************
    /*! BlockDisplayObject renders one Block as a Cube
    *
    * \ingroup gui
    *
    * Features: - Picking is enabled - connect to signal clicked()
    *           - three different materials: hoverMaterial is used on mouseOver,
    *                                        inactiveMaterial is used for not allocated blocks
    */
   //*******************************************************************************************************************
   class BlockDisplayObject: public QObject
   {

   Q_OBJECT

   public:

      explicit BlockDisplayObject( BlockForest * blockForest, const shared_ptr<blockforest::BlockID> & id  , QObject * parent = 0 );
      virtual ~BlockDisplayObject();

      void draw( QGLPainter *painter );

      shared_ptr<blockforest::BlockID> blockId() { return blockId_; }

      void setShrinkingFactor(qreal factor);


      //** Slice Indicators        *************************************************************************************
      /*! \name Slice Indicators  */
      //@{
      typedef ISliceChangeListener::SliceID SliceID;
      SliceID addSlice       ( int sliceDim, double position );
      void    removeSlice    ( SliceID id )                     { slices.remove(id); }
      void    setSliceActive ( SliceID id, bool active = true ) { slices[id].active = active; }
      //@}
      //****************************************************************************************************************

   signals:
      void hoverChanged();

   protected:
      bool event( QEvent *e );

   private:
      BlockForest * blockForest_;
      shared_ptr<blockforest::BlockID> blockId_;

      QGLSceneNode * mesh_;
      QGLMaterial  * blockMaterial_;
      QGLMaterial  * inactiveMaterial_;
      QGLMaterial  * hoverMaterial_;
      QGLMaterial  * activeSliceMaterial_;
      QGLMaterial  * inactiveSliceMaterial_;


      int  objectId_;

      bool hovering_;

      QMatrix4x4 modelViewMat_;

      /// Static counter, incremented for each object, is used as objectId_
      /// -> every cube gets its own objectId_ for picking
      static int colorId;


      // Slices
      struct Slice
      {
         bool   active;
         int    sliceDimension; // 0,1 or 2
         double position;
      };
      static int nextSliceId;
      QMap<SliceID, Slice > slices;

   };

} // namespace gui
} // namespace walberla


