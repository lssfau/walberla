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
//! \file BlockDisplayObject.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "BlockDisplayObject.h"
#include "qglbuilder.h"
#include "qglcube.h"
#include "qglview.h"

#include "blockforest/BlockForest.h"

#include "core/math/Vector3.h"


namespace walberla {
namespace gui {

   int BlockDisplayObject::colorId = 0;
   int BlockDisplayObject::nextSliceId = 0;


   BlockDisplayObject::BlockDisplayObject( BlockForest * blockForest,
                                           const shared_ptr<blockforest::BlockID> & id ,
                                           QObject * par )
         : QObject(par),
           blockForest_(blockForest),
           blockId_(id),
           blockMaterial_(NULL),
           hoverMaterial_(NULL),
           hovering_(false)
   {
      QGLBuilder builder;
      builder << QGL::Faceted;
      builder << QGLCube();

      objectId_ = colorId++;

      mesh_ = builder.finalizedSceneNode();

      blockMaterial_ = new QGLMaterial(this);
      blockMaterial_->setAmbientColor(QColor(0, 174, 255,210));
      blockMaterial_->setDiffuseColor(QColor(0, 174, 255,210));
      blockMaterial_->setSpecularColor(QColor(255, 255, 255,210));
      blockMaterial_->setShininess(128);

      inactiveMaterial_ = new QGLMaterial(this);
      inactiveMaterial_->setAmbientColor(QColor(240, 100, 100, 100));
      inactiveMaterial_->setDiffuseColor(QColor(240, 100, 100, 100));
      inactiveMaterial_->setSpecularColor(QColor(255, 255, 255, 120));
      inactiveMaterial_->setShininess(200);

      hoverMaterial_ = new QGLMaterial(this);
      hoverMaterial_->setAmbientColor(QColor(200, 200, 200));
      hoverMaterial_->setDiffuseColor(QColor(255, 255, 255));
      hoverMaterial_->setSpecularColor(QColor(255, 255, 255));
      hoverMaterial_->setShininess(128);

      activeSliceMaterial_ = new QGLMaterial(this);
      activeSliceMaterial_->setAmbientColor(QColor(255, 0, 0,240));
      activeSliceMaterial_->setDiffuseColor(QColor(255, 0, 0));
      activeSliceMaterial_->setSpecularColor(QColor(255, 0, 0));
      activeSliceMaterial_->setShininess(128);

      inactiveSliceMaterial_ = new QGLMaterial(this);
      inactiveSliceMaterial_->setAmbientColor(QColor(1, 255, 43, 200 ));
      inactiveSliceMaterial_->setDiffuseColor(QColor(255, 255, 255));
      inactiveSliceMaterial_->setSpecularColor(QColor(255, 255, 255));
      inactiveSliceMaterial_->setShininess(128);
   }

   BlockDisplayObject::~BlockDisplayObject()
   {
      delete mesh_;
   }


   BlockDisplayObject::SliceID BlockDisplayObject::addSlice( int sliceDim, double position )
   {
      WALBERLA_ASSERT_GREATER_EQUAL( sliceDim, 0 );
      WALBERLA_ASSERT_LESS_EQUAL( sliceDim, 2 );

      int sliceID = nextSliceId++;
      slices[ sliceID ].active= true;
      slices[ sliceID ].sliceDimension = sliceDim;
      slices[ sliceID ].position = position ;

      return sliceID;
   }



   void BlockDisplayObject::setShrinkingFactor(qreal factor)
   {
      AABB bb;
      blockForest_->getAABB( bb, *blockId_ );

      const qreal bbSize []= { bb.max(0) - bb.min(0),
                               bb.max(1) - bb.min(1),
                               bb.max(2) - bb.min(2) };

      modelViewMat_ = QMatrix4x4(); // reset to identity
      modelViewMat_.translate( bb.min(0) + factor/2 * bbSize[0],
                               bb.min(1) + factor/2 * bbSize[1],
                               bb.min(2) + factor/2 * bbSize[2]);

      modelViewMat_.scale( factor );
      modelViewMat_.scale( bbSize[0],bbSize[1],bbSize[2] );
   }



   void BlockDisplayObject::draw( QGLPainter *painter )
   {
      bool localBlock = blockForest_->blockExistsLocally( * blockId_ );

      QGLMaterial *curMat = 0;
      if (hovering_ )
         curMat = hoverMaterial_;
      else if ( !hovering_ && localBlock )
         curMat = blockMaterial_;
      else if ( !hovering_ && ! localBlock )
         curMat = inactiveMaterial_;

      WALBERLA_ASSERT( curMat );


      painter->setStandardEffect( QGL::LitMaterial );

      // Mark the object for object picking purposes.
      int prevObjectId = painter->objectPickId();
      painter->setObjectPickId(objectId_);

      painter->modelViewMatrix().push();
      (painter->modelViewMatrix() ) *= modelViewMat_;

      // Draw slices
      for( auto i = slices.begin(); i != slices.end(); ++i )
      {
         painter->modelViewMatrix().push();

         Slice & slice = i.value();
         QGLMaterial * sliceMat = slice.active ? activeSliceMaterial_ : inactiveSliceMaterial_;
         painter->setColor( sliceMat->diffuseColor() );
         painter->setFaceMaterial( QGL::AllFaces, sliceMat );

         Vector3<double> translateVec (0, 0, 0 );
         translateVec[ uint_c(slice.sliceDimension) ] = slice.position - 0.5;
         painter->modelViewMatrix().translate( translateVec[0], translateVec[1], translateVec[2] );

         Vector3<double> scaleVec ( 1.35, 1.35, 1.35 );
         scaleVec[ uint_c(slice.sliceDimension) ] = 0.01;
         painter->modelViewMatrix().scale( scaleVec[0], scaleVec[1], scaleVec[2] );
         mesh_->draw( painter );
         painter->modelViewMatrix().pop();
      }

      painter->setColor(curMat->diffuseColor());
      painter->setFaceMaterial(QGL::AllFaces, curMat);

      // Draw the block
      mesh_->draw(painter);

      // Revert to the previous object identifier.
      painter->setObjectPickId(prevObjectId);

      painter->modelViewMatrix().pop();
   }


   bool BlockDisplayObject::event( QEvent *e )
   {
      // Convert the raw event into a signal representing the user's action.
      if (e->type() == QEvent::Enter)
      {
         hovering_ = true;
         emit hoverChanged();
      }
      if (e->type() == QEvent::Leave )
      {
         hovering_ = false;
         emit hoverChanged();
      }

      return QObject::event(e);
   }


} // namespace gui
} // namespace walberla


