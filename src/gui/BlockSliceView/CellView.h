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
//! \file CellView.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "stencil/Directions.h"

#include <QGraphicsRectItem>


class QPainter;
class QStyleOptionGraphicsItem;
class QPointF;

namespace walberla {
namespace gui {




   //*******************************************************************************************************************
   /*! Class for drawing a single cell of a grid slice
   *
   * \ingroup gui
   *
   */
   //*******************************************************************************************************************
   class CellView : public QGraphicsRectItem
   {
   public:
      static const double BLOCK_PIXEL_SIZE;

      CellView(QGraphicsItem * par =0);

      void setText(const QString & text) { text_ = text; update(); }

      /// Enables arrow display and sets x and y value
      /// xVal and yVal should not form a vector with length > 1
      void setArrow(double xVal, double yVal);

      void reset();

      void addText(const QString & s) { texts.push_back(s); }

      void setPDF( stencil::Direction dir, const double & val);

   protected:
      void paintVector( QPainter * p );
      void paintTexts( QPainter * p );

      void paintPDFs( QPainter * p );

      void paintArrow( QPainter * p, const QPointF & src, const QPointF & dst, double arrowSize );

      virtual void paint( QPainter *p, const QStyleOptionGraphicsItem *option, QWidget *widget );

      QVector<QString> texts;

      bool pdfsEnabled_;
      double pdfs_[9];

      bool arrowEnabled_;
      double arrowX_;
      double arrowY_;

      double vectorArrowSize;
      double pdfArrowSize;
      QString text_;
   };




} // namespace gui
} // namespace walberla


