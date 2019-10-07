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
//! \file CellView.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "CellView.h"
#include "stencil/D2Q9.h"

#include <QPainter>
#include <QPen>


namespace walberla {
namespace gui {

   static const double Pi    = 3.14159265358979323846264338327950288419717;
   static const double TwoPi = 2.0 * Pi;

   const double CellView::BLOCK_PIXEL_SIZE = 100;

   CellView::CellView( QGraphicsItem * par )
      : QGraphicsRectItem(0,0,BLOCK_PIXEL_SIZE,BLOCK_PIXEL_SIZE,par),
        pdfsEnabled_(false),
        arrowEnabled_(false),
        arrowX_( 0.0 ),
        arrowY_( 0.0 ),
        vectorArrowSize(20),
        pdfArrowSize(3)
   {
      setPen( QPen(Qt::black) );

      for(int i=0; i<9; ++i)
         pdfs_[i] = 0;
   }

   void CellView::setPDF( stencil::Direction dir, const double & val)
   {
      pdfsEnabled_ = true;
      pdfs_ [ stencil::D2Q9::idx[dir] ]  = val;
      update();
   }



   void CellView::paint(QPainter *p, const QStyleOptionGraphicsItem *option, QWidget *widget)
   {
      QGraphicsRectItem::paint( p, option, widget );

      p->save();
      p->setRenderHint( QPainter::TextAntialiasing );

      p->setPen(Qt::black);
      p->setBrush(Qt::black);
      QFont f = p->font();
      f.setPixelSize(12);
      p->setFont(f);
      QFontMetrics fontMetrics = p->fontMetrics();

      QRectF textRect( QPointF(0,0), fontMetrics.size(0,text_) );
      textRect.moveCenter( QPointF(50,50) );
      p->drawText(textRect,text_);


      if(arrowEnabled_)
         paintVector(p);

      if(pdfsEnabled_)
         paintPDFs(p);

      paintTexts( p );

      p->restore();
   }


   void CellView::paintPDFs (QPainter * p)
   {
      const double halfCellSize = BLOCK_PIXEL_SIZE * 0.5;

      // Paint arrows
      p->save();
      p->setPen(QPen(QColor(145,192,255), 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
      p->setBrush(QColor(145,192,255));
      for ( auto d = stencil::D2Q9::begin(); d != stencil::D2Q9::end(); ++d )
      {
         const double dx = d.cx() * 0.9;
         const double dy = d.cy() * 0.9;

         QPointF cellMidPoint = QPointF(halfCellSize, halfCellSize);
         QPointF dst  = cellMidPoint + QPointF( dx * halfCellSize, - dy  * halfCellSize);
         paintArrow( p, cellMidPoint, dst, pdfArrowSize );

      }
      p->restore();


      // Paint text ( separate loop, to paint text over arrows )
      p->save();
      p->setRenderHint( QPainter::TextAntialiasing );

      for ( auto d = stencil::D2Q9::begin(); d != stencil::D2Q9::end(); ++d )
      {
         const double val = pdfs_[ d.toIdx() ];
         const double dx = d.cx() * 0.9;
         const double dy = d.cy() * 0.9;

         QPointF cellMidPoint = QPointF(halfCellSize, halfCellSize);
         QPointF dst  =  QPointF( dx * halfCellSize, - dy  * halfCellSize);

         QFont f = p->font();
         f.setPixelSize(5);
         p->setFont(f);
         QFontMetrics fontMetrics = p->fontMetrics();

         QPointF textPos =  cellMidPoint + 0.8 * dst;
         QString s =  QString("%1").arg(val);
         textPos.setX(textPos.x() - 0.5 * fontMetrics.width(s));
         p->drawText(textPos,s);
      }
      p->restore();

   }


   void CellView::paintArrow(QPainter * p, const QPointF & src, const QPointF & dst, double arrowSize)
   {
      QLineF line(src, dst);
      if (qFuzzyCompare(line.length(), qreal(0.)))
          return;

      // Draw the line itself
      p->drawLine(line);

      // Draw the arrows
      double angle = ::acos(line.dx() / line.length());
      if (line.dy() >= 0)
          angle = TwoPi - angle;

      QPointF destArrowP1 = dst + QPointF(sin(angle - Pi / 3) * arrowSize,
                                          cos(angle - Pi / 3) * arrowSize);
      QPointF destArrowP2 = dst + QPointF(sin(angle - Pi + Pi / 3) * arrowSize,
                                          cos(angle - Pi + Pi / 3) * arrowSize);


      p->drawPolygon(QPolygonF() << line.p2() << destArrowP1 << destArrowP2);
   }


   void CellView::paintVector(QPainter * p)
   {
      p->setRenderHint( QPainter::Antialiasing );
      const double halfCellSize = BLOCK_PIXEL_SIZE * 0.5;
      QPointF cellMidPoint = QPointF(halfCellSize, halfCellSize);

      QPointF sourcePoint = cellMidPoint + QPointF(-arrowX_* halfCellSize,  arrowY_ * halfCellSize);
      QPointF destPoint   = cellMidPoint + QPointF( arrowX_ *halfCellSize, -arrowY_ * halfCellSize);

      paintArrow(p,sourcePoint,destPoint, vectorArrowSize);
   }

   void CellView::paintTexts(QPainter * p)
   {
      const int TEXT_DIST = 8;

      QFont f = p->font();
      f.setPixelSize(7);
      p->setFont(f);
      p->setRenderHint( QPainter::TextAntialiasing );

      QFontMetrics fontMetrics = p->fontMetrics();

      QPoint curPos (5,fontMetrics.height()+5);

      foreach(const QString & s, texts)
      {
         p->drawText(curPos,s);
         curPos.setX( curPos.x()  + fontMetrics.size(0,s).width() + TEXT_DIST);
      }
   }

   void CellView::reset()
   {
      texts.clear();
      //icons.clear();

      pdfsEnabled_ = false;

      text_="";
      arrowEnabled_ = false;
      setBrush(QBrush(QColor(255,255,255)));
      update();
   }


   void CellView::setArrow(double xVal, double yVal)
   {
      arrowEnabled_ = true;
      arrowX_ = xVal;
      arrowY_ = yVal;
   }



} // namespace gui
} // namespace walberla


