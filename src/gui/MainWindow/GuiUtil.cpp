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
//! \file GuiUtil.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "GuiUtil.h"
#include "qtgradientdialog.h"

#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include <QLinearGradient>
#include <QMimeData>
#include <QPainter>
#include <QTreeWidgetItem>
#include <cassert>


namespace walberla {
namespace gui {

//======================================================================================================================
//
//  MIME DATA PACKING/UNPACKING USED FOR DRAG & DROP
//
//======================================================================================================================



QMimeData * createMimeDataFromPointer(void * p, const QString & role)
{
    QMimeData * data = new QMimeData();

    QByteArray d;
    QDataStream s( &d, QIODevice::Unbuffered | QIODevice::ReadWrite );


    int bytesWritten = s.writeRawData( ( const char * ) & p, int_c( sizeof( void * ) ) );
    WALBERLA_ASSERT_EQUAL( bytesWritten, sizeof( void*) );
    WALBERLA_UNUSED(bytesWritten);

    data->setData( role, d );

    return data;
}

void * getPointerFromMimeData(const QMimeData *data, const QString& type)
{
    if( !data || !data->hasFormat( type ) )
        return 0;

    QByteArray d( data->data( type ) );
    QDataStream s (d);

    void * ptr = 0;


    int bytesRead = s.readRawData( (char*)&ptr, int_c( sizeof( void * ) ) );
    WALBERLA_ASSERT_EQUAL( bytesRead, sizeof(void*) );
    WALBERLA_UNUSED(bytesRead);

    return ptr;
}




//======================================================================================================================
//
//  GRADIENT SELECTION
//
//======================================================================================================================

const QLinearGradient & getDefaultColormap()
{
   static bool firstRun=true;
   static QLinearGradient ret;
   if(!firstRun)
      return ret;


   firstRun = false;
   ret.setColorAt(1.0,QColor(0,0,255));
   ret.setColorAt(0.5,QColor(0,255,0));
   ret.setColorAt(0,QColor(255,255,255));

   return ret;
}

bool showColormapSelect(QLinearGradient & grad)
{
   bool ok=false;

   QGradient newColorMap = QtGradientDialog::getGradient(&ok,grad,0,"Select Color Map" );
   if(!ok)
      return false;

   QLinearGradient lastGrad;
   QLinearGradient colorMapGradient(0,0,100,0);
   for(int i=0; i< newColorMap.stops().size(); ++i)
      colorMapGradient.setColorAt(newColorMap.stops()[i].first, newColorMap.stops()[i].second );

   grad = colorMapGradient;


   return true;
}






} // namespace gui
} // namespace walberla
