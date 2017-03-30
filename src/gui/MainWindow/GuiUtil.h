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
//! \file GuiUtil.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Utility functions for Qt GUI
//
//======================================================================================================================

#pragma once

#include <QComboBox>
#include <QString>
#include <QVariant>


//Forward Declarations
class QMimeData;
class QLinearGradient;
class QTreeWidgetItem;


namespace walberla {
namespace gui {



   //===================================================================================================================
   //
   //  MIME DATA PACKING/UNPACKING USED FOR DRAG & DROP
   //
   //===================================================================================================================

   /// Application specific MIME type, encoding a blockforest::Block pointer
   const QString BLOCK_MIMETYPE = "application/BlockPointer";

   /// Packs an arbitrary pointer in a QMimeData object
   QMimeData * createMimeDataFromPointer(void * p, const QString & role);

   /// Retrieves pointer from QMimeData, which was packed using createMimeDataFromPointer()
   void *      getPointerFromMimeData   (const QMimeData *data, const QString& type);



   template<class T>
   class VPtr
   {
   public:
      static T* asPtr( QVariant v ) {
         return (T *) v.value<void *> ();
      }

      static QVariant asQVariant( T* ptr )  {
         return qVariantFromValue ( (void *) ptr );
      }
   };



   //===================================================================================================================
   //
   //  GRADIENT SELECTION
   //
   //===================================================================================================================


   //*******************************************************************************************************************
   /*!\brief Shows a gradient selection dialog
    *
    * \ingroup gui
    *
    * \param[out] grad  The selected gradient
    * \return           false if the user canceled the selection
    *******************************************************************************************************************/
   bool showColormapSelect(QLinearGradient & grad);
   const QLinearGradient & getDefaultColormap();








} // namespace gui
} // namespace walberla

