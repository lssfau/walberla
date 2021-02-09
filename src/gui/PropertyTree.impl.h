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
//! \file PropertyTree.impl.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================



#ifdef WALBERLA_ENABLE_GUI

#include <QStandardItemModel>

namespace walberla {
namespace gui {

   template <typename T>
   PropertyTree::ItemID PropertyTree::addItem(const std::string & name, const T & val, ItemID  parent)
   {
      if(!parent)
         parent = model_->invisibleRootItem();

      QStandardItem * firstCol  = new QStandardItem(QString::fromStdString(name));
      QStandardItem * secondCol = new QStandardItem(QString("%1").arg(val));

      parent->appendRow(QList<ItemID>() << firstCol << secondCol);

      return firstCol;
   }
   template <>
   PropertyTree::ItemID  PropertyTree::addItem(const std::string & name,const std::string & val, ItemID parent);


   template <typename T>
   void PropertyTree::modifyItem(const std::string & name, const T & val, ItemID id)
   {
      id->setText(QString::fromStdString(name));
      id->parent()->child(id->row(),1)->setText(QString("%1").arg(val));
   }
   template <>
   void PropertyTree::modifyItem(const std::string & name, const std::string & val, ItemID id);


   template <typename T>
   void PropertyTree::modifyItem(const T & val, ItemID id)
   {
      id->parent()->child(id->row(),1)->setText(QString("%1").arg(val));
   }
   template <>
   void PropertyTree::modifyItem(const std::string & val, ItemID id);


} // namespace gui
} // namespace walberla

#else

namespace walberla {
namespace gui  {

   struct QStandardItem {};
   struct QStandardItemModel {};

   template <typename T>
   PropertyTree::ItemID PropertyTree::addItem(const std::string &, const T &, ItemID)
   {
      return nullptr;
   }

   template <typename T>
   void PropertyTree::modifyItem(const std::string &, const T &, ItemID)
   {
   }

   template <typename T>
   void PropertyTree::modifyItem(const T &, ItemID)
   {
   }

} // namespace gui
} // namespace walberla



#endif




