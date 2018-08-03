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
//! \file PropertyTree.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "PropertyTree.h"
#include "core/debug/Debug.h"


#ifdef WALBERLA_ENABLE_GUI

namespace walberla {
namespace gui {


   PropertyTree::PropertyTree()
      : model_(0),lastBlock_(0)
   {
      model_ = new QStandardItemModel();
      model_->setColumnCount(2);
   }


   PropertyTree::~PropertyTree()
   {
      delete model_;
   }


   void PropertyTree::setBlock( IBlock * b)
   {
      if(b != lastBlock_)
      {
         model_->clear();
         this->clear();
      }
      lastBlock_ = b;
      fillTree(b);
   }

   PropertyTree::ItemID PropertyTree::addItem(const std::string & name, ItemID parent)
   {
      if(!parent)
         parent = model_->invisibleRootItem();

      WALBERLA_ASSERT_NOT_NULLPTR( parent );

      QStandardItem * firstCol  = new QStandardItem(QString::fromStdString(name) );
      QStandardItem * secondCol = new QStandardItem("");
      firstCol->setEditable(false);
      secondCol->setEditable(false);

      parent->appendRow(QList<ItemID>() << firstCol << secondCol);

      return firstCol;
   }

   template <>
   PropertyTree::ItemID  PropertyTree::addItem(const std::string & name,const std::string & val, ItemID parent)
   {
      if(!parent)
         parent = model_->invisibleRootItem();

      QStandardItem * firstCol  = new QStandardItem(QString::fromStdString(name));
      QStandardItem * secondCol = new QStandardItem(QString::fromStdString(val));
      firstCol->setEditable(false);
      secondCol->setEditable(false);
      parent->appendRow(QList<ItemID>() << firstCol << secondCol );

      return firstCol;
   }


   template <>
   void PropertyTree::modifyItem(const std::string & name, const std::string & val, ItemID id)
   {
      id->setText(QString::fromStdString(name));
      model_->item(id->row(),1)->setText(QString::fromStdString(val));
   }

   template <>
   void PropertyTree::modifyItem(const std::string & val, ItemID id)
   {
      model_->item(id->row(),1)->setText(QString::fromStdString(val));
   }


   void PropertyTree::removeItem(ItemID id)
   {
      QStandardItem * parent = id->parent();
      // due to a Qt Bug or (feature?!) parent is null when it should be invisibleRootItem
      if ( ! parent )
         parent = model_->invisibleRootItem();

      parent->removeRow(id->row() );
   }
} // namespace gui
} // namespace walberla

#else

namespace walberla {
namespace gui {
   PropertyTree::PropertyTree()
      : model_(nullptr), lastBlock_(nullptr)
   {}


   PropertyTree::~PropertyTree()
   = default;


   void PropertyTree::setBlock( IBlock * a)
   {
      //dummy to prevent lastBlock_ not used
      lastBlock_ = a;
   }

   void PropertyTree::removeItem(ItemID )
   {
   }

   PropertyTree::ItemID PropertyTree::addItem(const std::string & , ItemID )
   {
      return nullptr;
   }
} // namespace gui
} // namespace walberla

#endif




