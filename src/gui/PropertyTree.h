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
//! \file PropertyTree.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Interface for displaying tree based key-value information
//
//======================================================================================================================

#pragma once



class QStandardItemModel;
class QStandardItem;

#include "domain_decomposition/IBlock.h"
#include <string>

namespace walberla {
namespace gui {

   /**
    * Interface for providing Key-Value pairs for Qt Tree/List Views
    *
    * \ingroup gui
    *
    * Usage:
    *    - Subclass from PropertyTree and implement abstract methods
    *    - use fillTree() method to provide data using the protected methods addItem
    *    - the returned ItemID's can be used to modify, remove and add children to the item
    *      (modifyItem(), removeItem(), addItem() )
    *    - PropertyTree's can be registered in the gui::Manager for display
    */
   class PropertyTree
   {
      public:

#        ifdef WALBERLA_ENABLE_GUI
         typedef QStandardItem*  ItemID;
#        else
         using ItemID = void *;
#        endif

         PropertyTree();
         virtual ~PropertyTree();

         /************** @name Methods used by gui to retrieve data ********************************/
         //@{

         /**
          * Updates Model, using given block
          *
          * This method is called by gui when a different block should be displayed,
          * or the data inside a block has changed.
          * Should not be changed by subclasses
          */
         void setBlock( IBlock * b);

         /**
          * Retrieves the Qt Model Data structure
          *
          * used by gui, do not change in subclass
          */
         QStandardItemModel *      getModel()   { return model_; }

         /**
          * Returns name/type of data of this model
          *
          * When multiple PropertyTrees have been registered
          * the user can choose which data to display using this name
          */
         virtual const std::string getName() = 0;

         /**
          * Virtual constructor
          *
          * a new object has to be instantiated for each opened treeView
          * that displays this model
          * In the implementation just return "make_shared<YourDerivedClass>()"
          */
         virtual shared_ptr<PropertyTree> create() const = 0;

          //@}
      protected:

         /************** @name Abstract Methods   **************************************************/
         //@{

         /**
          * Method to fill the tree
          *
          * This method is called whenever the data in the block may have changed
          * or the block itself has changed. The first time fillTree is called add
          * your display data using addItem and store the ItemIDs.
          * On the next call, change the value of the items using modifyItem
          * If the block is not the same as on the last call, clear() is called
          * between the two fillTree invocations. In the clear method the saved ItemID's
          * can be deleted
          */
         virtual void fillTree( IBlock * b ) = 0;

         /**
          * Removes all data from the Property Tree
          *
          * This method should restore the object to its initial state
          * i.e. delete all stored data, usually the ItemID's
          * see also fillTree()
          */
         virtual void clear() = 0;

         //@}


         /************** @name Methods to add Data   ********************************************/
         //@{
         /**
          * Adds property value pair to the tree
          * @param name   the string that appears in the first column of the treeView
          * @param val    value that appears in the second column
          *               T can be every type that QString.arg() accepts
          *               i.e. all native data types and std::string;
          * @param parent key-value pairs can be organized hierarchically, here an ItemID
          *               of another item can be specified, which is used as parent
          * @return       identifier for this item which can be used to modify/remove this item or
          *               to add children
          */
         template <typename T>
         ItemID addItem(const std::string & name, const T & val, ItemID parent=nullptr);

         /**
          * Convenience method, behaves like addItem above with no value -> empty second column
          */
         ItemID addItem(const std::string & name, ItemID parent=nullptr);

         /**
          * Changes an existing item
          *
          * Store the ItemID's between fillTree calls, to modify items instead of
          * clearing and new adding, since otherwise the TreeView is reset after each timestep
          * and expansion states etc. are lost
          *
          * @param name new first column string for the item
          * @param val  new value for second column
          * @param id   id returned by addItem to specify which item is modified
          */
         template <typename T>
         void modifyItem(const std::string & name, const T & val, ItemID id );

         /**
          * Convenience method, behaves like modifyItem(string,T,ItemID) but does not change
          * the first column
          */
         template <typename T>
         void modifyItem(const T & val, ItemID id );


         /**
          * Removes an item from the tree
          *
          * To remove all items, call clear() instead of this method
          *
          * @param id  item to remove
          *
          */
         void removeItem(ItemID id);

         //@}

      private:
         QStandardItemModel * model_;
         IBlock    * lastBlock_;
   };


} // namespace gui
} // namespace walberla

#include "PropertyTree.impl.h"

