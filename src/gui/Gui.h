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
//! \file Gui.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "PropertyTree.h"
#include "blockforest/StructuredBlockForest.h"
#include "timeloop/ITimeloop.h"

#include <functional>
#include <vector>


namespace walberla {
namespace gui {


   struct GuiImpl; // pimpl pattern, to avoid the ifdef ENABLE_GUI in the header
   class  DisplayAdaptor;

   class GUI
   {
   public:

      GUI(timeloop::ITimeloop & timeloop, const shared_ptr<StructuredBlockForest> & blockForest, int& argc, char ** argv);
      ~GUI();


      //** Interface to User        ************************************************************************************
      /*! \name Interface to User  */
      //@{

      void run();

      void registerPropertyTree( const shared_ptr<PropertyTree>& propertyTree );

      using DisplayAdaptorCreatorFunc = std::function<DisplayAdaptor *(const IBlock &, ConstBlockDataID)>;
      void registerDisplayAdaptorCreator( const DisplayAdaptorCreatorFunc & creatorFunc );

      static void breakpoint( const std::string & comment, const std::string & file, int line );

      //@}
      //****************************************************************************************************************


      //** Interface to GUI        *************************************************************************************
      /*! \name Interface to GUI  */
      //@{

      void createAdaptors( std::vector<DisplayAdaptor*> & out ) const;

      const std::vector<shared_ptr<PropertyTree> > & getPropertyTrees() const;

      //@}
      //****************************************************************************************************************


   private:
      DisplayAdaptor * findDisplayAdaptorForBlockID ( ConstBlockDataID bdId ) const;

      std::vector<DisplayAdaptorCreatorFunc> displayAdaptorCreatorFuncs_;

      timeloop::ITimeloop   & timeloop_;
      StructuredBlockForest & blockForest_;

      GuiImpl * pImpl;

      static GUI * lastInstance_;
   };




} // namespace gui
} // namespace walberla



#ifdef WALBERLA_ENABLE_GUI

#define WALBERLA_GUI_BREAKPOINT(msg) {\
   std::ostringstream stringStream;\
   stringStream << msg;\
   ::walberla::gui::GUI::breakpoint(stringStream.str(), __FILE__, __LINE__ ); \
}
#else

#define WALBERLA_GUI_BREAKPOINT(msg)

#endif



//======================================================================================================================
//
//  EXPORTS
//
//======================================================================================================================

namespace walberla {
   using gui::GUI;
}
