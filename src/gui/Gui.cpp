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
//! \file Gui.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "Gui.h"
#include "waLBerlaDefinitions.h"


#ifdef WALBERLA_ENABLE_GUI

#include "MainWindow.h"
#include "ScalarFieldDisplayAdaptor.h"
#include "FlagFieldDisplayAdaptor.h"
#include "VectorFieldDisplayAdaptor.h"
#include "ScalarField3DisplayAdaptor.h"

#include "core/Abort.h"
#include "core/logging/Logging.h"

#include <QApplication>
#include <stdexcept>


namespace walberla {
namespace gui {

void throwOnAbort( const std::string & message )
{
   throw std::runtime_error( message );
}


// last instance pointer -> points to GUI that handles breakpoints
GUI * GUI::lastInstance_ = NULL;

struct GuiImpl
{
   QApplication * app;
   MainWindow   * mainWindow;
   std::vector< shared_ptr<PropertyTree>  > propertyTrees_;
};



GUI::GUI(timeloop::ITimeloop & timeloop, const shared_ptr<StructuredBlockForest> & blockForest, int& argc, char ** argv)
   :  timeloop_( timeloop ),
      blockForest_( *blockForest )
{
   QCoreApplication::setOrganizationName( "FAU"      );
   QCoreApplication::setApplicationName ( "waLBerla" );

   pImpl = new GuiImpl();
   pImpl->app = new QApplication( argc, argv );
   pImpl->mainWindow = new MainWindow( timeloop_, blockForest_, *this );

   if ( lastInstance_ != NULL ) {
      WALBERLA_LOG_WARNING( "More than one GUI may lead to problems!");
   }
   lastInstance_ = this;

   walberla::Abort::instance()->resetAbortFunction( &walberla::Abort::exceptionAbort );
}

GUI::~GUI()
{
   if ( lastInstance_ == this )
      lastInstance_ = NULL;

   delete pImpl->mainWindow;
   delete pImpl->app;
   delete pImpl;
}


void GUI::registerDisplayAdaptorCreator( const DisplayAdaptorCreatorFunc & creatorFunc )
{
   displayAdaptorCreatorFuncs_.push_back( creatorFunc );
}


DisplayAdaptor * GUI::findDisplayAdaptorForBlockID ( ConstBlockDataID bdId ) const
{
   for( auto blockIt = blockForest_.begin(); blockIt != blockForest_.end(); ++blockIt )
   {
      IBlock & block = *blockIt;
      if( ! block.isBlockDataAllocated(bdId) )
         continue;


      // User registered Adaptors
      for( auto it = displayAdaptorCreatorFuncs_.begin(); it != displayAdaptorCreatorFuncs_.end(); ++it )
      {
         DisplayAdaptor * result = (*it)( block, bdId);
         if ( result )
            return result;
      }

      // Flag fields
      if ( block.isDataOfType< FlagField< uint8_t> >( bdId) ) return new FlagFieldDisplayAdaptor<FlagField< uint8_t> >( bdId );
      if ( block.isDataOfType< FlagField<uint16_t> >( bdId) ) return new FlagFieldDisplayAdaptor<FlagField<uint16_t> >( bdId );
      if ( block.isDataOfType< FlagField<uint32_t> >( bdId) ) return new FlagFieldDisplayAdaptor<FlagField<uint32_t> >( bdId );
      if ( block.isDataOfType< FlagField<uint64_t> >( bdId) ) return new FlagFieldDisplayAdaptor<FlagField<uint64_t> >( bdId );

      // Ghost Layer fields
      if ( block.isDataOfType< GhostLayerField<uint8_t ,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<uint8_t ,1>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<uint16_t,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<uint16_t,1>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<uint32_t,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<uint32_t,1>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<uint64_t,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<uint64_t,1>   >( bdId );

      if ( block.isDataOfType< GhostLayerField<int8_t ,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<int8_t ,1>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<int16_t,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<int16_t,1>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<int32_t,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<int32_t,1>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<int64_t,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<int64_t,1>   >( bdId );

      if ( block.isDataOfType< GhostLayerField<float ,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<float, 1>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<double,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<double,1>   >( bdId );

      if ( block.isDataOfType< GhostLayerField<Vector3<real_t>,1> >( bdId) ) return new VectorFieldDisplayAdaptor<GhostLayerField<Vector3<real_t>,1> >( bdId );
      if ( block.isDataOfType< GhostLayerField<real_t,3> >         ( bdId) ) return new ScalarField3DisplayAdaptor<GhostLayerField<real_t,3> >        ( bdId );


      if ( block.isDataOfType< GhostLayerField<real_t,2>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<real_t, 2>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<real_t,3>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<real_t, 3>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<real_t,4>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<real_t, 4>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<real_t,5>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<real_t, 5>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<real_t,6>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<real_t, 6>   >( bdId );
      if ( block.isDataOfType< GhostLayerField<real_t,7>   >( bdId) ) return new ScalarFieldDisplayAdaptor<GhostLayerField<real_t, 7>   >( bdId );


      // Fields
      if ( block.isDataOfType< Field<uint8_t ,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<Field<uint8_t ,1>   >( bdId );
      if ( block.isDataOfType< Field<uint16_t,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<Field<uint16_t,1>   >( bdId );
      if ( block.isDataOfType< Field<uint32_t,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<Field<uint32_t,1>   >( bdId );
      if ( block.isDataOfType< Field<uint64_t,1>   >( bdId) ) return new ScalarFieldDisplayAdaptor<Field<uint64_t,1>   >( bdId );

      if ( block.isDataOfType< Field<int8_t ,1>   >( bdId) )  return new ScalarFieldDisplayAdaptor<Field<int8_t ,1>   >( bdId );
      if ( block.isDataOfType< Field<int16_t,1>   >( bdId) )  return new ScalarFieldDisplayAdaptor<Field<int16_t,1>   >( bdId );
      if ( block.isDataOfType< Field<int32_t,1>   >( bdId) )  return new ScalarFieldDisplayAdaptor<Field<int32_t,1>   >( bdId );
      if ( block.isDataOfType< Field<int64_t,1>   >( bdId) )  return new ScalarFieldDisplayAdaptor<Field<int64_t,1>   >( bdId );

      if ( block.isDataOfType< Field<float ,1>   >( bdId) )   return new ScalarFieldDisplayAdaptor<Field<float, 1>   >( bdId );
      if ( block.isDataOfType< Field<double,1>   >( bdId) )   return new ScalarFieldDisplayAdaptor<Field<double,1>   >( bdId );

      if ( block.isDataOfType< Field<Vector3<real_t>,1> >( bdId) ) return new VectorFieldDisplayAdaptor< Field<Vector3<real_t>,1> >( bdId );

   }

   return NULL;
}


void GUI::registerPropertyTree( const shared_ptr<PropertyTree>& propertyTree ) {
   pImpl->propertyTrees_.push_back( propertyTree );
}


const std::vector<shared_ptr<PropertyTree> > & GUI::getPropertyTrees() const {
   return pImpl->propertyTrees_;
}


void GUI::createAdaptors( std::vector<DisplayAdaptor*> & out ) const
{
   for( uint_t i = 0; i< blockForest_.numberOfBlockDataItems(); ++i )
   {
      DisplayAdaptor * da = findDisplayAdaptorForBlockID( ConstBlockDataID(i) );
      if( da )
         out.push_back( da );
   }
}


void GUI::run()
{

   pImpl->mainWindow->showMaximized();
   pImpl->app->exec();
}




void GUI::breakpoint( const std::string & comment, const std::string & file, int line )
{
   if ( GUI::lastInstance_ == NULL )
      return;

   GUI::lastInstance_->pImpl->mainWindow->breakpoint( QString::fromStdString(comment),
                                                    QString::fromStdString( file ),
                                                    line );
}



} // namespace gui
} // namespace walberla



#else //WALBERLA_ENABLE_GUI


namespace walberla {
namespace gui {

GUI * GUI::lastInstance_ = nullptr;

GUI::GUI(timeloop::ITimeloop & timeloop, const shared_ptr<StructuredBlockForest> & blockForest, int& , char ** )
   :  timeloop_(timeloop),
      blockForest_( *blockForest )
{
}

GUI::~GUI() = default;

void GUI::run() {
   timeloop_.run();
}

void GUI::registerPropertyTree( const shared_ptr<PropertyTree>&  ) {}

const std::vector<shared_ptr<PropertyTree> > & GUI::getPropertyTrees() const  {
   throw( "Should not happen!" );

}

void breakpoint( const std::string & , const std::string & , int  ){}


void GUI::registerDisplayAdaptorCreator( const DisplayAdaptorCreatorFunc &  )
{

}



} // namespace gui
} // namespace walberla


#endif //WALBERLA_ENABLE_GUI








