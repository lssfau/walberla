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
//! \file CommunicationExport.impl.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "python_coupling/PythonWrapper.h"


#ifdef WALBERLA_BUILD_WITH_PYTHON


#include "field/communication/PackInfo.h"
#include "field/communication/UniformMPIDatatypeInfo.h"

#include "python_coupling/helper/MplHelpers.h"
#include "python_coupling/helper/BoostPythonHelpers.h"
#include "python_coupling/helper/MplHelpers.h"


namespace walberla {
namespace field {


namespace internal {

   //===================================================================================================================
   //
   //  createPackInfo Export
   //
   //===================================================================================================================

   template< typename FieldType >
   boost::python::object createPackInfoToObject( BlockDataID bdId, uint_t numberOfGhostLayers )
   {
      typedef typename FieldType::value_type T;
      const uint_t F_SIZE = FieldType::F_SIZE;
      typedef GhostLayerField<T,F_SIZE> GlField_T;
      if ( numberOfGhostLayers > 0  )
         return boost::python::object( make_shared< field::communication::PackInfo<GlField_T> >( bdId, numberOfGhostLayers ) );
      else
         return boost::python::object( make_shared< field::communication::PackInfo<GlField_T> >( bdId ) );
   }

   FunctionExporterClass( createPackInfoToObject, boost::python::object( BlockDataID, uint_t  ) );

   template< typename FieldVector>
   boost::python::object createPackInfo( const shared_ptr<StructuredBlockStorage> & bs,
                                         const std::string & blockDataName, uint_t numberOfGhostLayers )
   {
      auto bdId = python_coupling::blockDataIDFromString( *bs, blockDataName );
      if ( bs->begin() == bs->end() ) {
         // if no blocks are on this field an arbitrary PackInfo can be returned
         return createPackInfoToObject< GhostLayerField<real_t,1> > ( bdId, numberOfGhostLayers );
      }

      IBlock * firstBlock =  & ( * bs->begin() );
      python_coupling::Dispatcher<FieldVector, Exporter_createPackInfoToObject > dispatcher( firstBlock );
      return dispatcher( bdId )( bdId, numberOfGhostLayers ) ;
   }


   //===================================================================================================================
   //
   //  createMPIDatatypeInfo
   //
   //===================================================================================================================


   template< typename FieldType >
   boost::python::object createMPIDatatypeInfoToObject( BlockDataID bdId, uint_t numberOfGhostLayers )
   {
      typedef typename FieldType::value_type T;
      const uint_t F_SIZE = FieldType::F_SIZE;
      typedef GhostLayerField<T,F_SIZE> GlField_T;
      using field::communication::UniformMPIDatatypeInfo;

      if ( numberOfGhostLayers > 0 )
         return boost::python::object( make_shared< UniformMPIDatatypeInfo<GlField_T> >( bdId, numberOfGhostLayers ) );
      else
         return boost::python::object( make_shared< UniformMPIDatatypeInfo<GlField_T> >( bdId ) );
   }

   FunctionExporterClass( createMPIDatatypeInfoToObject, boost::python::object( BlockDataID, uint_t  ) );

   template< typename FieldVector>
   boost::python::object createMPIDatatypeInfo( const shared_ptr<StructuredBlockStorage> & bs,
                                                const std::string & blockDataName,
                                                uint_t numberOfGhostLayers)
   {
      auto bdId = python_coupling::blockDataIDFromString( *bs, blockDataName );
      if ( bs->begin() == bs->end() ) {
         // if no blocks are on this field an arbitrary MPIDatatypeInfo can be returned
         return createMPIDatatypeInfoToObject< GhostLayerField<real_t,1> > ( bdId, numberOfGhostLayers );
      }

      IBlock * firstBlock =  & ( * bs->begin() );
      python_coupling::Dispatcher<FieldVector, Exporter_createMPIDatatypeInfoToObject > dispatcher( firstBlock );
      return dispatcher( bdId )( bdId, numberOfGhostLayers );
   }


} // namespace internal





template<typename FieldTypes>
void exportCommunicationClasses()
{
   using namespace boost::python;

   def( "createMPIDatatypeInfo",&internal::createMPIDatatypeInfo<FieldTypes>, ( arg("blocks"), arg("blockDataName"), arg("numberOfGhostLayers" ) =0 ) );
   def( "createPackInfo",       &internal::createPackInfo<FieldTypes>,        ( arg("blocks"), arg("blockDataName"), arg("numberOfGhostLayers" ) =0 ) );
}


} // namespace moduleName
} // namespace walberla




#endif // WALBERLA_BUILD_WITH_PYTHON
