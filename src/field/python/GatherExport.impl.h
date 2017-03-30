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
//! \file GatherExport.impl.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "field/Gather.h"
#include "python_coupling/helper/MplHelpers.h"
#include "python_coupling/helper/BlockStorageExportHelpers.h"
#include "python_coupling/helper/ModuleScope.h"
#include "python_coupling/helper/SliceToCellInterval.h"


namespace walberla {
namespace field {


namespace internal {

   //===================================================================================================================
   //
   //  Gather
   //
   //===================================================================================================================


   template<typename Field_T>
   boost::python::object gatherToObject( const shared_ptr<StructuredBlockStorage> & blocks, BlockDataID fieldID,
                                         CellInterval boundingBox = CellInterval(), int targetRank = 0 )
   {
      typedef Field< typename Field_T::value_type, Field_T::F_SIZE > ResultField;
      auto result = make_shared< ResultField > ( 0,0,0 );
      field::gather< Field_T, ResultField > ( *result, blocks, fieldID, boundingBox, targetRank, MPI_COMM_WORLD );

      if ( MPIManager::instance()->worldRank() == targetRank )
         return boost::python::object(result);
      else
         return boost::python::object();
   }

   FunctionExporterClass( gatherToObject,
                          boost::python::object( const shared_ptr<StructuredBlockStorage> &,
                                                 BlockDataID, CellInterval,int ) );

   template<typename FieldTypes>
   static boost::python::object gatherWrapper (  const shared_ptr<StructuredBlockStorage> & blocks, const std::string & blockDataStr,
                                                 const boost::python::tuple & slice,  int targetRank = 0 )
   {
      using namespace boost::python;

      auto fieldID = python_coupling::blockDataIDFromString( *blocks, blockDataStr );
      CellInterval boundingBox = python_coupling::globalPythonSliceToCellInterval( blocks, slice );

      if ( blocks->begin() == blocks->end() ) {
         // if no blocks are on this process the field::gather function can be called with any type
         // however we have to call it, otherwise a deadlock occurs
         gatherToObject< Field<real_t,1> > ( blocks, fieldID, boundingBox, targetRank );
         return object();
      }

      IBlock * firstBlock =  & ( * blocks->begin() );
      python_coupling::Dispatcher<FieldTypes, Exporter_gatherToObject > dispatcher( firstBlock );
      auto func = dispatcher( fieldID );
      if ( func.empty() )
      {
         PyErr_SetString( PyExc_RuntimeError, "This function cannot handle this type of block data.");
         throw error_already_set();
      }
      else
      {
         return func( blocks, fieldID, boundingBox, targetRank) ;
      }
   }

} // namespace internal



template<typename FieldTypes >
void exportGatherFunctions()
{
   using namespace boost::python;
   python_coupling::ModuleScope fieldModule( "field" );

   def( "gather",  &internal::gatherWrapper<FieldTypes>,  ( arg("blocks"), arg("blockDataName"), arg("slice"), arg("targetRank") = 0 ) );
}




} // namespace moduleName
} // namespace walberla


