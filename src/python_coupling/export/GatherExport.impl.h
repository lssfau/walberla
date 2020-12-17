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
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "field/Gather.h"
#include "python_coupling/helper/MplHelpers.h"
#include "python_coupling/helper/BlockStorageExportHelpers.h"
#include "python_coupling/helper/SliceToCellInterval.h"


namespace walberla {
namespace field {

//*******************************************************************************************************************
/*! Exports the gather functionality of waLberla
*
* With field.gather a corresponding field will the gathered to the specified process. This field can be viewed as a
* numpy array with field.toArrayOn all other porcesses an empty pybind11::object will be returned.
*
* \hint For large scale simulations it is also possible to provide a slice to keep the gathered data low!
*/
//*******************************************************************************************************************
namespace py = pybind11;
template<typename... FieldTypes >
void exportGatherFunctions(py::module_ &m);

namespace internal {
namespace py = pybind11;
   //===================================================================================================================
   //
   //  Gather
   //
   //===================================================================================================================

class GatherExporter
{
 public:
   GatherExporter(const shared_ptr<StructuredBlockForest> & blocks, ConstBlockDataID fieldId,
                  CellInterval boundingBox = CellInterval(),  int targetRank = 0 )
      : blocks_( blocks ), fieldId_(fieldId), boundingBox_( boundingBox ), targetRank_(targetRank)
   {}

   template< typename FieldType>
   void operator() ( python_coupling::NonCopyableWrap<FieldType> )
   {
      typedef Field< typename FieldType::value_type, FieldType::F_SIZE > ResultField;
      IBlock * firstBlock =  & ( * blocks_->begin() );
      if( firstBlock->isDataClassOrSubclassOf<FieldType>(fieldId_) )
      {
         auto result = make_shared< ResultField > ( 0,0,0 );
         field::gather< FieldType, ResultField > ( *result, blocks_, fieldId_, boundingBox_, targetRank_, MPI_COMM_WORLD );

         if ( MPIManager::instance()->worldRank() == targetRank_ )
            resultField_ = py::cast(result);
         else
            resultField_ = py::none();

      }
   }
   py::object getResultField()
   {
      return resultField_;
   }

 private:
   py::object resultField_;
   shared_ptr< StructuredBlockStorage > blocks_;
   ConstBlockDataID fieldId_;
   std::string vtkName_;
   CellInterval boundingBox_;
   int targetRank_ ;
};


template<typename... FieldTypes>
static py::object gatherWrapper(const shared_ptr<StructuredBlockForest> & blocks, const std::string & name,
                                const py::tuple & slice, int targetRank = 0 )
{
   BlockDataID fieldID = python_coupling::blockDataIDFromString( *blocks, name );
   CellInterval boundingBox = python_coupling::globalPythonSliceToCellInterval( blocks, slice );

   if ( blocks->begin() == blocks->end() ) {
      // if no blocks are on this process the field::gather function can be called with any type
      // however we have to call it, otherwise a deadlock occurs
      auto result = make_shared< Field<real_t, 1> > ( 0,0,0 );
      field::gather< Field<real_t, 1>, Field<real_t, 1> > ( *result, blocks, fieldID, boundingBox, targetRank, MPI_COMM_WORLD );
      return py::none();
   }

   GatherExporter exporter(blocks, fieldID, boundingBox, targetRank);
   python_coupling::for_each_noncopyable_type< FieldTypes... >  ( std::ref(exporter) );

   if ( ! exporter.getResultField() ) {
      throw py::value_error("Failed to gather Field");
   }
   else {
      return exporter.getResultField();
   }
}

} // namespace internal


namespace py = pybind11;
using namespace pybind11::literals;
template<typename... FieldTypes >
void exportGatherFunctions(py::module_ &m)
{
   py::module_ m2 = m.def_submodule("field", "Field Extension of the waLBerla python bindings");

   m2.def(
      "gather",
      [](const shared_ptr<StructuredBlockForest> & blocks, const std::string & name,
         const py::tuple & slice, int targetRank = 0 ) {
        return internal::gatherWrapper< FieldTypes... >(blocks, name, slice, targetRank);
      },
      "blocks"_a, "name"_a, "slice"_a, "targetRank"_a = uint_t(0));
}

} // namespace moduleName
} // namespace walberla


