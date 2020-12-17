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
//! \file FieldExport.impl.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "core/VectorTrait.h"
#include "core/logging/Logging.h"

#include "field/AddToStorage.h"
#include "field/Field.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "field/communication/UniformMPIDatatypeInfo.h"
#include "field/vtk/FlagFieldMapping.h"
#include "field/vtk/VTKWriter.h"

#include "python_coupling/PythonWrapper.h"
#include "python_coupling/helper/MplHelpers.h"
#include "python_coupling/helper/PybindHelper.h"

#include <iostream>
#include <type_traits>

#include "GatherExport.impl.h"
#include "pybind11/numpy.h"
#include <pybind11/stl.h>

namespace walberla
{
namespace field
{

//*******************************************************************************************************************
/*! Exports all Fields given in the Sequence
*
* Put only Fields in the sequence! The corresponding GhostLayerFields are exported automatically
*
* \warning Make sure that the same adaptor type is exported only once!
*/
//*******************************************************************************************************************
template<typename... FieldTypes >
void exportFields();



//*******************************************************************************************************************
/*! Exports all GhostLayerFieldAdaptors given in the Sequence
*
* \warning Make sure that the same adaptor type is exported only once!
*/
//*******************************************************************************************************************
template<typename... AdaptorTypes>
void exportGhostLayerFieldAdaptors();

template<typename AdaptorType>
void exportGhostLayerFieldAdaptor();

namespace internal
{
namespace py = pybind11;

template<class T> struct PythonFormatString                    { inline static char * get() { static char value [] = "B"; return value; } };

template<>        struct PythonFormatString<double>            { inline static char * get() { static char value [] = "d"; return value; } };
template<>        struct PythonFormatString<float>             { inline static char * get() { static char value [] = "f"; return value; } };
template<>        struct PythonFormatString<unsigned short>    { inline static char * get() { static char value [] = "H"; return value; } };
template<>        struct PythonFormatString<int>               { inline static char * get() { static char value [] = "i"; return value; } };
template<>        struct PythonFormatString<unsigned int>      { inline static char * get() { static char value [] = "I"; return value; } };
template<>        struct PythonFormatString<long>              { inline static char * get() { static char value [] = "l"; return value; } };
template<>        struct PythonFormatString<unsigned long>     { inline static char * get() { static char value [] = "L"; return value; } };
template<>        struct PythonFormatString<long long>         { inline static char * get() { static char value [] = "q"; return value; } };
template<>        struct PythonFormatString<unsigned long long>{ inline static char * get() { static char value [] = "Q"; return value; } };
template<>        struct PythonFormatString<int8_t>            { inline static char * get() { static char value [] = "c"; return value; } };
template<>        struct PythonFormatString<int16_t>           { inline static char * get() { static char value [] = "h"; return value; } };
template<>        struct PythonFormatString<uint8_t>           { inline static char * get() { static char value [] = "C"; return value; } };

//===================================================================================================================
//
//  Aligned Allocation
//
//===================================================================================================================

template< typename T >
shared_ptr< field::FieldAllocator< T > > getAllocator(uint_t alignment)
{
   if (alignment == 0)
      return shared_ptr< field::FieldAllocator< T > >(); // leave to default - auto-detection of alignment
   else if (alignment == 16)
      return make_shared< field::AllocateAligned< T, 16 > >();
   else if (alignment == 32)
      return make_shared< field::AllocateAligned< T, 32 > >();
   else if (alignment == 64)
      return make_shared< field::AllocateAligned< T, 64 > >();
   else if (alignment == 128)
      return make_shared< field::AllocateAligned< T, 128 > >();
   else
   {
      throw py::value_error("Alignment parameter has to be one of 0, 16, 32, 64, 128.");
      return shared_ptr< field::FieldAllocator< T > >();
   }
}

template< typename GhostLayerField_T >
class GhostLayerFieldDataHandling : public field::BlockDataHandling< GhostLayerField_T >
{
 public:
   typedef typename GhostLayerField_T::value_type Value_T;

   GhostLayerFieldDataHandling(const weak_ptr< StructuredBlockStorage >& blocks, const uint_t nrOfGhostLayers,
                               const Value_T& initValue, const Layout layout, uint_t alignment = 0)
      : blocks_(blocks), nrOfGhostLayers_(nrOfGhostLayers), initValue_(initValue), layout_(layout),
        alignment_(alignment)
   {}

   GhostLayerField_T* allocate(IBlock* const block)
   {
      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blocks, "Trying to access 'AlwaysInitializeBlockDataHandling' for a block "
                                         "storage object that doesn't exist anymore");
      GhostLayerField_T* field = new GhostLayerField_T(
         blocks->getNumberOfXCells(*block), blocks->getNumberOfYCells(*block), blocks->getNumberOfZCells(*block),
         nrOfGhostLayers_, initValue_, layout_, getAllocator< Value_T >(alignment_));
      return field;
   }

   GhostLayerField_T* reallocate(IBlock* const block) { return allocate(block); }

 private:
   weak_ptr< StructuredBlockStorage > blocks_;

   uint_t nrOfGhostLayers_;
   Value_T initValue_;
   Layout layout_;
   uint_t alignment_;
};

//===================================================================================================================
//
//  Field functions redefined for easier export
//
//===================================================================================================================

template< typename Field_T >
py::object field_size(const Field_T& field)
{
   return py::make_tuple(field.xSize(), field.ySize(), field.zSize(), field.fSize());
}

template< typename GlField_T >
py::tuple field_sizeWithGhostLayer(const GlField_T& field)
{
   return py::make_tuple(field.xSizeWithGhostLayer(), field.ySizeWithGhostLayer(), field.zSizeWithGhostLayer(),
                         field.fSize());
}

template< typename Field_T >
py::tuple field_allocSize(const Field_T& field)
{
   return py::make_tuple(field.xAllocSize(), field.yAllocSize(), field.zAllocSize(), field.fAllocSize());
}

template< typename Field_T >
py::tuple field_strides(const Field_T& field)
{
   return py::make_tuple(field.xStride(), field.yStride(), field.zStride(), field.fStride());
}

template< typename Field_T >
py::tuple field_offsets(const Field_T& field)
{
   return py::make_tuple(field.xOff(), field.yOff(), field.zOff());
}

template< typename Field_T >
py::object field_layout(const Field_T& f)
{
   return py::cast(f.layout());
}

template< typename Field_T >
void field_swapDataPointers(Field_T& f1, Field_T& f2)
{
   if (!f1.hasSameAllocSize(f2) || !f1.hasSameSize(f2) || f1.layout() != f2.layout())
   {
      throw py::value_error("The data of fields with different sizes or layout cannot be swapped");
   }
   f1.swapDataPointers(f2);
}

template< typename Field_T >
py::object copyAdaptorToField(const Field_T& f)
{
   typedef GhostLayerField< typename Field_T::value_type, Field_T::F_SIZE > ResField;
   auto res = make_shared< ResField >(f.xSize(), f.ySize(), f.zSize(), f.nrOfGhostLayers());

   auto srcIt = f.beginWithGhostLayerXYZ();
   auto dstIt = res->beginWithGhostLayerXYZ();
   while (srcIt != f.end())
   {
      for (cell_idx_t fCoord = 0; fCoord < cell_idx_c(Field_T::F_SIZE); ++fCoord)
         dstIt.getF(fCoord) = srcIt.getF(fCoord);

      ++srcIt;
      ++dstIt;
   }
   return py::cast(res);
}

//===================================================================================================================
//
//  Field export
//
//===================================================================================================================
template< typename Field_T >
py::array_t< typename Field_T::value_type > toNumpyArray(const Field_T& field)
{
   using T    = typename Field_T::value_type;
   const T* ptr = field.dataAt(0, 0, 0, 0);

   if (field.fSize() == 1)
   {
      return pybind11::array_t< T, 0 >({ field.xSize(), field.ySize(), field.zSize() },
                                       { static_cast< size_t >(field.xStride()) * sizeof(T),
                                         static_cast< size_t >(field.yStride()) * sizeof(T),
                                         static_cast< size_t >(field.zStride()) * sizeof(T) },
                                       ptr, py::cast(field));
   }
   else
   {
      return pybind11::array_t< T, 0 >(
         { field.xSize(), field.ySize(), field.zSize(), field.fSize() },
         { static_cast< size_t >(field.xStride()) * sizeof(T), static_cast< size_t >(field.yStride()) * sizeof(T),
           static_cast< size_t >(field.zStride()) * sizeof(T), static_cast< size_t >(field.fStride()) * sizeof(T) },
         ptr, py::cast(field));
   }
}

template< typename GlField_T >
py::array_t< typename GlField_T::value_type > toNumpyArrayWithGhostLayers(const GlField_T& field)
{
   using T    = typename GlField_T::value_type;
   const T* ptr     = field.dataAt(-static_cast< cell_idx_t >(field.nrOfGhostLayers()),
                                   -static_cast< cell_idx_t >(field.nrOfGhostLayers()),
                                   -static_cast< cell_idx_t >(field.nrOfGhostLayers()), 0);


   if (field.fSize() == 1)
   {
      return pybind11::array_t< T, 0 >({ field.xSizeWithGhostLayer(), field.ySizeWithGhostLayer(), field.zSizeWithGhostLayer() },
                                       { static_cast< size_t >(field.xStride()) * sizeof(T),
                                         static_cast< size_t >(field.yStride()) * sizeof(T),
                                         static_cast< size_t >(field.zStride()) * sizeof(T) },
                                       ptr, py::cast(field));
   }
   else
   {
      return pybind11::array_t< T, 0 >(
         { field.xSizeWithGhostLayer(), field.ySizeWithGhostLayer(), field.zSizeWithGhostLayer(), field.fSize() },
         { static_cast< size_t >(field.xStride()) * sizeof(T), static_cast< size_t >(field.yStride()) * sizeof(T),
           static_cast< size_t >(field.zStride()) * sizeof(T), static_cast< size_t >(field.fStride()) * sizeof(T) },
         ptr, py::cast(field));
   }
}


struct FieldExporter
{
   FieldExporter(py::module_& m) : m_(m) {}
   template< typename FieldType >
   void operator()(python_coupling::NonCopyableWrap< FieldType >) const
   {
      typedef typename FieldType::value_type T;
      const uint_t F_SIZE = FieldType::F_SIZE;
      typedef GhostLayerField< T, F_SIZE > GlField_T;
      typedef Field< T, F_SIZE > Field_T;

      std::string data_type_name = PythonFormatString<T>::get();

      std::string class_name = "Field_" + data_type_name + "_" + std::to_string(FieldType::F_SIZE);

      py::class_< Field_T, shared_ptr< Field_T > >(m_, class_name.c_str())
         .def_property_readonly("layout", &field_layout< Field_T >)
         .def_property_readonly("size", &field_size< Field_T >)
         .def_property_readonly("allocSize", &field_allocSize< Field_T >)
         .def_property_readonly("strides", &field_strides< Field_T >)
         .def_property_readonly("offsets", &field_offsets< Field_T >)
         .def("clone", &Field_T::clone, py::return_value_policy::copy)
         .def("cloneUninitialized", &Field_T::cloneUninitialized, py::return_value_policy::copy)
         .def("swapDataPointers", &field_swapDataPointers< Field_T >)
         .def("__getitem__",        [](const Field_T& self, const py::object& index) {
             return py::cast(self).attr("__array__")().attr("__getitem__")(index);
         } )
          .def("__setitem__",        [](const Field_T& self, const py::object& index,
                                        const typename Field_T::value_type& value) {
              py::cast(self).attr("__array__")().attr("__setitem__")(index, value);
          } )
         .def("__array__", &toNumpyArray< Field_T >);

      std::string class_nameGL =
         "GhostLayerField_" + data_type_name + "_" + std::to_string(FieldType::F_SIZE);

      py::class_< GlField_T, shared_ptr< GlField_T >, Field_T >(m_, class_nameGL.c_str())
         .def_property_readonly("sizeWithGhostLayer", &GlField_T::xSizeWithGhostLayer)
         .def_property_readonly("nrOfGhostLayers", &GlField_T::nrOfGhostLayers)
         .def("__array__", &toNumpyArrayWithGhostLayers< GlField_T >);

      using field::communication::PackInfo;
      std::string FieldPackInfo_name = "FieldPackInfo_" + data_type_name + "_" + std::to_string(FieldType::F_SIZE);
      py::class_< PackInfo< GlField_T >, shared_ptr< PackInfo< GlField_T > >, walberla::communication::UniformPackInfo >(m_, FieldPackInfo_name.c_str());

      using field::communication::UniformMPIDatatypeInfo;
      std::string FieldMPIDataTypeInfo_name = "FieldMPIDataTypeInfo_" + data_type_name + "_" + std::to_string(FieldType::F_SIZE);
      py::class_< UniformMPIDatatypeInfo< GlField_T >, shared_ptr< UniformMPIDatatypeInfo< GlField_T > >, walberla::communication::UniformMPIDatatypeInfo >(
         m_, FieldMPIDataTypeInfo_name.c_str());
   }
   const py::module_& m_;
};


struct FieldAllocatorExporter
{
   FieldAllocatorExporter(py::module_& m) : m_(m) {}
   template< typename T >
   void operator()(python_coupling::NonCopyableWrap< T >) const
   {
      std::string data_type_name = PythonFormatString<T>::get();
      std::string class_nameFieldAllocator = "FieldAllocator_" + data_type_name;
      py::class_< FieldAllocator< T >, shared_ptr< FieldAllocator< T > > >(m_, class_nameFieldAllocator.c_str())
         .def("incrementReferenceCount", &FieldAllocator< T >::incrementReferenceCount)
         .def("decrementReferenceCount", &FieldAllocator< T >::decrementReferenceCount);
   }
   const py::module_& m_;
};

//===================================================================================================================
//
//  addToStorage
//
//===================================================================================================================

class AddToStorageExporter
{
 public:
   AddToStorageExporter(const shared_ptr< StructuredBlockForest >& blocks, const std::string& name, py::object& dtype, uint_t fs,
                        uint_t gl, Layout layout, real_t initValue, uint_t alignment)
      : blocks_(blocks), name_(name), dtype_(dtype), fs_(fs), gl_(gl), layout_(layout), initValue_(initValue), alignment_(alignment), found_(false)
   {}

   template< typename FieldType >
   void operator()(python_coupling::NonCopyableWrap<FieldType>) const
   {
      using namespace py;
      typedef typename FieldType::value_type T;
      const uint_t F_SIZE = FieldType::F_SIZE;

      if (F_SIZE != fs_) return;
      if(python_coupling::isCppEqualToPythonType<T>(py::cast<std::string>(dtype_.attr("__name__"))))
      {
         typedef internal::GhostLayerFieldDataHandling< GhostLayerField< T, F_SIZE > > DataHandling;
         auto dataHandling = walberla::make_shared< DataHandling >(blocks_, gl_, initValue_, layout_, alignment_);
         blocks_->addBlockData(dataHandling, name_);
      }
      found_ = true;
   }

   bool successful() const { return found_; }

 private:
   shared_ptr< StructuredBlockStorage > blocks_;
   std::string name_;
   py::object dtype_;
   uint_t fs_;
   uint_t gl_;
   Layout layout_;
   real_t initValue_;
   uint_t alignment_;
   mutable bool found_;
};

template< typename... FieldTypes >
void addToStorage(const shared_ptr< StructuredBlockForest >& blocks, const std::string& name, py::object& dtype,
                  uint_t fs, uint_t gl, Layout layout, real_t initValue, uint_t alignment)
{
   using namespace py;

   auto result = make_shared< py::object >();
   AddToStorageExporter exporter(blocks, name, dtype, fs, gl, layout, initValue, alignment);
   python_coupling::for_each_noncopyable_type< FieldTypes... >(exporter);

   if (!exporter.successful())
   {
      throw py::value_error("Adding GhostLayerField failed. Maybe the data type and/or the fsize is not exported to python yet");
   }
}

inline void addFlagFieldToStorage(const shared_ptr< StructuredBlockStorage >& blocks, const std::string& name,
                                  uint_t nrOfBits, uint_t gl)
{
   if (nrOfBits == 8)
      field::addFlagFieldToStorage< FlagField< uint8_t > >(blocks, name, gl);
   else if (nrOfBits == 16)
      field::addFlagFieldToStorage< FlagField< uint16_t > >(blocks, name, gl);
   else if (nrOfBits == 32)
      field::addFlagFieldToStorage< FlagField< uint32_t > >(blocks, name, gl);
   else if (nrOfBits == 64)
      field::addFlagFieldToStorage< FlagField< uint64_t > >(blocks, name, gl);
   else
   {
      throw py::value_error("Allowed values for number of bits are: 8,16,32,64");
   }
}

//===================================================================================================================
//
//  createField
//
//===================================================================================================================

class CreateFieldExporter
{
 public:
   CreateFieldExporter( uint_t xs, uint_t ys, uint_t zs, uint_t fs, uint_t gl,
                        Layout layout, const py::object & dtype, uint_t alignment,
                        const shared_ptr<py::object> & resultPointer  )
      : xs_( xs ), ys_(ys), zs_(zs), fs_(fs), gl_(gl),
        layout_( layout),  dtype_( dtype ), alignment_(alignment), resultPointer_( resultPointer )
   {}

   template< typename FieldType>
   void operator() ( python_coupling::NonCopyableWrap<FieldType> ) const
   {
      typedef typename FieldType::value_type T;
      const uint_t F_SIZE = FieldType::F_SIZE;

      if( F_SIZE != fs_ )
         return;

      if(python_coupling::isCppEqualToPythonType<T>(py::cast<std::string>(dtype_.attr("__name__"))))
      {
         T initVal = T();
         *resultPointer_ = py::cast( make_shared< GhostLayerField<T, F_SIZE> >( xs_,ys_,zs_, gl_, initVal, layout_,
                                                                             getAllocator<T>(alignment_)));
      }
   }

 private:
   uint_t xs_;
   uint_t ys_;
   uint_t zs_;
   uint_t fs_;
   uint_t gl_;
   Layout layout_;
   py::object dtype_;
   uint_t alignment_;
   shared_ptr<py::object> resultPointer_;
};

template<typename... FieldTypes>
py::object createPythonField( std::array< uint_t, 4 > size,
                              py::object & dtype,
                              uint_t ghostLayers,
                              Layout layout,
                              uint_t alignment)
{
   uint_t xSize = size[0];
   uint_t ySize = size[1];
   uint_t zSize = size[2];
   uint_t fSize = size[3];

   auto result = make_shared<py::none>();
   CreateFieldExporter exporter( xSize,ySize, zSize, fSize, ghostLayers, layout, dtype, alignment, result );
   python_coupling::for_each_noncopyable_type< FieldTypes... >  ( exporter );

   return *result;
}

//===================================================================================================================
//
//  createVTKWriter
//
//===================================================================================================================

class CreateVTKWriterExporter
{
 public:
   CreateVTKWriterExporter( const shared_ptr<StructuredBlockForest> & blocks,
                            ConstBlockDataID fieldId, const std::string & vtkName)
      : blocks_( blocks ), fieldId_(fieldId), vtkName_( vtkName )
   {}

   template< typename FieldType>
   void operator() ( python_coupling::NonCopyableWrap<FieldType> )
   {
      IBlock * firstBlock =  & ( * blocks_->begin() );
      if( firstBlock->isDataClassOrSubclassOf<FieldType>(fieldId_) )
         writer_ = shared_ptr<field::VTKWriter<FieldType> >( new field::VTKWriter<FieldType>(fieldId_, vtkName_));
   }

   shared_ptr< vtk::BlockCellDataWriterInterface > getCreatedWriter() {
      return writer_;
   }

 private:
   shared_ptr< vtk::BlockCellDataWriterInterface > writer_;
   shared_ptr< StructuredBlockStorage > blocks_;
   ConstBlockDataID fieldId_;
   std::string vtkName_;
};


template<typename... FieldTypes>
inline shared_ptr<vtk::BlockCellDataWriterInterface> createVTKWriter(const shared_ptr<StructuredBlockForest> & blocks,
                                                                     const std::string & name,
                                                                     const std::string & nameInVtkOutput = "")
{
   std::string vtkName = nameInVtkOutput;
   if( vtkName.empty())
      vtkName = name;

   if ( blocks->begin() == blocks->end() )
      return shared_ptr<vtk::BlockCellDataWriterInterface>();
   auto fieldID = python_coupling::blockDataIDFromString( *blocks, name );

   CreateVTKWriterExporter exporter(blocks, fieldID, vtkName);
   python_coupling::for_each_noncopyable_type< FieldTypes... >  ( std::ref(exporter) );
   if ( ! exporter.getCreatedWriter() ) {
      throw py::value_error("Failed to create VTK writer");
   }
   else {
      return exporter.getCreatedWriter();
   }
}


//===================================================================================================================
//
//  createBinarizationFieldWriter
//
//===================================================================================================================

class CreateBinarizationVTKWriterExporter
{
 public:
   CreateBinarizationVTKWriterExporter( const shared_ptr<StructuredBlockStorage> & blocks,
                                        ConstBlockDataID fieldId, const std::string & vtkName, uint_t mask)
      : blocks_( blocks ), fieldId_(fieldId), vtkName_( vtkName ), mask_(mask)
   {}

   template< typename FieldType>
   void operator() ( python_coupling::NonCopyableWrap<FieldType> )
   {
      IBlock * firstBlock =  & ( * blocks_->begin() );
      if( firstBlock->isDataClassOrSubclassOf<FieldType>(fieldId_) )
      {
         typedef field::BinarizationFieldWriter< FieldType > Writer;
         writer_ = shared_ptr< Writer >(new Writer(fieldId_, vtkName_, static_cast< typename FieldType::value_type >(mask_)));
      }
   }

   shared_ptr< vtk::BlockCellDataWriterInterface > getCreatedWriter() {
      return writer_;
   }

 private:
   shared_ptr< vtk::BlockCellDataWriterInterface > writer_;
   shared_ptr< StructuredBlockStorage > blocks_;
   ConstBlockDataID fieldId_;
   std::string vtkName_;
   uint_t mask_;
};


template<typename... FieldTypes>
inline shared_ptr<vtk::BlockCellDataWriterInterface> createBinarizationVTKWriter(const shared_ptr<StructuredBlockStorage> & blocks,
                                                                                 const std::string & name,
                                                                                 uint_t mask,
                                                                                 const std::string & nameInVtkOutput = "")
{
   std::string vtkName = nameInVtkOutput;
   if( vtkName.empty())
      vtkName = name;

   if ( blocks->begin() == blocks->end() )
      return shared_ptr<vtk::BlockCellDataWriterInterface>();
   auto fieldID = python_coupling::blockDataIDFromString( *blocks, name );

   CreateBinarizationVTKWriterExporter exporter(blocks, fieldID, vtkName, mask);
   python_coupling::for_each_noncopyable_type< FieldTypes... >  ( std::ref(exporter) );
   if ( ! exporter.getCreatedWriter() ) {
      throw py::value_error("Failed to create binarization field writer");
   }
   else {
      return exporter.getCreatedWriter();
   }
}

} // namespace internal

namespace py = pybind11;
template< typename... FieldTypes >
void exportFields(py::module_& m)
{
   using namespace py;

   py::module_ m2 = m.def_submodule("field", "Field Extension of the waLBerla python bindings");

   py::enum_< Layout >(m2, "Layout").value("fzyx", fzyx).value("zyxf", zyxf).export_values();

   python_coupling::for_each_noncopyable_type< FieldTypes... >(internal::FieldExporter(m2));
   python_coupling::for_each_noncopyable_type< real_t, int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t >(internal::FieldAllocatorExporter(m2));

   m2.def(
      "createField",
      [](std::array< uint_t, 4 > size, py::object & dtype, uint_t ghostLayers, Layout layout, uint_t alignment) {
        return internal::createPythonField< FieldTypes... >(size, dtype, ghostLayers, layout, alignment);
      },
      "size"_a, "dtype"_a, "ghostLayers"_a = uint_t(1), "layout"_a = zyxf, "alignment"_a = 0);

   m2.def(
      "addToStorage",
      [](const shared_ptr< StructuredBlockForest > & blocks, const std::string & name, py::object &dtype, uint_t fSize,
         Layout layout, uint_t ghostLayers, real_t initValue, uint_t alignment) {
         return internal::addToStorage< FieldTypes... >(blocks, name, dtype, fSize, ghostLayers, layout, initValue, alignment);
      },
      "blocks"_a, "name"_a, "dtype"_a, "fSize"_a = 1, "layout"_a = zyxf, "ghostLayers"_a = uint_t(1), "initValue"_a = 0.0, "alignment"_a = 0);

   m2.def( "createVTKWriter",
           [](const shared_ptr<StructuredBlockForest> & blocks, const std::string & name,
              const std::string & nameInVtkOutput = ""){
              return internal::createVTKWriter< FieldTypes... >(blocks, name, nameInVtkOutput);
           },
       "blocks"_a, "name"_a, "nameInVtkOutput"_a="" );


   #define UintFields Field<uint8_t,1 >, Field<uint16_t, 1>, Field<uint32_t, 1>, Field<uint64_t, 1>
   m2.def( "createBinarizationVTKWriter",
           [](const shared_ptr<StructuredBlockForest> & blocks, const std::string & name,
              uint_t mask, const std::string & nameInVtkOutput = ""){
             return internal::createBinarizationVTKWriter< UintFields >(blocks, name, mask, nameInVtkOutput);
           },
           "blocks"_a, "name"_a, "mask"_a, "nameInVtkOutput"_a="" );


}

} // namespace field
} // namespace walberla
