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
//! \ingroup blockforest
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/communication/UniformDirectScheme.h"

#include "python_coupling/helper/MplHelpers.h"

#include <pybind11/pybind11.h>

namespace walberla
{
namespace blockforest
{
namespace py = pybind11;
namespace internal
{
//===================================================================================================================
//
//  UniformBufferedScheme
//
//===================================================================================================================

/// the purpose of this class could also be solved by adding return_internal_reference to "createUniformDirectScheme"
/// however this is not easily possible since it returns not a reference but an py::object
template< typename Stencil >
class UniformBufferedSchemeWrapper : public blockforest::communication::UniformBufferedScheme< Stencil >
{
 public:
   UniformBufferedSchemeWrapper(const shared_ptr< StructuredBlockForest >& bf, const int tag)
      : blockforest::communication::UniformBufferedScheme< Stencil >(bf, tag), blockforest_(bf)
   {}

 private:
   shared_ptr< StructuredBlockForest > blockforest_;
};

struct UniformBufferedSchemeExporter
{
   UniformBufferedSchemeExporter(py::module_& m) : m_(m) {}
   template< typename Stencil >
   void operator()(python_coupling::NonCopyableWrap< Stencil >) const
   {
      typedef UniformBufferedSchemeWrapper< Stencil > UBS;
      std::string class_name = "UniformBufferedScheme" + std::string(Stencil::NAME);

      py::class_< UBS, shared_ptr< UBS > >(m_, class_name.c_str())
         .def("__call__", &UBS::operator())
         .def("communicate", &UBS::communicate)
         .def("startCommunication", &UBS::startCommunication)
         .def("wait", &UBS::wait)
         .def("addPackInfo", &UBS::addPackInfo)
         .def("addDataToCommunicate", &UBS::addDataToCommunicate)
         .def("localMode", &UBS::localMode)
         .def("setLocalMode", &UBS::setLocalMode);
   }
   const py::module_& m_;
};

class UniformBufferedSchemeCreator
{
 public:
   UniformBufferedSchemeCreator( const shared_ptr<StructuredBlockForest> & bf,
                               const std::string & stencilName,
                               const int tag )
      : blockforest_( bf), stencilName_( stencilName ), tag_( tag )
   {}

   template<typename Stencil>
   void operator() ( python_coupling::NonCopyableWrap<Stencil> )
   {

      if ( std::string(Stencil::NAME) == stencilName_ ) {
         result_ = py::cast( make_shared< UniformBufferedSchemeWrapper<Stencil> > ( blockforest_, tag_ ) );
      }
   }

   py::object getResult() { return result_; }
 private:
   py::object result_;
   shared_ptr<StructuredBlockForest> blockforest_;
   std::string stencilName_;
   const int tag_;
};


template<typename... Stencils>
py::object createUniformBufferedScheme( const shared_ptr<StructuredBlockForest> & bf,
                                        const std::string & stencil, const int tag )
{
   UniformBufferedSchemeCreator creator( bf, stencil, tag );
   python_coupling::for_each_noncopyable_type< Stencils... >  ( std::ref(creator) );

   if ( !creator.getResult() )
   {
      throw py::value_error("Unknown stencil.");
   }
   return creator.getResult();
}

//===================================================================================================================
//
//  UniformDirectScheme
//
//===================================================================================================================

template< typename Stencil >
class UniformDirectSchemeWrapper : public blockforest::communication::UniformDirectScheme< Stencil >
{
 public:
   UniformDirectSchemeWrapper(const shared_ptr< StructuredBlockForest >& bf, const int tag)
      : blockforest::communication::UniformDirectScheme< Stencil >(
           bf, shared_ptr< walberla::communication::UniformMPIDatatypeInfo >(), tag),
        blockforest_(bf)
   {}

 private:
   shared_ptr< StructuredBlockForest > blockforest_;
};

struct UniformDirectSchemeExporter
{
   UniformDirectSchemeExporter(py::module_& m) : m_(m) {}
   template< typename Stencil >
   void operator()(python_coupling::NonCopyableWrap< Stencil >) const
   {
      typedef UniformDirectSchemeWrapper< Stencil > UDS;
      std::string class_name = "UniformDirectScheme_" + std::string(Stencil::NAME);

      py::class_< UDS, shared_ptr<UDS> >(m_, class_name.c_str() )
         .def("__call__", &UDS::operator())
         .def("communicate", &UDS::communicate)
         .def("startCommunication", &UDS::startCommunication)
         .def("wait", &UDS::wait)
         .def("addDataToCommunicate", &UDS::addDataToCommunicate);
   }
   const py::module_ m_;
};

class UniformDirectSchemeCreator
{
 public:
   UniformDirectSchemeCreator( const shared_ptr<StructuredBlockForest> & bf,
                               const std::string & stencilName,
                               const int tag )
      : blockforest_( bf), stencilName_( stencilName ), tag_( tag )
   {}

   template<typename Stencil>
   void operator() ( python_coupling::NonCopyableWrap<Stencil> )
   {

      if ( std::string(Stencil::NAME) == stencilName_ ) {
         result_ = py::cast( make_shared< UniformDirectSchemeWrapper<Stencil> > ( blockforest_, tag_ ) );
      }
   }

   py::object getResult() { return result_; }
 private:
   py::object result_;
   shared_ptr<StructuredBlockForest> blockforest_;
   std::string stencilName_;
   const int tag_;
};


template<typename... Stencils>
py::object createUniformDirectScheme( const shared_ptr<StructuredBlockForest> & bf,
                                      const std::string & stencil, const int tag )
{
   UniformDirectSchemeCreator creator( bf, stencil, tag );
   python_coupling::for_each_noncopyable_type< Stencils... >  ( std::ref(creator) );

   if ( !creator.getResult() )
   {
      throw py::value_error("Unknown stencil.");
   }
   return creator.getResult();
}




} // namespace internal

template< typename... Stencils >
void exportUniformDirectScheme(py::module_& m)
{
   using namespace py;

   python_coupling::for_each_noncopyable_type< Stencils... >(internal::UniformDirectSchemeExporter(m));
   m.def( "createUniformDirectScheme",
           [](const shared_ptr<StructuredBlockForest> & blocks, const std::string & stencil, const int tag)
            {
             return internal::createUniformDirectScheme< Stencils... >(blocks, stencil, tag);
            },
           "blocks"_a, "stencil"_a, "tag"_a=778 );


}

template< typename... Stencils >
void exportUniformBufferedScheme(py::module_& m)
{
   using namespace py;

   py::enum_< LocalCommunicationMode >(m, "LocalCommunicationMode")
      .value("START", START)
      .value("WAIT", WAIT)
      .value("BUFFER", BUFFER)
      .export_values();

   python_coupling::for_each_noncopyable_type< Stencils... >(internal::UniformBufferedSchemeExporter(m));
   m.def( "createUniformBufferedScheme",
          [](const shared_ptr<StructuredBlockForest> & blocks, const std::string & stencil, const int tag)
          {
            return internal::createUniformBufferedScheme< Stencils... >(blocks, stencil, tag);
          },
          "blocks"_a, "stencil"_a, "tag"_a=778 );
}

} // namespace blockforest
} // namespace walberla