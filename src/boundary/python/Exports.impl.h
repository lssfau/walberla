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
//! \file Exports.impl.h
//! \ingroup boundary
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "python_coupling/Manager.h"
#include "python_coupling/helper/SliceToCellInterval.h"
#include "python_coupling/helper/ConfigFromDict.h"
#include "boundary/Boundary.h"
#include "field/FlagField.h"


namespace walberla {
namespace boundary {

   namespace internal
   {
      using python_coupling::localPythonSliceToCellInterval;

      template<typename BH>
      shared_ptr<BoundaryConfiguration> boundaryConfFromDict( BH & h, const BoundaryUID & uid, boost::python::dict d ) {
         shared_ptr<Config> cfg = python_coupling::configFromPythonDict( d );
         return h.createBoundaryConfiguration( uid, cfg->getGlobalBlock() );
      }
      template<typename BH>
      void BH_setBoundary1( BH & h, const std::string & name, cell_idx_t x, cell_idx_t y , cell_idx_t z, boost::python::dict conf ) {
         h.setBoundary( name,x,y,z,  *boundaryConfFromDict( h, name, conf) );
      }
      template<typename BH>
      void BH_setBoundary2( BH & h, const std::string & name, const boost::python::tuple & index, boost::python::dict conf ) {
         h.setBoundary( name,
                        localPythonSliceToCellInterval( *h.getFlagField() , index ),
                        *boundaryConfFromDict( h, name, conf) );
      }
      template<typename BH>
      void BH_forceBoundary1( BH & h, const std::string & name, cell_idx_t x, cell_idx_t y , cell_idx_t z, boost::python::dict conf ) {
         h.forceBoundary( name,x,y,z,  *boundaryConfFromDict( h, name, conf) );
      }
      template<typename BH>
      void BH_forceBoundary2( BH & h, const std::string & name, const boost::python::tuple & index, boost::python::dict conf )  {
         h.forceBoundary( name,
                        localPythonSliceToCellInterval( *h.getFlagField() , index ),
                        *boundaryConfFromDict( h, name, conf) );
      }

      template<typename BH>
      void BH_forceBoundary3( BH & h, const GhostLayerField<int,1> & indexField , boost::python::dict boundaryInfo )
      {
         using namespace boost::python;
         list keys = boundaryInfo.keys();

         std::map<int, FlagUID > flagUIDs;
         std::map<int, shared_ptr<BoundaryConfiguration> > boundaryConfigs;

         for (int i = 0; i < len( keys ); ++i)
         {
            int key =  extract<int>( keys[i] );
            extract<std::string> extracted_str_val  ( boundaryInfo[key] );
            extract<dict       > extracted_dict_val ( boundaryInfo[key] );

            if ( extracted_str_val.check() )
            {
               std::string boundaryName = extracted_str_val;
               flagUIDs[key] = FlagUID ( boundaryName );
            }
            else if ( extracted_dict_val.check() )
            {
               dict info = extracted_dict_val;
               std::string boundaryName = extract<std::string>( info["name"] );

               dict configDict = extract<dict>( info["config"] );

               flagUIDs[key] = FlagUID ( boundaryName );
               boundaryConfigs[key] = boundaryConfFromDict( h, boundaryName, configDict);
            }
            else {
               PyErr_SetString( PyExc_ValueError, "Invalid parameter");
               throw error_already_set();
            }
         }

         if ( indexField.xyzSize() != h.getFlagField()->xyzSize() || indexField.nrOfGhostLayers() > h.getFlagField()->nrOfGhostLayers() ) {
            PyErr_SetString( PyExc_ValueError, "Index field has to have same size as flag field");
            throw error_already_set();
         }

         // iterate over flag field
         cell_idx_t gl = cell_idx_c( indexField.nrOfGhostLayers() );
         for( cell_idx_t z = -gl; z < cell_idx_c( indexField.zSize() ) + gl; ++z )
            for( cell_idx_t y = -gl; y < cell_idx_c( indexField.ySize() ) + gl; ++y )
               for( cell_idx_t x = -gl; x < cell_idx_c( indexField.xSize() ) + gl; ++x )
               {
                  int index = indexField(x,y,z);
                  if ( flagUIDs.find( index ) != flagUIDs.end() )
                  {
                     if ( boundaryConfigs.find( index ) != boundaryConfigs.end()  )
                        h.forceBoundary( flagUIDs[index],x,y,z, * boundaryConfigs[index] );
                     else
                        h.forceBoundary( flagUIDs[index],x,y,z );
                  }
               }
      }

      template<typename BH>
      void BH_setDomainSlice( BH & h, const boost::python::tuple & index ) {
         h.setDomain( localPythonSliceToCellInterval( *h.getFlagField() , index ) );
      }

      template<typename BH>
      void BH_forceDomainSlice( BH & h, const boost::python::tuple & index ) {
         h.forceDomain( localPythonSliceToCellInterval( *h.getFlagField() , index ) );
      }

      template<typename BH>
      void BH_fillDomainSlice( BH & h, const boost::python::tuple & index ) {
         h.fillWithDomain( localPythonSliceToCellInterval( *h.getFlagField() , index ) );
      }
      template<typename BH>
      void BH_removeDomainSlice( BH & h, const boost::python::tuple & index ) {
         h.removeDomain( localPythonSliceToCellInterval( *h.getFlagField() , index ) );
      }
      template<typename BH>
      void BH_removeBoundarySlice( BH & h, const boost::python::tuple & index ) {
         h.removeBoundary( localPythonSliceToCellInterval( *h.getFlagField() , index ) );
      }
      template<typename BH>
      void BH_clearSlice( BH & h, const boost::python::tuple & index ) {
         h.clear( localPythonSliceToCellInterval( *h.getFlagField() , index ) );
      }



      struct BoundaryHandlingExporter
      {
         template< typename BH>
         void operator() ( python_coupling::NonCopyableWrap<BH> )
         {
            using namespace boost::python;
            void ( BH::*p_exe1 )( uint_t )= &BH::operator();

            bool ( BH::*p_isEmpty        ) ( cell_idx_t, cell_idx_t, cell_idx_t ) const = &BH::isEmpty       ;
            bool ( BH::*p_isNearBoundary ) ( cell_idx_t, cell_idx_t, cell_idx_t ) const = &BH::isNearBoundary;
            bool ( BH::*p_isBoundary     ) ( cell_idx_t, cell_idx_t, cell_idx_t ) const = &BH::isBoundary    ;
            bool ( BH::*p_isDomain       ) ( cell_idx_t, cell_idx_t, cell_idx_t ) const = &BH::isDomain      ;

            void ( BH::*p_setDomain      ) ( cell_idx_t, cell_idx_t, cell_idx_t )= &BH::setDomain  ;
            void ( BH::*p_forceDomain    ) ( cell_idx_t, cell_idx_t, cell_idx_t )= &BH::forceDomain;
            void ( BH::*p_fillWithDomain1) ( const uint_t )                      = &BH::fillWithDomain;
            void ( BH::*p_fillWithDomain2) ( cell_idx_t, cell_idx_t, cell_idx_t )= &BH::fillWithDomain;

            void ( BH::*p_removeDomain1 ) ( const uint_t )                      = &BH::removeDomain;
            void ( BH::*p_removeDomain2 ) ( cell_idx_t, cell_idx_t, cell_idx_t )= &BH::removeDomain;

            void ( BH::*p_removeBoundary1 ) ( const uint_t )                      = &BH::removeBoundary;
            void ( BH::*p_removeBoundary2 ) ( cell_idx_t, cell_idx_t, cell_idx_t )= &BH::removeBoundary;

            void ( BH::*p_clear1 ) ( const uint_t )                      = &BH::clear;
            void ( BH::*p_clear2 ) ( cell_idx_t, cell_idx_t, cell_idx_t )= &BH::clear;

            typename BH::FlagField * ( BH::*p_getFlagField) () = &BH::getFlagField;

            typename BH::flag_t ( BH::*p_getNearBoundaryFlag ) () const  = &BH::getNearBoundaryFlag;
            typename BH::flag_t ( BH::*p_getBoundaryMask )     () const  = &BH::getBoundaryMask;
            typename BH::flag_t ( BH::*p_getDomainMask )       () const  = &BH::getDomainMask;

            class_< BH, boost::noncopyable > ( "BoundaryHandling", no_init )
               .def( "__call__",        p_exe1,                    ( arg("numberOfGhostLayersToInclude")=0 )  )
               .def( "isEmpty"       ,  p_isEmpty,                 ( arg("x"), arg("y"), arg("z") ) )
               .def( "isNearBoundary",  p_isNearBoundary,          ( arg("x"), arg("y"), arg("z") ) )
               .def( "isBoundary"    ,  p_isBoundary,              ( arg("x"), arg("y"), arg("z") ) )
               .def( "isDomain"      ,  p_isDomain,                ( arg("x"), arg("y"), arg("z") ) )
               .def( "setDomain"     ,  p_setDomain,               ( arg("x"), arg("y"), arg("z") ) )
               .def( "setDomain"     ,  BH_forceDomainSlice<BH>,   ( arg("slice") ) )
               .def( "forceDomain"   ,  p_forceDomain,             ( arg("x"), arg("y"), arg("z") ) )
               .def( "forceDomain"   ,  BH_forceDomainSlice<BH>,   ( arg("slice") ) )
               .def( "fillWithDomain",  p_fillWithDomain1,         ( arg("numberOfGhostLayersToInclude")=0 ) )
               .def( "fillWithDomain",  p_fillWithDomain2,         ( arg("x"), arg("y"), arg("z") ) )
               .def( "fillWithDomain",  BH_fillDomainSlice<BH>,    ( arg("slice") ) )
               .def( "setBoundary",     &BH_setBoundary1<BH>,      ( arg("name"), arg("x"), arg("y"), arg("z"), arg("boundaryParams")=dict() ) )
               .def( "setBoundary",     &BH_setBoundary2<BH>,      ( arg("name"), arg("slice"), arg("boundaryParams")=dict() ) )
               .def( "forceBoundary",   &BH_forceBoundary1<BH>,    ( arg("name"), arg("x"), arg("y"), arg("z"), arg("boundaryParams")=dict() ) )
               .def( "forceBoundary",   &BH_forceBoundary2<BH>,    ( arg("name"), arg("slice"), arg("boundaryParams")=dict() ) )
               .def( "forceBoundary",   &BH_forceBoundary3<BH>,    ( arg("indexField"), arg("boundaryInfo") ) )
               .def( "removeDomain",    p_removeDomain1,           ( arg("numberOfGhostLayersToInclude")=0 ) )
               .def( "removeDomain",    p_removeDomain2,           ( arg("x"), arg("y"), arg("z") ) )
               .def( "removeDomain",    BH_removeDomainSlice<BH>,  ( arg("slice") ) )
               .def( "removeBoundary",  p_removeBoundary1,         ( arg("numberOfGhostLayersToInclude")=0 ) )
               .def( "removeBoundary",  p_removeBoundary2,         ( arg("x"), arg("y"), arg("z") ) )
               .def( "removeBoundary",  BH_removeBoundarySlice<BH>,( arg("slice") ) )
               .def( "clear",           p_clear1,                  ( arg("numberOfGhostLayersToInclude")=0 ) )
               .def( "clear",           p_clear2,                  ( arg("x"), arg("y"), arg("z") ) )
               .def( "clear",           BH_clearSlice<BH>,         ( arg("slice") ) )
               .def( "getFlagField",    p_getFlagField,              return_internal_reference<>() )
               .def( "getNearBoundaryFlag", p_getNearBoundaryFlag )
               .def( "getBoundaryMask",     p_getBoundaryMask     )
               .def( "getDomainMask",       p_getDomainMask       )
            ;
         }
      };
   } // namespace internal

   template<typename BoundaryHandlings>
   void exportModuleToPython()
   {
      python_coupling::for_each_noncopyable_type< BoundaryHandlings > ( internal::BoundaryHandlingExporter() );

      auto pythonManager = python_coupling::Manager::instance();
      pythonManager->addBlockDataConversion< BoundaryHandlings> ();

   }


} // namespace boundary
} // namespace walberla


