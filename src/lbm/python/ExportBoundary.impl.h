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
//! \file ExportBoundary.impl.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "python_coupling/helper/MplHelpers.h"
#include "python_coupling/helper/BoostPythonHelpers.h"
#include "python_coupling/helper/ModuleScope.h"
#include "python_coupling/helper/BlockStorageExportHelpers.h"
#include "python_coupling/helper/ConfigFromDict.h"
#include "python_coupling/helper/SliceToCellInterval.h"
#include "python_coupling/Manager.h"

#include "core/logging/Logging.h"

#include "boundary/python/Exports.h"

#include "lbm/boundary/factories/ExtendedBoundaryHandlingFactory.h"

#include "field/python/FieldExport.h"
#include "field/adaptors/AdaptorCreators.h"

#include <boost/mpl/transform.hpp>
#include <boost/mpl/copy.hpp>

namespace walberla {
namespace lbm {


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

      // iterate over flag field
      if ( indexField.xyzSize() != h.getFlagField()->xyzSize() || indexField.nrOfGhostLayers() > h.getFlagField()->nrOfGhostLayers() ) {
         WALBERLA_LOG_DEVEL("indexField " << indexField.xyzSizeWithGhostLayer() );
         WALBERLA_LOG_DEVEL("flagField  " << h.getFlagField()->xyzSizeWithGhostLayer() );
         PyErr_SetString( PyExc_ValueError, "Index field has to have same size as flag field");
         throw error_already_set();
      }


      cell_idx_t gl = cell_idx_c( indexField.nrOfGhostLayers() );
      for( cell_idx_t z = -gl; z < cell_idx_c( indexField.zSizeWithGhostLayer() ); ++z )
         for( cell_idx_t y = -gl; y < cell_idx_c( indexField.ySizeWithGhostLayer() ); ++y )
            for( cell_idx_t x = -gl; x < cell_idx_c( indexField.xSizeWithGhostLayer() ); ++x )
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




   class ExtendedBoundaryHandlingCreator
   {
   public:

      ExtendedBoundaryHandlingCreator( const shared_ptr<StructuredBlockStorage> & blocks,
                                       const std::string & name, BlockDataID pdfFieldID, BlockDataID flagFieldID )
         : blocks_( blocks ), name_( name ),
           pdfFieldID_( pdfFieldID ), flagFieldID_( flagFieldID ), found_ ( false )
      {}


      template<typename LatticeModel>
      void operator()( python_coupling::NonCopyableWrap<LatticeModel>  )
      {
         using namespace boost::python;

         IBlock & firstBlock = *(  blocks_->begin() );
         if ( firstBlock.isDataClassOrSubclassOf< PdfField<LatticeModel> >( pdfFieldID_ ) )
         {
            WALBERLA_ASSERT( !found_ );
            ExtendedBoundaryHandlingFactory<LatticeModel, FlagField<uint8_t> >::addBoundaryHandlingToStorage(
                     blocks_, name_, flagFieldID_, pdfFieldID_, FlagUID("fluid") );
            found_ = true;
         }
      }

      bool successful() const { return found_; }

   private:
      shared_ptr<StructuredBlockStorage> blocks_;
      std::string name_;
      BlockDataID pdfFieldID_;
      BlockDataID flagFieldID_;

      bool found_;
   };


   template<typename LatticeModels>
   void addBoundaryHandlingToStorage( const shared_ptr<StructuredBlockStorage> & bs,
                                      const std::string & name,
                                      const std::string & pdfFieldStringID,
                                      const std::string & flagFieldStringID  )
   {
      using namespace boost::python;

      if ( bs->begin() == bs->end() )
         return;


      IBlock & firstBlock = *(  bs->begin() );

      BlockDataID pdfFieldID  = python_coupling::blockDataIDFromString( firstBlock, pdfFieldStringID  );
      BlockDataID flagFieldID = python_coupling::blockDataIDFromString( firstBlock, flagFieldStringID );

      if ( ! firstBlock.isDataClassOrSubclassOf< FlagField<uint8_t> >( flagFieldID ) )
      {
         PyErr_SetString( PyExc_ValueError, "Unknown FlagField type. Please provide a FlagField with 8 bit per cell");
         throw error_already_set();
      }

      ExtendedBoundaryHandlingCreator creator( bs, name, pdfFieldID, flagFieldID );
      python_coupling::for_each_noncopyable_type<LatticeModels>( std::ref( creator ) );

      if ( ! creator.successful() )
      {
         PyErr_SetString( PyExc_ValueError, "No Boundary Handling found for this lattice model");
         throw error_already_set();
      }
   }

   template< typename LatticeModel >
   struct ExtendedBoundaryHandlingFromLatticeModel
   {
      typedef typename ExtendedBoundaryHandlingFactory< LatticeModel, FlagField<uint8_t> >::BoundaryHandling type;
   };
 }


template< typename LatticeModels >
struct ExtendedBoundaryHandlingsFromLatticeModels
{
   typedef typename boost::mpl::transform< LatticeModels, internal::ExtendedBoundaryHandlingFromLatticeModel<boost::mpl::_1> >::type type;
};


template<typename LatticeModels, typename FlagFields>
void exportBoundary()
{
   using namespace boost::python;

   python_coupling::ModuleScope scope("lbm");

   // Extended Boundary Handling
   boundary::exportModuleToPython< typename ExtendedBoundaryHandlingsFromLatticeModels< LatticeModels >::type >();

   def ( "addBoundaryHandlingToStorage", &internal::addBoundaryHandlingToStorage<LatticeModels>,
            ( arg("blocks"),
              arg("name"),
              arg("pdfField"),
              arg("flagField") ) );

   auto pythonManager = python_coupling::Manager::instance();
   pythonManager->addBlockDataConversion< typename ExtendedBoundaryHandlingsFromLatticeModels< LatticeModels >::type > ();
}



} // namespace lbm
} // namespace walberla


