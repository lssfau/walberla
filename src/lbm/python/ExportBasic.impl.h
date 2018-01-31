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
//! \file ExportBasic.impl.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "python_coupling/helper/MplHelpers.h"
#include "python_coupling/helper/BoostPythonHelpers.h"
#include "python_coupling/helper/ModuleScope.h"
#include "python_coupling/helper/BlockStorageExportHelpers.h"

#include "core/logging/Logging.h"

#include "boundary/python/Exports.h"
#include "field/adaptors/AdaptorCreators.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/ForceModel.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/Adaptors.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/SplitPureSweep.h"

#include "field/python/FieldExport.h"

#include <boost/mpl/transform.hpp>
#include <boost/mpl/copy.hpp>


namespace walberla {
namespace lbm {


namespace internal
{
   //===================================================================================================================
   //
   //  LatticeModel
   //
   //===================================================================================================================

   class LatticeModelCreator
   {
   public:
      LatticeModelCreator( bool compressible, uint_t eqOrder, const std::string & stencil,
                           boost::python::object collisionModel, boost::python::object forceModel )
         : compressible_( compressible ), equilibriumAccuracyOrder_( eqOrder ), stencil_( stencil ),
           collisionModel_( collisionModel ), forceModel_( forceModel )
      {}

      template<typename LatticeModel_T>
      void operator() ( python_coupling::NonCopyableWrap<LatticeModel_T> )
      {
         using namespace boost::python;
         using namespace collision_model;

         typedef typename LatticeModel_T::Stencil Stencil_T;

         if ( compressible_              != LatticeModel_T::compressible             ) return;
         if ( equilibriumAccuracyOrder_  != LatticeModel_T::equilibriumAccuracyOrder ) return;
         if ( stencil_                   != Stencil_T::NAME                          ) return;

         if ( ! extract< typename LatticeModel_T::CollisionModel>(collisionModel_).check() ) return;
         if ( ! extract< typename LatticeModel_T::ForceModel    >(forceModel_    ).check() ) return;

         result_ = object( LatticeModel_T( extract< typename LatticeModel_T::CollisionModel>(collisionModel_),
                                           extract< typename LatticeModel_T::ForceModel    >(forceModel_    ) ) );
      }

      boost::python::object getResult() const { return result_; }

   private:
      bool compressible_;
      uint_t equilibriumAccuracyOrder_;
      std::string stencil_;
      boost::python::object collisionModel_;
      boost::python::object forceModel_;


      boost::python::object result_;
   };

   struct LatticeModelExporter
   {
      template<typename LatticeModel_T>
      static boost::python::list getDirections( LatticeModel_T & )
      {
         using namespace boost::python;
         list directionList;
         const int dimension = LatticeModel_T::Stencil::D;
         if( dimension < 1 || dimension > 3) {
            PyErr_SetString( PyExc_ValueError, "Only stencils with dimensions 1,2,3 supported");
            throw error_already_set();
         }
         for( auto it = LatticeModel_T::Stencil::begin(); it != LatticeModel_T::Stencil::end(); ++it )
         {
            if(dimension == 1) {
               directionList.append( make_tuple(it.cx() ) );
            }
            else if( dimension == 2) {
               directionList.append( make_tuple(it.cx(), it.cy()) );
            } else if (dimension == 3) {
               directionList.append( make_tuple(it.cx(), it.cy(), it.cz() ) );
            }
         }
         return directionList;
      }
      template<typename LatticeModel_T>
      static bool isCompressible(LatticeModel_T &) {
         return LatticeModel_T::compressible;
      }
      template<typename LatticeModel_T>
      static int equilibriumAccuracyOrder(LatticeModel_T &) {
         return LatticeModel_T::equilibriumAccuracyOrder;
      }
      template<typename LatticeModel_T>
      static boost::python::str stencilName(LatticeModel_T &) {
        return boost::python::str(LatticeModel_T::Stencil::NAME);
      }
      template<typename LatticeModel_T>
      static boost::python::str communicationStencilName(LatticeModel_T &) {
        return boost::python::str(LatticeModel_T::CommunicationStencil::NAME);
      }

      template<typename LatticeModel_T>
      void operator()( python_coupling::NonCopyableWrap<LatticeModel_T>  )
      {
         using namespace boost::python;

         typename LatticeModel_T::CollisionModel& ( LatticeModel_T::*p_collisionModel) () = &LatticeModel_T::collisionModel;
         typename LatticeModel_T::ForceModel    & ( LatticeModel_T::*p_forceModel)     () = &LatticeModel_T::forceModel;

         class_< LatticeModel_T > ( LatticeModel_T::NAME, no_init )
                  .add_property( "collisionModel", make_function( p_collisionModel, return_internal_reference<>() ) )
                  .add_property( "forceModel",     make_function( p_forceModel    , return_internal_reference<>() ) )
                  .add_property( "compressible",            &LatticeModelExporter::isCompressible<LatticeModel_T>)
                  .add_property( "equilibriumAccuracyOrder",&LatticeModelExporter::equilibriumAccuracyOrder<LatticeModel_T>)
                  .add_property( "stencilName",             &LatticeModelExporter::stencilName<LatticeModel_T>)
                  .add_property( "communicationStencilName",&LatticeModelExporter::communicationStencilName<LatticeModel_T>)
                  .add_property( "directions",              &LatticeModelExporter::getDirections<LatticeModel_T> )
         ;
      }
   };


   template<typename LatticeModels>
   boost::python::object makeLatticeModel( const std::string & stencil, boost::python::object collisionModel, boost::python::object forceModel,
                                           bool compressible, uint_t equilibriumAccuracyOrder )
   {
      using namespace boost::python;

      LatticeModelCreator creator( compressible, equilibriumAccuracyOrder, stencil, collisionModel, forceModel );
      python_coupling::for_each_noncopyable_type<LatticeModels>( boost::ref(creator) );

      if ( creator.getResult() == object() )
      {
         PyErr_SetString( PyExc_ValueError, "No such LatticeModel available.");
         throw error_already_set();
      }

      return creator.getResult();
   }

   //===================================================================================================================
   //
   //  PdfField
   //
   //===================================================================================================================

   template<typename LatticeModel_T>
   typename PdfField<LatticeModel_T>::iterator pythonSliceToFieldIterator( PdfField<LatticeModel_T> & field,
                                                                           boost::python::tuple pyIndex )
   {
      using python_coupling::localPythonSliceToCellInterval;
      return field.beginSliceXYZ( localPythonSliceToCellInterval(field, pyIndex) );
   }

   template<typename LatticeModel_T>
   void pdfField_setDensityAndVelocity( PdfField<LatticeModel_T> & field, boost::python::tuple pyIndex,
                                        const Vector3<real_t> & velocity, real_t rho )
   {
      typename PdfField<LatticeModel_T>::iterator beginIterator = pythonSliceToFieldIterator<LatticeModel_T>( field, pyIndex );
      DensityAndVelocityRange< LatticeModel_T, typename  PdfField<LatticeModel_T>::iterator >::set( beginIterator, field.end(), field.latticeModel(), velocity, rho );
   }

   template<typename LatticeModel_T>
   void pdfField_setToEquilibrium( PdfField<LatticeModel_T> & field, boost::python::tuple pyIndex,
                                   const Vector3<real_t> & velocity, real_t rho )
   {
      typename PdfField<LatticeModel_T>::iterator beginIterator = pythonSliceToFieldIterator<LatticeModel_T>( field, pyIndex );
      EquilibriumRange< LatticeModel_T, typename  PdfField<LatticeModel_T>::iterator >::set( beginIterator, field.end(), velocity, rho );
   }

   template<typename LatticeModel_T>
   boost::python::list pdfField_getPressureTensor( PdfField<LatticeModel_T> & field, cell_idx_t x, cell_idx_t y, cell_idx_t z )
   {
      using namespace boost::python;

      Matrix3<real_t> m = field.getPressureTensor( x,y,z );
      list result;
      for(uint_t i=0; i<3; ++i )
      {
         list row;
         for(uint_t j=0; j<3; ++j)
            row.append( m(i,j) );
         result.append(row);
      }
      return result;
   }


   struct PdfFieldExporter
   {
      template<typename LatticeModel_T>
      void operator()( python_coupling::NonCopyableWrap<LatticeModel_T>  )
      {
         using namespace boost::python;

         typedef PdfField<LatticeModel_T> PdfField_T;
         typedef GhostLayerField<real_t, LatticeModel_T::Stencil::Size >  Base;


         LatticeModel_T & ( PdfField_T::*ptr_latticeModel ) ()  = &PdfField_T::latticeModel;

         //real_t ( PdfField_T::*ptr_getShearRate                  )( cell_idx_t, cell_idx_t, cell_idx_t         ) const = &PdfField_T::getShearRate;

         real_t ( PdfField_T::*ptr_getDensity                    )( cell_idx_t, cell_idx_t, cell_idx_t         ) const = &PdfField_T::getDensity;
         real_t ( PdfField_T::*ptr_getDensitySI                  )( cell_idx_t, cell_idx_t, cell_idx_t, real_t ) const = &PdfField_T::getDensitySI;

         Vector3<real_t> ( PdfField_T::*ptr_getMomentumDensity            )( cell_idx_t, cell_idx_t, cell_idx_t         ) const = &PdfField_T::getMomentumDensity;
         Vector3<real_t> ( PdfField_T::*ptr_getEquilibriumMomentumDensity )( cell_idx_t, cell_idx_t, cell_idx_t         ) const = &PdfField_T::getEquilibriumMomentumDensity;
         real_t          ( PdfField_T::*ptr_getShearRate                  )( cell_idx_t, cell_idx_t, cell_idx_t                 ) const  = &PdfField_T::getShearRate;
         Vector3<real_t> ( PdfField_T::*ptr_getVelocity                   )( cell_idx_t, cell_idx_t, cell_idx_t                 ) const  = &PdfField_T::getVelocity;
         Vector3<real_t> ( PdfField_T::*ptr_getVelocitySI                 )( cell_idx_t, cell_idx_t, cell_idx_t, real_t, real_t ) const  = &PdfField_T::getVelocitySI;
         Vector3<real_t> ( PdfField_T::*ptr_getEquilibriumVelocity        )( cell_idx_t, cell_idx_t, cell_idx_t                 ) const  = &PdfField_T::getEquilibriumVelocity;

         try
         {
            class_< PdfField_T, shared_ptr<PdfField_T>, bases< Base >, boost::noncopyable >("PdfField", no_init)
               .add_property("latticeModel", make_function( ptr_latticeModel, return_internal_reference<>() ) )
               .def( "setDensityAndVelocity", pdfField_setDensityAndVelocity<LatticeModel_T>, ( args("slice"), args("velocity"), args("density") ) )
               .def( "setToEquilibrium",      pdfField_setToEquilibrium<LatticeModel_T>,      ( args("slice"), args("velocity"), args("density") ) )
               .def( "getShearRate",                   ptr_getShearRate                           , ( arg("x"), arg("y"), arg("z")  ) )
               .def( "getDensity",                     ptr_getDensity                             , ( arg("x"), arg("y"), arg("z") ) )
               .def( "getDensitySI",                   ptr_getDensitySI                           , ( arg("x"), arg("y"), arg("z"), arg("rho_SI") ) )
               .def( "getMomentumDensity"            , ptr_getMomentumDensity                     , ( arg("x"), arg("y"), arg("z") ) )
               .def( "getEquilibriumMomentumDensity" , ptr_getEquilibriumMomentumDensity          , ( arg("x"), arg("y"), arg("z") ) )
               .def( "getVelocity"                   , ptr_getVelocity                            , ( arg("x"), arg("y"), arg("z") ) )
               .def( "getVelocitySI"                 , ptr_getVelocitySI                          , ( arg("x"), arg("y"), arg("z"), arg("dx_SI"), arg("dy_SI") ) )
               .def( "getEquilibriumVelocity"        , ptr_getEquilibriumVelocity                 , ( arg("x"), arg("y"), arg("z") ) )
               .def( "getPressureTensor",              pdfField_getPressureTensor<LatticeModel_T> , ( arg("x"), arg("y"), arg("z") ) )
            ;
         }
         catch (...) {
            WALBERLA_LOG_WARNING( "Exporting PDFField failed. Did you forget to export the corresponding GhostLayerField?" );
         }
      }
   };

   class AddPdfFieldToStorageExporter
   {
   public:
      AddPdfFieldToStorageExporter( const shared_ptr<StructuredBlockStorage> & blocks,
                           const std::string & name,
                           boost::python::object latticeModel,
                           const Vector3<real_t> & initialVelocity, real_t initialDensity,
                           uint_t gl, field::Layout layout,
                           const std::string & densityAdaptorName, const std::string & velocityAdaptorName,
                           const std::string & shearRateAdaptorName )
         : blocks_( blocks ), name_( name ), latticeModel_(latticeModel),
           initialVelocity_( initialVelocity ), initialDensity_( initialDensity ),
           gl_(gl),layout_( layout),
           densityAdaptorName_( densityAdaptorName ), velocityAdaptorName_( velocityAdaptorName ),
           shearRateAdaptorName_( shearRateAdaptorName),
           found_(false)
      {}

      template< typename LatticeModel_T>
      void operator() ( python_coupling::NonCopyableWrap<LatticeModel_T> )
      {
         using namespace boost::python;

         if( ! extract<LatticeModel_T>( latticeModel_ ).check() )
            return;

         WALBERLA_ASSERT( !found_ );
         BlockDataID pdfFieldID = addPdfFieldToStorage<LatticeModel_T>( blocks_, name_, extract<LatticeModel_T>( latticeModel_ ),
                                                                        initialVelocity_, initialDensity_, gl_, layout_ );

         if ( densityAdaptorName_.size() > 0 ) {
            field::addFieldAdaptor< typename lbm::Adaptor<LatticeModel_T>::Density > ( blocks_, pdfFieldID, densityAdaptorName_ );
         }
         if ( velocityAdaptorName_.size() > 0 ) {
            field::addFieldAdaptor< typename lbm::Adaptor<LatticeModel_T>::Velocity >( blocks_, pdfFieldID, velocityAdaptorName_ );
         }
         if ( shearRateAdaptorName_.size() > 0 ) {
            field::addFieldAdaptor< typename lbm::Adaptor<LatticeModel_T>::ShearRate >( blocks_, pdfFieldID, shearRateAdaptorName_ );
         }
         found_ = true;
      }

      bool successful() { return found_; }

   private:
      shared_ptr< StructuredBlockStorage > blocks_;
      std::string name_;
      boost::python::object latticeModel_;
      Vector3<real_t> initialVelocity_;

      real_t initialDensity_;
      uint_t gl_;
      field::Layout layout_;

      std::string densityAdaptorName_;
      std::string velocityAdaptorName_;
      std::string shearRateAdaptorName_;

      bool found_;
   };

   template<typename LatticeModels>
   void addPdfFieldToStorage( const shared_ptr<StructuredBlockStorage> & blocks,
                              const std::string & identifier,
                              boost::python::object latticeModel,
                              const Vector3<real_t> & initialVelocity, real_t initialDensity,
                              uint_t gl, field::Layout layout,
                              const std::string & densityAdaptor, const std::string & velocityAdaptor,
                              const std::string & shearRateAdaptor )
   {
      using namespace boost::python;

      internal::AddPdfFieldToStorageExporter exporter( blocks, identifier, latticeModel, initialVelocity,
                                                       initialDensity, gl, layout, densityAdaptor, velocityAdaptor, shearRateAdaptor );
      python_coupling::for_each_noncopyable_type< LatticeModels >  ( boost::ref(exporter) );
      if ( ! exporter.successful() ) {
         PyErr_SetString( PyExc_ValueError, "Adding Pdf Field failed.");
         throw error_already_set();
      }
   }

   struct AdaptorExporter
   {
      template< typename LatticeModel_T>
      void operator() ( python_coupling::NonCopyableWrap<LatticeModel_T> )
      {
         field::exportGhostLayerFieldAdaptor< typename Adaptor< LatticeModel_T >::Density   > ( );
         field::exportGhostLayerFieldAdaptor< typename Adaptor< LatticeModel_T >::Velocity  > ( );
         field::exportGhostLayerFieldAdaptor< typename Adaptor< LatticeModel_T >::ShearRate > ( );
      }
   };

    template< typename LatticeModel >
    struct VelocityAdaptorFromLatticeModel
    {
        typedef typename lbm::Adaptor< LatticeModel>::Velocity type;
    };

    template< typename LatticeModel >
    struct DensityAdaptorFromLatticeModel
    {
        typedef typename lbm::Adaptor< LatticeModel>::Density type;
    };

    template< typename LatticeModel >
    struct ShearRateAdaptorFromLatticeModel
    {
        typedef typename lbm::Adaptor< LatticeModel>::ShearRate type;
    };

    class SweepWrapper
    {
    public:
        SweepWrapper()
        {}

        SweepWrapper( const std::function<void(IBlock*) > & sweepToWrap )
                : sweepToWrap_( sweepToWrap ) {}

        void operator() ( IBlock * block )
        {
           if ( sweepToWrap_ )
              sweepToWrap_( block );
        }
    protected:
        std::function<void(IBlock*) > sweepToWrap_;
    };


    inline shared_ptr<SweepWrapper> makeSplitPureSweep(const shared_ptr<StructuredBlockStorage> & blocks,
                                                       const std::string & pdfFieldIDStr)
    {
       if ( blocks->begin() == blocks->end() ) {
          PyErr_SetString( PyExc_RuntimeError, "No blocks on this process" );
          throw boost::python::error_already_set();
       }

       IBlock & firstBlock = *(  blocks->begin() );
       BlockDataID pdfFieldID = python_coupling::blockDataIDFromString( firstBlock, pdfFieldIDStr );

       typedef lbm::D3Q19< lbm::collision_model::SRT, true >  LM_SRT_Compressible;
       typedef lbm::D3Q19< lbm::collision_model::SRT, false > LM_SRT_Incompressible;
       typedef lbm::D3Q19< lbm::collision_model::TRT, true >  LM_TRT_Compressible;
       typedef lbm::D3Q19< lbm::collision_model::TRT, false > LM_TRT_Incompressible;


       if      (firstBlock.isDataOfType<PdfField<LM_SRT_Compressible>>(pdfFieldID)) {
          return make_shared<SweepWrapper>(SplitPureSweep<LM_SRT_Compressible>(pdfFieldID));
       }
       else if (firstBlock.isDataOfType<PdfField<LM_SRT_Incompressible>>(pdfFieldID)) {
          return make_shared<SweepWrapper>(SplitPureSweep<LM_SRT_Incompressible>(pdfFieldID));
       }
       else if (firstBlock.isDataOfType<PdfField<LM_TRT_Compressible>>(pdfFieldID)) {
          return make_shared<SweepWrapper>(SplitPureSweep<LM_TRT_Compressible>(pdfFieldID));
       }
       else if (firstBlock.isDataOfType<PdfField<LM_TRT_Incompressible>>(pdfFieldID)) {
          return make_shared<SweepWrapper>(SplitPureSweep<LM_TRT_Compressible>(pdfFieldID));
       }
       else {
          PyErr_SetString( PyExc_RuntimeError, "No split-pure sweep available for this lattice model" );
          throw boost::python::error_already_set();
       }
    }
}

template< typename LatticeModels >
struct VelocityAdaptorsFromLatticeModels
{
   typedef typename boost::mpl::transform< LatticeModels, internal::VelocityAdaptorFromLatticeModel<boost::mpl::_1 > >::type type;
};

template< typename LatticeModels >
struct DensityAdaptorsFromLatticeModels
{
   typedef typename boost::mpl::transform< LatticeModels, internal::DensityAdaptorFromLatticeModel<boost::mpl::_1> >::type type;
};

template< typename LatticeModels >
struct ShearRateAdaptorsFromLatticeModels
{
   typedef typename boost::mpl::transform< LatticeModels, internal::ShearRateAdaptorFromLatticeModel<boost::mpl::_1> >::type type;
};

template< typename LatticeModels >
struct AdaptorsFromLatticeModels
{
   typedef typename boost::mpl::copy< typename DensityAdaptorsFromLatticeModels<LatticeModels>::type,
                                      boost::mpl::front_inserter<  typename VelocityAdaptorsFromLatticeModels<LatticeModels>::type  > >
                                     ::type  tmp1;

   typedef typename boost::mpl::copy< typename ShearRateAdaptorsFromLatticeModels<LatticeModels>::type,
                                      boost::mpl::front_inserter<  tmp1  > >
                                     ::type  type;
};



template<typename LatticeModels, typename FlagFields>
void exportBasic()
{
   using namespace boost::python;

   python_coupling::ModuleScope scope("lbm");

   // Lattice Models
   exportForceModels();
   exportCollisionModels();
   python_coupling::for_each_noncopyable_type<LatticeModels>( internal::LatticeModelExporter() );
   python_coupling::for_each_noncopyable_type<LatticeModels>( internal::AdaptorExporter() );

   def( "makeLatticeModel", internal::makeLatticeModel<LatticeModels>,
        ( arg("stencil"), arg("collisionModel"),
          arg("forceModel")= force_model::None(),
          arg("compressible")=false, arg("equilibriumAccuracyOrder")=2)  );

   // PdfField
   python_coupling::for_each_noncopyable_type< LatticeModels > ( internal::PdfFieldExporter() );
   def( "addPdfFieldToStorage", &internal::addPdfFieldToStorage<LatticeModels>,
            (( arg("blocks")                            ),
             ( arg("name")                              ),
             ( arg("latticeModel")                      ),
             ( arg("initialVelocity") = Vector3<real_t>() ),
             ( arg("initialDensity")  = real_t(1)        ),
             ( arg("ghostlayers")     = uint_t(1)        ),
             ( arg("layout")          = field::zyxf      ),
             ( arg("densityAdaptor")  = std::string()    ),
             ( arg("velocityAdaptor") = std::string()    ),
             ( arg("shearRateAdaptor")= std::string()    )
             ) );

   // Split Sweeps
   def("makeSplitPureSweep", &internal::makeSplitPureSweep, (arg("blocks"), arg("pdfFieldId")));
   class_<internal::SweepWrapper, shared_ptr<internal::SweepWrapper>, boost::noncopyable>("SplitPureSweep", no_init)
           .def("__call__", &internal::SweepWrapper::operator());


   auto pythonManager = python_coupling::Manager::instance();
   pythonManager->addBlockDataConversion< typename AdaptorsFromLatticeModels                 < LatticeModels >::type > ();
}



} // namespace lbm
} // namespace walberla


