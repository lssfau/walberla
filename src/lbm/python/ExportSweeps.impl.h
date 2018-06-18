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
//! \file ExportSweeps.impl.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "python_coupling/helper/MplHelpers.h"
#include "python_coupling/helper/BoostPythonHelpers.h"
#include "python_coupling/helper/ModuleScope.h"
#include "python_coupling/helper/BlockStorageExportHelpers.h"

#include "python_coupling/Manager.h"

#include "core/logging/Logging.h"

#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "field/python/FieldExport.h"
#include "field/adaptors/AdaptorCreators.h"

namespace walberla {
namespace lbm {


namespace internal
{
   template< typename LM, typename Filter, typename In, typename Out >
   void addSweepMemberFunctions( boost::python::class_< CellwiseSweep<LM,Filter,In,Out>, shared_ptr<CellwiseSweep<LM,Filter,In,Out> > > & c )
   {
      using namespace boost::python;

      typedef CellwiseSweep<LM,Filter,In,Out> Sweep;
      c.def( "stream",       &Sweep::stream,        (arg("block"),arg("numberOfGhostLayersToInclude")=uint_t(0) ) )
       .def( "collide",      &Sweep::collide,       (arg("block"),arg("numberOfGhostLayersToInclude")=uint_t(0) ) )
       .def("streamCollide", &Sweep::streamCollide, (arg("block"),arg("numberOfGhostLayersToInclude")=uint_t(0) ) )
       .def("__call__",      &Sweep::operator(),    (arg("block"),arg("numberOfGhostLayersToInclude")=uint_t(0) ) )
      ;
   }

   struct CellwiseSweepExporterWithoutFlagField
   {
      template<typename LatticeModel_T>
      void operator()( python_coupling::NonCopyableWrap<LatticeModel_T>  )
      {
         using namespace boost::python;

         typedef GhostLayerField< real_t, 3 > VelocityField_T;

         typedef CellwiseSweep<LatticeModel_T,
                  field::DefaultEvaluationFilter,
                  DefaultDensityEquilibriumVelocityCalculation,
                  DefaultDensityVelocityCallback
                  > Sweep1;
         typedef CellwiseSweep<LatticeModel_T,
                  field::DefaultEvaluationFilter,
                  DefaultDensityEquilibriumVelocityCalculation,
                  VelocityCallback<VelocityField_T>
                  > Sweep2;

         auto sweep1 = class_< Sweep1, shared_ptr<Sweep1> > ( "CellwiseSweep", no_init );
         auto sweep2 = class_< Sweep2, shared_ptr<Sweep2> > ( "CellwiseSweep", no_init );
         addSweepMemberFunctions( sweep1 );
         addSweepMemberFunctions( sweep2 );
      }
   };

   struct CellwiseSweepExporterWithFlagField
   {
      template<typename LatticeModel_FlagField_Pair>
      void operator()( python_coupling::NonCopyableWrap<LatticeModel_FlagField_Pair>  )
      {
         using namespace boost::python;

         typedef typename LatticeModel_FlagField_Pair::first  LatticeModel_T;
         typedef typename LatticeModel_FlagField_Pair::second FlagField_T;

         typedef GhostLayerField< real_t,3> VelocityField_T;

         typedef CellwiseSweep<LatticeModel_T,
                  field::FlagFieldEvaluationFilter<FlagField_T>,
                  DefaultDensityEquilibriumVelocityCalculation,
                  DefaultDensityVelocityCallback
                  > Sweep1;
         typedef  CellwiseSweep<LatticeModel_T,
                  field::FlagFieldEvaluationFilter<FlagField_T>,
                  DefaultDensityEquilibriumVelocityCalculation,
                  VelocityCallback<VelocityField_T>
                  > Sweep2;

         auto sweep1 = class_< Sweep1, shared_ptr<Sweep1> > ( "CellwiseSweep", no_init );
         auto sweep2 = class_< Sweep2, shared_ptr<Sweep2> > ( "CellwiseSweep", no_init );
         addSweepMemberFunctions( sweep1 );
         addSweepMemberFunctions( sweep2 );
      }
   };


   class CellwiseSweepCreator
   {
   public:

      CellwiseSweepCreator( const shared_ptr<StructuredBlockStorage> & blocks, BlockDataID pdfFieldID,
                            const std::string & flagFieldStringID, const std::string & velocityFieldStringID,
                            const Set< FlagUID > & cellsToEvaluate )
         : blocks_( blocks ), pdfFieldID_( pdfFieldID ), flagFieldStringID_( flagFieldStringID ),
           velocityFieldStringID_( velocityFieldStringID ), cellsToEvaluate_( cellsToEvaluate )
      {}


      template<typename LatticeModel_FlagField_Pair>
      void operator()( python_coupling::NonCopyableWrap<LatticeModel_FlagField_Pair>  )
      {
         typedef GhostLayerField<real_t,3> VelocityField_T;

         using namespace boost::python;
         using python_coupling::blockDataIDFromString;


         typedef typename LatticeModel_FlagField_Pair::first  LatticeModel_T;
         typedef typename LatticeModel_FlagField_Pair::second FlagField_T;

         if ( blocks_->begin() == blocks_->end() )
            return;

         IBlock & firstBlock = *(  blocks_->begin() );

         if( ! firstBlock.isDataClassOrSubclassOf<PdfField<LatticeModel_T> >(pdfFieldID_) )
            return;

         if( flagFieldStringID_.size() > 0 )
         {
            // Use FlagFilter
            BlockDataID flagFieldID = blockDataIDFromString( firstBlock, flagFieldStringID_ );
            if ( ! firstBlock.isDataClassOrSubclassOf< FlagField_T>( flagFieldID ) )
               return;

            if ( velocityFieldStringID_.size() > 0 )
            {
               BlockDataID velocityFieldID = blockDataIDFromString( firstBlock, velocityFieldStringID_ );
               result_ = object( makeCellwiseSweep<LatticeModel_T, FlagField_T, VelocityField_T >( pdfFieldID_, flagFieldID,
                                                                                                   cellsToEvaluate_, velocityFieldID ) );
            }
            else
            {
               result_ = object( makeCellwiseSweep<LatticeModel_T, FlagField_T>( pdfFieldID_, flagFieldID, cellsToEvaluate_ ) );
            }
         }
         else
         {
            if ( velocityFieldStringID_.size() > 0 )
            {
               BlockDataID velocityFieldID = blockDataIDFromString( firstBlock, velocityFieldStringID_ );
               result_ = object ( makeCellwiseSweep<LatticeModel_T,
                                           field::DefaultEvaluationFilter,
                                           DefaultDensityEquilibriumVelocityCalculation,
                                           VelocityCallback<VelocityField_T>
                                           >
                                           ( pdfFieldID_, field::DefaultEvaluationFilter(),
                                             DefaultDensityEquilibriumVelocityCalculation(),
                                             VelocityCallback<VelocityField_T>( velocityFieldID ) ) );
            }
            else
            {
               result_ = object( makeCellwiseSweep<LatticeModel_T>( pdfFieldID_ ) );
            }
         }
      }

      boost::python::object getResult() const { return result_; }

   private:
      shared_ptr<StructuredBlockStorage> blocks_;
      BlockDataID pdfFieldID_;
      std::string flagFieldStringID_;
      std::string velocityFieldStringID_;
      Set< FlagUID >  cellsToEvaluate_;

      boost::python::object result_;
   };


   template<typename LatticeModel_FlagField_Pairs>
   boost::python::object makeCellwiseSweep( const shared_ptr<StructuredBlockStorage> & bs,
                                            const std::string & pdfFieldStringID,
                                            const std::string & flagFieldStringID, boost::python::list flagList,
                                            const std::string & velocityFieldStringID )
   {
      using namespace boost::python;

      if ( bs->begin() == bs->end() )
         return object();
      IBlock & firstBlock = *(  bs->begin() );

      if( flagFieldStringID.size() > 0 || len(flagList) > 0 )
      {
         if ( !( flagFieldStringID.size() > 0 && len(flagList) > 0 ) )
         {
            PyErr_SetString( PyExc_ValueError, "Pass flagFieldID and flagList or none of them.");
            throw error_already_set();
         }
      }

      BlockDataID pdfFieldID = python_coupling::blockDataIDFromString( firstBlock, pdfFieldStringID );

      auto flagUidSet = python_coupling::uidSetFromStringContainer< FlagUID >( flagList );

      CellwiseSweepCreator creator( bs, pdfFieldID, flagFieldStringID, velocityFieldStringID, flagUidSet );
      python_coupling::for_each_noncopyable_type<LatticeModel_FlagField_Pairs>( std::ref( creator ) );

      if ( creator.getResult() == object() )
      {
         PyErr_SetString( PyExc_ValueError, "No Sweep available with these options.");
         throw error_already_set();
      }
      return creator.getResult();
   }

}

template<typename LatticeModels, typename FlagFields>
void exportSweeps()
{
   using namespace boost::python;


   python_coupling::ModuleScope scope("lbm");

   // Cellwise Sweep
   typedef python_coupling::combine_vectors<LatticeModels,FlagFields> LatticeModel_FlagField_Pairs;
   python_coupling::for_each_noncopyable_type<LatticeModel_FlagField_Pairs>( internal::CellwiseSweepExporterWithFlagField() );
   python_coupling::for_each_noncopyable_type<LatticeModels>               ( internal::CellwiseSweepExporterWithoutFlagField() );
   def( "makeCellwiseSweep", internal::makeCellwiseSweep<LatticeModel_FlagField_Pairs>,
                  ( arg("blocks"), arg("pdfFieldID"),
                    arg("flagFieldID")=std::string(), arg("flagList")=boost::python::list(),
                    arg("velocityFieldID")=std::string() ) );

}



} // namespace lbm
} // namespace walberla


