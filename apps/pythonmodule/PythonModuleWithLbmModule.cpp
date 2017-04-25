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
//! \file PythonModule.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================
#include "python_coupling/PythonWrapper.h"
#include "waLBerlaDefinitions.h"
#include "blockforest/python/Exports.h"
#include "field/GhostLayerField.h"
#include "field/python/Exports.h"
#include "geometry/python/Exports.h"
#include "lbm/all.h"
#include "lbm/lattice_model/ForceModel.h"
#include "lbm/python/Exports.h"
#include "postprocessing/python/Exports.h"
#include "python_coupling/Manager.h"
#include "python_coupling/helper/ModuleInit.h"
#include "stencil/all.h"
#include "timeloop/python/Exports.h"
#include "vtk/python/Exports.h"

#include <boost/mpl/vector.hpp>
#include <boost/mpl/insert_range.hpp>


namespace bmpl = boost::mpl;
using namespace walberla;

typedef bmpl::vector<
            Field<walberla::real_t,1>,
            Field<walberla::real_t,2>,
            Field<walberla::real_t,3>,
            Field<walberla::real_t,4>,
            Field<walberla::real_t,5>,
            Field<walberla::real_t,9>,
            Field<walberla::real_t,19>,
            Field<walberla::real_t,27>,

            Field<walberla::int8_t,1>,
            Field<walberla::int16_t,1>,
            Field<walberla::int32_t,1>,

            Field<walberla::int64_t,1>,
            Field<walberla::int64_t,2>,
            Field<walberla::int64_t,3>,
            Field<walberla::int64_t,4>,
            Field<walberla::int64_t,5>,

            Field<walberla::uint8_t,1>,
            Field<walberla::uint16_t,1>,
            Field<walberla::uint32_t,1>
      > FieldTypes;


typedef bmpl::vector<
                      GhostLayerField<walberla::real_t,1>,
                      GhostLayerField<walberla::real_t,3>
                    >  FieldTypesForMeshGeneration;


typedef bmpl::vector< FlagField<walberla::uint8_t>,
                      FlagField<walberla::uint16_t> > FlagFieldTypes;

typedef bmpl::vector< stencil::D2Q5,
                      stencil::D2Q9,
                      stencil::D3Q7,
                      stencil::D3Q19,
                      stencil::D3Q27 > Stencils;

typedef lbm::collision_model::SRTField< GhostLayerField<walberla::real_t,1> > SRTField_T;

typedef GhostLayerField<walberla::real_t,3> VecField_T;


// -------------------------    D2Q9 -----------------------------------------------------------------------------

typedef bmpl::list< lbm::D2Q9  < lbm::collision_model::SRT, false, lbm::force_model::None, 2>,
                    lbm::D2Q9  < lbm::collision_model::SRT, false, lbm::force_model::SimpleConstant, 2>,
                    lbm::D2Q9  < lbm::collision_model::SRT, false, lbm::force_model::GuoConstant, 2>,
                    lbm::D2Q9  < lbm::collision_model::SRT, false, lbm::force_model::LuoConstant, 2>,
                    lbm::D2Q9  < lbm::collision_model::SRT, true, lbm::force_model::None, 2>,
                    lbm::D2Q9  < lbm::collision_model::SRT, true, lbm::force_model::SimpleConstant, 2>,
                    lbm::D2Q9  < lbm::collision_model::SRT, true, lbm::force_model::GuoConstant, 2>,
                    lbm::D2Q9  < lbm::collision_model::SRT, false, lbm::force_model::GuoField<VecField_T>, 2>,
                    lbm::D2Q9  < lbm::collision_model::SRT, true, lbm::force_model::LuoConstant, 2>,
                    lbm::D2Q9  < lbm::collision_model::SRT, true, lbm::force_model::None, 1>,
                    lbm::D2Q9  < lbm::collision_model::SRT, false, lbm::force_model::None, 1>
> LM_D2Q19_SRT;

typedef bmpl::list<   lbm::D2Q9  < lbm::collision_model::TRT, false, lbm::force_model::None, 2>,
                      lbm::D2Q9  < lbm::collision_model::TRT, false, lbm::force_model::SimpleConstant, 2>,
                      lbm::D2Q9  < lbm::collision_model::TRT, false, lbm::force_model::GuoConstant, 2>,
                      lbm::D2Q9  < lbm::collision_model::TRT, false, lbm::force_model::LuoConstant, 2>,

                      lbm::D2Q9  < lbm::collision_model::TRT, true, lbm::force_model::None, 2>,
                      lbm::D2Q9  < lbm::collision_model::TRT, true, lbm::force_model::SimpleConstant, 2>,
                      lbm::D2Q9  < lbm::collision_model::TRT, true, lbm::force_model::GuoConstant, 2>,
                      lbm::D2Q9  < lbm::collision_model::TRT, true, lbm::force_model::LuoConstant, 2>,
                      lbm::D2Q9  < lbm::collision_model::TRT, true, lbm::force_model::None, 1>,
                      lbm::D2Q9  < lbm::collision_model::TRT, false, lbm::force_model::None, 1>
> LM_D2Q19_TRT;

typedef lbm::collision_model::SRTField<GhostLayerField<real_t,1> > SRTFieldCollisionModel;
typedef bmpl::list<   lbm::D2Q9  < SRTFieldCollisionModel, false, lbm::force_model::None, 2>,
                      lbm::D2Q9  < SRTFieldCollisionModel, false, lbm::force_model::SimpleConstant, 2>,
                      lbm::D2Q9  < SRTFieldCollisionModel, false, lbm::force_model::GuoConstant, 2>,
                      lbm::D2Q9  < SRTFieldCollisionModel, false, lbm::force_model::LuoConstant, 2>,

                      lbm::D2Q9  < SRTFieldCollisionModel, true, lbm::force_model::None, 2>,
                      lbm::D2Q9  < SRTFieldCollisionModel, true, lbm::force_model::SimpleConstant, 2>,
                      lbm::D2Q9  < SRTFieldCollisionModel, true, lbm::force_model::GuoConstant, 2>,
                      lbm::D2Q9  < SRTFieldCollisionModel, true, lbm::force_model::LuoConstant, 2>
> LM_D2Q9_SRTField;

// -------------------------    D3Q19 ----------------------------------------------------------------------------


typedef bmpl::list<   lbm::D3Q19  < lbm::collision_model::SRT, false, lbm::force_model::None, 2>,
                      lbm::D3Q19  < lbm::collision_model::SRT, false, lbm::force_model::SimpleConstant, 2>,
                      lbm::D3Q19  < lbm::collision_model::SRT, false, lbm::force_model::GuoConstant, 2>,
                      lbm::D3Q19  < lbm::collision_model::SRT, false, lbm::force_model::LuoConstant, 2>,

                      lbm::D3Q19  < lbm::collision_model::SRT, true, lbm::force_model::None, 2>,
                      lbm::D3Q19  < lbm::collision_model::SRT, true, lbm::force_model::SimpleConstant, 2>,
                      lbm::D3Q19  < lbm::collision_model::SRT, true, lbm::force_model::GuoConstant, 2>,
                      lbm::D3Q19  < lbm::collision_model::SRT, true, lbm::force_model::LuoConstant, 2>
> LM_D3Q19_SRT;

typedef bmpl::list<   lbm::D3Q19  < lbm::collision_model::TRT, false, lbm::force_model::None, 2>,
                      lbm::D3Q19  < lbm::collision_model::TRT, false, lbm::force_model::SimpleConstant, 2>,
                      lbm::D3Q19  < lbm::collision_model::TRT, false, lbm::force_model::GuoConstant, 2>,
                      lbm::D3Q19  < lbm::collision_model::TRT, false, lbm::force_model::LuoConstant, 2>,

                      lbm::D3Q19  < lbm::collision_model::TRT, true, lbm::force_model::None, 2>,
                      lbm::D3Q19  < lbm::collision_model::TRT, true, lbm::force_model::SimpleConstant, 2>,
                      lbm::D3Q19  < lbm::collision_model::TRT, true, lbm::force_model::GuoConstant, 2>,
                      lbm::D3Q19  < lbm::collision_model::TRT, true, lbm::force_model::LuoConstant, 2>
> LM_D3Q19_TRT;


typedef bmpl::list<lbm::D3Q19 < lbm::collision_model::D3Q19MRT, false, lbm::force_model::None, 2>,
                   lbm::D3Q19 < lbm::collision_model::D3Q19MRT, false, lbm::force_model::SimpleConstant, 2>,
                   lbm::D3Q19 < lbm::collision_model::D3Q19MRT, false, lbm::force_model::GuoConstant, 2>,
                   lbm::D3Q19 < lbm::collision_model::D3Q19MRT, false, lbm::force_model::LuoConstant, 2>,

                   lbm::D3Q27 < lbm::collision_model::D3Q27Cumulant, true, lbm::force_model::None, 2>

> LM_D3Q19_Extra;



typedef lbm::collision_model::SRTField<GhostLayerField<real_t,1> > SRTFieldCollisionModel;
typedef bmpl::list<   lbm::D3Q19  < SRTFieldCollisionModel, false, lbm::force_model::None, 2>,
                      lbm::D3Q19  < SRTFieldCollisionModel, false, lbm::force_model::SimpleConstant, 2>,
                      lbm::D3Q19  < SRTFieldCollisionModel, false, lbm::force_model::GuoConstant, 2>,
                      lbm::D3Q19  < SRTFieldCollisionModel, false, lbm::force_model::LuoConstant, 2>,

                      lbm::D3Q19  < SRTFieldCollisionModel, true, lbm::force_model::None, 2>,
                      lbm::D3Q19  < SRTFieldCollisionModel, true, lbm::force_model::SimpleConstant, 2>,
                      lbm::D3Q19  < SRTFieldCollisionModel, true, lbm::force_model::GuoConstant, 2>,
                      lbm::D3Q19  < SRTFieldCollisionModel, true, lbm::force_model::LuoConstant, 2>
> LM_D3Q19_SRTField;


typedef bmpl::insert_range< LM_D2Q19_SRT, bmpl::end< LM_D2Q19_SRT >::type, LM_D2Q19_TRT>     ::type D2Q9_SRT_TRT;
typedef bmpl::insert_range< D2Q9_SRT_TRT, bmpl::end< D2Q9_SRT_TRT >::type, LM_D2Q9_SRTField>::type LatticeModels2D;


typedef bmpl::insert_range< LM_D3Q19_SRT,  bmpl::end< LM_D3Q19_SRT > ::type, LM_D3Q19_TRT>::type D3Q19_SRT_TRT;
typedef bmpl::insert_range< D3Q19_SRT_TRT, bmpl::end< D3Q19_SRT_TRT >::type, LM_D3Q19_Extra>::type D3Q19_SRT_TRT_MRT;
typedef bmpl::insert_range< D3Q19_SRT_TRT_MRT, bmpl::end< D3Q19_SRT_TRT_MRT >::type, LM_D3Q19_SRTField>::type LatticeModels3D;


typedef bmpl::insert_range< LatticeModels2D, bmpl::end< LatticeModels2D >::type, LatticeModels3D>::type LatticeModels;

typedef typename lbm::AdaptorsFromLatticeModels<LatticeModels>::type LBMAdaptors;



using namespace walberla;

struct InitObject
{
   InitObject()
   {
      namespace bmpl = boost::mpl;

      auto pythonManager = python_coupling::Manager::instance();

      // Field
      pythonManager->addExporterFunction( field::exportModuleToPython<FieldTypes> );
      pythonManager->addExporterFunction( field::exportGatherFunctions<bmpl::joint_view<FieldTypes,LBMAdaptors>::type > );
      pythonManager->addBlockDataConversion< bmpl::joint_view<FieldTypes,LBMAdaptors>::type >() ;

      // Blockforest
      pythonManager->addExporterFunction( blockforest::exportModuleToPython<Stencils> );

      // Geometry
      pythonManager->addExporterFunction( geometry::exportModuleToPython );

      // VTK
      pythonManager->addExporterFunction( vtk::exportModuleToPython );

      // LBM
      typedef bmpl::vector< FlagField<walberla::uint8_t> > FlagFieldTypesForLBM;
      pythonManager->addExporterFunction( lbm::exportBasic<LatticeModels, FlagFieldTypesForLBM> );
      pythonManager->addExporterFunction( lbm::exportBoundary<LatticeModels, FlagFieldTypesForLBM> );
      pythonManager->addExporterFunction( lbm::exportSweeps<LatticeModels, FlagFieldTypesForLBM> );

      // Postprocessing
      pythonManager->addExporterFunction( postprocessing::exportModuleToPython<FieldTypesForMeshGeneration, FlagFieldTypes> );

      // Timeloop
      pythonManager->addExporterFunction( timeloop::exportModuleToPython );

      python_coupling::initWalberlaForPythonModule();
   }
};
InitObject globalInitObject;

