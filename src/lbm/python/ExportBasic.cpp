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
//! \file PythonExports.cpp
//! \ingroup domain_decomposition
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "python_coupling/PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "ExportBasic.h"
#include "core/logging/Logging.h"
#include "python_coupling/Manager.h"
#include "python_coupling/helper/ModuleScope.h"

#include "domain_decomposition/BlockDataID.h"

#include "lbm/field/Adaptors.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/ForceModel.h"


using namespace boost::python;


namespace walberla {
namespace lbm {


namespace internal {

   using namespace boost::python;


   template<typename ForceModel_T> bool FM_getShiftMacVel( const ForceModel_T & ) { return ForceModel_T::shiftMacVel; }
   template<typename ForceModel_T> bool FM_getShiftEquVel( const ForceModel_T & ) { return ForceModel_T::shiftEquVel; }
   template<typename ForceModel_T> bool FM_getConstant   ( const ForceModel_T & ) { return ForceModel_T::constant;    }

   template<typename ForceModel_T>
   void addForceModelDefs( class_< ForceModel_T> & c )
   {
      c.add_property( "shiftMacVel",            FM_getShiftMacVel<ForceModel_T> )
       .add_property( "shiftEquVel",            FM_getShiftEquVel<ForceModel_T> )
       .add_property( "constant",               FM_getConstant   <ForceModel_T>)
       .def( "setConstantBodyForceIfPossible", &ForceModel_T::setConstantBodyForceIfPossible )
      ;
   }


   list getMRTRelaxationRates( const collision_model::D3Q19MRT & l )
   {
      list result;
      for(uint_t i=0; i<19; ++i )
         result.append( l.s(i) );
      return result;
   }

   list getCumulantRelaxationRates( const collision_model::D3Q27Cumulant & l )
   {
      list result;
      for(uint_t i=0; i<10; ++i )
         result.append( l.omega(i) );
      return result;
   }
}

void exportForceModels()
{
   static bool alreadyExported = false;
   if( alreadyExported )
      return;
   alreadyExported = true;

   using namespace internal;
   using namespace force_model;

   python_coupling::ModuleScope scope("forceModels");

   typedef GhostLayerField<Vector3<real_t>,1 > VecField;

   auto export1 = class_< None >("NoForce" ); // "None" is a Python keyword -> renamed to NoForce
   addForceModelDefs( export1 );

   auto export2 = class_< SimpleConstant >("SimpleConstant",
                  init<Vector3<real_t> , uint_t> ( (arg("force"), arg("level")=uint_t(0) ) ) );
   const Vector3<real_t>& (SimpleConstant::*p_SimpleConstantForce)() const = &SimpleConstant::force;
   export2.def("force", p_SimpleConstantForce, return_value_policy<copy_const_reference>());
   addForceModelDefs( export2 );

   auto export3 = class_< EDMField<VecField> > ("EDMField",
                  init< BlockDataID> ( (arg("forceFieldID") )) );
   addForceModelDefs( export3 );

   auto export4 = class_< LuoConstant > ("LuoConstant",
                   init<Vector3<real_t> , uint_t> ( (arg("force"), arg("level")=uint_t(0) ) ) );
   const Vector3<real_t>& (LuoConstant::*p_LuoConstantForce)() const = &LuoConstant::force;
   export4.def("force", p_LuoConstantForce, return_value_policy<copy_const_reference>() );
   addForceModelDefs( export4 );

   auto export5 = class_< LuoField<VecField> > ("LuoField",
                  init< BlockDataID> ( (arg("forceFieldID") )) );
   addForceModelDefs( export5 );

   auto export6 = class_< GuoConstant > ("GuoConstant",
                   init<Vector3<real_t> , uint_t> ( (arg("force"), arg("level")=uint_t(0) ) ) );
   const Vector3<real_t>& (GuoConstant::*p_GuoConstantForce)() const = &GuoConstant::force;
   export6.def("force", p_GuoConstantForce, return_value_policy<copy_const_reference>() );
   addForceModelDefs( export6 );

   auto export7 = class_< GuoField<VecField> > ("GuoField",
                  init< BlockDataID> ( (arg("forceFieldID") )) );
   addForceModelDefs( export7 );

   auto export8 = class_< Correction<VecField> > ("Correction",
                  init< BlockDataID> ( (arg("previousMomentumDensityFieldID") )) );
   addForceModelDefs( export8 );
}



shared_ptr< collision_model::SRTField<GhostLayerField<real_t,1> > > createSRTFieldLatticeModel( const shared_ptr<StructuredBlockStorage> & bs, const std::string & blockDataName, uint_t level )
{
   auto blockDataID = python_coupling::blockDataIDFromString( *bs, blockDataName );
   return make_shared< collision_model::SRTField<GhostLayerField<real_t,1> > >( blockDataID, level );
}

void exportCollisionModels()
{
   static bool alreadyExported = false;
   if( alreadyExported )
      return;
   alreadyExported = true;

   using namespace internal;

   python_coupling::ModuleScope scope("collisionModels");

   using collision_model::SRT;
   using collision_model::SRTField;
   using collision_model::TRT;
   using collision_model::D3Q19MRT;
   using collision_model::D3Q27Cumulant;

   def( "levelDependentRelaxationParameter", collision_model::levelDependentRelaxationParameter, (arg("targetLevel"), arg("parameterLevel"), arg("level")) );
   def( "viscosityFromOmega", collision_model::viscosityFromOmega, (arg("omega") ) );
   def( "omegaFromViscosity", collision_model::omegaFromViscosity, (arg("viscosity") ) );

   // SRT
   {
      real_t ( SRT::*ptr_omega     )() const = &SRT::omega;
      real_t ( SRT::*ptr_viscosity )() const = &SRT::viscosity;

      class_< SRT > ("SRT", init<real_t, uint_t>( ( arg("omega"), arg("level")=uint_t(0) )  ) )
          .add_property("omega", ptr_omega)
          .add_property("viscosity", ptr_viscosity )
          .add_property("level", &SRT::level )
          .def("reset", &SRT::reset, (arg("omega"), arg("level")=uint_t(0) ) )
       ;
   }

   // SRTField
   {
      typedef SRTField<GhostLayerField<real_t,1> > SRTField_T;
      class_< SRTField_T >( "SRTField", no_init )
          .def( "__init__", make_constructor( &createSRTFieldLatticeModel) )
       ;
   }

   // TRT
   {
      real_t ( TRT::*ptr_lambda_e )() const = &TRT::lambda_e;
      real_t ( TRT::*ptr_lambda_d )() const = &TRT::lambda_d;
      real_t ( TRT::*ptr_viscosity )() const = &TRT::viscosity;

      class_< TRT> ( "TRT", init<real_t, real_t, uint_t>( (arg("lambda_e"), arg("lambda_d"), arg("level")=uint_t(0) ) ) )
         .add_property("lambda_e", ptr_lambda_e )
         .add_property("lambda_d", ptr_lambda_d )
         .add_property("lambda_o", ptr_lambda_d )
         .add_property("viscosity", ptr_viscosity )
         .add_property("level", &TRT::level )
         .def( "reset", &TRT::reset, ( arg("lambda_e"), arg("lambda_d"), arg("level")=uint_t(0) ) )
         .def( "resetWithMagicNumber",     &TRT::resetWithMagicNumber,     (arg("omega"), arg("magicNumber") = TRT::threeSixteenth,arg("level")=uint_t(0) ) )
         .def( "constructWithMagicNumber", &TRT::constructWithMagicNumber, (arg("omega"), arg("magicNumber") = TRT::threeSixteenth, arg("level")=uint_t(0) ) )
         .staticmethod("constructWithMagicNumber")
      ;
   }

   // MRT
   {
      real_t ( D3Q19MRT::*ptr_viscosity )() const = &D3Q19MRT::viscosity;

      class_< D3Q19MRT >( "D3Q19MRT", init<real_t, real_t, real_t, real_t, real_t, real_t, uint_t>(
                                         ( arg("s1"), arg("s2"), arg("s4"), arg("s9"), arg("s10"), arg("s16"), arg("level")=uint_t(0) ) ))
          .add_property( "relaxationRates", getMRTRelaxationRates )
          .add_property( "viscosity", ptr_viscosity )
          .def( "constructTRT", &D3Q19MRT::constructTRT, (arg("lambda_e"), arg("lambda_d"), arg("level")=uint_t(0) ) )
          .staticmethod("constructTRT")
          .def( "constructTRTWithMagicNumber", &D3Q19MRT::constructTRTWithMagicNumber, (arg("omega"), arg("magicNumber")=D3Q19MRT::threeSixteenth, arg("level")=uint_t(0) ) )
          .staticmethod("constructTRTWithMagicNumber")
          .def( "constructPan", &D3Q19MRT::constructPan, (arg("lambda_e"), arg("lambda_d"), arg("level")=uint_t(0) ) )
          .staticmethod("constructPan")
          .def( "constructPanWithMagicNumber", &D3Q19MRT::constructPanWithMagicNumber, (arg("omega"), arg("magicNumber")=D3Q19MRT::threeSixteenth, arg("level")=uint_t(0) ) )
          .staticmethod("constructPanWithMagicNumber")
       ;
   }

   // Cumulant
   {
      real_t ( D3Q27Cumulant::*ptr_viscosity )() const = &D3Q27Cumulant::viscosity;

      class_< D3Q27Cumulant >( "D3Q27Cumulant", init<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t,real_t, uint_t>(
                                         ( arg("omega1"), arg("omega2")=real_t(1), arg("omega3")=real_t(1), arg("omega4")=real_t(1), arg("omega5")=real_t(1), arg("omega6")=real_t(1), arg("omega7")=real_t(1), arg("omega8")=real_t(1), arg("omega9")=real_t(1), arg("omega10")=real_t(1), arg("level")=uint_t(0) ) ))
          .add_property( "relaxationRates", getCumulantRelaxationRates )
          .add_property( "viscosity", ptr_viscosity )
       ;

   }

}



} // namespace lbm
} // namespace walberla


#endif //WALBERLA_BUILD_WITH_PYTHON
