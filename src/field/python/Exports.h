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
//! \file Exports.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"


#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "python_coupling/helper/ModuleScope.h"
#include "FieldExport.h"
#include "GatherExport.h"
#include "CommunicationExport.h"

namespace walberla {
namespace field {


   //*******************************************************************************************************************
   /*! Exports the field types and corresponding function given in the type sequence to python
   *
   * Automatically exports the Field, GhostLayerField and (if possible) the FlagField version,
   * so it is enough to add the basic Field<> types to the sequence
   * For fields that store and unsigned integer type and have an fSize of one, the FlagField is also exported
   *
   * Example:
    \code
      typedef boost::mpl::vector< Field<real_t,1>,
                                  Field<uint16_t, 1>
                                  Field<Vector3<real_t>, 1 >  FieldVector;
    \endcode

   *  This exports the following types:
   *     - Field<real_t,1>, GhostLayerField<real_t,1>
   *     - Field<uint16_t,1>, GhostLayerField<uint16_t,1>, FlagField<uint16_t>
   *     - Field< Vector3<real_t>, 1 > , GhostLayerField< Vector3<real_t>, 1 >
   *
   *  Additionally the following free functions are exported
   *     - field.createField
   *     - field.createFlagField
   *     - field.addToStorage
   *     - field.createPackInfo
   *     - field.createMPIDatatypeInfo
   *
   *  See Python documentation of waLBerla.field module for details
   *
   * \warning Make sure that the same field type is exported only once!
   */
   //*******************************************************************************************************************
   template<typename FieldTypes>
   void exportModuleToPython()
   {
      python_coupling::ModuleScope fieldModule( "field" );

      exportFields<FieldTypes>();

      exportCommunicationClasses<FieldTypes>();
   }


} // namespace field
} // namespace walberla


#endif // WALBERLA_BUILD_WITH_PYTHON
