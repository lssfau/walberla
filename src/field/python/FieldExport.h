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
//! \file FieldExport.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#ifdef WALBERLA_BUILD_WITH_PYTHON


#include <string>

namespace walberla {
namespace field {


   //*******************************************************************************************************************
   /*! Exports all Fields given in the Sequence
   *
   * Put only Fields in the sequence! The corresponding GhostLayerFields and FlagFields are exported automatically
   *
   * \warning Make sure that the same adaptor type is exported only once!
   */
   //*******************************************************************************************************************
   template<typename FieldTypes >
   void exportFields();



   //*******************************************************************************************************************
   /*! Exports all GhostLayerFieldAdaptors given in the Sequence
   *
   * \warning Make sure that the same adaptor type is exported only once!
   */
   //*******************************************************************************************************************
   template<typename AdaptorTypes>
   void exportGhostLayerFieldAdaptors();

   template<typename AdaptorType>
   void exportGhostLayerFieldAdaptor();


} // namespace field
} // namespace walberla

#include "FieldExport.impl.h"


#endif //WALBERLA_BUILD_WITH_PYTHON
