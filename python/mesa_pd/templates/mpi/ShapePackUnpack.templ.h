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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/shape/BaseShape.h>
{%- for shape in particle.shapes %}
#include <mesa_pd/data/shape/{{shape}}.h>
{%- endfor %}

#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

#include <memory>

namespace walberla {
namespace mpi {
   template< typename T,    // Element type of SendBuffer
             typename G>    // Growth policy of SendBuffer
   mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf,
                                            const std::shared_ptr<mesa_pd::data::BaseShape>& bs )
   {
      buf.addDebugMarker( "up" );
      buf << bs->getShapeType();
      bs->pack(buf);
      return buf;
   }

   template< typename T>    // Element type  of RecvBuffer
   mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf,
                                          std::shared_ptr<mesa_pd::data::BaseShape>& bs )
   {
      using namespace mesa_pd::data;

      buf.readDebugMarker( "up" );

      mesa_pd::data::BaseShape::ShapeTypeT shapeType;
      buf >> shapeType;
      switch (shapeType)
      {
         case BaseShape::INVALID_SHAPE :
            bs = std::make_unique<mesa_pd::data::BaseShape>();
            bs->unpack(buf);
            break;
         {%- for shape in particle.shapes %}
         case {{shape}}::SHAPE_TYPE :
            bs = std::make_unique<mesa_pd::data::{{shape}}>();
            bs->unpack(buf);
            break;
         {%- endfor %}
         default : WALBERLA_ABORT("Shape type (" << shapeType << ") could not be determined!");
      }
      return buf;
   }
} //namespace mpi
} //namespace walberla
