//======================================================================================================================
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
//! \\file {{class_name}}.h
//! \\author lbmpy
//======================================================================================================================

#pragma once

#include "stencil/Directions.h"
#include "core/cell/CellInterval.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "domain_decomposition/IBlock.h"
#include "field/GhostLayerField.h"
{% if target is equalto 'gpu' -%}
#include "gpu/GPUWrapper.h"
#include "gpu/GPUField.h"
{%- endif %}

{% if target is equalto 'cpu' -%}
#define FUNC_PREFIX
{%- elif target is equalto 'gpu' -%}
#define FUNC_PREFIX __global__
{%- endif %}

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif

namespace walberla{
namespace {{namespace}} {

class {{class_name}} {

public:
   using PdfField_T = {{src_field | field_type(is_gpu=is_gpu)}};
   using value_type = typename PdfField_T::value_type;

   static const bool inplace = {% if inplace -%} true {%- else -%} false {%- endif -%};

   /**
    * Packs all pdfs from the given cell interval to the send buffer.
    */
   void packAll(
      {{- [ "PdfField_T * " + src_field.name, "CellInterval & ci",
            "unsigned char * outBuffer", kernels['packAll'].kernel_selection_parameters,
             ["gpuStream_t stream = nullptr"] if is_gpu else []]
          | type_identifier_list -}}
   ) const;

   /**
    * Unpacks all pdfs from the send buffer to the given cell interval.
    */
   void unpackAll(
      {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
            "unsigned char * inBuffer", kernels['unpackAll'].kernel_selection_parameters,
             ["gpuStream_t stream = nullptr"] if is_gpu else []]
          | type_identifier_list -}}
   ) const;

   /**
    * Copies data between two blocks on the same process.
    * All pdfs from the sending interval are copied onto the receiving interval.
    */
   void localCopyAll(
      {{- [ "PdfField_T * " + src_field.name, "CellInterval & srcInterval",
            "PdfField_T * " + dst_field.name, "CellInterval & dstInterval",
          kernels['localCopyAll'].kernel_selection_parameters,
             ["gpuStream_t stream = nullptr"] if is_gpu else []]
      | type_identifier_list -}}
   ) const;

   /**
    * Packs only those populations streaming in directions aligned with the sending direction dir from the given
    * cell interval.
    * For example, in 2D, if dir == N, the pdfs streaming in directions NW, N, NE are packed.
    */
   void packDirection(
      {{- [ "PdfField_T * " + src_field.name, "CellInterval & ci",
            "unsigned char * outBuffer", kernels['packDirection'].kernel_selection_parameters,
             ["gpuStream_t stream = nullptr"] if is_gpu else []]
          | type_identifier_list -}}
   ) const;

   /**
    * Unpacks only those populations streaming in directions aligned with the sending direction dir to the given
    * cell interval.
    * For example, in 2D, if dir == N, the pdfs streaming in directions NW, N, NE are unpacked.
    */
   void unpackDirection(
      {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
            "unsigned char * inBuffer", kernels['unpackDirection'].kernel_selection_parameters,
             ["gpuStream_t stream = nullptr"] if is_gpu else []]
          | type_identifier_list -}}
   ) const;

   /**
    * Copies data between two blocks on the same process.
    * PDFs streaming aligned with the direction dir are copied from the sending interval
    * onto the receiving interval.
    */
   void localCopyDirection(
      {{- [ "PdfField_T * " + src_field.name, "CellInterval & srcInterval",
            "PdfField_T * " + dst_field.name, "CellInterval & dstInterval",
            kernels['localCopyDirection'].kernel_selection_parameters,
             ["gpuStream_t stream = nullptr"] if is_gpu else []]
          | type_identifier_list -}}
   ) const;

   /**
    * Returns the number of bytes that will be packed from / unpacked to the cell interval
    * when using packDirection / unpackDirection
    * @param ci  The cell interval
    * @param dir The communication direction
    * @return    The required size of the buffer, in bytes
    */
   uint_t size (CellInterval & ci, stencil::Direction dir) const {
      return ci.numCells() * sizes[dir] * sizeof(value_type);
   }

   /**
    * Returns the number of bytes that will be packed from / unpacked to the cell interval
    * when using packAll / unpackAll
    * @param ci  The cell interval
    * @return    The required size of the buffer, in bytes
    */
   uint_t size (CellInterval & ci) const {
      return ci.numCells() * {{stencil_size}} * sizeof(value_type);
   }

   {% block AdditionalPublicDeclarations %}
   {% endblock %}

 private:
   const uint_t sizes[{{direction_sizes|length}}] { {{ direction_sizes | join(', ') }} };
};

}  // namespace {{namespace}}
}  // namespace walberla
