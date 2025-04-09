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
#include "core/DataTypes.h"
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

   class {{class_name}}
   {
    public:
      using Field_T    = {{src_field | field_type(is_gpu=is_gpu)}};
      using value_type = typename Field_T::value_type;

      void packEqual({{-[
         "Field_T * " + src_field.name, "CellInterval & ci", "unsigned char * outBuffer",
         kernels['packEqual'].kernel_selection_parameters, ["gpuStream_t stream = nullptr"] if is_gpu else[]]
                         | type_identifier_list -}}) const;

      void unpackEqual({{-[
         "Field_T * " + dst_field.name, "CellInterval & ci", "unsigned char * inBuffer",
         kernels['unpackEqual'].kernel_selection_parameters, ["gpuStream_t stream = nullptr"] if is_gpu else[]
      ] | type_identifier_list -}}) const;

      void localCopyEqual(
         {{- [ "Field_T * " + src_field.name, "CellInterval & srcInterval",
                "Field_T * " + dst_field.name, "CellInterval & dstInterval",
                kernels['localCopyEqual'].kernel_selection_parameters,
                ["gpuStream_t stream = nullptr"] if is_gpu else []]
             | type_identifier_list -}}) const;

      void packFineToCoarse({{-[
         "Field_T * " + src_field.name, "CellInterval & ci", "unsigned char * outBuffer",
         kernels['packFineToCoarse'].kernel_selection_parameters,
         ["gpuStream_t stream = nullptr"] if is_gpu else[]
      ] | type_identifier_list -}}) const;

      void unpackCoarseToFine({{-[
         "Field_T * " + dst_field.name, "CellInterval & ci", "unsigned char * inBuffer",
         kernels['unpackCoarseToFine'].kernel_selection_parameters, ["gpuStream_t stream = nullptr"] if is_gpu else[]
      ] | type_identifier_list -}}) const;

   };
} // end namespace
} // end walberla