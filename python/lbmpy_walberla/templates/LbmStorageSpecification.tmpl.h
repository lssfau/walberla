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

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"

#include "domain_decomposition/IBlock.h"
#include "field/GhostLayerField.h"

#include "stencil/{{stencil_name}}.h"
#include "stencil/Directions.h"

{% if target is equalto 'cpu' -%}
#define FUNC_PREFIX
{%- elif target is equalto 'gpu' -%}
#define FUNC_PREFIX __global__
#include "gpu/GPUWrapper.h"
#include "gpu/GPUField.h"
{%- endif %}

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#if defined WALBERLA_CXX_COMPILER_IS_GNU || defined WALBERLA_CXX_COMPILER_IS_CLANG
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

namespace walberla
{
namespace {{namespace}}{

class {{class_name}}
{
 public:
   // Used lattice stencil
   using Stencil = stencil::{{stencil_name}};
   // Lattice stencil used for the communication (should be used to define which block directions need to be communicated)
   using CommunicationStencil = stencil::{{communication_stencil_name}};
   // If false used correction: Lattice Boltzmann Model for the Incompressible Navier–Stokes Equation, He 1997
   static const bool compressible = {% if compressible %}true{% else %}false{% endif %};
   // Cut off for the lattice Boltzmann equilibrium
   static const int equilibriumAccuracyOrder = {{equilibrium_accuracy_order}};
   // If true the equilibrium is computed in regard to "delta_rho" and not the actual density "rho"
   static const bool equilibriumDeviationOnly = {% if equilibrium_deviation_only -%} true {%- else -%} false {%- endif -%};
   // If streaming pattern is inplace (esotwist, aa, ...) or not (pull, push)
   static const bool inplace = {% if inplace -%} true {%- else -%} false {%- endif -%};
   // If true the background deviation (rho_0 = 1) is subtracted for the collision step.
   static const bool zeroCenteredPDFs = {% if zero_centered -%} true {%- else -%} false {%- endif -%};
   // Lattice weights
   static constexpr {{dtype}} w[{{stencil_size}}] = { {{weights}} };
   // Inverse lattice weights
   static constexpr {{dtype}} wInv[{{stencil_size}}] = { {{inverse_weights}} };

   struct AccessorEVEN
   {
      static constexpr cell_idx_t readX[{{stencil_size}}] = { {{even_read[0]}} };
      static constexpr cell_idx_t readY[{{stencil_size}}] = { {{even_read[1]}} };
      static constexpr cell_idx_t readZ[{{stencil_size}}] = { {{even_read[2]}} };
      static constexpr cell_idx_t readD[{{stencil_size}}] = { {{even_read[3]}} };

      static constexpr cell_idx_t writeX[{{stencil_size}}] = { {{even_write[0]}} };
      static constexpr cell_idx_t writeY[{{stencil_size}}] = { {{even_write[1]}} };
      static constexpr cell_idx_t writeZ[{{stencil_size}}] = { {{even_write[2]}} };
      static constexpr cell_idx_t writeD[{{stencil_size}}] = { {{even_write[3]}} };
   };

   struct AccessorODD
   {
      static constexpr cell_idx_t readX[{{stencil_size}}] = { {{odd_read[0]}} };
      static constexpr cell_idx_t readY[{{stencil_size}}] = { {{odd_read[1]}} };
      static constexpr cell_idx_t readZ[{{stencil_size}}] = { {{odd_read[2]}} };
      static constexpr cell_idx_t readD[{{stencil_size}}] = { {{odd_read[3]}} };

      static constexpr cell_idx_t writeX[{{stencil_size}}] = { {{odd_write[0]}} };
      static constexpr cell_idx_t writeY[{{stencil_size}}] = { {{odd_write[1]}} };
      static constexpr cell_idx_t writeZ[{{stencil_size}}] = { {{odd_write[2]}} };
      static constexpr cell_idx_t writeD[{{stencil_size}}] = { {{odd_write[3]}} };
   };

   // Compute kernels to pack and unpack MPI buffers
   class PackKernels {

    public:
      using PdfField_T = {{src_field | field_type(is_gpu=is_gpu)}};
      using value_type = typename PdfField_T::value_type;

      {% if nonuniform -%}
      {% if target is equalto 'cpu' -%}
      using MaskField_T = GhostLayerField< uint32_t, 1 >;
      {%- elif target is equalto 'gpu' -%}
      using MaskField_T = gpu::GPUField< uint32_t >;
      {%- endif %}
      {%- endif %}

      static const bool inplace = {% if inplace -%} true {%- else -%} false {%- endif -%};

      /**
      * Packs all pdfs from the given cell interval to the send buffer.
      * */
      void packAll(
         {{- [ "PdfField_T * " + src_field.name, "CellInterval & ci",
                "unsigned char * outBuffer", kernels['packAll'].kernel_selection_parameters,
                ["gpuStream_t stream = nullptr"] if is_gpu else []]
             | type_identifier_list -}}
      ) const;

      /**
       * Unpacks all pdfs from the send buffer to the given cell interval.
       * */
      void unpackAll(
         {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
                "unsigned char * inBuffer", kernels['unpackAll'].kernel_selection_parameters,
                ["gpuStream_t stream = nullptr"] if is_gpu else []]
             | type_identifier_list -}}
      ) const;

      /**
       * Copies data between two blocks on the same process.
       * All pdfs from the sending interval are copied onto the receiving interval.
       * */
      void localCopyAll(
         {{- [ "PdfField_T * " + src_field.name, "CellInterval & srcInterval",
                "PdfField_T * " + dst_field.name, "CellInterval & dstInterval",
                kernels['localCopyAll'].kernel_selection_parameters,
                ["gpuStream_t stream = nullptr"] if is_gpu else []]
             | type_identifier_list -}}
      ) const;

      /**
       * Packs only those populations streaming in directions aligned with the sending direction dir from the given cell interval.
       * For example, in 2D, if dir == N, the pdfs streaming in directions NW, N, NE are packed.
       * */
      void packDirection(
         {{- [ "PdfField_T * " + src_field.name, "CellInterval & ci",
                "unsigned char * outBuffer", kernels['packDirection'].kernel_selection_parameters,
                ["gpuStream_t stream = nullptr"] if is_gpu else []]
             | type_identifier_list -}}
      ) const;

      /**
       * Unpacks only those populations streaming in directions aligned with the sending direction dir to the given cell interval.
       * For example, in 2D, if dir == N, the pdfs streaming in directions NW, N, NE are unpacked.
       * */
      void unpackDirection(
         {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
                "unsigned char * inBuffer", kernels['unpackDirection'].kernel_selection_parameters,
                ["gpuStream_t stream = nullptr"] if is_gpu else []]
             | type_identifier_list -}}
      ) const;

      /** Copies data between two blocks on the same process.
        * PDFs streaming aligned with the direction dir are copied from the sending interval onto the receiving interval.
        * */
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
       * */
      uint_t size (CellInterval & ci, stencil::Direction dir) const {
         return ci.numCells() * sizes[dir] * uint_c(sizeof(value_type));
      }

      /**
       * Returns the number of bytes that will be packed from / unpacked to the cell interval
       * when using packAll / unpackAll
       * @param ci  The cell interval
       * @return    The required size of the buffer, in bytes
       * */
      uint_t size (CellInterval & ci) const {
         return ci.numCells() * {{stencil_size}} * uint_c(sizeof(value_type));
      }

      {% if nonuniform -%}

      /**
       * Unpacks and uniformly redistributes populations coming from a coarse block onto the fine grid.
       * */
      void unpackRedistribute(
         {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
                "unsigned char * inBuffer", kernels['unpackRedistribute'].kernel_selection_parameters,
                ["gpuStream_t stream = nullptr"] if is_gpu else []]
             | type_identifier_list -}}
      ) const;

      /**
       * Partially coalesces and packs populations streaming from a fine block into a coarse block
       * */
      void packPartialCoalescence(
         {{- [ "PdfField_T * " + src_field.name, "MaskField_T * " + mask_field.name, "CellInterval & ci",
                "unsigned char * outBuffer", kernels['packPartialCoalescence'].kernel_selection_parameters,
                ["gpuStream_t stream = nullptr"] if is_gpu else []]
             | type_identifier_list -}}
      ) const;

      /**
       * Prepares a coarse block for coalescence by setting every population that must be coalesced from fine blocks to zero.
       * */
      void zeroCoalescenceRegion(
         {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
                kernels['zeroCoalescenceRegion'].kernel_selection_parameters,
                ["gpuStream_t stream = nullptr"] if is_gpu else []]
             | type_identifier_list -}}
      ) const;

      /**
       * Unpacks and coalesces populations coming from a fine block onto the fine grid
       * */
      void unpackCoalescence(
         {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
                "unsigned char * inBuffer", kernels['unpackCoalescence'].kernel_selection_parameters,
                ["gpuStream_t stream = nullptr"] if is_gpu else []]
             | type_identifier_list -}}
      ) const;

      /**
       * Returns the number of bytes that will be unpacked to the cell interval
       * when using unpackRedistribute. This is 2^{-d} of the data that would be
       * unpacked during same-level communication.
       * @param ci  The cell interval
       * @return    The required size of the buffer, in bytes
       * */
      uint_t redistributeSize(CellInterval & ci) const {
         return size(ci) >> {{dimension}};
      }

      /**
       * Returns the number of bytes that will be packed from the cell interval
       * when using packPartialCoalescence.
       * @param ci  The cell interval
       * @param dir The communication direction
       * @return    The required size of the buffer, in bytes
       * */
      uint_t partialCoalescenceSize(CellInterval & ci, stencil::Direction dir) const {
         return size(ci, dir) >> {{dimension}};
      }

      {%- endif %}

    private:
      const uint_t sizes[{{direction_sizes|length}}] { {{ direction_sizes | join(', ') }} };
   };

   using value_type = PackKernels::value_type;

};

}} //{{namespace}}/walberla