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
//! \\file D3Q27StorageSpecification.h
//! \\author lbmpy
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"

#include "domain_decomposition/IBlock.h"
#include "field/GhostLayerField.h"

#include "stencil/D3Q27.h"
#include "stencil/Directions.h"

#define FUNC_PREFIX

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
namespace lbm{

class D3Q27StorageSpecification
{
 public:
   // Used lattice stencil
   using Stencil = stencil::D3Q27;
   // Lattice stencil used for the communication (should be used to define which block directions need to be communicated)
   using CommunicationStencil = stencil::D3Q27;
   // If false used correction: Lattice Boltzmann Model for the Incompressible Navier–Stokes Equation, He 1997
   static const bool compressible = false;
   // Cut off for the lattice Boltzmann equilibrium
   static const int equilibriumAccuracyOrder = 2;
   // If true the equilibrium is computed in regard to "delta_rho" and not the actual density "rho"
   static const bool equilibriumDeviationOnly = true;
   // If streaming pattern is inplace (esotwist, aa, ...) or not (pull, push)
   static const bool inplace = false;
   // If true the background deviation (rho_0 = 1) is subtracted for the collision step.
   static const bool zeroCenteredPDFs = true;
   // Lattice weights
   static constexpr double w[27] = { double(0.296296296296296), double(0.0740740740740741), double(0.0740740740740741), double(0.0740740740740741), double(0.0740740740740741), double(0.0740740740740741), double(0.0740740740740741), double(0.0185185185185185), double(0.0185185185185185), double(0.0185185185185185), double(0.0185185185185185), double(0.0185185185185185), double(0.0185185185185185), double(0.0185185185185185), double(0.0185185185185185), double(0.0185185185185185), double(0.0185185185185185), double(0.0185185185185185), double(0.0185185185185185), double(0.00462962962962963), double(0.00462962962962963), double(0.00462962962962963), double(0.00462962962962963), double(0.00462962962962963), double(0.00462962962962963), double(0.00462962962962963), double(0.00462962962962963) };
   // Inverse lattice weights
   static constexpr double wInv[27] = { double(3.37500000000000), double(13.5000000000000), double(13.5000000000000), double(13.5000000000000), double(13.5000000000000), double(13.5000000000000), double(13.5000000000000), double(54.0000000000000), double(54.0000000000000), double(54.0000000000000), double(54.0000000000000), double(54.0000000000000), double(54.0000000000000), double(54.0000000000000), double(54.0000000000000), double(54.0000000000000), double(54.0000000000000), double(54.0000000000000), double(54.0000000000000), double(216.000000000000), double(216.000000000000), double(216.000000000000), double(216.000000000000), double(216.000000000000), double(216.000000000000), double(216.000000000000), double(216.000000000000) };

   struct AccessorEVEN
   {
      static constexpr cell_idx_t readX[27] = { 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1 };
      static constexpr cell_idx_t readY[27] = { 0, -1, 1, 0, 0, 0, 0, -1, -1, 1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, -1, 1, 1, -1, -1, 1, 1 };
      static constexpr cell_idx_t readZ[27] = { 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1 };
      static constexpr cell_idx_t readD[27] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 };

      static constexpr cell_idx_t writeX[27] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      static constexpr cell_idx_t writeY[27] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      static constexpr cell_idx_t writeZ[27] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      static constexpr cell_idx_t writeD[27] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 };
   };

   struct AccessorODD
   {
      static constexpr cell_idx_t readX[27] = { 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1 };
      static constexpr cell_idx_t readY[27] = { 0, -1, 1, 0, 0, 0, 0, -1, -1, 1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, -1, 1, 1, -1, -1, 1, 1 };
      static constexpr cell_idx_t readZ[27] = { 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1 };
      static constexpr cell_idx_t readD[27] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 };

      static constexpr cell_idx_t writeX[27] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      static constexpr cell_idx_t writeY[27] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      static constexpr cell_idx_t writeZ[27] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      static constexpr cell_idx_t writeD[27] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 };
   };

   // Compute kernels to pack and unpack MPI buffers
   class PackKernels {

    public:
      using PdfField_T = field::GhostLayerField<double, 27>;
      using value_type = typename PdfField_T::value_type;

      

      static const bool inplace = false;

      /**
      * Packs all pdfs from the given cell interval to the send buffer.
      * */
      void packAll(PdfField_T * pdfs_src, CellInterval & ci, unsigned char * outBuffer) const;

      /**
       * Unpacks all pdfs from the send buffer to the given cell interval.
       * */
      void unpackAll(PdfField_T * pdfs_dst, CellInterval & ci, unsigned char * inBuffer) const;

      /**
       * Copies data between two blocks on the same process.
       * All pdfs from the sending interval are copied onto the receiving interval.
       * */
      void localCopyAll(PdfField_T * pdfs_src, CellInterval & srcInterval, PdfField_T * pdfs_dst, CellInterval & dstInterval) const;

      /**
       * Packs only those populations streaming in directions aligned with the sending direction dir from the given cell interval.
       * For example, in 2D, if dir == N, the pdfs streaming in directions NW, N, NE are packed.
       * */
      void packDirection(PdfField_T * pdfs_src, CellInterval & ci, unsigned char * outBuffer, stencil::Direction dir) const;

      /**
       * Unpacks only those populations streaming in directions aligned with the sending direction dir to the given cell interval.
       * For example, in 2D, if dir == N, the pdfs streaming in directions NW, N, NE are unpacked.
       * */
      void unpackDirection(PdfField_T * pdfs_dst, CellInterval & ci, unsigned char * inBuffer, stencil::Direction dir) const;

      /** Copies data between two blocks on the same process.
        * PDFs streaming aligned with the direction dir are copied from the sending interval onto the receiving interval.
        * */
      void localCopyDirection(PdfField_T * pdfs_src, CellInterval & srcInterval, PdfField_T * pdfs_dst, CellInterval & dstInterval, stencil::Direction dir) const;

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
         return ci.numCells() * 27 * uint_c(sizeof(value_type));
      }

      

    private:
      const uint_t sizes[27] { 0, 9, 9, 9, 9, 9, 9, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1 };
   };

   using value_type = PackKernels::value_type;

};

}} //lbm/walberla