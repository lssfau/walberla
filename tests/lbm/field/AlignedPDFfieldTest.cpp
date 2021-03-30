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
//! \file AlignedPDFfieldTest.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/all.h"

#include "field/allocation/FieldAllocator.h"

#include "lbm/all.h"

namespace walberla
{
using LatticeModel_T = lbm::D3Q19< lbm::collision_model::SRT, false >;
typedef lbm::PdfField< LatticeModel_T > PdfField_T;

int main(int argc, char** argv)
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv(argc, argv);

   auto blocks = blockforest::createUniformBlockGridFromConfig(walberlaEnv.config());

   auto parameters    = walberlaEnv.config()->getOneBlock("Parameters");
   const real_t omega = parameters.getParameter< real_t >("omega", real_c(1.4));

   LatticeModel_T latticeModel = LatticeModel_T(omega);

   shared_ptr< field::FieldAllocator< real_t > > alloc_32 = make_shared< field::AllocateAligned< real_t, 32 > >();
   BlockDataID pdfFieldId_32 = lbm::addPdfFieldToStorage(blocks, "pdf field", latticeModel, field::fzyx, Set<SUID>::emptySet(), Set<SUID>::emptySet(), alloc_32);

   shared_ptr< field::FieldAllocator< real_t > > alloc_64 = make_shared< field::AllocateAligned< real_t, 64 > >();
   BlockDataID pdfFieldId_64 = lbm::addPdfFieldToStorage(blocks, "pdf field", latticeModel, field::fzyx, Set<SUID>::emptySet(), Set<SUID>::emptySet(), alloc_64);

   shared_ptr< field::FieldAllocator< real_t > > alloc_128 = make_shared< field::AllocateAligned< real_t, 128 > >();
   BlockDataID pdfFieldId_128 = lbm::addPdfFieldToStorage(blocks, "pdf field", latticeModel, field::fzyx, Set<SUID>::emptySet(), Set<SUID>::emptySet(), alloc_128);

   for (auto& block : *blocks)
   {
      auto pdfField_32 = block.getData< PdfField_T >(pdfFieldId_32);
      WALBERLA_CHECK_EQUAL((size_t) pdfField_32->dataAt(0, 0, 0, 0) % 32, 0)
      WALBERLA_CHECK_EQUAL((size_t) pdfField_32->dataAt(0, 0, 0, 1) % 32, 0)
      WALBERLA_CHECK_EQUAL((size_t) pdfField_32->dataAt(0, 0, 1, 0) % 32, 0)
      WALBERLA_CHECK_EQUAL((size_t) pdfField_32->dataAt(0, 1, 0, 0) % 32, 0)

      auto pdfField_64 = block.getData< PdfField_T >(pdfFieldId_64);
      WALBERLA_CHECK_EQUAL((size_t) pdfField_64->dataAt(0, 0, 0, 0) % 64, 0)
      WALBERLA_CHECK_EQUAL((size_t) pdfField_64->dataAt(0, 0, 0, 1) % 64, 0)
      WALBERLA_CHECK_EQUAL((size_t) pdfField_64->dataAt(0, 0, 1, 0) % 64, 0)
      WALBERLA_CHECK_EQUAL((size_t) pdfField_64->dataAt(0, 1, 0, 0) % 64, 0)

      auto pdfField_128 = block.getData< PdfField_T >(pdfFieldId_128);
      WALBERLA_CHECK_EQUAL((size_t) pdfField_128->dataAt(0, 0, 0, 0) % 128, 0)
      WALBERLA_CHECK_EQUAL((size_t) pdfField_128->dataAt(0, 0, 0, 1) % 128, 0)
      WALBERLA_CHECK_EQUAL((size_t) pdfField_128->dataAt(0, 0, 1, 0) % 128, 0)
      WALBERLA_CHECK_EQUAL((size_t) pdfField_128->dataAt(0, 1, 0, 0) % 128, 0)
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }