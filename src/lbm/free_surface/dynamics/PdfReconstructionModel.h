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
//! \file PdfReconstructionModel.h
//! \ingroup dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Class that specifies the number of reconstructed PDFs at the free surface interface.
//
//======================================================================================================================

#pragma once

#include "core/StringUtility.h"
#include "core/stringToNum.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Class that specifies the number of reconstructed PDFs at the free surface interface. PDFs need to be reconstructed
 * as they might be missing, i.e., PDFs streaming from gas to interface are not available inherently.
 *
 * Available models:
 * - NormalBasedKeepCenter: reconstruct all PDFs for which n * c_i >= 0 (approach by Koerner et al., 2005); some
 *                           already available PDFs will be overwritten
 *
 * - NormalBasedReconstructCenter: reconstruct all PDFs for which n * c_i >= 0 (including the center PDF); some already
 *                                 available PDFs coming from liquid will be overwritten
 *
 * - OnlyMissing: reconstruct only missing PDFs (no already available PDF gets overwritten)
 *
 * - All: reconstruct all PDFs (any already available PDF is overwritten)
 *
 * - OnlyMissingMin-N-largest: Reconstruct only missing PDFs but at least N (and therefore potentially overwrite
 *                             available PDFs); "smallest" or "largest" specifies whether PDFs with smallest or largest
 *                             n * c_i get overwritten first. This model is motivated by the dissertation of Simon
 *                             Bogner, 2017, section 4.2.1, where it is argued that at least 3 PDFs must be
 *                             reconstructed, as otherwise the free surface boundary condition is claimed to be
 *                             underdetermined. However, a mathematical proof for this statement is not given.
 *
 * - OnlyMissingMin-N-smallest: see comment at "OnlyMissingMin-N-largest"
 *
 * - OnlyMissingMin-N-normalBasedKeepCenter: see comment at "OnlyMissingMin-N-largest"; if less than N PDFs are unknown,
 *                                           reconstruct according to the model "NormalBasedKeepCenter"
 * ********************************************************************************************************************/
class PdfReconstructionModel
{
 public:
   enum class ReconstructionModel {
      NormalBasedReconstructCenter,
      NormalBasedKeepCenter,
      OnlyMissing,
      All,
      OnlyMissingMin,
   };

   enum class FallbackModel {
      Largest,
      Smallest,
      NormalBasedKeepCenter,
   };

   PdfReconstructionModel(const std::string& modelName) : modelName_(modelName), modelType_(chooseType(modelName))
   {
      if (modelType_ == ReconstructionModel::OnlyMissingMin)
      {
         const std::vector< std::string > substrings = string_split(modelName, "-");

         modelName_         = substrings[0];                        // "OnlyMissingMin"
         numMinReconstruct_ = stringToNum< uint_t >(substrings[1]); // N
         fallbackModelName_ = substrings[2]; // "smallest" or "largest" or "normalBasedKeepCenter"
         fallbackModel_     = chooseFallbackModel(fallbackModelName_);

         if (fallbackModel_ != FallbackModel::Largest && fallbackModel_ != FallbackModel::Smallest &&
             fallbackModel_ != FallbackModel::NormalBasedKeepCenter)
         {
            WALBERLA_ABORT("The specified PDF reconstruction fallback-model " << modelName << " is not available.");
         }
      }
   }

   inline ReconstructionModel getModelType() const { return modelType_; }
   inline std::string getModelName() const { return modelName_; }
   inline uint_t getNumMinReconstruct() const { return numMinReconstruct_; }
   inline FallbackModel getFallbackModel() const { return fallbackModel_; }
   inline std::string getFallbackModelName() const { return fallbackModelName_; }
   inline std::string getFullModelSpecification() const
   {
      if (modelType_ == ReconstructionModel::OnlyMissingMin)
      {
         return modelName_ + "-" + std::to_string(numMinReconstruct_) + "-" + fallbackModelName_;
      }
      else { return modelName_; }
   }

 private:
   ReconstructionModel chooseType(const std::string& modelName)
   {
      if (!string_icompare(modelName, "NormalBasedReconstructCenter"))
      {
         return ReconstructionModel::NormalBasedReconstructCenter;
      }

      else
      {
         if (!string_icompare(modelName, "NormalBasedKeepCenter"))
         {
            return ReconstructionModel::NormalBasedKeepCenter;
         }
         else
         {
            if (!string_icompare(modelName, "OnlyMissing")) { return ReconstructionModel::OnlyMissing; }
            else
            {
               if (!string_icompare(modelName, "All")) { return ReconstructionModel::All; }
               else
               {
                  if (!string_icompare(string_split(modelName, "-")[0], "OnlyMissingMin"))
                  {
                     return ReconstructionModel::OnlyMissingMin;
                  }
                  else
                  {
                     WALBERLA_ABORT("The specified PDF reconstruction model " << modelName << " is not available.");
                  }
               }
            }
         }
      }
   }

   FallbackModel chooseFallbackModel(const std::string& fallbackModelName)
   {
      if (!string_icompare(fallbackModelName, "largest")) { return FallbackModel::Largest; }
      else
      {
         if (!string_icompare(fallbackModelName, "smallest")) { return FallbackModel::Smallest; }
         else
         {
            if (!string_icompare(fallbackModelName, "normalBasedKeepCenter"))
            {
               return FallbackModel::NormalBasedKeepCenter;
            }
            else
            {
               WALBERLA_ABORT("The specified PDF reconstruction fallback-model " << fallbackModelName
                                                                                 << " is not available.");
            }
         }
      }
   }

   std::string modelName_;
   ReconstructionModel modelType_;
   uint_t numMinReconstruct_;
   std::string fallbackModelName_;
   FallbackModel fallbackModel_;
}; // class PdfReconstructionModel
} // namespace free_surface
} // namespace walberla