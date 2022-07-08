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
//! \file PdfRefillingModel.h
//! \ingroup dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \author Michael Zikeli
//! \brief Defines how cells are refilled (i.e. PDFs reinitialized) after the cell was converted from gas to interface.
//
//======================================================================================================================

#pragma once

#include "core/StringUtility.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 *   Class that specifies how PDFs are reinitialized in cells that are converted from gas to interface.
 *
 *   Available options are:
 *       - EquilibriumRefilling:
 *              initialize PDFs with equilibrium using average density and velocity from neighboring cells; default
 *              approach used in any known publication with free surface LBM
 *
 *       - AverageRefilling:
 *              initialize PDFs with average PDFs (in the respective directions) of neighboring cells
 *
 *       - EquilibriumAndNonEquilibriumRefilling:
 *              initialize PDFs with EquilibriumRefilling and add the non-equilibrium contribution of neighboring cells
 *
 *       - ExtrapolationRefilling:
 *              initialize PDFs with PDFs extrapolated (in surface normal direction) from neighboring cells
 *
 *       - GradsMomentsRefilling:
 *              initialize PDFs with EquilibriumRefilling and add the contribution of the non-equilibrium pressure
 *              tensor
 *
 * See src/lbm/free_surface/dynamics/PdfRefillingSweep.h for a detailed description of the models.
 *
 * The models and their implementation are inspired by the equivalent functionality of the lbm-particle coupling, see
 * src/lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h.
 **********************************************************************************************************************/
class PdfRefillingModel
{
 public:
   enum class RefillingModel {
      EquilibriumRefilling,
      AverageRefilling,
      EquilibriumAndNonEquilibriumRefilling,
      ExtrapolationRefilling,
      GradsMomentsRefilling
   };

   PdfRefillingModel(const std::string& modelName) : modelName_(modelName), modelType_(chooseType(modelName)) {}

   PdfRefillingModel(const RefillingModel& modelType) : modelName_(chooseName(modelType)), modelType_(modelType)
   {
      switch (modelType_)
      {
      case RefillingModel::EquilibriumRefilling:
         break;
      case RefillingModel::AverageRefilling:
         break;
      case RefillingModel::EquilibriumAndNonEquilibriumRefilling:
         break;
      case RefillingModel::ExtrapolationRefilling:
         break;
      case RefillingModel::GradsMomentsRefilling:
         break;
      }
   }

   inline RefillingModel getModelType() const { return modelType_; }
   inline std::string getModelName() const { return modelName_; }
   inline std::string getFullModelSpecification() const { return getModelName(); }

   static inline std::initializer_list< const RefillingModel > getTypeIterator() { return listOfAllEnums; }

 private:
   RefillingModel chooseType(const std::string& modelName)
   {
      if (!string_icompare(modelName, "EquilibriumRefilling")) { return RefillingModel::EquilibriumRefilling; }

      if (!string_icompare(modelName, "AverageRefilling")) { return RefillingModel::AverageRefilling; }

      if (!string_icompare(modelName, "EquilibriumAndNonEquilibriumRefilling"))
      {
         return RefillingModel::EquilibriumAndNonEquilibriumRefilling;
      }

      if (!string_icompare(modelName, "ExtrapolationRefilling")) { return RefillingModel::ExtrapolationRefilling; }

      if (!string_icompare(modelName, "GradsMomentsRefilling")) { return RefillingModel::GradsMomentsRefilling; }

      WALBERLA_ABORT("The specified PDF reinitialization model " << modelName << " is not available.");
   }

   std::string chooseName(RefillingModel const& modelType) const
   {
      std::string modelName;
      switch (modelType)
      {
      case RefillingModel::EquilibriumRefilling:
         modelName = "EquilibriumRefilling";
         break;
      case RefillingModel::AverageRefilling:
         modelName = "AverageRefilling";
         break;
      case RefillingModel::EquilibriumAndNonEquilibriumRefilling:
         modelName = "EquilibriumAndNonEquilibriumRefilling";
         break;
      case RefillingModel::ExtrapolationRefilling:
         modelName = "ExtrapolationRefilling";
         break;
      case RefillingModel::GradsMomentsRefilling:
         modelName = "GradsMomentsRefilling";
         break;
      }
      return modelName;
   }

   std::string modelName_;
   RefillingModel modelType_;
   static constexpr std::initializer_list< const RefillingModel > listOfAllEnums = {
      RefillingModel::EquilibriumRefilling, RefillingModel::AverageRefilling,
      RefillingModel::EquilibriumAndNonEquilibriumRefilling, RefillingModel::ExtrapolationRefilling,
      RefillingModel::GradsMomentsRefilling
   };

}; // class PdfRefillingModel
} // namespace free_surface
} // namespace walberla