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
//! \file ExcessMassDistributionModel.h
//! \ingroup dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Class that specifies how excessive mass is distributed.
//
//======================================================================================================================

#pragma once

#include "core/StringUtility.h"
#include "core/stringToNum.h"

#include <string>

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Class that specifies how excessive mass is distributed after cell conversions from interface to liquid or interface
 * to gas.
 * For example, when converting an interface cell with fill level 1.1 to liquid with fill level 1,0, an excessive mass
 * corresponding to the fill level 0.1 must be distributed.
 *
 * Available models:
 *  - EvenlyAllInterface:
 *       Excess mass is distributed evenly among all neighboring interface cells (see dissertations of T. Pohl, S.
 *       Donath, S. Bogner).
 *
 *  - EvenlyNewInterface:
 *       Excess mass is distributed evenly among newly converted neighboring interface cells (see Koerner et al., 2005).
 *       Falls back to EvenlyAllInterface if not applicable.
 *
 *  - EvenlyOldInterface:
 *       Excess mass is distributed evenly among old neighboring interface cells, i.e., cells that are non-newly
 *       converted to interface. Falls back to EvenlyAllInterface if not applicable.
 *
 *  - WeightedAllInterface:
 *       Excess mass is distributed weighted with the direction of the interface normal among all neighboring interface
 *       cells (see dissertation of N. Thuerey, 2007). Falls back to EvenlyAllInterface if not applicable.
 *
 *  - WeightedNewInterface:
 *       Excess mass is distributed weighted with the direction of the interface normal among newly converted
 *       neighboring interface cells. Falls back to WeightedAllInterface if not applicable.
 *
 *  - WeightedOldInterface:
 *       Excess mass is distributed weighted with the direction of the interface normal among old neighboring interface
 *       cells, i.e., cells that are non-newly converted to interface. Falls back to WeightedAllInterface if not
 * applicable.
 *
 *  - EvenlyLiquidAndAllInterface:
 *      Excess mass is distributed evenly among all neighboring interface and liquid cells (see p.47 in master thesis of
 *      M. Lehmann, 2019). The excess mass distributed to liquid cells does neither modify the cell's density nor fill
 *      level. Instead, it is stored in an additional excess mass field. Therefore, not only the converted interface
 *      cells' excess mass is distributed, but also the excess mass of liquid cells stored in this additional field.
 *
 *  - EvenlyLiquidAndAllInterfacePreferInterface:
 *      Similar to EvenlyLiquidAndAllInterface, however, excess mass is preferably distributed to interface cells. It is
 *      distributed to liquid cells only if there are no neighboring interface cells available.
 * ********************************************************************************************************************/
class ExcessMassDistributionModel
{
 public:
   enum class ExcessMassModel {
      EvenlyAllInterface,
      EvenlyNewInterface,
      EvenlyOldInterface,
      WeightedAllInterface,
      WeightedNewInterface,
      WeightedOldInterface,
      EvenlyLiquidAndAllInterface,
      EvenlyLiquidAndAllInterfacePreferInterface
   };

   ExcessMassDistributionModel(const std::string& modelName) : modelName_(modelName), modelType_(chooseType(modelName))
   {}

   ExcessMassDistributionModel(const ExcessMassModel& modelType)
      : modelName_(chooseName(modelType)), modelType_(modelType)
   {
      switch (modelType_)
      {
      case ExcessMassModel::EvenlyAllInterface:
         break;
      case ExcessMassModel::EvenlyNewInterface:
         break;
      case ExcessMassModel::EvenlyOldInterface:
         break;
      case ExcessMassModel::WeightedAllInterface:
         break;
      case ExcessMassModel::WeightedNewInterface:
         break;
      case ExcessMassModel::WeightedOldInterface:
         break;
      case ExcessMassModel::EvenlyLiquidAndAllInterface:
         break;
      case ExcessMassModel::EvenlyLiquidAndAllInterfacePreferInterface:
         break;
      }
   }

   inline ExcessMassModel getModelType() const { return modelType_; }
   inline std::string getModelName() const { return modelName_; }
   inline std::string getFullModelSpecification() const { return getModelName(); }

   inline bool isEvenlyType() const
   {
      return modelType_ == ExcessMassModel::EvenlyAllInterface || modelType_ == ExcessMassModel::EvenlyNewInterface ||
             modelType_ == ExcessMassModel::EvenlyOldInterface;
   }

   inline bool isWeightedType() const
   {
      return modelType_ == ExcessMassModel::WeightedAllInterface ||
             modelType_ == ExcessMassModel::WeightedNewInterface || modelType_ == ExcessMassModel::WeightedOldInterface;
   }

   inline bool isEvenlyLiquidAndAllInterfacePreferInterfaceType() const
   {
      return modelType_ == ExcessMassModel::EvenlyLiquidAndAllInterface ||
             modelType_ == ExcessMassModel::EvenlyLiquidAndAllInterfacePreferInterface;
   }

   static inline std::initializer_list< const ExcessMassModel > getTypeIterator() { return listOfAllEnums; }

 private:
   ExcessMassModel chooseType(const std::string& modelName)
   {
      if (!string_icompare(modelName, "EvenlyAllInterface")) { return ExcessMassModel::EvenlyAllInterface; }

      if (!string_icompare(modelName, "EvenlyNewInterface")) { return ExcessMassModel::EvenlyNewInterface; }

      if (!string_icompare(modelName, "EvenlyOldInterface")) { return ExcessMassModel::EvenlyOldInterface; }

      if (!string_icompare(modelName, "WeightedAllInterface")) { return ExcessMassModel::WeightedAllInterface; }

      if (!string_icompare(modelName, "WeightedNewInterface")) { return ExcessMassModel::WeightedNewInterface; }

      if (!string_icompare(modelName, "WeightedOldInterface")) { return ExcessMassModel::WeightedOldInterface; }

      if (!string_icompare(modelName, "EvenlyLiquidAndAllInterface"))
      {
         return ExcessMassModel::EvenlyLiquidAndAllInterface;
      }

      if (!string_icompare(modelName, "EvenlyLiquidAndAllInterfacePreferInterface"))
      {
         return ExcessMassModel::EvenlyLiquidAndAllInterfacePreferInterface;
      }

      WALBERLA_ABORT("The specified PDF reinitialization model " << modelName << " is not available.");
   }

   std::string chooseName(ExcessMassModel const& modelType) const
   {
      std::string modelName;
      switch (modelType)
      {
      case ExcessMassModel::EvenlyAllInterface:
         modelName = "EvenlyAllInterface";
         break;
      case ExcessMassModel::EvenlyNewInterface:
         modelName = "EvenlyNewInterface";
         break;
      case ExcessMassModel::EvenlyOldInterface:
         modelName = "EvenlyOldInterface";
         break;
      case ExcessMassModel::WeightedAllInterface:
         modelName = "WeightedAllInterface";
         break;
      case ExcessMassModel::WeightedNewInterface:
         modelName = "WeightedNewInterface";
         break;
      case ExcessMassModel::WeightedOldInterface:
         modelName = "WeightedOldInterface";
         break;

      case ExcessMassModel::EvenlyLiquidAndAllInterface:
         modelName = "EvenlyLiquidAndAllInterface";
         break;
      case ExcessMassModel::EvenlyLiquidAndAllInterfacePreferInterface:
         modelName = "EvenlyLiquidAndAllInterfacePreferInterface";
         break;
      }
      return modelName;
   }

   std::string modelName_;
   ExcessMassModel modelType_;
   static constexpr std::initializer_list< const ExcessMassModel > listOfAllEnums = {
      ExcessMassModel::EvenlyAllInterface,          ExcessMassModel::EvenlyNewInterface,
      ExcessMassModel::EvenlyOldInterface,          ExcessMassModel::WeightedAllInterface,
      ExcessMassModel::WeightedNewInterface,        ExcessMassModel::WeightedOldInterface,
      ExcessMassModel::EvenlyLiquidAndAllInterface, ExcessMassModel::EvenlyLiquidAndAllInterfacePreferInterface
   };

}; // class ExcessMassDistributionModel
} // namespace free_surface
} // namespace walberla