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
//! \file Iterator.h
//! \ingroup config
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

namespace walberla {
namespace config {


   struct ConfigGenerator
   {
      virtual ~ConfigGenerator() = default;
      virtual shared_ptr<Config> next() = 0;
   };


   class Iterator
   {
   public:
      Iterator() = default;
      Iterator( const shared_ptr<ConfigGenerator> & configGenerator )
          : generator_ ( configGenerator )
       {
          currentConfig_ = generator_->next();
       }

      Iterator & operator++ () {
          currentConfig_ = generator_->next();
          return *this;
       }

       inline bool operator==( const Iterator& it ) const { return it.currentConfig_ == currentConfig_; }
       inline bool operator!=( const Iterator& it ) const { return it.currentConfig_ != currentConfig_; }

       inline shared_ptr<Config> & operator*()  { return  currentConfig_; }
       inline shared_ptr<Config> * operator->() { return &currentConfig_; }

    private:
       shared_ptr<ConfigGenerator> generator_;
       shared_ptr<Config>          currentConfig_;
   };



} // namespace config
} // namespace walberla


