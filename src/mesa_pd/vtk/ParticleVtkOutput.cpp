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
//! \file ParticleVtkOutput.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "ParticleVtkOutput.h"

namespace walberla {
namespace mesa_pd {
namespace vtk {

std::vector< ParticleVtkOutput::Attributes > ParticleVtkOutput::getAttributes() const
{
   std::vector< Attributes > attributes;
   for (const auto& s : selectors_)
   {
      attributes.emplace_back( s.second->type_string, s.first, s.second->components );
   }
   return attributes;
}

std::vector< math::Vector3< real_t > > ParticleVtkOutput::getPoints()
{
   std::vector< math::Vector3< real_t > > result;
   result.reserve(ps_->size());

   particleIndices_.clear();
   particleIndices_.reserve(ps_->size());

   for (auto pIt = ps_->begin(); pIt != ps_->end(); ++pIt)
   {
      if (particleSelector_(pIt))
      {
         particleIndices_.emplace_back(pIt.getIdx());
         result.emplace_back(pIt->getPosition());
      }
   }

   return result;
}

void ParticleVtkOutput::push( std::ostream& os, const uint_t data, const uint_t point, const uint_t component )
{
   WALBERLA_ASSERT_LESS( data, selectors_.size() );
   WALBERLA_ASSERT_LESS( point, particleIndices_.size() );
   WALBERLA_ASSERT_LESS( particleIndices_[point], ps_->size() );

   selectors_[data].second->push(os, *(*ps_)[particleIndices_[point]], component);
}

void ParticleVtkOutput::push( walberla::vtk::Base64Writer& b64, const uint_t data, const uint_t point, const uint_t component )
{
   WALBERLA_ASSERT_LESS( data, selectors_.size() );
   WALBERLA_ASSERT_LESS( point, particleIndices_.size() );
   WALBERLA_ASSERT_LESS( particleIndices_[point], ps_->size() );

   selectors_[data].second->push(b64, *(*ps_)[particleIndices_[point]], component);
}

void ParticleVtkOutput::addOutput(const std::string& name, std::shared_ptr<IOutputSelector> selector)
{
   if ( std::find_if(selectors_.begin(), selectors_.end(), [&name](const auto& item){return item.first==name;} ) !=
        selectors_.end() )
   {
      WALBERLA_LOG_WARNING("Output " << name << " already registered!");
      return;
   }
   selectors_.emplace_back(name, selector);
}

} // namespace vtk
} // namespace pe
} // namespace walberla

