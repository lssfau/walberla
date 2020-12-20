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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/ParticleStorage.h>
#include "OutputSelector.h"

#include "core/Set.h"
#include "core/uid/SUID.h"

#include "vtk/Base64Writer.h"
#include "vtk/PointDataSource.h"
#include "vtk/UtilityFunctions.h"

#include <memory>
#include <string>
#include <vector>

namespace walberla {
namespace mesa_pd {
namespace vtk {

class ParticleVtkOutput : public walberla::vtk::PointDataSource
{
public:
   using ParticleSelectorFunc = std::function<bool (const data::ParticleStorage::iterator& pIt)>;

   ParticleVtkOutput( const std::shared_ptr<data::ParticleStorage>& ps )
      : ps_(ps) { }

   std::vector< Attributes > getAttributes() const override;
   void configure() override {}
   std::vector< Vector3< real_t > > getPoints() override;
   void push( std::ostream& os , const uint_t data, const uint_t point, const uint_t component ) override;
   void push( walberla::vtk::Base64Writer& b64, const uint_t data, const uint_t point, const uint_t component ) override;

   template <typename T>
   void addOutput(const std::string& name) { addOutput(name, std::make_shared<OutputSelector<T>>(T())); }

   void addOutput(const std::string& name, const std::shared_ptr<IOutputSelector>& selector);

   ///sets a function which decides which particles should be written to file
   void setParticleSelector( const ParticleSelectorFunc& func) {particleSelector_ = func;}

   ///returns the number of particles written during the last write
   size_t getParticlesWritten() const {return particleIndices_.size();}

private:
   const std::shared_ptr<data::ParticleStorage>& ps_;
   std::vector<size_t> particleIndices_; ///< stores indices of particles which got selected by particleSelector_
   std::vector<std::pair<std::string, std::shared_ptr<IOutputSelector>>> selectors_;
   ParticleSelectorFunc particleSelector_ = [](const data::ParticleStorage::iterator& /*pIt*/){return true;};
};

} // namespace vtk
} // namespace pe
} // namespace walberla

