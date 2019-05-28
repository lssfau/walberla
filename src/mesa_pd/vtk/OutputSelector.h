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
//! \file OutputSelector.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ParticleStorage.h>
#include "WriteOutput.h"

#include <vtk/Base64Writer.h>
#include <vtk/VTKTrait.h>

#include <ostream>
#include <string>

namespace walberla {
namespace mesa_pd {
namespace vtk {

class IOutputSelector
{
public:
   IOutputSelector( char const * const ts, const uint_t c) : type_string(ts), components(c) {}
   virtual ~IOutputSelector() {}
   virtual void push( std::ostream& os , const data::Particle&& p, const uint_t component ) = 0;
   virtual void push( walberla::vtk::Base64Writer& b64, const data::Particle&& p, const uint_t component ) = 0;

   const std::string type_string;
   const uint_t components;
};

template <typename Selector>
class OutputSelector : public IOutputSelector
{
public:
   OutputSelector(Selector s) :
      IOutputSelector(walberla::vtk::VTKTrait<typename Selector::return_type>::type_string,
                      walberla::vtk::VTKTrait<typename Selector::return_type>::components),
      selector_(s) {}
   inline void push( std::ostream& os , const data::Particle&& p, const uint_t component ) override
   {
      writeOutput( os, selector_(p), component );
   }
   inline void push( walberla::vtk::Base64Writer& b64, const data::Particle&& p, const uint_t component ) override
   {
      writeOutput( b64, selector_(p), component );
   }
private:
   Selector selector_;
};

} // namespace vtk
} // namespace pe
} // namespace walberla

