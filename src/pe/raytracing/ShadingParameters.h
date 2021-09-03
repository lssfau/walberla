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
//! \file ShadingParameters.h
//! \author Lukas Werner
//
//======================================================================================================================

#pragma once

#include <pe/basic.h>
#include <pe/Types.h>
#include <pe/raytracing/Color.h>

namespace walberla {
namespace pe {
namespace raytracing {

struct ShadingParameters {
   Color diffuseColor;  //!< Primary color of the material.
   Color ambientColor;  //!< Color the material has even when its not directly lit.
   Color specularColor; //!< Color the specular highlight has.
   real_t shininess;    //!< How shiny a material is (approximate range is between 1 and 100).
   
   /*!\brief Instantiation constructor for the Shading struct.
    */
   ShadingParameters () = default;
   
   /*!\brief Instantiation constructor for the Shading struct.
    * \param _diffuseColor Primary color of the material.
    * \param _ambientColor Color the material has even when its not directly lit.
    * \param _specularColor Color this material contributes to on its specular highlights.
    * \param _shininess Shininess of the material.
    */
   ShadingParameters (const Color& _diffuseColor, const Color& _ambientColor, const Color& _specularColor, real_t _shininess)
   : diffuseColor(_diffuseColor), ambientColor(_ambientColor), specularColor(_specularColor), shininess(_shininess) {
      
   }
   
   /*!\brief Instantiation constructor for the Shading struct.
    * \param config Config handle.
    *
    * Config has to contain diffuseColor (Color), ambientColor (Color), specularColor (Color), shininess (real_t).
    * Colors are Vec3's with values from 0 to 1.
    */
   ShadingParameters (const Config::BlockHandle& config) {
      diffuseColor = config.getParameter<Color>("diffuseColor");
      ambientColor = config.getParameter<Color>("ambientColor");
      specularColor = config.getParameter<Color>("specularColor");
      shininess = config.getParameter<real_t>("shininess");
   }
   
   /*!\brief Makes a rendered object shiny by setting the shininess and adjusting the specularColor.
    * \param _shininess Shininess
    */
   ShadingParameters& makeGlossy(real_t _shininess = 30) {
      shininess = _shininess;
      specularColor.set(real_t(1), real_t(1), real_t(1));
      return *this;
   }
   
   /*!\brief Makes the rendered object matte by setting the shininess attribute to zero and adjusting the specularColor.
    */
   ShadingParameters& makeMatte() {
      shininess = 0;
      specularColor.set(real_t(0.1), real_t(0.1), real_t(0.1));
      return *this;
   }
};

} //namespace raytracing
} //namespace pe
} //namespace walberla

