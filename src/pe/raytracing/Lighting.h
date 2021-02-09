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
//! \file Lighting.h
//! \author Lukas Werner
//
//======================================================================================================================

#pragma once

#include <pe/basic.h>
#include <pe/Types.h>
#include <core/math/Vector3.h>
#include <pe/raytracing/Color.h>

namespace walberla {
namespace pe {
namespace raytracing {

/*!\brief The Lighting struct defines the properties of a point light in the scene.
 */
struct Lighting {
   Vec3 pointLightOrigin;
   Color diffuseColor;
   Color specularColor;
   Color ambientColor;
   
   /*!\brief Instantiation constructor for the Lighting struct.
    */
   Lighting () = default;
   
   /*!\brief Instantiation constructor for the Lighting struct.
    * \param pointLightOrigin Origin of the point light.
    * \param diffuseColor Diffuse color (base color of the light).
    * \param specularColor Specular color (color of light refractions on an objects surface).
    * \param ambientColor Color of the ambient light in the scene.
    */
   Lighting (const Vec3& _pointLightOrigin,
             const Color& _diffuseColor, const Color& _specularColor, const Color& _ambientColor)
   : pointLightOrigin(_pointLightOrigin),
   diffuseColor(_diffuseColor), specularColor(_specularColor), ambientColor(_ambientColor) {
      
   }
   
   /*!\brief Instantiation constructor for the Lighting struct.
    * \param config Config handle.
    *
    * The config block has to contain a pointLightOrigin parameter (Vec3).
    * Optional are ambientColor (Vec3), diffuseColor (Vec3), specularColor (Vec3).
    * Colors are Vec3's with values from 0 to 1.
    */
   Lighting (const Config::BlockHandle& config) {
      WALBERLA_CHECK(config.isValid(), "No valid config passed to raytracer lighting.");

      pointLightOrigin = config.getParameter<Vec3>("pointLightOrigin");
      diffuseColor = config.getParameter<Color>("diffuseColor", Color(1,1,1));
      specularColor = config.getParameter<Color>("specularColor", Color(1,1,1));
      ambientColor = config.getParameter<Color>("ambientColor", Color(0.5,0.5,0.5));
   }
};

} //namespace raytracing
} //namespace pe
} //namespace walberla
