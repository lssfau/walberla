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
//! \file PointLight.h
//! \author Lukas Werner
//
//======================================================================================================================

#pragma once

#include <pe/basic.h>
#include <pe/Types.h>
#include <core/math/Vector3.h>

namespace walberla {
namespace pe {
namespace raytracing {
struct Lighting {
   Vec3 pointLightOrigin;
   
   Vec3 ambientLight;
   
   Vec3 diffuseColor;
   real_t diffusePower;
   
   Vec3 specularColor;
   real_t specularPower;
   
   /*!\brief Instantiation constructor for the Lighting struct.
    */
   Lighting () {
      
   }
   
   /*!\brief Instantiation constructor for the Lighting struct.
    * \param pointLightOrigin Origin of the point light.
    * \param ambientLight Color of the ambient light.
    * \param diffuseColor Diffuse color.
    * \param diffusePower Diffuse color power.
    * \param specularColor Specular color.
    * \param specularPower Specular color power.
    */
   Lighting (const Vec3& _pointLightOrigin, const Vec3& _ambientLight,
             const Vec3& _diffuseColor, real_t _diffusePower,
             const Vec3& _specularColor, real_t _specularPower)
   : pointLightOrigin(_pointLightOrigin), ambientLight(_ambientLight),
   diffuseColor(_diffuseColor), diffusePower(_diffusePower),
   specularColor(_specularColor), specularPower(_specularPower) {
      
   }
   
   /*!\brief Instantiation constructor for the Lighting struct.
    * \param config Config handle.
    *
    * The config block has to contain a pointLightOrigin parameter (Vec3).
    * Optional are ambientLight (Vec3), for diffuse coloring diffuseColor (Vec3) and diffusePower (real) and
    * for specular color specularColor (Vec3) and specularPower (real).
    * Colors are Vec3's with values from 0 to 1.
    */
   Lighting (const Config::BlockHandle& config) {
      WALBERLA_CHECK(config.isValid(), "No valid config passed to raytracer lighting.");

      pointLightOrigin = config.getParameter<Vec3>("pointLightOrigin"),
      ambientLight = config.getParameter<Vec3>("ambientLight", Vec3(0,0,0)),
      diffuseColor = config.getParameter<Vec3>("diffuseColor", Vec3(0,0,0));
      diffusePower = config.getParameter<real_t>("diffusePower", real_t(0));
      specularColor = config.getParameter<Vec3>("specularColor", Vec3(0,0,0));
      specularPower = config.getParameter<real_t>("specularPower", real_t(0));
   }
};
}
}
}
