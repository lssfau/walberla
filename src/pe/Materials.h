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
//! \file Materials.h
//! \author Klaus Iglberger
//! \brief Header file for materials
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iostream>
#include <string>
#include <vector>
#include <pe/Types.h>
#include <core/DataTypes.h>
#include <core/debug/Debug.h>
#include <core/math/Shims.h>
#include <core/math/MatrixMxN.h>

namespace walberla {
namespace pe {

//=================================================================================================
//
//  CLASS MATERIAL
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Rigid body material.
 *
 * A material specifies the properties of a rigid body: the density of the body, the coefficient
 * of restitution and the coefficients of static and dynamic friction.\n
 * The \b pe module provides several predefined materials that can be directly used:
 *
 * - iron
 * - copper
 * - granite
 * - oak
 * - fir
 *
 * In order to create a new custom material use the createMaterial() function:

   \code
   // Creating a new material using the following material properties:
   // - name/identifier: myMaterial
   // - density: 2.54
   // - coefficient of restitution: 0.8
   // - coefficient of static friction: 0.1
   // - coefficient of dynamic friction: 0.05
   // - Poisson's ratio: 0.2
   // - Young's modulus: 80.0
   // - Contact stiffness: 100
   // - dampingN: 10
   // - dampingT: 11
   MaterialID myMaterial = createMaterial( "myMaterial", 2.54, 0.8, 0.1, 0.05, 0.2, 80, 100, 10, 11 );
   \endcode

 * The following functions can be used to acquire a specific MaterialID or to get a specific
 * property of a material:

   \code
   // Searching a material
   MaterialID myMaterial = Material::find( "myMaterial" );

   // Getting the density, coefficient of restitution, coefficient of static and
   // dynamic friction, Poisson's ratio and Young's modulus of the material
   real_t density = Material::getDensity( myMaterial );
   real_t cor     = Material::getRestitution( myMaterial );
   real_t csf     = Material::getStaticFriction( myMaterial );
   real_t cdf     = Material::getDynamicFriction( myMaterial );
   real_t poisson = Material::getPoissonRatio( myMaterial );
   real_t young   = Material::getYoungModulus( myMaterial ):
   \endcode
 */
class WALBERLA_PUBLIC Material
{
private:
   //**Type definitions****************************************************************************
   using SizeType = Materials::size_type;  //!< Size type of the Material class.
   //**********************************************************************************************

public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline Material( const std::string& name, real_t density, real_t cor,
                             real_t csf, real_t cdf, real_t poisson, real_t young, real_t stiffness,
                             real_t dampingN, real_t dampingT );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Copy assignment operator********************************************************************
   // No explicitly declared copy assignment operator.
   //**********************************************************************************************

   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline const std::string&   getName()            const;
   inline real_t               getDensity()         const;
   inline real_t               getRestitution()     const;
   inline real_t               getStaticFriction()  const;
   inline real_t               getDynamicFriction() const;
   inline real_t               getPoissonRatio()    const;
   inline real_t               getYoungModulus()    const;
   inline real_t               getStiffness()       const;
   inline real_t               getDampingN()        const;
   inline real_t               getDampingT()        const;

   static        MaterialID         find( const std::string& name );
   static        std::vector<MaterialID> findPrefix( const std::string& prefix );
   static inline const std::string&   getName( MaterialID material );
   static inline real_t               getDensity( MaterialID material );
   static inline real_t               getRestitution( MaterialID material );
   static inline real_t               getRestitution( MaterialID material1, MaterialID material2 );
   static inline real_t               getStaticFriction( MaterialID material );
   static inline real_t               getStaticFriction( MaterialID material1, MaterialID material2 );
   static inline real_t               getDynamicFriction( MaterialID material );
   static inline real_t               getDynamicFriction( MaterialID material1, MaterialID material2 );
   static inline real_t               getPoissonRatio( MaterialID material );
   static inline real_t               getYoungModulus( MaterialID material );
   static inline real_t               getYoungModulus( MaterialID material1, MaterialID material2 );
   static inline real_t               getStiffness( MaterialID material );
   static inline real_t               getStiffness( MaterialID material1, MaterialID material2 );
   static inline real_t               getDampingN( MaterialID material );
   static inline real_t               getDampingN( MaterialID material1, MaterialID material2 );
   static inline real_t               getDampingT( MaterialID material );
   static inline real_t               getDampingT( MaterialID material1, MaterialID material2 );
   static std::string                 toString( const MaterialID& v );
   //@}
   //**********************************************************************************************

   //**Set functions*******************************************************************************
   /*!\name Set functions */
   //@{
   static inline void setRestitution( MaterialID material1, MaterialID material2, real_t cor );
   static inline void setStaticFriction( MaterialID material1, MaterialID material2, real_t csf );
   static inline void setDynamicFriction( MaterialID material1, MaterialID material2, real_t cdf );
   //@}
   //**********************************************************************************************

private:
   //**Setup functions*****************************************************************************
   /*!\name Setup functions */
   //@{
   static bool activateMaterials();
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string name_;  //!< The name of the material.
   real_t density_;      //!< The density of the material.
   real_t restitution_;  //!< The coefficient of restitution (COR) of a self-similar collision \f$ [0..1] \f$.
                       /*!< The COR represents the energy dissipated during a collision between
                            self-similar bodies, that is bodies with similar materials. A value of
                            0 corresponds to completely inelastic collision where all energy is
                            dissipated, a value of 1 corresponds to a completely elastic collision
                            where no energy is lost. The COR is assumed to be rate-independent. The
                            COR is often determined experimentally by measuring the pre- and
                            post-impact relative velocities:
                            \f[ C_R = \frac{V_{2,after}-V_{1,after}}{V_{2,before}-V_{1,before}} \f]
                            During a collision, the COR values of the two colliding
                            rigid bodies can be used by the collision response mechanism to
                            determine the restitution factor of the contact point. */
   real_t static_;       //!< The coefficient of static friction (CSF) \f$ [0..\infty) \f$.
                       /*!< The CSF is a dimensionless, non-negative quantity representing the
                            amount of static friction between two touching rigid bodies. Static
                            friction occurs in case the relative tangential velocity between the
                            two bodies is 0. Then the force magnitudes of the normal and friction
                            force are related by an inequality:
                            \f[ |\vec{f_t}| \leq \mu_s |\vec{f_n}| \f]
                            The direction of the friction must oppose acceleration if sliding is
                            imminent and is unresticted otherwise. */
   real_t dynamic_;      //!< The coefficient of dynamic friction (CDF) \f$ [0..\infty) \f$.
                       /*!< The CDF is a dimensionless, non-negative quantity representing the
                            amount of dynamic friction between two touching rigid bodies. Dynamic
                            friction occurs in case the relative tangential velocity between the
                            two bodies is greater than 0. Then the force magnitudes of the normal
                            and friction force are related by an inequality:
                            \f[ |\vec{f_t}| = -\mu_d |\vec{f_n}| \frac{\vec{v_t}}{|\vec{v_t}|} \f] */
   real_t poisson_;      //!< The Poisson's ratio for the material \f$ [-1..0.5] \f$.
                       /*!< When a material is compressed in one direction, it usually tends to
                            expand in the other two directions perpendicular to the direction of
                            compression. This effect is called Poisson effect. In this context, the
                            Poisson's ratio is the ratio of the contraction or transverse strain
                            (perpendicular to the applied load) to the extension or axial strain
                            (in the direction of the applied load). For stable, isotropic, linear
                            elastic materials this ratio cannot be less than -1.0 nor greater than
                            0.5 due to the requirement that Young's modulus has positive values. */
   real_t young_;        //!< The Young's modulus for the material \f$ (0..\infty) \f$.
                       /*!< The Young's modulus is a measure for the stiffness of an isotropic
                            elastic material. It is defined as the ratio of the uniaxial stress
                            over the uniaxial strain in the range of stress in which Hooke's law
                            holds. The SI unit for Young's modulus is \f$ Pa \f$ or \f$ N/m^2 \f$. */
   real_t stiffness_;    //!< The stiffness of the contact region \f$ (0..\infty) \f$.
                       /*!< Rigid body theory assumes that the deformation during contact is
                            localized to the contact region. This local compliance can be modelled
                            simplified as a spring-damper. The spring constant corresponds to this
                            parameter. */
   real_t dampingN_;     //!< The damping at the contact region in normal direction \f$ [0..\infty) \f$.
                       /*!< Rigid body theory assumes that the deformation during contact is
                            localized to the contact region. This local compliance in normal
                            direction can be modelled simplified as a spring-damper. The viscous
                            damping coefficient corresponds to this parameter. */
   real_t dampingT_;     //!< The damping at the contact region in tangential direction \f$ [0..\infty) \f$.
                       /*!< Friction counteracts the tangential relative velocity and thus can be
                            modelled as a viscous damper with a limited damping force. The viscous
                            damping coefficient corresponds to this parameter.*/

   static Materials materials_;      //!< Vector for the registered materials.
   static MatN corTable_;            //!< Table for the coefficients of restitution.
   static MatN csfTable_;            //!< Table for the coefficients of static friction.
   static MatN cdfTable_;            //!< Table for the coefficients of dynamic friction.
   static bool materialsActivated_;  //!< Helper variable for the automatic registration process.
   static unsigned int anonymousMaterials_;   //!< Counter for the amount of anonymous materials.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond internal */
   friend MaterialID createMaterial( const std::string& name, real_t density, real_t cor,
                                     real_t csf, real_t cdf, real_t poisson, real_t young,
                                     real_t stiffness, real_t dampingN, real_t dampingT );
   friend MaterialID createMaterial( real_t density, real_t cor, real_t csf, real_t cdf,
                                     real_t poisson, real_t young,
                                     real_t stiffness, real_t dampingN, real_t dampingT );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constructor of the Material class.
 *
 * \param name The name of the material.
 * \param density The density of the material \f$ (0..\infty) \f$.
 * \param cor The coefficient of restitution (COR) of the material \f$ [0..1] \f$.
 * \param csf The coefficient of static friction (CSF) of the material \f$ [0..\infty) \f$.
 * \param cdf The coefficient of dynamic friction (CDF) of the material \f$ [0..\infty) \f$.
 * \param poisson The Poisson's ratio of the material \f$ [-1..0.5] \f$.
 * \param young The Young's modulus of the material \f$ (0..\infty) \f$.
 * \param stiffness The stiffness in normal direction of the material's contact region.
 * \param dampingN The damping coefficient in normal direction of the material's contact region.
 * \param dampingT The damping coefficient in tangential direction of the material's contact region.
 */
inline Material::Material( const std::string& name, real_t density, real_t cor,
                           real_t csf, real_t cdf, real_t poisson, real_t young,
                           real_t stiffness, real_t dampingN, real_t dampingT )
   : name_       ( name )       // The name of the material
   , density_    ( density )    // The density of the material
   , restitution_( cor )        // The coefficient of restitution of the material
   , static_     ( csf )        // The coefficient of static friction of the material
   , dynamic_    ( cdf )        // The coefficient of dynamic friction of the material
   , poisson_    ( poisson )    // The Poisson's ratio for the material
   , young_      ( young )      // The Young's modulus for the material
   , stiffness_  ( stiffness )  // The stiffness in normal direction of the material's contact region.
   , dampingN_   ( dampingN )   // The damping coefficient in normal direction of the material's contact region.
   , dampingT_   ( dampingT )   // The damping coefficient in tangential direction of the material's contact region.
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the name of the material.
 *
 * \return The name of the material.
 */
inline const std::string& Material::getName() const
{
   return name_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the density of the material.
 *
 * \return The density of the material.
 */
inline real_t Material::getDensity() const
{
   return density_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the coefficient of restitution of the material.
 *
 * \return The coefficient of restitution of the material.
 */
inline real_t Material::getRestitution() const
{
   return restitution_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the coefficient of static friction of the material.
 *
 * \return The coefficient of static friction of the material.
 */
inline real_t Material::getStaticFriction() const
{
   return static_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the coefficient of dynamic friction of the material.
 *
 * \return The coefficient of dynamic friction of the material.
 */
inline real_t Material::getDynamicFriction() const
{
   return dynamic_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the Poisson's ratio of the material.
 *
 * \return The Poisson's ratio of the material.
 */
inline real_t Material::getPoissonRatio() const
{
   return poisson_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the Young's modulus of the material.
 *
 * \return The Young's modulus of the material.
 */
inline real_t Material::getYoungModulus() const
{
   return young_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the stiffness in normal direction of the material's contact region.
 *
 * \return The stiffness in normal direction of the material's contact region.
 */
inline real_t Material::getStiffness() const
{
   return stiffness_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the damping coefficient in normal direction of the material's contact region.
 *
 * \return The damping coefficient in normal direction of the material's contact region.
 */
inline real_t Material::getDampingN() const
{
   return dampingN_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the damping coefficient in tangential direction of the material's contact region.
 *
 * \return The damping coefficient in tangential direction of the material's contact region.
 */
inline real_t Material::getDampingT() const
{
   return dampingT_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the name of the given material.
 *
 * \param material The material to be queried.
 * \return The name of the given material.
 */
inline const std::string& Material::getName( MaterialID material )
{
   WALBERLA_ASSERT( material < materials_.size(), "Invalid material ID" );
   return materials_[material].getName();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the density of the given material.
 *
 * \param material The material to be queried.
 * \return The density of the given material.
 */
inline real_t Material::getDensity( MaterialID material )
{
   WALBERLA_ASSERT( material < materials_.size(), "Invalid material ID" );
   return materials_[material].getDensity();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the coefficient of restitution of the given material.
 *
 * \param material The material to be queried.
 * \return The coefficient of restitution of the given material.
 */
inline real_t Material::getRestitution( MaterialID material )
{
   WALBERLA_ASSERT( material < materials_.size(), "Invalid material ID" );
   return materials_[material].getRestitution();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the composite coefficient of restitution for a collision between two rigid bodies.
 *
 * \param material1 The material of the first colliding rigid body.
 * \param material2 The material of the second colliding rigid body.
 * \return The resulting composite coefficient of restitution of the collision.
 */
inline real_t Material::getRestitution( MaterialID material1, MaterialID material2 )
{
   WALBERLA_ASSERT( material1 < materials_.size(), "Invalid material ID" );
   WALBERLA_ASSERT( material2 < materials_.size(), "Invalid material ID" );
   return corTable_( material1, material2 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the coefficient of static friction of the given material.
 *
 * \param material The material to be queried.
 * \return The coefficient of static friction of the given material.
 */
inline real_t Material::getStaticFriction( MaterialID material )
{
   WALBERLA_ASSERT( material < materials_.size(), "Invalid material ID" );
   return materials_[material].getStaticFriction();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the coefficient of static friction for a collision between two rigid bodies.
 *
 * \param material1 The material of the first colliding rigid body.
 * \param material2 The material of the second colliding rigid body.
 * \return The resulting coefficient of static friction of the collision.
 */
inline real_t Material::getStaticFriction( MaterialID material1, MaterialID material2 )
{
   WALBERLA_ASSERT( material1 < materials_.size(), "Invalid material ID" );
   WALBERLA_ASSERT( material2 < materials_.size(), "Invalid material ID" );
   return csfTable_( material1, material2 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the coefficient of dynamic friction of the given material.
 *
 * \param material The material to be queried.
 * \return The coefficient of dynamic friction of the given material.
 */
inline real_t Material::getDynamicFriction( MaterialID material )
{
   WALBERLA_ASSERT( material < materials_.size(), "Invalid material ID" );
   return materials_[material].getDynamicFriction();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the coefficient of dynamic friction for a collision between two rigid bodies.
 *
 * \param material1 The material of the first colliding rigid body.
 * \param material2 The material of the second colliding rigid body.
 * \return The resulting coefficient of dynamic friction of the collision.
 */
inline real_t Material::getDynamicFriction( MaterialID material1, MaterialID material2 )
{
   WALBERLA_ASSERT( material1 < materials_.size(), "Invalid material ID" );
   WALBERLA_ASSERT( material2 < materials_.size(), "Invalid material ID" );
   return cdfTable_( material1, material2 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the Poisson's ratio of the given material.
 *
 * \param material The material to be queried.
 * \return The Poisson's ratio of the given material.
 */
inline real_t Material::getPoissonRatio( MaterialID material )
{
   WALBERLA_ASSERT( material < materials_.size(), "Invalid material ID" );
   return materials_[material].getPoissonRatio();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the Young's modulus of the given material.
 *
 * \param material The material to be queried.
 * \return The Young's modulus of the given material.
 */
inline real_t Material::getYoungModulus( MaterialID material )
{
   WALBERLA_ASSERT( material < materials_.size(), "Invalid material ID" );
   return materials_[material].getYoungModulus();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the (effective) Young's modulus for a collision between two rigid bodies.
 *
 * \param material1 The material of the first colliding rigid body.
 * \param material2 The material of the second colliding rigid body.
 * \return The resulting (effective) Young's modulus of the collision.
 *
 * This function returns the effective Young's modulus for a collision between two rigid bodies.
 * The effective Young's modulus is calculated as

          \f[ \frac{1}{E_{eff}} = \frac{1 - \nu_1^2}{E_1} + \frac{1 - \nu_2^2}{E_2}, \f]

 * where \f$ E_1 \f$ and \f$ E_2 \f$ are the Young's modulus for the first and second material,
 * respectively, and \f$ \nu_1 \f$ and \f$ \nu_2 \f$ are the Poisson's ratio for the materials.
 */
inline real_t Material::getYoungModulus( MaterialID material1, MaterialID material2 )
{
   WALBERLA_ASSERT( material1 < materials_.size(), "Invalid material ID" );
   WALBERLA_ASSERT( material2 < materials_.size(), "Invalid material ID" );

   const real_t nu1( getPoissonRatio( material1 ) );
   const real_t nu2( getPoissonRatio( material2 ) );
   const real_t y1 ( getYoungModulus( material1 ) );
   const real_t y2 ( getYoungModulus( material2 ) );

   const real_t tmp1( y2 * ( real_c(1) - nu1*nu1 ) );
   const real_t tmp2( y1 * ( real_c(1) - nu2*nu2 ) );

   return ( ( y1*y2 ) / ( tmp1 + tmp2 ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the stiffness in normal direction of the material's contact region.
 *
 * \param material The material to be queried.
 * \return The stiffness in normal direction of the contact region of the given material.
 */
inline real_t Material::getStiffness( MaterialID material )
{
   WALBERLA_ASSERT( material < materials_.size(), "Invalid material ID" );
   return materials_[material].getStiffness();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the stiffness in normal direction of the contact between two materials.
 *
 * \param material1 The material of the first colliding rigid body.
 * \param material2 The material of the second colliding rigid body.
 * \return The stiffness in normal direction of the contact between two materials.
 *
 * Rigid body theory assumes that deformation during contact is localized to the contact region.
 * Therefore the contact region is often modelled simplified as a spring-damper. When two bodies
 * are in contact the spring-dampers are serially connected and thus the contact stiffness can
 * be expressed as the series connection of two springs: \f$ k_*^{-1} = k_1^{-1} + k_2^{-1}\f$.
 */
inline real_t Material::getStiffness( MaterialID material1, MaterialID material2 )
{
   WALBERLA_ASSERT( material1 < materials_.size(), "Invalid material ID" );
   WALBERLA_ASSERT( material2 < materials_.size(), "Invalid material ID" );

   return math::inv( math::inv( getStiffness( material1 ) ) + math::inv( getStiffness( material2 ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the damping coefficient in normal direction of the material's contact region.
 *
 * \param material The material to be queried.
 * \return The damping in normal direction of the contact region of the given material.
 */
inline real_t Material::getDampingN( MaterialID material )
{
   WALBERLA_ASSERT( material < materials_.size(), "Invalid material ID" );
   return materials_[material].getDampingN();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the damping in normal direction of the contact between two materials.
 *
 * \param material1 The material of the first colliding rigid body.
 * \param material2 The material of the second colliding rigid body.
 * \return The damping in normal direction of the contact between two materials.
 *
 * Rigid body theory assumes that deformation during contact is localized to the contact region.
 * Therefore the contact region is often modelled simplified as a spring-damper. When two bodies
 * are in contact the spring-dampers are serially connected and thus the contact damping can
 * be expressed as the series connection of two viscous dampers: \f$ c_*^{-1} = c_1^{-1} + c_2^{-1}\f$.
 */
inline real_t Material::getDampingN( MaterialID material1, MaterialID material2 )
{
   WALBERLA_ASSERT( material1 < materials_.size(), "Invalid material ID" );
   WALBERLA_ASSERT( material2 < materials_.size(), "Invalid material ID" );

   return math::inv( math::inv( getDampingN( material1 ) ) + math::inv( getDampingN( material2 ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the damping coefficient in tangential direction of the material's contact region.
 *
 * \param material The material to be queried.
 * \return The damping in tangential direction of the contact region of the given material.
 */
inline real_t Material::getDampingT( MaterialID material )
{
   WALBERLA_ASSERT( material < materials_.size(), "Invalid material ID" );
   return materials_[material].getDampingT();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the damping in tangential direction of the contact between two materials.
 *
 * \param material1 The material of the first colliding rigid body.
 * \param material2 The material of the second colliding rigid body.
 * \return The damping in tangential direction of the contact between two materials.
 *
 * Rigid body theory assumes that deformation during contact is localized to the contact region.
 * Therefore the contact region is often modelled simplified as a spring-damper. When two bodies
 * are in contact the spring-dampers are serially connected and thus the contact damping can
 * be expressed as the series connection of two viscous dampers: \f$ c_*^{-1} = c_1^{-1} + c_2^{-1}\f$.
 */
inline real_t Material::getDampingT( MaterialID material1, MaterialID material2 )
{
   WALBERLA_ASSERT( material1 < materials_.size(), "Invalid material ID" );
   WALBERLA_ASSERT( material2 < materials_.size(), "Invalid material ID" );

   return math::inv( math::inv( getDampingT( material1 ) ) + math::inv( getDampingT( material2 ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the coefficient of restitution between material \a material1 and \a material2.
 *
 * \param material1 The material of the first colliding rigid body.
 * \param material2 The material of the second colliding rigid body.
 * \param cor The coefficient of restitution between \a material1 and \a material2.
 * \return void
 */
inline void Material::setRestitution( MaterialID material1, MaterialID material2, real_t cor )
{
   WALBERLA_ASSERT( material1 < materials_.size(), "Invalid material ID" );
   WALBERLA_ASSERT( material2 < materials_.size(), "Invalid material ID" );
   corTable_( material1, material2 ) = cor;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the coefficient of static friction between material \a material1 and \a material2.
 *
 * \param material1 The material of the first colliding rigid body.
 * \param material2 The material of the second colliding rigid body.
 * \param csf The coefficient of static friction between \a material1 and \a material2.
 * \return void
 */
inline void Material::setStaticFriction( MaterialID material1, MaterialID material2, real_t csf )
{
   WALBERLA_ASSERT( material1 < materials_.size(), "Invalid material ID" );
   WALBERLA_ASSERT( material2 < materials_.size(), "Invalid material ID" );
   csfTable_( material1, material2 ) = csf;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the coefficient of dynamic friction between material \a material1 and \a material2.
 *
 * \param material1 The material of the first colliding rigid body.
 * \param material2 The material of the second colliding rigid body.
 * \param cdf The coefficient of dynamic friction between \a material1 and \a material2.
 * \return void
 */
inline void Material::setDynamicFriction( MaterialID material1, MaterialID material2, real_t cdf )
{
   WALBERLA_ASSERT( material1 < materials_.size(), "Invalid material ID" );
   WALBERLA_ASSERT( material2 < materials_.size(), "Invalid material ID" );
   cdfTable_( material1, material2 ) = cdf;
}
//*************************************************************************************************




//=================================================================================================
//
//  MATERIAL IRON
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specification of the material iron.
 *
 * The Iron class represents the material iron. It is implemented as a veneer class for the
 * Material base class to set the properties of iron:
 *
 * - Name: "iron"
 * - Density: \f$ 7874 \frac{kg}{m^3} \f$
 * - Coefficient of restitution: 0.5
 * - Coefficient of static friction: 0.1
 * - Coefficient of dynamic friction: 0.1
 * - Poisson's ratio: 0.24
 * - Young's modulus: 200 GPa
 * - Stiffness: \f$ ~200 \frac{N}{m} \f$
 * - Normal Damping: \f$ 0 \frac{Ns}{m} \f$
 * - Tangential Damping: \f$ 0 \frac{Ns}{m} \f$
 *
 * Since several parameters are not unitless they might not match the scaling of the simulation.
 * In that case custom materials must be created. Also even though the stiffness is proportional
 * to Young's modulus the proportionality constant depends on other parameters such as the shape of
 * the contact region or the radii of the objects. Thus if the simulation does rely on the value of
 * the stiffness the user must supply an appropriate stiffness coefficient. Since no published
 * values were available for the damping coefficients they are deactivated.
 *
 * The iron material is automatically registered and can be directly used:
   \code
   MaterialID iron = Material::find("iron");
   \endcode
 */
class WALBERLA_PUBLIC Iron : public Material
{
public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline Iron();
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default constructor for the Iron class.
 */
inline Iron::Iron()
   : Material( "iron", real_c( 7.874 ), real_c( 0.5 ), real_c( 0.1 ), real_c( 0.1 ), real_c( 0.24 ), real_c( 200 ), real_c( 200 ), real_c( 0 ), real_c( 0 ) )  // Initialization of the Material base class
{

}
//*************************************************************************************************




//=================================================================================================
//
//  MATERIAL COPPER
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specification of the material copper.
 *
 * The Copper class represents the material copper. It is implemented as a veneer class for
 * the Material base class to set the properties of iron:
 *
 * - Name: "copper"
 * - Density: \f$ 8920 \frac{kg}{m^3} \f$
 * - Coefficient of restitution: 0.5
 * - Coefficient of static friction: 0.1
 * - Coefficient of dynamic friction: 0.1
 * - Poisson's ratio: 0.33
 * - Young's modulus: 117 GPa
 * - Stiffness: \f$ ~117 \frac{N}{m} \f$
 * - Normal Damping: \f$ 0 \frac{Ns}{m} \f$
 * - Tangential Damping: \f$ 0 \frac{Ns}{m} \f$
 *
 * Since several parameters are not unitless they might not match the scaling of the simulation.
 * In that case custom materials must be created. Also even though the stiffness is proportional
 * to Young's modulus the proportionality constant depends on other parameters such as the shape of
 * the contact region or the radii of the objects. Thus if the simulation does rely on the value of
 * the stiffness the user must supply an appropriate stiffness coefficient. Since no published
 * values were available for the damping coefficients they are deactivated.
 *
 * The copper material is automatically registered and can be directly used:
   \code
   MaterialID iron = Material::find("copper");
   \endcode
 */
class WALBERLA_PUBLIC Copper : public Material
{
public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline Copper();
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default constructor for the Copper class.
 */
inline Copper::Copper()
   : Material( "copper", real_c( 8920 ), real_c( 0.5 ), real_c( 0.1 ), real_c( 0.1 ), real_c( 0.33 ), real_c( 117 ), real_c( 117 ), real_c( 0 ), real_c( 0 ) )  // Initialization of the Material base class
{

}
//*************************************************************************************************




//=================================================================================================
//
//  MATERIAL GRANITE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specification of the material granite.
 *
 * The Granite class represents the material granite. It is implemented as a veneer class for
 * the Material base class to set the properties of granite:
 *
 * - Name: "granite"
 * - Density: \f$ 2800 \frac{kg}{dm^3} \f$
 * - Coefficient of restitution: 0.5
 * - Coefficient of static friction: 0.1
 * - Coefficient of dynamic friction: 0.1
 * - Poisson's ratio: 0.25
 * - Young's modulus: 55 GPa
 * - Stiffness: \f$ ~55 \frac{N}{m} \f$
 * - Normal Damping: \f$ 0 \frac{Ns}{m} \f$
 * - Tangential Damping: \f$ 0 \frac{Ns}{m} \f$
 *
 * Since several parameters are not unitless they might not match the scaling of the simulation.
 * In that case custom materials must be created. Also even though the stiffness is proportional
 * to Young's modulus the proportionality constant depends on other parameters such as the shape of
 * the contact region or the radii of the objects. Thus if the simulation does rely on the value of
 * the stiffness the user must supply an appropriate stiffness coefficient. Since no published
 * values were available for the damping coefficients they are deactivated.
 *
 * The granite material is automatically registered and can be directly used:
   \code
   MaterialID iron = Material::find("granite");
   \endcode
 */
class WALBERLA_PUBLIC Granite : public Material
{
public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline Granite();
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default constructor for the Granite class.
 */
inline Granite::Granite()
   : Material( "granite", real_c( 2800 ), real_c( 0.5 ), real_c( 0.1 ), real_c( 0.1 ), real_c( 0.25 ), real_c( 55 ), real_c( 55 ), real_c( 0 ), real_c( 0 ) )  // Initialization of the Material base class
{

}
//*************************************************************************************************




//=================================================================================================
//
//  MATERIAL OAK
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specification of the material oak.
 *
 * The Oak class represents the material oak wood. It is implemented as a veneer class for the
 * Material base class to set the properties of oak wood:
 *
 * - Name: "oak"
 * - Density: \f$ 800 \frac{kg}{m^3} \f$
 * - Coefficient of restitution: 0.5
 * - Coefficient of static friction: 0.1
 * - Coefficient of dynamic friction: 0.1
 * - Poisson's ratio: 0.35
 * - Young's modulus: 11 GPa
 * - Stiffness: \f$ ~11 \frac{N}{m} \f$
 * - Normal Damping: \f$ 0 \frac{Ns}{m} \f$
 * - Tangential Damping: \f$ 0 \frac{Ns}{m} \f$
 *
 * Since several parameters are not unitless they might not match the scaling of the simulation.
 * In that case custom materials must be created. Also even though the stiffness is proportional
 * to Young's modulus the proportionality constant depends on other parameters such as the shape of
 * the contact region or the radii of the objects. Thus if the simulation does rely on the value of
 * the stiffness the user must supply an appropriate stiffness coefficient. Since no published
 * values were available for the damping coefficients they are deactivated.
 *
 * The oak material is automatically registered and can be directly used:
   \code
   MaterialID iron = Material::find("oak");
   \endcode
 */
class WALBERLA_PUBLIC Oak : public Material
{
public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline Oak();
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default constructor for the Oak class.
 */
inline Oak::Oak()
   : Material( "oak", real_c( 800 ), real_c( 0.5 ), real_c( 0.1 ), real_c( 0.1 ), real_c( 0.35 ), real_c( 11 ), real_c( 11 ), real_c( 0 ), real_c( 0 ) )  // Initialization of the Material base class
{

}
//*************************************************************************************************




//=================================================================================================
//
//  MATERIAL FIR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specification of the material fir.
 *
 * The Fir class represents the material fir wood. It is implemented as a veneer class for the
 * Material base class to set the properties of fir wood:
 *
 * - Name: "fir"
 * - Density: \f$ 500 \frac{kg}{m^3} \f$
 * - Coefficient of restitution: 0.5
 * - Coefficient of static friction: 0.1
 * - Coefficient of dynamic friction: 0.1
 * - Poisson's ratio: 0.34
 * - Young's modulus: 13 GPa
 * - Stiffness: \f$ ~13 \frac{N}{m} \f$
 * - Normal Damping: \f$ 0 \frac{Ns}{m} \f$
 * - Tangential Damping: \f$ 0 \frac{Ns}{m} \f$
 *
 * Since several parameters are not unitless they might not match the scaling of the simulation.
 * In that case custom materials must be created. Also even though the stiffness is proportional
 * to Young's modulus the proportionality constant depends on other parameters such as the shape of
 * the contact region or the radii of the objects. Thus if the simulation does rely on the value of
 * the stiffness the user must supply an appropriate stiffness coefficient. Since no published
 * values were available for the damping coefficients they are deactivated.
 *
 * The fir material is automatically registered and can be directly used:
   \code
   MaterialID iron = Material::find("fir");
   \endcode
 */
class WALBERLA_PUBLIC Fir : public Material
{
public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline Fir();
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default constructor for the Fir class.
 */
inline Fir::Fir()
   : Material( "fir", real_c( 500 ), real_c( 0.5 ), real_c( 0.1 ), real_c( 0.1 ), real_c( 0.34 ), real_c( 13 ), real_c( 13 ), real_c( 0 ), real_c( 0 ) )  // Initialization of the Material base class
{

}
//*************************************************************************************************




//=================================================================================================
//
//  MATERIAL CONSTANTS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief ID for an invalid material.
 *
 * This MaterialID is returned by the getMaterial() function in case no material with the
 * specified name is returned. This value should not be used to create rigid bodies or in
 * any other function!
 */
const MaterialID invalid_material = static_cast<MaterialID>( -1 );
//*************************************************************************************************




//=================================================================================================
//
//  MATERIAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Material functions */
//@{
WALBERLA_PUBLIC MaterialID createMaterial( const std::string& name, real_t density, real_t cor,
                           real_t csf, real_t cdf, real_t poisson, real_t young, real_t stiffness, real_t dampingN, real_t dampingT );
WALBERLA_PUBLIC MaterialID createMaterial( real_t density, real_t cor, real_t csf, real_t cdf,
                           real_t poisson, real_t young, real_t stiffness, real_t dampingN, real_t dampingT );
//@}
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla

