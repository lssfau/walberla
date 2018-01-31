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
//! \file FPClassify.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include <cmath>


namespace walberla {
namespace math {


/*******************************************************************************************************************//**
 * \brief   Query if 'x' is NaN.
 *
 * \param   x  Value to be checked.
 *
 * \return  true if x is NaN, else false.
 **********************************************************************************************************************/
template<typename T>
inline bool isnan(T x)
{
   return (std::isnan)(x);
}


 /******************************************************************************************************************//**
 * \brief   Query if 'x' is infinite.
 *
 * \param   x  Value to be checked.
 *
 * \return  true if x is infinite, else false.
 **********************************************************************************************************************/
template<typename T>
inline bool isinf(T x)
{
   return (std::isinf)(x);
}


 /******************************************************************************************************************//**
 * \brief   Query if 'x' is finite.
 *
 * \param   x  Value to be checked.
 *
 * \return  true if x is not NaN and not infinite, else false.
 **********************************************************************************************************************/
template<typename T>
bool finite(T x)
{
    return (std::isfinite)(x);
}



} // namespace math
} // namespace walberla
