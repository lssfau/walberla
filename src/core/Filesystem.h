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
//! \file Filesystem.h
//! \ingroup core
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once


#if defined(WALBERLA_USE_STD_FILESYSTEM)
#include <filesystem>
#elif defined(WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM)
#define _LIBCPP_NO_EXPERIMENTAL_DEPRECATION_WARNING_FILESYSTEM
#include <experimental/filesystem>
#else
#include <boost/filesystem.hpp>
#endif



namespace walberla {
namespace filesystem {

#if defined(WALBERLA_USE_STD_FILESYSTEM)
using namespace std::filesystem;
#elif defined(WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM)
using namespace std::experimental::filesystem;
#else
using namespace boost::filesystem;
#endif

}
}
