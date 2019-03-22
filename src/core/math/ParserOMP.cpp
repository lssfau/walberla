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
//! \file ParserOMP.cpp
//! \ingroup core
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#include "ParserOMP.h"
#include "core/debug/Debug.h"
#include "core/DataTypes.h"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num()  0
#endif


namespace walberla {
namespace math {



FunctionParserOMP::FunctionParserOMP()
: parser_(new FunctionParser[uint_c(omp_get_max_threads())])
#ifndef NDEBUG
   , num_parsers_(omp_get_max_threads())
#endif
{
}

void FunctionParserOMP::parse(const std::string & eq)
{
   WALBERLA_ASSERT_EQUAL(omp_get_max_threads(), num_parsers_);

   #ifdef _OPENMP
   #pragma omp parallel for schedule(static)
   #endif
   for (int t = 0; t < omp_get_max_threads(); ++t)
   {
      parser_[uint_c(t)].parse(eq);
   }
}

double FunctionParserOMP::evaluate(const std::map<std::string,double> & symbolTable) const
{
   WALBERLA_ASSERT_EQUAL(omp_get_max_threads(), num_parsers_);

   return parser_[uint_c(omp_get_thread_num())].evaluate(symbolTable);
}



} // namespace math
} // namespace walberla
