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
//! \file MetisWrapper.h
//! \ingroup blockforest
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/mpi/MPIWrapper.h"

namespace walberla {
namespace core {

int METIS_PartGraphKway( int64_t *nvtxs, int64_t *ncon, int64_t *xadj, int64_t *adjncy, int64_t *vwgt, int64_t *vsize, int64_t *adjwgt,
                         int64_t *nparts, double *tpwgts, double *ubvec, int64_t *options, int64_t *edgecut, int64_t *part );

int METIS_PartGraphRecursive( int64_t *nvtxs, int64_t *ncon, int64_t *xadj, int64_t *adjncy, int64_t *vwgt, int64_t *vsize, int64_t *adjwgt,
                              int64_t *nparts, double *tpwgts, double *ubvec, int64_t *options, int64_t *edgecut, int64_t *part );

int METIS_SetDefaultOptions(int64_t *options );

#ifndef METIS_NOPTIONS
#  define METIS_NOPTIONS 40
#endif

extern const int METIS_OK;
extern const int METIS_ERROR;
extern const int METIS_ERROR_INPUT;
extern const int METIS_ERROR_MEMORY;

extern const int METIS_OPTION_PTYPE;
extern const int METIS_OPTION_OBJTYPE;
extern const int METIS_OPTION_CTYPE;
extern const int METIS_OPTION_IPTYPE;
extern const int METIS_OPTION_RTYPE;
extern const int METIS_OPTION_DBGLVL;
extern const int METIS_OPTION_NITER;
extern const int METIS_OPTION_NCUTS;
extern const int METIS_OPTION_SEED;
extern const int METIS_OPTION_NO2HOP;
extern const int METIS_OPTION_MINCONN;
extern const int METIS_OPTION_CONTIG;
extern const int METIS_OPTION_COMPRESS;
extern const int METIS_OPTION_CCORDER;
extern const int METIS_OPTION_PFACTOR;
extern const int METIS_OPTION_NSEPS;
extern const int METIS_OPTION_UFACTOR;
extern const int METIS_OPTION_NUMBERING;
extern const int METIS_OPTION_HELP;
extern const int METIS_OPTION_TPWGTS;
extern const int METIS_OPTION_NCOMMON;
extern const int METIS_OPTION_NOOUTPUT;
extern const int METIS_OPTION_BALANCE;
extern const int METIS_OPTION_GTYPE;
extern const int METIS_OPTION_UBVEC;

} // namespace core
} // namespace walberla
