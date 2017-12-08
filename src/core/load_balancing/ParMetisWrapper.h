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
//! \file ParMetisWrapper.h
//! \ingroup blockforest
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "MetisWrapper.h"

#include "core/DataTypes.h"
#include "core/mpi/MPIWrapper.h"

namespace walberla {
namespace core {

int ParMETIS_V3_AdaptiveRepart(
   int64_t *vtxdist, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
   int64_t *vsize, int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *ncon,
   int64_t *nparts, double *tpwgts, double *ubvec, double *ipc2redist,
   int64_t *options, int64_t *edgecut, int64_t *part, MPI_Comm *comm );

int ParMETIS_V3_PartKway(
   int64_t *vtxdist, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
   int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *ncon, int64_t *nparts,
   double *tpwgts, double *ubvec, int64_t *options, int64_t *edgecut, int64_t *part,
   MPI_Comm *comm );

int ParMETIS_V3_PartGeomKway(
   int64_t *vtxdist, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
   int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *ndims, double *xyz,
   int64_t *ncon, int64_t *nparts, double *tpwgts, double *ubvec, int64_t *options,
   int64_t *edgecut, int64_t *part, MPI_Comm *comm );

int ParMETIS_V3_PartGeom(
   int64_t *vtxdist, int64_t *ndims, double *xyz, int64_t *part, MPI_Comm *comm );

int ParMETIS_V3_RefineKway(
   int64_t *vtxdist, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
   int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *ncon, int64_t *nparts,
   double *tpwgts, double *ubvec, int64_t *options, int64_t *edgecut,
   int64_t *part, MPI_Comm *comm );

} // namespace core
} // namespace walberla
