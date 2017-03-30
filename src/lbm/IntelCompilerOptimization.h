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
//! \file IntelCompilerOptimization.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/perf_analysis/extern/iacaMarks.h"

// Blocked x-loop:
#ifdef __INTEL_COMPILER
#define INTEL_COMPILER_XLOOP_OPT
#endif


#ifdef INTEL_COMPILER_XLOOP_OPT


#define X_LOOP(loopBody) {\
                              const cell_idx_t unroll = 4;\
                              cell_idx_t x = 0;\
                              cell_idx_t outerCounter = 0;\
                              const cell_idx_t xSizeMinusUnroll = xSize-unroll;\
                              for ( outerCounter = 0; outerCounter <= xSizeMinusUnroll; outerCounter+=unroll )\
                              {\
                                 _Pragma("unroll")\
                                 _Pragma("vector always")\
                                 _Pragma("ivdep")\
                                 /*_Pragma("vector aligned")*/\
                                 x = outerCounter;\
                                 for( cell_idx_t innerCounter = 0; innerCounter != unroll; ++innerCounter )\
                                 {\
                                    loopBody\
                                    ++x;\
                                 }\
                              }\
                              for( ; x < xSize; ++x ) {\
                                 loopBody\
                              }\
                           }


#define X_LOOP_IACA(loopBody) {\
                              const cell_idx_t unroll = 4;\
                              cell_idx_t x = 0;\
                              cell_idx_t outerCounter = 0;\
                              const cell_idx_t xSizeMinusUnroll = xSize-unroll;\
                              for ( outerCounter = 0; outerCounter <= xSizeMinusUnroll; outerCounter+=unroll )\
                              {\
                                 _Pragma("unroll")\
                                 _Pragma("vector always")\
                                 _Pragma("ivdep")\
                                 /*_Pragma("vector aligned")*/\
                                 x = outerCounter;\
                                 for( cell_idx_t innerCounter = 0; innerCounter != unroll; ++innerCounter )\
                                 {\
                                    IACA_START\
                                    loopBody\
                                    ++x;\
                                 }\
                                 IACA_END\
                              }\
                              for( ; x < xSize; ++x ) {\
                                 loopBody\
                              }\
                           }


#else


#define X_LOOP(loopBody) for( cell_idx_t x = 0; x != xSize; ++x ) { loopBody }


#define X_LOOP_IACA(loopBody) for( cell_idx_t x = 0; x != xSize; ++x ) { IACA_START loopBody } IACA_END


#endif
