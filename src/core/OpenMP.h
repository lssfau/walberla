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
//! \file OpenMP.h
//! \brief Guarded OpenMP include
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"
#include "waLBerlaDefinitions.h"

// MPI SECTION //

#ifdef WALBERLA_BUILD_WITH_OPENMP

#define     WALBERLA_OPENMP_SECTION() if(true)
#define WALBERLA_NON_OPENMP_SECTION() if(false)

#else

#define     WALBERLA_OPENMP_SECTION() if(false)
#define WALBERLA_NON_OPENMP_SECTION() if(true)

#endif

#ifdef _OPENMP

#include <omp.h>

#else

namespace walberla {

#define WALBERLA_OPENMP_FUNCTION_ERROR WALBERLA_ABORT( "Invalid OpenMP function call! In case of compiling without OpenMP, OpenMP functions are not available and shouldn't be called!" );

/* schedule kind constants */
using omp_sched_t = enum omp_sched_t {
omp_sched_static  = 1,
omp_sched_dynamic = 2,
omp_sched_guided  = 3,
omp_sched_auto    = 4
};

/* set API functions */
inline void    omp_set_num_threads (int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void    omp_set_dynamic     (int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void    omp_set_nested      (int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void    omp_set_max_active_levels (int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void    omp_set_schedule          (omp_sched_t, int) { WALBERLA_OPENMP_FUNCTION_ERROR }

/* query API functions */
inline int     omp_get_num_threads  (void) { return 1; }
inline int     omp_get_dynamic      (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_get_nested       (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_get_max_threads  (void) { return 1; }
inline int     omp_get_thread_num   (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_get_num_procs    (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_in_parallel      (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_in_final         (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_get_active_level        (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_get_level               (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_get_ancestor_thread_num (int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_get_team_size           (int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_get_thread_limit        (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_get_max_active_levels   (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void    omp_get_schedule            (omp_sched_t *, int *) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_get_max_task_priority   (void) { WALBERLA_OPENMP_FUNCTION_ERROR }

/* lock API functions */
using omp_lock_t = struct omp_lock_t {
    void * _lk;
};

inline void    omp_init_lock    (omp_lock_t *) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void    omp_set_lock     (omp_lock_t *) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void    omp_unset_lock   (omp_lock_t *) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void    omp_destroy_lock (omp_lock_t *) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_test_lock    (omp_lock_t *) { WALBERLA_OPENMP_FUNCTION_ERROR }

/* nested lock API functions */
using omp_nest_lock_t = struct omp_nest_lock_t {
    void * _lk;
};

inline void    omp_init_nest_lock    (omp_nest_lock_t *) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void    omp_set_nest_lock     (omp_nest_lock_t *) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void    omp_unset_nest_lock   (omp_nest_lock_t *) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void    omp_destroy_nest_lock (omp_nest_lock_t *) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int     omp_test_nest_lock    (omp_nest_lock_t *) { WALBERLA_OPENMP_FUNCTION_ERROR }

/* lock hint type for dynamic user lock */
using omp_lock_hint_t = enum omp_lock_hint_t {
    omp_lock_hint_none           = 0,
    omp_lock_hint_uncontended    = 1,
    omp_lock_hint_contended      = (1<<1 ),
    omp_lock_hint_nonspeculative = (1<<2 ),
    omp_lock_hint_speculative    = (1<<3 ),
    kmp_lock_hint_hle            = (1<<16),
    kmp_lock_hint_rtm            = (1<<17),
    kmp_lock_hint_adaptive       = (1<<18)
};

/* hinted lock initializers */
inline void omp_init_lock_with_hint(omp_lock_t *, omp_lock_hint_t) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void omp_init_nest_lock_with_hint(omp_nest_lock_t *, omp_lock_hint_t) { WALBERLA_OPENMP_FUNCTION_ERROR }

/* time API functions */
inline double  omp_get_wtime (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline double  omp_get_wtick (void) { WALBERLA_OPENMP_FUNCTION_ERROR }

/* OpenMP 4.0 */
inline int   omp_get_default_device (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void  omp_set_default_device (int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int   omp_is_initial_device (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int   omp_get_num_devices (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int   omp_get_num_teams (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int   omp_get_team_num (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int   omp_get_cancellation (void) { WALBERLA_OPENMP_FUNCTION_ERROR }

#   include <stdlib.h>
/* OpenMP 4.5 */
inline int    omp_get_initial_device (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void*  omp_target_alloc(size_t, int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void   omp_target_free(void *, int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int    omp_target_is_present(void *, int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int    omp_target_memcpy(void *, void *, size_t, size_t, size_t, int, int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int    omp_target_memcpy_rect(void *, void *, size_t, int, const size_t *,
                                        const size_t *, const size_t *, const size_t *, const size_t *, int, int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int    omp_target_associate_ptr(void *, void *, size_t, size_t, int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int    omp_target_disassociate_ptr(void *, int) { WALBERLA_OPENMP_FUNCTION_ERROR }

/* OpenMP 4.0 affinity API */
using omp_proc_bind_t = enum omp_proc_bind_t {
    omp_proc_bind_false = 0,
    omp_proc_bind_true = 1,
    omp_proc_bind_master = 2,
    omp_proc_bind_close = 3,
    omp_proc_bind_spread = 4
};

inline omp_proc_bind_t omp_get_proc_bind (void) { WALBERLA_OPENMP_FUNCTION_ERROR }

/* OpenMP 4.5 affinity API */
inline int  omp_get_num_places (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int  omp_get_place_num_procs (int) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void omp_get_place_proc_ids (int, int *) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int  omp_get_place_num (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline int  omp_get_partition_num_places (void) { WALBERLA_OPENMP_FUNCTION_ERROR }
inline void omp_get_partition_place_nums (int *) { WALBERLA_OPENMP_FUNCTION_ERROR }

} //namespace walberla

#endif


