# encoding: utf-8

import math
import operator
from functools import reduce
import waLBerla as wlb

# Constants that define the size of blocks that are used in the benchmarks
MIN_CELLS_PER_BLOCK = 16
MAX_CELLS_PER_BLOCK = 256
INC_CELLS_PER_BLOCK = 16
# Amount of cells per block
cells_per_block_interval = range(MIN_CELLS_PER_BLOCK, MAX_CELLS_PER_BLOCK + 1, INC_CELLS_PER_BLOCK)
# Blocks with size in [16, 32, 64, 128, 256]
cells_per_block = [num_cells for num_cells in cells_per_block_interval]
# Number of active MPI processes
num_processes = wlb.mpi.numProcesses()
# Whether to overlap computation with communication
overlap_communication = [False, True]
# Whether MPI supports buffers in GPU memory
cuda_enabled_mpi = [False, True]
# Supported communication schemes
communication_schemes = ['GPUPackInfo_Streams', 'UniformGPUScheme_Baseline', 'UniformGPUScheme_Memcpy']


def calculate_time_steps(runtime, expected_mlups, domain_size):
    cells = reduce(operator.mul, domain_size, 1)
    time_steps_per_second = expected_mlups * 1e6 / cells
    return int(time_steps_per_second * runtime)


def side_length_to_fill_memory(memory_fill_percentage, memory_in_gb):
    bytes_per_cell = 19 * 2 * 8
    max_cells = memory_in_gb * 1e9 / bytes_per_cell * memory_fill_percentage
    return int(max_cells**(1/3))


def get_block_decomposition(block_decomposition, num_processes):
    bx = by = bz = 1
    blocks_per_axis = int(math.log(num_processes, 2))
    for i in range(blocks_per_axis):
        decomposition_axis = block_decomposition[i % len(block_decomposition)]
        if decomposition_axis == 'y':
            by *= 2
        elif decomposition_axis == 'z':
            bz *= 2
        elif decomposition_axis == 'x':
            bx *= 2

    assert (bx * by * bz) == num_processes
    return bx, by, bz
