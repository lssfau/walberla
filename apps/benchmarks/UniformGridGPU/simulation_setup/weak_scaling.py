# encoding: utf-8

import itertools
import waLBerla as wlb
from base import get_block_decomposition, communication_schemes, overlap_communication, \
                 cuda_enabled_mpi, num_processes, calculate_time_steps, side_length_to_fill_memory
from benchmark import BenchmarkScenario, CommunicationSchemeType


# Stores the scenarios for the current simulation
scenarios = wlb.ScenarioManager()

# Generates all block decompositions of xyz, 2 directions at a time
#block_decompositions = itertools.combinations_with_replacement('xyz', r=2)
block_decompositions = ['xyz', 'yzx', 'zyx', 'yxz']

# compute number of cells depending on GPU memory i.e. by specifying the percentage of GPU memory to fill
gpu_memory_gb = 16
cells_per_block = [side_length_to_fill_memory(pc, gpu_memory_gb) for pc in (0.8, 0.5, 0.05)]

expected_mlups = 200  # to compute how many time steps have to be done
time_per_scenarios = 5  # benchmark time in seconds

fully_periodic = [False, True]

if num_processes == 1:
    scenario_generator = itertools.product(communication_schemes, [False, ], [False, ],
                                           block_decompositions, cells_per_block, fully_periodic)
else:
    scenario_generator = itertools.product(communication_schemes, [True],
                                           cuda_enabled_mpi, block_decompositions, cells_per_block, fully_periodic)

testcase_name = "weak-scaling"

for scenario_params in scenario_generator:
    # Extract parameters from tuple
    comm_scheme, is_communication_overlapped, is_cuda_enabled_mpi, decomposition_axes, num_cells_per_block, fully_periodic = scenario_params
    if comm_scheme != 'UniformGPUScheme_Baseline' and is_cuda_enabled_mpi is True:
        # Skip CUDA enabled MPI tests for GPUPackInfo tests
        continue
    elif comm_scheme == 'GPUPackInfo_Baseline' and is_communication_overlapped is True:
        # Skip communication overlap tests for GPUPackInfo without streams 
        continue

    # Convert the axes decompositions to string
    decomposition_axes_str = ''.join(decomposition_axes)
    # Compute block decomposition based on the specified axes and the number of processes
    blocks = get_block_decomposition(decomposition_axes, num_processes)
    # Estimate number of time steps
    time_steps = max(50, calculate_time_steps(time_per_scenarios, expected_mlups, 3 * (num_cells_per_block,)))
    # Create a benchmark scenario
    scenario = BenchmarkScenario(testcase=testcase_name, decomposition_axes=decomposition_axes_str,
                                 time_steps=time_steps, fully_periodic=fully_periodic)
    # Domain Setup parameters
    domain_setup = scenario.scenario_config['DomainSetup']
    domain_setup['cellsPerBlock'] = 3 * (num_cells_per_block,)
    domain_setup['nrOfProcesses'] = blocks
    domain_setup['blocks'] = blocks
    # Additional parameters for benchmarking
    params = scenario.scenario_config['Parameters']
    params['cudaEnabledMPI'] = is_cuda_enabled_mpi
    params['overlapCommunication'] = is_communication_overlapped
    params['communicationScheme'] = CommunicationSchemeType[comm_scheme]
    # Add scenario for execution
    scenarios.add(scenario)
