# encoding: utf-8

import itertools
import waLBerla as wlb
from base import get_block_decomposition, communication_schemes, overlap_communication, \
                 cuda_enabled_mpi, num_processes
from benchmark import BenchmarkScenario, CommunicationSchemeType


# Stores the scenarios for the current simulation
scenarios = wlb.ScenarioManager()

# Generates all block decompositions of xyz, 2 directions at a time
#block_decompositions = itertools.combinations_with_replacement('xyz', r=2)
block_decompositions = ['xyz', 'yzx', 'yxz', 'zyx']

cells_per_block = [256,]

if num_processes == 1:
    scenario_generator = itertools.product(communication_schemes, [False,], [False,],
                                           block_decompositions, cells_per_block)
else:
    scenario_generator = itertools.product(communication_schemes, overlap_communication, 
                                           cuda_enabled_mpi, block_decompositions, cells_per_block)

testcase_name = "strong-scaling"

for scenario_params in scenario_generator:
    # Extract parameters from tuple
    comm_scheme, is_communication_overlapped, is_cuda_enabled_mpi, decomposition_axes, num_cells_per_block = scenario_params
    if comm_scheme != 'UniformGPUScheme_Baseline' and is_cuda_enabled_mpi is True:
        # Skip CUDA enabled MPI tests for GPUPackInfo tests
        continue
    elif comm_scheme == 'GPUPackInfo_Baseline' and is_communication_overlapped is True:
        # Skip communication overlap tests for GPUPackInfo baseline
        continue

    # Convert the axes decompositions to string
    decomposition_axes_str = ''.join(decomposition_axes)
    # Compute block decomposition based on the specified axes and the number of processes
    blocks = get_block_decomposition(decomposition_axes, num_processes)
    # Create a benchmark scenario
    scenario = BenchmarkScenario(testcase=testcase_name, decomposition_axes=decomposition_axes_str)
    # Domain Setup parameters
    domain_setup = scenario.scenario_config['DomainSetup']
    domain_setup['cellsPerBlock'] = tuple(num_cells_per_block // block for block in blocks)
    domain_setup['nrOfProcesses'] = blocks
    domain_setup['blocks'] = blocks
    # Additional parameters for benchmarking
    params = scenario.scenario_config['Parameters']
    params['cudaEnabledMPI'] = is_cuda_enabled_mpi
    params['overlapCommunication'] = is_communication_overlapped
    params['communicationScheme'] = CommunicationSchemeType[comm_scheme]
    # Add scenario for execution
    scenarios.add(scenario)

