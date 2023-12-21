import os
import waLBerla as wlb
from waLBerla.tools.config import block_decomposition
from waLBerla.tools.sqlitedb import sequenceValuesToScalars, checkAndUpdateSchema, storeSingle
import sys
import sqlite3
from math import prod

# Number of time steps run for a workload of 128^3 per GPU
# if double as many cells are on the GPU, half as many time steps are run etc.
# increase this to get more reliable measurements
TIME_STEPS_FOR_128_BLOCK = 1000
DB_FILE = os.environ.get('DB_FILE', "gpu_benchmark.sqlite3")

BASE_CONFIG = {
    'DomainSetup': {
        'cellsPerBlock': (256, 128, 128),
        'periodic': (1, 1, 1),
    },
    'Parameters': {
        'omega': 1.8,
        'cudaEnabledMPI': False,
        'warmupSteps': 5,
        'outerIterations': 3,
    }
}

ldc_setup = {'Border': [
    {'direction': 'N', 'walldistance': -1, 'flag': 'UBB'},
    {'direction': 'W', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'E', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'B', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'T', 'walldistance': -1, 'flag': 'NoSlip'},

]}


def num_time_steps(block_size, time_steps_for_128_block=200):
    cells = block_size[0] * block_size[1] * block_size[2]
    time_steps = (128 ** 3 / cells) * time_steps_for_128_block
    return int(time_steps)


def cuda_block_size_ok(block_size, regs_per_threads=168):
    """Checks if a given CUDA block size does not exceed the SM register limit.
    168 registers per thread was obtained using cuobjdump on both SRT and Cumulant
    kernels. You might want to validate that for your own kernels."""

    return prod(block_size) * regs_per_threads < 64 * (2 ** 10)


def domain_block_size_ok(block_size, total_mem, gls=1, q=27, size_per_value=8):
    """Checks if a single block of given size fits into GPU memory"""
    return prod(b + 2 * gls for b in block_size) * q * size_per_value < total_mem


class Scenario:
    def __init__(self, cells_per_block=(256, 128, 128), periodic=(1, 1, 1), cuda_blocks=(128, 1, 1),
                 timesteps=None, time_step_strategy="normal", omega=1.8, cuda_enabled_mpi=False,
                 inner_outer_split=(1, 1, 1), warmup_steps=5, outer_iterations=3,
                 init_shear_flow=False, boundary_setup=False,
                 vtk_write_frequency=0, remaining_time_logger_frequency=-1,
                 additional_info=None):

        if boundary_setup:
            init_shear_flow = False
            periodic = (0, 0, 0)

        self.blocks = block_decomposition(wlb.mpi.numProcesses())

        self.cells_per_block = cells_per_block
        self.periodic = periodic

        self.time_step_strategy = time_step_strategy
        self.omega = omega
        self.timesteps = timesteps if timesteps else num_time_steps(cells_per_block)
        self.cuda_enabled_mpi = cuda_enabled_mpi
        self.inner_outer_split = inner_outer_split
        self.init_shear_flow = init_shear_flow
        self.boundary_setup = boundary_setup
        self.warmup_steps = warmup_steps
        self.outer_iterations = outer_iterations
        self.cuda_blocks = cuda_blocks

        self.vtk_write_frequency = vtk_write_frequency
        self.remaining_time_logger_frequency = remaining_time_logger_frequency

        self.config_dict = self.config(print_dict=False)
        self.additional_info = additional_info

    @wlb.member_callback
    def config(self, print_dict=True):
        from pprint import pformat
        config_dict = {
            'DomainSetup': {
                'blocks': self.blocks,
                'cellsPerBlock': self.cells_per_block,
                'periodic': self.periodic,
                'oneBlockPerProcess': True
            },
            'Parameters': {
                'omega': self.omega,
                'cudaEnabledMPI': self.cuda_enabled_mpi,
                'warmupSteps': self.warmup_steps,
                'outerIterations': self.outer_iterations,
                'timeStepStrategy': self.time_step_strategy,
                'timesteps': self.timesteps,
                'initShearFlow': self.init_shear_flow,
                'gpuBlockSize': self.cuda_blocks,
                'innerOuterSplit': self.inner_outer_split,
                'vtkWriteFrequency': self.vtk_write_frequency,
                'remainingTimeLoggerFrequency': self.remaining_time_logger_frequency
            },
            'Logging': {
                'logLevel': 'info',  # info progress detail tracing
            }

        }
        if self.boundary_setup:
            config_dict["Boundaries"] = ldc_setup

        if print_dict:
            wlb.log_info_on_root("Scenario:\n" + pformat(config_dict))
            if self.additional_info:
                wlb.log_info_on_root("Additional Info:\n" + pformat(self.additional_info))
        return config_dict

    @wlb.member_callback
    def results_callback(self, **kwargs):
        data = {}
        data.update(self.config_dict['Parameters'])
        data.update(self.config_dict['DomainSetup'])
        data.update(kwargs)

        if self.additional_info is not None:
            data.update(self.additional_info)

        data['executable'] = sys.argv[0]
        data['compile_flags'] = wlb.build_info.compiler_flags
        data['walberla_version'] = wlb.build_info.version
        data['build_machine'] = wlb.build_info.build_machine
        sequenceValuesToScalars(data)

        result = data
        sequenceValuesToScalars(result)
        num_tries = 4
        # check multiple times e.g. may fail when multiple benchmark processes are running
        table_name = f"runs_{data['stencil']}_{data['streamingPattern']}_{data['collisionSetup']}_{prod(self.blocks)}"
        table_name = table_name.replace("-", "_")  # - not allowed for table name would lead to syntax error
        for num_try in range(num_tries):
            try:
                checkAndUpdateSchema(result, table_name, DB_FILE)
                storeSingle(result, table_name, DB_FILE)
                break
            except sqlite3.OperationalError as e:
                wlb.log_warning(f"Sqlite DB writing failed: try {num_try + 1}/{num_tries}  {str(e)}")


# -------------------------------------- Profiling -----------------------------------
def profiling():
    """Tests different communication overlapping strategies"""
    wlb.log_info_on_root("Running 2 timesteps for profiling")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()
    cells = (128, 128, 128)
    cuda_enabled_mpi = False

    scenarios.add(Scenario(cells_per_block=cells, time_step_strategy='kernelOnly',
                           inner_outer_split=(1, 1, 1), timesteps=2, cuda_enabled_mpi=cuda_enabled_mpi,
                           outer_iterations=1, warmup_steps=0))


# -------------------------------------- Functions trying different parameter sets -----------------------------------


def overlap_benchmark():
    """Tests different communication overlapping strategies"""
    wlb.log_info_on_root("Running different communication overlap strategies")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()
    cuda_enabled_mpi = False
    inner_outer_splits = [(1, 1, 1), (4, 1, 1), (8, 1, 1), (16, 1, 1), (32, 1, 1),
                          (4, 4, 1), (8, 8, 1), (16, 16, 1), (32, 32, 1),
                          (4, 4, 4), (8, 8, 8), (16, 16, 16), (32, 32, 32)]

    # no overlap
    scenarios.add(Scenario(time_step_strategy='noOverlap',
                           inner_outer_split=(1, 1, 1),
                           cuda_enabled_mpi=cuda_enabled_mpi,
                           outer_iterations=1))

    for inner_outer_split in inner_outer_splits:
        scenario = Scenario(time_step_strategy='simpleOverlap',
                            inner_outer_split=inner_outer_split,
                            cuda_enabled_mpi=cuda_enabled_mpi,
                            outer_iterations=1)
        scenarios.add(scenario)


def single_gpu_benchmark():
    """Benchmarks only the LBM compute kernel"""
    wlb.log_info_on_root("Running single GPU benchmarks")
    wlb.log_info_on_root("")

    gpu_mem_gb = int(os.environ.get('GPU_MEMORY_GB', 8))
    gpu_mem = gpu_mem_gb * (2 ** 30)
    gpu_type = os.environ.get('GPU_TYPE')

    additional_info = {}
    if gpu_type is not None:
        additional_info['gpu_type'] = gpu_type

    scenarios = wlb.ScenarioManager()
    block_sizes = [(i, i, i) for i in (32, 64, 128, 256)]
    cuda_blocks = [(32, 1, 1), (64, 1, 1), (128, 1, 1), (256, 1, 1), (512, 1, 1),
                   (32, 2, 1), (64, 2, 1), (128, 2, 1), (256, 2, 1),
                   (32, 4, 1), (64, 4, 1), (128, 4, 1),
                   (32, 8, 1), (64, 8, 1),
                   (32, 16, 1)]
    for block_size in block_sizes:
        for cuda_block_size in cuda_blocks:
            # cuda_block_size = (256, 1, 1) and block_size = (64, 64, 64) would be cut to cuda_block_size = (64, 1, 1)
            if cuda_block_size > block_size:
                continue
            if not cuda_block_size_ok(cuda_block_size):
                wlb.log_info_on_root(f"Cuda block size {cuda_block_size} would exceed register limit. Skipping.")
                continue
            if not domain_block_size_ok(block_size, gpu_mem):
                wlb.log_info_on_root(f"Block size {block_size} would exceed GPU memory. Skipping.")
                continue
            scenario = Scenario(cells_per_block=block_size,
                                cuda_blocks=cuda_block_size,
                                time_step_strategy='kernelOnly',
                                timesteps=num_time_steps(block_size, 2000),
                                outer_iterations=1,
                                additional_info=additional_info)
            scenarios.add(scenario)


def validation_run():
    """Run with full periodic shear flow or boundary scenario (ldc) to check if the code works"""
    wlb.log_info_on_root("Validation run")
    wlb.log_info_on_root("")

    time_step_strategy = "noOverlap"  # "simpleOverlap"

    scenarios = wlb.ScenarioManager()
    scenario = Scenario(cells_per_block=(128, 128, 128),
                        time_step_strategy=time_step_strategy,
                        timesteps=10001,
                        outer_iterations=1,
                        warmup_steps=0,
                        init_shear_flow=False,
                        boundary_setup=True,
                        vtk_write_frequency=5000,
                        remaining_time_logger_frequency=30)
    scenarios.add(scenario)


wlb.log_info_on_root(f"Batch run of benchmark scenarios, saving result to {DB_FILE}")
# Select the benchmark you want to run
# single_gpu_benchmark()  # benchmarks different CUDA block sizes and domain sizes and measures single GPU
# performance of compute kernel (no communication)
# overlap_benchmark()  # benchmarks different communication overlap options
# profiling()  # run only two timesteps on a smaller domain for profiling only
validation_run()
