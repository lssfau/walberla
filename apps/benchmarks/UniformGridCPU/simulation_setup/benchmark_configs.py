import os
import waLBerla as wlb
from waLBerla.tools.config import block_decomposition
from waLBerla.tools.sqlitedb import sequenceValuesToScalars, checkAndUpdateSchema, storeSingle
import sys
import sqlite3
from math import prod

# Number of time steps run for a workload of 128^3 per process
# if double as many cells are on the process, half as many time steps are run etc.
# increase this to get more reliable measurements
TIME_STEPS_FOR_128_BLOCK = 5
DB_FILE = os.environ.get('DB_FILE', "cpu_benchmark.sqlite3")


def num_time_steps(block_size, time_steps_for_128_block=TIME_STEPS_FOR_128_BLOCK):
    cells = block_size[0] * block_size[1] * block_size[2]
    time_steps = (128 ** 3 / cells) * time_steps_for_128_block
    return int(time_steps)


ldc_setup = {'Border': [
    {'direction': 'N', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'N', 'walldistance': -1, 'flag': 'UBB'},
    {'direction': 'E', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'T', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'B', 'walldistance': -1, 'flag': 'NoSlip'},
]}


class Scenario:
    def __init__(self, cells_per_block=(128, 128, 128), periodic=(1, 1, 1),
                 timesteps=None, time_step_strategy="normal", omega=1.8, inner_outer_split=(1, 1, 1),
                 warmup_steps=2, outer_iterations=3, init_shear_flow=False, boundary_setup=False,
                 vtk_write_frequency=0, remaining_time_logger_frequency=-1):

        if boundary_setup:
            init_shear_flow = False
            periodic = (0, 0, 0)

        self.blocks = block_decomposition(wlb.mpi.numProcesses())

        self.cells_per_block = cells_per_block
        self.periodic = periodic

        self.time_step_strategy = time_step_strategy
        self.omega = omega
        self.timesteps = timesteps if timesteps else num_time_steps(cells_per_block)
        self.inner_outer_split = inner_outer_split
        self.init_shear_flow = init_shear_flow
        self.boundary_setup = boundary_setup
        self.warmup_steps = warmup_steps
        self.outer_iterations = outer_iterations

        self.vtk_write_frequency = vtk_write_frequency
        self.remaining_time_logger_frequency = remaining_time_logger_frequency

        self.config_dict = self.config(print_dict=False)

    @wlb.member_callback
    def config(self, print_dict=True):
        from pprint import pformat
        config_dict = {
            'DomainSetup': {
                'blocks': self.blocks,
                'cellsPerBlock': self.cells_per_block,
                'periodic': self.periodic,
            },
            'Parameters': {
                'omega': self.omega,
                'warmupSteps': self.warmup_steps,
                'outerIterations': self.outer_iterations,
                'timeStepStrategy': self.time_step_strategy,
                'timesteps': self.timesteps,
                'initShearFlow': self.init_shear_flow,
                'innerOuterSplit': self.inner_outer_split,
                'vtkWriteFrequency': self.vtk_write_frequency,
                'remainingTimeLoggerFrequency': self.remaining_time_logger_frequency
            }
        }
        if self.boundary_setup:
            config_dict["Boundaries"] = ldc_setup

        if print_dict:
            wlb.log_info_on_root("Scenario:\n" + pformat(config_dict))
        return config_dict

    @wlb.member_callback
    def results_callback(self, **kwargs):
        data = {}
        data.update(self.config_dict['Parameters'])
        data.update(self.config_dict['DomainSetup'])
        data.update(kwargs)

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
        table_name = table_name.replace("-", "_")
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

    scenarios.add(Scenario(cells_per_block=cells, time_step_strategy='kernelOnly',
                           inner_outer_split=(1, 1, 1), timesteps=2,
                           outer_iterations=1, warmup_steps=0))


# -------------------------------------- Functions trying different parameter sets -----------------------------------


def overlap_benchmark():
    """Tests different communication overlapping strategies"""
    wlb.log_info_on_root("Running different communication overlap strategies")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()
    inner_outer_splits = [(1, 1, 1), (4, 1, 1)]
    cells_per_block = [(i, i, i) for i in (16, 32, 64, 128)]

    for cell_per_block in cells_per_block:
        scenarios.add(Scenario(time_step_strategy='noOverlap',
                               inner_outer_split=(1, 1, 1),
                               cells_per_block=cell_per_block))

        for inner_outer_split in inner_outer_splits:
            scenario = Scenario(time_step_strategy='simpleOverlap',
                                inner_outer_split=inner_outer_split,
                                cells_per_block=cell_per_block)
            scenarios.add(scenario)


def scaling_benchmark():
    """Tests different communication overlapping strategies"""
    wlb.log_info_on_root("Running scaling benchmark")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()
    cells_per_block = [(32, 32, 32), (128, 128, 128)]

    for cell_per_block in cells_per_block:
        scenarios.add(Scenario(time_step_strategy='noOverlap',
                               inner_outer_split=(1, 1, 1),
                               cells_per_block=cell_per_block))


def single_node_benchmark():
    """Benchmarks only the LBM compute kernel"""
    wlb.log_info_on_root("Running single Node benchmarks")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()
    block_sizes = [(i, i, i) for i in (8, 16, 32, 64, 128)]
    for block_size in block_sizes:
        scenario = Scenario(cells_per_block=block_size,
                            time_step_strategy='kernelOnly',
                            timesteps=num_time_steps(block_size))
        scenarios.add(scenario)


def validation_run():
    """Run with full periodic shear flow or boundary scenario (ldc) to check if the code works"""
    wlb.log_info_on_root("Validation run")
    wlb.log_info_on_root("")

    time_step_strategy = 'simpleOverlap'  # 'noOverlap'

    scenarios = wlb.ScenarioManager()
    scenario = Scenario(cells_per_block=(64, 64, 64),
                        time_step_strategy=time_step_strategy,
                        timesteps=101,
                        outer_iterations=1,
                        warmup_steps=0,
                        init_shear_flow=True,
                        boundary_setup=False,
                        vtk_write_frequency=100,
                        remaining_time_logger_frequency=10)
    scenarios.add(scenario)


wlb.log_info_on_root(f"Batch run of benchmark scenarios, saving result to {DB_FILE}")
# Select the benchmark you want to run
single_node_benchmark()  # benchmarks different CUDA block sizes and domain sizes and measures single GPU
# performance of compute kernel (no communication)
# overlap_benchmark()  # benchmarks different communication overlap options
# profiling()  # run only two timesteps on a smaller domain for profiling only
# validation_run()
# scaling_benchmark()
