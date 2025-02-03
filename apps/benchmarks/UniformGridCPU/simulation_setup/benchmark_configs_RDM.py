import os
import waLBerla as wlb
from waLBerla.tools.config import block_decomposition
from waLBerla.tools.sqlitedb import sequenceValuesToScalars, checkAndUpdateSchema, storeSingle
import sys
import sqlite3
from pprint import pformat

try:
    import machinestate as ms
except ImportError:
    ms = None

# Number of time steps run for a workload of 128^3 per process
# if double as many cells are on the process, half as many time steps are run etc.
# increase this to get more reliable measurements
TIME_STEPS_FOR_128_BLOCK = int(os.environ.get('TIME_STEPS_FOR_128_BLOCK', 100))
DB_FILE = os.environ.get('DB_FILE', "cpu_benchmark.sqlite3")
BENCHMARK = int(os.environ.get('BENCHMARK', 0))

WeakX = int(os.environ.get('WeakX', 128))
WeakY = int(os.environ.get('WeakY', 128))
WeakZ = int(os.environ.get('WeakZ', 128))

StrongX = int(os.environ.get('StrongX', 128))
StrongY = int(os.environ.get('StrongY', 128))
StrongZ = int(os.environ.get('StrongZ', 128))


def num_time_steps(block_size, time_steps_for_128_block=TIME_STEPS_FOR_128_BLOCK):
    """
    Calculate the number of time steps based on the block size.

    This function computes the number of time steps required for a given block size
    by scaling the time steps that could be executed on one process within one second 
    for a 128x128x128 cells_per_block to the given cells_per_block size.

    Parameters:
    block_size (tuple): A tuple of three integers representing the dimensions of the cells_per_block (x, y, z).
    time_steps_for_128_block (int, optional): The number of time steps for a 128x128x128 block. Default is 100.

    Returns:
    int: The calculated number of time steps, with a minimum value of 5.
    """
    cells = block_size[0] * block_size[1] * block_size[2]
    time_steps = (128 ** 3 / cells) * time_steps_for_128_block
    if time_steps < 5:
        time_steps = 5
    return int(time_steps)


ldc_setup = {'Border': [
    {'direction': 'N', 'walldistance': -1, 'flag': 'UBB'},
    {'direction': 'W', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'E', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'B', 'walldistance': -1, 'flag': 'NoSlip'},
    {'direction': 'T', 'walldistance': -1, 'flag': 'NoSlip'},
]}


class Scenario:
    def __init__(self, cells_per_block=(128, 128, 128), periodic=(1, 1, 1), blocks_per_process=1,
                 timesteps=None, time_step_strategy="normal", omega=1.8, inner_outer_split=(1, 1, 1),
                 warmup_steps=2, outer_iterations=3, init_shear_flow=False, boundary_setup=False,
                 vtk_write_frequency=0, remaining_time_logger_frequency=-1, db_file_name=None):

        if boundary_setup:
            init_shear_flow = False
            periodic = (0, 0, 0)

        self.blocks_per_process = blocks_per_process
        self.blocks = block_decomposition(self.blocks_per_process * wlb.mpi.numProcesses())

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
        self.db_file_name = DB_FILE if db_file_name is None else db_file_name

        self.config_dict = self.config(print_dict=False)

    @wlb.member_callback
    def config(self, print_dict=True):
        config_dict = {
            'DomainSetup': {
                'blocks': self.blocks,
                'cellsPerBlock': self.cells_per_block,
                'periodic': self.periodic,
                'cartesianSetup': (self.blocks_per_process == 1)
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

        if ms:
            state = ms.MachineState(extended=False, anonymous=True)
            state.generate()  # generate subclasses
            state.update()  # read information
            data["MachineState"] = str(state.get())
        else:
            print("MachineState module is not available. MachineState was not saved")

        sequenceValuesToScalars(data)

        result = data
        sequenceValuesToScalars(result)
        num_tries = 4
        # check multiple times e.g. may fail when multiple benchmark processes are running
        table_name = "runs"
        table_name = table_name.replace("-", "_")
        for num_try in range(num_tries):
            try:
                checkAndUpdateSchema(result, table_name, self.db_file_name)
                storeSingle(result, table_name, self.db_file_name)
                break
            except sqlite3.OperationalError as e:
                wlb.log_warning(f"Sqlite DB writing failed: try {num_try + 1}/{num_tries}  {str(e)}")


# -------------------------------------- Functions trying different parameter sets -----------------------------------


def weak_scaling_benchmark():
    wlb.log_info_on_root("Running weak scaling benchmark with one block per proc")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()

    for t in ["simpleOverlap"]:
        scenarios.add(Scenario(time_step_strategy=t,
                               inner_outer_split=(1, 1, 1),
                               cells_per_block=(WeakX, WeakY, WeakZ),
                               boundary_setup=True,
                               outer_iterations=1,
                               db_file_name="weakScalingUniformGridOneBlock.sqlite3"))


def strong_scaling_benchmark():
    wlb.log_info_on_root("Running strong scaling benchmark with one block per proc")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()

    domain_size = (StrongX, StrongY, StrongZ)
    blocks = block_decomposition(wlb.mpi.numProcesses())
    cells_per_block = tuple([d // b for d, b in zip(domain_size, reversed(blocks))])

    for t in ["simpleOverlap"]:
        scenarios.add(Scenario(cells_per_block=cells_per_block,
                               time_step_strategy=t,
                               outer_iterations=1,
                               timesteps=num_time_steps(cells_per_block),
                               boundary_setup=True,
                               db_file_name="strongScalingUniformGridOneBlock.sqlite3"))


def single_node_benchmark():
    """Benchmarks only the LBM compute kernel"""
    wlb.log_info_on_root("Running single Node benchmarks")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()
    scenario = Scenario(cells_per_block=(128, 128, 128),
                        time_step_strategy='kernelOnly',
                        outer_iterations=1,
                        timesteps=10)
    scenarios.add(scenario)


def validation_run():
    """Run with full periodic shear flow or boundary scenario (ldc) to check if the code works"""
    wlb.log_info_on_root("Validation run")
    wlb.log_info_on_root("")

    time_step_strategy = "noOverlap"  # "noOverlap"

    scenarios = wlb.ScenarioManager()
    scenario = Scenario(cells_per_block=(64, 64, 64),
                        time_step_strategy=time_step_strategy,
                        timesteps=201,
                        outer_iterations=1,
                        warmup_steps=0,
                        init_shear_flow=False,
                        boundary_setup=True,
                        vtk_write_frequency=50,
                        remaining_time_logger_frequency=10)
    scenarios.add(scenario)


wlb.log_info_on_root(f"Batch run of benchmark scenarios, saving result to {DB_FILE}")
# Select the benchmark you want to run
# single_node_benchmark()  # benchmarks different CUDA block sizes and domain sizes and measures single GPU
# performance of compute kernel (no communication)
# overlap_benchmark()  # benchmarks different communication overlap options
# profiling()  # run only two timesteps on a smaller domain for profiling only
# validation_run()
# scaling_benchmark()

if BENCHMARK == 0:
    single_node_benchmark()
elif BENCHMARK == 1:
    weak_scaling_benchmark()
elif BENCHMARK == 2:
    strong_scaling_benchmark()
else:
    validation_run()
