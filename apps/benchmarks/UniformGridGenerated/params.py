import math
import os
import operator
import waLBerla as wlb
from waLBerla.tools.sqlitedb import *
from functools import reduce


def prod(seq):
    return reduce(operator.mul, seq, 1)


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


def calculate_time_steps(runtime, expected_mlups, domain_size):
    cells = prod(domain_size)
    time_steps_per_second = expected_mlups * 1e6 / cells
    return int(time_steps_per_second * runtime)


class BenchmarkScenario:
    def __init__(self, block_size=(256, 128, 128), direct_comm=True, time_step_mode='aa', db_file_name='uniform_grid_gen.sqlite'):
        self.block_size = block_size
        self.direct_comm = direct_comm
        self.time_step_mode = time_step_mode
        self.threads = int(os.environ['OMP_NUM_THREADS'])
        self.processes = wlb.mpi.numProcesses()
        self.db_file_name = db_file_name

    @wlb.member_callback
    def config(self, **kwargs):
        time_steps_for_128_cubed = 50
        time_steps = int(128**3 / prod(self.block_size) * time_steps_for_128_cubed)
        time_steps = max(10, time_steps)
        cfg = {
            'DomainSetup': {
                'blocks': (1, 1, self.processes),
                'cellsPerBlock': (self.block_size[0], self.block_size[1], self.block_size[2] * self.threads),
                'periodic': (1, 1, 1),
            },
            'Parameters': {
                'timesteps': time_steps,
                'warmupSteps': 6,
                'outerIterations': 3,
                'vtkWriteFrequency': 0,
                'remainingTimeLoggerFrequency': 0,
                'omega': 1.6,
                'timeStepMode': self.time_step_mode,
                'directComm': self.direct_comm,
            }
        }
        return cfg

    @wlb.member_callback
    def results_callback(self, mlupsPerProcess, optimizations, **kwargs):
        cfg = self.config()
        result = {
            'block_size': self.block_size,
            'mlups_per_core': mlupsPerProcess / self.threads,
            'threads': self.threads,
            'processes': self.processes,
            'time_step_mode': self.time_step_mode,
            'direct_comm': self.direct_comm,
            'time_steps': cfg['Parameters']['timesteps'],
            'I_MPI_PIN_CELL': os.environ.get('I_MPI_PIN_CELL', ''),
            'I_MPI_PIN_DOMAIN': os.environ.get('I_MPI_PIN_CELL', ''),
        }

        optimizations = eval(optimizations)
        result.update(optimizations)
        result.update(kwargs)
        sequenceValuesToScalars(result)
        checkAndUpdateSchema(result, "runs", self.db_file_name)
        storeSingle(result, "runs", self.db_file_name)


def benchmark():
    scenarios = wlb.ScenarioManager()
    for block_size in [(128, 128, 128), (128, 64, 64), (64, 64, 128), (64, 64, 64), (64, 32, 32), (32, 32, 32), (16, 16, 16), (256, 128, 64), (512, 128, 32)]:
        for direct_comm in (True, False):
            for time_step_mode in ['aa', 'aaKernelOnly', 'twoField']:
                sc = BenchmarkScenario(block_size=block_size, direct_comm=direct_comm, time_step_mode=time_step_mode)
                scenarios.add(sc)

benchmark()

