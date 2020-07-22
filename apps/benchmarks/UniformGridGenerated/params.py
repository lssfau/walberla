import math
import os
import operator
import waLBerla as wlb
from waLBerla.tools.sqlitedb import sequenceValuesToScalars, checkAndUpdateSchema, storeSingle
from waLBerla.tools.config import block_decomposition
from functools import reduce
import sqlite3


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


def domain_decomposition_func_z(processes, threads, block_size):
    return {
        'blocks': (1, 1, processes),
        'cellsPerBlock': (block_size[0], block_size[1], block_size[2] * threads)
    }


def domain_decomposition_func_full(processes, threads, block_size):
    return {
        'blocks': block_decomposition(processes),
        'cellsPerBlock': (block_size[0], block_size[1], block_size[2] * threads)
    }


class BenchmarkScenario:
    def __init__(self, block_size=(256, 128, 128), direct_comm=True,
                 time_step_mode='aa', two_field_kernel_type='generated',
                 domain_decomposition_func=domain_decomposition_func_z,
                 db_file_name='uniform_grid_gen.sqlite'):
        self.block_size = block_size
        self.direct_comm = direct_comm
        self.time_step_mode = time_step_mode
        self.two_field_kernel_type = two_field_kernel_type
        self.domain_decomposition_func = domain_decomposition_func
        self.threads = int(os.environ['OMP_NUM_THREADS'])
        self.processes = wlb.mpi.numProcesses()
        self.db_file_name = db_file_name

    @wlb.member_callback
    def config(self, **kwargs):
        time_steps_for_128_cubed = 10
        time_steps = int(128**3 / prod(self.block_size) * time_steps_for_128_cubed)
        time_steps = max(10, time_steps)
        cfg = {
            'DomainSetup': {
                'periodic': (1, 1, 1),
            },
            'Parameters': {
                'timesteps': time_steps,
                'warmupSteps': 2,
                'outerIterations': 4,
                'vtkWriteFrequency': 0,
                'remainingTimeLoggerFrequency': 0,
                'omega': 1.6,
                'timeStepMode': self.time_step_mode,
                'twoFieldKernelType': self.two_field_kernel_type,
                'directComm': self.direct_comm,
            }
        }
        cfg['DomainSetup'].update(self.domain_decomposition_func(self.processes, self.threads, self.block_size))
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
        num_tries = 4
        for num_try in range(num_tries):
            try:
                checkAndUpdateSchema(result, "runs", self.db_file_name)
                storeSingle(result, "runs", self.db_file_name)
                break
            except sqlite3.OperationalError as e:
                wlb.log_warning(f"Sqlite DB writing failed: try {num_try + 1}/4 " + str(e))


def block_size_ok(sc):
    block_size = sc.config()['DomainSetup']['cellsPerBlock']
    return prod(block_size) * 19 < 2 ** 31


def single_node_benchmark():
    scenarios = wlb.ScenarioManager()
    for block_size in [(128, 128, 128), (128, 64, 64), (64, 64, 128), (64, 128, 64), (64, 64, 64),
                       (1024, 64, 32), (2048, 64, 16),
                       (64, 32, 32), (32, 32, 32), (16, 16, 16), (256, 128, 64), (512, 128, 32)]:
        for direct_comm in (True, False):
            for time_step_mode in ['aa', 'aaKernelOnly', 'twoField']:
                if time_step_mode == 'twoField':
                    for kt in ['generated', 'manualGeneric', 'manualD3Q19']:
                        sc = BenchmarkScenario(block_size=block_size, direct_comm=direct_comm,
                                               time_step_mode=time_step_mode, two_field_kernel_type=kt,
                                               domain_decomposition_func=domain_decomposition_func_z
                                               )
                        if not block_size_ok(sc):
                            continue
                        scenarios.add(sc)
                else:
                    sc = BenchmarkScenario(block_size=block_size, direct_comm=direct_comm,
                                           domain_decomposition_func=domain_decomposition_func_z,
                                           time_step_mode=time_step_mode)
                    if not block_size_ok(sc):
                        continue
                        # scenarios.add(sc)


def single_node_benchmark_small():
    scenarios = wlb.ScenarioManager()
    for block_size in [(128, 128, 128), (128, 64, 64), (64, 64, 128), (64, 128, 64), (64, 64, 64),
                       (1024, 64, 32), (2048, 64, 16), (64, 32, 32), (32, 32, 32), (16, 16, 16),
                       (256, 128, 64), (512, 128, 32)]:
        for direct_comm in (True, False):
            for time_step_mode in ['aa', 'aaKernelOnly', 'twoField']:
                sc = BenchmarkScenario(block_size=block_size, direct_comm=direct_comm, time_step_mode=time_step_mode)
                if not block_size_ok(sc):
                    continue
                scenarios.add(sc)


def weak_scaling():
    scenarios = wlb.ScenarioManager()
    for block_size in [(128, 128, 128), (128, 64, 64), (64, 64, 128), (64, 128, 64), (64, 64, 64),
                       (1024, 64, 32), (2048, 64, 16), (64, 32, 32), (32, 32, 32), (16, 16, 16),
                       (256, 128, 64), (512, 128, 32)]:
        for direct_comm in (True, False):
            sc = BenchmarkScenario(block_size=block_size, direct_comm=direct_comm,
                                   domain_decomposition_func=domain_decomposition_func_full)
            if not block_size_ok(sc):
                continue
            scenarios.add(sc)


single_node_benchmark()
