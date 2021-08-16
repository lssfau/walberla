import os
import waLBerla as wlb
import pandas as pd

from waLBerla.tools.sqlitedb import sequenceValuesToScalars
from waLBerla.tools.config import block_decomposition

import sys


class Scenario:
    def __init__(self, time_step_strategy, cells_per_block=(256, 256, 256),
                 cuda_enabled_mpi=False):
        # output frequencies
        self.vtkWriteFrequency = 0

        # simulation parameters
        self.timesteps = 101
        self.cells_per_block = cells_per_block
        self.blocks = block_decomposition(wlb.mpi.numProcesses())
        self.periodic = (1, 1, 1)
        self.size = (self.cells_per_block[0] * self.blocks[0],
                     self.cells_per_block[1] * self.blocks[1],
                     self.cells_per_block[2] * self.blocks[2])

        self.timeStepStrategy = time_step_strategy
        self.warmupSteps = 10

        self.cudaEnabledMpi = cuda_enabled_mpi

        # bubble parameters
        self.bubbleRadius = min(self.size) // 4
        self.bubbleMidPoint = (self.size[0] / 2, self.size[1] / 2, self.size[2] / 2)

        self.config_dict = self.config()

        self.csv_file = "benchmark.csv"

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                'blocks': self.blocks,
                'cellsPerBlock': self.cells_per_block,
                'periodic': self.periodic,
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'remainingTimeLoggerFrequency': -1,
                'timeStepStrategy': self.timeStepStrategy,
                'warmupSteps': self.warmupSteps,
                'scenario': 1,
                'cudaEnabledMpi': self.cudaEnabledMpi
            },
            'Boundaries_GPU': {
                'Border': []
            },
            'Boundaries_CPU': {
                'Border': []
            },
            'Bubble': {
                'bubbleMidPoint': self.bubbleMidPoint,
                'bubbleRadius': self.bubbleRadius,
            },
        }

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

        df = pd.DataFrame.from_records([data])
        if not os.path.isfile(self.csv_file):
            df.to_csv(self.csv_file, index=False)
        else:
            df.to_csv(self.csv_file, index=False, mode='a', header=False)


def benchmark():
    scenarios = wlb.ScenarioManager()
    block_size = (64, 64, 64)

    scenarios.add(Scenario(time_step_strategy='normal', cells_per_block=block_size))


def kernel_benchmark():
    scenarios = wlb.ScenarioManager()
    block_sizes = [(i, i, i) for i in (8, 16, 32, 64, 128)]

    for time_step_strategy in ['phase_only', 'hydro_only', 'kernel_only', 'normal']:
        for block_size in block_sizes:
            scenario = Scenario(time_step_strategy=time_step_strategy,
                                cells_per_block=block_size)
            scenarios.add(scenario)


# benchmark()
kernel_benchmark()
