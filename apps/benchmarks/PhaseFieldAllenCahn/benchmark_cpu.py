import os
import waLBerla as wlb
import pandas as pd

from waLBerla.tools.sqlitedb import sequenceValuesToScalars

import sys


class Scenario:
    def __init__(self, timeStepStrategy, overlappingWidth):
        # output frequencies
        self.vtkWriteFrequency = -1

        # simulation parameters
        self.timesteps = 500
        edge_size = 50
        self.cells = (edge_size, edge_size, edge_size)
        self.blocks = (1, 1, 1)
        self.periodic = (1, 1, 1)
        self.size = (self.cells[0] * self.blocks[0],
                     self.cells[1] * self.blocks[1],
                     self.cells[2] * self.blocks[2])

        self.timeStepStrategy = timeStepStrategy
        self.overlappingWidth = overlappingWidth
        self.cuda_block_size = (-1, -1, -1)
        self.warmupSteps = 10

        # bubble parameters
        self.bubbleRadius = min(self.size) // 4
        self.bubbleMidPoint = (self.size[0] / 2, self.size[1] / 2, self.size[2] / 2)

        self.scenario = 1  # 1 rising bubble, 2 RTI
        self.config_dict = self.config()

        self.csv_file = "benchmark_cpu.csv"

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                'blocks': self.blocks,
                'cellsPerBlock': self.cells,
                'periodic': self.periodic,
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'useGui': 0,
                'remainingTimeLoggerFrequency': 10.0,
                'timeStepStrategy': self.timeStepStrategy,
                'overlappingWidth': self.overlappingWidth,
                'gpuBlockSize': self.cuda_block_size,
                'warmupSteps': self.warmupSteps,
                'scenario': self.scenario,
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


def overlap_benchmark():
    scenarios = wlb.ScenarioManager()
    overlappingWidths = [(1, 1, 1), (4, 1, 1), (8, 1, 1), (16, 1, 1), (32, 1, 1),
                         (4, 4, 1), (8, 8, 1), (16, 16, 1), (32, 32, 1),
                         (4, 4, 4), (8, 8, 8), (16, 16, 16), (32, 32, 32),
                         (64, 32, 32), (64, 64, 32), (64, 64, 64)]
    # no overlap
    scenarios.add(Scenario(timeStepStrategy='normal', overlappingWidth=(1, 1, 1)))

    # overlap
    for overlap_strategy in ['overlap']:
        for overlappingWidth in overlappingWidths:
            scenario = Scenario(timeStepStrategy=overlap_strategy,
                                overlappingWidth=overlappingWidth)
            scenarios.add(scenario)


def kernel_benchmark():
    scenarios = wlb.ScenarioManager()

    # overlap
    scenario = Scenario(timeStepStrategy='overlap',
                        overlappingWidth=(8, 8, 8))
    scenarios.add(scenario)


# overlap_benchmark()
kernel_benchmark()
