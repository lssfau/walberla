import os
import waLBerla as wlb
import pandas as pd

from waLBerla.tools.sqlitedb import sequenceValuesToScalars

import sys


class Scenario:
    def __init__(self, timeStepStrategy, overlappingWidth, cuda_block_size):
        # output frequencies
        self.vtkWriteFrequency = 0

        # simulation parameters
        self.timesteps = 101
        edge_size = 32
        self.cells = (edge_size, edge_size, edge_size)
        self.blocks = (1, 1, 1)
        self.periodic = (1, 1, 1)
        self.size = (self.cells[0] * self.blocks[0],
                     self.cells[1] * self.blocks[1],
                     self.cells[2] * self.blocks[2])

        self.timeStepStrategy = timeStepStrategy
        self.overlappingWidth = overlappingWidth
        self.cuda_block_size = cuda_block_size
        self.warmupSteps = 10

        # bubble parameters
        self.bubbleRadius = min(self.size) // 4
        self.bubbleMidPoint = (self.size[0] / 2, self.size[1] / 2, self.size[2] / 2)

        self.scenario = 1  # 1 rising bubble, 2 RTI
        self.config_dict = self.config()

        self.csv_file = "benchmark.csv"

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
                'remainingTimeLoggerFrequency': -1,
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
                         (4, 4, 4), (8, 8, 8), (16, 16, 16), (32, 32, 32)]

    cuda_blocks = [(32, 1, 1), (64, 1, 1), (128, 1, 1), (256, 1, 1), (512, 1, 1),
                   (32, 2, 1), (64, 2, 1), (128, 2, 1), (256, 2, 1),
                   (32, 4, 1), (64, 4, 1), (128, 4, 1),
                   (32, 4, 2), (64, 4, 2), (128, 2, 2),
                   (32, 8, 1), (64, 8, 1), (64, 4, 2),
                   (32, 16, 1), (16, 16, 1), (16, 16, 2)]

    # no overlap
    scenarios.add(Scenario(timeStepStrategy='normal', overlappingWidth=(1, 1, 1), cuda_block_size=(16, 16, 1)))

    # overlap
    for overlap_strategy in ['overlap']:
        for overlappingWidth in overlappingWidths:
            for cuda_block in cuda_blocks:
                scenario = Scenario(timeStepStrategy=overlap_strategy,
                                    overlappingWidth=overlappingWidth,
                                    cuda_block_size=cuda_block)
                scenarios.add(scenario)


def kernel_benchmark():
    scenarios = wlb.ScenarioManager()

    # overlap
    # for overlap_strategy in ['phase_only', 'hydro_only', 'kernel_only']:
    for overlap_strategy in ['overlap']:
        scenario = Scenario(timeStepStrategy=overlap_strategy,
                            overlappingWidth=(1, 1, 1),
                            cuda_block_size=(128, 1, 1))
        scenarios.add(scenario)


# overlap_benchmark()
kernel_benchmark()
