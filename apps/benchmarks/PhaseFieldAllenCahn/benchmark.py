import os
import waLBerla as wlb
import pandas as pd

from waLBerla.tools.sqlitedb import sequenceValuesToScalars
from waLBerla.tools.config import block_decomposition

import sys
from math import prod

try:
    import machinestate as ms
except ImportError:
    ms = None


def domain_block_size_ok(block_size, total_mem, gls=1, q_phase=15, q_hydro=27, size_per_value=8):
    """Checks if a single block of given size fits into GPU memory"""

    cells = prod(b + 2 * gls for b in block_size)
    # 3 values for the velocity and two for the phase field and the temporary phase field
    values_per_cell = 2 * q_phase + 2 * q_hydro + 3 + 2
    needed_memory = values_per_cell * cells * size_per_value
    return needed_memory < total_mem


class Scenario:
    def __init__(self, time_step_strategy,
                 cuda_block_size,
                 cells_per_block=(256, 256, 256),
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
        self.cuda_block_size = cuda_block_size
        self.warmupSteps = 10

        self.cudaEnabledMpi = cuda_enabled_mpi

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
                'cellsPerBlock': self.cells_per_block,
                'periodic': self.periodic,
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'remainingTimeLoggerFrequency': -1,
                'timeStepStrategy': self.timeStepStrategy,
                'gpuBlockSize': self.cuda_block_size,
                'warmupSteps': self.warmupSteps,
                'scenario': self.scenario,
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
        if ms:
            state = ms.MachineState(extended=False, anonymous=True)
            state.generate()                        # generate subclasses
            state.update()                          # read information
            data["MachineState"] = str(state.get())
        else:
            print("MachineState module is not available. MachineState was not saved")

        sequenceValuesToScalars(data)

        df = pd.DataFrame.from_records([data])
        if not os.path.isfile(self.csv_file):
            df.to_csv(self.csv_file, index=False)
        else:
            df.to_csv(self.csv_file, index=False, mode='a', header=False)


def benchmark():
    scenarios = wlb.ScenarioManager()

    gpu_mem_gb = int(os.environ.get('GPU_MEMORY_GB', 40))
    gpu_mem = gpu_mem_gb * (2 ** 30)

    block_size = (320, 320, 320)
    cuda_enabled_mpi = True

    if not domain_block_size_ok(block_size, gpu_mem):
        wlb.log_info_on_root(f"Block size {block_size} would exceed GPU memory. Skipping.")
    else:
        scenarios.add(Scenario(time_step_strategy='normal',
                               cuda_block_size=(128, 1, 1),
                               cells_per_block=block_size,
                               cuda_enabled_mpi=cuda_enabled_mpi))


benchmark()
