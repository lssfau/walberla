import waLBerla as wlb
from waLBerla.tools.config import block_decomposition
from waLBerla.tools.sqlitedb import sequenceValuesToScalars, checkAndUpdateSchema, storeSingle
import sqlite3
import os
import sys

try:
    import machinestate as ms
except ImportError:
    ms = None

DB_FILE = os.environ.get('DB_FILE', "gpu_benchmark.sqlite3")
BENCHMARK = int(os.environ.get('BENCHMARK', 0))

WeakX = int(os.environ.get('WeakX', 128))
WeakY = int(os.environ.get('WeakY', 128))
WeakZ = int(os.environ.get('WeakZ', 128))

StrongX = int(os.environ.get('StrongX', 128))
StrongY = int(os.environ.get('StrongY', 128))
StrongZ = int(os.environ.get('StrongZ', 128))


class Scenario:
    def __init__(self,
                 domain_size=(64, 64, 64),
                 root_blocks=(2, 2, 2),
                 num_processes=1,
                 refinement_depth=0,
                 cells_per_block=(32, 32, 32),
                 timesteps=101,
                 gpu_enabled_mpi=False,
                 vtk_write_frequency=0,
                 logger_frequency=30,
                 blockforest_filestem="blockforest",
                 write_setup_vtk=True,
                 db_file_name=None):

        self.domain_size = domain_size
        self.root_blocks = root_blocks
        self.cells_per_block = cells_per_block
        self.periodic = (0, 0, 1)

        self.refinement_depth = refinement_depth
        self.num_processes = num_processes
        self.bfs_filestem = blockforest_filestem
        self.write_setup_vtk = write_setup_vtk

        self.timesteps = timesteps
        self.gpu_enabled_mpi = gpu_enabled_mpi
        self.vtk_write_frequency = vtk_write_frequency
        self.logger_frequency = logger_frequency

        self.db_file_name = DB_FILE if db_file_name is None else db_file_name

        self.config_dict = self.config(print_dict=False)

    @wlb.member_callback
    def config(self, print_dict=True):
        from pprint import pformat
        config_dict = {
            'DomainSetup': {
                'domainSize': self.domain_size,
                'rootBlocks': self.root_blocks,
                'cellsPerBlock': self.cells_per_block,
                'periodic': self.periodic,
            },
            'SetupBlockForest': {
                'refinementDepth': self.refinement_depth,
                'numProcesses': self.num_processes,
                'blockForestFilestem': self.bfs_filestem,
                'writeVtk': self.write_setup_vtk,
                'outputStatistics': True,
                'writeSetupForestAndReturn': True,
            },
            'Parameters': {
                'omega': 1.95,
                'timesteps': self.timesteps,
                'remainingTimeLoggerFrequency': self.logger_frequency,
                'vtkWriteFrequency': self.vtk_write_frequency,
                'useVTKAMRWriter': True,
                'oneFilePerProcess': False,
                'writeOnlySlice': False,
                'gpuEnabledMPI': self.gpu_enabled_mpi,
                'gpuBlockSize': (128, 1, 1),
            },
            'Logging': {
                'logLevel': "info",
            }
        }

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
            state.generate()                        # generate subclasses
            state.update()                          # read information
            data["MachineState"] = str(state.get())
        else:
            print("MachineState module is not available. MachineState was not saved")

        sequenceValuesToScalars(data)
        result = data
        sequenceValuesToScalars(result)
        num_tries = 4
        # check multiple times e.g. may fail when multiple benchmark processes are running
        table_name = f"runs"
        table_name = table_name.replace("-", "_")
        for num_try in range(num_tries):
            try:
                checkAndUpdateSchema(result, table_name, self.db_file_name)
                storeSingle(result, table_name, self.db_file_name)
                break
            except sqlite3.OperationalError as e:
                wlb.log_warning(f"Sqlite DB writing failed: try {num_try + 1}/{num_tries}  {str(e)}")


def validation_run():
    """Run with full periodic shear flow or boundary scenario (ldc) to check if the code works"""
    wlb.log_info_on_root("Validation run")

    domain_size = (192, 192, 64)
    cells_per_block = (64, 64, 64)

    root_blocks = tuple([d // c for d, c in zip(domain_size, cells_per_block)])

    scenarios = wlb.ScenarioManager()
    scenario = Scenario(domain_size=domain_size,
                        root_blocks=root_blocks,
                        cells_per_block=cells_per_block,
                        timesteps=0,
                        vtk_write_frequency=0,
                        refinement_depth=3,
                        gpu_enabled_mpi=False)
    scenarios.add(scenario)


def weak_scaling_ldc(num_proc, gpu_enabled_mpi=False, uniform=True):
    wlb.log_info_on_root("Running weak scaling benchmark...")

    # This benchmark must run from 16 GPUs onwards
    if wlb.mpi.numProcesses() > 1:
        num_proc = wlb.mpi.numProcesses()

    if uniform:
        factor = 3 * num_proc
        name = "uniform"
    else:
        if num_proc % 16 != 0:
            raise RuntimeError("Number of processes must be dividable by 16")
        factor = int(num_proc // 16)
        name = "nonuniform"

    cells_per_block = (WeakX, WeakY, WeakZ)
    domain_size = (cells_per_block[0] * 3, cells_per_block[1] * 3, cells_per_block[2] * factor)

    root_blocks = tuple([d // c for d, c in zip(domain_size, cells_per_block)])

    scenarios = wlb.ScenarioManager()
    scenario = Scenario(blockforest_filestem=f"blockforest_{name}_{num_proc}",
                        domain_size=domain_size,
                        root_blocks=root_blocks,
                        num_processes=num_proc,
                        cells_per_block=cells_per_block,
                        refinement_depth=0 if uniform else 3,
                        timesteps=10,
                        gpu_enabled_mpi=gpu_enabled_mpi,
                        db_file_name=f"weakScalingGPU{name}LDC.sqlite3")
    scenarios.add(scenario)


def strong_scaling_ldc(num_proc, gpu_enabled_mpi=False, uniform=True):
    wlb.log_info_on_root("Running strong scaling benchmark...")

    # This benchmark must run from 64 GPUs onwards
    if wlb.mpi.numProcesses() > 1:
        num_proc = wlb.mpi.numProcesses()

    if num_proc % 64 != 0:
        raise RuntimeError("Number of processes must be dividable by 64")

    cells_per_block = (StrongX, StrongY, StrongZ)

    if uniform:
        domain_size = (cells_per_block[0] * 2, cells_per_block[1] * 2, cells_per_block[2] * 16)
        name = "uniform"
    else:
        factor = int(num_proc / 64)
        blocks64 = block_decomposition(factor)
        cells_per_block = tuple([int(c / b) for c, b in zip(cells_per_block, reversed(blocks64))])
        domain_size = (cells_per_block[0] * 3, cells_per_block[1] * 3, cells_per_block[2] * factor)
        name = "nonuniform"

    root_blocks = tuple([d // c for d, c in zip(domain_size, cells_per_block)])

    scenarios = wlb.ScenarioManager()
    scenario = Scenario(blockforest_filestem=f"blockforest_{name}_{num_proc}",
                        domain_size=domain_size,
                        root_blocks=root_blocks,
                        num_processes=num_proc,
                        cells_per_block=cells_per_block,
                        refinement_depth=0 if uniform else 3,
                        timesteps=10,
                        gpu_enabled_mpi=gpu_enabled_mpi,
                        db_file_name=f"strongScalingGPU{name}LDC.sqlite3")
    scenarios.add(scenario)


if BENCHMARK == 0:
    validation_run()
elif BENCHMARK == 1:
    weak_scaling_ldc(1, True, False)
elif BENCHMARK == 2:
    strong_scaling_ldc(1, True, False)
else:
    print(f"Invalid benchmark case {BENCHMARK}")


