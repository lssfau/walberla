import os
import sys
import pandas as pd

import waLBerla as wlb
from waLBerla.tools.config import block_decomposition
from waLBerla.tools.sqlitedb import sequenceValuesToScalars


def num_time_steps(block_size, time_steps_for_256_block=50):
    # Number of time steps run for a workload of 256^3 cells per process
    # if double as many cells are on the process, half as many time steps are run etc.
    # increase this to get more reliable measurements
    cells = block_size[0] * block_size[1] * block_size[2]
    time_steps = (256 ** 3 / cells) * time_steps_for_256_block
    if time_steps < 10:
        time_steps = 10
    return int(time_steps)


class Scenario:
    def __init__(self, cells=(256, 256, 256), heat_solver_rk_or_lbm=False,
                 weak_scaling=False, time_step_strategy='NormalTimestep',
                 cuda_enabled_mpi=False, cuda_blocks=(256, 1, 1), timesteps=None, rk_solver_order=2):
        # simulation parameters
        self.timesteps = timesteps if timesteps else num_time_steps(cells)

        self.cells = cells
        self.blocks = block_decomposition(wlb.mpi.numProcesses())
        self.periodic = (1, 0, 0)

        self.heat_solver_rk_or_lbm = heat_solver_rk_or_lbm
        self.weak_scaling = weak_scaling
        self.time_step_strategy = time_step_strategy
        self.rk_solver_order = rk_solver_order

        # GPU specific parameters will be neglected for benchmarks on CPU
        self.cuda_enabled_mpi = cuda_enabled_mpi
        self.cuda_blocks = cuda_blocks

        self.config_dict = self.config()

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                'blocks': self.blocks,
                'cellsPerBlock': self.cells,
                'periodic': self.periodic,
                'weakScaling': self.weak_scaling,
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'pre_thermal_timesteps': 0,
                'HeatSolverRKOrLBM': self.heat_solver_rk_or_lbm,
                'orderRKSolver': self.rk_solver_order,
                'vtkWriteFrequency': 0,
                'dbWriteFrequency': 0,
                'gpuBlockSize': self.cuda_blocks,
                'gpuEnabledMpi': self.cuda_enabled_mpi
            },
            'PhysicalParameters': {
                'density_liquid': 1.0,
                'density_gas': 1.0,
                'sigma_ref': 0.025,
                'sigma_t': -5e-4,
                'mobility': 0.05,
                'temperature_ref': 10,
                'heat_conductivity_liquid': 0.2,
                'heat_conductivity_gas': 0.2,
                'relaxation_time_liquid': 1.1,
                'relaxation_time_gas': 1.1,
                'interface_thickness': 4,
                'velocity_wall': 0
            },
            'BoundariesAllenCahn': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'T', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'B', 'walldistance': -1, 'flag': 'NoSlip'},
                ],
            },
            'BoundariesHydro': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'T', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'B', 'walldistance': -1, 'flag': 'NoSlip'},
                ],
            },
            'BoundariesThermal': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'DiffusionDirichletStatic'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'DiffusionDirichletDynamic'},
                    {'direction': 'T', 'walldistance': -1, 'flag': 'NeumannByCopy'},
                    {'direction': 'B', 'walldistance': -1, 'flag': 'NeumannByCopy'},
                ],
            },
            'Benchmark': {
                'timeStepStrategy': self.time_step_strategy
            }
        }

    @wlb.member_callback
    def write_benchmark_results(self, **kwargs):
        data = {'timesteps': self.timesteps,
                'cells_per_block_0': self.cells[0],
                'cells_per_block_1': self.cells[1],
                'cells_per_block_2': self.cells[2],
                'blocks_0': self.blocks[0],
                'blocks_1': self.blocks[1],
                'blocks_2': self.blocks[2],
                'cuda_blocks_0': self.cuda_blocks[0],
                'cuda_blocks_1': self.cuda_blocks[1],
                'cuda_blocks_2': self.cuda_blocks[2],
                'stencil_phase': kwargs.get('stencil_phase'),
                'stencil_hydro': kwargs.get('stencil_hydro'),
                'stencil_thermal': kwargs.get('stencil_thermal'),
                'collision_space_phase': kwargs.get('collision_space_phase'),
                'collision_space_hydro': kwargs.get('collision_space_hydro'),
                'collision_space_thermal': kwargs.get('collision_space_thermal'),
                'field_type': kwargs.get('field_type'),
                'field_type_pdfs': kwargs.get('field_type_pdfs'),
                'number_of_processes': kwargs.get('number_of_processes'),
                'threads': kwargs.get('threads'),
                'MLUPS': kwargs.get('MLUPS'),
                'MLUPS_process': kwargs.get('MLUPS_process'),
                'timeStepsPerSecond': kwargs.get('timeStepsPerSecond'),
                'timeStepStrategy': self.time_step_strategy,
                'thermalSolver': "Runge Kutta" if self.heat_solver_rk_or_lbm else "lattice Boltzmann",
                'RKSolverOrder': self.rk_solver_order}

        data.update(self.config_dict['PhysicalParameters'])
        data['executable'] = sys.argv[0]
        data['compile_flags'] = wlb.build_info.compiler_flags
        data['walberla_version'] = wlb.build_info.version
        data['build_machine'] = wlb.build_info.build_machine

        sequenceValuesToScalars(data)

        csv_file = f"thermocapillary_benchmark.csv"

        df = pd.DataFrame.from_records([data])
        if not os.path.isfile(csv_file):
            df.to_csv(csv_file, index=False, sep=';')
        else:
            df.to_csv(csv_file, index=False, mode='a', header=False, sep=';')


def profiling():
    scenarios = wlb.ScenarioManager()
    cuda_enabled_mpi = False
    time_step_strategy = "PhaseOnly"  # , "HydroOnly", "ThermalOnly", "KernelOnly"
    heat_solver_rk_or_lbm = True

    cells_per_block = (256, 256, 256)
    cuda_blocks = (128, 1, 1)
    scenarios.add(Scenario(cells=cells_per_block,
                           heat_solver_rk_or_lbm=heat_solver_rk_or_lbm,
                           time_step_strategy=time_step_strategy,
                           cuda_enabled_mpi=cuda_enabled_mpi,
                           cuda_blocks=cuda_blocks,
                           timesteps=2))


def benchmark():
    scenarios = wlb.ScenarioManager()
    cuda_enabled_mpi = True
    cells = (512, 256, 256)

    for i in range(3):
        scenarios.add(Scenario(cells=cells,
                               heat_solver_rk_or_lbm=False,
                               rk_solver_order=2,
                               time_step_strategy="NormalTimestep",
                               cuda_enabled_mpi=cuda_enabled_mpi,
                               cuda_blocks=(64, 2, 1),
                               timesteps=100))

    # for rk_solver_order in (2, 4):
    #     scenarios.add(Scenario(cells=cells,
    #                            heat_solver_rk_or_lbm=True,
    #                            rk_solver_order=rk_solver_order,
    #                            time_step_strategy="NormalTimestep",
    #                            cuda_enabled_mpi=cuda_enabled_mpi,
    #                            cuda_blocks=(64, 2, 1),
    #                            timesteps=100))


def single_kernel_benchmark():
    scenarios = wlb.ScenarioManager()
    cuda_enabled_mpi = False

    cells_per_block = [(i, i, i) for i in (384, )]
    cuda_blocks = [(64, 1, 1), (128, 1, 1), (256, 1, 1),
                   (64, 2, 1), (128, 2, 1),
                   (64, 4, 1)]
    for heat_solver_rk_or_lbm in (True, ):
        for time_step_strategy in ("PhaseOnly", "HydroOnly", "ThermalOnly", "KernelOnly"):
            for cells in cells_per_block:
                for cuda_block_size in cuda_blocks:
                    scenarios.add(Scenario(cells=cells,
                                           heat_solver_rk_or_lbm=heat_solver_rk_or_lbm,
                                           rk_solver_order=2,
                                           time_step_strategy=time_step_strategy,
                                           cuda_enabled_mpi=cuda_enabled_mpi,
                                           cuda_blocks=cuda_block_size))


# profiling()
benchmark()
# single_kernel_benchmark()
