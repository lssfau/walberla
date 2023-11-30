import waLBerla as wlb
from waLBerla.tools.sqlitedb import sequenceValuesToScalars, checkAndUpdateSchema, storeSingle
import sqlite3
import os
import sys

DB_FILE = os.environ.get('DB_FILE', "cpu_benchmark.sqlite3")

class Scenario:
    def __init__(self,
                 domain_size=(64, 64, 64),
                 root_blocks=(2, 2, 2),
                 num_processes=1,
                 refinement_depth=0,
                 cells_per_block=(32, 32, 32),
                 timesteps=101,
                 vtk_write_frequency=0,
                 logger_frequency=0,
                 blockforest_filestem="blockforest",
                 write_setup_vtk=False):

        self.domain_size = domain_size
        self.root_blocks = root_blocks
        self.cells_per_block = cells_per_block
        self.periodic = (0, 0, 0)

        self.refinement_depth = refinement_depth
        self.num_processes = num_processes
        self.bfs_filestem = blockforest_filestem
        self.write_setup_vtk = write_setup_vtk

        self.timesteps = timesteps
        self.vtk_write_frequency = vtk_write_frequency
        self.logger_frequency = logger_frequency

        self.config_dict = self.config(print_dict=False)

    @wlb.member_callback
    def config(self, print_dict=True):
        from pprint import pformat
        config_dict = {
            'DomainSetup': {
                'domainSize': self.domain_size,
                'rootBlocks': self.root_blocks,
                'cellsPerBlock': self.cells_per_block,
                'periodic': self.periodic
            },
            'SetupBlockForest': {
                'refinementDepth': self.refinement_depth,
                'numProcesses': self.num_processes,
                'blockForestFilestem': self.bfs_filestem,
                'writeVtk': self.write_setup_vtk,
                'outputStatistics': False
            },
            'Parameters': {
                'omega': 1.95,
                'timesteps': self.timesteps,
                'remainingTimeLoggerFrequency': self.logger_frequency,
                'vtkWriteFrequency': self.vtk_write_frequency,
            },
            'Logging': {
                'logLevel': "info",
            }
        }

        if(print_dict):
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
        table_name = f"runs"
        table_name = table_name.replace("-", "_")
        for num_try in range(num_tries):
            try:
                checkAndUpdateSchema(result, table_name, DB_FILE)
                storeSingle(result, table_name, DB_FILE)
                break
            except sqlite3.OperationalError as e:
                wlb.log_warning(f"Sqlite DB writing failed: try {num_try + 1}/{num_tries}  {str(e)}")


def validation_run():
    """Run with full periodic shear flow or boundary scenario (ldc) to check if the code works"""
    wlb.log_info_on_root("Validation run")

    domain_size = (96, 96, 96)
    cells_per_block = (32, 32, 32)

    root_blocks = tuple([d // c for d, c in zip(domain_size, cells_per_block)])

    scenarios = wlb.ScenarioManager()
    scenario = Scenario(domain_size=domain_size,
                        root_blocks=root_blocks,
                        num_processes=1,
                        refinement_depth=1,
                        cells_per_block=cells_per_block,
                        timesteps=201,
                        vtk_write_frequency=100,
                        logger_frequency=5,
                        write_setup_vtk=True)
    scenarios.add(scenario)

def scaling():
    wlb.log_info_on_root("Running scaling benchmark...")

    numProc = wlb.mpi.numProcesses()

    domain_size = (256, 256, 128 * numProc)
    cells_per_block = (64, 64, 64)
    root_blocks = tuple([d // c for d, c in zip(domain_size, cells_per_block)])

    scenarios = wlb.ScenarioManager()
    scenario = Scenario(domain_size=domain_size,
                        root_blocks=root_blocks,
                        cells_per_block=cells_per_block,
                        refinement_depth=2,
                        timesteps=10)
    scenarios.add(scenario)

validation_run()
# scaling()
