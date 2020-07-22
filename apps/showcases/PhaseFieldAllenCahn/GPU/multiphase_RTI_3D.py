import waLBerla as wlb
import waLBerla.tools.sqlitedb as wlbSqlite
from waLBerla.core_extension import makeSlice

import numpy as np
from lbmpy.phasefield_allen_cahn.parameter_calculation import calculate_parameters_rti


class Scenario:
    def __init__(self):
        # output frequencies
        self.vtkWriteFrequency = 1000
        self.dbWriteFrequency = 200

        # simulation parameters
        self.timesteps = 27001

        self.cells = (64, 256, 64)
        self.blocks = (1, 1, 1)
        self.periodic = (1, 0, 1)
        self.size = (self.cells[0] * self.blocks[0],
                     self.cells[1] * self.blocks[1],
                     self.cells[2] * self.blocks[2])

        # physical parameters
        self.density_heavy = 1.0
        self.reference_time = 6000
        self.parameters = calculate_parameters_rti(reference_length=128,
                                                   reference_time=self.reference_time,
                                                   density_heavy=self.density_heavy,
                                                   capillary_number=9.1,
                                                   reynolds_number=128,
                                                   atwood_number=1.0,
                                                   peclet_number=140,
                                                   density_ratio=3,
                                                   viscosity_ratio=3)

        # everything else
        self.dbFile = "RTI.db"

        self.scenario = 2  # 1 rising bubble, 2 RTI

        self.counter = 0
        self.yPositions = []

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
                'dbWriteFrequency': self.dbWriteFrequency,
                'useGui': 0,
                'remainingTimeLoggerFrequency': 10.0,
                'scenario': self.scenario,
            },
            'PhysicalParameters': {
                'density_liquid': self.density_heavy,
                'density_gas': self.parameters["density_light"],
                'surface_tension': self.parameters["surface_tension"],
                'mobility': self.parameters.get("mobility", 0.1),
                'gravitational_acceleration': self.parameters["gravitational_acceleration"],
                'relaxation_time_liquid': self.parameters.get("relaxation_time_heavy"),
                'relaxation_time_gas': self.parameters.get("relaxation_time_light"),
            },
            'Boundaries': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
                ]
            },
        }

    @wlb.member_callback
    def at_end_of_time_step(self, blocks, time_loop):
        t = time_loop.getCurrentTimeStep()
        ny = self.size[1]
        l0 = self.size[0]
        if t % self.dbWriteFrequency == 0:
            location_of_spike = -100
            location_of_bubble = -100
            location_of_saddle = -100
            mass = -100
            spike_data = wlb.field.gather(blocks, 'phase', makeSlice[self.size[0] // 2, :, self.size[2] // 2])
            if spike_data:
                spike_field = np.asarray(spike_data.buffer()).squeeze()
                location_of_spike = (np.argmax(spike_field > 0.5) - ny // 2) / l0

            bubble_data = wlb.field.gather(blocks, 'phase', makeSlice[0, :, 0])
            if bubble_data:
                bubble_field = np.asarray(bubble_data.buffer()).squeeze()
                location_of_bubble = (np.argmax(bubble_field > 0.5) - ny // 2) / l0

            saddle_data = wlb.field.gather(blocks, 'phase', makeSlice[0, :, self.size[2] // 2])
            if saddle_data:
                saddle_field = np.asarray(saddle_data.buffer()).squeeze()
                location_of_saddle = (np.argmax(saddle_field > 0.5) - ny // 2) / l0

            phase = wlb.field.gather(blocks, 'phase', makeSlice[:, :, :])
            if phase:
                phase_field = np.asarray(phase.buffer()).squeeze()
                mass = np.sum(phase_field)

            self.write_result_to_database(t, location_of_spike, location_of_bubble, location_of_saddle, mass)

    def write_result_to_database(self, t, spike, bubble, saddle, mass):
        """Writes the simulation result stored in the global variables shapeFactors and angles to
               an sqlite3 database, and resets the global variables."""
        result = {'waLBerlaVersion': wlb.build_info.version,
                  'xCells': self.size[0],
                  'yCells': self.size[1],
                  'zCells': self.size[2],
                  'spike': spike,
                  'bubble': bubble,
                  'saddle': saddle,
                  'mass': mass,
                  'timesteps': t,
                  'normalized_time': t / self.reference_time,
                  'processes': wlb.mpi.numProcesses(),
                  }
        try:
            wlbSqlite.checkAndUpdateSchema(result, 'interface_location', self.dbFile)
            wlbSqlite.storeSingle(result, 'interface_location', self.dbFile)
        except Exception as e:
            wlb.log_warning("Failed to store run in database " + str(e) + "\n" + str(result))


scenarios = wlb.ScenarioManager()
scenarios.add(Scenario())
