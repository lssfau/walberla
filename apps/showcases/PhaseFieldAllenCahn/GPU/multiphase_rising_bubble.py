import waLBerla as wlb
import waLBerla.tools.sqlitedb as wlbSqlite
import numpy as np

from waLBerla.core_extension import makeSlice
from lbmpy.phasefield_allen_cahn.parameter_calculation import calculate_dimensionless_rising_bubble


class Scenario:
    def __init__(self):
        # output frequencies
        self.vtkWriteFrequency = 1000
        self.dbWriteFrequency = 200

        # simulation parameters
        self.timesteps = 10000
        self.cells = (64, 128, 64)
        self.blocks = (1, 1, 1)
        self.periodic = (0, 0, 0)
        self.size = (self.cells[0] * self.blocks[0],
                     self.cells[1] * self.blocks[1],
                     self.cells[2] * self.blocks[2])

        self.overlappingWidth = (8, 1, 1)
        self.timeStepStrategy = 'normal'

        # physical parameters
        self.bubbleRadius = 16
        self.bubbleMidPoint = (self.size[0] / 2, self.bubbleRadius + 10, self.size[2] / 2)

        # physical parameters
        self.density_heavy = 1.0
        self.reference_time = 18000
        self.parameters = calculate_dimensionless_rising_bubble(reference_time=self.reference_time,
                                                                density_heavy=self.density_heavy,
                                                                bubble_radius=self.bubbleRadius,
                                                                bond_number=1,
                                                                reynolds_number=40,
                                                                density_ratio=1000,
                                                                viscosity_ratio=100)

        # everything else
        self.dbFile = "risingBubble3D.db"

        self.scenario = 1  # 1 rising bubble, 2 RTI

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
                'timeStepStrategy': self.timeStepStrategy,
                'overlappingWidth': self.overlappingWidth,
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
                    {'direction': 'W', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'E', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'T', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'B', 'walldistance': -1, 'flag': 'NoSlip'},
                ]
            },
            'Bubble': {
                'bubbleMidPoint': self.bubbleMidPoint,
                'bubbleRadius': self.bubbleRadius,
                'bubble': True
            },
        }

    @wlb.member_callback
    def at_end_of_time_step(self, blocks, time_loop):
        t = time_loop.getCurrentTimeStep()
        if t % self.dbWriteFrequency == 0:
            wlb_field = wlb.field.gather(blocks, 'phase', makeSlice[:, :, self.size[2] // 2])
            if wlb_field:
                phase_field = np.asarray(wlb_field.buffer()).squeeze()

                location_of_gas = np.where(phase_field < 0.5)
                cov = np.cov(location_of_gas)
                eig = np.linalg.eig(cov)
                axis_of_the_bubble = np.sqrt(eig[0])

                center_of_mass = np.mean(location_of_gas, axis=1)
                self.yPositions.append(center_of_mass[1])
                if len(self.yPositions) > 1:
                    speed = self.yPositions[-1] - self.yPositions[-2]
                    self.write_result_to_database(t, speed, axis_of_the_bubble, center_of_mass)

                self.counter += 1

    def write_result_to_database(self, t, speed, axis_of_the_bubble, center_of_mass):
        """Writes the simulation result stored in the global variables shapeFactors and angles to
               an sqlite3 database, and resets the global variables."""
        result = {'waLBerlaVersion': wlb.build_info.version,
                  'xCells': self.size[0],
                  'yCells': self.size[1],
                  'zCells': self.size[2],
                  'bubbleDiameter': self.bubbleRadius * 2.0,
                  'bubble_axis_x': axis_of_the_bubble[0],
                  'bubble_axis_y': axis_of_the_bubble[1],
                  'center_of_mass_x': center_of_mass[0],
                  'center_of_mass_y': center_of_mass[1],
                  'rising_speed': speed,
                  'timesteps': t,
                  'processes': wlb.mpi.numProcesses(),
                  }
        try:
            wlbSqlite.checkAndUpdateSchema(result, 'data', self.dbFile)
            wlbSqlite.storeSingle(result, 'data', self.dbFile)
        except Exception as e:
            wlb.log_warning("Failed to store run in database " + str(e) + "\n" + str(result))


scenarios = wlb.ScenarioManager()
scenarios.add(Scenario())
