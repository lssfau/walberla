import waLBerla as wlb
import math


class Scenario:
    def __init__(self):
        # physical parameters
        self.drop_diameter = 20
        self.pool_depth = self.drop_diameter * 0.5

        self.bond_number = 3.18
        self.weber_number = 2010
        self.ohnesorge_number = 0.0384

        self.relaxation_rate_heavy = 1.988
        self.density_heavy = 1
        self.density_ratio = 1000
        self.dynamic_viscosity_ratio = 100

        self.impact_angle_degree = 0    # drop impact angle in degree

        self.interface_thickness = 4
        self.mobility = 0.03

        self.relaxation_time_heavy = 1 / self.relaxation_rate_heavy
        self.kinematic_viscosity_heavy = 1 / 3 * (self.relaxation_time_heavy - 0.5)
        self.dynamic_viscosity_heavy = self.kinematic_viscosity_heavy * self.density_heavy

        self.density_light = self.density_heavy / self.density_ratio
        self.dynamic_viscosity_light = self.dynamic_viscosity_heavy / self.dynamic_viscosity_ratio
        self.kinematic_viscosity_light = self.dynamic_viscosity_light / self.density_light
        self.relaxation_time_light = 3 * self.kinematic_viscosity_light + 0.5

        self.surface_tension = ((self.dynamic_viscosity_heavy / self.ohnesorge_number)**2 / self.drop_diameter
                                / self.density_heavy)

        self.gravitational_acceleration = (- self.bond_number * self.surface_tension / self.drop_diameter**2 /
                                           self.density_heavy)

        self.impact_velocity_magnitude = (self.weber_number * self.surface_tension / self.drop_diameter /
                                          self.density_heavy)**0.5

        self.impact_velocity = (self.impact_velocity_magnitude * math.sin(self.impact_angle_degree * math.pi / 180),
                                -self.impact_velocity_magnitude * math.cos(self.impact_angle_degree * math.pi / 180),
                                0)

        self.reference_time = self.drop_diameter / self.impact_velocity_magnitude
        self.timesteps = 15 * self.reference_time
        self.vtkWriteFrequency = self.reference_time
        self.meshWriteFrequency = self.reference_time

        # domain parameters
        self.size = (self.drop_diameter * 10,
                     self.drop_diameter * 5,
                     self.drop_diameter * 10)
        self.blocks = (1, 1, 1)
        self.periodic = (1, 0, 1)
        self.cells = (self.size[0] // self.blocks[0], self.size[1] // self.blocks[1], self.size[2] // self.blocks[2])
        self.drop_mid_point = (self.size[0] / 2, self.pool_depth + self.drop_diameter / 2, self.size[2] / 2)

        self.scenario = 3   # 1 rising bubble, 2 RTI, 3 drop, 4 taylor bubble set up

        self.counter = 0
        self.yPositions = []

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                'blocks': self.blocks,
                'domainSize': self.size,
                'cellsPerBlock': self.cells,
                'periodic': self.periodic,
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'meshWriteFrequency': self.meshWriteFrequency,
                'remainingTimeLoggerFrequency': 10.0,
                'scenario': self.scenario,
            },
            'PhysicalParameters': {
                'density_liquid': self.density_heavy,
                'density_gas': self.density_light,
                'surface_tension': self.surface_tension,
                'mobility': self.mobility,
                'gravitational_acceleration': self.gravitational_acceleration,
                'relaxation_time_liquid': self.relaxation_time_heavy - 0.5,
                'relaxation_time_gas': self.relaxation_time_light - 0.5,
                'interface_thickness': self.interface_thickness,
            },
            'Boundaries': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
                ]
            },
            'Drop': {
                'drop_mid_point': self.drop_mid_point,
                'drop_radius': self.drop_diameter / 2,
                'pool_depth': self.pool_depth,
                'impact_velocity': self.impact_velocity,
            },
        }


scenarios = wlb.ScenarioManager()
scenarios.add(Scenario())
