import math
import os
import pickle as pl
import waLBerla as wlb
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from waLBerla.core_extension import makeSlice
from waLBerla.tools.sqlitedb import sequenceValuesToScalars
from scipy.ndimage.filters import gaussian_filter
import scipy.spatial as spatial

from lbmpy.phasefield_allen_cahn.parameter_calculation import calculate_parameters_taylor_bubble


def intersection(points1, points2, eps):
    tree = spatial.KDTree(points1)
    distances, indices = tree.query(points2, k=1, distance_upper_bound=eps)
    intersection_points = tree.data[indices[np.isfinite(distances)]]
    return intersection_points


def contour_points(contour, steps=1):
    return np.row_stack([path.interpolated(steps).vertices
                         for linecol in contour.collections
                         for path in linecol.get_paths()])


def test_line(center_x, center_z, angle, contour, tol):
    contact = -1

    line_size = 200

    points1 = contour_points(contour)
    points2 = np.zeros((line_size, 2))
    points2[:, 0] = center_x + np.linspace(0, center_x, line_size) * np.cos(np.radians(angle))
    points2[:, 1] = center_z + np.linspace(0, center_z, line_size) * np.sin(np.radians(angle))

    intersection_points = intersection(points1, points2, tol)

    if len(intersection_points) != 0:
        contact = 1

    return contact


def get_circle(midpoint, radius):
    theta = np.linspace(0, 2 * np.pi, 100)

    a = midpoint[0] + radius * np.cos(theta)
    b = midpoint[1] + radius * np.sin(theta)

    return a, b


class Scenario:
    def __init__(self, cuda_enabled_mpi=False):
        self.density_liquid = 1.0
        self.reference_time = 16000
        self.reference_length = 128
        d1 = 0.0254  # (0.0508, 0.0381, 0.0254)
        d2 = 0.0127  # (0.0254, 0.0127, 0.0127)
        self.interface_width = 4
        self.mobility = 0.05

        num_processes = wlb.mpi.numProcesses()

        # output frequencies
        self.vtkWriteFrequency = self.reference_time
        self.dbWriteFrequency = self.reference_time // 25
        self.meshWriteFrequency = self.reference_time
        self.pngWriteFrequency = self.reference_time

        # simulation parameters
        self.diameter = self.reference_length
        self.timesteps = self.reference_time * 15 + 1
        self.cells = (self.diameter, (self.reference_length * 15) // num_processes, self.diameter)
        self.blocks = (1, num_processes, 1)
        self.periodic = (0, 0, 0)
        self.size = (self.cells[0] * self.blocks[0],
                     self.cells[1] * self.blocks[1],
                     self.cells[2] * self.blocks[2])
        self.inner_radius = self.diameter // 4

        self.center_x = self.size[0] / 2
        self.center_y = self.size[1] / 2
        self.center_z = self.size[2] / 2

        self.scenario = 4  # 1 rising bubble or droplet, 2 RTI, 3 bubble field, 4 taylor bubble set up

        self.counter = 0
        self.yPositions = []

        self.eccentricity_or_pipe_ratio = False  # if True eccentricity is conducted otherwise pipe ratio
        self.ratio = 0

        self.start_transition = (self.size[1] // 2) - 2 * self.diameter
        self.length_transition = 4 * self.diameter

        setup = "eccentricity" if self.eccentricity_or_pipe_ratio else "ratio"

        self.csv_file = f"Taylor_bubble_D_{self.diameter}_DasC_{setup}_{self.ratio}_W_" \
                        f"{self.interface_width}_M_{self.mobility}.csv"

        d = self.diameter / 2
        dh = self.diameter - d

        resolution = self.diameter / 128

        self.Donut_D = 0.1 * self.diameter / resolution
        self.Donut_h = dh / 6
        self.DonutTime = 0.5 * (self.diameter + d) / 2

        parameters = calculate_parameters_taylor_bubble(reference_length=self.reference_length,
                                                        reference_time=self.reference_time,
                                                        density_heavy=self.density_liquid,
                                                        diameter_outer=d1,
                                                        diameter_inner=d2)

        self.density_gas = parameters["density_light"]
        self.surface_tension = parameters["surface_tension"]

        self.gravitational_acceleration = parameters["gravitational_acceleration"]

        self.relaxation_time_liquid = parameters.get("relaxation_time_heavy")
        self.relaxation_time_gas = parameters.get("relaxation_time_light")

        self.cudaEnabledMpi = cuda_enabled_mpi
        self.cuda_blocks = (64, 2, 2)

        self.config_dict = self.config()

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                'blocks': self.blocks,
                'cellsPerBlock': self.cells,
                'periodic': self.periodic,
                'inner_radius': self.inner_radius,
                'ratio': self.ratio,
                'start_transition': self.start_transition,
                'length_transition': self.length_transition,
                'eccentricity_or_pipe_ration': self.eccentricity_or_pipe_ratio,
                'tube': True
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'dbWriteFrequency': self.dbWriteFrequency,
                'meshWriteFrequency': self.meshWriteFrequency,
                'remainingTimeLoggerFrequency': 60.0,
                'scenario': self.scenario,
                'cudaEnabledMpi': self.cudaEnabledMpi,
                'gpuBlockSize': self.cuda_blocks
            },
            'PhysicalParameters': {
                'density_liquid': self.density_liquid,
                'density_gas': self.density_gas,
                'surface_tension': self.surface_tension,
                'mobility': self.mobility,
                'gravitational_acceleration': self.gravitational_acceleration,
                'relaxation_time_liquid': self.relaxation_time_liquid,
                'relaxation_time_gas': self.relaxation_time_gas,
                'interface_thickness': self.interface_width
            },
            'Boundaries': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'W', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'E', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'T', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'B', 'walldistance': -1, 'flag': 'NoSlip'},
                ],
            },
            'Torus': {
                'Donut_midpoint': 5 * self.diameter,
                'Donut_h': self.Donut_h,
                'Donut_D': self.Donut_D,
                'DonutTime': self.DonutTime
            }
        }

    @wlb.member_callback
    def interface_diffusion(self, blocks):
        for block in blocks:
            phase_field_array = wlb.field.toArray(block['phase'])
            phase_field_array[:, :, :] = gaussian_filter(phase_field_array[:, :, :], sigma=2)

    @wlb.member_callback
    def at_end_of_time_step(self, blocks, **kwargs):
        t = kwargs["timeStep"]
        target_rank = kwargs["target_rank"]
        bounding_box_min = kwargs["bounding_box_min"]
        bounding_box_max = kwargs["bounding_box_max"]
        center_of_mass = [kwargs["center_of_mass_X"], kwargs["center_of_mass_Y"], kwargs["center_of_mass_Z"]]
        if t % self.dbWriteFrequency == 0:
            wlb_field = wlb.field.gather(blocks, 'phase', makeSlice[:, bounding_box_min:bounding_box_max, :],
                                         target_rank)
            if wlb_field:
                data = {'timestep': t}
                data.update(self.config_dict['Parameters'])
                data.update(self.config_dict['DomainSetup'])
                data.update(self.config_dict['Torus'])
                data.update(kwargs)
                del data["bounding_box_min"]
                del data["bounding_box_max"]
                data['total_velocity_Y / sum_inv_phi'] = kwargs['total_velocity_Y'] / kwargs['sum_inv_phi']

                phase_field = np.asarray(wlb_field).squeeze()
                assert np.isfinite(np.sum(phase_field)), "NaN detected in bounding Box"
                location_of_gas = np.where(phase_field < 0.5)
                bubble_tip = np.max(location_of_gas, axis=1)

                fig_handle = plt.figure()

                plt.axis('equal')
                ax = plt.gca()
                ax.set(ylim=(0, self.diameter + 1))
                my_contour = plt.contour(np.rot90(phase_field[:, math.floor(center_of_mass[1]) - bounding_box_min, :]),
                                         [0.5])

                # For eccentricity test cases
                center_x = self.center_x
                center_z = self.center_z
                new_radius = self.inner_radius

                if self.eccentricity_or_pipe_ratio:
                    shift = self.ratio * self.center_x / 2

                    if center_of_mass[1] < self.start_transition:
                        center_x = self.center_x
                    elif self.start_transition < center_of_mass[1] < self.start_transition + self.length_transition:
                        tmp = math.pi * (center_of_mass[1] - self.start_transition) / self.length_transition
                        shift_tmp = shift * 0.5 * (1 - math.cos(tmp))
                        center_x = self.center_x + shift_tmp
                    else:
                        center_x = self.center_x + shift

                else:
                    shift = self.ratio * self.center_x / 2

                    if center_of_mass[1] < self.start_transition:
                        new_radius = self.inner_radius
                    elif self.start_transition < center_of_mass[1] < self.start_transition + self.length_transition:
                        tmp = math.pi * (center_of_mass[1] - self.start_transition) / self.length_transition
                        shift_tmp = shift * 0.5 * (1 - math.cos(tmp))
                        new_radius = self.inner_radius + shift_tmp
                    else:
                        new_radius = self.inner_radius + shift

                start_angle = 0
                tol = 0.5

                # Search for two lines where one intersects and one does not:
                contact = test_line(center_x, center_z, start_angle, my_contour, tol)
                angle = 0
                angle1 = 0

                if contact == -1:
                    num = np.linspace(0, 180, 500)
                    for i in range(0, len(num)):
                        test = test_line(center_x, center_z, num[i], my_contour, tol)
                        if test != -1:
                            angle = num[i]
                            break

                    num = np.linspace(0, -180, 500)
                    for i in range(0, len(num)):
                        test = test_line(center_x, center_z, num[i], my_contour, tol)
                        if test != -1:
                            angle1 = num[i]
                            break

                    theta = 360 - (angle - angle1)
                else:
                    theta = 360

                if self.pngWriteFrequency > 0 and t % self.pngWriteFrequency == 0:
                    plt.plot(center_x + np.linspace(0, center_x) * np.cos(np.radians(angle)),
                             center_z + np.linspace(0, center_z) * np.sin(np.radians(angle)), 'b-')
                    plt.plot(center_x + np.linspace(0, center_x) * np.cos(np.radians(angle1)),
                             center_z + np.linspace(0, center_z) * np.sin(np.radians(angle1)), 'r-')

                    radius = self.diameter // 2
                    circle1 = get_circle([radius, radius], radius)
                    plt.plot(circle1[0], circle1[1], 'k--', linewidth=2)

                    circle2 = get_circle([center_x, center_z], new_radius)
                    plt.fill(circle2[0], circle2[1], 'lightgrey')

                    plt.savefig(f"angle_measurement_ratio_{self.ratio}_{self.counter:06d}.png", dpi=600)

                    outfile = open(f"angle_measurement_ratio_{self.ratio}_{self.counter:06d}.pkl", 'wb')

                    pl.dump(fig_handle, outfile)
                    outfile.close()
                    self.counter += 1
                plt.cla()

                data['bubble_tip_y'] = bubble_tip[1]
                data['center_of_mass_x'] = center_of_mass[0]
                data['center_of_mass_y'] = center_of_mass[1]
                data['center_of_mass_z'] = center_of_mass[2]

                data['xCells'] = self.size[0]
                data['yCells'] = self.size[1]
                data['zCells'] = self.size[2]

                data['theta'] = theta
                if self.eccentricity_or_pipe_ratio:
                    data['eccentricity'] = self.ratio
                else:
                    data['pipe_ratio'] = self.ratio

                sequenceValuesToScalars(data)

                df = pd.DataFrame.from_records([data])
                if not os.path.isfile(self.csv_file):
                    df.to_csv(self.csv_file, index=False)
                else:
                    df.to_csv(self.csv_file, index=False, mode='a', header=False)


scenarios = wlb.ScenarioManager()
scenarios.add(Scenario())
