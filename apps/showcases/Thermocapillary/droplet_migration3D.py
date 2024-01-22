import os
import shutil

import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
import pandas as pd

from lbmpy.phasefield_allen_cahn.parameter_calculation import calculate_parameters_droplet_migration

import waLBerla as wlb
from waLBerla.core_extension import makeSlice
from waLBerla.tools.sqlitedb import sequenceValuesToScalars


class Scenario:
    def __init__(self):
        self.dropletRadius = 32

        self.cells = (6 * self.dropletRadius, 2 * self.dropletRadius, 4 * self.dropletRadius)
        self.blocks = (2, 1, 2)
        self.periodic = (1, 0, 0)
        self.size = (self.cells[0] * self.blocks[0],
                     self.cells[1] * self.blocks[1],
                     self.cells[2] * self.blocks[2])

        self.dropletMidPoint = (65, 0, self.size[2] // 2)

        self.initial_temperature = 0

        self.contact_angle = os.environ.get('CONTACT_ANGLE', 90)
        self.maximum_heat_flux = float(os.environ.get('MAXIMUM_HEAT_FLUX', 0.2))
        single_heat_source = True

        self.heat_source_z_shift = 0  # int(self.dropletRadius * (2.0/3.0))
        self.maximum_heat_flux_1 = self.maximum_heat_flux
        self.maximum_heat_flux_2 = 0 if single_heat_source else self.maximum_heat_flux

        self.capillary_number = 0.01
        self.reynolds_number = 0.16
        self.marangoni_number = 0.08
        self.peclet_number = 1.0

        self.sigma_ref = 5e-3
        self.heat_ratio = 1
        self.interface_width = 4
        self.viscosity_ratio = 1

        self.parameters = calculate_parameters_droplet_migration(radius=self.dropletRadius,
                                                                 reynolds_number=self.reynolds_number,
                                                                 capillary_number=self.capillary_number,
                                                                 marangoni_number=self.marangoni_number,
                                                                 peclet_number=self.peclet_number,
                                                                 viscosity_ratio=self.viscosity_ratio,
                                                                 heat_ratio=self.heat_ratio,
                                                                 interface_width=self.interface_width,
                                                                 reference_surface_tension=self.sigma_ref)

        self.timesteps = 1 + self.parameters.reference_time * 100
        self.pre_thermal_timesteps = self.parameters.reference_time
        self.vtkWriteFrequency = 0
        self.dbWriteFrequency = self.parameters.reference_time

        # use either finite difference (4th order RK scheme) or lattice Boltzmann to solve the heat equation
        self.heat_solver_fd = True  # "lattice Boltzmann"
        self.config_dict = self.config()

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
                'pre_thermal_timesteps': self.pre_thermal_timesteps,
                'heat_solver_fd': self.heat_solver_fd,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'dbWriteFrequency': self.dbWriteFrequency,
                'meshWriteFrequency': 0,
                'remainingTimeLoggerFrequency': 60.0,
                'gpuEnabledMpi': False
            },
            'PhysicalParameters': {
                'density_liquid': self.parameters.density_heavy,
                'density_gas': self.parameters.density_light,
                'sigma_ref': self.parameters.sigma_ref,
                'sigma_t': self.parameters.sigma_t,
                'mobility': self.parameters.mobility,
                'temperature_ref': self.parameters.tmp_ref,
                'heat_conductivity_liquid': self.parameters.heat_conductivity_heavy,
                'heat_conductivity_gas': self.parameters.heat_conductivity_light,
                'relaxation_time_liquid': self.parameters.relaxation_time_heavy,
                'relaxation_time_gas': self.parameters.relaxation_time_light,
                'interface_thickness': self.parameters.interface_thickness,
                'velocity_wall': 0.001,
                'contact_angle': self.contact_angle
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
                    {'direction': 'N', 'walldistance': -1, 'flag': 'UBB'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'T', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'B', 'walldistance': -1, 'flag': 'NoSlip'},
                ],
            },
            'BoundariesThermal': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'DiffusionDirichlet'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'DiffusionDirichlet'},
                    {'direction': 'W', 'walldistance': -1, 'flag': 'NeumannByCopy'},
                    {'direction': 'E', 'walldistance': -1, 'flag': 'NeumannByCopy'},
                    {'direction': 'T', 'walldistance': -1, 'flag': 'DiffusionDirichlet'},
                    {'direction': 'B', 'walldistance': -1, 'flag': 'DiffusionDirichlet'},
                ],
            },
            'Droplet': {
                'dropletRadius': self.dropletRadius,
                'dropletMidPoint': self.dropletMidPoint,
            },
            'HeatSource': {
                'heatSourceMidPointOne': (181, 21, self.size[2] // 2 + self.heat_source_z_shift),
                'heatSourceMidPointTwo': (181, 21, self.size[2] // 2 - self.heat_source_z_shift),
                'ws': 6,
                'sizeDiffusedHotspot': 8,
                'maximumHeatFluxOne': self.maximum_heat_flux_1,
                'maximumHeatFluxTwo': self.maximum_heat_flux_2
            },
        }

    @wlb.member_callback
    def at_end_of_time_step(self, blocks, **kwargs):
        t = kwargs['timeStep']
        velocity_field_wlb = wlb.field.gather(blocks, 'vel', makeSlice[:, :, self.size[2] // 2])
        temperature_field_wlb_xy = wlb.field.gather(blocks, 'temperature', makeSlice[:, :, self.size[2] // 2])
        temperature_field_wlb_xz = wlb.field.gather(blocks, 'temperature', makeSlice[:, 21, :])
        phase_field_wlb = wlb.field.gather(blocks, 'phase', makeSlice[:, :, :])

        if temperature_field_wlb_xy:
            velocity_field = np.asarray(velocity_field_wlb).squeeze()
            temperature_field_xy = np.asarray(temperature_field_wlb_xy).squeeze()
            temperature_field_xz = np.asarray(temperature_field_wlb_xz).squeeze()
            phase_field = np.asarray(phase_field_wlb).squeeze()
            
            var = int(100 * self.maximum_heat_flux)
            path = f"results_contact_angle_{int(self.contact_angle)}_degree_heat_flux_{var}"

            if t == 0:
                if not os.path.exists(path):
                    os.makedirs(path)
                else:
                    shutil.rmtree(path)
                    os.makedirs(path)

            location_of_droplet = np.where(phase_field > 0.5)
            center_of_mass = np.mean(location_of_droplet, axis=1)

            xx, yy = np.meshgrid(np.arange(phase_field.shape[0]), np.arange(phase_field.shape[1]))

            intx = integrate.cumtrapz(velocity_field[:, :, 1].T, xx, axis=1, initial=0)[0]
            inty = integrate.cumtrapz(velocity_field[:, :, 0].T, yy, axis=0, initial=0)
            psi1 = intx-inty
            intx = integrate.cumtrapz(velocity_field[:, :, 1].T, xx, axis=1, initial=0)
            inty = integrate.cumtrapz(velocity_field[:, :, 0].T, yy, axis=0, initial=0)[:, 0][:, None]
            psi2 = intx - inty

            psi = 0.5 * (psi1 + psi2)

            fig, ax = plt.subplots()
            fig.set_figheight(4)
            fig.set_figwidth(32)
            ax.contour(xx, yy, psi, levels=50, cmap="RdBu")
            ax.contour(xx, yy, phase_field[:, :, self.size[2] // 2].T, levels=[0.5, ], colors=['k'], linewidths=[4, ])
            ax.contour(xx, yy, temperature_field_xy[:, :].T, colors=['g'])
            plt.grid()
            plt.xticks(fontsize=40)
            plt.yticks(fontsize=40)
            plt.ylim((0, phase_field.shape[1]))
            plt.xlim((0, phase_field.shape[0]))

            plt.savefig(f'{path}/dropletXY{t}.png', dpi=300)
            plt.close('all')

            xx, zz = np.meshgrid(np.arange(phase_field.shape[0]), np.arange(phase_field.shape[2]))

            fig, ax = plt.subplots()
            fig.set_figheight(16)
            fig.set_figwidth(32)
            ax.contour(xx, zz, phase_field[:, int(center_of_mass[1]), :].T, levels=[0.5, ], colors=['k'], linewidths=[4, ])
            ax.contour(xx, zz, temperature_field_xz[:, :].T, colors=['g'])
            plt.grid()
            plt.xticks(fontsize=40)
            plt.yticks(fontsize=40)
            plt.ylim((0, phase_field.shape[2]))
            plt.xlim((0, phase_field.shape[0]))

            plt.savefig(f'{path}/dropletXZ{t}.png', dpi=300)
            plt.close('all')

            data = {'timestep': t}
            data.update({'total_timesteps': self.timesteps})
            data.update({'pre_thermal_timesteps': self.pre_thermal_timesteps})
            data.update({'center_of_mass_x': center_of_mass[0]})
            data.update({'center_of_mass_y': center_of_mass[1]})
            data.update({'center_of_mass_z': center_of_mass[2]})
            data.update({'capillary_number': self.capillary_number})
            data.update({'reynolds_number': self.reynolds_number})
            data.update({'marangoni_number': self.marangoni_number})
            data.update({'peclet_number': self.peclet_number})
            data.update({'viscosity_ratio': self.viscosity_ratio})
            data.update(self.config_dict['PhysicalParameters'])
            data.update(self.config_dict['Droplet'])
            data.update(self.config_dict['HeatSource'])

            stencil_phase = kwargs.get('stencil_phase')
            stencil_hydro = kwargs.get('stencil_hydro')
            stencil_thermal = kwargs.get('stencil_thermal')
            collision_space_phase = kwargs.get('collision_space_phase')
            collision_space_hydro = kwargs.get('collision_space_hydro')
            collision_space_thermal = kwargs.get('collision_space_thermal')

            data.update({'stencil_phase': stencil_phase})
            data.update({'stencil_hydro': stencil_hydro})
            data.update({'stencil_thermal': stencil_thermal})
            data.update({'collision_space_phase': collision_space_phase})
            data.update({'collision_space_hydro': collision_space_hydro})
            data.update({'collision_space_thermal': collision_space_thermal})
            data.update({'contact_angle': self.contact_angle})

            sequenceValuesToScalars(data)

            csv_file = f"{path}/results.csv"

            df = pd.DataFrame.from_records([data])
            if not os.path.isfile(csv_file):
                df.to_csv(csv_file, index=False, sep=';')
            else:
                df.to_csv(csv_file, index=False, mode='a', header=False, sep=';')


scenarios = wlb.ScenarioManager()
scenarios.add(Scenario())
