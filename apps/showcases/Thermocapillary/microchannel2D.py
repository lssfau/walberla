import os
import shutil

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from lbmpy.phasefield_allen_cahn.analytical import analytical_solution_microchannel

import waLBerla as wlb
from waLBerla.core_extension import makeSlice
from waLBerla.tools.sqlitedb import sequenceValuesToScalars


class Scenario:
    def __init__(self, heat_solver_rk_or_lbm, case):
        self.reference_length = 256

        # simulation parameters
        self.timesteps = 1
        self.pre_thermal_timesteps = 1
        self.vtkWriteFrequency = 0
        self.dbWriteFrequency = 10000

        self.cells = (2 * self.reference_length, self.reference_length, 1)
        self.blocks = (1, 1, 1)
        self.domain_size = tuple([c * b for c, b in zip(self.cells, self.blocks)])
        self.periodic = (1, 0, 0)

        self.heat_solver_rk_or_lbm = heat_solver_rk_or_lbm
        self.order_rk_solver = 4

        self.case = case

        if self.case == 1:
            self.rho_h = 1.0
            self.rho_l = 1.0
            self.mu_l = 0.2
            self.mu_h = 0.2

            self.sigma_ref = 0.025
            self.sigma_t = -5e-4

            self.kappa_h = 0.2
            self.kappa_l = 0.2
            self.temperature_ref = 10
            self.temperature_h = 20
            self.temperature_l = 4

            self.interface_thickness = 5
            self.mobility = 0.05
            self.velocity_wall = 0
        else:
            self.rho_h = 1.0
            self.rho_l = 1.0
            self.mu_l = 0.2
            self.mu_h = 0.2

            self.sigma_ref = 0.025
            self.sigma_t = -5e-4

            self.kappa_h = 0.04
            self.kappa_l = 0.2
            self.temperature_ref = 10
            self.temperature_h = 20
            self.temperature_l = 4

            self.interface_thickness = 5
            self.mobility = 0.05
            self.velocity_wall = 0

        x, y, u_x, u_y, t_a = analytical_solution_microchannel(self.reference_length,
                                                               self.domain_size[0], self.domain_size[1],
                                                               self.kappa_h, self.kappa_l,
                                                               self.temperature_h, self.temperature_ref,
                                                               self.temperature_l,
                                                               self.sigma_t, self.mu_l)
        self.x = x
        self.y = y
        self.u_x = u_x
        self.u_y = u_y
        self.t_a = t_a
        l0 = self.reference_length
        self.XX, self.YY = np.meshgrid(np.arange(self.domain_size[0]) - l0, np.arange(self.domain_size[1]) - l0 // 2)

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
                'HeatSolverRKOrLBM': self.heat_solver_rk_or_lbm,
                'orderRKSolver': self.order_rk_solver,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'dbWriteFrequency': self.dbWriteFrequency,
                'remainingTimeLoggerFrequency': 30.0,
                'gpuBlockSize': (128, 1, 1),
                'gpuEnabledMpi': False
            },
            'PhysicalParameters': {
                'density_liquid': self.rho_h,
                'density_gas': self.rho_l,
                'sigma_ref': self.sigma_ref,
                'sigma_t': self.sigma_t,
                'mobility': self.mobility,
                'temperature_ref': self.temperature_ref,
                'heat_conductivity_liquid': self.kappa_h,
                'heat_conductivity_gas': self.kappa_l,
                'relaxation_time_liquid': 3 * self.mu_h,
                'relaxation_time_gas': 3 * self.mu_l,
                'interface_thickness': self.interface_thickness,
                'velocity_wall': self.velocity_wall
            },
            'BoundariesAllenCahn': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
                ],
            },
            'BoundariesHydro': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
                ],
            },
            'BoundariesThermal': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'DiffusionDirichletStatic'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'DiffusionDirichletDynamic'},
                ],
            },
            'MicroChannel': {
                'Th': self.temperature_h,
                'T0': self.temperature_l,
            }
        }

    @wlb.member_callback
    def at_end_of_time_step(self, blocks, **kwargs):
        t = kwargs['timeStep']
        velocity_field_wlb = wlb.field.gather(blocks, 'vel', makeSlice[:, :, :])
        temperature_field_wlb = wlb.field.gather(blocks, 'temperature', makeSlice[:, :, :])

        if temperature_field_wlb:
            stencil_phase = kwargs.get('stencil_phase')
            stencil_hydro = kwargs.get('stencil_hydro')
            stencil_thermal = kwargs.get('stencil_thermal')
            collision_space_phase = kwargs.get('collision_space_phase')
            collision_space_hydro = kwargs.get('collision_space_hydro')
            collision_space_thermal = kwargs.get('collision_space_thermal')
            field_type = kwargs.get('field_type')
            field_type_pdfs = kwargs.get('field_type_pdfs')
            if self.heat_solver_rk_or_lbm:
                path = f"results_{stencil_phase}_{stencil_hydro}_{collision_space_phase}_{collision_space_hydro}_RK{self.order_rk_solver}_{field_type}_{field_type_pdfs}_case_{self.case}"
            else:
                path = f"results_{stencil_phase}_{stencil_hydro}_{stencil_thermal}_" \
                       f"{collision_space_phase}_{collision_space_hydro}_{collision_space_thermal}_{field_type}_{field_type_pdfs}_case_{self.case}"

            if t == 0:
                if not os.path.exists(path):
                    os.makedirs(path)
                else:
                    shutil.rmtree(path)
                    os.makedirs(path)

            velocity_field = np.asarray(velocity_field_wlb).squeeze()
            temperature_field = np.asarray(temperature_field_wlb).squeeze()
            l_2_temp_num = 0.0
            l_2_temp_den = 0.0
            l_2_u_num = 0.0
            l_2_u_den = 0.0
            from math import sqrt
            for x in range(self.domain_size[0]):
                for y in range(self.domain_size[1]):
                    u = sqrt(velocity_field[x, y, 0] ** 2 + velocity_field[x, y, 1] ** 2)
                    u_den = sqrt(self.u_x.T[x, y] ** 2 + self.u_y.T[x, y] ** 2)
                    l_2_temp_num += ((temperature_field[x, y] - self.t_a.T[x, y]) ** 2)
                    l_2_temp_den += ((self.t_a.T[x, y]) ** 2)
                    l_2_u_num += ((u - u_den) ** 2)
                    l_2_u_den += (u_den ** 2)
            l_2_t = sqrt(l_2_temp_num / l_2_temp_den)
            l_2_u = sqrt(l_2_u_num / l_2_u_den)

            fig, ax = plt.subplots()
            fig.set_figheight(5)
            fig.set_figwidth(10)
            levels = range(11, 24)
            plt.contour(self.x, self.y, self.t_a, linestyles='dashed', levels=levels, colors=['k'])
            plt.grid()
            contour_2 = plt.contour(self.XX, self.YY, temperature_field.T, levels=levels, colors=['k'])
            clabels = ax.clabel(contour_2, inline=1, fontsize=10, fmt='%2.0lf')
            [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabels]
            plt.ylim((-128, 128))
            plt.xlim((-256, 256))
            plt.savefig(f'{path}/temperature_{t}.png', dpi=300)
            plt.close('all')

            fig1, ax = plt.subplots()
            fig1.set_figheight(5)
            fig1.set_figwidth(10)
            n = 15
            ax.quiver(self.x[::n, ::n] + 1.1, self.y[::n, ::n] - 2.5, self.u_x.T[::n, ::n].T, self.u_y.T[::n, ::n].T,
                      angles='xy', scale_units='xy', scale=0.00001, color='r')
            # c = np.sqrt(velocity_field[::n, ::n, 0].T * velocity_field[::n, ::n, 0].T +
            #             velocity_field[::n, ::n, 1].T * velocity_field[::n, ::n, 1].T)
            ax.quiver(self.XX[::n, ::n] + 1.1, self.YY[::n, ::n] - 2,
                      velocity_field[::n, ::n, 0].T, velocity_field[::n, ::n, 1].T,
                      angles='xy', scale_units='xy', scale=0.00001, color='b')
            plt.grid()
            plt.ylim((-128, 128))
            plt.xlim((-256, 256))
            plt.savefig(f'{path}/velocity_{t}.png', dpi=300)
            plt.close('all')

            data = {'timestep': t}
            data.update({"L2_T": l_2_t})
            data.update({"L2_U": l_2_u})
            data.update({'total_timesteps': self.timesteps})
            data.update({'pre_thermal_timesteps': self.pre_thermal_timesteps})
            data.update({'stencil_phase': stencil_phase})
            data.update({'stencil_hydro': stencil_hydro})
            if not self.heat_solver_rk_or_lbm:
                data.update({'stencil_thermal': stencil_thermal})
            data.update({'collision_space_phase': collision_space_phase})
            data.update({'collision_space_hydro': collision_space_hydro})
            if not self.heat_solver_rk_or_lbm:
                data.update({'collision_space_thermal': collision_space_thermal})
            data.update(self.config_dict['PhysicalParameters'])

            sequenceValuesToScalars(data)

            csv_file = f"{path}/results.csv"

            df = pd.DataFrame.from_records([data])
            if not os.path.isfile(csv_file):
                df.to_csv(csv_file, index=False, sep=';')
            else:
                df.to_csv(csv_file, index=False, mode='a', header=False, sep=';')


scenarios = wlb.ScenarioManager()
# scenarios.add(Scenario(heat_solver_rk_or_lbm=False, case=1))
scenarios.add(Scenario(heat_solver_rk_or_lbm=True, case=1))
# scenarios.add(Scenario(heat_solver_rk_or_lbm=False, case=2))
scenarios.add(Scenario(heat_solver_rk_or_lbm=True, case=2))
