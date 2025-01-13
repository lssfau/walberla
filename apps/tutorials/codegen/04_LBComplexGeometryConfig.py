import os
import waLBerla as wlb
from waLBerla.tools.config import block_decomposition
from waLBerla.tools.sqlitedb import sequenceValuesToScalars, checkAndUpdateSchema, storeSingle
import sys

from math import prod

default_boundary_conditions = {'Border': [
                                    {'direction': 'W', 'flag': 'UBB'},
                                    {'direction': 'E', 'flag': 'Outflow'},
                                    {'direction': 'S', 'flag': 'NoSlip'},
                                    {'direction': 'N', 'flag': 'FreeSlip'},
                               ]}

class Scenario:
    def __init__(self, simTime=100, uLB=0.02, velocity_SI=1, viscosity_SI=1.562e-5, Href_SI=0.2,
                 boundary_conditions=None, domainScaling=(3, 1, 1), cells_per_block=(256, 64, 64), 
                 periodic=(0, 0, 1), dx_SI=0.1, meshFile=str("bunny.obj"), cuda_blocks=(64, 1, 1),
                 cuda_enabled_mpi=False, inner_outer_split=(1, 1, 1), vtk_write_frequency=0, 
                 vtk_output_str=str("Fluid"), vtk_flag_output_str=str("Flags"), 
                 remaining_time_logger_frequency=-1):
        
        self.simulationTime = simTime
        self.uLB = uLB

        self.velocity_SI  = velocity_SI
        self.viscosity_SI = viscosity_SI

        self.Href_SI = Href_SI              # Reference height of bunny mm

        self.boundary_conditions = boundary_conditions if boundary_conditions else default_boundary_conditions

        self.domainScaling = domainScaling
        self.cells_per_block = cells_per_block
        self.periodic = periodic
        self.dx_SI = dx_SI

        self.meshFile = meshFile

        self.cuda_blocks = cuda_blocks
        self.cuda_enabled_mpi = cuda_enabled_mpi
        self.inner_outer_split = inner_outer_split

        self.vtk_write_frequency = vtk_write_frequency
        self.vtk_output_str = vtk_output_str
        self.vtk_flag_output_str = vtk_flag_output_str
        self.remaining_time_logger_frequency = remaining_time_logger_frequency

        self.config_dict = self.config(print_dict=False)

    @wlb.member_callback
    def config(self, print_dict=True):
        from pprint import pformat
        config_dict = {
            'DomainSetup': {
                'domainScaling': self.domainScaling,
                'cellsPerBlock': self.cells_per_block,
                'periodic': self.periodic,
                'dx_SI': self.dx_SI
            },
            'Geometry':{
                'meshFile': self.meshFile
            },
            'Parameters': {
                'simTime': self.simulationTime,
                'uLB': self.uLB,
                'velocity_SI': self.velocity_SI,
                'viscosity_SI':self.viscosity_SI,
                'Href_SI': self.Href_SI,
                'cudaEnabledMPI': self.cuda_enabled_mpi,
                'innerOuterSplit': self.inner_outer_split,
                'gpuBlockSize': self.cuda_blocks,
                'vtkWriteFrequency': self.vtk_write_frequency,
                'vtkOutputString': self.vtk_output_str,
                'vtkFlagOutputString':self.vtk_flag_output_str,
                'remainingTimeLoggerFrequency': self.remaining_time_logger_frequency
            },
            'Logging': {
                'logLevel': 'detail',  # info progress detail tracing
            }
        }

        config_dict["Boundaries"] = self.boundary_conditions

        if print_dict:
            wlb.log_info_on_root("Scenario:\n" + pformat(config_dict))

        return config_dict


def run():
    # BC options are: NoSlip, NoSlipQBB, ObjNoSlipQBB, FreeSlip, UBB, Outflow (defined for east boundary only)
    border_bc_dir = {'Border': [    {'direction': 'N', 'flag': 'NoSlipQBB'},
                                    {'direction': 'S', 'flag': 'NoSlipQBB'},
                                    {'direction': 'W', 'flag': 'UBB'},
                                    {'direction': 'E', 'flag': 'Outflow'},
                                ],
                    }

    object_bc_dir = [  {'flag': 'ObjNoSlipQBB'},
                       {'flag': 'NoSlip'}
                    ]
    
    scenarios = wlb.ScenarioManager()

    for object_bc in object_bc_dir:
        boundary_conditions = border_bc_dir.copy()
        boundary_conditions['ObjectBC'] = object_bc

        scenario = Scenario(    simTime=100,
                                uLB=0.08,
                                boundary_conditions=boundary_conditions,
                                cells_per_block=(128, 32, 32),
                                dx_SI=0.03, 
                                meshFile=str("cube.obj"),
                                domainScaling=(6, 2, 2),
                                periodic=(0, 0, 1), 
                                cuda_enabled_mpi=False,
                                vtk_write_frequency=1000,
                                vtk_output_str=f"Fluid_{object_bc['flag']}",
                                vtk_flag_output_str=f"Flags_{object_bc['flag']}",
                                remaining_time_logger_frequency=5
                            )

        scenarios.add(scenario)

#---------------------------------------------------------------------------------------------------------------------
run()
