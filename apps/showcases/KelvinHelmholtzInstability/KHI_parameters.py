import os
import waLBerla as wlb

class Scenario:
    def __init__( self, reNumber:float, machNumber:float=0.1,  cellsPerLength:int=64, simulationPeriods:float=5.0, 
                  viscosity_SI:float=1e-5, Href_SI:float=1.0, momentumThickness_SI:float=3e-4, gpu_blocks=(128, 1, 1), 
                  gpu_enabled_mpi:bool=False, usedPitchedMem:bool=False, vtk_periods:float=0.1, vtk_start:int=0, plot_periods:float=0.1, 
                  plot_start:int=0, stabilityCheckFrequency:int = 10000, additional_info=None ):

        self.reNumber = reNumber
        self.machNumber = machNumber

        self.cellsPerLength = cellsPerLength

        self.simulationPeriods = simulationPeriods

        self.viscosity_SI = viscosity_SI
        self.Href_SI = Href_SI
        self.momentumThickness_SI = momentumThickness_SI

        self.gpu_blocks = gpu_blocks
        self.gpu_enabled_mpi = gpu_enabled_mpi
        self.usedPitchedMem = usedPitchedMem

        self.vtk_periods = vtk_periods
        self.vtk_start   = vtk_start

        self.plot_periods = plot_periods
        self.plot_start   = plot_start
        
        self.checkFrequency = stabilityCheckFrequency

        self.config_dict = self.config(print_dict=False)
        self.additional_info = additional_info

    @wlb.member_callback
    def config(self, print_dict=True):
        from pprint import pformat
        config_dict = {
            'Parameters': {
                'reNumber': self.reNumber,
                'machNumber': self.machNumber,
                'cellsPerLength': self.cellsPerLength,
                'simulationPeriods': self.simulationPeriods,
                'viscosity_SI': self.viscosity_SI,
                'Href_SI': self.Href_SI,
                'momentumThickness_SI': self.momentumThickness_SI,
                'gpuBlockSize': self.gpu_blocks,
                'gpuEnabledMPI':self.gpu_enabled_mpi,
                'usedPitchedMem':self.usedPitchedMem,
                'vtk_periods': self.vtk_periods,
                'vtk_start': self.vtk_start,
                'plot_periods': self.plot_periods,
                'plot_start': self.plot_start
            },
            'StabilityChecker': {
                'checkFrequency': self.checkFrequency,
                'streamOutput':   False,
                'vtkOutput':      True
            },
            'Logging': {
                'logLevel': 'info',  # info progress detail tracing
            }
        }

        if print_dict:
            wlb.log_info_on_root("Scenario:\n" + pformat(config_dict))
            if self.additional_info:
                wlb.log_info_on_root("Additional Info:\n" + pformat(self.additional_info))

        return config_dict

def main():
    wlb.log_info_on_root("Method Comparison")
    wlb.log_info_on_root("")
    
    scenarios = wlb.ScenarioManager()

    cases = [   ([1600], [8, 16, 32, 64, 128, 256]),
                ([1e4],  [8, 16, 32, 64, 128, 256, 512]),
                ([1e6],  [8, 16, 32, 64, 128, 256, 512]) 
            ]

    for reNumbers, resolutions in cases:
        for re in reNumbers:
            for nCells in resolutions:
                scenario = Scenario(reNumber=re, simulationPeriods=3, momentumThickness_SI=3e-6, cellsPerLength=nCells, additional_info=None)
                scenarios.add(scenario)

main()



