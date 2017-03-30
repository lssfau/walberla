from __future__ import print_function, absolute_import, division, unicode_literals

import os
import math


def createJobscript( wall_time = None, nodes = None, cores = None, job_class = None,
                      initial_dir = '~', job_name="waLBerla", 
                      exe_name = None, parameter_files = [], commands = [], hyperthreading=2, 
                      bg_connectivity = "Torus", 
                      output_file=None, error_file=None, **kwargs ):
    
    CORES_PER_NODE = 16 * hyperthreading
    
    
    validNodeCountSmaller512 = [32,64,128,256]
    def nodeCountValid( n ):
        return ( n in validNodeCountSmaller512 ) or ( n % 512 == 0 )
    
    def sanitizeNodes( requestedNodes ):
        """Allowed node counts on JUEQUEEN are 32,64,128,256 and multiples of 512"""
        if requestedNodes in validNodeCountSmaller512: return requestedNodes
        if requestedNodes % 512 == 0: return requestedNodes
        
        for limit in validNodeCountSmaller512:
            if requestedNodes < limit: return limit
        
        # round up to next multiple of 512
        return int ( (requestedNodes / 512 + 1) * 512 )
    
    
    if nodes is not None and cores is not None:
        raise ValueError("You can either specify nodes or cores - not both.")
    
    if hyperthreading not in [1,2,4]:
        raise ValueError("JUQUEEN hyperthreading has to be 1,2 or 4 (requested %d)" %(hyperthreading,) )
    
    if nodes is None and cores is None:
        raise ValueError('Specify either cores or nodes.')
    
    if cores > CORES_PER_NODE and cores % CORES_PER_NODE != 0:
        raise ValueError("When using more than one node, the number of cores has to be a multiple of 16")
    
    if nodes is None:
        nodes = sanitizeNodes( int(math.ceil( cores / CORES_PER_NODE )) )
    if cores is None:
        if not nodeCountValid( nodes ):
            raise ValueError("Allowed node counts are 32,64,128,256 and multiples of 512")
        cores = nodes * CORES_PER_NODE

    if not output_file: output_file = job_name
    if not error_file:  error_file  = job_name
    
    
    assert( nodeCountValid(nodes) )
    
    tasks_per_node = min( CORES_PER_NODE, cores )
    
    template_file = os.path.join(  os.path.dirname( os.path.realpath(__file__)  ), "juqueen.job" )
    
    
    result = open(template_file).read().format( cores = cores, 
                                                nodes = nodes, 
                                                initial_dir = initial_dir, 
                                                tasks_per_node = tasks_per_node,
                                                job_name = job_name,
                                                bg_connectivity = bg_connectivity,
                                                output_file = output_file,
                                                error_file = error_file,
                                                wall_time = wall_time )
    
    exec_line = "runjob --np %d --ranks-per-node %d : %s %s\n" 
    
    if exe_name is not None:
        for param_file in parameter_files:
            result += exec_line %( cores, tasks_per_node, exe_name, param_file )
    
    for exe_paramfile_pair in commands:
        if type(exe_paramfile_pair) is not tuple:
            result += exe_paramfile_pair + "\n"
        else:
            result += exec_line %( cores, tasks_per_node, exe_paramfile_pair[0], exe_paramfile_pair[1] )
    
    
    return result

