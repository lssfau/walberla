from os import getcwd
from waLBerla.tools.jobscripts import createJobscript
from datetime import timedelta


for node_count in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2400]:
    with open("job_weak_scaling_{:04d}.sh".format(node_count), 'w') as f:
        js = createJobscript(nodes=node_count,
                             output_file='out_lbm_bench_%j.txt',
                             error_file='err_lbm_bench_%j.txt',
                             initial_dir=getcwd(),
                             exe_name='UniformGridBenchmarkGPU',
                             parameter_files=['weak_scaling.py'],
                             wall_time=timedelta(minutes=25),
                             machine='pizdaint_hybrid'
                             )
        f.write(js)
