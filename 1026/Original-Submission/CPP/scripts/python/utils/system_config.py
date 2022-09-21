import copy
import logging

logger = logging.getLogger('config')

input_path = 'data/'
lut_output_path = input_path
plot_output_path = 'data/figures/'


def generate_flops_per_sec(flops_per_cycle, cpu_ghz):
    flops_per_sec = copy.deepcopy(cpu_ghz)
    for cores, ghz in flops_per_sec.items():
        flops_per_sec[cores] = ghz * 1e9 * flops_per_cycle * int(cores)
    return flops_per_sec


arch = 'skx'

if arch == 'skx':

    # CLAIX18 (Turbo - SKYLAKE (2x) Xeon Platinum 8160)
    CPU_GHZ = {'1': 3.5, '12': 2.6, '24': 2}
    CPU_FPC = 32
    CPU_FPS = generate_flops_per_sec(CPU_FPC, CPU_GHZ)
    GPU_FPS = 7e12  # V100

    # TODO: Have the script auto generate tables with this data instead of manually extracting it.
    # Statistics for GEMMs of 3 matrices, whose elements sum up to the elements of a tensor (absolute best).
    GEMM = {'100-100-100': {'1': {'mean': 0.78456, 'median': 0.81856, 'std': 0.0895013},
                            '12': {'mean': 0.676944, 'median': 0.687896, 'std': 0.0320906},
                            '24': {'mean': 0.554468, 'median': 0.566198, 'std': 0.062533}},
            '200-200-200': {'1': {'mean': 0.857106, 'median': 0.892538, 'std': 0.0810871},
                            '12': {'mean': 0.748761, 'median': 0.748748, 'std': 0.0593422},
                            '24': {'mean': 0.744228, 'median': 0.745608, 'std': 0.0111339}},
            '300-300-300': {'1': {'mean': 0.863431, 'median': 0.894378, 'std': 0.0775792},
                            '12': {'mean': 0.857263, 'median': 0.883854, 'std': 0.0517553},
                            '24': {'mean': 0.833132, 'median': 0.834759, 'std': 0.0134416}},
            '299-301-41': {'1': {'mean': 0.857106, 'median': 0.892538, 'std': 0.0810871},
                           '12': {'mean': 0.748761, 'median': 0.748748, 'std': 0.0593422},
                           '24': {'mean': 0.744228, 'median': 0.745608, 'std': 0.0111339}}}
elif arch == 'haswell':

    # Node 72 (No Turbo - HASWELL (2x) E5-2680 v3)
    CPU_GHZ = {'1': 2.5, '12': 2.5, '24': 2.5}
    CPU_FPC = 16
    CPU_FPS = generate_flops_per_sec(CPU_FPC, CPU_GHZ)

elif arch == 'sandy':

    # Rubik (No Turbo - SANDYBRIDGE (2x) E5-2670 v2)
    CPU_GHZ = {'1': 2.5, '10': 2.5, '20': 2.5}
    CPU_FPC = 8
    CPU_FPS = generate_flops_per_sec(CPU_FPC, CPU_GHZ)

else:
    logger.error('Architecture {} not supported.'.format(arch))
    exit()
