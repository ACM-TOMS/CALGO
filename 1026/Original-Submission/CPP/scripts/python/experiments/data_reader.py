import logging
import os.path

import numpy as np
import pandas as pd

from scripts.python.utils.system_config import input_path
from scripts.python.utils.utils import modes_string

logger = logging.getLogger('Data reader')


def read_cuda_csv(file_name):
    if os.path.isfile(file_name):
        return pd.read_csv(file_name, delimiter=';')
    else:
        return []


def read_ttb(file_name):
    if os.path.isfile(file_name):
        x = []
        file = open(file_name, 'r')
        for y in file.readlines():
            if y:
                x.append(float(y))
        return x
    else:
        return []


def read_data(backend, threads, modes, suffix=None):
    results_dict = {}
    folder = input_path + '{}/'.format(backend)
    mode_name = modes_string(modes)
    if suffix is None:
        suffix = ''
    else:
        suffix = '_{}'.format(suffix)

    als_file_name = '{}ALS_{}_{}_{}{}.csv'.format(folder, backend, mode_name, threads, suffix)
    als_omp_file_name = '{}ALS_OMP_{}_{}_{}{}.csv'.format(folder, backend, mode_name, threads, suffix)
    als_cuda_file_name = '{}ALS_CUDA_{}_{}_{}{}.csv'.format(folder, backend, mode_name, threads, suffix)
    als_omp_cuda_file_name = '{}ALS_OMP_CUDA_{}_{}_{}{}.csv'.format(folder, backend, mode_name, threads, suffix)
    cals_file_name = '{}CALS_{}_{}_{}{}.csv'.format(folder, backend, mode_name, threads, suffix)
    cals_cuda_file_name = '{}CALS_CUDA_{}_{}_{}{}.csv'.format(folder, backend, mode_name, threads, suffix)
    ttb_file_name = '{}TTB_{}_{}.csv'.format(folder, modes_string(modes), threads)

    if not os.path.isfile(als_file_name):
        logger.warning('File: {} not found.'.format(als_file_name))
        alsdata = []
    else:
        alsdata = pd.read_csv(als_file_name, delimiter=';')
    if not os.path.isfile(als_omp_file_name):
        logger.warning('File: {} not found.'.format(als_omp_file_name))
        alsompdata = []
    else:
        alsompdata = pd.read_csv(als_omp_file_name, delimiter=';')
    if not os.path.isfile(als_cuda_file_name):
        logger.warning('File: {} not found.'.format(als_cuda_file_name))
        alscudadata = []
    else:
        alscudadata = pd.read_csv(als_cuda_file_name, delimiter=';')
    if not os.path.isfile(als_omp_cuda_file_name):
        logger.warning('File: {} not found.'.format(als_omp_cuda_file_name))
        alsompcudadata = []
    else:
        alsompcudadata = pd.read_csv(als_omp_cuda_file_name, delimiter=';')
    if not os.path.isfile(cals_file_name):
        logger.warning('File: {} not found.'.format(cals_file_name))
        calsdata = []
    else:
        calsdata = pd.read_csv(cals_file_name, delimiter=';')
    if not os.path.isfile(cals_cuda_file_name):
        logger.warning('File: {} not found.'.format(cals_cuda_file_name))
        calscudadata = []
    else:
        calscudadata = read_cuda_csv(cals_cuda_file_name)
    if not os.path.isfile(ttb_file_name):
        logger.warning('File: {} not found.'.format(ttb_file_name))
        ttbdata = []
    else:
        ttbdata = read_ttb(ttb_file_name)

    results_dict.update({'backend': backend,
                         'threads': threads,
                         'modes': modes,
                         'alsdata': alsdata,
                         'alsompdata': alsompdata,
                         'alscudadata': alscudadata,
                         'alsompcudadata': alsompcudadata,
                         'calsdata': calsdata,
                         'calscudadata': calscudadata,
                         'ttbdata': ttbdata})
    return results_dict


def extract_als_data(df):
    total_time = []
    comp_time = []
    flops = []
    ranks = []
    it_time = []

    for _, row in df.iterrows():
        total_time.append(row['TOTAL'])
        flops.append(row['FLOPS'] * row['ITER'])
        ranks.append(row['KTENSOR_RANK'])
        it_time.append(row['ITERATION'] * row['ITER'])

        comp_sum = 0.0
        for m in range(len(row['TENSOR_MODES'].split('-'))):
            comp_sum += row['MODE_{}_TOTAL_MTTKRP'.format(m)] + row['MODE_{}_UPDATE'.format(m)]
        comp_sum += row['ERROR']
        comp_sum *= row['ITER']
        comp_time.append(comp_sum)

    time_diff = (np.array(total_time) - np.array(comp_time)) / total_time
    # print('------- Begin ALS reading ------------------------');
    # print('ALS comp VS total min: {0:.0%}'.format(np.min(time_diff)))
    # print('ALS comp VS total mean: {0:.0%}'.format(np.mean(time_diff)))
    # print('ALS comp VS total max: {0:.0%}'.format(np.max(time_diff)))
    # print('ALS comp VS total std: {0:.0%}'.format(np.std(time_diff)))
    # print('------- End ALS reading ------------------------');
    return np.array(total_time), np.array(comp_time), np.array(flops), np.array(ranks), np.array(it_time)
