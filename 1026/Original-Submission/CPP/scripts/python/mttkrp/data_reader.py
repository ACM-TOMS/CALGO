import logging
import os.path

import pandas as pd

from scripts.python.utils.system_config import input_path
from scripts.python.utils.utils import modes_string

logger = logging.getLogger('Reader')


def read_bench(backend, threads, modes, mttkrp_method, gpu=False):
    bench_dict = {}
    folder = input_path + '{}/benchmark/'.format(backend)
    mode_name = modes_string(modes)
    if gpu:
        cuda = "_CUDA"
    else:
        cuda = ""
    bench_file_name = '{}benchmark{}_{}_{}_{}_{}.csv'.format(folder, cuda, backend, mode_name, threads, mttkrp_method)
    if os.path.isfile(bench_file_name):
        logger.info('Found file: {}'.format(bench_file_name))
        bench_dict.update({'backend': backend,
                           'modes': modes,
                           'threads': threads,
                           'mttkrp_method': mttkrp_method,
                           'data': pd.read_csv(bench_file_name, delimiter=';')})
    else:
        logger.warning('File: {} Not Found!'.format(bench_file_name))
    return bench_dict
