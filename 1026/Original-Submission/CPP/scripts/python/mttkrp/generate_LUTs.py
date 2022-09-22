import logging

logging.basicConfig(level=logging.INFO)
from pathlib import Path

from scripts.python.mttkrp.data_reader import read_bench
from scripts.python.utils.system_config import lut_output_path
from scripts.python.utils.utils import MttkrpMethod, modes_string

logger = logging.getLogger('LUT generator')


def best_method_per_rank(backend, modes, threads, gpu=False, export=False):
    mttkrp_methods = [MttkrpMethod.MTTKRP, MttkrpMethod.TWOSTEP0, MttkrpMethod.TWOSTEP1]
    best = []

    for method in mttkrp_methods:
        bench_dict = read_bench(backend, threads, modes, method, gpu)
        get_best_method(best, bench_dict)

    if export:
        export_to_file(backend, modes, threads, best, gpu)

    return best


def get_best_method(best, bench_dict):
    method = bench_dict['mttkrp_method']
    for midx, mode in enumerate(bench_dict['modes']):

        df = bench_dict['data']
        dm = df.loc[df['MODE'] == midx]
        if len(best) < midx + 1:
            best.append([(e['RANK'], e['FLOPS'], e['TIME'], int(method)) for i, e in dm.iterrows()])
        else:
            for i, e in dm.reset_index().iterrows():
                (rank, best_flops, best_time, best_method) = best[midx][i]
                if not rank == e['RANK']:
                    logger.error('Rank misalignment. Benchmarks were not performed for the same ranks.')
                    exit()
                if best_flops / best_time < e['FLOPS'] / e['TIME']:
                    best[midx][i] = (e['RANK'], e['FLOPS'], e['TIME'], int(method))


def export_to_file(backend, modes, threads, best, gpu):
    if gpu:
        path_to_files = lut_output_path + "GPU" + '/lookup_tables/' + modes_string(modes)
    else:
        path_to_files = lut_output_path + backend + '/lookup_tables/' + modes_string(modes) + '/' + str(threads)
    logger.info('Path to files: {}'.format(path_to_files))

    Path(path_to_files).mkdir(parents=True, exist_ok=True)

    for idx, mode in enumerate(modes):
        f = open(path_to_files + '/' + str(idx), "w")
        for (rank, flops, time, method) in best[idx]:
            f.write('{} {}\n'.format(rank, method))


if __name__ == '__main__':

    backend = 'MKL'
    gpu = True
    modes_l = []
    # modes_l.append((117, 18, 702))
    # modes_l.append((299, 301, 41))
    # modes_l.append((405, 136, 19))
    modes_l.append((255, 281, 25))
    # modes_l.append((100, 100, 100))
    # modes_l.append((200, 200, 200))
    # modes_l.append((300, 300, 300))
    for modes in modes_l:
        for threads in [1, 24]:
            if not gpu:
                best_method_per_rank(backend, modes, threads, gpu=gpu, export=True)
            elif threads == 24:  # make sure that if gpu is enabled, only the 24 thread case is used.
                best_method_per_rank(backend, modes, threads, gpu=gpu, export=True)
