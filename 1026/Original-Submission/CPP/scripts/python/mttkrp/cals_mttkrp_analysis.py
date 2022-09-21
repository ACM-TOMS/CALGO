import logging

import matplotlib.pyplot as plt
import numpy as np

from scripts.python.mttkrp.data_reader import read_bench
from scripts.python.mttkrp.generate_LUTs import best_method_per_rank
from scripts.python.utils.figure_config import markersize, linewidth, fig_format, fig_size_in
from scripts.python.utils.system_config import CPU_FPS, GPU_FPS, plot_output_path
from scripts.python.utils.utils import modes_string, modes_title_string, MttkrpMethod

logger = logging.getLogger('CALS MTTKRP Analyzer')


def mttkrp_performance_analysis(backend, modes, threads, gpu=False):
    fig, ax = plt.subplots(len(modes), 1, sharex='all')
    fig.set_size_inches(w=fig_size_in['width'], h=5)

    if gpu:
        mfps = GPU_FPS
    else:
        mfps = CPU_FPS[str(threads)]
    mttkrp_configurations = [MttkrpMethod.MTTKRP, MttkrpMethod.TWOSTEP0, MttkrpMethod.TWOSTEP1]

    for midx, mode in enumerate(modes):
        for m in mttkrp_configurations:
            dic = read_bench(backend, threads, modes, m, gpu)
            df = dic['data']
            dm = df.loc[df['MODE'] == midx]
            y = dm['FLOPS'] / dm['TIME'] / mfps
            label, color = get_label_and_color(dic['mttkrp_method'], midx, modes)
            ax[midx].plot(dm['RANK'], y, '-o', label=label, c=color, markersize=markersize, linewidth=linewidth)
            # ax[midx].set_xticks(np.arange(dm['RANK'].min(), dm['RANK'].max(), step=49))
            ax[midx].set_ylim([0, 1])
            ax[midx].set_yticks(np.arange(0, 1.1, step=0.2))
            ax[midx].set_ylabel('Efficiency')
            ax[midx].grid(True)
            ax[midx].set_xscale('log')
            ax[midx].set_title('Mode {}'.format(midx))
            ax[midx].tick_params(labelright=True)
            ax[midx].legend(loc='upper left')
    ax[-1].set_xlabel('Components')
    ax[-1].set_xticks([1, 10, 100, 1000])
    ax[-1].set_xticklabels(['1', '10', '100', '1000'])
    fig.suptitle('CALS Implementations (BLAS: ' + backend
                 + ', Threads: ' + str(threads)
                 + ', Modes: ' + modes_title_string(modes) + ')')

    plt.tight_layout()
    fig.savefig(plot_output_path
                + 'MTTKRP_Perf_' + backend
                + '_threads_' + str(threads)
                + '_modes_' + modes_string(modes)
                + fig_format)
    return fig


def get_label_and_color(method, mode_index, modes):
    simple = True

    I = modes[0]
    J = modes[1]
    K = modes[2]
    label, color = '', ''
    if method == MttkrpMethod.MTTKRP:
        if simple:
            label = 'KRP + MTTKRP'
        else:
            if mode_index == 0:
                label = '({} x {}) * ({} x R) - {} block'.format(I, J * K, J * K, 1)
            elif mode_index == 1:
                label = '({} x {}) * ({} x R) - {} blocks'.format(J, I, I, K)
            elif mode_index == 2:
                label = '({} x {}) * ({} x R) - {} block'.format(K, I * J, I * J, 1)
            else:
                logger.error('midx exceeds 3D size ({}). This analyzer is designed for 3D tensors.'.format(mode_index))
                exit()
        color = 'r'
    elif method == MttkrpMethod.TWOSTEP0:
        if simple:
            label = 'TTM + TTV'
        else:
            if mode_index == 0:
                label = '({} x {}) * ({} x R) - {} block'.format(I * J, K, K, 1)
            elif mode_index == 1:
                label = '({} x {}) * ({} x R) - {} block'.format(I * J, K, K, 1)
            elif mode_index == 2:
                label = '({} x {}) * ({} x R) - {} block'.format(I * K, J, J, 1)
            else:
                logger.error('midx exceeds 3D size ({}). This analyzer is designed for 3D tensors.'.format(mode_index))
                exit()
        color = 'b'
    elif method == MttkrpMethod.TWOSTEP1:
        if simple:
            label = 'TTM + TTV'
        else:
            if mode_index == 0:
                label = '({} x {}) * ({} x R) - {} blocks'.format(I * K, J, J, K)
            elif mode_index == 1:
                label = '({} x {}) * ({} x R) - {} block'.format(J * K, I, I, 1)
            elif mode_index == 2:
                label = '({} x {}) * ({} x R) - {} blocks'.format(J * K, I, I, K)
            else:
                logger.error('midx exceeds 3D size ({}). This analyzer is designed for 3D tensors.'.format(mode_index))
                exit()
        color = 'g'
    else:
        logger.error('MTTKRP method not found: {}'.format(method))
        exit()
    return label, color


def plot_best_mttkrp(backend, modes, threads, best, ax=None, c=None, f_ylabel=True, f_xlabel=True, f_legend=True, gpu=False):

    mfps = CPU_FPS[str(threads)]
    if gpu:
        mfps = 7e12

    # Accumulate the FLOPS and the TIME per mode, in order to calculate the efficiency of MTTKRP for all modes
    x = np.array(0)
    y_time = np.array(0)
    y_flops = np.array(0)
    for idx, data in enumerate(best):
        if idx == 0:
            x = np.array([rank for rank, flops, time, method in data])
            y_time = np.array([time for rank, flops, time, method in data])
            y_flops = np.array([flops for rank, flops, time, method in data])
        else:
            y_time += np.array([time for rank, flops, time, method in data])
            y_flops += np.array([flops for rank, flops, time, method in data])
    y = y_flops / y_time / mfps

    # Remove rank values from 11 to 19 (to conform to log scale)
    x, y = list(x), list(y)
    for rank in range(11, 20):
        if rank in x:
            ind = x.index(rank)
            x.pop(ind)
            y.pop(ind)
    x, y = np.array(x), np.array(y)

    if threads == 1:
        label = '1 thread'
    else:
        label = '{} threads'.format(threads)

    fig = None
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    ax.plot(x, y, '-o', color=c, label=label, markersize=markersize, linewidth=linewidth)

    ###################
    # Figure formatting
    ###################

    # y axis
    ax.set_ylim([0, 1])
    ax.set_yticks(np.arange(0, 1.1, step=0.2))
    if f_ylabel:
        ax.set_ylabel('Efficiency')

    # x axis
    # if modes == (300, 300, 300) or modes == (299, 301, 41):
    ax.set_xscale('log')
    ax.set_xticks([1, 10, 100, 1000])

    if f_xlabel:
        ax.set_xticklabels([str(i) for i in [1, 10, 100, 1000]])
        ax.set_xlabel('Components')

    # rest
    ax.tick_params(labelleft=True, left=True)
    ax.grid(b=True, axis='y')

    if f_legend:
        ax.legend(loc='upper left')

    ax.set_title(modes_title_string(modes))
    if fig:
        fig.savefig(plot_output_path
                    + 'MTTKRP_Best_Benchmark_'
                    + backend
                    + '_threads_' + str(threads)
                    + '_modes_' + modes_string(modes)
                    + fig_format)
    return ax


if __name__ == '__main__':

    backend = 'MKL'
    gpu = False
    modes_l = []
    # modes_l.append((117, 18, 702))
    # modes_l.append((299, 301, 41))
    # modes_l.append((405, 136, 19))
    # modes_l.append((100, 100, 100))
    modes_l.append((200, 200, 200))
    # modes_l.append((300, 300, 300))
    for modes in modes_l:
        for threads in [1, 24]:
            mttkrp_performance_analysis(backend, modes, threads, gpu)
            best = best_method_per_rank(backend, modes, threads, gpu)
            plot_best_mttkrp(backend, modes, threads, best, gpu=gpu)
    plt.show()
