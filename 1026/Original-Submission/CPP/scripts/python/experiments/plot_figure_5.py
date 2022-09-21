import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.python.experiments.data_reader import read_data, extract_als_data
from scripts.python.utils.figure_config import fig_format, fig_size_in, markersize, linewidth
from scripts.python.utils.system_config import plot_output_path, CPU_FPS, GEMM
from scripts.python.utils.utils import modes_string, modes_title_string, plot_gemm


def plot_x_ranks(ax, als_df):
    ttime, ctime, flops, ranks, ittime = extract_als_data(als_df)
    flop_cumsum_als = np.cumsum(flops)
    ranks = list(ranks)
    x = []
    for i in range(1, 21):
        x.append(flop_cumsum_als[ranks.index(i)])
    for el in x:
        ax.axvline(el, ymax=0.04, c='g', linewidth=linewidth)


def performance_plot_both(dic, ax=None, print_all=False):
    als_df = dic['alsdata']
    als_omp = dic['alsompdata']
    cals_df = dic['calsdata']
    ccals_df = dic['ccalsdata']
    ttb_l = dic['ttbdata']

    if ax is None:
        fig, ax = plt.subplots(1, 1)
        fig.set_size_inches(w=4.68596, h=3.5)

    threads = str(cals_df['NUM_THREADS'][0])

    mfps = CPU_FPS[threads]
    gemm = GEMM[modes_string(dic['modes'])][str(threads)]

    ttime, ctime, flops, ranks, ittime = extract_als_data(als_df)
    flop_cumsum_als = np.cumsum(flops)

    ttime_omp, ctime_omp, flops_omp, ranks_omp, ittime_omp = extract_als_data(als_omp)
    flop_cumsum_als_omp = np.cumsum(flops_omp)

    flop_cumsum_cals = cals_df['FLOPS'].cumsum()

    print()
    print('CALS Flops: {:>14}, Total: {:>8.2f}, Iteration sum: {:>8.2f}'.format(list(flop_cumsum_cals)[-1],
                                                                      cals_df['ITERATION'].sum(),
                                                                      cals_df['TOTAL'].max()))
    print('OALS Flops: {:>14}, Total: {:>8.2f}'.format(list(flop_cumsum_als_omp)[-1],
                                                ttime_omp.max()))
    print(' ALS Flops: {:>14}, Total: {:>8.2f}, Iteration sum: {:>8.2f}'.format(list(flop_cumsum_als)[-1],
                                                                      ittime.sum(),
                                                                      ttime.sum()))
    print()

    ax.step(flop_cumsum_cals, cals_df['FLOPS'] / cals_df['ITERATION'] / mfps,
            '-', label='CALS', color='C0', markersize=markersize, linewidth=linewidth)

    print('{} {} {} {}'.format(flops_omp[-1], ttime_omp.max(), mfps, flops_omp[-1] / ttime_omp.max() / mfps))

    if threads != '1':
        val = flop_cumsum_als_omp[-1] / ttime_omp.max() / mfps
        ax.step([flop_cumsum_als_omp[0], flop_cumsum_als_omp[-1]],
                [val, val],
                '-', label='OMP ALS', color='C6', markersize=markersize, linewidth=linewidth)

    ax.step(flop_cumsum_als, flops / ttime / mfps,
            '-', label='ALS', color='C1', markersize=markersize, linewidth=linewidth)

    if ttb_l:
        ax.step(flop_cumsum_als, flops / np.array(ttb_l) / mfps,
                '-', label='TTB', color='C4', markersize=markersize, linewidth=linewidth)

    plot_gemm(gemm, ax, flop_cumsum_als)

    # Plot the CALS buffer size as xticks
    # xticks = np.arange(1, cals_df['COLS'].count(), step=3)
    # plt.xticks(ticks=xticks, labels=np.array(cals_df['COLS'])[xticks - 1], rotation=45, fontsize=3)

    # Plot the ALS ranks as xticks
    # xticks = np.arange(1, len(ranks), step=1)
    # plt.xticks(ticks=xticks, labels=ranks[xticks - 1], rotation=45, fontsize=3)

    # Plot total distance as xticks
    flop_cumsum_cals = np.array(flop_cumsum_cals)
    ax.set_xticks([0, 0.33 * flop_cumsum_cals[-1], 0.66 * flop_cumsum_cals[-1], flop_cumsum_cals[-1]])

    if threads == '24' or print_all:
        ax.set_xticklabels(['0', '.33', '.66', '1'])
    # xticks = np.arange(1, len(ranks), step=1)
    # plt.xticks(ticks=xticks, labels=ranks[xticks - 1], rotation=45, fontsize=3)

    # if (dic['modes'] == (200, 200, 200)) and (threads == '1'):
    plot_x_ranks(ax, als_df)

    if ((dic['modes'] == (100, 100, 100) or dic['modes'] == (299, 301, 41)) and (threads == '12')) or print_all:
        ax.legend()

    if ((dic['modes'] == (100, 100, 100)) or (dic['modes'] == (100, 100, 100) and threads == '1')) or print_all:
        ax.set_ylabel('Efficiency (Threads: {})'.format(threads))
    else:
        ax.tick_params(labelleft=False, left=True)

    if dic['modes'] == (299, 301, 41):
        ax.set_ylabel('Efficiency (Threads: {})'.format(threads))
        ax.tick_params(labelleft=True, left=True)
    else:
        if threads == '24':
            ax.set_xlabel('Total computation')

    # if threads_on_title:
    #     ax.set_title('Threads: {}'.format(threads))
    # else:
    #     ax.set_title(mode_string_title(dic['modes']))

    if threads == "1" or print_all:
        ax.set_title(modes_title_string(dic['modes']))

    ax.set_xlim([-0.02 * flop_cumsum_cals[-1], flop_cumsum_cals[-1] + 0.02 * flop_cumsum_cals[-1]])

    ax.set_ylim([0, 1])
    ax.set_yticks(ticks=np.arange(0, 1.1, step=0.1))

    ax.grid(True, axis='y')
    # plt.tight_layout()
    if ax is None:
        plt.savefig(plot_output_path
                    + 'ALS_v_CALS_'
                    + dic['backend']
                    + '_modes_' + modes_string(dic['modes'])
                    + '_threads_' + str(dic['threads'])
                    + fig_format)


if __name__ == '__main__':
    backend = 'MKL'

    fig, ax = plt.subplots(3, 3, sharey='all', sharex='col')
    fig.set_size_inches(w=fig_size_in['width'], h=6)
    for im, modes in enumerate([(100, 100, 100), (200, 200, 200), (300, 300, 300)]):
        for it, threads in enumerate([1, 12, 24]):
            results_dict = read_data(backend, threads, modes)
            performance_plot_both(results_dict, ax[it][im])
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.11)
    plt.savefig(plot_output_path
                + 'ALS_v_CALS_'
                + backend
                + '_efficiency'
                + fig_format)
    plt.show()
