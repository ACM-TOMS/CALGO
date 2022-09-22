import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.python.experiments.data_reader import read_data, extract_als_data
from scripts.python.utils.figure_config import markersize, linewidth, colors, fig_size_in, fig_format
from scripts.python.utils.system_config import plot_output_path
from scripts.python.utils.utils import modes_title_string

logger = logging.getLogger('Figure 3')

logging.basicConfig(level=logging.INFO)

def speedup_plot(backend, threads, modes, ax, count):
    x = np.arange(1, 21, 1)
    y = []
    yc = []
    for r in range(1, 21):
        dic = read_data(backend, threads, modes, 'speedup_{}'.format(r))

        als = dic['alsdata']
        cals = dic['calsdata']
        cals_cuda = dic['calscudadata']

        tcals = cals['ITERATION'].sum()
        logger.info("Total CALS time: {}".format(tcals))

        ttime, _, _, _, _ = extract_als_data(als)
        tals = np.sum(ttime)
        logger.info("Total ALS time: {}".format(tals))

        y.append(tals / tcals)

        if threads == 24 and isinstance(cals_cuda, pd.DataFrame):
            tccals = cals_cuda['TOTAL'].max()
            yc.append(tals / tccals)

    if threads == 1:
        label = '1 thread'
    else:
        label = '{} threads'.format(threads)

    yticks = [1]
    if not yc:
        yticks.extend(list(np.arange(3, np.max(np.array(y)) + 3, 2)))
    else:
        yticks.extend(list(np.arange(3, np.max(np.array(yc)) + 3, 2)))
    yticks = np.array(yticks)

    if not yc:
        max_y = np.max(np.array(ax.get_yticks()))
        if np.max(np.array(yticks)) > max_y:
            ax.set_yticks(yticks)
            ax.set_yticklabels([str(i) for i in list(yticks)])
            ax.set_ylim([0, np.max(np.array(y)) + 0.1 * np.max(np.array(y))])
            ax.set_ylabel('Speedup')
    else:
        ax.set_yscale('log')
        yticks = [1, 10, 100]
        ax.set_yticks(yticks)
        ax.set_yticklabels([str(i) for i in yticks])
        ax.set_ylim([0.9, 110])
        ax.set_ylabel('Speedup')

    if modes == (300, 300, 300):
        xticks = [1, 5, 10, 15, 20]
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(i) for i in list(xticks)])
        ax.set_xlabel('Components')

    ax.set_title(modes_title_string(dic['modes']))

    ax.grid(b=True, which='both', axis='y')
    ax.plot(x, y, '-o', color=colors[count], label=label, markersize=markersize, linewidth=linewidth)
    if yc:
        ax.plot(x, yc, '-o', color='C2', label='CUDA', markersize=markersize, linewidth=linewidth)

    if modes == (100, 100, 100):
        ax.legend(ncol=2)

    return ax


if __name__ == '__main__':

    backend = 'MKL'

    fig, ax = plt.subplots(3, 1, sharex='all')
    fig.set_size_inches(w=fig_size_in['width'], h=5)
    for im, modes in enumerate([(100, 100, 100), (200, 200, 200), (300, 300, 300)]):
        for i, threads in enumerate([1, 24]):
            speedup_plot(backend, threads, modes, ax[im], i)
        ax[im].axhline(1, ls='-', color='red', markersize=markersize, linewidth=linewidth)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.23)
    fig.savefig(plot_output_path
                + 'Speedup_'
                + backend
                + fig_format)
    plt.show()
