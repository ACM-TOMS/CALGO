import matplotlib.pyplot as plt

from scripts.python.mttkrp.cals_mttkrp_analysis import plot_best_mttkrp
from scripts.python.mttkrp.generate_LUTs import best_method_per_rank
from scripts.python.utils.figure_config import fig_format, fig_size_in, f_latex
from scripts.python.utils.system_config import plot_output_path

if __name__ == '__main__':

    backend = 'MKL'

    fig, ax = plt.subplots(3, 1, sharex='all')
    fig.set_size_inches(w=fig_size_in['width'] - 0.2 * fig_size_in['width'], h=5)
    f_xlabel = False
    f_ylabel = True
    f_legend = False
    for im, modes in enumerate([(100, 100, 100), (200, 200, 200), (300, 300, 300)]):

        if modes == (300, 300, 300):
            f_xlabel = True
        else:
            f_xlabel = False

        if modes == (200, 200, 200):
            f_legend = True
        else:
            f_legend = False

        threads = 1
        best = best_method_per_rank(backend, modes, threads)
        ax[im] = plot_best_mttkrp(backend, modes, threads, best, ax[im], 'C0',
                                  f_ylabel=f_ylabel, f_xlabel=f_xlabel, f_legend=f_legend)

        threads = 12
        best = best_method_per_rank(backend, modes, threads)
        ax[im] = plot_best_mttkrp(backend, modes, threads, best, ax[im], 'C1',
                                  f_ylabel=f_ylabel, f_xlabel=f_xlabel, f_legend=f_legend)

        threads = 24
        best = best_method_per_rank(backend, modes, threads)
        ax[im] = plot_best_mttkrp(backend, modes, threads, best, ax[im], 'C3',
                                  f_ylabel=f_ylabel, f_xlabel=f_xlabel, f_legend=f_legend)

    fig.tight_layout()
    plt.subplots_adjust(wspace=0.23)
    fig.savefig(plot_output_path
                + 'MTTKRP_Best_Benchmark_'
                + backend
                + fig_format)

    plt.show()
