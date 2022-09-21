import matplotlib.pyplot as plt

from scripts.python.mttkrp.cals_mttkrp_analysis import plot_best_mttkrp
from scripts.python.mttkrp.generate_LUTs import best_method_per_rank
from scripts.python.utils.figure_config import fig_format, fig_size_in, f_latex
from scripts.python.utils.system_config import plot_output_path
from scripts.python.utils.utils import modes_string

if __name__ == '__main__':

    backend = 'MKL'
    modes = (300, 300, 300)

    f_xlabel = True
    f_ylabel = True
    f_legend = True

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(w=fig_size_in['width'] - 0.2 * fig_size_in['width'], h=3)

    threads = 1
    best = best_method_per_rank(backend, modes, threads)
    ax = plot_best_mttkrp(backend, modes, threads, best, ax, 'C0',
                          f_ylabel=f_ylabel, f_xlabel=f_xlabel, f_legend=f_legend)
    threads = 12
    best = best_method_per_rank(backend, modes, threads)
    ax = plot_best_mttkrp(backend, modes, threads, best, ax, 'C1',
                          f_ylabel=f_ylabel, f_xlabel=f_xlabel, f_legend=f_legend)
    threads = 24
    best = best_method_per_rank(backend, modes, threads)
    ax = plot_best_mttkrp(backend, modes, threads, best, ax, 'C3',
                          f_ylabel=f_ylabel, f_xlabel=f_xlabel, f_legend=f_legend)

    fig.tight_layout()
    fig.savefig(plot_output_path
                + 'MTTKRP_Best_Benchmark_'
                + backend
                + '_modes_' + modes_string(modes)
                + '_sec_3'
                + fig_format)

    plt.show()
