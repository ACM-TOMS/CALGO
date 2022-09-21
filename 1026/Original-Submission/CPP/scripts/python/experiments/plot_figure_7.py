import matplotlib.pyplot as plt
import pandas as pd

from scripts.python.experiments.data_reader import read_data
from scripts.python.utils.figure_config import fig_size_in, fig_format
from scripts.python.utils.system_config import plot_output_path
from scripts.python.utils.utils import modes_title_string

if __name__ == '__main__':
    backend = 'MKL'

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(w=fig_size_in['width'], h=3)
    modes_list = [(100, 100, 100), (200, 200, 200), (300, 300, 300)]
    threads = [24]
    # results = {'ALS': [], 'OMP ALS': [], 'CALS': [], 'ALSC': [], 'OMP ALSC': [], 'CALSC': []}
    results = {'ALS CUDA': [], 'OMP ALS CUDA': [], 'CALS CUDA': []}
    for modes in modes_list:
        for th in threads:
            dic = read_data(backend, th, modes)
            als = dic['alsdata']
            alsomp = dic['alsompdata']
            alscuda = dic['alscudadata']
            alsompcuda = dic['alsompcudadata']
            cals = dic['calsdata']
            ttb = dic['ttbdata']
            cuda = dic['calscudadata']

            # results['ALS'].extend([als['TOTAL'].sum()])
            # results['OMP ALS'].extend([alsomp['TOTAL'].max()])
            # results['CALS'].extend([cals['ITERATION'].sum()])

            results['ALS CUDA'].extend([alscuda['TOTAL'].sum()])
            results['OMP ALS CUDA'].extend([alsompcuda['TOTAL'].max()])
            if th == 24 and isinstance(cuda, pd.DataFrame):
                print('CUDA Iteration - CUDA Total: {:0.3f}'.format(cuda['ITERATION'].sum() - cuda['TOTAL'].max()))
                results['CALS CUDA'].extend([cuda['TOTAL'].max()])

    index = [modes_title_string(i) for i in modes_list]
    df = pd.DataFrame(results, index=index)
    # df.plot.bar(ax=ax, color=['C1', 'C6', 'C0', '#5fd35f', '#2ca02c', '#165016'], rot=0)
    df.plot.bar(ax=ax, color=['#5fd35f', '#2ca02c', '#165016'], rot=0)
    ax.set_ylabel('Time in seconds')
    old_lim = list(ax.get_ylim())
    old_lim[1] += 0.05 * old_lim[1]
    ax.set_ylim(old_lim)
    for p in ax.patches:
        height = round(p.get_height(), 1)
        if not p.get_height() == 0:
            ax.annotate(str(height),
                        xy=(p.get_x() + p.get_width() / 2, height),
                        xytext=(0, 1),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(plot_output_path
                + 'CUDA_v_CALS_'
                + backend
                + fig_format)
    plt.show()
