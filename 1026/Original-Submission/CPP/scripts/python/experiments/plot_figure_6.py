import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.python.experiments.data_reader import extract_als_data, read_data
from scripts.python.utils.figure_config import fig_size_in, fig_format
from scripts.python.utils.system_config import plot_output_path
from scripts.python.utils.utils import modes_string, modes_title_string

if __name__ == '__main__':

    modes_v = []
    modes_v.append((299, 301, 41))
    modes_v.append((405, 136, 19))
    modes_v.append((255, 281, 25))

    for modes in modes_v:
        backend = 'MKL'

        fig, ax = plt.subplots(1, 1)
        fig.set_size_inches(w=fig_size_in['width'], h=3)

        index = [1, 24]
        columns = {'TTB': [], 'CP-ALS': [], 'OMP ALS': [], 'CALS': []}
        df = pd.DataFrame(index=index, columns=columns)
        for th in index:
            dic = read_data(backend, th, modes)

            als = dic['alsdata']
            als_cuda = dic['alscudadata']
            als_omp = dic['alsompdata']
            als_omp_cuda = dic['alsompcudadata']
            cals = dic['calsdata']
            cals_cuda = dic['calscudadata']
            ttb = dic['ttbdata']

            df.at[th, 'TTB'] = np.sum(ttb)

            ttime, _, _, _, _ = extract_als_data(als)
            df.at[th, 'CP-ALS'] = np.sum(ttime)

            tcals = cals['TOTAL'].max()
            df.at[th, 'CALS'] = tcals


            if th == 24:
                ttime_omp, _, _, _, _ = extract_als_data(als_omp)
                df.at[th, 'OMP ALS'] = np.max(ttime_omp)

            if th == 24 and isinstance(als_cuda, pd.DataFrame):
                ttime_cu, _, _, _, _ = extract_als_data(als_cuda)
                df.at['CUDA', 'CP-ALS'] = np.sum(ttime_cu)

            if th == 24 and isinstance(als_omp_cuda, pd.DataFrame):
                ttime_omp_cu, _, _, _, _ = extract_als_data(als_omp_cuda)
                df.at['CUDA', 'OMP ALS'] = np.max(ttime_omp_cu)

            if th == 24 and isinstance(cals_cuda, pd.DataFrame):
                print('CUDA Iteration - CUDA Total: {:0.3f}'.format(cals_cuda['ITERATION'].sum() - cals_cuda['TOTAL'].max()))
                df.at['CUDA', 'CALS'] = cals_cuda['ITERATION'].sum()

        index = ['1 thread', '24 threads', 'CUDA']
        df.index = index

        print(df.to_latex(float_format="{:0.2f}".format, na_rep='-'))
        df.plot.bar(ax=ax, color=['C4', 'C1', 'C6', 'C0', 'C3'], rot=0)
        ax.set_title(modes_title_string(modes))
        ax.set_ylabel('Time in seconds')
        # ax.set_yscale('log')
        # yticks = [1, 10, 100, 1000, 10000]
        # ax.set_yticks(yticks)
        # ax.set_yticklabels([str(i) for i in yticks])
        # ax.set_ylim([0.9, 110])

        old_lim = list(ax.get_ylim())
        old_lim[1] += 0.05 * old_lim[1]
        ax.set_ylim(old_lim)
        for p in ax.patches:
            height = int(round(p.get_height()))
            if not p.get_height() == 0:
                ax.annotate(str(height),
                            xy=(p.get_x() + p.get_width() / 2, height),
                            xytext=(0, 1),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom')
        plt.tight_layout()
        plt.savefig(plot_output_path
                    + 'ALS_v_CALS_'
                    + backend
                    + '_real_'
                    + modes_string(modes)
                    + fig_format)
    plt.show()
