import numpy as np
import pandas as pd

from scripts.python.experiments.data_reader import extract_als_data


def speedup(dic):
    als = dic['alsdata']
    cals = dic['calsdata']
    ccals = dic['ccalsdata']

    tcals = cals['TOTAL'].max()
    icals = cals['ITERATION'].sum()
    if isinstance(ccals, pd.DataFrame):
        tccals = cals['TOTAL'].max()
        iccals = cals['ITERATION'].sum()

    ttime, ctime, flops, ranks = extract_als_data(als)
    tals = np.sum(ttime)

    print('ALS time: {}'.format(tals))
    print('CALS time: It: {} To:{}'.format(icals, tcals))
    print('ALS v CALS(t) Speedup: ', '{:.2f}'.format(tals / tcals))
    print('ALS v CALS(i) Speedup: ', '{:.2f}'.format(tals / icals))

    if isinstance(ccals, pd.DataFrame):
        print('CUDA time: It: {} To:{}'.format(iccals, tccals))
        print('CALS v CUDA(t) Speedup: ', '{:.2f}'.format(tcals / tccals))
        print('CALS v CUDA(i) Speedup: ', '{:.2f}'.format(icals / iccals))
