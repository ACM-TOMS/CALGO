import logging

from scripts.python.experiments.data_reader import read_data

logger = logging.getLogger('Letter defrag')

logging.basicConfig(level=logging.INFO)


def defrag_plot(backend, threads, modes):
    dic = read_data(backend, threads, modes, 'defrag')

    cals = dic['calsdata']

    cals['P'] = cals['DEFRAGMENTATION'] / cals['ITERATION']
    logger.info("Percent mean: {}".format(cals['P'].mean()))
    logger.info("Percent median: {}".format(cals['P'].median()))
    logger.info("Percent max: {}".format(cals['P'].max()))
    logger.info("Percent min: {}".format(cals['P'].min()))


if __name__ == '__main__':

    backend = 'MKL'

    for im, modes in enumerate([(200, 200, 200)]):
        for i, threads in enumerate([1, 24]):
            defrag_plot(backend, threads, modes)
