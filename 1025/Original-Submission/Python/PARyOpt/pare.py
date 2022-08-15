"""
---
Copyright (c) 2018 Baskar Ganapathysubramanian, Balaji Sesha Sarath Pokuri

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
---

PARE: Parallel Asynchronous Remote Evaluator

Requires python>=3.5, preferably 3.6
Requires PARyOpt>=1.0.1 -- perform a `pip install paryopt ` for ensuring all dependencies

Highly recommended to use venv.

Execute the script as `python pare.py locations_file.txt`

This script reads a file named 'locations_file.txt' which contains the locations to be
evaluated (each location per line), creates respective folders, transfers to the remote
machine, submits a job, waits for completion and finally downloads the
results from the remote machine to the local machine.

"""
import getpass
from typing import List

import numpy as np

from PARyOpt.evaluators import AsyncSbatchEvaluator
import sys
import os
import logging
import time

from PARyOpt.evaluators.async_sbatch import VALUE_FROM_FILE
from PARyOpt.evaluators.connection import Connection, Host


def folder_generator(directory: str, x: np.array) -> None:
    # write config file and put relevant files in the folder "directory".
    # The current location of control is in "directory"
    with open(os.path.join(directory, 'config.txt'), 'w') as f:
        f.write('{}'.format(x[0]*100.))


# for AsyncSbatchEvaluator
def sbatch_job_script_gen(lcl_dir: str, x: np.array):
    """
    Job script that needs to be submitted using command sbatch. The whole script should be returned as a string.
    The user needs to modify this as needed
    :param lcl_dir:
    :param x:
    :return:
    """
    # make changes
    return """#!/bin/bash
#SBATCH --job-name='pare_test'
#SBATCH --output='output.txt'
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:00:01

echo "Done {}"
""".format(' '.join([str(s) for s in x]))


def result_parser(local_dir: str, remote_dir: str, conn: Connection, x: np.array) -> float:
    """
    Parse results from a file downloaded from cluster. The file will be located in remote_dir, which is an extension \
    of the files in local_dir.
    Currently, it will just download the remote folder. The user  has to write the script to parse the results
    and return the appropriate value.
    :param local_dir:
    :param remote_dir:
    :param conn:
    :param x:
    :return:
    """
    # filename = 'myresults.txt'
    # fl = io.BytesIO()
    # conn.sftp().getfo(os.path.join(remote_dir, filename), fl)
    # return float(fl.getvalue().decode(encoding='utf-8').strip())
    conn.get_dir(remote_dir=remote_dir, local_dir=local_dir)
    # with open('result_file.plt', 'r') as rf:
    #     for line in rf:
    #             result_line = abs(float(line.strip()))
    # return result_line
    return 0.0


def results_filename(loc_file: str) -> str:
    """
    create a results file based on the locations file
    :param loc_file: string
    :return: modified string
    """
    # make changes
    res_file = loc_file.split('.')[0] + '_results'
    if '.' in loc_file:
        res_file += '.' + loc_file.split('.')[-1]
    return res_file


def parse_locations_file(loc_file: str, delimiter=None) -> List[np.array]:
    """
    parses the given filename to return a list of evaluation locations
    :param delimiter:
    :param loc_file: filename
    :return: list of arrays
    """
    # make changes
    loc_list = []
    if delimiter is None:
        delimiter = ' '
    with open(loc_file, 'r') as lf:
        for line in lf:
            loc_list.append(np.asarray([float(a) for a in line.split(delimiter)]))

    return loc_list


def write_results(pare_results: list, locations: list, res_file: str) -> None:
    """
    writes the results in the order of the way locations were asked for
    :param pare_results: results from PARE
    :param locations: requested order of results
    :param res_file:
    :return: None
    """
    # make changes
    completed_dict = dict([(tuple(x.tolist()), y) for x,y in pare_results])
    loc_list = [tuple(a.tolist()) for a in locations]
    with open(res_file, 'w') as rf:
        for loc in loc_list:
            rf.write('{}\n'.format(completed_dict[loc]))


if __name__ == "__main__":
    # logging setup
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    async_fraction = 0.9
    jobs_per_batch = 2

    fh = logging.FileHandler('pare_frac_{}_jobs_{}.log'.format(time.strftime("%Y.%m.%d-%H%M%S"), async_fraction,
                                                               jobs_per_batch), mode='a')

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    locations_file = ""
    if len(sys.argv) != 2:
        logger.info('Number of arguments to PARE are incorrect. Please check. Exiting!')
        logger.info('Execute the script as `python3.6 pare.py locations_file.txt`')
        exit(code=192)
    else:
        locations_file = sys.argv[1]
        logger.info('Starting evaluations from file: {}'.format(locations_file))

    next_query_points = parse_locations_file(locations_file)            # type: List[np.array]
    results_file = results_filename(locations_file)                     # type: str

    # parallel, asynchronous, on Condo
    cost_function = AsyncSbatchEvaluator(Host(getpass.getuser(), 'cluster_login.address'),
                                         job_generator=folder_generator, job_script=sbatch_job_script_gen,
                                         lcl_parse_result=result_parser,
                                         remote_jobs_dir='Projects/pare',   # give absolute path on cluster
                                         required_fraction=async_fraction,
                                         max_pending=jobs_per_batch)

    completed, pending_query_points, failed = cost_function(next_query_points, [])

    assert len(failed) == 0, "There are failed query points at {}... Exiting!".format(failed)
    assert len(pending_query_points) == 0,  "There are incomplete query points at {}... Exiting!".format(pending_query_points)

    write_results(completed, next_query_points, results_file)

    logger.info('Results written to {}'.format(results_file))
    logger.removeHandler(fh)

