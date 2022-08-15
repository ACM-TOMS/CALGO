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
"""

## --- end license text --- ##
import subprocess

"""
SLURM scheduler asynchronous evaluator sub-class
"""
from datetime import datetime, timedelta
import re
import tempfile
import numpy as np
import os
import io
from .connection import Connection, Host
from typing import Callable, Union

from .paryopt_async import AsyncFunctionEvaluator, EvaluationFailed, ValueNotReady, EvaluateAgain


class AsyncLocalSbatchEvaluator(AsyncFunctionEvaluator):
    """
    Class for cost functions that evaluated by launching a job on the cluster running the SLURM job scheduler. The
    optimization routine runs on a single core on the HPC while it launches multiple jobs on the same HPC.

    :param job_generator: callable that sets up the run directory for a given x (by e.g. writing config files). \
    It will be passed two arguments: the job directory and the point to evaluate at (x).
    :param job_script: either a string (for a fixed job script), or a callable that returns the job script string. \
    In the latter case, job_script will be passed two arguments: the job directory and the point to evaluate at (x).
    :param parse_result: callable that returns the cost function evaluated at X. \
    It is passed two arguments: the job dir, remote job dir, and X. \
    It is executed on the local machine. This requires the sbatch machine to have Python installed.
    :param jobs_dir: optional base directory to generate jobs in - default is $PWD/opt_jobs.
    :param squeue_update_rate: minimum time between squeue calls. Lower for better job latency, higher to be more \
    polite
    :param required_fraction: fraction of points which must complete before continuing to the next iteration see \
    AsyncEvaluator for more info and implementation
    :param max_pending: maximum simultaneous queued jobs, defaults to 25, see AsyncEvaluator for implementation \
    """

    def __init__(self,
                 host: Host,
                 job_generator: Callable[[str, np.array], None],
                 job_script: Union[str, Callable[[str, np.array], str]],
                 parse_result: Callable[[str, np.array], float] = None,
                 jobs_dir: str = os.path.join(os.getcwd(), 'opt_jobs'),
                 squeue_update_rate: timedelta = timedelta(seconds=30),
                 required_fraction=1.0,
                 max_pending=25):
        super().__init__(required_fraction=required_fraction, max_pending=max_pending)

        self.host = host
        self.job_generator = job_generator
        self.job_script = job_script
        self.parse_result = parse_result
        self.jobs_dir = jobs_dir

        self.connection = Connection()
        self.connection.connect(self.host)

        self._last_squeue = ""
        self._last_squeue_upd_time = datetime.min
        self._squeue_upd_rate = squeue_update_rate
        self._total_jobs_submitted = 0

    def _update_squeue(self):
        proc = subprocess.Popen(['squeue', '-u', '`whoami'], stdout=subprocess.STDOUT,
                                stdin=subprocess.DEVNULL, start_new_session=True)
        self._last_squeue = proc.stdout
        self._last_squeue_upd_time = datetime.now()

    def squeue(self):
        if (datetime.now() - self._last_squeue_upd_time) >= self._squeue_upd_rate:
            self._update_squeue()
        return self._last_squeue

    def start(self, x: np.array) -> (str, int, datetime):
        """
        Generate job directory and fill in data related to the job like directory, job id and submit time
        """
        if not os.path.exists(self.jobs_dir):
            os.mkdir(self.jobs_dir)

        directory = tempfile.mkdtemp(prefix='job_' + str(self._total_jobs_submitted) + '_', dir=self.jobs_dir)
        os.mkdir(directory)

        self.job_generator(directory, x)

        # generate job script in directory/job.sh
        job_script_path = os.path.join(directory, 'job.sh')
        with open(job_script_path, 'w') as f:
            job_script = self.job_script
            if callable(job_script):
                job_script = job_script(directory, x)
            f.write(job_script)

        # run sbatch from inside the directory:
        run_cmd = ['sbatch', 'job.sh']
        proc = subprocess.Popen(run_cmd, stdout=open(os.path.join(directory, 'OUT.LOG'), 'wb'),
                                stdin=subprocess.DEVNULL,
                                shell=(run_cmd is str), cwd=directory, start_new_session=True)

        self._total_jobs_submitted += 1
        # grab the job ID from the output (assuming it was submitted successfully)
        match = re.match(r'Submitted batch job (\d+)', proc.stdout)
        if not match:
            raise RuntimeError('sbatch in ' + directory + ' failed: ' + proc.stderr)

        job_id = match.group(1)
        print('Job ID: {} for {}, in directory:{}'.format(job_id, x, directory))

        return directory, job_id, datetime.now()

    def check_for_result(self, x: np.array, data: (str, int, datetime)) -> \
            Union[ValueNotReady, EvaluateAgain, EvaluationFailed, float]:
        """
        checks for result if the jobid is complete and ping time is after update rate

        :param x: location of evaluation
        :param data: data related to the location. Typically this is directory information, job id and submit time
        :return: one of ValueNotReady, EvaluateAgain, EvaluationFailed or float
        """
        directory = data[0]
        job_id = data[1]
        submit_time = data[2]

        # if squeue hasn't had a chance to update since submission, wait until it does
        # (this effectively enforces a minimum time for jobs)
        if (datetime.now() - submit_time) < self._squeue_upd_rate:
            return ValueNotReady()

        queue = self.squeue()
        status = re.search(r'\s*' + re.escape(job_id) + r'\s+\S+\s+\S+\s+\S+\s+(\S+)\s+', queue)
        if status and status.group(1) != 'CG':
            print('Waiting for ' + directory + ' to complete (status: ' + status.group(1) + ')')
            return ValueNotReady()

        print('  ' + directory + ' completed')

        # it's either no longer in the queue or giving a 'complete' status, so it's done
        try:
            return self.parse_result(directory, x)
        except Exception as err:
            if __debug__:
                raise
            else:
                return EvaluationFailed(err)
