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


# helper local parse result functions
def VALUE_FROM_FILE(filename):
    def value_from_file(local_dir: str, remote_dir: str, conn: Connection, x: np.array) -> float:
        fl = io.BytesIO()
        conn.sftp().getfo(os.path.join(remote_dir, filename), fl)
        return float(fl.getvalue().decode(encoding='utf-8').strip())
    return value_from_file


class AsyncSbatchEvaluator(AsyncFunctionEvaluator):
    """
    Class for cost functions that evaluated by launching a job on a remote machine running the SLURM job scheduler.

    :param host: Host object containing the credentials for the server to connect to
    :param job_generator: callable that sets up the run directory for a given x (by e.g. writing config files). \
    It will be passed two arguments: the job directory and the point to evaluate at (x).
    :param job_script: either a string (for a fixed job script), or a callable that returns the job script string. \
    In the latter case, job_script will be passed two arguments: the job directory and the point to evaluate at (x).
    :param remote_parse_result: callable that returns the cost function evaluated at x. \
    It will be passed two arguments: the job directory and the point to evaluate at (x). \
    This will be called after the command returned by run_cmd_generator has terminated (gracefully or otherwise). \
    If the process did not terminate successfully or the result is otherwise unavailable, parse_result \
    should raise any exception. This will signal the optimization routine to not try this point again. \
    This function will be executed *on the remote host*. This requires the remote host to have a matching version \
    of Python installed and the dill module.
    :param lcl_parse_result: callable that returns the cost function evaluated at X. \
    It is passed three arguments: the local job dir, remote job dir, the Connection object to the remote, and X. \
    It is executed on the local machine. This does not require the remote to have Python installed.
    :param lcl_jobs_dir: optional base directory to generate jobs in - default is $PWD/opt_jobs.
    :param remote_jobs_dir: optional base directory to upload jobs to - default is $HOME/paryopt_jobs.
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
                 lcl_parse_result: Callable[[str, str, Connection, np.array], float] = None,
                 remote_parse_result: Callable[[str, np.array], float] = None,
                 lcl_jobs_dir: str = os.path.join(os.getcwd(), 'opt_jobs'),
                 squeue_update_rate: timedelta = timedelta(seconds=30),
                 remote_jobs_dir: str = 'paryopt_jobs',
                 required_fraction=1.0,
                 max_pending=25):
        super().__init__(required_fraction=required_fraction, max_pending=max_pending)

        if (lcl_parse_result and remote_parse_result) or (not lcl_parse_result and not remote_parse_result):
            raise Exception('You must specify either a local parse result function or remote parse result function.')

        self.host = host
        self.job_generator = job_generator
        self.job_script = job_script
        self.local_parse_result = lcl_parse_result
        self.remote_parse_result = remote_parse_result
        self.lcl_jobs_dir = lcl_jobs_dir
        self.remote_jobs_dir = remote_jobs_dir

        self.connection = Connection()
        self.connection.connect(self.host)

        self._last_squeue = ""
        self._last_squeue_upd_time = datetime.min
        self._squeue_upd_rate = squeue_update_rate

    def _update_squeue(self):
        stdout, stderr = self.connection.exec_command('squeue -u `whoami`')
        self._last_squeue = stdout
        self._last_squeue_upd_time = datetime.now()

    def squeue(self):
        if (datetime.now() - self._last_squeue_upd_time) >= self._squeue_upd_rate:
            self._update_squeue()
        return self._last_squeue

    def start(self, x: np.array) -> (str, str, int, datetime):
        """
        Generate job directory on local machine, fill in data related to the job like directory, job id and submit time
        """
        if not os.path.exists(self.lcl_jobs_dir):
            os.mkdir(self.lcl_jobs_dir)

        prefix = 'job_' + ('_'.join([str(v) for v in x])) + '_'
        directory = tempfile.mkdtemp(prefix=prefix, dir=self.lcl_jobs_dir)

        self.job_generator(directory, x)

        # generate job script in directory/job.sh
        job_script_path = os.path.join(directory, 'job.sh')
        with open(job_script_path, 'w') as f:
            job_script = self.job_script
            if callable(job_script):
                job_script = job_script(directory, x)
            f.write(job_script)

        # copy it to the remote machine
        remote_dir = os.path.join(self.remote_jobs_dir, os.path.basename(os.path.normpath(directory)))
        self.connection.put_dir(directory, remote_dir)

        # run sbatch from inside the directory
        stdout, stderr = self.connection.exec_command('sbatch job.sh', cwd=remote_dir)

        # grab the job ID from the output (assuming it was submitted successfully)
        match = re.match(r'Submitted batch job (\d+)', stdout)
        if not match:
            raise RuntimeError('sbatch in ' + remote_dir + ' failed: ' + stderr)

        job_id = match.group(1)

        return directory, remote_dir, job_id, datetime.now()

    def check_for_result(self, x: np.array, data: (str, str, int, datetime)) -> \
            Union[ValueNotReady, EvaluateAgain, EvaluationFailed, float]:
        """
        checks for result if the jobid is complete and ping time is after update rate

        :param x: location of evaluation
        :param data: data related to the location. Typically this is directory information, job id and submit time
        :return: one of ValueNotReady, EvaluateAgain, EvaluationFailed or float
        """
        lcl_dir = data[0]
        remote_dir = data[1]
        job_id = data[2]
        submit_time = data[3]

        # if squeue hasn't had a chance to update since submission, wait until it does
        # (this effectively enforces a minimum time for jobs)
        if (datetime.now() - submit_time) < self._squeue_upd_rate:
            return ValueNotReady()

        queue = self.squeue()
        status = re.search(r'\s*' + re.escape(job_id) + r'\s+\S+\s+\S+\s+\S+\s+(\S+)\s+', queue)
        if status and status.group(1) != 'CG':
            print('Waiting for ' + lcl_dir + ' to complete (status: ' + status.group(1) + ')')
            return ValueNotReady()

        print('  ' + lcl_dir + ' completed')

        # it's either no longer in the queue or giving a 'complete' status, so it's done
        try:
            if self.local_parse_result:
                return self.local_parse_result(lcl_dir, remote_dir, self.connection, x)
            elif self.remote_parse_result:
                return self.connection.call_on_remote(self.remote_parse_result, remote_dir, x, remote_cwd=remote_dir)
        except Exception as err:
            if __debug__:
                raise
            else:
                return EvaluationFailed(err)
