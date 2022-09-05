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
Local asynchronous evaluator sub-class
"""
from typing import Union, Callable, List, Any
import signal
import tempfile
import os
import errno
import subprocess
import numpy as np
from multiprocessing import cpu_count

from .paryopt_async import AsyncFunctionEvaluator, ValueNotReady, EvaluationFailed, EvaluateAgain


def _pid_is_running(pid: int) -> bool:
    """
    Adapted from this post: https://stackoverflow.com/a/6940314

    :param pid: process ID to check for
    :return: True if a process with pid is running
    """
    if pid <= 0:
        raise ValueError('Invalid PID: ' + str(pid))

    try:
        # Try to send "signal 0" to pid.
        # Signal 0 is not a real signal. It is a special value that doesn't actually send a signal, but still performs
        # error checking (i.e. kill will only return successfully if we are able to send signals to pid, which means
        # it must be alive).
        os.kill(pid, 0)
    except OSError as err:
        if err.errno == errno.ESRCH:  # no such process
            return False
        elif err.errno == errno.EPERM:
            raise SystemError('Permission failure when checking if PID is running - cannot detect job completion')
        else:
            raise  # EINVAL

    # no errors, so the process must be alive (since we can send signals to it)
    return True


class AsyncLocalEvaluator(AsyncFunctionEvaluator):
    """
    Class for cost functions that evaluated by launching a long-running process on the local machine.

    :param job_generator: callable that sets up the run directory for a given x (by e.g. writing config files). \
    It will be passed two arguments: the job directory and the point to evaluate at (x).
    :param run_cmd_generator: callable that returns the command to run the job. \
    It will be passed two arguments: the job directory and the point to evaluate at (x). \
    If run_cmd_generator returns a string, the string will be run by the default shell (typically /bin/sh) \
    via Popen with shell=True. If run_cmd_generator returns a list, it will be passed to Popen. \
    In both cases, the CWD is set to the job directory.
    :param parse_result: callable that returns the cost function evaluated at x. \
    It will be passed two arguments: the job directory and the point to evaluate a t (x).\
    This will be called after the command returned by run_cmd_generator has terminated (gracefully or otherwise). \
    If the process did not terminate successfully or the result is otherwise unavailable, parse_result \
    should raise any exception. This will signal the optimization routine to not try this point again.
    :param jobs_dir: optional base directory to run jobs in - default is $PWD/opt_jobs.
    :param required_fraction: fraction of points which must complete before continuing to the next iteration see \
    AsyncEvaluator for more info and implementation
    :param max_pending: maximum simultaneous processes, defaults to multiprocessing.cpu_count() see AsyncEvaluator for \
    implementation
    """
    def __init__(self,
                 job_generator: Callable[[str, np.array], None],
                 run_cmd_generator: Callable[[str, np.array], Union[str, List[Any]]],
                 parse_result: Callable[[str, np.array], float],
                 jobs_dir: str = os.path.join(os.getcwd(), 'opt_jobs'),
                 required_fraction=1.0,
                 max_pending=cpu_count()):
        super().__init__(required_fraction=required_fraction, max_pending=max_pending)
        self.job_generator = job_generator
        self.run_cmd_generator = run_cmd_generator
        self.parse_result = parse_result
        self.jobs_dir = jobs_dir

        # ignore SIGCHLD so we don't have to deal with zombie processes
        # https://stackoverflow.com/questions/16807603/python-non-blocking-non-defunct-process
        signal.signal(signal.SIGCHLD, signal.SIG_IGN)

    def start(self, x: np.array) -> (str, int):
        if not os.path.exists(self.jobs_dir):
            os.mkdir(self.jobs_dir)
        directory = tempfile.mkdtemp(prefix='job_', dir=self.jobs_dir)

        self.job_generator(directory, x)
        run_cmd = self.run_cmd_generator(directory, x)
        if type(run_cmd) is list:
            run_cmd = [str(s) for s in run_cmd]

        proc = subprocess.Popen(run_cmd, stdout=open(os.path.join(directory, 'OUT.LOG'), 'wb'),
                                stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                shell=(run_cmd is str), cwd=directory, start_new_session=True)
        pid = proc.pid
        return directory, pid

    def check_for_result(self, x: np.array, data: (str, int)) -> Union[ValueNotReady,
                                                                       EvaluationFailed,
                                                                       EvaluateAgain,
                                                                       float]:
        """
        Checks the status of pid, in data, and calls parse_result if the job is done.

        :param x: location of function evaluation
        :param data: list of directory and pid
        :return: either ValueNotReady float or EvaluationFailed()
        """
        pid = data[1]
        if _pid_is_running(pid):
            return ValueNotReady()

        try:
            directory = data[0]
            return self.parse_result(directory, x)
        except Exception as err:
            return EvaluationFailed(err)
