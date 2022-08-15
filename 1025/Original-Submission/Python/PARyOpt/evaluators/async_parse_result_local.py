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
AsyncLocalParseResultEvaluator
An evaluator that only parses files periodically for function evaluations and
does not actually execute a script. It periodically checks for the file "job_folder/if_parse.txt" 
whether the function evaluation is complete or not. If the file contains anything other than 'False', 
it assumes the external function evaluation is completed and hence checks for "job_folder/y.txt" 
for the cost function value.
"""
import tempfile
from datetime import datetime, timedelta
import os

import numpy as np
from typing import Union, Callable

from PARyOpt.evaluators.paryopt_async import AsyncFunctionEvaluator, ValueNotReady, EvaluationFailed, EvaluateAgain


class AsyncLocalParseResultEvaluator(AsyncFunctionEvaluator):
    """
    Fills files in a set of folders and periodically checks if the external evaluator has finished \
    evaluation. Supports asynchronous evaluations. Only on a local machine.

    :var folder_num: folder into which the location values are written into
    :var job_generator: callable that sets up the run directory for a given x (by e.g. writing config files). \
    It will be passed two arguments: the job directory and the point to evaluate at (x).
    :var jobs_dir: base directory where the jobs are written into
    :var parse_result: function to parse result from directory , if not specified, it will search in \
    jobs_dir/folder_<folder_num>/if_parse.txt and y.txt
    :var total_folders: total number of folders to write into
    :var wait_time: min time to wait before parsing the results folder
    """

    def __init__(self,
                 parse_result: Callable[[str, np.array], float] = None,
                 job_generator: Callable[[str, np.array], None] = None,
                 jobs_dir: str = os.path.join(os.getcwd(), 'opt_jobs'),
                 wait_time: timedelta = timedelta(minutes=1),
                 total_folders: int = 16,
                 max_pending: int = 0,
                 required_fraction: float = 1.0):
        super().__init__(required_fraction=required_fraction, max_pending=max_pending)
        self.job_generator = job_generator
        self.parse_result = parse_result
        self.jobs_dir = jobs_dir
        self.wait_time = wait_time
        self.folder_num = 0
        self.total_folders = total_folders

    def start(self, x: np.array) -> (str, datetime):
        """
        function to start a cost function evaluation (write to file in this case) given the location \
        of evaluation

        :param x: location of evaluation
        :return: folder name in which it was submitted
        """
        if not os.path.exists(self.jobs_dir):
            os.mkdir(self.jobs_dir)
        folder_name = tempfile.mkdtemp(prefix='folder_', dir=self.jobs_dir)
        if callable(self.job_generator):
            self.job_generator(folder_name, x)
        else:
            # default method to fill folder with simulation related files
            with open(folder_name + '/x.txt', 'a') as f:
                f.write('{}\n'.format(x))
            self.folder_num += 1
            # cycle the folder numbers so that the folders are reused
            self.folder_num %= self.total_folders
        return folder_name, datetime.now()

    def check_for_result(self, x: np.array, data: (str, datetime)) -> Union[ValueNotReady,
                                                                            EvaluationFailed,
                                                                            EvaluateAgain,
                                                                            float]:
        """
        Function to check if a location has been evaluated or not (ValueNotReady). If it completes, categorize the \
        result as one of a float , EvaluationFailed or EvaluateAgain

        :param x: location to evaluate
        :param data: directory
        :return: Union[ValueNotReady, EvaluationFailed, EvaluateAgain, float]
        """
        job_folder = data[0]
        submit_time = data[1]

        if datetime.now() - submit_time < self.wait_time:
            return ValueNotReady()

        if self.parse_result is Callable:
            return self.parse_result(x, data)

        else:
            # default parsing method
            if_parse = "False"
            y_val = None
            # check if the folder has to be parsed for cost function value
            with open(job_folder + '/if_parse.txt', 'r') as f_parse:
                if_parse = f_parse.readline().strip()
            if if_parse == "False":
                # cost function evaluation not yet done
                return ValueNotReady()
            else:
                # cost function evaluation is done, parse y.txt
                with open(job_folder + '/y.txt', 'r') as f:
                    y_vals = f.readlines()
                    # read all the lines in y.txt and
                    # get the latest function value, which corresponds to the latest x.txt
                    y_val = float(y_vals[-1])
                if abs(y_val) > 1e5:
                    # evaluation has failed
                    return EvaluationFailed('parsed value is very large')
                else:
                    return float(y_val)
