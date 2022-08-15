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
Asynchronous evaluator super class
"""
import signal
from typing import Any, Union, Tuple, List

import numpy as np
import time
import pickle
import os
import math
import sys
import logging

# logging
logger = logging.getLogger(__name__)


class ValueNotReady:
    """
    Indicates a function value is not ready yet.
    """
    pass


class EvaluationFailed:
    """
    Indicates evaluation was not able to complete successfully, with an error value (i.e. an exception).
    """
    def __init__(self, reason: str):
        self.reason = reason

    def __str__(self):
        return str(self.reason)


class EvaluateAgain:
    """
    Indicates evaluation needs to be done again, due to some reason (for eg., during hardware failures)
    """
    def __init__(self, reason: str):
        self.reason = reason

    def __str__(self):
        return str(self.reason)


class _CriticalSection:
    """
    Context manager for temporarily disabling SIGINT (Ctrl-C) interrupts.
    """
    def __init__(self, msg: str = 'Program will exit after critical section.'):
        self.msg = msg
        self.interrupted = False
        self._prev_handler = None

    def _handler(self, signum, stack):
        self.interrupted = True
        if self.msg:
            print(self.msg, file=sys.stderr)

    def __enter__(self):
        self._prev_handler = signal.getsignal(signal.SIGINT)
        signal.signal(signal.SIGINT, self._handler)

    def __exit__(self, exc_type, exc_val, exc_tb):
        signal.signal(signal.SIGINT, self._prev_handler)
        self._prev_handler = None

        if self.interrupted:
            raise KeyboardInterrupt()


class AsyncFunctionEvaluator:
    """
    Abstract base class for long-running cost functions (e.g. external simulations).
    Must be subclassed. Subclasses should fill in start() and check_for_results().
    Automatically saves state as jobs are submitted.
    """
    def __init__(self, required_fraction: float=1.0, max_pending: int=0):
        """

        :param required_fraction: percentage of points that must have completed (or failed)
                                    before evaluate_population will return.
        :param max_pending: maximum simultaneous pending points (cost function evaluations)
                            this is an inclusive upper bound (e.g. max_pending=8 means 8 points can eval simultaneously)
                            use 0 for no limit
        """
        # history is a dictionary of x -> (result, data) pairs.
        # res is the result returned by check_for_result, data is the data returned by
        # start(). The start() data is only kept for logging purposes.
        self.history = {}

        # pending is a dictionary of x -> data.
        # data is whatever user data is returned by start().
        # By having an entry in pending, it is assumed that an evaluation for x has been
        # started, and a result will eventually be returned by check_for_result().
        self.pending = {}

        self.save_path = "async_state.pickle"
        if os.path.exists(self.save_path):
            self._load()

        self._req_frac = required_fraction
        assert (self._req_frac > 0.0) and (self._req_frac <= 1.0)

        self._max_pending = max_pending

    def __call__(self, xs: List[np.array], old_xs: List[np.array] = list()):
        return self.evaluate_population(xs, old_xs)

    def evaluate_population(self, xs: List[np.array], if_ready_xs: List[np.array] = list())\
            -> Tuple[List[Tuple[np.array, float]], List[np.array], List[np.array]]:
        """
        Evaluates a population of x values, encoded as a list of 1D numpy arrays.
        Returns a tuple containing three lists:

        - Completed values: [ (x1, y1), (x2, y2), ... ]
        - Pending values - evaluation is in progress, but not complete: [ x1, x2, ... ]
        - Failed  values - evaluation completed unsuccessfully: [ x1, x2, ... ]

        The union of completed, failed, and pending is equal to the union of xs and if_ready_xs.
        The __init__ parameter required_fraction tunes how many completed/pending values are returned.

        :param xs: list of new points to check
        :param if_ready_xs: List of points to include in the return tuple if they available by the time we
                            evaluate the minimum required percentage of xs. These points do not count towards
                            the minimum required completed points.
        :return: ( [(x, y), ...] completed, [x, ...] pending, [x, ...] failed )
        """
        completed = []
        pending = [(True, x) for x in xs] + [(False, x) for x in if_ready_xs]
        failed = []
        n_completed_new = 0
        n_req_pts = math.ceil(self._req_frac * len(xs))
        check_iter = 0
        while True:
            check_iter += 1
            still_pending = []
            for is_new_pt, x in pending:
                y = self._get_point(x)
                if type(y) is float:
                    logger.info("X completed: {}, value: {}".format(x, y))
                    completed.append((x, y))
                    if is_new_pt:
                        n_completed_new += 1
                elif type(y) is EvaluationFailed:
                    logger.info("X failed: {} because: {}".format(x, str(y)))
                    failed.append(x)
                    if is_new_pt:
                        n_completed_new += 1
                elif type(y) is ValueNotReady:
                    logger.debug("X not ready: {}".format(x))
                    still_pending.append((is_new_pt, x))
                elif type(y) is EvaluateAgain:
                    logger.info("Re-evaluating : {}".format(x))
                    self._remove_point(x)
                    still_pending.append((is_new_pt, x))
                else:
                    raise TypeError("y: expected float, EvaluationFailed, or ValueNotReady - got " + str(type(y)))
            pending = still_pending

            # sleep if we still aren't done
            # sleep time increases with check iteration and eventually asymptotes at 10 sec
            if n_completed_new < n_req_pts:
                time.sleep(9 * np.tanh(check_iter) + 1.)
            else:
                break

        with _CriticalSection('The program will exit once the async evaluator finishes saving.'):
            self._save()

        return completed, [x[1] for x in pending], failed

    def _get_point(self, x: np.array) -> Union[ValueNotReady,
                                               EvaluationFailed,
                                               EvaluateAgain,
                                               float]:

        # convert np.array to tuple so it can be used as a dict key
        # only used for dictionaries - we still pass the np.array to the child functions
        x_key = tuple(x.tolist())

        if x_key in self.history:
            return self.history[x_key][0]
        elif x_key in self.pending:
            res = self.check_for_result(x, self.pending[x_key])
            if type(res) is float:
                data = self.pending.pop(x_key)
                self.history[x_key] = (res, data)
            return res
        else:
            if 0 < self._max_pending <= len(self.pending):
                return ValueNotReady()

            with _CriticalSection('The program will exit once job submission completes.'):
                self.pending[x_key] = self.start(x)
                self._save()
            return self._get_point(x)  # recurse once, will take pending case

    def _remove_point(self, x: np.array):
        """
        Removes a point from the history and pending dictionaries, if it exists, then saves the new state to disk.
        It is safe to call this function even if x has not been submitted.

        :param x: point to remove
        """
        x_key = tuple(x.tolist())

        self.history.pop(x_key, None)
        self.pending.pop(x_key, None)
        with _CriticalSection('The program will exit once the async evaluator finishes saving.'):
            self._save()

    def _save(self):
        data = {
            "pending": self.pending,
            "history": self.history,
        }
        with open(self.save_path, 'wb') as f:
            pickle.dump(data, f)

    def _load(self):
        with open(self.save_path, 'rb') as f:
            data = pickle.load(f)
        self.pending = data["pending"]
        self.history = data["history"]

        # remove failed entries, in case the user has fixed things since the last run
        # kept_history = {}
        # for k, v in self.history:
        #     if type(v) is EvaluationFailed:
        #         kept_history[k] = v
        #     else:
        #         print('Removing previously failed point "' + str(k) + '" from history (to allow re-evaluation)')

    def start(self, x: np.array) -> Any:
        """
        Start a cost function evaluation at the given x location.
        This method may return anything - the data will be passed on to check_for_result().
        The only restriction is that the return value should be pickle-able to enable restart support.

        :param x: point to begin evaluation at
        :return: user data that will be fed into check_for_result()
        """
        raise NotImplementedError()

    def check_for_result(self, x: np.array, data: Any) -> Union[ValueNotReady, EvaluationFailed, EvaluateAgain, float]:
        """
        Returns the cost function evaluation at x, if the value is available. This method is only called after start.

        :param x: the point to evaluate the cost function at
        :param data: user data returned by start(x)
        :return: an instance of ValueNotReady if such, EvaluationFailed(reason), or the cost function value at x (float)
        """
        raise NotImplementedError()
