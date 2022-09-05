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
from typing import Callable, List, Tuple
import numpy as np

from . import async_local
from . import paryopt_async
from . import async_parse_result_local
from . import async_sbatch
from . import connection

from .async_local import AsyncLocalEvaluator
try:
    from .async_sbatch import AsyncSbatchEvaluator
except ImportError:
    AsyncSbatchEvaluator = None
from .async_parse_result_local import AsyncLocalParseResultEvaluator


class FunctionEvaluator:
    """
    The simplest function evaluator - evaluates a Python function for each point.

    :param func: cost function to evaluate at each x
    """
    def __init__(self, func: Callable[[np.array], float]):
        self.evaluate = func

    def __call__(self, xs: List[np.array], if_ready_xs: List[np.array] = list()):
        return self.evaluate_population(xs, if_ready_xs)

    def evaluate_population(self, xs: List[np.array], old_xs: List[np.array] = list()) \
            -> Tuple[List[Tuple[np.array, float]], List[np.array], List[np.array]]:
        xs += old_xs
        return [(x, self.evaluate(x)) for x in xs], [], []
