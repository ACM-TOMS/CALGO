"""
Simple setup of a time-multigrid hierarchy for a problem.

Creates a time-multigrid hierarchy using the finest
problem, the number of levels, and the coarsening factor.
"""

import copy
import warnings
from typing import List

from pymgrit.core.application import Application


def simple_setup_problem(problem: Application, level: int, coarsening: int) -> List[Application]:
    """
    Simple setup of a time-multigrid hierarchy for a problem.

    Creates a time-multigrid hierarchy using the finest
    problem, the number of levels, and the coarsening factor.

    :param problem: Application problem on the finest grid
    :param level: Number of time-grid levels
    :param coarsening: Coarsening factor to be used for all levels
    :return: List of application problems; one application problem per time-grid level
    """
    problem_structure = [problem]

    if len(problem.t[::coarsening * level]) == 1:
        warnings.warn(
            "This choice leads to a coarsest grid with only one time point, which is the initial point. "
            "It is recommended to choose a structure with at least two points on the coarsest grid.")

    for i in range(level - 1):
        problem_tmp = copy.deepcopy(problem)
        tmp_t = problem_structure[-1].t[::coarsening]
        problem_tmp.t_start = tmp_t[0]
        problem_tmp.t_end = tmp_t[-1]
        problem_tmp.t = tmp_t
        problem_tmp.nt = len(tmp_t)
        problem_structure.append(problem_tmp)

    return problem_structure
