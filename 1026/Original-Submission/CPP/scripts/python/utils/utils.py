import copy
from enum import IntEnum

from scripts.python.utils.figure_config import markersize, linewidth, f_latex


class MttkrpMethod(IntEnum):
    MTTKRP = 0
    TWOSTEP0 = 1
    TWOSTEP1 = 2
    AUTO = 3
    CTF = 4
    PLANC = 5


def modes_string(modes, symbol='-'):
    mode_name = ""
    for i, m in enumerate(modes):
        if i != len(modes) - 1:
            mode_name += str(m) + symbol
        else:
            mode_name += str(m)
    return mode_name


def modes_title_string(modes):
    if f_latex:
        return modes_string(modes, symbol='$\\times$')
    else:
        return modes_string(modes, symbol='x')


def plot_gemm(gemm, ax, x, color='C3'):
    # Adjust x to start from little bit left and right of actual values of x
    x_copy = copy.deepcopy(x)
    x_copy[0] -= 0.1 * x_copy[-1]
    x_copy[-1] += 0.1 * x_copy[-1]

    # Plot gemm
    ax.axhline(gemm['mean'], ls='-', color=color, markersize=markersize, linewidth=linewidth, label='GEMM')
    ax.axhline(gemm['median'], ls='--', color=color, markersize=markersize, linewidth=linewidth)
    ax.fill_between(x_copy, gemm['mean'] - gemm['std'], gemm['mean'] + gemm['std'], alpha=0.2, color=color)
