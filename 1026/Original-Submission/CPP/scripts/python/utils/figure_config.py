import matplotlib
import matplotlib.pyplot as plt
plt.switch_backend('Qt5Agg')

# General configurations
linewidth = 1
markersize = 2
colors = ['C0', 'C1', 'C3']

fig_size_in = {'width': 5.47807, 'height': 6}
fig_format = '.png'

# f_latex = True
f_latex = False
if f_latex:
    matplotlib.use("pgf")
    matplotlib.rcParams.update({'text.usetex': True})
    matplotlib.rcParams.update({'pgf.rcfonts': False})
    matplotlib.rcParams.update({'font.family': 'serif'})
    matplotlib.rcParams.update({'pgf.texsystem': 'xelatex'})
    # matplotlib.rcParams.update({'font.size': 22})
    # matplotlib.rcParams.update({'figure.dpi': 300})

    fig_size_in = {'width': 5.47807, 'height': 3}
    fig_format = '.pgf'