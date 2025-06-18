import matplotlib.pyplot as plt
from plotnine import theme_set, theme_matplotlib
from kimmdy_paper_theme import plot_colors, auto_init, single_column, double_column, convert_to_rgb

# https://www.nature.com/documents/NRJs-guide-to-preparing-final-artwork.pdf

auto_init()

DPI = 300
rc = plt.rcParams
default_size = rc["figure.figsize"]
plot_colors = plot_colors

rc["ytick.right"] = False
rc["xtick.top"] = False
rc["axes.ymargin"] = 0


# rc["savefig.facecolor"] = "white"
# rc["figure.facecolor"] = "white"
# rc["axes.facecolor"] = "white"

theme_set(theme_matplotlib(rc=rc))

HITS_DARKBLUE = '#003063'
HITS_SIGNALBLUE = '#004f9f'
HITS_CYAN = '#0088c2'
HITS_CYAN_LIGHT = '#55b4dc'
HITS_GREEN = '#019050'
HITS_GREEN_LIGHT = '#89b77a'
HITS_MAGENTA = '#c3006b'
HITS_MAGENTA_LIGHT = '#da7da6'
HITS_YELLOW = '#ffcc00'
HITS_YELLOW_LIGHT = '#ffe07d'
MPI_GREEN = '#006c66'
MPI_GREEN_SECONDARY = '#055'
MPG_grey_dark = '#777777'
experiment = '#0088c2'
experiment_light = '#55b4dc'
kimmdy = '#c3006b'
kimmdy_light = '#da7da6'
