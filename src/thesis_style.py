from kimmdy_paper_theme.kimmdy_paper_theme import mm_to_inch
import matplotlib.pyplot as plt
from plotnine import theme_set, theme_matplotlib
from kimmdy_paper_theme import (
    auto_init,
    convert_to_rgb,
)
from plotnine import * # type: ignore

# https://www.nature.com/documents/NRJs-guide-to-preparing-final-artwork.pdf

# A4 = 210 x 297 mm = 8.3 x 11.7 in
DPI = 300
double_column = mm_to_inch(180)
single_column = mm_to_inch(80)

auto_init()

HITS_DARKBLUE = "#003063"
HITS_SIGNALBLUE = "#004f9f"
HITS_CYAN = "#0088c2"
HITS_CYAN_LIGHT = "#55b4dc"
HITS_GREEN = "#019050"
HITS_GREEN_LIGHT = "#89b77a"
HITS_MAGENTA = "#c3006b"
HITS_MAGENTA_LIGHT = "#da7da6"
HITS_YELLOW = "#ffcc00"
HITS_YELLOW_LIGHT = "#ffe07d"

RED="#cc0130"
DARK_RED="#ad0028"


# https://www.mpikg.mpg.de/6339023/Logo-Guide-Print-_-Max-Planck-Gesellschaft.pdf
MPG_grey_dark = "#777777"
MPI_GREEN = "#006c66"
MPI_GREEN_SECONDARY = "#055000"
MPI_DARK_GREY = "#777777"

experiment = "#0088c2"
experiment_light = "#55b4dc"
kimmdy = "#c3006b"
kimmdy_light = "#da7da6"

rc = plt.rcParams

rc["ytick.right"] = False
rc["xtick.top"] = False
rc["axes.ymargin"] = 0

# rc["savefig.facecolor"] = "white"
# rc["figure.facecolor"] = "white"
# rc["axes.facecolor"] = "white"

default_colors=[
    HITS_MAGENTA,
    HITS_DARKBLUE,
    MPI_GREEN,
    MPI_GREEN_SECONDARY,
    HITS_SIGNALBLUE,
    HITS_CYAN,
    HITS_GREEN,
    HITS_YELLOW,
]

class custom_theme(theme_matplotlib):
    def __init__(self, rc=rc):
        super().__init__(rc=rc)
        self += theme(
            strip_background=element_rect(fill="white", color="black", size=1),
            text = element_text(color="black", size=7, family="Roboto")
        )

theme_set(custom_theme(rc=rc))

