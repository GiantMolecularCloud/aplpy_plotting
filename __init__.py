#####################################################################
#                          APLPY PLOTTING                           #
#####################################################################

"""
These functions will produce plots of channel maps, moment maps,
pV diagrams, ... in a quality that (hopefully) allows publishing.
"""

###################################################################################################


import os
import aplpy
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import Angle
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text',usetex=True)
from matplotlib.cbook import MatplotlibDeprecationWarning
import warnings
warnings.simplefilter('ignore', MatplotlibDeprecationWarning)

###################################################################################################

# General style settings
########################

# This is defined here to ensure consistency of the output over subsequent calls, i.e. I want my 
# plots to look the same. Some things must be adapted to the dataset and should thus be changeable 
# before calling the plot function. Global variables do exactly that.

# settings that are dataset dependend, i.e. need to be changed for a new dataset
tick_label_xformat = 'hh:mm:ss.ss'
tick_label_yformat = 'dd:mm:ss.ss'
ticks_xspacing = Angle('0 1 0', unit='hourangle')
ticks_yspacing = 1.0*u.arcmin
ticks_minor_frequency = 5



# fixed settings, i.e. settings that chould be the same for each and every plot
velo_fontsize         = 10.0		# unit: point
colorbar_fontsize     = 10.0		# unit: point
colorbar_width        = 0.15	    # relative to panel size
scalebar_frame        = False
scalebar_linestyle    = 'solid'		# or any other plt.plot linestyle
scalebar_linewidth    = 2			# unit: points
scalebar_color        = 'red'		# any named color or mpl.color instance
scalebar_fontsize     = 10.0    	# only used in channel map to prevent bar sliding over map
beam_frame            = False
beam_color            = 'black'
ticks_color           = 'black'		# this setting overrules the matplotlibrc defaults



# define new viridis colormap with less dark blue
viridis_cropped = colors.ListedColormap(mpl.cm.viridis(np.linspace(0.1,1.0,100)))


###################################################################################################

# import plotting functions
from .aplpy_plot import aplpy_plot
from .aplpy_plot_pv import aplpy_plot_pv
from .aplpy_channel_map import aplpy_channel_map
from .aplpy_map_grid import aplpy_map_grid