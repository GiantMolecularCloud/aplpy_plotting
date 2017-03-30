#####################################################################
#                          APLPY PLOTTING                           #
#####################################################################

"""
These functions will produce plots of channel maps, moment maps,
pV diagrams, ... in a quality that (hopefully) allows publishing.
"""

###################################################################################################


import os as __os__
import aplpy as __aplpy__
import numpy as __np__
from astropy.coordinates import SkyCoord as __SkyCoord__
from astropy import units as __u__
from astropy.coordinates import Angle as __Angle__
import matplotlib as __mpl__
import matplotlib.colors as __colors__
import matplotlib.pyplot as __plt__
from matplotlib import rc as __rc__
__rc__('text',usetex=True)
from matplotlib.cbook import MatplotlibDeprecationWarning as __MatplotlibDeprecationWarning__
import warnings as __warnings__
__warnings__.simplefilter('ignore', __MatplotlibDeprecationWarning__)

###################################################################################################

# General style settings
########################

# This is defined here to ensure consistency of the output over subsequent calls, i.e. I want my 
# plots to look the same. Some things must be adapted to the dataset and should thus be changeable 
# before calling the plot function. Global variables do exactly that.

# settings that are dataset dependend, i.e. need to be changed for a new dataset
global tick_label_xformat; tick_label_xformat = 'hh:mm:ss.s'
global tick_label_yformat; tick_label_yformat = 'dd:mm:ss'
global ticks_xspacing; ticks_xspacing = __Angle__('0 1 0', unit='hourangle')
global ticks_yspacing; ticks_yspacing = 10.0*__u__.arcsec
global ticks_minor_frequency; ticks_minor_frequency = 5



# fixed settings, i.e. settings that should be the same for each and every plot
_velo_fontsize         = 10.0		# unit: point
_colorbar_fontsize     = 10.0		# unit: point
#_colorbar_width        = 0.15	    # relative to panel size
_scalebar_frame        = False
_scalebar_linestyle    = 'solid'		# or any other plt.plot linestyle
_scalebar_linewidth    = 2			# unit: points
_scalebar_color        = 'red'		# any named color or mpl.color instance
_scalebar_fontsize     = 10.0    	# only used in channel map to prevent bar sliding over map
_beam_frame            = False
_beam_color            = 'black'
_ticks_color           = 'black'		# this setting overrules the matplotlibrc defaults



# define new viridis colormap with less dark blue
viridis_cropped = __colors__.ListedColormap(__mpl__.cm.viridis(__np__.linspace(0.1,1.0,100)))


###################################################################################################

# import plotting functions
from .aplpy_plot import aplpy_plot
from .aplpy_plot_pv import aplpy_plot_pv
from .aplpy_channel_map import aplpy_channel_map
from .aplpy_map_grid import aplpy_map_grid