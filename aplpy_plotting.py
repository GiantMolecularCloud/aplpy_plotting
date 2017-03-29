#####################################################################
#                          APLPY PLOTTING                           #
#####################################################################
# These functions will produce plots of channel maps, moment maps,  #
# pV diagrams, ... in a quality that (hopefully) allows publishing. #
#####################################################################

# I know that I shouldn't use dozens of if-else statements but rather try-except to get correct 
# handling of exceptions. When I started wrinting this code I didn't know what exception could
# do for my case and right now I don't have time to change this script.

###################################################################################################

from __future__ import division
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

# General plotting settings (common to all plots)
# This is defined here to ensure consistency of the output over subsequent calls, i.e. I want my 
# plots to look the same. Some things must be adapted to the dataset and should thus be changeable 
# when calling the plot function. Global variables do exactly that.

# fixed settings
velo_fontsize         = 10.0
colorbar_fontsize     = 10.0
colorbar_width        = 0.15    # relative to panel size
scalebar_frame        = False
scalebar_linestyle    = 'solid'
scalebar_linewidth    = 2
scalebar_color        = 'red'
scalebar_fontsize     = 10.0    # only used in channel map to prevent bar sliding over map
beam_frame            = False
beam_color            = 'black'
ticks_color           = 'black'

# user-adjustable settings come with a default value
global tick_label_xformat; tick_label_xformat = 'hh:mm:ss.ss'
global tick_label_yformat; tick_label_yformat = 'dd:mm:ss.ss'
global ticks_xspacing; ticks_xspacing = Angle('0 1 0', unit='hourangle')
global ticks_yspacing; ticks_yspacing = 1.0*u.arcmin
global ticks_minor_frequency; ticks_minor_frequency = 5

# define new viridis colormap with less dark blue
global viridis_cropped; viridis_cropped = colors.ListedColormap(mpl.cm.viridis(np.linspace(0.1,1.0,100)))

###################################################################################################

def aplpy_plot(fitsfile, **kwargs):

    """
    aplpy_plot (main plotting function)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Provides a simple interface to APLpy plots. A single filename is enough to get a simple plot.

    Mandatory unnamed arguments:
        fitsfile    Path and file name of the fits image to be plotted
        

    Optional arguments:
        out         Path and file name of the created plot.
                    If not specified, the plot will be saved where the input image
                    is located. Default format: png
        
        cmap        Colormap to plot the image.
                    If not specified, the matplotlib default will be used, usually
                    this is viridis.
                    Every named matplotlib colormap can be used or any matplotlib 
                    colormap object.
                    For grayscale use cmap='grayscale'
                    
        vmin        Minimum value for colormap normalization.
        vmax        Maximum value for colormap normalization.
        
        stretch     Colormap strech, e.g. 'linear', 'log', 'sqrt', 'arcsinh'.
                    Linear and log scaling are fully implemented, other scaling can
                    cause errors when unusual argument combinations are chosen.
        
        invert      Invert the colormap if True.
        
        recenter    Center the image on this location and set image width/height.
                    You either have to specify radius or width and height.
                    Must be an array containing an astropy.SkyCoord object plus one 
                    or two angular distances: [SkyCoord(...), 1*u.arcmin]
                    
        contour     List of contour elements.
                    Each contour element must be a list of 'image file', list of 
                    contour levels and list of colors for each contour. If only one
                    color is given, it will be used for all contours.
        clabel      Draw labels on the contours? Needs contours to be present.
                    Either set to True to use the matplotlib defaults or give a 
                    dictionary with allowed axes.clabel arguments.
                    Use clabel={'fmt': '%i'} to get labels without trailing zeros.
                    
        colorbar_location   As the name says.
                            Can be, e.g. 'bottom', 'right', ...
        colorbar_label      Can be specified only when colorbar_location is given.
                            String containing the label.
                            
        scalebar_length     As the name says. Must be an angular astropy.units 
                            object, e.g. 10.0*u.arcmin
        scalebar_label      String containing the label to be plotted below the 
                            scalebar.
        scalebar_corner     Where should the scalebar and label be plotted?
                            E.g. 'bottom left'
        
        beam_corner Plot a beam ellipse as given in the fits header in this corner 
                    of the plot. Does not have to be a corner. 'bottom' also works.
        
        overlay     List of overlay elements.
                    Each overlay element is a list of shape (plt.scatter type), 
                    astropy.SkyCoord object specifying the position, size as 
                    angular astropy.unit and further plt.scatter kwargs (can be 
                    empty).

    General style settings
        Settings that do not have to be changed for each plot but maybe once per 
        script or once per project. Often used ones are tick_label_xformat, 
        ticks_xspacing and the corresponding settings for y.
            Can be accessed via
        import aplpy_plotting as ap
        ap.setting = ...
        See the aplpy_plotting.py main file the exact setting names if you need to 
        change them.


    example:

    aplpy_plot('abc.fits', 
            out  = 'abc.png',       
            cmap = 'jet', 
            vmin = 0, 
            vmax = 1, 
            stretch  = 'linear', 
            recenter = [SkyCoord('01h23m45.6s 12d34m45.6s'), 1.0*u.armin], 
            contour  = ['contour.fits', [1,2,3], ['white', 'grey', 'black']], 
            colorbar_location = 'right', 
            colorbar_label  = 'Jy/beam', 
            scalebar_length = 1.0, 
            scalebar_label  = 'string', 
            scalebar_corner = 'bottom', 
            beam_corner     = 'bottom left',
            overlay  = [['circle', SkyCoord('01h23m45.6s 12d34m45.6s'), 1.0*u.arcmin, {'linewidth': 1.0}]]
            )
    """

    
    print "--> plotting map "+fitsfile
    
    fig = aplpy.FITSFigure(fitsfile)
    
    if 'cmap' in kwargs:
        if kwargs['cmap'] == 'grayscale':
            fig.show_grayscale()
            if 'stretch' in kwargs:
                fig.show_grayscale(stretch=kwargs['stretch'])
                if 'vmin' and 'vmax' in kwargs:
                    fig.show_grayscale(vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'])
                    if 'invert' in kwargs:
                        fig.show_grayscale(vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'], invert=kwargs['invert'])
        else:
            fig.show_colorscale(cmap=kwargs['cmap'])
            if 'stretch' in kwargs:
                fig.show_colorscale(cmap=kwargs['cmap'], stretch=kwargs['stretch'])
                if 'vmin' and 'vmax' in kwargs:
                    fig.show_colorscale(cmap=kwargs['cmap'], vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'])
                    if 'invert' in kwargs:
                        fig.show_colorscale(cmap=kwargs['cmap'], vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'], invert=kwargs['invert'])
    else:
        fig.show_colorscale()
    
    # recenter image
    if 'recenter' in kwargs:
        if (len(kwargs['recenter']) == 2):
            fig.recenter(kwargs['recenter'][0].ra.degree, kwargs['recenter'][0].dec.degree, radius=kwargs['recenter'][1].to(u.degree).value)
        elif (len(kwargs['recenter']) == 3):
            fig.recenter(kwargs['recenter'][0].ra.degree, kwargs['recenter'][0].dec.degree, width=kwargs['recenter'][1].to(u.degree).value, height=kwargs['recenter'][2].to(u.degree).value)
        else:
            print "--> specify SkyCoord(x,y) and either radius or width, height. not recentering"
    
    # contours
    if 'contour' in kwargs:
        for cont_i in np.arange(len(kwargs['contour'])):
            if len(kwargs['contour'][cont_i]) == 3:
                fig.show_contour(data=kwargs['contour'][cont_i][0], levels=kwargs['contour'][cont_i][1], colors=kwargs['contour'][cont_i][2])
            else:
                print "--> wrong number or format of contour parameters in image "+str(cont_i)+". not plotting contours"
            if 'clabel' in kwargs:
                if kwargs['clabel'] == True:
                    fig._layers['contour_set_'+str(cont_i+1)].clabel()
                if isinstance(kwargs['clabel'],dict):
                    fig._layers['contour_set_'+str(cont_i+1)].clabel(**kwargs['clabel'])
    
    # colorbar settings
    if 'colorbar_location' in kwargs:
        fig.add_colorbar()
        fig.colorbar.show()
        fig.colorbar.set_location(kwargs['colorbar_location'])
        if 'colorbar_label' in kwargs:
            fig.colorbar.set_axis_label_text(kwargs['colorbar_label'])
        if 'stretch' in kwargs:
            if (kwargs['stretch'] == 'log'):
                log_ticks = [float('{:.2f}'.format(round(x,int(-1*np.log10(kwargs['vmin']))))) for x in np.logspace(np.log10(kwargs['vmin']),np.log10(kwargs['vmax']),num=10, endpoint=True)]
                fig.colorbar.set_ticks(log_ticks)

    # scale bar
    if 'scalebar_length' and 'scalebar_label' and 'scalebar_corner' in kwargs:
        fig.add_scalebar(length=kwargs['scalebar_length'].to(u.degree).value, label=kwargs['scalebar_label'], corner=kwargs['scalebar_corner'], frame=scalebar_frame)
        fig.scalebar.set_linestyle(scalebar_linestyle)
        fig.scalebar.set_linewidth(scalebar_linewidth)
        fig.scalebar.set_color(scalebar_color)

    # beam settings
    if 'beam_corner' in kwargs:
        fig.add_beam()
        fig.beam.show()
        fig.beam.set_corner(kwargs['beam_corner'])
        fig.beam.set_frame(beam_frame)
        fig.beam.set_color(beam_color)
    
    # data set overlay
    if 'label_text' in kwargs:
        fig.add_label(0.5, 0.9, kwargs['label_text'][i].replace('_','$\_$'), color='black', relative=True, size=velo_fontsize)
        if 'label_kwargs' in kwargs:
            fig.add_label(0.5, 0.9, kwargs['label_text'][i].replace('_','$\_$'), color='black', relative=True, size=velo_fontsize, **kwargs['label_kwargs'])
    
    if 'overlay' in kwargs:
        for olay in np.arange(len(kwargs['overlay'])):
            if (kwargs['overlay'][olay][0] == 'circle'):
                fig.show_circles(xw=kwargs['overlay'][olay][1].ra.degree, yw=kwargs['overlay'][olay][1].dec.degree, radius=kwargs['overlay'][olay][2].to(u.degree).value, **kwargs['overlay'][olay][3])
            else:
                fig.show_markers(xw=kwargs['overlay'][olay][1].ra.degree, yw=kwargs['overlay'][olay][1].dec.degree, marker=kwargs['overlay'][olay][0], s=kwargs['overlay'][olay][2], **kwargs['overlay'][olay][3])

    # ticks + labels
    fig.tick_labels.show()
    fig.tick_labels.set_xformat(tick_label_xformat)
    fig.tick_labels.set_yformat(tick_label_yformat)
    fig.ticks.show()
    fig.ticks.set_xspacing(ticks_xspacing.to(u.degree).value)
    fig.ticks.set_yspacing(ticks_yspacing.to(u.degree).value)
    fig.ticks.set_minor_frequency(ticks_minor_frequency)
    fig.ticks.set_color(ticks_color)

    if 'out' in kwargs:
        fig.save(kwargs['out'], dpi=300, transparent=True)
        print "--> saved file as "+kwargs['out']
    else:
        fig.save(os.path.splitext(fitsfile)[0]+'.png', dpi=300, transparent=True)
        print "--> saved plot as "+os.path.splitext(fitsfile)[0]+'.png'



###################################################################################################

def aplpy_channel_map(fitscube, ncols, nrows, chan_start, chan_iter, **kwargs):
    
    """
    aplpy_channel: channel map plotting function
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Plot a grid of channel maps, i.e. a regular grid of plots displaying every 
    n-th channel of a three dimensional fits cube.
    Axis labels are plotted on the bottom left panel only. The same applies
    to beam and scalebar (if given). The last panel always contains a  
    horizontal colorbar.

    Mandatory unnamed, ordered arguments:
        fitsfile    Path and file name of the fits cube to be plotted. Remove
                    degenerate axes such as the Stokes axis.
        
        ncols       number of columns of the grid
        nrows       number of rows of the grid
        
        chan_start  Channel number of first channel to plot in top left corner.
        
        chan_iter   Plot every chan_iter channel.
        

    Optional arguments:
        chan_width  Channel width in km/s.
        chan_velo0  Velocity of first channel in the cube.
                    These arguments are used to calculate thevelocity of each 
                    channel map. Will be replaced in the future by the information 
                    giben in the fits header but has to be specified manually at
                    the moment.
    
        out         Path and file name of the created plot.
                    If not specified, the plot will be saved where the input image
                    is located. Default format: png
                    
        figsize     Fiugre size as tpuel in inches. Default: (8.27, 11.69) A4 size
        
        cmap        Colormap to plot the image.
                    If not specified, the matplotlib default will be used, usually
                    this is viridis.
                    Every named matplotlib colormap can be used or any matplotlib 
                    colormap object.
                    For grayscale use cmap='grayscale'
                    
        vmin        Minimum value for colormap normalization.
        vmax        Maximum value for colormap normalization.
        
        stretch     Colormap strech, e.g. 'linear', 'log', 'sqrt', 'arcsinh'.
                    Linear and log scaling are fully implemented, other scaling can
                    cause errors when unusual argument combinations are chosen.
        
        recenter    Center the image on this location and set image width/height.
                    You either have to specify radius or width and height.
                    Must be an array containing an astropy.SkyCoord object plus one 
                    or two angular distances: [SkyCoord(...), 1*u.arcmin]
                    
        contour     List of contour elements.
                    Each contour element must be a list of mask type, image file, list
                    of contour levels and list of colors for each contour. The mask 
                    type can be either 'cubemask' or 'pixelmask'. A cube mask is a fits 
                    cube with the channel corresponding to the image channel is plotted
                    as contour, i.e. the contours can be different for each channel.
                    Pixel mask specifies a single plane image, i.e. the contours are 
                    identical in each panel. If only one color is given, it will be 
                    used for all contours.
                    Might be replaced be automatic detection if cube or single plane
                    in the future.
                    
        colorbar_location   As the name says.
                            Can be, e.g. 'bottom', 'right', ...
        colorbar_cmap       The colorbar is generated via matplotlib and thus
                            specified manually at the moment. Make sure you use
                            the same colorbar as in the cmap argument.
                            Will be updated in a way that cmap is enough and this 
                            argument will be not necessary anymore.
        colorbar_label      Can be specified only when colorbar_location is given.
                            String containing the label.
                            
        scalebar_length     As the name says. Must be an angular astropy.units 
                            object, e.g. 10.0*u.arcmin
        scalebar_label      String containing the label to be plotted below the 
                            scalebar.
        scalebar_corner     Where should the scalebar and label be plotted?
                            E.g. 'bottom left'
        
        beam_corner Plot a beam ellipse as given in the fits header in this corner 
                    of the plot. Does not have to be a corner. 'bottom' also works.
    
    
    General style settings
        Settings that do not have to be changed for each plot but maybe once per 
        script or once per project. Often used ones are tick_label_xformat, 
        ticks_xspacing and the corresponding settings for y.
            Can be accessed via
        import aplpy_plotting as ap
        ap.setting = ...
        See the aplpy_plotting.py main file the exact setting names if you need to 
        change them.


    example:
    
    aplpy_channel_map(fitscube, 
                      ncols, 
                      nrows, 
                      chan_start, 
                      chan_iter, 
                      chan_width=1.0, 
                      chan_vel0=0.0, 
                      cmap='jet', 
                      vmin=-0.05, 
                      vmax=3.5, 
                      stretch='linear', 
                      recenter = [SkyCoord('01h23m45.6s 12d34m45.6s'), 10.0*u.arcsec], 
                      contour=['pixelmask', 'contour.fits', [1,2,3], ['white', 'grey', 'black']], 
                      beam_corner='bottom left', 
                      colorbar_cmap=mpl.cm.jet, 
                      colorbar_label='flux density [Jy/beam]', 
                      scalebar_length=1.0, 
                      scalebar_label='string', 
                      scalebar_corner='bottom',
                      out='cube.png'
                      )
    """
    
    print "--> plotting channel map "+fitscube
    if not ('chan_width' and 'chan_vel0' in kwargs):
        print "--> channel width and velocity of first channel not given. will not calculate velocity"
    
    if 'figsize' in kwargs:
        main_figsize = kwargs['figsize']
    else:
        main_figsize = (8.27, 11.69)    # A4 in inches
    
    main_fig = plt.figure(figsize=main_figsize)
    for i in np.arange(nrows*ncols):
        
        # get subplot specific info
        subplt_size = [0.05+(i%ncols)*0.9/ncols, 0.95-np.ceil((i+1)/ncols)*0.9/nrows, 0.9/ncols, 0.9/nrows]
        #chn_slice   = [chan_start+i*chan_iter,0]   # if 4th axis is present
        chn_slice   = [chan_start+i*chan_iter]  # if no 4th axis available

        # introduce a check here to assure that the cube has only three axis (no forth degenerate axis, like Stokes)
        
        if 'chan_width' and 'chan_vel0' in kwargs:
            chn_velo    = (i*chan_iter+chan_start)*kwargs['chan_width']+kwargs['chan_vel0']
        
        # plot channel map if not last panel
        if ( i < nrows*ncols-1 ):
            print "--> panel "+str(i+1)+" of "+str(nrows*ncols)
            
            fig = aplpy.FITSFigure(fitscube, figure=main_fig, subplot=subplt_size, dimensions=[0,1], slices=chn_slice)
            if 'vmin' and 'vmax' in kwargs:
                if 'stretch' in kwargs:
                    if 'cmap' in kwargs:
                        if kwargs['cmap'] == 'grayscale':
                            fig.show_grayscale(vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'])
                        else:
                            fig.show_colorscale(cmap=kwargs['cmap'], vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'])
                    else:
                        fig.show_colorscale(vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'])
                else:
                    fig.show_colorscale(vmin=kwargs['vmax'], vmax=kwargs['vmax'])
            else:
                fig.show_colorscale()
            
            # recenter image
            if 'recenter' in kwargs:
                if (len(kwargs['recenter']) == 2):
                    fig.recenter(kwargs['recenter'][0].ra.degree, kwargs['recenter'][0].dec.degree, radius=kwargs['recenter'][1].to(u.degree).value)
                elif (len(kwargs['recenter']) == 3):
                    fig.recenter(kwargs['recenter'][0].ra.degree, kwargs['recenter'][0].dec.degree, width=kwargs['recenter'][1].to(u.degree).value, height=kwargs['recenter'][2].to(u.degree).value)
                else:
                    print "--> specify SkyCoord(x,y) and either radius or width, height. not recentering"
            
            # contours?
            if 'contour' in kwargs:
                if len(kwargs['contour']) == 4:
                    if kwargs['contour'][0] == 'pixelmask':
                        fig.show_contour(data=kwargs['contour'][1], levels=kwargs['contour'][2], colors=kwargs['contour'][3])
                    elif kwargs['contour'][0] == 'cubemask':
                        fig.show_contour(data=kwargs['contour'][1], dimensions=[0,1], slices=chn_slice, levels=kwargs['contour'][2], colors=kwargs['contour'][3])
                    else:
                        print "--> wrong mask type specification"
                else:
                    print "--> wrong number or format of contour parameters. not plotting contours"
            
            
            # ticks + labels
            fig.axis_labels.hide()
            fig.tick_labels.hide()
            fig.ticks.show()
            fig.ticks.set_xspacing(ticks_xspacing.to(u.degree).value)
            fig.ticks.set_yspacing(ticks_yspacing.to(u.degree).value)
            fig.ticks.set_minor_frequency(ticks_minor_frequency)
            fig.ticks.set_color(ticks_color)
            
            # velocity or channel number overlay
            if 'chan_width' and 'chan_vel0' in kwargs:
                fig.add_label(0.8, 0.8, str('{:3.1f}'.format(chn_velo))+r'\,km/s', color='black', relative=True, size=velo_fontsize)
            else:
                fig.add_label(0.8, 0.8, 'channel '+str(chn_slice), color='black', relative=True, size=velo_fontsize)
            
            # beam settings
            if 'beam_corner' in kwargs:
                fig.add_beam()
                fig.beam.show()
                fig.beam.set_corner(kwargs['beam_corner'])
                fig.beam.set_frame(beam_frame)
                fig.beam.set_color(beam_color)
            
        # colorbar settings
        if (i == ncols*nrows-1):
            if 'colorbar_cmap' and 'colorbar_label' in kwargs:
                print "--> panel "+str(i+1)+" of "+str(ncols*nrows)+": colorbar"
                ax1 = main_fig.add_axes([0.05+(i%ncols)*0.9/ncols+0.05/ncols, 0.95-np.ceil((i+1)/ncols)*0.9/nrows+0.5*0.9/nrows, 0.9*0.9/ncols, colorbar_width*0.9/nrows])
                if 'stretch' in kwargs:
                    if (kwargs['stretch'] == 'linear'):
                        colorbar = mpl.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=mpl.colors.Normalize(vmin=kwargs['vmin'], vmax=kwargs['vmax']), orientation='horizontal')
                        colorbar.set_label(kwargs['colorbar_label'])
                    elif (kwargs['stretch'] == 'log'):
                        log_ticks = [float('{:.2f}'.format(x)) for x in np.logspace(np.log10(kwargs['vmin']),np.log10(kwargs['vmax']),num=5, endpoint=True)]
                        colorbar = mpl.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=mpl.colors.LogNorm(vmin=kwargs['vmin'], vmax=kwargs['vmax']), ticks=log_ticks, orientation='horizontal')
                        colorbar.set_label(kwargs['colorbar_label'])
                        colorbar.set_ticks(log_ticks)
                        colorbar.set_ticklabels(['{:.2f}'.format(x) for x in log_ticks])
                    else:
                        print "--> only linear and log stretch supported!"
                else:
                    colorbar = mpl.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=mpl.colors.Normalize(vmin=kwargs['vmin'], vmax=kwargs['vmax']), orientation='horizontal')
                    colorbar.set_label(kwargs['colorbar_label'])
            else:
                print "--> you need to specify colorbar_cmap and colorbar_label to plot a colorbar"
            
        # add axis label and scale bar if bottom left plot
        if ( i == (nrows-1)*ncols ):
            fig.axis_labels.show()
            fig.tick_labels.show()
            fig.tick_labels.set_xformat(tick_label_xformat)
            fig.tick_labels.set_yformat(tick_label_yformat)
            fig.ticks.show()
            fig.ticks.set_xspacing(ticks_xspacing.to(u.degree).value)
            fig.ticks.set_yspacing(ticks_yspacing.to(u.degree).value)
            fig.ticks.set_minor_frequency(ticks_minor_frequency)
            fig.ticks.set_color(ticks_color)

#           fig.remove_scalebar()
            if 'scalebar_length' and 'scalebar_label' and 'scalebar_corner' in kwargs:
                fig.add_scalebar(length=kwargs['scalebar_length'].to(u.degree).value, label=kwargs['scalebar_label'], corner=kwargs['scalebar_corner'], frame=scalebar_frame)
                fig.scalebar.set_font(size=scalebar_fontsize)
                fig.scalebar.set_linestyle(scalebar_linestyle)
                fig.scalebar.set_linewidth(scalebar_linewidth)
                fig.scalebar.set_color(scalebar_color)
    
    if 'out' in kwargs:
        fig.save(kwargs['out'], dpi=300, transparent=True)
        print "--> saved file as "+kwargs['out']
    else:
        fig.save(os.path.splitext(fitscube)[0]+'.png', dpi=300, transparent=True)
        print "--> saved file as "+os.path.splitext(fitscube)[0]+'.png'


###################################################################################################

def aplpy_plot_pv(fitspv, **kwargs):
    
    """
    aplpy_plot_pv: position-velocity slice plotting
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Mandatory unnamed arguments:
        fitspv      Path and file name of the fits image to be plotted. This must
                    have one velocity axis and one spatial (offset) axis, such as
                    generated by CASA impv().
        

    Optional arguments:
        out         Path and file name of the created plot.
                    If not specified, the plot will be saved where the input image
                    is located. Default format: png
        
        figsize     Fiugre size as tpuel in inches. No default.
        
        xlabel      Label of the xaxis. As the PV slice can be calculated in 
                    arbitrary direction, this can be just an offset or along a
                    coordniate axis (e.g. RA).
        ylabel      Label of the yaxis. Both labels or none must be given.
                    If neither xlabel, nor ylabel is specified, the header 
                    information is used.
        
        cmap        Colormap to plot the image.
                    If not specified, the matplotlib default will be used, usually
                    this is viridis.
                    Every named matplotlib colormap can be used or any matplotlib 
                    colormap object.
                    For grayscale use cmap='grayscale'
                    
        vmin        Minimum value for colormap normalization.
        vmax        Maximum value for colormap normalization.
        
        stretch     Colormap strech, e.g. 'linear', 'log', 'sqrt', 'arcsinh'.
                    Linear and log scaling are fully implemented, other scaling can
                    cause errors when unusual argument combinations are chosen.
        
        recenter    Center the image on this location and set image width/height.
                    You either have to specify radius or width and height.
                    Must be an array containing an astropy.SkyCoord object plus one 
                    or two angular distances: [SkyCoord(...), 1*u.arcmin]
                    
        contour     List of contour elements.
                    Each contour element must be a list of 'image file', list of 
                    contour levels and list of colors for each contour. If only one
                    color is given, it will be used for all contours.
                    
        colorbar_location   As the name says.
                            Can be, e.g. 'bottom', 'right', ...
        colorbar_label      Can be specified only when colorbar_location is given.
                            String containing the label.
    
    
    General style settings
        Settings that do not have to be changed for each plot but maybe once per 
        script or once per project. Often used ones are tick_label_xformat, 
        ticks_xspacing and the corresponding settings for y.
            Can be accessed via
        import aplpy_plotting as ap
        ap.setting = ...
        See the aplpy_plotting.py main file the exact setting names if you need to 
        change them.


    example:
    
    aplpy_plot_pv('abc.fits', 
                  figsize = (8.27, 11.69),
                  vmin    = 0,
                  vmax    = 100,
                  stretch = 'linear',
                  cmap    = ap.viridis_cropped,
                  contour = [['xyz.fits', [1,2,3], 'black']], 
                  colorbar_location = 'right', 
                  colorbar_label = 'intensity [Jy\,beam$^{-1}$]', 
                  xlabel = 'offset [arcsec]', 
                  ylabel = 'velocity [km\,s$^{-1}$', 
                  out = 'abc.png'
                  )
    """
    
    print "--> plotting map "+fitspv
    
    if 'figsize' in kwargs:
        fig = aplpy.FITSFigure(fitspv, figsize=kwargs['figsize'])
    else:
        fig = aplpy.FITSFigure(fitspv)
    
    if 'vmin' and 'vmax' in kwargs:
        if 'stretch' in kwargs:
            if 'cmap' in kwargs:
                if kwargs['cmap'] == 'grayscale':
                    fig.show_grayscale(vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'], aspect='auto')
                else:
                    fig.show_colorscale(cmap=kwargs['cmap'], vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'], aspect='auto')
            else:
                fig.show_colorscale(vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'], aspect='auto')
        else:
            fig.show_colorscale(vmin=kwargs['vmax'], vmax=kwargs['vmax'], aspect='auto')
    else:
        fig.show_colorscale(aspect='auto')
        
    # recenter image
    if 'recenter' in kwargs:
        print "--> recenter not implemented yet"
    
    # contours?
    if 'contour' in kwargs:
        for cont_i in np.arange(len(kwargs['contour'])):
            if len(kwargs['contour'][cont_i]) == 3:
                fig.show_contour(data=kwargs['contour'][cont_i][0], levels=kwargs['contour'][cont_i][1], colors=kwargs['contour'][cont_i][2])
            else:
                print "--> wrong number or format of contour parameters in image "+str(cont_i)+". not plotting contours"

    # colorbar settings
    if 'colorbar_location' in kwargs:
        fig.add_colorbar()
        fig.colorbar.show()
        fig.colorbar.set_location(kwargs['colorbar_location'])
        if 'colorbar_label' in kwargs:
            fig.colorbar.set_axis_label_text(kwargs['colorbar_label'])
        if 'stretch' in kwargs:
            if (kwargs['stretch'] == 'log'):
                log_ticks = [float('{:.2f}'.format(round(x,int(-1*np.log10(kwargs['vmin']))))) for x in np.logspace(np.log10(kwargs['vmin']),np.log10(kwargs['vmax']),num=10, endpoint=True)]
                fig.colorbar.set_ticks(log_ticks)

    # scale bar
    # not possible with CTYPE='OFFSET'

    # ticks + labels
    fig.tick_labels.show()
#   fig.tick_labels.set_xformat(tick_label_xformat_pv)
#   fig.tick_labels.set_yformat(tick_label_yformat_pv)
    fig.ticks.show()
#   fig.ticks.set_xspacing(ticks_xspacing.to(u.degree).value)
#   fig.ticks.set_yspacing(ticks_yspacing.to(u.degree).value)
    fig.ticks.set_minor_frequency(ticks_minor_frequency)
    fig.ticks.set_color(ticks_color)
    
    # axis labels
    if 'xlabel' and 'ylabel' in kwargs:
        fig.set_axis_labels(kwargs['xlabel'],kwargs['ylabel'])
    else:
        print "--> you need to give both labels"

    if 'out' in kwargs:
        fig.save(kwargs['out'], dpi=300, transparent=True)
        print "--> saved file as "+kwargs['out']
    else:
        fig.save(os.path.splitext(fitspv)[0]+'.png', dpi=300, transparent=True)
        print "--> saved plot as "+os.path.splitext(fitspv)[0]+'.png'



###################################################################################################

def aplpy_map_grid(fitsimages, ncols, nrows, **kwargs):
    
    """
    aplpy_map_grid: grid plotting function
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Plot a grid of images. Each panel will have an overlay with the file name
    (can be changed with label_text). Axis labels are plotted on the bottom left 
    panel only. The same applies to beam and scalebar (if given). The last panel 
    always contains a horizontal colorbar. Each image can have its own contour.
    

    Mandatory unnamed, ordered arguments:
        fitsfile    List of fits image to be plotted.
        
        ncols       number of columns of the grid
        nrows       number of rows of the grid
        

    Optional arguments:
        out         Path and file name of the created plot.
                    If not specified, the plot will be saved where the input image
                    is located. Default format: png
                    
        figsize     Fiugre size as tpuel in inches. Default: (8.27, 11.69) A4 size
        
        cmap        Colormap to plot the image.
                    If not specified, the matplotlib default will be used, usually
                    this is viridis.
                    Every named matplotlib colormap can be used or any matplotlib 
                    colormap object.
                    For grayscale use cmap='grayscale'
                    
        vmin        Minimum value for colormap normalization.
        vmax        Maximum value for colormap normalization.
        
        stretch     Colormap strech, e.g. 'linear', 'log', 'sqrt', 'arcsinh'.
                    Linear and log scaling are fully implemented, other scaling can
                    cause errors when unusual argument combinations are chosen.
        
        recenter    Center the image on this location and set image width/height.
                    You either have to specify radius or width and height.
                    Must be an array containing an astropy.SkyCoord object plus one 
                    or two angular distances: [SkyCoord(...), 1*u.arcmin]
                    
        contour     List of contour elements in the same order as the input images.
                    The first contour element will be overplotted the first image, 
                    the second contour element over the second image and so on.
                    Each contour element must be a list of image file, list of 
                    contour levels and list of colors for each contour. If only one
                    color is given, it will be used for all contours.
                    If you want to have the same contour on each image, you have 
                    specify the contour as often as there are input image. (This will
                    probably change in the near future that you can give just one 
                    contour to be plotted on all images.)
                    
        label_text          List of labels to label to individual panels. Default is
                            to get a label from the file name without extension and
                            path.
        label_kwargs        Keyword arguments to format the label. All plt.text
                            kwargs are allowed, especially bbox=props can be used.
                            See example below.
                    
        colorbar_location   As the name says.
                            Can be, e.g. 'bottom', 'right', ...
        colorbar_cmap       The colorbar is generated via matplotlib and thus
                            specified manually at the moment. Make sure you use
                            the same colorbar as in the cmap argument.
                            Will be updated in a way that cmap is enough and this 
                            argument will be not necessary anymore.
        colorbar_label      Can be specified only when colorbar_location is given.
                            String containing the label.
                            
        scalebar_length     As the name says. Must be an angular astropy.units 
                            object, e.g. 10.0*u.arcmin
        scalebar_label      String containing the label to be plotted below the 
                            scalebar.
        scalebar_corner     Where should the scalebar and label be plotted?
                            E.g. 'bottom left'
        
        beam_corner Plot a beam ellipse as given in the fits header in this corner 
                    of the plot. Does not have to be a corner. 'bottom' also works.
                    Beams are plotted in each panel as they can be different for 
                    different images.
                    
        overlay     List of overlay elements.
                    Each overlay element is a list of shape (plt.scatter type), 
                    astropy.SkyCoord object specifying the position, size as 
                    angular astropy.unit and further plt.scatter kwargs (can be 
                    empty).
    
    
    General style settings
        Settings that do not have to be changed for each plot but maybe once per 
        script or once per project. Often used ones are tick_label_xformat, 
        ticks_xspacing and the corresponding settings for y.
            Can be accessed via
        import aplpy_plotting as ap
        ap.setting = ...
        See the aplpy_plotting.py main file the exact setting names if you need to 
        change them.


    example:
    
    props = {'boxstyle': "round", 'facecolor': "w", 'edgecolor': "black", 'linewidth': 0.5, 'alpha': 0.8}
    aplpy_map_grid([image1.fits,image2.fits], 
                   ncols, 
                   nrows, 
                   figsize  = (8.27, 11.69), 
                   cmap     = 'jet', 
                   vmin     = -0.05, 
                   vmax     = 3.5, 
                   stretch  = 'linear', 
                   recenter = SkyCoord('01h23m45.6s 12d34m45.6s'), 
                   contour  = [['contour1.fits', [1,2,3], ['white', 'grey', 'black']],['contour2.fits',[9,8,7],'blue']], 
                   label_text      = [image1_label, image2_label], 
                   label_kwargs    = {'bbox': props}, 
                   overlay         = [['circle', SkyCoord('01h23m45.6s 12d34m45.6s'), 1.0*u.arcmin, {'linewidth': 1.0}]],
                   beam_corner     = 'bottom left', 
                   colorbar_cmap   = mpl.cm.jet, 
                   colorbar_label  = 'flux density [Jy/beam]', 
                   scalebar_length = 1.0, 
                   scalebar_label  = 'string', 
                   scalebar_corner = 'bottom'
                   )
    """
    
    print "--> plotting map grid of these maps: ", fitsimages
    
    if 'figsize' in kwargs:
        main_figsize = kwargs['figsize']
    else:
        main_figsize = (8.27, 11.69)    # A4 in inches
    
    main_fig = plt.figure(figsize = main_figsize)
    for i in np.arange(len(fitsimages)):
        
        # get subplot specific info
        subplt_size = [0.05+(i%ncols)*0.9/ncols, 0.95-np.ceil((i+1)/ncols)*0.9/nrows, 0.9/ncols, 0.9/nrows]
        
        # plot channel map if not last panel
        if ( i < nrows*ncols-1 ):
            print "--> panel "+str(i+1)+" of "+str(nrows*ncols)
            
            fig = aplpy.FITSFigure(fitsimages[i], figure=main_fig, subplot=subplt_size)
            if 'vmin' and 'vmax' in kwargs:
                if 'stretch' in kwargs:
                    if 'cmap' in kwargs:
                        if kwargs['cmap'] == 'grayscale':
                            fig.show_grayscale(vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'])
                        else:
                            fig.show_colorscale(cmap=kwargs['cmap'], vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'])
                    else:
                        fig.show_colorscale(vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'])
                else:
                    fig.show_colorscale(vmin=kwargs['vmax'], vmax=kwargs['vmax'])
            else:
                fig.show_colorscale()
            
            # force subfigure to be only as wide as the image
#           fig.image.axes.set_adjustable('box-forced')
            
            # recenter image
            if 'recenter' in kwargs:
                if (len(kwargs['recenter']) == 2):
                    fig.recenter(kwargs['recenter'][0].ra.degree, kwargs['recenter'][0].dec.degree, radius=kwargs['recenter'][1].to(u.degree).value)
                elif (len(kwargs['recenter']) == 3):
                    fig.recenter(kwargs['recenter'][0].ra.degree, kwargs['recenter'][0].dec.degree, width=kwargs['recenter'][1].to(u.degree).value, height=kwargs['recenter'][2].to(u.degree).value)
                else:
                    print "--> specify SkyCoord(x,y) and either radius or width, height. not recentering"
            
            # contours?
            if 'contour' in kwargs:
                for cont_i in np.arange(len(kwargs['contour'])):
                    if len(kwargs['contour'][cont_i]) == 3:
                        fig.show_contour(data=kwargs['contour'][cont_i][0], levels=kwargs['contour'][cont_i][1], colors=kwargs['contour'][cont_i][2])
                    else:
                        print "--> wrong number or format of contour parameters in image "+str(cont_i)+". not plotting contours"
            
            # ticks + labels
            fig.axis_labels.hide()
            fig.tick_labels.hide()
            fig.ticks.show()
            fig.ticks.set_xspacing(ticks_xspacing.to(u.degree).value)
            fig.ticks.set_yspacing(ticks_yspacing.to(u.degree).value)
            fig.ticks.set_minor_frequency(ticks_minor_frequency)
            fig.ticks.set_color(ticks_color)
            
            # data set overlay
            if 'label_text' in kwargs:
                fig.add_label(0.5, 0.9, kwargs['label_text'][i].replace('_','$\_$'), color='black', relative=True, size=velo_fontsize)
                if 'label_kwargs' in kwargs:
                    fig.add_label(0.5, 0.9, kwargs['label_text'][i].replace('_','$\_$'), color='black', relative=True, size=velo_fontsize, **kwargs['label_kwargs'])
            else:
                fig.add_label(0.5, 0.9, fitsimages[i][:-5].replace('_','$\_$'), color='black', relative=True, size=velo_fontsize)
            
            if 'overlay' in kwargs:
                for olay in np.arange(len(kwargs['overlay'])):
                    if (kwargs['overlay'][olay][0] == 'circle'):
                        fig.show_circles(xw=kwargs['overlay'][olay][1].ra.degree, yw=kwargs['overlay'][olay][1].dec.degree, radius=kwargs['overlay'][olay][2].to(u.degree).value, **kwargs['overlay'][olay][3])
                    else:
                        fig.show_markers(xw=kwargs['overlay'][olay][1].ra.degree, yw=kwargs['overlay'][olay][1].dec.degree, marker=kwargs['overlay'][olay][0], s=kwargs['overlay'][olay][2], **kwargs['overlay'][olay][3])
            
            # beam settings
            if 'beam_corner' in kwargs:
                fig.add_beam()
                fig.beam.show()
                fig.beam.set_corner(kwargs['beam_corner'])
                fig.beam.set_frame(beam_frame)
                fig.beam.set_color(beam_color)
            
        # add axis label and scale bar if bottom left plot
        if ( i == (nrows-1)*ncols ):
            fig.axis_labels.show()
            fig.tick_labels.show()
            fig.tick_labels.set_xformat(tick_label_xformat)
            fig.tick_labels.set_yformat(tick_label_yformat)
            fig.ticks.show()
            fig.ticks.set_xspacing(ticks_xspacing.to(u.degree).value)
            fig.ticks.set_yspacing(ticks_yspacing.to(u.degree).value)
            fig.ticks.set_minor_frequency(ticks_minor_frequency)
            fig.ticks.set_color(ticks_color)

    #       fig.remove_scalebar()
            if 'scalebar_length' and 'scalebar_label' and 'scalebar_corner' in kwargs:
                fig.add_scalebar(length=kwargs['scalebar_length'].to(u.degree).value, label=kwargs['scalebar_label'], corner=kwargs['scalebar_corner'], frame=scalebar_frame)
                fig.scalebar.set_font(size=scalebar_fontsize)
                fig.scalebar.set_linestyle(scalebar_linestyle)
                fig.scalebar.set_linewidth(scalebar_linewidth)
                fig.scalebar.set_color(scalebar_color)
        
        # colorbar settings
        if 'colorbar_cmap' and 'colorbar_label' in kwargs:
            ax1 = main_fig.add_axes([0.05+(ncols-1+0.05)*0.9/ncols, 0.05+0.5*0.9/nrows, 0.9*0.9/ncols, colorbar_width*0.9/nrows])
            if 'vmin' and 'vmax' in kwargs:
                colorbar = mpl.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=mpl.colors.Normalize(vmin=kwargs['vmin'], vmax=kwargs['vmax']), orientation='horizontal')
            else:
                colorbar = mpl.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=mpl.colors.Normalize(vmin=0.0, vmax=1.0), orientation='horizontal')
            colorbar.ax.tick_params(labelsize = colorbar_fontsize)
            colorbar.ax.set_xticklabels([item.get_text() for item in colorbar.ax.get_xticklabels()], rotation=90)
            colorbar.set_label(kwargs['colorbar_label'], size = colorbar_fontsize)
        else:
            print "--> you need to define both colorbar_location and colorbar_label to plot a colorbar"
    
    if 'out' in kwargs:
        fig.save(kwargs['out'], dpi=300, transparent=True)
        print "--> saved file as "+kwargs['out']
    else:
        fig.save('grid_plot.png', dpi=300, transparent=True)
        print "--> saved file as grid_plot.png"
