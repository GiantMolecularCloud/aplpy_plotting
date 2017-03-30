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

import aplpy_plotting as ap

import os as __os__
import aplpy as __aplpy__
import numpy as __np__
from astropy import units as __u__
import matplotlib as __mpl__
import matplotlib.pyplot as __plt__
from matplotlib import rc as __rc__
__rc__('text',usetex=True)
from matplotlib.cbook import MatplotlibDeprecationWarning as __MatplotlibDeprecationWarning__
import warnings as __warnings__
__warnings__.simplefilter('ignore', __MatplotlibDeprecationWarning__)

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
        See the aplpy_plotting.py main file for the exact setting names if you need
        to change them.


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
    
    print("--> plotting channel map "+fitscube)
    if not ('chan_width' and 'chan_vel0' in kwargs):
        print("--> channel width and velocity of first channel not given. will not calculate velocity")
    
    if 'figsize' in kwargs:
        main_figsize = kwargs['figsize']
    else:
        main_figsize = (8.27, 11.69)    # A4 in inches
    
    main_fig = __plt__.figure(figsize=main_figsize)

    # convert to float for python 2 compatibility
    ncols_f = float(ncols)
    nrows_f = flaot(nrows)

    for i in __np__.arange(nrows*ncols):
        
        # get subplot specific info
        subplt_size = [0.05+(i%ncols_f)*0.9/ncols_f, 0.95-__np__.ceil((i+1)/ncols_f)*0.9/nrows_f, 0.9/ncols_f, 0.9/nrows_f]
        #chn_slice   = [chan_start+i*chan_iter,0]   # if 4th axis is present
        chn_slice   = [chan_start+i*chan_iter]  # if no 4th axis available

        # introduce a check here to assure that the cube has only three axis (no forth degenerate axis, like Stokes)
        
        if 'chan_width' and 'chan_vel0' in kwargs:
            chn_velo    = (i*chan_iter+chan_start)*kwargs['chan_width']+kwargs['chan_vel0']
        
        # plot channel map if not last panel
        if ( i < nrows*ncols-1 ):
            print("--> panel "+str(i+1)+" of "+str(nrows*ncols))
            
            fig = __aplpy__.FITSFigure(fitscube, figure=main_fig, subplot=subplt_size, dimensions=[0,1], slices=chn_slice)
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
                    fig.recenter(kwargs['recenter'][0].ra.degree, kwargs['recenter'][0].dec.degree, radius=kwargs['recenter'][1].to(__u__.degree).value)
                elif (len(kwargs['recenter']) == 3):
                    fig.recenter(kwargs['recenter'][0].ra.degree, kwargs['recenter'][0].dec.degree, width=kwargs['recenter'][1].to(__u__.degree).value, height=kwargs['recenter'][2].to(__u__.degree).value)
                else:
                    print("--> specify SkyCoord(x,y) and either radius or width, height. not recentering")
            
            # contours?
            if 'contour' in kwargs:
                if len(kwargs['contour']) == 4:
                    if kwargs['contour'][0] == 'pixelmask':
                        fig.show_contour(data=kwargs['contour'][1], levels=kwargs['contour'][2], colors=kwargs['contour'][3])
                    elif kwargs['contour'][0] == 'cubemask':
                        fig.show_contour(data=kwargs['contour'][1], dimensions=[0,1], slices=chn_slice, levels=kwargs['contour'][2], colors=kwargs['contour'][3])
                    else:
                        print("--> wrong mask type specification")
                else:
                    print("--> wrong number or format of contour parameters. not plotting contours")
            
            
            # ticks + labels
            fig.axis_labels.hide()
            fig.tick_labels.hide()
            fig.ticks.show()
            fig.ticks.set_xspacing(ap.ticks_xspacing.to(__u__.degree).value)
            fig.ticks.set_yspacing(ap.ticks_yspacing.to(__u__.degree).value)
            fig.ticks.set_minor_frequency(ap.ticks_minor_frequency)
            fig.ticks.set_color(ap.ticks_color)
            
            # velocity or channel number overlay
            if 'chan_width' and 'chan_vel0' in kwargs:
                fig.add_label(0.8, 0.8, str('{:3.1f}'.format(chn_velo))+r'\,km\,s$^{-1}$', color='black', relative=True, size=ap._velo_fontsize)
            else:
                fig.add_label(0.8, 0.8, 'channel '+str(chn_slice), color='black', relative=True, size=ap._velo_fontsize)
            
            # beam settings
            if 'beam_corner' in kwargs:
                fig.add_beam()
                fig.beam.show()
                fig.beam.set_corner(kwargs['beam_corner'])
                fig.beam.set_frame(ap._beam_frame)
                fig.beam.set_color(ap._beam_color)
            
        # colorbar settings
        if (i == ncols*nrows-1):
            if 'colorbar_cmap' and 'colorbar_label' in kwargs:
                print("--> panel "+str(i+1)+" of "+str(ncols*nrows)+": colorbar")
                ax1 = main_fig.add_axes([0.05+(i%ncols_f)*0.9/ncols_f+0.05/ncols_f, 0.95-__np__.ceil((i+1)/ncols_f)*0.9/nrows_f+0.5*0.9/nrows_f, 0.9*0.9/ncols_f, colorbar_width*0.9/nrows_f])
                if 'stretch' in kwargs:
                    if (kwargs['stretch'] == 'linear'):
                        colorbar = __mpl__.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=__mpl__.colors.Normalize(vmin=kwargs['vmin'], vmax=kwargs['vmax']), orientation='horizontal')
                        colorbar.set_label(kwargs['colorbar_label'])
                    elif (kwargs['stretch'] == 'log'):
                        log_ticks = [float('{:.2f}'.format(x)) for x in __np__.logspace(__np__.log10(kwargs['vmin']),__np__.log10(kwargs['vmax']),num=5, endpoint=True)]
                        colorbar = __mpl__.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=__mpl__.colors.LogNorm(vmin=kwargs['vmin'], vmax=kwargs['vmax']), ticks=log_ticks, orientation='horizontal')
                        colorbar.set_label(kwargs['colorbar_label'])
                        colorbar.set_ticks(log_ticks)
                        colorbar.set_ticklabels(['{:.2f}'.format(x) for x in log_ticks])
                    else:
                        print("--> only linear and log stretch supported!")
                else:
                    colorbar = __mpl__.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=__mpl__.colors.Normalize(vmin=kwargs['vmin'], vmax=kwargs['vmax']), orientation='horizontal')
                    colorbar.ax.tick_params(labelsize = ap._colorbar_fontsize)
                    colorbar.set_label(kwargs['colorbar_label'], size = ap._colorbar_fontsize)
            else:
                print("--> you need to specify colorbar_cmap and colorbar_label to plot a colorbar")
            
        # add axis label and scale bar if bottom left plot
        if ( i == (nrows-1)*ncols ):
            fig.axis_labels.show()
            fig.tick_labels.show()
            fig.tick_labels.set_xformat(ap.tick_label_xformat)
            fig.tick_labels.set_yformat(ap.tick_label_yformat)
            fig.ticks.show()
            fig.ticks.set_xspacing(ap.ticks_xspacing.to(__u__.degree).value)
            fig.ticks.set_yspacing(ap.ticks_yspacing.to(__u__.degree).value)
            fig.ticks.set_minor_frequency(ap.ticks_minor_frequency)
            fig.ticks.set_color(ap._ticks_color)

#           fig.remove_scalebar()
            if 'scalebar_length' and 'scalebar_label' and 'scalebar_corner' in kwargs:
                fig.add_scalebar(length=kwargs['scalebar_length'].to(__u__.degree).value, label=kwargs['scalebar_label'], corner=kwargs['scalebar_corner'], frame=_scalebar_frame)
                fig.scalebar.set_font(size=ap._scalebar_fontsize)
                fig.scalebar.set_linestyle(ap._scalebar_linestyle)
                fig.scalebar.set_linewidth(ap._scalebar_linewidth)
                fig.scalebar.set_color(ap._scalebar_color)
    
    if 'out' in kwargs:
        fig.save(kwargs['out'], dpi=300, transparent=True)
        print("--> saved file as "+kwargs['out'])
    else:
        fig.save(__os__.path.splitext(fitscube)[0]+'.png', dpi=300, transparent=True)
        print("--> saved file as "+__os__.path.splitext(fitscube)[0]+'.png')


###################################################################################################
