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
from matplotlib import rc as __rc__
__rc__('text',usetex=True)
from matplotlib.cbook import MatplotlibDeprecationWarning as __MatplotlibDeprecationWarning__
import warnings as __warnings__
__warnings__.simplefilter('ignore', __MatplotlibDeprecationWarning__)

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

        figsize     Fiugre size as tupel in inches. No default.

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
        legend      Add a legend that lists the contour objects. The default is a
                    legend with two columns over the image listing the countour file
                    names.
        legend_kwargs       Specify the legend details as in ax.legend().

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

        execute_code        This option allows to pass arbitrary code that is executed
                            just before saving the figure and can be used to access
                            aplpy functionality that is not mapped by aplpy_plotting.
                            The code must be given in a list of strings. Note that the
                            correct namespaces need to be given, e.g. access numpy with
                            __np__ instead of np. The figure objects are called fig, or
                            main_fig in plots with multiple figures. Example:
                            execute_code = ["fig.show_lines([0,10],[2,5],color='k'"]

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


    print("--> plotting map "+fitsfile)

    if 'figsize' in kwargs:
        main_figsize = kwargs['figsize']
        fig = __aplpy__.FITSfigure(fitsfile, figsize=main_figsize)
    else:
        fig = __aplpy__.FITSFigure(fitsfile)

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
            if 'vmin' and 'vmax' in kwargs:
                fig.show_colorscale(cmap=kwargs['cmap'], vmin=kwargs['vmin'], vmax=kwargs['vmax'])
                if 'stretch' in kwargs:
                    fig.show_colorscale(cmap=kwargs['cmap'], stretch=kwargs['stretch'], vmin=kwargs['vmin'], vmax=kwargs['vmax'])
                    if 'invert' in kwargs:
                        fig.show_colorscale(cmap=kwargs['cmap'], vmin=kwargs['vmin'], vmax=kwargs['vmax'], stretch=kwargs['stretch'], invert=kwargs['invert'])
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

    # contours
    if 'contour' in kwargs:
        contournum = 0  # counting variable for number of total contours plotted
        for cont_i in __np__.arange(len(kwargs['contour'])):
            if len(kwargs['contour'][cont_i]) == 3:
                fig.show_contour(data=kwargs['contour'][cont_i][0], levels=kwargs['contour'][cont_i][1], colors=kwargs['contour'][cont_i][2])
                if 'legend' in kwargs:
                    if (kwargs['legend'] == True):
                        fig._ax1.collections[contournum].set_label(kwargs['contour'][cont_i][0].replace('_','$\_$'))
                    elif (isinstance(kwargs['legend'], (list,tuple))):
                        fig._ax1.collections[contournum].set_label(kwargs['legend'][cont_i])
                    else:
                        print("--> wrong format for legend: either True or list of names for each contour.")
                    contournum += len(kwargs['contour'][cont_i][1])     # count up plotted contours
            elif len(kwargs['contour'][cont_i]) == 4:
                fig.show_contour(data=kwargs['contour'][cont_i][0], slices=[kwargs['contour'][1]], levels=kwargs['contour'][cont_i][2], colors=kwargs['contour'][cont_i][3])
                if 'legend' in kwargs:
                    if (kwargs['legend'] == True):
                        fig._ax1.collections[contournum].set_label(kwargs['contour'][cont_i][0].replace('_','$\_$'))
                    elif (isinstance(kwargs['legend'], (list,tuple))):
                        fig._ax1.collections[contournum].set_label(kwargs['legend'][cont_i])
                    else:
                        print("--> wrong format for legend: either True or list of names for each contour.")
                    contournum += len(kwargs['contour'][cont_i][1])     # count up plotted contours
            else:
                print("--> wrong number or format of contour parameters in image "+str(cont_i)+". not plotting contours")
            if 'clabel' in kwargs:
                if kwargs['clabel'] == True:
                    fig._layers['contour_set_'+str(cont_i+1)].clabel()
                if isinstance(kwargs['clabel'],dict):
                    fig._layers['contour_set_'+str(cont_i+1)].clabel(**kwargs['clabel'])

    # regions
    if 'regions' in kwargs:
        # add a test if pyregions is installed
        for reg_i in __np__.arange(len(kwargs['regions'])):
            fig.show_regions(kwargs['regions'][reg_i])

    # colorbar settings
    if 'colorbar_location' in kwargs:
        fig.add_colorbar()
        fig.colorbar.show()
        fig.colorbar.set_location(kwargs['colorbar_location'])
        if 'colorbar_label' in kwargs:
            fig.colorbar.set_axis_label_text(kwargs['colorbar_label'])
        if 'stretch' in kwargs:
            if (kwargs['stretch'] == 'log'):
                log_ticks = [float('{:.2f}'.format(round(x,int(-1*__np__.log10(kwargs['vmin']))))) for x in __np__.logspace(__np__.log10(kwargs['vmin']),__np__.log10(kwargs['vmax']),num=10, endpoint=True)]
                fig.colorbar.set_ticks(log_ticks)
        fig.colorbar.set_font(size=ap._colorbar_fontsize)
        fig.colorbar.set_axis_label_font(size=ap._colorbar_fontsize)
        fig.colorbar.set_frame_color(ap._frame_color)

    # scale bar
    if 'scalebar_length' and 'scalebar_label' and 'scalebar_corner' in kwargs:
        fig.add_scalebar(length=kwargs['scalebar_length'].to(__u__.degree).value, label=kwargs['scalebar_label'], corner=kwargs['scalebar_corner'], frame=ap._scalebar_frame)
        fig.scalebar.set_font(size=ap._scalebar_fontsize)
        fig.scalebar.set_linestyle(ap._scalebar_linestyle)
        fig.scalebar.set_linewidth(ap._scalebar_linewidth)
        fig.scalebar.set_color(ap._scalebar_color)

    # beam settings
    if 'beam_corner' in kwargs:
        fig.add_beam()
        fig.beam.show()
        fig.beam.set_corner(kwargs['beam_corner'])
        fig.beam.set_frame(ap._beam_frame)
        fig.beam.set_color(ap._beam_color)

    # data set overlay
    if 'label_text' in kwargs:
        if isinstance(kwargs['label_text'], str):
            if 'label_kwargs' in kwargs:
                fig.add_label(0.5, 0.9, kwargs['label_text'].replace('_','$\_$'), color='black', relative=True, size=ap._velo_fontsize, **kwargs['label_kwargs'])
            else:
                fig.add_label(0.5, 0.9, kwargs['label_text'].replace('_','$\_$'), color='black', relative=True, size=ap._velo_fontsize)
        elif isinstance(kwargs['label_text'], list):
            if (len(kwargs['label_text']) == 3):
                if 'label_kwargs' in kwargs:
                    fig.add_label(kwargs['label_text'][0], kwargs['label_text'][1], kwargs['label_text'][2].replace('_','$\_$'), color='black', relative=True, size=ap._velo_fontsize, **kwargs['label_kwargs'])
                else:
                    fig.add_label(kwargs['label_text'][0], kwargs['label_text'][1], kwargs['label_text'][2].replace('_','$\_$'), color='black', relative=True, size=ap._velo_fontsize)
            else:
                raise ValueError('Need xpos, ypos and string to plot label.')
        else:
            raise TypeError('Unknown type for label_text. Either string or list with three entries are supported.')

    if 'overlay' in kwargs:
        for olay in __np__.arange(len(kwargs['overlay'])):
            if (kwargs['overlay'][olay][0] == 'circle'):
                fig.show_circles(xw=kwargs['overlay'][olay][1].ra.degree, yw=kwargs['overlay'][olay][1].dec.degree, radius=kwargs['overlay'][olay][2].to(__u__.degree).value, **kwargs['overlay'][olay][3])
            else:
                fig.show_markers(xw=kwargs['overlay'][olay][1].ra.degree, yw=kwargs['overlay'][olay][1].dec.degree, marker=kwargs['overlay'][olay][0], s=kwargs['overlay'][olay][2], **kwargs['overlay'][olay][3])

    # ticks + labels
    fig.tick_labels.show()
    fig.tick_labels.set_xformat(ap.tick_label_xformat)
    fig.tick_labels.set_yformat(ap.tick_label_yformat)
    fig.tick_labels.set_font(size=ap._tick_label_fontsize)
    #fig.tick_labels.hide_x()
    #fig.tick_labels.hide_y()
    fig.ticks.show()
    fig.ticks.set_xspacing(ap.ticks_xspacing.to(__u__.degree).value)
    fig.ticks.set_yspacing(ap.ticks_yspacing.to(__u__.degree).value)
    fig.ticks.set_minor_frequency(ap.ticks_minor_frequency)
    fig.ticks.set_color(ap._ticks_color)
    fig.frame.set_color(ap._frame_color)
    fig.axis_labels.set_font(size=ap._tick_label_fontsize)

    # add legend
    adjbbox = True
    if 'legend' in kwargs:
        adjbbox = False
        if 'legend_kwargs' in kwargs:
            fig._ax1.legend(**kwargs['legend_kwargs'])
        else:
            fig._ax1.legend(loc=3, bbox_to_anchor=(0.,1.02,1.,0.1), ncol=2, mode='expand', borderaxespad=0., fontsize=ap._colorbar_fontsize)

    # execute additional code passed by the user
    if 'execute_code' in kwargs:
        if (isinstance(kwargs['execute_code'], (list,tuple))):
            for codes in kwargs['execute_code']:
                exec(codes)
        else:
            print("--> Code to execute must be given in a list of strings")

    if 'out' in kwargs:
        fig.save(kwargs['out'], dpi=300, transparent=True, adjust_bbox=adjbbox)
        print("--> saved file as "+kwargs['out'])
    else:
        fig.save(__os__.path.splitext(fitsfile)[0]+'.png', dpi=300, transparent=True, adjust_bbox=adjbbox)
        print("--> saved plot as "+__os__.path.splitext(fitsfile)[0]+'.png')



###################################################################################################
