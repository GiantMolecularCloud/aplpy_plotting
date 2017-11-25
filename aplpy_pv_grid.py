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

def aplpy_pv_grid(fitsimages, ncols, nrows, **kwargs):

    """
    aplpy_map_grid: grid plotting function
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Plot a grid of position-velocity diagrams. Each panel will have an overlay with
    the file name (can be changed with label_text). Axis labels are plotted on the
    bottom left panel only. The same applies to beam and scalebar (if given). The
    last panel always contains a horizontal colorbar. Each image can have its own
    contour.
    (Basically this is just a copy of aplpy_map_grid with adjustments for plotting
    position-velocity diagrams.)


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

        xlabel      Label of the xaxis. As the PV slice can be calculated in
                    arbitrary direction, this can be just an offset or along a
                    coordniate axis (e.g. RA).
        ylabel      Label of the yaxis. Both labels or none must be given.
                    If neither xlabel, nor ylabel is specified, the header
                    information is used.

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

        overlay     List of overlay elements.
                    Each overlay element is a list of shape (plt.scatter type),
                    astropy.SkyCoord object specifying the position, size as
                    angular astropy.unit and further plt.scatter kwargs (can be
                    empty).

        execute_code        This option allows to pass arbitrary code that is executed
                            for each panel just before moving on to the next panel. It
                            can be used to access aplpy functionality that is not
                            mapped by aplpy_plotting. The code must be given in a list
                            of strings. Note that the correct namespaces need to be
                            given, e.g. access numpy with __np__ instead of np. The
                            figure objects are called fig, or main_fig in plots with
                            multiple figures. Example:
                            execute_code = [["fig.show_lines([0,10],[2,5],color='k'"]]
        execute_once        Execute arbitrary code passed as list of strings. This is
                            executed just once before saving the figure. The intended
                            use is to create a single legend for all the panels. The
                            same comments as in execute_code apply.

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

    props = {'boxstyle': "round", 'facecolor': "w", 'edgecolor': "black", 'linewidth': 0.5, 'alpha': 0.8}
    aplpy_map_grid([image1.fits,image2.fits],
                   ncols,
                   nrows,
                   figsize  = (8.27, 11.69),
                   cmap     = 'jet',
                   vmin     = -0.05,
                   vmax     = 3.5,
                   stretch  = 'linear',
                   xlabel   = 'offset ["]',
                   ylabel   = 'v$_{rad}$ [km\,s^{-1}$]',
                   recenter = SkyCoord('01h23m45.6s 12d34m45.6s'),
                   contour  = [['contour1.fits', [1,2,3], ['white', 'grey', 'black']],['contour2.fits',[9,8,7],'blue']],
                   label_text      = [image1_label, image2_label],
                   label_kwargs    = {'bbox': props},
                   overlay         = [['circle', SkyCoord('01h23m45.6s 12d34m45.6s'), 1.0*u.arcmin, {'linewidth': 1.0}]],
                   colorbar_cmap   = mpl.cm.jet,
                   colorbar_label  = 'flux density [Jy/beam]',
                   scalebar_length = 1.0,
                   scalebar_label  = 'string',
                   scalebar_corner = 'bottom'
                   )
    """

    print("--> plotting map grid of these maps: ", fitsimages)

    if 'figsize' in kwargs:
        main_figsize = kwargs['figsize']
    else:
        main_figsize = (8.27, 11.69)    # A4 in inches

    main_fig = __plt__.figure(figsize = main_figsize)

    # convert to float for python 2 compatibility
    ncols_f = float(ncols)
    nrows_f = float(nrows)

    for i in __np__.arange(len(fitsimages)):

        # get subplot specific info
        subplt_size = [0.05+(i%ncols_f)*0.9/ncols_f, 0.95-__np__.ceil((i+1)/ncols_f)*0.9/nrows_f, 0.9/ncols_f, 0.9/nrows_f]

        # plot channel map if not last panel
        if ( i < nrows*ncols-1 ):
            print("--> panel "+str(i+1)+" of "+str(nrows*ncols))

            # add option to use single vmin/vmax or list (one for each panel)
            fig = __aplpy__.FITSFigure(fitsimages[i], figure=main_fig, subplot=subplt_size)
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

            # force subfigure to be only as wide as the image
#           fig.image.axes.set_adjustable('box-forced')

            # recenter image
            if 'recenter' in kwargs:
                if (len(kwargs['recenter']) == 4):
                    fig.recenter(kwargs['recenter'][0].value, kwargs['recenter'][1].value, width=kwargs['recenter'][2].value, height=kwargs['recenter'][3].value)
                else:
                    print("--> specify center (x), center (y), width and height as astropy.units object. not recentering")

            # contours?
            if 'contour' in kwargs:
                if len(kwargs['contour'][i]) == 3:
                    fig.show_contour(data=kwargs['contour'][i][0], levels=kwargs['contour'][i][1], colors=kwargs['contour'][i][2])
                elif len(kwargs['contour'][i]) == 4:
                    fig.show_contour(data=kwargs['contour'][i][0], slices=[kwargs['contour'][i][1]], levels=kwargs['contour'][i][2], colors=kwargs['contour'][i][3])
                else:
                    print("--> wrong number or format of contour parameters in image "+str(i)+". not plotting contours")

            # ticks + labels
            fig.axis_labels.hide()
            fig.tick_labels.hide()
            fig.ticks.show()
            fig.ticks.set_minor_frequency(ap.ticks_minor_frequency)
            fig.ticks.set_color(ap._ticks_color)
            fig.frame.set_color(ap._frame_color)

            # data set overlay
            if 'label_text' in kwargs:
                fig.add_label(0.5, 0.9, kwargs['label_text'][i].replace('_','$\_$'), color='black', relative=True, size=ap._velo_fontsize)
                if 'label_kwargs' in kwargs:
                    fig.add_label(0.5, 0.9, kwargs['label_text'][i].replace('_','$\_$'), color='black', relative=True, size=ap._velo_fontsize, **kwargs['label_kwargs'])
            else:
                fig.add_label(0.5, 0.9, fitsimages[i][:-5].replace('_','$\_$'), color='black', relative=True, size=ap._velo_fontsize)

            if 'overlay' in kwargs:
                for olay in __np__.arange(len(kwargs['overlay'])):
                    if (kwargs['overlay'][olay][0] == 'circle'):
                        fig.show_circles(xw=kwargs['overlay'][olay][1].ra.degree, yw=kwargs['overlay'][olay][1].dec.degree, radius=kwargs['overlay'][olay][2].to(__u__.degree).value, **kwargs['overlay'][olay][3])
                    else:
                        fig.show_markers(xw=kwargs['overlay'][olay][1].ra.degree, yw=kwargs['overlay'][olay][1].dec.degree, marker=kwargs['overlay'][olay][0], s=kwargs['overlay'][olay][2], **kwargs['overlay'][olay][3])

        # add axis label and scale bar if bottom left plot
        if ( i == (nrows-1)*ncols ):
            fig.axis_labels.show()
            fig.tick_labels.show()
            fig.tick_labels.set_font(size=ap._tick_label_fontsize)
            fig.ticks.set_minor_frequency(ap.ticks_minor_frequency)
            fig.axis_labels.set_font(size=ap._tick_label_fontsize)

            # axis labels
            if 'xlabel' and 'ylabel' in kwargs:
                fig.set_axis_labels(kwargs['xlabel'],kwargs['ylabel'])
            else:
                print("--> you need to give both labels")

    #       fig.remove_scalebar()
            if 'scalebar_length' and 'scalebar_label' and 'scalebar_corner' in kwargs:
                fig.add_scalebar(length=kwargs['scalebar_length'].to(__u__.degree).value, label=kwargs['scalebar_label'], corner=kwargs['scalebar_corner'], frame=ap._scalebar_frame)
                fig.scalebar.set_font(size=ap._scalebar_fontsize)
                fig.scalebar.set_linestyle(ap._scalebar_linestyle)
                fig.scalebar.set_linewidth(ap._scalebar_linewidth)
                fig.scalebar.set_color(ap._scalebar_color)

        # execute additional code passed by the user
        if 'execute_code' in kwargs:
            if (isinstance(kwargs['execute_code'], (list,tuple))):
                if (len(kwargs['execute_code']) == len(fitsimages)):
                    for codes in kwargs['execute_code'][i]:
                        exec(codes)
                else:
                    print("Please give exactly one list element for each panel. Can also be empty or list of multiple commands.")
            else:
                print("Code to execute must be given in a list of list of strings")

    # colorbar settings
    if 'colorbar_cmap' and 'colorbar_label' in kwargs:
        i = ncols*nrows-1
        print("--> panel "+str(i+1)+" of "+str(ncols*nrows)+": colorbar")
        ax1 = main_fig.add_axes([0.05+(ncols_f-1+0.05)*0.9/ncols_f, 0.05+0.5*0.9/nrows_f, 0.9*0.9/ncols_f, ap._colorbar_width*0.9/nrows_f])

        if 'stretch' in kwargs:
            if (kwargs['stretch'] == 'linear'):
                colorbar = __mpl__.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=__mpl__.colors.Normalize(vmin=kwargs['vmin'], vmax=kwargs['vmax']), orientation='horizontal')
            elif (kwargs['stretch'] == 'log'):
                log_ticks = [float('{:.2f}'.format(x)) for x in __np__.logspace(__np__.log10(kwargs['vmin']),__np__.log10(kwargs['vmax']),num=5, endpoint=True)]
                colorbar = __mpl__.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=__mpl__.colors.LogNorm(vmin=kwargs['vmin'], vmax=kwargs['vmax']), ticks=log_ticks, orientation='horizontal')
                colorbar.set_label(kwargs['colorbar_label'])
                colorbar.set_ticks(log_ticks)
                colorbar.set_ticklabels(['{:.2f}'.format(x) for x in log_ticks])
                colorbar.outline.set_edgecolor(ap._frame_color)
                #colorbar.dividers.set_color(ap._frame_color)
            else:
                print("--> only linear and log stretch supported!")
        else:
            colorbar = __mpl__.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=__mpl__.colors.Normalize(vmin=kwargs['vmin'], vmax=kwargs['vmax']), orientation='horizontal')
        colorbar.set_label(kwargs['colorbar_label'], size = ap._colorbar_fontsize)
        colorbar.ax.tick_params(labelsize = ap._colorbar_fontsize)
        #colorbar.ax.set_xticklabels([item.get_text() for item in colorbar.ax.get_xticklabels()], rotation=90)
        colorbar.outline.set_edgecolor(ap._frame_color)
        #colorbar.dividers.set_color(ap._frame_color)       # possibly broken in mpl 2.0.0
    else:
        print("--> you need to define both colorbar_location and colorbar_label to plot a colorbar")


    if 'execute_once' in kwargs:
        if (isinstance(kwargs['execute_once'], (list,tuple))):
            for code in kwargs['execute_once']:
                exec(code)
        else:
            print("Please give a list of commands to execute.")

    if 'out' in kwargs:
        fig.save(kwargs['out'], dpi=300, transparent=True)
        print("--> saved file as "+kwargs['out'])
    else:
        fig.save('grid_plot.png', dpi=300, transparent=True)
        print("--> saved file as grid_plot.png")

###################################################################################################
