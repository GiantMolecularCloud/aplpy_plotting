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
            
            fig = __aplpy__.FITSFigure(fitsimages[i], figure=main_fig, subplot=subplt_size)
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
                    fig.recenter(kwargs['recenter'][0].ra.degree, kwargs['recenter'][0].dec.degree, radius=kwargs['recenter'][1].to(__u__.degree).value)
                elif (len(kwargs['recenter']) == 3):
                    fig.recenter(kwargs['recenter'][0].ra.degree, kwargs['recenter'][0].dec.degree, width=kwargs['recenter'][1].to(__u__.degree).value, height=kwargs['recenter'][2].to(__u__.degree).value)
                else:
                    print("--> specify SkyCoord(x,y) and either radius or width, height. not recentering")
            
            # contours?
            if 'contour' in kwargs:
                for cont_i in __np__.arange(len(kwargs['contour'])):
                    if len(kwargs['contour'][cont_i]) == 3:
                        fig.show_contour(data=kwargs['contour'][cont_i][0], levels=kwargs['contour'][cont_i][1], colors=kwargs['contour'][cont_i][2])
                    else:
                        print("--> wrong number or format of contour parameters in image "+str(cont_i)+". not plotting contours")
            
            # ticks + labels
            fig.axis_labels.hide()
            fig.tick_labels.hide()
            fig.ticks.show()
            fig.ticks.set_xspacing(ap.ticks_xspacing.to(__u__.degree).value)
            fig.ticks.set_yspacing(ap.ticks_yspacing.to(__u__.degree).value)
            fig.ticks.set_minor_frequency(ap.ticks_minor_frequency)
            fig.ticks.set_color(ap.ticks_color)
            
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
            
            # beam settings
            if 'beam_corner' in kwargs:
                fig.add_beam()
                fig.beam.show()
                fig.beam.set_corner(kwargs['beam_corner'])
                fig.beam.set_frame(_beam_frame)
                fig.beam.set_color(_beam_color)
            
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

    #       fig.remove_scalebar()
            if 'scalebar_length' and 'scalebar_label' and 'scalebar_corner' in kwargs:
                fig.add_scalebar(length=kwargs['scalebar_length'].to(__u__.degree).value, label=kwargs['scalebar_label'], corner=kwargs['scalebar_corner'], frame=_scalebar_frame)
                fig.scalebar.set_font(size=ap._scalebar_fontsize)
                fig.scalebar.set_linestyle(ap._scalebar_linestyle)
                fig.scalebar.set_linewidth(ap._scalebar_linewidth)
                fig.scalebar.set_color(ap._scalebar_color)
        
        # colorbar settings
        if 'colorbar_cmap' and 'colorbar_label' in kwargs:
            ax1 = main_fig.add_axes([0.05+(ncols_f-1+0.05)*0.9/ncols_f, 0.05+0.5*0.9/nrows_f, 0.9*0.9/ncols_f, ap._colorbar_width*0.9/nrows_f])
            if 'vmin' and 'vmax' in kwargs:
                colorbar = __mpl__.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=__mpl__.colors.Normalize(vmin=kwargs['vmin'], vmax=kwargs['vmax']), orientation='horizontal')
            else:
                colorbar = __mpl__.colorbar.ColorbarBase(ax1, cmap=kwargs['colorbar_cmap'], norm=__mpl__.colors.Normalize(vmin=0.0, vmax=1.0), orientation='horizontal')
            colorbar.ax.tick_params(labelsize = ap._colorbar_fontsize)
            colorbar.ax.set_xticklabels([item.get_text() for item in colorbar.ax.get_xticklabels()], rotation=90)
            colorbar.set_label(kwargs['colorbar_label'], size = ap._colorbar_fontsize)
        else:
            print("--> you need to define both colorbar_location and colorbar_label to plot a colorbar")
    
    if 'out' in kwargs:
        fig.save(kwargs['out'], dpi=300, transparent=True)
        print("--> saved file as "+kwargs['out'])
    else:
        fig.save('grid_plot.png', dpi=300, transparent=True)
        print("--> saved file as grid_plot.png")

###################################################################################################
