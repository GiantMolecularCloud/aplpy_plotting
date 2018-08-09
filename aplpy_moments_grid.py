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
try:
    from matplotlib.cbook import MatplotlibDeprecationWarning as __MatplotlibDeprecationWarning__
    import warnings as __warnings__
    __warnings__.simplefilter('ignore', __MatplotlibDeprecationWarning__)
except:
    pass

###################################################################################################

def aplpy_moments_grid(fitsimages, scaling, **kwargs):

    """
    aplpy_moments_grid: plot moment maps side by side
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This is being tested at the moment. No description yet.

    """

    print("--> plotting map grid of these maps: ", fitsimages)


    # get keyword arguments
    main_figsize = kwargs.get('figsize', (8.27, 11.69))         # default: A4
    recenter     = kwargs.get('recenter', None)                 # default: do not recenter
    contours     = kwargs.get('contour', ['' for im in fitsimages])     # default: no contours
    contours     = kwargs.get('contours', ['' for im in fitsimages])    # default: no contours
    labels       = kwargs.get('labels', [im[:-5].replace('_','$\_$') for im in fitsimages])      # default: file name without ending
    label_kwargs = kwargs.get('label_kwargs', {})               # default: no kwargs
    scalebar     = kwargs.get('scalebar', None)                 # default: not displayed
    beam         = kwargs.get('beam', None)                     # default: not displayed
    out          = kwargs.get('out', 'moments.pdf')             # default: generic figure name


    # check inputs
    if not ( len(contours) == len(fitsimages) ):
        raise ValueError("Number of contours dows not match number of images.")


    # set up main figure containing all subplots and calculate number of rows
    main_fig = __plt__.figure(figsize = main_figsize)
    ncols = 3.0
    nrows = len(fitsimages)/3


    # loop over rows if more than 3 images are given
    for row in __np__.arange(len(fitsimages)/3):

        this_row_fits     = fitsimages[3*row:3*row+3]
        this_row_contours = contours[3*row:3*row+3]
        this_row_label    = labels[3*row:3*row+3]


        # loop over moments 0 to 2
        for mom in [0,1,2]:
            this_panel_fits     = this_row_fits[mom]
            this_panel_contours = this_row_contours[mom]
            this_panel_label    = this_row_label[mom]
            this_panel_scaling  = scaling[mom]
            this_panel_size     = (nrows,ncols,row*ncols+mom+1)


            # plot image
            fig = __aplpy__.FITSFigure(this_panel_fits,
                                       figure  = main_fig,
                                       subplot = this_panel_size
                                      )
            fig.show_colorscale(cmap    = this_panel_scaling[4],
                                vmin    = this_panel_scaling[0],
                                vmax    = this_panel_scaling[1],
                                stretch = this_panel_scaling[3]
                               )


            # recenter image
            if recenter is not None:
                if (len(recenter) == 2):
                    fig.recenter(recenter[0].ra.degree, recenter[0].dec.degree, radius=recenter[1].to(__u__.degree).value)
                elif (len(recenter) == 3):
                    fig.recenter(recenter[0].ra.degree, recenter[0].dec.degree, width=recenter[1].to(__u__.degree).value, height=recenter[2].to(__u__.degree).value)
                else:
                    raise SyntaxWarning("--> Specify center as SkyCoord(x,y) and either radius or width, height. Not recentering")


            # overplot contours
            if not this_panel_contours is '':
                if ( len(this_panel_contours) == 3 ):
                    fig.show_contour(data   = this_panel_contours[0],
                                     levels = this_panel_contours[1],
                                     colors = this_panel_contours[2]
                                    )
                else:
                    raise SyntaxWarning("--> Wrong number or format of contour parameters. Not plotting contours")


            # panel label
            if not this_panel_label is None:
                fig.add_label(0.5,0.9, this_panel_label.replace('_','$\_$'), color='black', relative=True, size=ap._velo_fontsize, **label_kwargs)


            # colorbar settings
            # show colorbars only in top images
            if ( row == 0 ):
                fig.add_colorbar()
                fig.colorbar.set_location('top')
                fig.colorbar.set_width(0.2)
                fig.colorbar.set_axis_label_text(this_panel_scaling[2])
                fig.colorbar.set_axis_label_font(size=ap._colorbar_fontsize)
                fig.colorbar.set_font(size=ap._colorbar_fontsize)
                fig.colorbar.set_frame_color(ap._frame_color)


            # set up panel ticks + labels
            fig.axis_labels.hide()
            fig.tick_labels.hide()
            fig.ticks.show()
            fig.ticks.set_xspacing(ap.ticks_xspacing.to(__u__.degree).value)
            fig.ticks.set_yspacing(ap.ticks_yspacing.to(__u__.degree).value)
            fig.ticks.set_minor_frequency(ap.ticks_minor_frequency)
            fig.ticks.set_color(ap._ticks_color)
            fig.frame.set_color(ap._frame_color)


            # add axis label and scale bar if bottom left plot
            if ( row == nrows-1 ) and ( mom == 0 ):
                fig.axis_labels.show()
                fig.tick_labels.show()
                fig.tick_labels.set_xformat(ap.tick_label_xformat)
                fig.tick_labels.set_yformat(ap.tick_label_yformat)
                fig.tick_labels.set_font(size=ap._tick_label_fontsize)
                fig.ticks.show()
                fig.ticks.set_xspacing(ap.ticks_xspacing.to(__u__.degree).value)
                fig.ticks.set_yspacing(ap.ticks_yspacing.to(__u__.degree).value)
                fig.ticks.set_minor_frequency(ap.ticks_minor_frequency)
                fig.ticks.set_color(ap._ticks_color)
                fig.axis_labels.set_font(size=ap._tick_label_fontsize)

                # add beam
                if not beam is None:
                    fig.add_beam()
                    fig.beam.show()
                    fig.beam.set_corner(beam)
                    fig.beam.set_frame(ap._beam_frame)
                    fig.beam.set_color(ap._beam_color)

                # add scalebar
                if not scalebar is None:
                    fig.add_scalebar(length = scalebar[0].to(__u__.degree).value,
                                     label  = scalebar[1],
                                     corner = scalebar[2],
                                     frame  = ap._scalebar_frame
                                    )
                    fig.scalebar.set_font(size = ap._scalebar_fontsize)
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

    __mpl__.pyplot.subplots_adjust(wspace=0.001, hspace=0.001)
#    __mpl__.pyplot.tight_layout()
    fig.save(out, dpi=300, transparent=True)
    print("--> saved file as "+out)

###################################################################################################
