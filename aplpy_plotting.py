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
colorbar_width        = 0.15	# relative to panel size
scalebar_frame        = False
scalebar_linestyle    = 'solid'
scalebar_linewidth    = 2
scalebar_color        = 'red'
scalebar_fontsize     = 10.0	# only used in channel map to prevent bar sliding over map
beam_frame            = False
beam_color            = 'black'
ticks_color           = 'black'

# user-adjustable settings come with a default value
global tick_label_xformat; tick_label_xformat = 'hh:mm:ss.ss'
global tick_label_yformat; tick_label_yformat = 'dd:mm:ss.ss'
global ticks_xspacing; ticks_xspacing = 0.2		# unit: degrees!
global ticks_yspacing; ticks_yspacing = 0.2 		# unit: degrees!
global ticks_minor_frequency; ticks_minor_frequency = 5

# define new viridis colormap with less dark blue
global viridis_cropped; viridis_cropped = colors.ListedColormap(mpl.cm.viridis(np.linspace(0.1,1.0,100)))

###################################################################################################


# main plotting function
#
# aplpy_plot('abc.fits', 
#			 out  = 'abc.png',       
#			 cmap = 'jet', 
#			 vmin = 0, 
#			 vmax = 1, 
#			 stretch  = 'linear', 
#			 recenter = SkyCoord('01h23m45.6s 12d34m45.6s'), 
#			 contour  = ['contour.fits', [1,2,3], ['white', 'grey', 'black']], 
#			 colorbar_location = 'right', 
#			 colorbar_label  = 'Jy/beam', 
#			 scalebar_length = 1.0, 
#			 scalebar_label  = 'string', 
#			 scalebar_corner = 'bottom', 
#			 beam_corner     = 'bottom left',
#			 overlay  = [['circle', SkyCoord('01h23m45.6s 12d34m45.6s'), 1.0*u.arcmin, {'linewidth': 1.0}]
#			 )
#
# the only mandatory argument is fitsfile
# grayscale plotting is accessible via cmap='grayscale'
# contours can be overplotted with contour=[...] where 
# [...] is a list of file name, list of contour levels 
# and corresponding colors

def aplpy_plot(fitsfile, **kwargs):
	
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
	if 'scalebar_length' and 'scalebar_label' and 'scalebar_corner' in kwargs:
		fig.add_scalebar(length=kwargs['scalebar_length'], label=kwargs['scalebar_label'], corner=kwargs['scalebar_corner'], frame=scalebar_frame)
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
	fig.ticks.set_xspacing(ticks_xspacing)
	fig.ticks.set_yspacing(ticks_yspacing)
	fig.ticks.set_minor_frequency(ticks_minor_frequency)
	fig.ticks.set_color(ticks_color)

	if 'out' in kwargs:
		fig.save(kwargs['out'], dpi=300, transparent=True)
		print "--> saved file as "+kwargs['out']
	else:
		fig.save(os.path.splitext(fitsfile)[0]+'.png', dpi=300, transparent=True)
		print "--> saved plot as "+os.path.splitext(fitsfile)[0]+'.png'



###################################################################################################

# channel map plotting function
#
# aplpy_channel_map(fitscube, 
#					ncols, 
#					nrows, 
#					chan_start, 
#					chan_iter, 
#					chan_width=1.0, 
#					chan_vel0=0.0, 
#					cmap='jet', 
#					vmin=-0.05, 
#					vmax=3.5, 
#					stretch='linear', 
#					recenter = SkyCoord('01h23m45.6s 12d34m45.6s'), 
#					contour=['contour.fits', [1,2,3], ['white', 'grey', 'black']], 
#					beam_corner='bottom left', 
#					colorbar_cmap=mpl.cm.jet, 
#					colorbar_label='flux density [Jy/beam]', 
#					scalebar_length=1.0, 
#					scalebar_label='string', 
#					scalebar_corner='bottom',
#					out='cube.png'
#					)
#
# mandatory (ordered!) arguments are: fitscube, ncols, nrows, chan_start, chan_iter to set up a basic 
# channel map 'fitscube' in ncols columns and nrows rows starting at channel chan_start at steps of 
# chan_iter
# grayscale plotting is accessible via cmap='grayscale'

def aplpy_channel_map(fitscube, ncols, nrows, chan_start, chan_iter, **kwargs):
	
	print "--> plotting channel map "+fitscube
	if not ('chan_width' and 'chan_vel0' in kwargs):
		print "--> channel width and velocity of first channel not given. will not calculate velocity"
	
	if 'figsize' in kwargs:
		main_figsize = kwargs['figsize']
	else:
		main_figsize = (8.27, 11.69)	# A4 in inches
	
	main_fig = plt.figure(figsize=main_figsize)
	for i in np.arange(nrows*ncols):
		
		# get subplot specific info
		subplt_size = [0.05+(i%ncols)*0.9/ncols, 0.95-np.ceil((i+1)/ncols)*0.9/nrows, 0.9/ncols, 0.9/nrows]
		#chn_slice   = [chan_start+i*chan_iter,0]	# if 4th axis is present
		chn_slice   = [chan_start+i*chan_iter]	# if no 4th axis available

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
			fig.ticks.set_xspacing(ticks_xspacing)
			fig.ticks.set_yspacing(ticks_yspacing)
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
			fig.ticks.set_xspacing(ticks_xspacing)
			fig.ticks.set_yspacing(ticks_yspacing)
			fig.ticks.set_minor_frequency(ticks_minor_frequency)
			fig.ticks.set_color(ticks_color)

#			fig.remove_scalebar()
			if 'scalebar_length' and 'scalebar_label' and 'scalebar_corner' in kwargs:
				fig.add_scalebar(length=kwargs['scalebar_length'], label=kwargs['scalebar_label'], corner=kwargs['scalebar_corner'], frame=scalebar_frame)
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

# position-velocity slice plotting
#
# aplpy_plot_pv('abc.fits', 
#				figsize = (8.27, 11.69),
#				vmin    = 0,
#				vmax    = 100,
#				stretch = 'linear',
#				cmap    = ap.viridis_cropped,
#				contour = [['xyz.fits', [1,2,3], 'black']], 
#				colorbar_location = 'right', 
#				colorbar_label = 'intensity [Jy\,beam$^{-1}$]', 
#				xlabel = 'offset [arcsec]', 
#				ylabel = 'velocity [km\,s$^{-1}$', 
#				out = 'abc.png'
#				)

def aplpy_plot_pv(fitspv, **kwargs):
	
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
#	fig.tick_labels.set_xformat(tick_label_xformat_pv)
#	fig.tick_labels.set_yformat(tick_label_yformat_pv)
	fig.ticks.show()
#	fig.ticks.set_xspacing(ticks_xspacing)
#	fig.ticks.set_yspacing(ticks_yspacing)
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

# grid plotting function
#
# props = {'boxstyle': "round", 'facecolor': "w", 'edgecolor': "black", 'linewidth': 0.5, 'alpha': 0.8}
# aplpy_map_grid([image1.fits,image2.fits], 
#				 ncols, 
#				 nrows, 
#				 figsize  = (8.27, 11.69), 
#				 cmap     = 'jet', 
#				 vmin     = -0.05, 
#				 vmax     = 3.5, 
#				 stretch  = 'linear', 
#				 recenter = SkyCoord('01h23m45.6s 12d34m45.6s'), 
#				 contour  = [['contour1.fits', [1,2,3], ['white', 'grey', 'black']],['contour2.fits',[9,8,7],'blue']], 
#				 label_text      = [image1_label, image2_label], 
#				 label_kwargs    = {'bbox': props}, 
#				 overlay         = [['circle', SkyCoord('01h23m45.6s 12d34m45.6s'), 1.0*u.arcmin, {'linewidth': 1.0}]],
#				 beam_corner     = 'bottom left', 
#				 colorbar_cmap   = mpl.cm.jet, 
#				 colorbar_label  = 'flux density [Jy/beam]', 
#				 scalebar_length = 1.0, 
#				 scalebar_label  = 'string', 
#				 scalebar_corner = 'bottom'
#				 )
#
# mandatory (ordered!) arguments are: fitsimages, ncols, nrows
# to set up a grid of single plane maps 'fitsimages' in ncols columns and nrows rows 
# grayscale plotting is accessible via cmap='grayscale'
# each image can have its own contour, give contour options in the same order as fitsimages

def aplpy_map_grid(fitsimages, ncols, nrows, **kwargs):
	
	print "--> plotting map grid of these maps: ", fitsimages
	
	if 'figsize' in kwargs:
		main_figsize = kwargs['figsize']
	else:
		main_figsize = (8.27, 11.69)	# A4 in inches
	
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
#			fig.image.axes.set_adjustable('box-forced')
			
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
			fig.ticks.set_xspacing(ticks_xspacing)
			fig.ticks.set_yspacing(ticks_yspacing)
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
						print "--> other shapes are not yet implemented"
			
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
			fig.ticks.set_xspacing(ticks_xspacing)
			fig.ticks.set_yspacing(ticks_yspacing)
			fig.ticks.set_minor_frequency(ticks_minor_frequency)
			fig.ticks.set_color(ticks_color)

	#		fig.remove_scalebar()
			if 'scalebar_length' and 'scalebar_label' and 'scalebar_corner' in kwargs:
				fig.add_scalebar(length=kwargs['scalebar_length'], label=kwargs['scalebar_label'], corner=kwargs['scalebar_corner'], frame=scalebar_frame)
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
