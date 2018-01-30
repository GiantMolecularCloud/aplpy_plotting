# aplpy_plotting
repo for some wrapper functions to make plotting with APLpy more convenient

for APLpy see https://github.com/aplpy/aplpy


## Installation:
- clone this repository or download the scripts
- add the path where the file are located to your python path
    ```
    import sys
    sys.path.append('/your/path/to/aplpy_plotting')
    ```
- import the wrapper functions
    ```
    import aplpy_plotting as ap
    ```
- call the plotting wrappers with `ap.xxxxx`

## known issues:
- not all possible combinations of keyword arguments work (e.g. giving just vmin and vmax but not cmap).
- channel maps: first channel to plot and step need to be given as channel number, giving velocity or frequency is not supported yet.
- The functions sometimes offer more parameters than listed in the inline help (`?ap.aplpy_channel_map`). These are options that I introduced for specific plots but didn't have the time yet to update the docstring.

## some examples:

### A simple figure
```
aplpy.plot('map.fits')
```
![moment map](http://www2.mpia-hd.mpg.de/homes/krieger/images/SWAG_moment_map.png)

### A channel map
It's also very easy to get a channel map with the function aplpy_channel_map. No need to type ~100 lines of code anymore.
```
ap.aplpy_channel_map('datacube.fits',
    2,      # how many columns?
    5,      # how many rows?
    65,     # first channel to plot
    15,     # plot every n-th channel
    cmap            = 'viridis',
    contour         = ['pixelmask', 'mask.fits', [0.5], 'grey'],
    scalebar_length = 0.3452,
    scalebar_label  = '50\,pc',
    scalebar_corner = 'bottom'
    )
```
![channel map](http://www2.mpia-hd.mpg.de/homes/krieger/images/SWAG_channelmap.png)

### A position velocity diagram
It gets a bit more tricky when not both axis are spatial coordinates but aplpy_plot_pv does the job for you.
```
aplpy_plot_pv('pv_file.fits',
    xlabel = 'Galactic Longitude',
    ylabel = 'optical velocity [km\,s$^{-1}$]'
    )
```
![pV diagram](http://www2.mpia-hd.mpg.de/homes/krieger/images/SWAG_pV-l.png)
