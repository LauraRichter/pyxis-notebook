
import subprocess

import matplotlib.pylab as plt
# maybe eventually don't need to import these here as they will be imported in the notebook?
import Pyxis
import mqt,stefcal,imager,lsm,std,ms
import pyfits
import numpy as np

# dictionary for plotting scale limits
# display dictionary
DISPSIGMA = {}
DISPSIGMA['68'] = 1.
DISPSIGMA['95'] = 2.
DISPSIGMA['99.7'] = 3.
DISPSIGMA['99.99'] = 4.
DISPSIGMA['99.9999'] = 5.

# NB_MODE options:
# 'browser' - plots figures inline or with qt, depending how notebook is set up
# 'tigger' plots figures with tigger
NB_MODE='browser'

def inline_plot(image,sigma,mu=0.0,df='99.7'):
    d_sig = display_frac[df]
    imshow(i,vmin=mu-sigma*d_sig,vmax=mu+sigma*d_sig)

def get_im_stats(data,sigma_clip=5.0,tolerance=0.01):
    """Compute the variance by iteratively removing outliers greater than a given sigma
    until the mean changes by no more than tolerance.

    Inputs
    ------
    data - 1d numpy array of data to compute variance
    sigma_clip - the amount of sigma to clip the data before the next iteration
    tolerance - the fractional change in the mean to stop iterating

    Outputs
    -------
    variance - the final background variance in the sigma clipped image
    """

    #Initialise diff and data_clip and mean and std
    diff = 1
    mean = data.mean()
    data_clip = data
    while diff > tolerance:
        data_clip = data_clip[np.abs(data_clip)<mean+sigma_clip*data_clip.std()]
        newmean = data_clip.mean()
        diff = np.abs(mean-newmean)/(mean+newmean)
        mean = newmean
    return np.mean(data_clip), np.sqrt(np.var(data_clip))

def plotscaled(im,ax,df='99.7'):
    # plot a greyscale figure on the given axes,
    # using the given color scaling factor
    
    # get estimates of the figure mu and sigma
    mu, sigma = get_im_stats(im)
    # plot greyscale plot
    d_sig = DISPSIGMA[df]
    ax.imshow(im,vmin=mu-sigma*d_sig,vmax=mu+sigma*d_sig)

def nbPlotFits(imlist,df='99.7'):
    # set up figure size and number of subplots within the image
    fig, axes = plt.subplots(1,len(imlist),figsize=(16,5))
    # plot each image in its subplot location
    for i, im in enumerate(imlist):
        # open the fits file
        f=pyfits.open(getattr(imager, im))
        # TODO for later: add cabability to plot multiple stokes
        #print np.squeeze(f[0].data).shape
        plotim = np.squeeze(f[0].data) 

        axes[i].set_title(im)       
        plotscaled(plotim,axes[i],df)

def tiggerPlotFits(imlist):
    # plot all of the images in tigger
    plotlist = ["tigger"]
    for im in imlist: plotlist.append(getattr(imager, im))
    subprocess.call(plotlist)

def plotfits(imlist,df='99.7',nbmode=NB_MODE):
    if nbmode == 'browser':
        nbPlotFits(imlist,df='99.7') 
    elif nbmode == 'tigger':
        tiggerPlotFits(imlist)
    else:
        print "TODO - sort out batch mode!"   

#if __name__ == '__main__':
#    NB_MODE = sys.argv



















