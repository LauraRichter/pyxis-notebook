
import subprocess

import matplotlib.pylab as plt
# maybe eventually don't need to import these here as they will be imported in the notebook?
import Pyxis
import mqt,stefcal,imager,lsm,std,ms
from Pyxis.Commands import II

import pyfits
import numpy as np



# for plotting gains in cpickles
import cPickle
import Timba.Meq.meq

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

def plotFits(imlist,df='99.7',nbmode=NB_MODE):
    if nbmode == 'browser':
        nbPlotFits(imlist,df='99.7') 
    elif nbmode == 'tigger':
        tiggerPlotFits(imlist)
    else:
        print "TODO - sort out batch mode!"
        
def imagelocation(imlist):
    print 'Image locations\n---------------'
    for im in imlist:
        print im, ': ', getattr(imager, im)   
    print

def wrap_angle(angle, period=2.0 * np.pi):
    """Wrap angle into interval centred on zero.

    This wraps the *angle* into the interval -*period* / 2 ... *period* / 2.

    Stoken from ephem_extra.py

    """
    return (angle + 0.5 * period) % period - 0.5 * period

def getCPdims(ms_name=None, g=None):

    base_name = II(ms_name+"$SUFFIX")
    if not g: 
        cpfile = II("$OUTDIR/plots-"+base_name+"/"+base_name+"-s$STEP.gain-spw${ms.DDID}.cp")
        g = cPickle.load(file(cpfile))
    
    # determine data characteristics
    arb_index = next(g['gains'].iterkeys())   
    d = g['gains'][arb_index]['solutions']
    full_ant_list = d.keys()
    NA = len(d)
    arb_index = next(d.iterkeys()) 
    d = np.array(d[arb_index]) # at this point in the dictionary it is data, in list-of-lists form
    NP = d.shape[0]
    NT = d.shape[1]
    NF = d.shape[2]

    return NA, NP, NT, NF, full_ant_list

def listSolPerMS(ms_name=None, g=None):

    base_name = II(ms_name+"$SUFFIX")
    if not g: 
        cpfile = II("$OUTDIR/plots-"+base_name+"/"+base_name+"-s$STEP.gain-spw${ms.DDID}.cp")
        g = cPickle.load(file(cpfile))
    
    print 'Gain solutions\n--------------'
    
    #print 'Solution file: ', cpfile
    print 'Config and SPW: ', base_name
    print 'Solutions: ', g['gains'].keys()
    
    # determine data characteristics
    NA, NP, NT, NF, full_ant_list = getCPdims(ms_name, g)
    
    print 'Number of polarisations :', NP
    print 'Number of times :', NT
    print 'Number of frequencies :', NF           
    print 'Number of antennas: ', NA
    print 'Antenas: ', [a for a in full_ant_list]

def listSol():

    for sub_ms in eval(II("$MS_List")):
        ms_name = sub_ms.split('/')[-1].split('.')[0]
        listSolPerMS(ms_name)
        print

def plotGainsPerMS(ms_name=None,ants=None):

    base_name = II(ms_name+"$SUFFIX")
    cpfile = II("$OUTDIR/plots-"+base_name+"/"+base_name+"-s$STEP.gain-spw${ms.DDID}.cp")
    g = cPickle.load(file(cpfile))

    NA, NP, NT, NF, tmp = getCPdims(ms_name, g)

    for cal in g['gains']:
        solns = g['gains'][cal]['solutions']        
        # get time indices
        xt = range(NT)
        
        # for all antennas and pols, plot amplitude
        fig, axes = plt.subplots(len(ants),NP,figsize=(17,5))
        figTitle = cal+' solutions for '+base_name
        fig.suptitle(figTitle, fontsize=20, horizontalalignment='center')
        
        for i, ant in enumerate(ants): #solns.keys(): 
            axes[i,0].set_ylabel('Amplitude')
            
            for pol in range(NP):
                # plot amplitude
                y0 = np.abs(solns[ant][pol][:,0]) 
                yn = np.abs(solns[ant][pol][:,NF-1]) 
                axes[i,pol].plot(xt,np.abs(solns[ant][pol][:,NF/2]))
                axes[i,pol].fill_between(xt, y0, yn, facecolor='gray', alpha=0.5, interpolate=True)
                # label
                axlabel = 'ANT '+ant+' POL '+str(pol)
                axes[i,pol].annotate(axlabel,xy=(0.93,0.9), xycoords='axes fraction', fontsize=10,
                    horizontalalignment='right', verticalalignment='bottom', color='red') 

        # for all antennas and pols, plot phase
        fig, axes = plt.subplots(len(ants),NP,figsize=(17,5))
        for i, ant in enumerate(ants):
            axes[i,0].set_ylabel('Phase')
            
            for pol in range(NP):
                # plot phase
                y0 = np.angle(solns[ant][pol][:,0]) 
                yn = np.angle(solns[ant][pol][:,NF-1]) 
                axes[i,pol].plot(xt, np.angle(solns[ant][pol][:,NF/2]))
                axes[i,pol].fill_between(xt, y0, yn, facecolor='gray', alpha=0.5, interpolate=True)
                
                #label
                axlabel = 'ANT '+ant+' POL '+str(pol)
                axes[i,pol].annotate(axlabel,xy=(0.93,0.9), xycoords='axes fraction', fontsize=10,
                    horizontalalignment='right', verticalalignment='bottom', color='red') 

def plotGains(ants=None):

    for sub_ms in eval(II("$MS_List")):
        ms_name = sub_ms.split('/')[-1].split('.')[0]
        plotGainsPerMS(ms_name,ants)
        print

#if __name__ == '__main__':
#    NB_MODE = sys.argv



















