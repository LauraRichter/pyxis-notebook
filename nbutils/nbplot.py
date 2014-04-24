
import subprocess

import matplotlib.pylab as plt
# maybe eventually don't need to import these here as they will be imported in the notebook?
import Pyxis
import mqt,stefcal,imager,lsm,std,ms
from Pyxis.Commands import II

import pyfits
import numpy as np

import os

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

def is_iterable(x):
    """Checks if object is iterable (but not a string or 0-dimensional array)."""
    return hasattr(x, '__iter__') and not (getattr(x, 'shape', None) == ())

def to_list(x):
    if is_iterable(x): 
        return x
    else:
        return [x]

def inline_plot(image,sigma,mu=0.0,df='99.7'):
    d_sig = display_frac[df]
    imshow(i,vmin=mu-sigma*d_sig,vmax=mu+sigma*d_sig)

def get_im_stats(data, sigma_clip=5.0, tolerance=0.01):
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
    """Plot a color greyscale image of an image, using 
    the intensity range specified.

    Inputs
    ------
    im      - image in 2 x 2 numpy array format
    ax      - axes on which to plot the image
    df      - intensity level of the plot
    """ 
    # plot a greyscale figure on the given axes,
    # using the given color scaling factor
    
    # get estimates of the figure mu and sigma
    mu, sigma = get_im_stats(im)
    # plot greyscale plot
    d_sig = DISPSIGMA[df]
    ax.imshow(im,vmin=mu-sigma*d_sig,vmax=mu+sigma*d_sig)

def stokes_audit(hdu0):
    """Determine stokes panels present in fits file.

    Inputs
    ------
    imlist  - fits file hdu[0]

    Outputs
    -------
    stokes  - string of stokes panels in the fits file
    """ 
    #print hdu0.header['NAXIS'],
    #print hdu0.header['NAXIS3'],
    #print hdu0.header['CRPIX3'], 
    #print hdu0.header['CRVAL3'], 
    #print hdu0.header['CDELT3']

    if hdu0.header['NAXIS3']==1 and hdu0.header['CDELT3']==1.0:
        return 'I'

    if hdu0.header['NAXIS3']==2 and hdu0.header['CDELT3']==3.0:
        return 'IV'

    if hdu0.header['NAXIS3']==2 and hdu0.header['CDELT3']==1.0:
        return 'IQ'

    if hdu0.header['NAXIS3']==4 and hdu0.header['CDELT3']==1.0:
        return 'IQUV'

def nbPlotFits(imlist,stokes=None,df='99.7'):
    """Plots list of images inline.

    Inputs
    ------
    imlist  - list of image types to plot, as used in pyxis:
              MODEL_IMAGE, RESTORED_IMAGE, RESIDUAL_IMAGE
    stokes  - stokes panel(s) to plot, e.g. 'IQ'
              Default: all panels present in the fits file.
    df      - intensity level to plot
    """ 
    # fits image dimensions:
    #    channels x stokes x ra x dec 

    # plot each image in its subplot location
    for im in to_list(imlist):
        # open the fits file
        f=pyfits.open(getattr(imager, im))

        # determine stokes panels in the fits image 
        fits_stokes = stokes_audit(f[0])
        # if stokes parameter is not set, plot all stokes present
        if not stokes:
            stokes = fits_stokes
        nstokes = len(stokes)

        # plot each specified stokes image
        if nstokes > 1:

            # set up figure size and number of subplots within the image
            fig, axes = plt.subplots(1,nstokes,figsize=(16,5))

            for s in stokes:
                ind = fits_stokes.index(s)
                plotim = f[0].data[0,ind,:,:] 
                axes[ind].set_title(im+'    STOKES '+stokes[ind])       
                plotscaled(plotim,axes[ind],df=df)

        else:
            ind = fits_stokes.index(stokes)
            plotim = f[0].data[0,ind,:,:] 
            fig, axes = plt.subplots(1,1,figsize=(6,5))
            axes.set_title(im+'    STOKES '+stokes)       
            plotscaled(plotim,axes,df=df)    

def tiggerPlotFits(imlist):
    """Plots list of images in tigger.

    Inputs
    ------
    imlist  - list of image types to plot, as used in pyxis:
              MODEL_IMAGE, RESTORED_IMAGE, RESIDUAL_IMAGE
    """ 
    plotlist = ["tigger"]
    for im in to_list(imlist): plotlist.append(getattr(imager, im))
    subprocess.call(plotlist)

def plotFits(imlist,stokes=None,df='99.7',nbmode=NB_MODE):
    """Sends list of images to plotter, in specified mode.

    Inputs
    ------
    imlist  - list of image types to plot, as used in pyxis:
              MODEL_IMAGE, RESTORED_IMAGE, RESIDUAL_IMAGE
    stokes  - stokes panel(s) to plot.
              Default: all panels present in the fits file.
    df      - intensity level to plot
    nbmode  - plotting mode:
              'browser' - inline ipython notebook plotting
              'tigger' - interactive plotting with tigger
    """ 

    if nbmode == 'browser':
        nbPlotFits(imlist,stokes=stokes,df=df) 
    elif nbmode == 'tigger':
        tiggerPlotFits(imlist)
    else:
        print "TODO - sort out batch mode!"
        
def imagelocation(imlist):
    """Determine and print the location of an image of a specified type,
    created most recently with a pyxis pipeline.
    """
    print 'Image locations\n---------------'
    for im in to_list(imlist):
        print im, ': ', getattr(imager, im)   
    print

def wrap_angle(angle, period=2.0 * np.pi):
    """Wrap angle into interval centred on zero.

    This wraps the *angle* into the interval -*period* / 2 ... *period* / 2.

    Stolen from ephem_extra.py
    """
    return (angle + 0.5 * period) % period - 0.5 * period

def getCPdims(ms_name=None, g=None):
    """Determine the dimensions of a solution cpickle file.

    Inputs
    ------
    ms_name - measurement set
    g       - loaded cpickle file (optional)

    Outputs
    -------
    NA      - number of antennas
    NP      - number of polarisations
    NT      - number of time points
    NF      - number of frequencies
    full_ant_list - list of all antennas present
    """ 

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
    """Lists solutions for specified MS.

    Inputs
    ------
    ms_name - measurement set
    ants    - list of antennas, e.g. [0,1,2]
              Default: plots all antennas
    g       - loaded cpickle file (optional)
    """ 

    base_name = II(ms_name+"$SUFFIX")

    sol_types = ['gain','diffgain']

    for st in sol_types:
        cpfile = II("$OUTDIR/plots-"+base_name+"/"+base_name+"-s$STEP."+st+"-spw${ms.DDID}.cp")
        if os.path.isfile(cpfile): listSolPerType(ms_name, st, g) 

def listSolPerType(ms_name=None, soltype=None, g=None): 
    """Lists solutions for specified MS and solution type.

    Inputs
    ------
    ms_name - measurement set
    ants    - list of antennas, e.g. [0,1,2]
              Default: plots all antennas
    soltype - solution to list: 'gain', 'diffgain' or ['gain,'diffgain']
              Default: ['gain,'diffgain']
    g       - loaded cpickle file 
    """   

    base_name = II(ms_name+"$SUFFIX")
    if not g: 
        cpfile = II("$OUTDIR/plots-"+base_name+"/"+base_name+"-s$STEP."+soltype+"-spw${ms.DDID}.cp")
        g = cPickle.load(file(cpfile))
    
    print soltype, 'solutions'
    print '-'*(len(soltype)+10)
    
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
    print

def listSol():
    """Lists solutions produced in previous pyxis pipeline step.
    """

    print 'Pipeline step', II('$STEP'), '\n'

    for sub_ms in eval(II("$MS_List")):
        ms_name = sub_ms.split('/')[-1].split('.')[0]
        listSolPerMS(ms_name)
        print

def plotGainsPerType(ms_name=None, ants=None, pols=None, soltype=None):
    """Plots solutions produced by pyxis, for specified MS and soltype.

    Inputs
    ------
    ms_name - measurement set
    ants    - list of antennas, e.g. ['0','1','2']
              Default: plots all antennas
    pols    - list of polarisations, e.g. ['0','1','2']
              Default: plots all polarisations
    soltype - solution to plot: 'gain', 'diffgain' or ['gain,'diffgain']
              Default: ['gain,'diffgain']
    """

    base_name = II(ms_name+"$SUFFIX")
    cpfile = II("$OUTDIR/plots-"+base_name+"/"+base_name+"-s$STEP."+soltype+"-spw${ms.DDID}.cp")
    g = cPickle.load(file(cpfile))

    NA, NP, NT, NF, full_ant_list = getCPdims(ms_name, g)

    # format antenna lisr
    if not ants: ants = full_ant_list
    ants = to_list(ants)

    # format polarisation list
    if not pols: pols = range(NP)
    pols = [int(p) for p in pols]

    #ants = full_ant_list

    for cal in g['gains']:
        solns = g['gains'][cal]['solutions']     

        # check that all of the requested antennas are present in the solutions file
        # if they aren't, leave them out of the list
        ants = [ a for a in ants if a in solns ]
   
        # get time indices
        xt = range(NT)

        # set up figure
        ncols = 2*len(pols)
        nrows = len(ants)

        fig, axes = plt.subplots(nrows,ncols,figsize=(5.0*ncols,1.5*nrows),squeeze=False,sharex=True)
        figTitle = '\n'+cal+' solutions for '+base_name
        fig.suptitle(figTitle, fontsize=16, horizontalalignment='center',verticalalignment='baseline')
        
        # for all antennas and pols, plot amplitude
        for a, ant in enumerate(ants): 
            
            for p, pol in enumerate(pols):
                # plot amplitude
                y_amp = np.abs(solns[ant][pol][:,:]) 
                y0 = np.min(y_amp,axis=1)
                yn = np.max(y_amp,axis=1)
                axes[a,2*p].plot(xt,np.abs(solns[ant][pol][:,NF/2]))
                axes[a,2*p].fill_between(xt, y0, yn, facecolor='gray', alpha=0.5, interpolate=True)
                # label
                axlabel = 'ANT '+ant+' POL '+str(pol)+' AMP'
                axes[a,2*p].annotate(axlabel,xy=(0.97,0.85), xycoords='axes fraction', fontsize=12,
                    horizontalalignment='right', verticalalignment='bottom', color='red',fontweight='normal') 
                axes[a,2*p].set_title(axlabel)
                axes[a,2*p].set_xlim([np.min(xt),np.max(xt)])

        # for all antennas and pols, plot phase
        for a, ant in enumerate(ants):
            
            for p, pol in enumerate(pols):
                # plot phase
                y_phase = wrap_angle(np.angle(solns[ant][pol][:]))
                y0 = np.min(y_phase,axis=1) 
                yn = np.max(y_phase,axis=1) 
                axes[a,2*p+1].plot(xt, y_phase[:,NF/2])
                axes[a,2*p+1].fill_between(xt, y0, yn, facecolor='gray', alpha=0.5, interpolate=True)
                
                #label
                axlabel = 'ANT '+ant+' POL '+str(pol)+' PHS'
                axes[a,2*p+1].annotate(axlabel,xy=(0.97,0.85), xycoords='axes fraction', fontsize=12,
                    horizontalalignment='right', verticalalignment='bottom', color='green',fontweight='normal') 
                axes[a,2*p+1].set_title(axlabel)
                axes[a,2*p+1].set_xlim([np.min(xt),np.max(xt)])
                # if the y limits are close to plus/minus pi, set to pi
                if np.abs(np.max(yn) - np.pi) < 0.1:
                    axes[a,2*p+1].set_ylim([-np.pi,np.pi])

def plotGainsPerMS(ms_name=None, ants=None, pols=None, soltype=None):
    """Plots solutions produced by pyxis, for specified MS.

    Inputs
    ------
    ms_name - measurement set
    ants    - list of antennas, e.g. ['0','1','2']
              Default: plots all antennas
    pols    - list of polarisations, e.g. ['0','1','2']
              Default: plots all polarisations
    soltype - solution to plot: 'gain', 'diffgain' or ['gain,'diffgain']
              Default: ['gain,'diffgain']
    """

    base_name = II(ms_name+"$SUFFIX")

    if not soltype: soltype = ['gain','diffgain']
    soltype=to_list(soltype)

    for st in soltype:
        cpfile = II("$OUTDIR/plots-"+base_name+"/"+base_name+"-s$STEP."+st+"-spw${ms.DDID}.cp")
        if os.path.isfile(cpfile): plotGainsPerType(ms_name, ants=ants, pols=pols, soltype=st) 

def plotGains(ants=None,pols=None,soltype=None):
    """Plots solutions produced by pyxis.

    Inputs
    ------
    ants    - comma separated string of antennas, e.g. '0,1,2'
              Default: plots all antennas
    pols    - comma separated string of polarisations, e.g. '0,1,2'
              Default: plots all polarisations
    soltype - solution to plot: 'gain', 'diffgain' or ['gain,'diffgain']
              Default: ['gain,'diffgain']
    """

    if not ants is None: ants = ants.split(',')
    if not pols is None: pols = pols.split(',')    

    for sub_ms in eval(II("$MS_List")):
        ms_name = sub_ms.split('/')[-1].split('.')[0]
        # check if soltype is a list. If not, make it a list.
        plotGainsPerMS(ms_name,ants=ants,pols=pols,soltype=soltype)
        print

#if __name__ == '__main__':
#    NB_MODE = sys.argv



















