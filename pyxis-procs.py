
import mqt,stefcal,imager,lsm,std,ms

## 2. Procedures
# Procedures are invoked from the command line (i.e. "pyxis runcal" or "pyxis per_ms[runcal]").
# Think of them as recipes or something like that.
# I've tried to keep this one simple and linear, with everything determined by the variables set above.
# The net result of this is that any processing stage can be re-created interactively in ipython, by 
# simply typing
# : import Pyxis
# : LSM=LSM1
# : stefcal.stefcal(restore=True);
# You can also restart the calibration at a specific step by supplying a goto_step>0 here.
def runcal (goto_step=1):
  ## initial calibration
  if goto_step > 1:
    info("########## restarting calibration from step $goto_step");
  # Calibration step -- this is just a user-defined label used in filenames (etc. "blahblah_s1"), which serves to keep the output from each step
  # of a pipeline separate. If this is numeric, then functions such as stefcal.stefcal() will increment it automatically first thing. Otherwise you can set
  # it yourself to some more fancy label. Here we also provide a way to hop to  particular step via goto_step
  v.STEP = goto_step-1;
  
  if goto_step < 2:
    # set the superglobal LSM
    v.LSM = LSM0
    info("########## solving for G with initial LSM");
    # no w-proj for dirty map to save time
    stefcal.stefcal(stefcal_reset_ifr_gains=True,dirty=dict(wprojplanes=0),restore=True);
    info("########## running source finder and updating model");
    ## now run pybdsm on restored image, output LSM will be given by variable cal.PYBDSM_OUTPUT
    lsm.pybdsm_search(threshold=7);
    ### merge new sources into sky model, give it a new name ($LSM1)
    lsm.tigger_convert("$LSM -a ${lsm.PYBDSM_OUTPUT} $LSM1 --rename -f");
  
  if goto_step < 3:
    info("########## solving for G with updated LSM (initial+pybdsm)");
    v.LSM = LSM1
    stefcal.stefcal(dirty=dict(wprojplanes=0));
    
  if goto_step < 4:
    info("########## re-solving for G to apply IFR solutions");
    stefcal.stefcal(dirty=dict(wprojplanes=0),restore=True);
    
    info("########## adding clean components to LSM");
    CCMODEL = II("ccmodel-ddid${ms.DDID}.fits");  # note the per-style variable interpolation done by the II() function
    ff = pyfits.open(imager.MODEL_IMAGE);
    dd = ff[0].data;
    dd *= 1.0769     # scale up to compensate for selfcal flux suppression
    # dd[dd<0] = 0;  # remove negative components
    ff.writeto(CCMODEL,clobber=True);
    # add model image to LSM
    lsm.tigger_convert("$LSM $LSM2 --add-brick=ccmodel:$CCMODEL:2 -f");

  if goto_step < 5:          
    info("########## solving for G with updated LSM (inital+pybdsm+cc)");
    v.LSM = LSM2
    stefcal.stefcal(dirty=dict(wprojplanes=0));
    
  if goto_step < 6:
    info("########## running DD solutions");
    v.LSM = LSM2
    # now, set dE tags on sources
    lsm.transfer_tags(LSMREF,LSM,tags="dE",tolerance=45*ARCSEC);
  
    # make final image
    stefcal.stefcal(dirty=dict(wprojplanes=0),diffgains=True,restore=True,label="dE"); 

def c_cal (goto_step=1):
  """Calibration for C-config data"""
  ## initial calibration
  if goto_step > 1:
    info("########## restarting calibration from step $goto_step");
  # Calibration step -- this is just a user-defined label used in filenames (etc. "blahblah_s1"), which serves to keep the output from each step
  # of a pipeline separate. If this is numeric, then functions such as stefcal.stefcal() will increment it automatically first thing. Otherwise you can set
  # it yourself to some more fancy label. Here we also provide a way to hop to  particular step via goto_step
  v.STEP = goto_step-1;
  
  if goto_step < 2:
    # set the superglobal LSM
    v.LSM = LSM0
    info("########## solving for G with initial LSM");
    # no w-proj for dirty map to save time
    stefcal.stefcal(stefcal_reset_all=True,
        dirty=dict(wprojplanes=0,npix=NPIX),
        restore=dict(npix=NPIX,threshold=CLEAN_THRESH[0],wprojplanes=128));
    info("########## running source finder and updating model");
    ## now run pybdsm on restored image, output LSM will be given by variable cal.PYBDSM_OUTPUT
    ### NB: select on radius to exclude any artefacts picked up around 3C147 itself
    lsm.pybdsm_search(thresh_pix=THRESH_PIX,thresh_isl=THRESH_ISL,select="r.gt.30s");
    ### merge new sources into sky model, give it a new name ($LSM1)
    lsm.tigger_convert("$LSM -a ${lsm.PYBDSM_OUTPUT} $LSM1 --rename -f");
  
  if goto_step < 3:
    info("########## solving for G+dE with updated LSM (initial+pybdsm)");
    v.LSM = LSM1
    # now, set dE tags on sources
    lsm.transfer_tags(LSMREF,LSM,tags="dE",tolerance=45*ARCSEC);
    stefcal.stefcal(stefcal_reset_all=True,diffgains=True,dirty=dict(wprojplanes=0,npix=NPIX));
    
  if goto_step < 4:
    info("########## re-solving for G to apply IFR solutions");
    v.LSM = LSM1
    stefcal.stefcal(diffgains=True,diffgain_apply_only=True,
      dirty=dict(wprojplanes=0,npix=NPIX),
      restore=dict(npix=NPIX,threshold=CLEAN_THRESH[1]));
    
    info("########## adding clean components to LSM");
    ff = pyfits.open(imager.MODEL_IMAGE);
    dd = ff[0].data;
    dd *= 1.0769     # scale up to compensate for selfcal flux suppression
    # dd[dd<0] = 0;  # remove negative components
    ff.writeto(LSM_CCMODEL,clobber=True);
    # add model image to LSM
    lsm.tigger_convert("$LSM $LSM2 --add-brick=ccmodel:$LSM_CCMODEL:2 -f");

  if goto_step < 5:
    info("########## re-running DD solutions");
    v.LSM = LSM2
    # make final image
    stefcal.stefcal(dirty=dict(wprojplanes=0,npix=NPIX),diffgains=True,
      restore=dict(npix=NPIX,threshold=CLEAN_THRESH[2]),
      label="dE"); 

  # make per-channel cube
  makecube(NPIX);
  # make noise images     
  makenoise();
  
def makenoise ():  
  # make noise images     
  addnoise();
  imager.make_image(channelize=1,dirty_image="$OUTFILE.noisecube.fits",npix=256,wprojplanes=0,stokes="I",column="MODEL_DATA");
  imager.make_image(dirty_image="$OUTFILE.noise.fits",npix=256,wprojplanes=0,stokes="I",column="MODEL_DATA");
  noise = pyfits.open(II("$OUTFILE.noise.fits"))[0].data.std();
  info(">>> maximum noise value is %.2f uJy"%(noise*1e+6));
  
def saveconf ():
  if OUTDIR and OUTDIR != ".":
    x.sh("cp pyxis-*.py pyxis-*.conf tdlconf.profiles $OUTDIR");

def makecube (npix=512,stokes="I"):
  imager.make_image(channelize=1,dirty_image="$OUTFILE.cube.fits",npix=npix,wprojplanes=0,stokes=stokes);
  
def swapfields (f1,f2):
  """Swaps two fields in an MS"""
  info("swapping FIELDs $f1 and $f2 in $MS");
  field = ms.msw(subtable="FIELD");
  for name in field.colnames():
    info("swapping column $name");
    col = field.getcol(name);
    col[f1],col[f2] = col[f2],col[f1];
    field.putcol(name,col);
  field.close();
  tab = ms.msw();
  fcol = tab.getcol("FIELD_ID");
  r1 = (fcol==f1)
  r2 = (fcol==f2)
  fcol[r1] = f2
  fcol[r2] = f1
  tab.putcol("FIELD_ID",fcol);
  tab.close();
  
def compute_vis_noise (noise=0):
  tab = ms.ms().query("FIELD_ID==%d"%ms.FIELD);
  spwtab = ms.ms(subtable="SPECTRAL_WINDOW");
  freq0 = spwtab.getcol("CHAN_FREQ")[ms.SPWID,0];
  global WAVELENGTH
  WAVELENGTH = 300e+6/freq0
  bw = spwtab.getcol("CHAN_WIDTH")[ms.SPWID,0];
  dt = INTEGRATION or tab.getcol("EXPOSURE",0,1)[0];
  dtf = (tab.getcol("TIME",tab.nrows()-1,1)-tab.getcol("TIME",0,1))[0]
  # close tables properly, else the calls below will hang waiting for a lock...
  tab.close();
  spwtab.close();
  info(">>> $MS freq %.2f MHz (lambda=%.2fm), bandwidth %.2g kHz, %.2fs integrations, %.2fh synthesis"%(freq0*1e-6,WAVELENGTH,bw*1e-3,dt,dtf/3600));
  if not noise:
    noise = SEFD/math.sqrt(2*bw*dt);
    info(">>> SEFD of %.2f Jy gives per-visibility noise of %.2f mJy"%(SEFD,noise*1000));
  else:
    info(">>> using per-visibility noise of %.2f mJy"%(noise*1000));
  return noise;

def addnoise (noise=0,rowchunk=100000):
  """adds noise to MODEL_DATA, writes to CORRECTED_DATA""";
  # compute expected noise
  noise = compute_vis_noise(noise);
  # fill MS with noise
    # setup stefcal options and run 
  info("Running turbo-sim to add noise to data");
  # setup args
  args = [ """${ms.MS_TDL} ${ms.CHAN_TDL} ms_sel.ms_ifr_subset_str=${ms.IFRS} noise_stddev=%g"""%noise ];
  mqt.run("${mqt.CATTERY}/Siamese/turbo-sim.py","simulate",section="addnoise",args=args);

  
