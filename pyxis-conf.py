import mqt,stefcal,imager,lsm,std,ms

_ms_config = [];
_ms_config_lastms = None;

def msconfig (pattern,*args,**kw):
  """Adds MS-specific configuration function or variable assignemnys.
  'pattern', e.g. "3C147-C-*" is matched (using shell wildcards) against the MS name.
  'args' have two forms: either callable functions (which are expected to do configurations
  and variable asignment, or they should come as 'name',value pairs (or keyword arguments)
  Whenever runcal() below is called, it will see which patterns the MS matches (in the order
  that they were added), and perform configuration according to this.
  
  Example: msconfig("3C147-C*",cconf,'imager.cellsize','2arcsec',LSMBASE='3C147-C.LSM');
  This will call cconf(), and set imager.cellsize='2arcsec' and v.LSMBASE='3C147-C.LSM' when
  the MS name matches the pattern"""
  if not isinstance(pattern,str):
    abort("invalid msconfig() pattern '$pattern'");
  cmdlist = [];
  while args:
    if callable(args[0]):
      cmdlist.append(args[0]);
      args = args[1:];
    elif len(args)<2 or not isinstance(args[0],str):
      abort("invalid msconfig() argument '%s'"%str(args[0]));
    else:
      cmdlist.append(args[0:2]);
      args = args[2:];
  cmdlist += list(kw.iteritems());
  _ms_config.append((pattern,cmdlist));
  
def _MSCONFIG_Template ():
  import fnmatch
  global _ms_config_lastms;
  if MS != _ms_config_lastms:
    _ms_config_lastms = MS;
    for pattern,cmdlist in _ms_config:
      if fnmatch.fnmatch(MS,pattern):
        info("$MS matches msconfig pattern $pattern:");
        for cmd in cmdlist:
          if callable(cmd):
            info("  calling %s()"%cmd.__name__);
            cmd();
          else:
            info("  assigning %s=%s"%cmd);
            assign(*cmd);
  return MS

def dconf ():
  global NPIX,CLEAN_THRESH,THRESH_ISL,THRESH_PIX
  """Sets config and imaging options for VLA-D"""
  imager.npix = NPIX = 2048
  imager.cellsize = "8arcsec"
  imager.wprojplanes = 0
  imager.CLEAN_ALGORITHM = "clark"
  v.LSMREF = "${MS:BASE}.refmodel.lsm.html"
  THRESH_PIX,THRESH_ISL = 50,15
  CLEAN_THRESH = ".5mJy",".12mJy",".06mJy"

def cconf ():
  global NPIX,CLEAN_THRESH,THRESH_ISL,THRESH_PIX
  """Sets config and imaging options for VLA-C"""
  imager.npix = NPIX = 4096
  imager.cellsize = "2arcsec"
  imager.wprojplanes = 128
  imager.CLEAN_ALGORITHM = "csclean"
  v.LSMREF = "${MS:BASE}.refmodel.lsm.html"
  THRESH_PIX,THRESH_ISL = 50,15
  CLEAN_THRESH = ".4mJy",".1mJy",".05mJy"
  
# Now for things specific to this script here.
dconf()


