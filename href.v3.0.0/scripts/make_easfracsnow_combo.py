# Author: B Blake / J Carley (EMC), T Alcott / C Alexander / I Jankov (GSD)
#
# Script to calculate and plot 6-h HREFv2 PQPF using the EAS Fractional Coverage method
# Apply variable radii to probability field
#
# Script History Log:
# 2017-03-02  B Blake   new script for HREFv2
# 2017-03-22  B Blake   6-hr time lagged NAM member added
# 2017-05-31  T Alcott  calculate_eas_probability function added
#                       computational resources reduced from 16 cores to 1 core
# 2017-06-01  B Blake   modified script containing loops over each ensemble member type
# 2018-01-05  M Pyle    3-hr only version, no bias correction
# 2018-01-18  M Pyle    Make a version for weasd
# 2018-02-08  M Pyle    Attempt to unify all EAS snow scripts into one
# 2018-03-21  M Pyle    pygrib free version

import os, sys, time
import numpy as np
import math as m
from datetime import datetime, timedelta
from scipy import ndimage, optimize, signal
# from scipy.stats import threshold
# from netCDF4 import Dataset
import fortranfile as F

starttime = time.time()

print 'Processing probabilistic snow'

try:
  os.environ["WGRIB2"]
except KeyError:
  print "NEED module loaded to define WGRIB2"
  exit(1)

WGRIB2=os.environ.get('WGRIB2','trash')
print 'found WGRIB2 as ', WGRIB2



try:
  os.environ["HOMEhref"]
except KeyError:
  print "NEED TO DEFINE HOMEhref"
  exit(1)
HOMEhref=os.environ.get('HOMEhref','trash')
print 'found HOMEhref as ', HOMEhref

try:
  os.environ["COMINhiresw"]
except KeyError:
  print "NEED TO DEFINE COMINhiresw"
  exit(1)
COMINhiresw=os.environ.get('COMINhiresw','trash')
print 'found COMINhiresw as ', COMINhiresw

try:
  os.environ["COMINfv3"]
except KeyError:
  print "NEED TO DEFINE COMINfv3"
  exit(1)
COMINfv3=os.environ.get('COMINfv3','trash')
print 'found COMINfv3 as ', COMINfv3

try:
  os.environ["COMOUT"]
except KeyError:
  print "NEED TO DEFINE COMOUT"
  exit(1)
COMOUT=os.environ.get('COMOUT','trash')
print 'found COMOUT as ', COMOUT

try:
  os.environ["PDY"]
except KeyError:
  print "NEED TO DEFINE PDY"
  exit(1)
PDY=os.environ.get('PDY','trash')
print 'found PDY as ', PDY

try:
  os.environ["cyc"]
except KeyError:
  print "NEED TO DEFINE cyc"
  exit(1)
cyc=os.environ.get('cyc','trash')
print 'found cyc as ', cyc

try:
  os.environ["dom"]
except KeyError:
  print "NEED TO DEFINE dom"
  exit(1)
dom=os.environ.get('dom','trash')
print 'found dom as: ', dom


try:
  os.environ["DATA"]
except KeyError:
  print "NEED TO DEFINE DATA"
  exit(1)
DATA=os.environ.get('DATA','trash')
print 'found DATA as ', DATA

# input directory and config file
sys.path.append(HOMEhref)
from eas_config import *
# grib-2 template 
template = HOMEhref + '/fix/weasd_rrfs' + dom + 'template.grib2'

record = 1		# SNOW from SREF pgrb212 file

# get latest run times and create output directory for plots if it doesn't already exist
if not os.path.exists(COMOUT):
  os.system("mkdir -p " + COMOUT)

# accumulation interval (hours)
fcst_hour = int(sys.argv[1])
qpf_interval = int(sys.argv[2])

DATArun=DATA+'/snow_'+str(fcst_hour)
if not os.path.exists(DATArun):
  os.system("mkdir -p "+DATArun)
os.system("cd "+DATArun)

start_hour = int(fcst_hour - qpf_interval)

fhr_range=str(start_hour)+'-'+str(fcst_hour)


print 'fcst_hour, qpf_interval, start_hour: ', fcst_hour, qpf_interval, start_hour

# maximum radius (km)
slim = max(rlist)
alpha = 0.5




#------------------------------------------------------------------------------------------
# read in calibration coefficients

if dom == 'conus':
  members = ['rrfs01','rrfs02','rrfs03','rrfs04','rrfs05','rrfs06','rrfs07','rrfs08','rrfs09']
elif dom == 'ak':
  members = ['arw','fv3nc','arw2','hrrrak']
else:
  members = ['arw','fv3nc','arw2']

pqpf_6h_calibrate = 'no'
pqpf_3h_calibrate = 'no'
# print 'Calibration coefficients for '+mem+' members not found. Using raw model QPF'
   
#--------------------------------------------------------------------------------
#### FUNCTIONS AND ROUTINES ####

def simplewgrib2(txtfile):
  tmps= []
  with open(txtfile) as F1:
    i=1
    nx,ny=[int(x) for x in next(F1).split()]
    ilim=nx*ny
    while i <= ilim:
      tmp=[float(x) for x in next(F1).split()]
      tmps.append(tmp)
      i=i+1
    array2d = np.asarray(tmps,dtype=np.float32)
    array2d.shape = (ny,nx)
    return array2d,nx,ny
  F1.close()


def process_nam_qpf(file3,file4,fhr):


    fhour=fhr
    if fhr%3 is 1:
      shour=fhour-1
      print 'process_nam_qpf remainder 1'
      os.system(WGRIB2+' '+file3+' -match "WEASD:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf.txt ')
      qpf1,nx,ny=simplewgrib2('qpf.txt')

    if fhr%3 is 2:
      shour1=fhour-2
      shour2=fhour-1
      print 'process_nam_qpf remainder 2 - f02 minus f01'
      os.system(WGRIB2+' '+file3+' -match "WEASD:surface:%i'%shour1+'-%i'%fhour+'" -end -text qpf.txt')
      qpfa,nx,ny=simplewgrib2('qpf.txt')
      os.system('rm -f qpf.txt')

      print 'shour1: ', shour1
      print 'shour2: ', shour2
      print 'file4: ', file4
      os.system(WGRIB2+' '+file4+' -match "WEASD:surface:%i'%shour1+'-%i'%shour2+'" -end -text qpf.txt')
      qpfb,nx,ny=simplewgrib2('qpf.txt')
      os.system('rm -f qpf.txt')

      qpf1=qpfa-qpfb

    if fhr%3 is 0:
      shour1=fhour-3
      shour2=fhour-2
      fhourm1=fhour-1
      print 'process_nam_qpf remainder 3 - f03 minus f02'
      os.system(WGRIB2+' '+file3+' -match "WEASD:surface:%i'%shour1+'-%i'%fhour+'" -end -text qpf.txt')
      qpfa,nx,ny=simplewgrib2('qpf.txt')
      os.system('rm -f qpf.txt')

      os.system(WGRIB2+' '+file4+' -match "WEASD:surface:%i'%shour1+'-%i'%fhourm1+'" -end -text qpf.txt')
      qpfb,nx,ny=simplewgrib2('qpf.txt')
      os.system('rm -f qpf.txt')

      qpf1=qpfa-qpfb

    return qpf1


# calculate footprint routine
def get_footprint(r):
    footprint = (np.ones(((r/dx)*2+1,(r/dx)*2+1))).astype(int)
    footprint[int(m.ceil(r/dx)),int(m.ceil(r/dx))] = 0
    dist = ndimage.distance_transform_edt(footprint,sampling=[dx,dx])
    footprint = np.where(np.greater(dist,r),0,1)
    return footprint

def get_footprint_flexi(r,i,j,nx,ny):

    footprint = (np.ones(((r/dx)*2+1,(r/dx)*2+1))).astype(int)
    footprint[int(m.ceil(r/dx)),int(m.ceil(r/dx))] = 0
    dist = ndimage.distance_transform_edt(footprint,sampling=[dx,dx])
    footprint = np.where(np.greater(dist,r),0,1)

    rdx2 = int(r/dx)
    nx1=nx-1
    ny1=ny-1

# if test is on domain grid indices.  Indices for footprint_flexi are footprint relative

    if i < rdx2:
      footprint[0:rdx2-i,:]=0
    if j < rdx2:
      footprint[:,0:rdx2-j]=0
    if i > nx1-rdx2:
      footprint[nx1-rdx2-i:,:]=0
    if j > ny1-rdx2:
      footprint[:,ny1-rdx2-j:]=0

    return footprint


# ENSEMBLE AGREEMENT SCALES PROBABILITY FUNCTION
#
#   ensemble_qpf: 3-D QPF array (members x lat x lon)
#   t:            QPF threshold in inches
#   rlist:        list of potential EAS radii (km)
#   alpha:        similarity criteria parameter
#   dx:           model grid spacing (km)
#   p_smooth:     width of Gaussian filter for smoothing radius field (grid points)
#
def calculate_eas_probability(ensemble_qpf,t,rlist,alpha,dx,p_smooth):
    exceed3d = np.where(np.greater_equal(ensemble_qpf/1.0,t),1,0)
    nm, isize, jsize = np.shape(exceed3d)
    print 'np.shape: ', nm, isize, jsize
    pr1, pr2 = [], []  # create lists of ens member pairs
    for m1 in range(nm-1):
      for m2 in range(m1+1,nm):
        pr1.append(m1)
        pr2.append(m2)
    optrad = np.zeros((isize,jsize)).astype(float) + float(max(rlist)) # set initial radius to max for better smoothing
    rlist = sorted(rlist,reverse=True)
    for i in range(len(rlist)):  # loop through rlist to find smallest radius over which member pairs meet similarity criteria
      dcrit = alpha
#      dcrit = alpha + ((1 - alpha) * (rlist[i] / float(max(rlist))))
      footprint = get_footprint(rlist[i])
      dijmean = np.zeros((isize,jsize)).astype(float)
      frac = np.zeros((nm,isize,jsize)).astype(float)
      for mem in range(nm):
        frac[mem,:,:] = np.around(signal.fftconvolve(exceed3d[mem,:,:],footprint,mode='same'))/float(np.sum(footprint))
      if i == 0:  # identify points where QPF threshold is not met within the maximum radius in any ens member
        fracsum = np.sum(frac,axis=0)
      for pair in range(len(pr1)):
        dijmean = dijmean + np.where(np.logical_and(np.greater(frac[pr1[pair]],0),np.greater(frac[pr2[pair]],0)),(frac[pr1[pair]]-frac[pr2[pair]])**2/(frac[pr1[pair]]**2+frac[pr2[pair]]**2),1)/float(len(pr1))
      optrad=np.where(np.less_equal(dijmean,dcrit),rlist[i],optrad)
#      p = np.where(np.less_equal(dijmean,dcrit),100.0*np.sum(frac,axis=0)/float(nm),p)
# smooth radius grid and zero out any dry areas
#    p = np.where(np.equal(fracsum,0),0,ndimage.filters.gaussian_filter(p,p_smooth))
    optrad = np.where(np.equal(fracsum,0),slim+5,ndimage.filters.gaussian_filter(optrad,p_smooth))
    return optrad

def calculate_pnt_probability(ensemble_qpf,t,p_smooth):
    exceed3d = np.where(np.greater_equal(ensemble_qpf/25.4,t),1,0)
    p_smooth_loc=p_smooth+2

    nm_use, isize, jsize = np.shape(exceed3d)
    pnt_prob = np.zeros((isize,jsize)).astype(float)

    for mem in range(nm_use):
        pnt_prob[:,:] = pnt_prob[:,:]+(exceed3d[mem,:,:]/float(nm_use))

    pnt_prob = 100.0 * pnt_prob
    pnt_prob = ndimage.filters.gaussian_filter(pnt_prob,p_smooth_loc)

    return pnt_prob



#--------------------------------------------------------------------------------
#### START OF SCRIPT ####

cy, cm, cd, ch = int(PDY[0:4]), int(PDY[4:6]), int(PDY[6:8]), int(cyc[0:2])
d0 = datetime(cy,cm,cd,ch,0)
starttime = d0+timedelta(start_hour/24.0)
endtime = d0+timedelta((start_hour+qpf_interval)/24.0)
memfiles = {}
memfiles4 = {}
latency4 = {}
itimes = []
fhours = []

latency = min_latency
stop = max_latency
stopnam = max_latency_nam

wgribdate=PDY+cyc

# get dimensions and message from template file
# grbs = pygrib.open(template)
# grbtmp = grbs[record]
# lats, lons = grbtmp.latlons()
# grbs.close()

os.system(WGRIB2+' '+template+' -rpn rcl_lat -text lat.txt  -rpn rcl_lon -text lon.txt')

lons,nx,ny=simplewgrib2('lon.txt')
lats,nx,ny=simplewgrib2('lat.txt')

nlats, nlons = np.shape(lats)

# define mask - NAM nest grid interpolated to grid 227 has undefined values
print 'dom here at decision point: ', dom
if dom == 'conus_avoid':
  print 'defining conus maskfile'
  maskfile = HOMEhref + '/fix/nam_mask.grib2'

print 'dom here at decision point2: ', dom
if dom == 'ak':
  print 'defining ak maskfile'
  maskfile = HOMEhref + '/fix/akhref_mask.grib2'

print 'dom here at decision point3: ', dom
if dom == 'conus_avoid' or dom == 'ak':
  print 'opening the maskregion stuff'
  print 'maskfile is: ', maskfile

  os.system(WGRIB2+' '+maskfile+'  -text mask.txt ')
  undefmask,nx,ny=simplewgrib2('mask.txt')
  undefmask=np.ma.masked_greater(undefmask,9.0e+20)
  maskregion = np.ma.filled(undefmask,-9999)


if qpf_interval == 1:
  print 'defined 1 h outbase'
  outbase = 'href.t'+cyc[0:2]+'z.'+dom+'.snow01_easfrac.f%02d'%(start_hour+qpf_interval)+'.grib2'
  incr =  1
if qpf_interval == 3:
  print 'defined 3 h outbase'
  outbase = 'href.t'+cyc[0:2]+'z.'+dom+'.snow03_easfrac.f%02d'%(start_hour+qpf_interval)+'.grib2'
  incr =  3
if qpf_interval == 6:
  print 'defined 6 h outbase'
  outbase = 'href.t'+cyc[0:2]+'z.'+dom+'.snow06_easfrac.f%02d'%(start_hour+qpf_interval)+'.grib2'
  incr =  3

outfile = DATA  + '/' + outbase

if os.path.exists(outfile):
  os.system('rm -f '+outfile)

prob = {}
qpf = {}
memcount = 0
check = 0
check_hrrr = 0

print 'members: ', members

for mem in members:
  while (len(itimes) < memcount+2) and (latency <= stop) and (start_hour+qpf_interval+latency <= 60):
    memname=mem[0:4]
    memnum=mem[4:6]
    print len(itimes), memcount+2
    itime = starttime-timedelta((start_hour+latency)/24.0)
    if memname == 'rrfs':
      file3 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day+'/fv3s.t%02d'%itime.hour+'z.m'+memnum+'.f%02d'%(start_hour+latency+incr)+'.grib2'
      file6 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day+'/fv3s.t%02d'%itime.hour+'z.m'+memnum+'.f%02d'%(start_hour+latency+incr+incr)+'.grib2'

    elif memname == 'nmmb':
      file3 = COMINhiresw + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/hiresw.t%02d'%itime.hour+'z.nmmb_5km.f%02d'%(start_hour+latency+incr)+'.'+dom+'.grib2'
      file6 = COMINhiresw + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/hiresw.t%02d'%itime.hour+'z.nmmb_5km.f%02d'%(start_hour+latency+incr+incr)+'.'+dom+'.grib2'

    elif memname == 'arw2':
      file3 = COMINhiresw + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/hiresw.t%02d'%itime.hour+'z.arw_5km.f%02d'%(start_hour+latency+incr)+'.'+dom+'mem2.grib2'
      file6 = COMINhiresw + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/hiresw.t%02d'%itime.hour+'z.arw_5km.f%02d'%(start_hour+latency+incr+incr)+'.'+dom+'mem2.grib2'

# have the FV3 preprocessing generate 3 h WEASD fields....handle like Hiresw?
    elif memname == 'fv3s':
      file3 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day+'/fv3s.t%02d'%itime.hour+'z.'+dom+'.f%02d'%(start_hour+latency+incr)+'.grib2'
      file6 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day+'/fv3s.t%02d'%itime.hour+'z.'+dom+'.f%02d'%(start_hour+latency+incr+incr)+'.grib2'


    if qpf_interval != 6:
      if os.path.exists(file3):
        print 'Found:',itime,'forecast hour',(start_hour+qpf_interval+latency)
        fhours.append(start_hour+qpf_interval+latency)
        itimes.append(itime)
        memfiles[itime] = file3
        print 'utilizing file3 in memfiles: ', file3
      else:
        print 'Missing:',itime,'forecast hour',(start_hour+qpf_interval+latency)
#        if os.path.exists(file3alt):
#          print 'alt Found:',itime_alt,'forecast hour',(start_hour+qpf_interval+latency)
#          fhours.append(start_hour+qpf_interval+latency+6)
#          print 'defined file3 fhours: ', start_hour+qpf_interval+latency+6
#          itimes.append(itime_alt)
#          memfiles[itime_alt] = file3alt
#          print 'using file3alt which is: ', file3alt
#        else:
#          print 'Completely missing:',itime,'forecast hour',(start_hour+qpf_interval+latency)
    else:
      print 'using 6 h block portion'
      if os.path.exists(file3) and os.path.exists(file6):
        print 'Found:',itime,'forecast hour',(start_hour+qpf_interval+latency)
        itimes.append(itime)
        fhours.append(start_hour+qpf_interval+latency)
        print 'defined file3 for 6h fhours: ', start_hour+qpf_interval+latency
        memfiles[itime] = [file3,file6]
      else:
        print 'Missing:',itime,'forecast hour',(start_hour+qpf_interval+latency)
#        if os.path.exists(file3alt) and os.path.exists(file6alt):
#          print 'alt Found:',itime_alt,'forecast hour',(start_hour+qpf_interval+latency)
#          itimes.append(itime_alt)
#          fhours.append(start_hour+qpf_interval+latency+6)
#          print 'defined fhours in alt 6h block: ', start_hour+qpf_interval+latency+6
#          memfiles[itime_alt] = [file3alt,file6alt]
#          print 'using file3alt which is: ', file3alt
#          print 'using file6alt which is: ', file6alt
#        else:
#          print 'Even alt is missing:',itime,'forecast hour',(start_hour+qpf_interval+latency)
#          print 'file3alt: ', file3alt
#          print 'file6alt: ', file6alt

    if mem == 'nam' or mem == 'hrrr' or mem == 'hrrrak' :
      latency = latency + 6
    else:
      print '12 hour latency member'
      latency = latency + 12

  if len(itimes) == (memcount+1):
    print 'Found 1 '+mem+' members valid:',starttime,'-',endtime
  else:
    print 'Could not find 2 '+mem+' members valid for start hour',start_hour
#    sys.exit(1)

#### READ IN QPF ####
  print 'here a memcount: ', memcount
  for itime in itimes[memcount:memcount+2]:	# Hi Res Window ARW members
    if qpf_interval == 6:
      file3,file6 = memfiles[itime]
    else:
      file3 = memfiles[itime]

    if dom == 'conus':
      print 'Processing member',(1+memcount),'of',nm
    if dom != 'conus':
      print 'Processing member',(1+memcount),'of',nm_nonconus
    print 'file3 near read is ', file3

    if qpf_interval != 1:
      if qpf_interval == 3:
        fhour=fhours[memcount]
        shour=fhour-3
      if qpf_interval == 6:
        fhour=fhours[memcount]-3
        shour=fhour-3
      print 'for file3 shour fhour: ', shour,fhour
      os.system(WGRIB2+' '+file3+' -match "WEASD:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf.txt')
      qpf3,nx,ny=simplewgrib2('qpf.txt')

      undefsnow,nx,ny=simplewgrib2('qpf.txt')
      undefsnow=np.ma.masked_greater(undefsnow,9.0e+20)
      snowmaskregion = np.ma.filled(undefsnow,-9999)
      qpf3 = np.where(np.equal(snowmaskregion,-9999),0,qpf3)

    else:

### HERE 

      if mem == 'nam':
        fcst_hour=fhours[memcount]
        file4 = memfiles4[itime]

        print 'call process_nam_qpf with file3: ', file3
        print 'call process_nam_qpf with file4: ', file4
        print 'call process_nam_qpf with fcst_hour: ', fcst_hour
        qpf1=process_nam_qpf(file3,file4,fcst_hour)

      else:
        print 'ready from file3 in 1h: ', file3
        fhour=fhours[memcount]
        shour=fhour-1
        print 'WEASD accum from : ', shour, ' to : ', fhour
        os.system(WGRIB2+' '+file3+' -match "WEASD:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf.txt')
        qpf1,nx,ny=simplewgrib2('qpf.txt')

        undefsnow,nx,ny=simplewgrib2('qpf.txt')
        undefsnow=np.ma.masked_greater(undefsnow,9.0e+20)
        snowmaskregion = np.ma.filled(undefsnow,-9999)
        qpf1 = np.where(np.equal(snowmaskregion,-9999),0,qpf1)

      print 'defining qpf from qpf1'
      qpf[itime] = qpf1*0.39370079

    if qpf_interval == 6:
      fhour=fhours[memcount]
      shour=fhour-3
      print '6h shour fhour: ', shour,fhour
      os.system(WGRIB2+' '+file6+' -match "WEASD:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf.txt')
      qpf6,nx,ny=simplewgrib2('qpf.txt')
      undefsnow,nx,ny=simplewgrib2('qpf.txt')
      undefsnow=np.ma.masked_greater(undefsnow,9.0e+20)
      snowmaskregion = np.ma.filled(undefsnow,-9999)
      qpf6 = np.where(np.equal(snowmaskregion,-9999),0,qpf6)
      qpf[itime] = (qpf3 + qpf6) * 0.39370079
    if qpf_interval == 3:
      qpf[itime] = qpf3*0.39370079   # Apply 10:1 SLR, convert to inches

    
    print 'max of qpf as snow: ', np.max(qpf[itime])
    print 'dom is: ', dom
    if dom == 'conus_avoid' or dom == 'ak':
      print 'qpf shape: ', np.shape(qpf)
      print 'maskregion shape: ', np.shape(maskregion)
      qpf[itime] = np.where(np.equal(maskregion,-9999),0,qpf[itime])

    # calculate threshold exceedance
    if qpf_interval == 1:
      thresh_use=snow_1h_thresh
    if qpf_interval == 3:
      thresh_use=snow_3h_thresh
    if qpf_interval == 6:
      thresh_use=snow_6h_thresh

    for t in thresh_use:
      for size in rlist:
        if memcount == 0:
          prob[t,size] = np.zeros((nlats,nlons))
        exceed = np.where(np.greater_equal(qpf[itime],t),1,0)
        filter_footprint = get_footprint(size)
        prob[t,size] = prob[t,size] + signal.fftconvolve(exceed,filter_footprint,mode='same')

    qpf[memcount]=qpf[itime]
    memcount = memcount + 1
  latency = min_latency


#-------------------------------------------------------------------#

# Get final probabilities

# redefine number of members (in case a TL member is missing)
nm = len(itimes)
nm_use=nm
print 'nm members is: ', nm
ensemble_qpf = np.zeros((nm,nlats,nlons)).astype(float)
for mem in range(0,len(itimes)):
  print 'mem in range is: ', mem
  print 'max(qpf(itimes), max(qpf(mem): ', np.max(qpf[itimes[mem]]), np.max(qpf[mem])
#  ensemble_qpf[mem,:,:] = qpf[itimes[mem]]
  ensemble_qpf[mem,:,:] = qpf[mem]

# Get final probabilities
probfinal = np.zeros((nlats,nlons))

filter_footprint_10 = get_footprint(10)
filter_footprint_25 = get_footprint(25)
filter_footprint_40 = get_footprint(40)
filter_footprint_55 = get_footprint(55)
filter_footprint_70 = get_footprint(70)
filter_footprint_85 = get_footprint(85)
filter_footprint_100 = get_footprint(100)

for t in thresh_use:
  t3 = time.time()
  optrad = calculate_eas_probability(ensemble_qpf,t,rlist,alpha,dx,p_smooth)
  t4 = time.time()
  print 'Time for optrad routine:', t4-t3
  pnt_prob = calculate_pnt_probability (ensemble_qpf, t, p_smooth)


  for row in range((slim/dx),nlats - (slim/dx)):
    for column in range((slim/dx),nlons - (slim/dx)):
      rad = (optrad[row,column]).astype(int)

      if (2.5 <= rad < 17.5):
        footprint_use  =  float(np.sum(filter_footprint_10))
      elif (17.5 <= rad < 32.5):
        footprint_use  =  float(np.sum(filter_footprint_25))
      elif (32.5 <= rad < 47.5):
        footprint_use  =  float(np.sum(filter_footprint_40))
      elif (47.5 <= rad < 62.5):
        footprint_use  =  float(np.sum(filter_footprint_55))
      elif (62.5 <= rad < 77.5):
        footprint_use  =  float(np.sum(filter_footprint_70))
      elif (77.5 <= rad < 92.5):
        footprint_use  =  float(np.sum(filter_footprint_85))
      elif (92.5 <= rad <= 100):
        footprint_use  =  float(np.sum(filter_footprint_100))
      elif (rad > 100):
        footprint_use  =  float(np.sum(filter_footprint_100))


      rdx=int(rad/dx)
# now insert something to compute the smaller footprint_use if in proper row/colum using flexi?

      if column < rdx:
         footprint_orig=np.sum(footprint_use)
         footprint_use = np.sum(get_footprint_flexi(rad,column,row,nlons,nlats))
         if float(footprint_use)/float(footprint_orig) < 0.5:
#            print 'W bound less than 0.5'
            if row > rdx and row < nlats-1-rdx:
#              print 'column, row, reduced, orig: ', column, row, footprint_use, footprint_orig
              footprint_use=int(0.51*footprint_orig)
              print 'revised column, row, reduced, orig: ', column, row, footprint_use, footprint_orig
         if float(footprint_use)/float(footprint_orig) < 0.25:
              print 'corner boost'
              footprint_use=int(0.26*footprint_orig)

      if row < rdx:
         footprint_orig=np.sum(footprint_use)
         footprint_use = np.sum(get_footprint_flexi(rad,column,row,nlons,nlats))
         if float(footprint_use)/float(footprint_orig) < 0.5:
#            print 'S bound less than 0.5'
            if column > rdx and column < nlons-1-rdx:
#              print 'column, row, reduced, orig: ', column, row, footprint_use, footprint_orig
              footprint_use=int(0.51*footprint_orig)
              print 'revised column, row, reduced, orig: ', column, row, footprint_use, footprint_orig
         if float(footprint_use)/float(footprint_orig) < 0.25:
              print 'corner boost'
              footprint_use=int(0.26*footprint_orig)

      if column > nlons-1-rdx:
         footprint_orig=np.sum(footprint_use)
         footprint_use = np.sum(get_footprint_flexi(rad,column,row,nlons,nlats))
         if float(footprint_use)/float(footprint_orig) < 0.5:
#            print 'E bound less than 0.5'
            if row > rdx and row < nlats-1-rdx:
#              print 'column, row, reduced, orig: ', column, row, footprint_use, footprint_orig
              footprint_use=int(0.51*footprint_orig)
              print 'revised column, row, reduced, orig: ', column, row, footprint_use, footprint_orig
         if float(footprint_use)/float(footprint_orig) < 0.25:
              print 'corner boost'
              footprint_use=int(0.26*footprint_orig)

      if row > nlats-1-rdx:
         footprint_orig=np.sum(footprint_use)
         footprint_use = np.sum(get_footprint_flexi(rad,column,row,nlons,nlats))
         if float(footprint_use)/float(footprint_orig) < 0.5:
#            print 'N bound less than 0.5'
            if column > rdx and column < nlons-1-rdx:
#              print 'column, row, reduced, orig: ', column, row, footprint_use, footprint_orig
              footprint_use=int(0.51*footprint_orig)
              print 'revised column, row, reduced, orig: ', column, row, footprint_use, footprint_orig
         if float(footprint_use)/float(footprint_orig) < 0.25:
              print 'corner boost'
              footprint_use=int(0.26*footprint_orig)





      if (2.5 <= rad < 17.5):
        probfinal[row,column] = prob[t,10][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
      elif (17.5 <= rad < 32.5):
        probfinal[row,column] = prob[t,25][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
      elif (32.5 <= rad < 47.5):
        probfinal[row,column] = prob[t,40][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
      elif (47.5 <= rad < 62.5):
        probfinal[row,column] = prob[t,55][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
      elif (62.5 <= rad < 77.5):
        probfinal[row,column] = prob[t,70][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
      elif (77.5 <= rad < 92.5):
        probfinal[row,column] = prob[t,85][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
      elif (92.5 <= rad <= 100):
        probfinal[row,column] = prob[t,100][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
      elif (rad > 100):
        probfinal[row,column] = prob[t,100][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
        optrad[row,column] = 0

# slight smoothing of probfinal
  probfinal = ndimage.filters.gaussian_filter(probfinal,1)

  if dom == 'conus_avoid' or dom == 'ak':
    print 'working final probability with mask'
    probfinal = np.where(np.equal(maskregion,-9999),0,probfinal)  # set to 0 for mask 
  t5 = time.time()
  print 'Time for get final probability routine for ',t, 'inch threshold: ',t5-t4
  print 'max of probfinal: ', np.max(probfinal)


  
  probstr=str(t*2.54)
  print 'probstr is: ', probstr

  byte=int(t*2.54*1000)
  byte44=0
  byte45=int(byte/65536)
  byte45rem=byte%65536
  byte46=int(byte45rem/256)
  byte47=byte45rem%256

# write binary file out of probfinal array, then import it into grib file using WGRIB2

  myfort = F.FortranFile('record_out.bin',mode='w')
  myfort.writeReals(probfinal)
  myfort.close()

  string="0:0:d="+wgribdate+":WEASD:surface:"+fhr_range+" hour acc fcst:prob > "+probstr
  print 'string used: ', string
  os.system(WGRIB2+' '+template+' -import_bin record_out.bin -set_metadata_str "'+string+'" -set_grib_type c3 -grib_out premod.grb')
  os.system(WGRIB2+' premod.grb -set_byte 4 12 197 -set_byte 4 17 0 -set_byte 4 24:35 0:0:0:0:0:255:0:0:0:0:0:0 -set_byte 4 36 '+str(nm_use)+' -set_byte 4 38:42 0:0:0:0:0 -set_byte 4 43 3 -set_byte 4 44 0 -set_byte 4 45 '+str(byte45)+' -set_byte 4 46 '+str(byte46)+' -set_byte 4 47 '+str(byte47)+' -set_byte 4 67 1 -append  -set_grib_type c3 -grib_out '+outfile)

  print 'byte, byte45, byte46, byte47: ', byte, byte45, byte46, byte47

  os.system('rm record_out.bin')
  os.system('rm premod.grb')

  print 'Wrote PSNOW for ', qpf_interval, ' to:',outfile, 'for ',t, 'inch threshold'
