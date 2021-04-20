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
# 2018-03-20  M Pyle    replace pygrib with wgrib2

import os, sys, time
import numpy as np
import math as m
from datetime import datetime, timedelta
from scipy import ndimage, optimize, signal
import fortranfile as F

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

starttime = time.time()

print 'Processing 24-h probabilistic QPF'

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
print 'found dom as ', dom

try:
  os.environ["DATA"]
except KeyError:
  print "NEED TO DEFINE DATA"
  exit(1)
DATA=os.environ.get('DATA','trash')
print 'found DATA as ', DATA

sys.path.append(HOMEhref)
from eas_config import *

# output directory
# grib-2 template 
template = HOMEhref + '/fix/pqpf_rrfs'+dom+'template.grib2'
record = 1		# PQPF from SREF pgrb212 file

# get latest run times and create output directory if it doesn't already exist

# accumulation interval (hours)
fcst_hour = int(sys.argv[1])
qpf_interval = int(sys.argv[2])
start_hour = int(fcst_hour - qpf_interval)

DATArun=DATA+'/qpf_'+str(fcst_hour)
if not os.path.exists(DATArun):
  os.system("mkdir -p " + DATArun)
os.system("cd "+DATArun)

fhr=fcst_hour
fhr_range=str(start_hour)+'-'+str(fcst_hour)


# maximum radius (km)
slim = max(rlist)
alpha = 0.5

os.system(WGRIB2+' '+template+' -rpn rcl_lat -text lat.txt  -rpn rcl_lon -text lon.txt')

lons,nx,ny=simplewgrib2('lon.txt')
lats,nx,ny=simplewgrib2('lat.txt')

# get dimensions and message from template file
# grbs = pygrib.open(template)
# grbtmp = grbs[record]
# lats, lons = grbtmp.latlons()
# grbs.close()

nlats, nlons = np.shape(lats)
# nlons, nlats = np.shape(lats)

print 'nlons, nlats: ', nlons, nlats

# define mask - NAM nest grid interpolated to grid 227 has undefined values
if dom == 'conus':
  maskfile = HOMEhref + '/fix/nam_mask.grib2'
if dom == 'ak':
  maskfile = HOMEhref + '/fix/akhref_mask.grib2'

if dom == 'avoidconus' or dom == 'ak':
  os.system(WGRIB2+' '+maskfile+' -text mask.txt ')
#  grbs2 = pygrib.open(maskfile)
  undefmask,nx,ny=simplewgrib2('mask.txt')
#  undefmask = grbs2[1].values

  undefmask=np.ma.masked_greater(undefmask,9.0e+20)
  maskregion = np.ma.filled(undefmask,-9999)
  print 'maskregion defined'
#  grbs2.close()

if not os.path.exists(COMOUT):
  os.system("mkdir -p " + COMOUT)

#------------------------------------------------------------------------------------------

if dom == 'conus':
  nm_use = nm_v3
  members = ['rrfs01','rrfs02','rrfs03','rrfs04','rrfs05','rrfs06','rrfs07','rrfs08','rrfs09']
elif dom == 'ak':
  nm_use = nm_ak
  members = ['arw','fv3nc','arw2','hrrrak']
else:
  nm_use = nm_nonconus
  members = ['arw','fv3nc','arw2']


for mem in members:
  if mem == 'arw':
    coeffs_arw = {}
  elif mem == 'nmmb':
    coeffs_nmmb = {}
  elif mem == 'fv3s':
    coeffs_fv3 = {}
  elif mem == 'arw2':
    coeffs_arw2 = {}
  elif mem == 'nam':
    coeffs_nam = {}
  elif mem == 'hrrr':
    coeffs_hrrr = {}

  pqpf_6h_calibrate = 'no'
  print 'Calibration coefficients for '+mem+' members not found. Using raw model QPF'
   

#--------------------------------------------------------------------------------
#### FUNCTIONS AND ROUTINES ####



def process_nam_qpf(file3,file4,fhr):


    fhour=fhr
    if fhr%3 is 1:
      shour=fhour-1
      print 'process_nam_qpf remainder 1'
      os.system(WGRIB2+' '+file3+' -match "APCP:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf.txt ')
      qpf1,nx,ny=simplewgrib2('qpf.txt')
      os.system('rm -f qpf.txt')


    if fhr%3 is 2:
      shour1=fhour-2
      shour2=fhour-1
      print 'process_nam_qpf remainder 2 - f02 minus f01'
      os.system(WGRIB2+' '+file3+' -match "APCP:surface:%i'%shour1+'-%i'%fhour+'" -end -text qpf.txt')
      qpfa,nx,ny=simplewgrib2('qpf.txt')
      os.system('rm -f qpf.txt')

      os.system(WGRIB2+' '+file4+' -match "APCP:surface:%i'%shour1+'-%i'%shour2+'" -end -text qpf.txt')
      qpfb,nx,ny=simplewgrib2('qpf.txt')
      os.system('rm -f qpf.txt')

      qpf1=qpfa-qpfb

    if fhr%3 is 0:
      shour1=fhour-3
      shour2=fhour-2
      fhourm1=fhour-1
      print 'process_nam_qpf remainder 3 - f03 minus f02'
      os.system(WGRIB2+' '+file3+' -match "APCP:surface:%i'%shour1+'-%i'%fhour+'" -end -text qpf.txt')
      qpfa,nx,ny=simplewgrib2('qpf.txt')
      os.system('rm -f qpf.txt')

      os.system(WGRIB2+' '+file4+' -match "APCP:surface:%i'%shour1+'-%i'%fhourm1+'" -end -text qpf.txt')
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
    exceed3d = np.where(np.greater_equal(ensemble_qpf/25.4,t),1,0)
    nm_use, isize, jsize = np.shape(exceed3d)
    print 'nm_use within get_footprint: ', nm_use
    print 'isize, jsize within get_footprint: ', isize, jsize
    pr1, pr2 = [], []  # create lists of ens member pairs
    for m1 in range(nm_use-1):
      for m2 in range(m1+1,nm_use):
        pr1.append(m1)
        pr2.append(m2)
    optrad = np.zeros((isize,jsize)).astype(float) + float(max(rlist)) # set initial radius to max for better smoothing
    rlist = sorted(rlist,reverse=True)
    for i in range(len(rlist)):  # loop through rlist to find smallest radius over which member pairs meet similarity criteria
      dcrit = alpha
#      dcrit = alpha + ((1 - alpha) * (rlist[i] / float(max(rlist))))
      footprint = get_footprint(rlist[i])
      dijmean = np.zeros((isize,jsize)).astype(float)
      frac = np.zeros((nm_use,isize,jsize)).astype(float)
      for mem in range(nm_use):
        frac[mem,:,:] = np.around(signal.fftconvolve(exceed3d[mem,:,:],footprint,mode='same'))/float(np.sum(footprint))
      if i == 0:  # identify points where QPF threshold is not met within the maximum radius in any ens member
        fracsum = np.sum(frac,axis=0)
      for pair in range(len(pr1)):
        dijmean = dijmean + np.where(np.logical_and(np.greater(frac[pr1[pair]],0),np.greater(frac[pr2[pair]],0)),(frac[pr1[pair]]-frac[pr2[pair]])**2/(frac[pr1[pair]]**2+frac[pr2[pair]]**2),1)/float(len(pr1))
      optrad=np.where(np.less_equal(dijmean,dcrit),rlist[i],optrad)
#      p = np.where(np.less_equal(dijmean,dcrit),100.0*np.sum(frac,axis=0)/float(nm_use),p)
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
itimes = []
fhours = []
latency = min_latency
stop = max_latency
stopnam = max_latency_nam

# create grib messages from template (only need to do this once)

wgribdate=PDY+cyc

# grbtmp['dataDate']=int('%i'%d0.year+'%02d'%d0.month+'%02d'%d0.day)
# grbtmp['dataTime']=int('%02d'%d0.hour+'00')
# grbtmp['startStep']=int(start_hour)
# grbtmp['endStep']=int(start_hour+qpf_interval)
# grbtmp['yearOfEndOfOverallTimeInterval']=endtime.year
# grbtmp['monthOfEndOfOverallTimeInterval']=endtime.month
# grbtmp['dayOfEndOfOverallTimeInterval']=endtime.day
# grbtmp['hourOfEndOfOverallTimeInterval']=endtime.hour
# grbtmp['scaleFactorOfUpperLimit']=3

if qpf_interval == 1:
  outbase = 'href.t'+cyc[0:2]+'z.'+dom+'.pqpf01_easfrac.f%02d'%(start_hour+qpf_interval)+'.grib2'
  incr = 1
  thresh_use=pqpf_1h_thresh
if qpf_interval == 3:
  outbase = 'href.t'+cyc[0:2]+'z.'+dom+'.pqpf03_easfrac.f%02d'%(start_hour+qpf_interval)+'.grib2'
  incr = 3
  thresh_use=pqpf_3h_thresh
if qpf_interval == 6:
  outbase = 'href.t'+cyc[0:2]+'z.'+dom+'.pqpf06_easfrac.f%02d'%(start_hour+qpf_interval)+'.grib2'
  incr = 3
  thresh_use=pqpf_6h_thresh
if qpf_interval == 12:
  outbase = 'href.t'+cyc[0:2]+'z.'+dom+'.pqpf12_easfrac.f%02d'%(start_hour+qpf_interval)+'.grib2'
  incr = 3
  thresh_use=pqpf_12h_thresh
if qpf_interval == 24:
  outbase = 'href.t'+cyc[0:2]+'z.'+dom+'.pqpf24_easfrac.f%02d'%(start_hour+qpf_interval)+'.grib2'
  incr = 3
  thresh_use=pqpf_24h_thresh

outfile = DATA + '/' + outbase

if os.path.exists(outfile):
  os.system('rm -f '+outfile)


# qpf is a dictionary - need unique key.  itime values are repeated, so is a bad choice.

prob = {}
qpf = {}
memcount = 0

## change 48 to 60...need to set different limits for different sources?

for mem in members:

  memname=mem[0:4]
  memnum=mem[4:6]
  print 'memname, memnum: ', memname, memnum

  while (len(itimes) < memcount+2) and (latency <= stop) and (start_hour+qpf_interval+latency <= 60):
    print 'len(itimes), memcount+2: ', len(itimes), memcount+2
    itime = starttime-timedelta((start_hour+latency)/24.0)
    print 'itime for this member: ', itime
    if memname == 'rrfs':
      file0 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/fv3s.t%02d'%itime.hour+'z.m'+memnum+'.f%02d'%(start_hour+latency)+'.grib2'
      print 'file0 is: ', file0
      file1 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/fv3s.t%02d'%itime.hour+'z.m'+memnum+'.f%02d'%(start_hour+latency+incr)+'.grib2'
      file2 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/fv3s.t%02d'%itime.hour+'z.m'+memnum+'.f%02d'%(start_hour+latency+2*incr)+'.grib2'
      file3 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/fv3s.t%02d'%itime.hour+'z.m'+memnum+'.f%02d'%(start_hour+latency+3*incr)+'.grib2'
      file4 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/fv3s.t%02d'%itime.hour+'z.m'+memnum+'.f%02d'%(start_hour+latency+4*incr)+'.grib2'
      file5 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/fv3s.t%02d'%itime.hour+'z.m'+memnum+'.f%02d'%(start_hour+latency+5*incr)+'.grib2'
      file6 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/fv3s.t%02d'%itime.hour+'z.m'+memnum+'.f%02d'%(start_hour+latency+6*incr)+'.grib2'
      file7 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/fv3s.t%02d'%itime.hour+'z.m'+memnum+'.f%02d'%(start_hour+latency+7*incr)+'.grib2'
      file8 = COMINfv3 + '.%02d'%itime.year+'%02d'%itime.month+'%02d'%itime.day + '/fv3s.t%02d'%itime.hour+'z.m'+memnum+'.f%02d'%(start_hour+latency+8*incr)+'.grib2'

    elif memname == 'fv3n':
      print 'no fv3nc here'
    elif memname == 'fv3s':
      print 'no fv3s here'
    elif memname == 'arw2':
      print 'no arw2 here'


    if qpf_interval == 1:
      if os.path.exists(file1):
        print 'Found:',itime,'forecast hour',(start_hour+1*incr+latency)
        fhours.append(start_hour+1*incr+latency)
        itimes.append(itime)
# define fully just in case
        memfiles[itime] = [file0,file1,file2,file3,file4,file5,file6,file7,file8]
      else:
        print 'did not find file1: ', file1
#        print 'trying to work file1alt: ', file1alt
#        if (os.path.exists(file1alt)):
#          fhours.append(start_hour+1*incr+latency+alt_fhrinc)
#          itimes.append(itime_alt)
#          memfiles[itime_alt] = [file0alt,file1alt,file2alt,file3alt,file4alt,file5alt,file6alt,file7alt,file8alt]
#        else:
#          print 'Alt cycle is missing as well'

    if qpf_interval == 3:
      if os.path.exists(file1):
        print 'Found:',itime,'forecast hour',(start_hour+1*incr+latency)
        fhours.append(start_hour+1*incr+latency)
        itimes.append(itime)
# define fully just in case
        memfiles[itime] = [file1,file2,file3,file4,file5,file6,file7,file8]
      else:
        print 'did not find file1: ', file1
#        print 'trying to work file1alt: ', file1alt
#        if (os.path.exists(file1alt)):
#          fhours.append(start_hour+1*incr+latency+alt_fhrinc)
#          itimes.append(itime_alt)
#          memfiles[itime_alt] = [file1alt,file2alt,file3alt,file4alt,file5alt,file6alt,file7alt,file8alt]
#        else:
#          print 'Alt cycle is missing as well'

    if qpf_interval == 6:
      if os.path.exists(file2):
        print 'Found:',itime,'forecast hour',(start_hour+2*incr+latency)
        itimes.append(itime)
        fhours.append(start_hour+1*incr+latency)
# define fully just in case
        memfiles[itime] = [file1,file2,file3,file4,file5,file6,file7,file8]
      else:
        print 'did not find file2 for qpf_6: ', file2
#        if (os.path.exists(file2alt)):
#          fhours.append(start_hour+1*incr+latency+alt_fhrinc)
#          itimes.append(itime_alt)
#          memfiles[itime_alt] = [file1alt,file2alt,file3alt,file4alt,file5alt,file6alt,file7alt,file8alt]
#        else:
#          print 'Alt cycle is missing as well'
#
    if qpf_interval == 12:
      if os.path.exists(file4):
        print 'Found:',itime,'forecast hour',(start_hour+4*incr+latency)
        itimes.append(itime)
        fhours.append(start_hour+1*incr+latency)
# define fully just in case
        memfiles[itime] = [file1,file2,file3,file4,file5,file6,file7,file8]
      else:
        print 'did not find file4 for qpf_12: ', file4
#        if (os.path.exists(file4alt)):
#          itimes.append(itime_alt)
#          fhours.append(start_hour+1*incr+latency+alt_fhrinc)
#          memfiles[itime_alt] = [file1alt,file2alt,file3alt,file4alt,file5alt,file6alt,file7alt,file8alt]
#        else:
#          print 'Alt cycle is missing as well'



    if qpf_interval == 24:
      if os.path.exists(file8):
        print 'Found:',itime,'forecast hour',(start_hour+8*incr+latency)
        itimes.append(itime)
        fhours.append(start_hour+1*incr+latency)
        memfiles[itime] = [file1,file2,file3,file4,file5,file6,file7,file8]
      else:
        print 'did not find file8 for qpf_24: ', file8
#        if (os.path.exists(file8alt)):
#          itimes.append(itime_alt)
#          fhours.append(start_hour+1*incr+latency+alt_fhrinc)
#          memfiles[itime_alt] = [file1alt,file2alt,file3alt,file4alt,file5alt,file6alt,file7alt,file8alt]
#        else:
#          print 'Alt cycle is missing as well'



    if mem == 'nam' or mem == 'hrrr' or mem == 'hrrrak':
      latency = latency + 6
    else:
      latency = latency + 12

  if len(itimes) == (memcount+1):
    print 'Found 1 '+mem+' members valid:',starttime,'-',endtime
  else:
    print 'Could not find 2 '+mem+' members valid for start hour',start_hour
#    sys.exit(1)

#### READ IN QPF ####
  for itime in itimes[memcount:memcount+2]:
    if qpf_interval == 24 or qpf_interval == 6 or qpf_interval == 12:

      file1,file2,file3,file4,file5,file6,file7,file8 = memfiles[itime]

# Process first 6 hours
      print 'Processing member',(1+memcount),'of',nm_use

      fhour=fhours[memcount]
      shour=fhour-3
      os.system(WGRIB2+' '+file1+' -match "APCP:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf1.txt ')
      qpf1,nx,ny=simplewgrib2('qpf1.txt')

## how get the proper hours for each file??
      fhour=fhours[memcount]+incr
      shour=fhour-3
      os.system(WGRIB2+' '+file2+' -match "APCP:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf2.txt ')
      qpf2,nx,ny=simplewgrib2('qpf2.txt')

      qpf12 = qpf1 + qpf2    

      print 'max of qpf12: ', np.max(qpf12)

      if qpf_interval == 6:
         print 'defined qpf[itime] from qpf12'
         qpf[itime]=qpf12
         print 'itime assigned: ', itime
         print 'used memcount instead: ', memcount

    if qpf_interval == 24 or qpf_interval == 12 :

## figure out fhour for two pieces here

      fhour=fhr+6
      fhour=fhours[memcount]+incr*2
      shour=fhour-3
      os.system(WGRIB2+' '+file3+' -match "APCP:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf3.txt ')
      qpf3,nx,ny=simplewgrib2('qpf3.txt')

      fhour=fhr+9
      fhour=fhours[memcount]+incr*3
      shour=fhour-3
      os.system(WGRIB2+' '+file4+' -match "APCP:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf4.txt ')
      qpf4,nx,ny=simplewgrib2('qpf4.txt')

# Process second 6 hours

      qpf34 = qpf3 + qpf4    

#### 12 h sum of the first two 6 h periods
      if qpf_interval == 12 :
        qpf[itime]= qpf12 + qpf34

    if qpf_interval == 24 :

# Process third 6 hours
## figure out fhour for two pieces here
      fhour=fhr+12
      fhour=fhours[memcount]+incr*4
      shour=fhour-3
      os.system(WGRIB2+' '+file5+' -match "APCP:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf5.txt ')
      qpf5,nx,ny=simplewgrib2('qpf5.txt')

      fhour=fhr+15
      fhour=fhours[memcount]+incr*5
      shour=fhour-3
      os.system(WGRIB2+' '+file6+' -match "APCP:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf6.txt ')
      qpf6,nx,ny=simplewgrib2('qpf6.txt')

      qpf56 = qpf5 + qpf6


# Process last 6 hours
## figure out fhour for two pieces here

      fhour=fhr+18
      fhour=fhours[memcount]+incr*6
      shour=fhour-3
      os.system(WGRIB2+' '+file7+' -match "APCP:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf7.txt ')
      qpf7,nx,ny=simplewgrib2('qpf7.txt')

      fhour=fhr+21
      fhour=fhours[memcount]+incr*7
      shour=fhour-3
      os.system(WGRIB2+' '+file8+' -match "APCP:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf8.txt ')
      qpf8,nx,ny=simplewgrib2('qpf8.txt')

      qpf78 = qpf7 + qpf8

      qpf[itime] = qpf12 + qpf34 + qpf56 + qpf78

######### 3 h APCP

    if qpf_interval == 3:

      file1,file2,file3,file4,file5,file6,file7,file8 = memfiles[itime]
      print 'from memfiles file1 for 3 h qpf: ', file1

# Process first 3 hours
      print 'Processing member',(1+memcount),'of',nm_use
      print 'fhours of mem: ', fhours[memcount]
      fhour=fhours[memcount]
      shour=fhour-3

      print 'shour: ', shour
      print 'fhour: ', fhour

      print 'nx, ny: ', nx, ny
      print 'for file1: ', file1
      os.system(WGRIB2+' '+file1+' -match "APCP:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf.txt ')
      qpfhere,nx,ny=simplewgrib2('qpf.txt')
      os.system('rm -f qpf.txt')
      qpf[itime]=qpfhere
      print 'from file1: ', file1
      print 'max 3 h qpf: ', np.max(qpf[itime])
      print 'mean 3 h qpf: ', np.mean(qpf[itime])
#      idx.close()

######### 1 h APCP

    if qpf_interval == 1:

      file0,file1,file2,file3,file4,file5,file6,file7,file8 = memfiles[itime]

# Process first 1 hour
      print 'Processing member',(1+memcount),'of',nm_use
      print 'fhours of mem: ', fhours[memcount]
      fhour=fhours[memcount]

      if mem == 'nam' :
        print 'call process_nam_qpf with: '
        print 'file1: ', file1
        print 'file0: ', file0
        print 'fcst_hour: ', fcst_hour
        fcst_hour=fhours[memcount]
        qpf[itime]=process_nam_qpf(file1,file0,fcst_hour)

      else:

        fhour=fhours[memcount]
        shour=fhour-1
        print 'from memfiles file1 for 1 h qpf: ', file1
#        idx = pygrib.index(file1,'name','lengthOfTimeRange')
        os.system(WGRIB2+' '+file1+' -match "APCP:surface:%i'%shour+'-%i'%fhour+'" -end -text qpf.txt ')
        qpf1,nx,ny=simplewgrib2('qpf.txt')
        qpf[itime] = qpf1

      print 'max of 1 h APCP: ', np.max(qpf[itime])
      print 'mean of 1 h APCP: ', np.mean(qpf[itime])
      print 'median of 1 h APCP: ', np.median(qpf[itime])




###############

    print 'here past qpf_interval tests'
    print 'here with itime: ', itime
    if dom == 'conusavoid' or dom == 'ak':
      qpf[itime] = np.where(np.equal(maskregion,-9999),0,qpf[itime])

# prob is okay as is based on local qpf[itime]

# but save to qpf[memcount] for use below

    qpf[memcount]=qpf[itime]

    # calculate threshold exceedance on QPF
    for t in thresh_use:
      for size in rlist:
        if memcount == 0:
          prob[t,size] = np.zeros((nlats,nlons))
        exceed = np.where(np.greater_equal(qpf[itime]/25.4,t),1,0)
        dontexceed = np.where(np.less(qpf[itime]/25.4,t),1,0)
        filter_footprint = get_footprint(size)

#        print 'np.shape(filter_footprint): ', np.shape(filter_footprint)
#        print 'np.shape(exceed): ', np.shape(exceed)
## exceed varies by member, so the fftconvolve will be sensitive to the order in which the members are done
        prob[t,size] = prob[t,size] + signal.fftconvolve(exceed,filter_footprint,mode='same')
        print 't, size, sum(exceed),sum(dontexceed), mean prob: ', t, size, np.sum(exceed),np.sum(dontexceed), np.mean(prob[t,size])

# possible to compute number of valid points here?

    print 'working memcount: ', memcount

    memcount = memcount + 1

  latency = min_latency


#-------------------------------------------------------------------#

# Get final probabilities
# redefine number of members (in case a TL member is missing)
nm_use = len(itimes)

print 'settled on nm_use for final probabilities: ', nm_use

ensemble_qpf = np.zeros((nm_use,nlats,nlons)).astype(float)

for mem in range(0,len(itimes)):
  print 'mem, itimes[mem],max(qpf): ',mem, itimes[mem], np.max(qpf[mem])
  ensemble_qpf[mem,:,:] = qpf[mem]
  print 'max of member: ', np.max(ensemble_qpf[mem,:,:])

# Get final probabilities
probfinal = np.zeros((nlats,nlons))
filter_footprint_10 = get_footprint(10)
filter_footprint_25 = get_footprint(25)
filter_footprint_40 = get_footprint(40)
filter_footprint_55 = get_footprint(55)
filter_footprint_70 = get_footprint(70)
filter_footprint_85 = get_footprint(85)
filter_footprint_100 = get_footprint(100)

print 'sums of the filters (10,25,40): ', np.sum(filter_footprint_10),np.sum(filter_footprint_25),np.sum(filter_footprint_40)
print 'sums of the filters (55,70,85,100): ', np.sum(filter_footprint_55),np.sum(filter_footprint_70),np.sum(filter_footprint_85),np.sum(filter_footprint_100)

for t in thresh_use:
  t3 = time.time()
  optrad = calculate_eas_probability(ensemble_qpf,t,rlist,alpha,dx,p_smooth)
  t4 = time.time()
  print 'Time for optrad routine:', t4-t3
  pnt_prob = calculate_pnt_probability (ensemble_qpf, t, p_smooth)


# how know the filter_footprint size for near boundary points?


## need something to account for reduced portion of filter_footprint actually within domain
#
#
  print 'nm_use for final prob: ', nm_use
  for row in range(nlats):
    for column in range(nlons):
      rad = (optrad[row,column]).astype(int)


# assign the default radius footprint size to footprint_use

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
              print 'column, row, reduced, orig: ', column, row, footprint_use, footprint_orig
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
        if (probfinal[row,column] > 100.1): 
           print 'row, column, probfinal[row,column]: ', row, column, probfinal[row,column]
           print 'prob[t,10][row,column], footprint_use: ', prob[t,10][row,column], footprint_use
      elif (17.5 <= rad < 32.5):
        probfinal[row,column] = prob[t,25][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
        if (probfinal[row,column] > 100.1): 
           print 'row, column, probfinal[row,column]: ', row, column, probfinal[row,column]
           print 'prob[t,25][row,column], footprint_use: ', prob[t,25][row,column], footprint_use,np.sum(filter_footprint_25)
      elif (32.5 <= rad < 47.5):
        probfinal[row,column] = prob[t,40][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
        if (probfinal[row,column] > 100.1): 
           print 'row, column, probfinal[row,column]: ', row, column, probfinal[row,column]
           print 'prob[t,40][row,column], footprint_use: ', prob[t,40][row,column], footprint_use
      elif (47.5 <= rad < 62.5):
        probfinal[row,column] = prob[t,55][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
        if (probfinal[row,column] > 100.1): 
           print 'row, column, probfinal[row,column]: ', row, column, probfinal[row,column]
           print 'prob[t,55][row,column], footprint_use: ', prob[t,55][row,column], footprint_use
      elif (62.5 <= rad < 77.5):
        probfinal[row,column] = prob[t,70][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
        if (probfinal[row,column] > 100.1): 
           print 'row, column, probfinal[row,column]: ', row, column, probfinal[row,column]
           print 'prob[t,70][row,column], footprint_use: ', prob[t,70][row,column], footprint_use
      elif (77.5 <= rad < 92.5):
        probfinal[row,column] = prob[t,85][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
        if (probfinal[row,column] > 100.1): 
           print 'row, column, probfinal[row,column]: ', row, column, probfinal[row,column]
           print 'prob[t,85][row,column], footprint_use: ', prob[t,85][row,column], footprint_use
      elif (92.5 <= rad <= 100):
        probfinal[row,column] = prob[t,100][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
        if (probfinal[row,column] > 100.1): 
           print 'row, column, probfinal[row,column]: ', row, column, probfinal[row,column]
           print 'prob[t,100][row,column], footprint_use: ', prob[t,100][row,column], footprint_use
      elif (rad > 100):
        probfinal[row,column] = prob[t,100][row,column]
        probfinal[row,column] = 100.0*probfinal[row,column] / float(footprint_use*nm_use)
        optrad[row,column] = 0

# slight smoothing of probfinal
  probfinal = ndimage.filters.gaussian_filter(probfinal,1)

  if dom == 'conusavoid' or dom == 'ak':
    probfinal = np.where(np.equal(maskregion,-9999),0,probfinal)  # set to 0 for mask 

  t5 = time.time()

  print 'Time for get final probability routine for ',t, 'inch threshold: ',t5-t4
  print 'max of probfinal: ', np.max(probfinal)
  print 'mean of probfinal: ', np.mean(probfinal)
  print 'probfinal dims: ', np.shape(probfinal)


  probstr=str(t*25.4)
  byte=int(m.ceil(t*25.4*1000))
  byte44=0
  byte45=int(byte/65536)
  byte45rem=byte%65536
  byte46=int(byte45rem/256)
  byte47=byte45rem%256

  myfort = F.FortranFile('record_out.bin',mode='w')
  myfort.writeReals(probfinal)
  myfort.close()

  string="0:0:d="+wgribdate+":APCP:surface:"+fhr_range+" hour acc fcst:prob >"+probstr+":"

  os.system(WGRIB2+' '+template+' -import_bin record_out.bin -set_metadata_str "'+string+'" -set_grib_type c3 -grib_out premod.grb')
  os.system(WGRIB2+' premod.grb -set_byte 4 12 197 -set_byte 4 17 0 -set_byte 4 24:35 0:0:0:0:0:255:0:0:0:0:0:0 -set_byte 4 36 '+str(nm_use)+' -set_byte 4 38:42 0:0:0:0:0 -set_byte 4 43 3 -set_byte 4 44 0 -set_byte 4 45 '+str(byte45)+' -set_byte 4 46 '+str(byte46)+' -set_byte 4 47 '+str(byte47)+' -set_byte 4 67 1 -append  -set_grib_type c3 -grib_out '+outfile)
  print 'Wrote ', qpf_interval, ' PQPF to:',outfile, 'for ',t, 'inch threshold'
  os.system('rm record_out.bin')
  os.system('rm premod.grb')
