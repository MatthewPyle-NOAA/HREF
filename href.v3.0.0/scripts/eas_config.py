# SSEOx CONFIGURATION FILE

###########################
# SETTINGS FOR ALL FIELDS #
###########################

# number of ens members
nm = 9

nm_v3 = 10
nm_ak = 8
nm_nonconus = 6
# model resolution (km)
dx = 3.0
# min latency (hours)
min_latency = 0
# max latency (hours)
max_latency = 12
max_latency_nam = 6

#################
# PQPF SETTINGS #
#################

# exceedance thresholds (in)

# possible optimal set?
pqpf_1h_thresh = [.01,0.25,0.5]
pqpf_3h_thresh = [.01,0.25,0.5]
pqpf_6h_thresh = [.01,0.25,0.5,1.0]
pqpf_12h_thresh = [0.1,0.25,0.5,1.0,2.0]
pqpf_24h_thresh = [0.1,0.25,0.5,1.0,2.0,3.0]

# neighborhood size (km)
pqpf_neighborhood = 40

# spatial filter radius (km)
pqpf_1h_filter = 60
pqpf_3h_filter = 100
pqpf_6h_filter = 100
rlist = [10,25,40,55,70,85,100]
rlist_hi = [10,25,40,55,70]

# EAS settings
alpha = 0.1	# similarity criteria parameter
p_smooth_low = 3	# width of Gaussian filter for smoothing radius field (grid points) 
p_smooth = 3	# width of Gaussian filter for smoothing radius field (grid points) 
p_smooth_high = 3	# width of Gaussian filter for smoothing radius field (grid points) 

# forecast hour weighting
pqpf_1h_weight1 = 1.0  # hour before
pqpf_1h_weight2 = 1.0  # valid hour
pqpf_1h_weight3 = 1.0  # hour after

# PQPF calibration
pqpf_6h_calibrate_arw = 'no'
pqpf_6h_calibrate_nmmb = 'no'
pqpf_6h_calibrate_fv3 = 'no'
pqpf_6h_calibrate_arw2 = 'no'
pqpf_6h_calibrate_nam = 'no'
pqpf_6h_calibrate_hrrr = 'no'
pqpf_3h_calibrate_arw = 'no'
pqpf_3h_calibrate_nmmb = 'no'
pqpf_3h_calibrate_arw2 = 'no'
pqpf_3h_calibrate_nam = 'no'
pqpf_1h_calibrate_arw = 'no'
pqpf_1h_calibrate_nmmb = 'no'
pqpf_1h_calibrate_arw2 = 'no'
pqpf_1h_calibrate_nam = 'no'
# minimum quantile for curve fitting
pqpf_6h_minpct = 75
pqpf_3h_minpct = 75
pqpf_1h_minpct = 75


########################
# FLASH-FLOOD SETTINGS #
########################

# ARI uses spatial filter from 6-h PQPF
# both products use pqpf neighborhood

# exceedance thresholds (in)
runoff_3h_thresh = [1.0,2.0,3.0]

# recurrence interval exceedance threshold
ari_years = 100

# percent of this exceedance threshold for a "hit"
ari_pct = 100

# spatial filter radius (km)
runoff_3h_filter = 100


#################
# SNOW SETTINGS #
#################

# exceedance thresholds (in)

 # do 1h 3" as neighborhood max
snow_1h_thresh = [1.0,3.0]
 # do 3h 3" and 6" as neighborhood max
snow_3h_thresh = [1.0,3.0]
 # do 6h 6" and 12" as neighborhood max
snow_6h_thresh = [1.0,3.0,6.0]

# neighborhood size (km)
snow_neighborhood = 40

# spatial filter radius (km)
snow_1h_filter = 60
snow_3h_filter = 100
snow_6h_filter = 100

# forecast hour weighting
snow_1h_weight1 = 1.0  # hour before
snow_1h_weight2 = 1.0  # valid hour
snow_1h_weight3 = 1.0  # hour after


