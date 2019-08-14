"""This is code to generate a land mask for the CORE forcing files
It also zonally averages and hemispherically symmetrizes the files for use in ocean-only aquaplanet simulations."""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc4
import scipy.interpolate as sci
import shutil
import os
import glob

'''Copy the original CORE forcing files from source directory, CORE_cnyf, to source directory, CORE_new'''

source_dir = '/Users/sragen/PycharmProjects/AQUA/DATA/CORE_cnyf'
dest_dir = '/Users/sragen/PycharmProjects/AQUA/DATA/CORE_new'
for filename in glob.glob(os.path.join(source_dir, '*.*')):
    shutil.copy(filename, dest_dir)

'''First, we use some model output to create a land mask. Then we interpolate the model mask and the CORE forcing data 
to create a new land mask that has the same dimensions as the CORE forcing files.
This new CORE mask can then be used to determine where there is land and where there is ocean in all of the CORE files 
used to force an ocean-only model.'''

# load data to use as land mask
f = xr.open_dataset('CORE_new/newCO2_control_2d.2399-2499.ocean_month.clim.nc', decode_times=False)
sst = np.array(f['sst'][:])
y_mask = np.array(f['yt_ocean'][:])
x_mask = np.array(f['xt_ocean'][:])
# sort the data so the longitude matches that of the CORE files
inds = np.where(x_mask[x_mask < 0])[0]
x_mask[inds] = x_mask[inds] + 360
i = np.argsort(x_mask)
x_mask = x_mask[i]
mask = sst[0,:,i].T

# load CORE temperature data
f0 = xr.open_dataset('CORE_new/t_10.15JUNE2009.nc')
x = np.array(f0['LON'][:])
y = np.array(f0['LAT'][:])
T_10 = np.array(f0['T_10'][:])

# Now, interpolate the mask grid onto the CORE grid
X_mask, Y_mask = np.meshgrid(x_mask, y_mask)
X, Y = np.meshgrid(x, y)
core_mask = sci.griddata((np.ravel(X_mask),np.ravel(Y_mask)), np.ravel(mask), (X,Y), method='linear')
######################################

'''Now that we have a new CORE land mask, we use it to ignore points over land in the CORE files to create new CORE 
forcing files. We don't want to see the direct influence of land in our forcing for ocean-only model runs.
We will also symmetrize the southern hemisphere zonally and mirror it across the equator to replace the northern 
hemisphere. Thus we will have hemispherically and zonally symmetric forcings.'''

#######################################
# load files:
f1 = nc4.Dataset('CORE_new/ncar_precip.15JUNE2009.nc', 'r+')
f2 = nc4.Dataset('CORE_new/ncar_rad.15JUNE2009.nc', 'r+')
f3 = nc4.Dataset('CORE_new/q_10.15JUNE2009.nc', 'r+')
f4 = nc4.Dataset('CORE_new/slp.15JUNE2009.nc', 'r+')
f5 = nc4.Dataset('CORE_new/t_10.15JUNE2009.nc', 'r+')
f6 = nc4.Dataset('CORE_new/u_10.15JUNE2009.nc', 'r+')
f7 = nc4.Dataset('CORE_new/v_10.15JUNE2009.nc', 'r+')
######################################

# File 1: ncar_precip.15JUNE2009.nc
glob_prec_adj1 = np.array(f1.variables['GLOB_PREC_ADJ1'][:])
glob_prec_adj1_mask = np.where(~np.isnan(core_mask), glob_prec_adj1, np.nan)
shem = np.nanmean(glob_prec_adj1_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
glob_prec_adj1_mask_sym = np.concatenate((shem,nhem), axis=1)
GLOB_PREC_ADJ1_mask_sym = np.tile(glob_prec_adj1_mask_sym[:,:, np.newaxis], (1, 1, 192))
f1.variables['GLOB_PREC_ADJ1'][:] = GLOB_PREC_ADJ1_mask_sym

prc_mod = np.array(f1.variables['PRC_MOD'][:])
prc_mod_mask = np.where(~np.isnan(core_mask), prc_mod, np.nan)
shem = np.nanmean(prc_mod_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
prc_mod_mask_sym = np.concatenate((shem,nhem), axis=1)
PRC_MOD_mask_sym = np.tile(prc_mod_mask_sym[:,:, np.newaxis], (1, 1, 192))
f1.variables['PRC_MOD'][:] = PRC_MOD_mask_sym

prc_mod1 = np.array(f1.variables['PRC_MOD1'][:])
prc_mod1_mask = np.where(~np.isnan(core_mask), prc_mod1, np.nan)
shem = np.nanmean(prc_mod1_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
prc_mod1_mask_sym = np.concatenate((shem,nhem), axis=1)
PRC_MOD1_mask_sym = np.tile(prc_mod1_mask_sym[:,:, np.newaxis], (1, 1, 192))
f1.variables['PRC_MOD1'][:] = PRC_MOD1_mask_sym

prec_raw = np.array(f1.variables['PREC_RAW'][:])
prec_raw_mask = np.where(~np.isnan(core_mask), prec_raw, np.nan)
shem = np.nanmean(prec_raw_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
prec_raw_mask_sym = np.concatenate((shem,nhem), axis=1)
PREC_RAW_mask_sym = np.tile(prec_raw_mask_sym[:,:, np.newaxis], (1, 1, 192))
f1.variables['PREC_RAW'][:] = PREC_RAW_mask_sym

rain = np.array(f1.variables['RAIN'][:])
rain_mask = np.where(~np.isnan(core_mask), rain, np.nan)
shem = np.nanmean(rain_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
rain_mask_sym = np.concatenate((shem,nhem), axis=1)
RAIN_mask_sym = np.tile(rain_mask_sym[:,:, np.newaxis], (1, 1, 192))
f1.variables['RAIN'][:] = RAIN_mask_sym

snow = np.array(f1.variables['SNOW'][:])
snow_mask = np.where(~np.isnan(core_mask), snow, np.nan)
shem = np.nanmean(snow_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
snow_mask_sym = np.concatenate((shem,nhem), axis=1)
SNOW_mask_sym = np.tile(snow_mask_sym[:,:, np.newaxis], (1, 1, 192))
f1.variables['SNOW'][:] = SNOW_mask_sym

# File 2: ncar_rad.15JUNE2009.nc
lwdn = np.array(f2.variables['LWDN'][:])
lwdn_mask = np.where(~np.isnan(core_mask), lwdn, np.nan)
shem = np.nanmean(lwdn_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
lwdn_mask_sym = np.concatenate((shem,nhem), axis=1)
LWDN_mask_sym = np.tile(lwdn_mask_sym[:,:, np.newaxis], (1, 1, 192))
f2.variables['LWDN'][:] = LWDN_mask_sym

lwdn_mod = np.array(f2.variables['LWDN_MOD'][:])
lwdn_mod_mask = np.where(~np.isnan(core_mask), lwdn_mod, np.nan)
shem = np.nanmean(lwdn_mod_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
lwdn_mod_mask_sym = np.concatenate((shem,nhem), axis=1)
LWDN_MOD_mask_sym = np.tile(lwdn_mod_mask_sym[:,:, np.newaxis], (1, 1, 192))
f2.variables['LWDN_MOD'][:] = LWDN_MOD_mask_sym

swdn = np.array(f2.variables['SWDN'][:])
swdn_mask = np.where(~np.isnan(core_mask), swdn, np.nan)
shem = np.nanmean(swdn_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
swdn_mask_sym = np.concatenate((shem,nhem), axis=1)
SWDN_mask_sym = np.tile(swdn_mask_sym[:,:, np.newaxis], (1, 1, 192))
f2.variables['SWDN'][:] = SWDN_mask_sym

swdn_mod = np.array(f2.variables['SWDN_MOD'][:])
swdn_mod_mask = np.where(~np.isnan(core_mask), swdn_mod, np.nan)
shem = np.nanmean(swdn_mod_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
swdn_mod_mask_sym = np.concatenate((shem,nhem), axis=1)
SWDN_MOD_mask_sym = np.tile(swdn_mod_mask_sym[:,:, np.newaxis], (1, 1, 192))
f2.variables['SWDN_MOD'][:] = SWDN_MOD_mask_sym

# File 3: q_10.15JUNE2009.nc
q_10 = np.array(f3.variables['Q_10'][:])
q_10_mask = np.where(~np.isnan(core_mask), q_10, np.nan)
shem = np.nanmean(q_10_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
q_10_mask_sym = np.concatenate((shem,nhem), axis=1)
Q_10_mask_sym = np.tile(q_10_mask_sym[:,:, np.newaxis], (1, 1, 192))
f3.variables['Q_10'][:] = Q_10_mask_sym

q_10_mod = np.array(f3.variables['Q_10_MOD'][:])
q_10_mod_mask = np.where(~np.isnan(core_mask), q_10_mod, np.nan)
shem = np.nanmean(q_10_mod_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
q_10_mod_mask_sym = np.concatenate((shem,nhem), axis=1)
Q_10_MOD_mask_sym = np.tile(q_10_mod_mask_sym[:,:, np.newaxis], (1, 1, 192))
f3.variables['Q_10_MOD'][:] = Q_10_MOD_mask_sym

# File 4: slp.15JUNE2009.nc
# Just set sea level pressure to a constant value of 101000 everywhere
slp = np.array(f4.variables['SLP'][:])*0
SLP_mask_sym = slp + 101000
f4.variables['SLP'][:] = SLP_mask_sym

# File 5: t_10.15JUNE2009.nc
t_10 = np.array(f5.variables['T_10'][:])
t_10_mask = np.where(~np.isnan(core_mask), t_10, np.nan)
shem = np.nanmean(t_10_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
t_10_mask_sym = np.concatenate((shem,nhem), axis=1)
T_10_mask_sym = np.tile(t_10_mask_sym[:,:, np.newaxis], (1, 1, 192))
f5.variables['T_10'][:] = T_10_mask_sym

t_10_mod = np.array(f5.variables['T_10_MOD'][:])
t_10_mod_mask = np.where(~np.isnan(core_mask), t_10_mod, np.nan)
shem = np.nanmean(t_10_mod_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
t_10_mod_mask_sym = np.concatenate((shem,nhem), axis=1)
T_10_MOD_mask_sym = np.tile(t_10_mod_mask_sym[:,:, np.newaxis], (1, 1, 192))
f5.variables['T_10_MOD'][:] = T_10_MOD_mask_sym

# File 6: u_10.15JUNE2009.nc
u_10 = np.array(f6.variables['U_10'][:])
u_10_mask = np.where(~np.isnan(core_mask), u_10, np.nan)
shem = np.nanmean(u_10_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
u_10_mask_sym = np.concatenate((shem,nhem), axis=1)
U_10_mask_sym = np.tile(u_10_mask_sym[:,:, np.newaxis], (1, 1, 192))
f6.variables['U_10'][:] = U_10_mask_sym

u_10_mod = np.array(f6.variables['U_10_MOD'][:])
u_10_mod_mask = np.where(~np.isnan(core_mask), u_10_mod, np.nan)
shem = np.nanmean(u_10_mod_mask[:,:47,:], axis=2)
tmp = shem[:,::-1]
mid = int(len(shem[:,0])/2)
nhem = np.concatenate((tmp[mid:,:],tmp[:mid,:]), axis=0)
u_10_mod_mask_sym = np.concatenate((shem,nhem), axis=1)
U_10_MOD_mask_sym = np.tile(u_10_mod_mask_sym[:,:, np.newaxis], (1, 1, 192))
f6.variables['U_10_MOD'][:] = U_10_MOD_mask_sym

# File 7: v_10.15JUNE2009.nc
# Just set meridional wind, v, to be 0 everywhere
v_10 = np.array(f7.variables['V_10'][:])*0 
V_10_mask_sym = v_10
f7.variables['V_10'][:] = V_10_mask_sym

v_10_mod = np.array(f7['V_10_MOD'][:])*0
V_10_MOD_mask_sym = v_10_mod
f7.variables['V_10_MOD'][:] = V_10_MOD_mask_sym

###############################
# Uncomment this block if you want to see plots of the masked and symmetrized T_10 variable.
'''cmap = plt.get_cmap('magma')
for i in range(0, len(V_10_mask_sym[:,0,0]), 200):
    fig, ax = plt.subplots()
    plt.contourf(f5.variables['LON'][:], f5.variables['LAT'][:], T_10_mask_sym[i,:,:], cmap=cmap, extend='both')
    ax.set_title('T 10 {0}'.format(i))
    cbar = plt.colorbar()
    plt.show(block=False)'''

###############################
'''Now, we need to extend the variables to cover the entire latitudinal extent of the model.'''

# File 1: ncar_precip.15JUNE2009.nc
glob_prec_adj1 = np.array(f1.variables['GLOB_PREC_ADJ1'][:])
data = glob_prec_adj1[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
glob_prec_adj1 = np.concatenate((conc1,data,conc2), axis=1)
f1.variables['GLOB_PREC_ADJ1'][:] = glob_prec_adj1

prc_mod = np.array(f1.variables['PRC_MOD'][:])
data = prc_mod[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
prc_mod = np.concatenate((conc1,data,conc2), axis=1)
f1.variables['PRC_MOD'][:] = prc_mod

prc_mod1 = np.array(f1.variables['PRC_MOD1'][:])
data = prc_mod1[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
prc_mod1 = np.concatenate((conc1,data,conc2), axis=1)
f1.variables['PRC_MOD1'][:] = prc_mod1

prec_raw = np.array(f1.variables['PREC_RAW'][:])
data = prec_raw[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
prec_raw = np.concatenate((conc1,data,conc2), axis=1)
f1.variables['PREC_RAW'][:] = prec_raw

rain = np.array(f1.variables['RAIN'][:])
data = rain[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
rain = np.concatenate((conc1,data,conc2), axis=1)
f1.variables['RAIN'][:] = rain

snow = np.array(f1.variables['SNOW'][:])
data = snow[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
snow = np.concatenate((conc1,data,conc2), axis=1)
f1.variables['SNOW'][:] = snow

# File 2: ncar_rad.15JUNE2009.nc
lwdn = np.array(f2.variables['LWDN'][:])
data = lwdn[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
lwdn = np.concatenate((conc1,data,conc2), axis=1)
f2.variables['LWDN'][:] = lwdn

lwdn_mod = np.array(f2.variables['LWDN_MOD'][:])
data = lwdn_mod[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
lwdn_mod = np.concatenate((conc1,data,conc2), axis=1)
f2.variables['LWDN_MOD'][:] = lwdn_mod

swdn = np.array(f2.variables['SWDN'][:])
data = swdn[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
swdn = np.concatenate((conc1,data,conc2), axis=1)
f2.variables['SWDN'][:] = swdn

swdn_mod = np.array(f2.variables['SWDN_MOD'][:])
data = swdn_mod[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
swdn_mod = np.concatenate((conc1,data,conc2), axis=1)
f2.variables['SWDN_MOD'][:] = swdn_mod

# File 3: q_10.15JUNE2009.nc
q_10 = np.array(f3.variables['Q_10'][:])
data = q_10[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
q_10 = np.concatenate((conc1,data,conc2), axis=1)
f3.variables['Q_10'][:] = q_10

q_10_mod = np.array(f3.variables['Q_10_MOD'][:])
data = q_10_mod[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
q_10_mod = np.concatenate((conc1,data,conc2), axis=1)
f3.variables['Q_10_MOD'][:] = q_10_mod

# File 5: t_10.15JUNE2009.nc
t_10 = np.array(f5.variables['T_10'][:])
data = t_10[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
t_10 = np.concatenate((conc1,data,conc2), axis=1)
f5.variables['T_10'][:] = t_10

t_10_mod = np.array(f5.variables['T_10_MOD'][:])
data = t_10_mod
data = data[:,7:-7,:]
conc = np.expand_dims(data[:,0,:], axis=1)
conc = np.tile(conc, (1,7,1))
t_10_mod = np.concatenate((conc,data,conc), axis=1)
f5.variables['T_10_MOD'][:] = t_10_mod

# File 6: u_10.15JUNE2009.nc
u_10 = np.array(f6.variables['U_10'][:])
data = u_10[:,7:-7,:]
conc1 = np.expand_dims(data[:,0,:], axis=1)
conc1 = np.tile(conc1, (1,7,1))
conc2 = np.expand_dims(data[:,-1,:], axis=1)
conc2 = np.tile(conc2, (1,7,1))
u_10 = np.concatenate((conc1,data,conc2), axis=1)
f6.variables['U_10'][:] = u_10

u_10_mod = np.array(f6.variables['U_10_MOD'][:])
data = u_10_mod
data = data[:,7:-7,:]
conc = np.expand_dims(data[:,0,:], axis=1)
conc = np.tile(conc, (1,7,1))
u_10_mod = np.concatenate((conc,data,conc), axis=1)
f6.variables['U_10_MOD'][:] = u_10_mod

###############################
# close files
f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()
f7.close()
################################

f5 = xr.open_dataset('CORE_new/t_10.15JUNE2009.nc')
lon = f5['LON']
lat = f5['LAT']
time = f5['TIME']
t_10 = f5['T_10']

b = np.linspace(240, 310, 11)
t = np.linspace(240, 310, 21)
cmap = plt.get_cmap('magma')
for i in range(0, len(time), 200):
    fig, ax = plt.subplots()
    plt.contourf(lon, lat, t_10[i,:,:], t, cmap=cmap, extend='both')
    plt.title('T 10 {0}'.format(i))
    cbar = plt.colorbar(ticks=b, boundaries=b, spacing='uniform')
    plt.show(block=False)

fig, ax = plt.subplots()
for i in range(0, len(time), 200):
    plt.plot(lat, t_10[i,:,0])
    plt.title('T_10')
    plt.show(block=False)
