#playing around with Ocean Parcels
#before running, make sure parcels is activated using 'conda activate py3_parcels'

import os
import numpy as np
from datetime import timedelta
import random
import glob
import datetime
import netCDF4
from parcels import (FieldSet, Field, NestedField, ParticleSet, JITParticle,
    AdvectionRK4_3D, ErrorCode, plotTrajectoriesFile, Variable)
from lanaau_build_spawn import getSpawningSites, getSpawningEvents
from lanaau_functions import (HandleBoundaries, AdvectionRK4_3D_Lanaau,
    RWDiffusion_Lanaau, deleteParticle, submergeParticle,
    Buoyancy, KillFish, SettleFish, GrowFish)
from LarvalParticle import define_Larval_Class

def setFields_nonest(dt):

    dim = {'lon':'Longitude_u', 'lat':'Latitude_t', 'depth':'Depth_w'} #u, t works for 1km; t,t works for 4km
    ts = np.expand_dims(np.arange(0, 73526400, 86400), axis = 1) #set up timestamps (corresponds to day of model)
        
    #combine everything into nested fieldset
    fieldset_nested = FieldSet.from_mitgcm(glob.glob("/10tb_abyss/emily/biophys_modeling/current_products/new_1km/*cdf"),
    variables = {'U':'u', 'V':'v', 'W':'w'}, dimensions=dim, timestamps=ts, deferred_load=True)

    #add fields for coastlines
    fieldset_nested.add_constant("dt", dt)
    northCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/northcoasts_1km.nc", 
            variable=("F_N", "northcoast_1km"), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    southCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/southcoasts_1km.nc", 
            variable=("F_S", "southcoast_1km"), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    eastCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/eastcoasts_1km.nc", 
            variable=("F_E", "eastcoast_1km"), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    westCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/westcoasts_1km.nc", 
            variable=("F_W", "westcoast_1km"), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    
    for coast in [northCoasts_1km, southCoasts_1km, eastCoasts_1km, westCoasts_1km]:
        fieldset_nested.add_field(coast)
        
    #set up field to check interpolation
    F = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_1km/1km_mask.nc", 
            variable=('F','F1'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    fieldset_nested.add_field(F)
    
    #set up bathymetry field
    bathy = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_1km/1km_bathy.nc", 
            variable=('seafloor','depth1km'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    fieldset_nested.add_field(bathy)
    fieldset_nested.seafloor.interp_method = "nearest"

    return(fieldset_nested)

def setFields_from_fieldset(dt):
    #trying nested interpolation fields from fieldsets
    #extra step but reading in just fieldsets (setFields) doesn't seem to work
    dim_u = {'lon':'Longitude_u', 'lat':'Latitude_t', 'depth':'Depth_t'}
    dim_v = {'lon':'Longitude_t', 'lat':'Latitude_v', 'depth':'Depth_t'}
    dim_w = {'lon':'Longitude_t', 'lat':'Latitude_t', 'depth':'Depth_w'}
    ts = np.expand_dims(np.arange(0, 73526400, 86400), axis = 1) #set up timestamps (corresponds to day of model in seconds)

    fieldset_1km = FieldSet.from_mitgcm(glob.glob("/10tb_abyss/emily/biophys_modeling/current_products/new_1km/*cdf"),
    variables = {'U1':'u', 'V1':'v', 'W1':'w'}, dimensions = {'U1':dim_u, 'V1':dim_v, 'W1':dim_w}, 
    timestamps=ts, deferred_load=True)
    
    fieldset_4km = FieldSet.from_mitgcm(glob.glob("/10tb_abyss/emily/biophys_modeling/current_products/final_4km/*cdf"),
    variables = {'U2':'u', 'V2':'v', 'W2':'w'}, dimensions = {'U2':dim_u, 'V2':dim_v, 'W2':dim_w}, 
    timestamps=ts, deferred_load=True)
    
    U_nested = NestedField('U', [fieldset_1km.U1, fieldset_4km.U2])
    V_nested = NestedField('V', [fieldset_1km.V1, fieldset_4km.V2])
    W_nested = NestedField('W', [fieldset_1km.W1, fieldset_4km.W2])
    
    #combine everything into nested fieldset
    fieldset_nested = fieldset_4km
    fieldset_nested.add_field(U_nested)
    fieldset_nested.add_field(V_nested)
    fieldset_nested.add_field(W_nested)
    
    #add separate fine fields for silliness near coasts
    fieldset_nested.add_field(fieldset_1km.U1)
    fieldset_nested.add_field(fieldset_1km.V1)
    fieldset_nested.add_field(fieldset_1km.W1)
    
    #set up field to check interpolation
    F1 = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_1km/1km_mask.nc", 
            variable=('F1','F1'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    F2 = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_4km/4km_mask.nc", 
            variable=('F2','F2'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    F = NestedField('F', [F1, F2])
    fieldset_nested.add_field(F)
    
    #set up bathymetry field
    B1 = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_1km/1km_bathy.nc", 
            variable='depth1km', dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    B2 = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_4km/4km_bathy.nc", 
            variable='depth4km', dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    bathy = NestedField('seafloor', [B1, B2])
    fieldset_nested.add_field(bathy)
    fieldset_nested.seafloor.interp_method = "nearest"

    #add fields for coastlines
    fieldset_nested.add_constant("dt", dt)
    northCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/northcoasts_1km.nc", 
            variable="northcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    southCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/southcoasts_1km.nc", 
            variable="southcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    eastCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/eastcoasts_1km.nc", 
            variable="eastcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    westCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/westcoasts_1km.nc", 
            variable="westcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    
    northCoasts_4km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/northcoasts_4km.nc", 
            variable="northcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    southCoasts_4km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/southcoasts_4km.nc", 
            variable="southcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    eastCoasts_4km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/eastcoasts_4km.nc", 
            variable="eastcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    westCoasts_4km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/westcoasts_4km.nc", 
            variable="westcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)

    fieldset_nested.add_field(NestedField('F_N', [northCoasts_1km, northCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_S', [southCoasts_1km, southCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_E', [eastCoasts_1km, eastCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_W', [westCoasts_1km, westCoasts_4km]))

    return(fieldset_nested)
    
def setFields(dt):
    #set up nested interpolation fields
    #bounds of 4km field: 190, 209.5141; 15.42, 27
    #bounds of 1km field: 198.695, 205.4807; 17.495, 21.68512

    #var and dim the same for both fieldsets
    dim_u = {'lon':'Longitude_u', 'lat':'Latitude_t', 'depth':'Depth_w'}
    dim_v = {'lon':'Longitude_t', 'lat':'Latitude_v', 'depth':'Depth_w'}
    dim_w = {'lon':'Longitude_t', 'lat':'Latitude_t', 'depth':'Depth_w'}
    ts = np.expand_dims(np.arange(0, 73526400, 86400), axis = 1) #set up timestamps (corresponds to day of model in seconds)
        
    U_fine = Field.from_netcdf(glob.glob("/10tb_abyss/emily/biophys_modeling/current_products/new_1km/*cdf"), 
            variable=('U1','u'), dimensions=dim_u, timestamps=ts, deferred_load=True, gridindexingtype = "mitgcm", fieldtype = "U")
    V_fine = Field.from_netcdf(glob.glob("/10tb_abyss/emily/biophys_modeling/current_products/new_1km/*cdf"),
            variable=('V1','v'), dimensions=dim_v, timestamps=ts, deferred_load=True, gridindexingtype = "mitgcm", fieldtype = "V") 
    W_fine = Field.from_netcdf(glob.glob("/10tb_abyss/emily/biophys_modeling/current_products/new_1km/*cdf"),
            variable=('W1','w'), dimensions=dim_w, timestamps=ts, deferred_load=True, gridindexingtype = "mitgcm")
            
    U_coarse = Field.from_netcdf(glob.glob("/10tb_abyss/emily/biophys_modeling/current_products/final_4km/*cdf"),
            variable=('U2','u'), dimensions=dim_u,  timestamps=ts, deferred_load=True, gridindexingtype = "mitgcm", fieldtype = "U") 
    V_coarse = Field.from_netcdf(glob.glob("/10tb_abyss/emily/biophys_modeling/current_products/final_4km/*cdf"),
            variable=('V2','v'), dimensions=dim_v, timestamps=ts, deferred_load=True, gridindexingtype = "mitgcm", fieldtype = "V")
    W_coarse = Field.from_netcdf(glob.glob("/10tb_abyss/emily/biophys_modeling/current_products/final_4km/*cdf"),
            variable=('W2','w'), dimensions=dim_w, timestamps=ts, deferred_load=True, gridindexingtype = "mitgcm")

    #nest fields
    U_nested = NestedField('U', [U_fine, U_coarse])
    V_nested = NestedField('V', [V_fine, V_coarse])
    W_nested = NestedField('W', [W_fine, W_coarse])

    #combine everything into nested fieldset
    fieldset_nested = FieldSet(U_nested, V_nested)
    fieldset_nested.add_field(W_nested)

    #add separate fine fields for silliness near coasts
    fieldset_nested.add_field(U_fine)
    fieldset_nested.add_field(V_fine)
    fieldset_nested.add_field(W_fine)
    
    #set up field to check interpolation
    F1 = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_1km/1km_mask.nc", 
            variable=('F1','F1'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    F2 = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_4km/4km_mask.nc", 
            variable=('F2','F2'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    F = NestedField('F', [F1, F2])
    fieldset_nested.add_field(F)
    
    #set up bathymetry field
    B1 = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_1km/1km_bathy.nc", 
            variable='depth1km', dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    B2 = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_4km/4km_bathy.nc", 
            variable='depth4km', dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    bathy = NestedField('seafloor', [B1, B2])
    fieldset_nested.add_field(bathy)
    fieldset_nested.seafloor.interp_method = "nearest"

    #add fields for coastlines
    fieldset_nested.add_constant("dt", dt)
    northCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/northcoasts_1km.nc", 
            variable="northcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    southCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/southcoasts_1km.nc", 
            variable="southcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    eastCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/eastcoasts_1km.nc", 
            variable="eastcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    westCoasts_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/westcoasts_1km.nc", 
            variable="westcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    
    northCoasts_4km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/northcoasts_4km.nc", 
            variable="northcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    southCoasts_4km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/southcoasts_4km.nc", 
            variable="southcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    eastCoasts_4km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/eastcoasts_4km.nc", 
            variable="eastcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    westCoasts_4km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/coast_files/westcoasts_4km.nc", 
            variable="westcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)

    fieldset_nested.add_field(NestedField('F_N', [northCoasts_1km, northCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_S', [southCoasts_1km, southCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_E', [eastCoasts_1km, eastCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_W', [westCoasts_1km, westCoasts_4km]))

    return(fieldset_nested)

def main():
    
    #set random seed
    #666, 808
    random.seed(808)
    np.random.seed(808)
    
    #set up timestep - at 12 minutes we don't skip over any cells at max velocity
    dt = 12 #in minutes
    
    fieldset = setFields(dt) #get fieldset

    #set up species-specific parameters and add to fieldset
    typeSpawner = 0 #broadcast:0, brood:1
    monthSpawner = True #year-round (False) or seasonal reproduction (True)
    spawnMonths = [4,5]
    timedSpawner = True #spawning during a certain time of day (True) or not (False)
    timeSpawnBegin = datetime.timedelta(hours = 7)
    timeSpawnEnd = datetime.timedelta(hours = 8)
    peakSpawner = False #peak spawning during a certain time of year
    peakSpawnMonths = "NA"
    moonSpawner = True
    spawnPhase = ["Waning gibbous"]
    
    hatchTime = 1 #time til hatching in days
    flexionTime = 1 #start of swimming behavior in days
    
    #pelagic larval duration in days
    minPLD = 3
    maxPLD = 100
    
    tsMortality = (1.0/maxPLD) / 120 #assumes that by max PLD, all larvae are dead
    eggMort = tsMortality
    larvaeOneMort = tsMortality
    larvaeTwoMort = tsMortality
    
    #add buoyancy based on life stage, in ascent speed in m/s (negative towards surface)
    buoyant = True
    eggBuoy = 0
    larvaeOneBuoy = 0
    larvaeTwoBuoy = 0
    
    maxAdultDepth = 27
    maxLarvalDepth = 100
    
    #TODO: swimming?
    #behaviorStart = -1

    #add global variables to fieldset
    fieldset.add_constant("hatchTime", hatchTime) 
    #fieldset.add_constant("behaviorStart", behaviorStart) 
    fieldset.add_constant("flexionTime", flexionTime)
    fieldset.add_constant("minPLD", minPLD)
    fieldset.add_constant("maxPLD", maxPLD)
    fieldset.add_constant("maxDepth", maxLarvalDepth)
    fieldset.add_constant("eggMort", eggMort)
    fieldset.add_constant("larvaeOneMort", larvaeOneMort)
    fieldset.add_constant("larvaeTwoMort", larvaeTwoMort)
    fieldset.add_constant("eggBuoy", eggBuoy)
    fieldset.add_constant("larvaeOneBuoy", larvaeOneBuoy)
    fieldset.add_constant("larvaeTwoBuoy", larvaeTwoBuoy)
    fieldset.add_constant("lonmin", min(fieldset.U[1].grid.lon))
    fieldset.add_constant("lonmax", max(fieldset.U[1].grid.lon))
    fieldset.add_constant("latmin", min(fieldset.V[1].grid.lat))
    fieldset.add_constant("latmax", max(fieldset.V[1].grid.lat))
    
    #add settlement boundary to fieldset
    settle_1km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_1km/1km_p_meandrina_settle_final.nc", 
            variable=('km1settle','1kmsettle'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    settle_4km = Field.from_netcdf("/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_4km/4km_p_meandrina_final.nc", 
            variable=('km4settle','4kmspawnsettle'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    #fieldset.add_field(settle_1km)
    #fieldset.add_field(settle_4km)
    #fieldset.km1settle.interp_method = "nearest"
    #fieldset.km4settle.interp_method = "nearest"
    
    settle = NestedField('settle', [settle_1km, settle_4km])
    fieldset.add_field(settle)
    fieldset.settle.interp_method = "nearest"
    
    #save spawning files
    spawn_1km_file = "/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_1km/1km_p_meandrina_spawn_final.nc"
    spawn_4km_file = "/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_4km/4km_p_meandrina_final.nc"
    depth_1km_file = "/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_1km/1km_p_meandrina_depth.nc"
    depth_4km_file = "/home/emily/emily_work_dirs/2021_parcels_MHI/spawn_files/2021_final_4km/4km_p_meandrina_depth.nc"
    
    #set up spawning sites
    lat, lon, dep = getSpawningSites(spawn_1km_file, 
                                     spawn_4km_file, 
                                     depth_1km_file, 
                                     depth_4km_file, 
                                     maxAdultDepth)

    #set up spawning events (FOR TIMED + SEASONAL SPAWNER)
    starttime, events = getSpawningEvents(monthSpawner, 
                                          spawnMonths, timedSpawner, 
                                          timeSpawnBegin, timeSpawnEnd, peakSpawner, 
                                          peakSpawnMonths, moonSpawner, spawnPhase)

    #TESTING - temporary
    #idx = random.sample(range(len(lon)), 10)
    #lon = [lon[i] for i in idx]
    #lat = [lat[i] for i in idx]
    #dep = [dep[i] for i in idx]
    
    #set up particle set
    LarvalParticle = define_Larval_Class(typeSpawner)
    pset = ParticleSet.from_list(fieldset = fieldset,
                                 pclass = LarvalParticle,
                                 lon = lon, lat = lat,
                                 depth = dep, time = starttime,
                                 repeatdt = timedelta(minutes=dt).total_seconds(),
                                 siteID = range(1, (len(lon) + 1)))

    print("Number of spawning sites:")
    print(len(lon))

    #functions to kernels
    #kernels = pset.Kernel(HandleBoundaries) + pset.Kernel(AdvectionRK4_3D_Lanaau) + pset.Kernel(RWDiffusion_Lanaau) + pset.Kernel(Buoyancy) + pset.Kernel(KillFish) + pset.Kernel(SettleFish) + pset.Kernel(GrowFish)
    kernels = pset.Kernel(KillFish) + pset.Kernel(SettleFish) + pset.Kernel(AdvectionRK4_3D_Lanaau) + pset.Kernel(RWDiffusion_Lanaau) + pset.Kernel(Buoyancy) + pset.Kernel(HandleBoundaries) + pset.Kernel(GrowFish)
    output_file = pset.ParticleFile(name="/10tb_abyss/emily/mhi_model_output/pocillopora-meandrina-nobuoy-run2.nc", outputdt=timedelta(hours = 24))
    counter = 1
    n_events = len(events)

    for duration in events:
        
         print("On " + str(counter) + " out of " + str(n_events) + " steps")

         #if weʻre at an EVEN number, we donʻt spawn
         if (counter % 2) == 0:
             pset.repeatdt = None

         #if weʻre at an ODD number, we do spawn
         else:
             pset.repeatdt = timedelta(minutes=dt).total_seconds()
             #pset.repeat_starttime = runningTime

         #pset.execute(kernels, runtime=duration, dt=timedelta(minutes=dt), output_file=output_file)

         pset.execute(kernels, runtime=duration, dt=timedelta(minutes=dt), output_file=output_file,
                      recovery={ErrorCode.ErrorThroughSurface: submergeParticle,
                                ErrorCode.ErrorOutOfBounds: deleteParticle})

         #pset.execute(kernels, runtime=duration, dt=timedelta(minutes=dt), output_file=output_file,
         #             recovery={ErrorCode.ErrorOutOfBounds: deleteParticle})
         counter += 1

main()
