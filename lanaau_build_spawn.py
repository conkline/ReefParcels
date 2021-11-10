import netCDF4 as nc
import numpy as np
import random
import datetime
from lanaau_temporal import getModelDay, getDate, getMoon, getTimeOfDay

def getSpawningSites(spawn_1km_file, spawn_4km_file, depth_1km_file, depth_4km_file, maxDepth):
  #read in files
  spawn_1km = nc.Dataset(spawn_1km_file)
  spawn_4km = nc.Dataset(spawn_4km_file)
  depth_1km = nc.Dataset(depth_1km_file)
  depth_4km = nc.Dataset(depth_4km_file)
  
  #set how many times we release a particle from a cell per time step
  n_particles = 10
              
  #generate random spawn sites within each habitat cell (1km)
  coordx = []; coordy = []
  spawn_array_1km = spawn_1km["1kmspawn"]
  for x in range(0, (spawn_array_1km.shape[1])):
    for y in range(0, (spawn_array_1km.shape[0])):
      if spawn_array_1km[y,x] != "--":
        for i in range(n_particles):
          coordx.append(x)
          coordy.append(y)
    
  n = len(coordx)    
  dx = (max(spawn_1km["lon"]) - min(spawn_1km["lon"])) / (spawn_array_1km.shape[1] * 2)
  dy = (max(spawn_1km["lat"]) - min(spawn_1km["lat"])) / (spawn_array_1km.shape[0] * 2)
  lon = spawn_1km["lon"][coordx].tolist()
  lat = spawn_1km["lat"][coordy].tolist()
  newlat_1km = (lat + np.random.uniform(low = -dy, high = dy, size = n)).tolist()
  newlon_1km = (lon + np.random.uniform(low = -dx, high = dx, size = n)).tolist()

  #randomly assign starting depths (1km)
  rand_depths_1km = []
  for i in range(n):
    try:
      rd = random.randrange(0, maxDepth+1, 1)
    except ValueError:
      rd = maxDepth #for opihi at surface
      
    depth = depth_1km["1kmdepth"][coordy[i], coordx[i]].tolist()
    if rd > depth:
      rd = depth
    rand_depths_1km.append(rd)
    
  #generate random spawn sites within each habitat cell (4km)
  coordx = []; coordy = []
  spawn_array_4km = spawn_4km["4kmspawnsettle"]
  lon_all = spawn_4km["lon"][0:spawn_array_4km.shape[1]].tolist()
  lat_all = spawn_4km["lat"][0:spawn_array_4km.shape[0]].tolist()
  max_y = max(spawn_1km["lat"]).tolist()
  min_x = min(spawn_1km["lon"]).tolist()
  
  for x in range(0, (spawn_array_4km.shape[1])):
    for y in range(0, (spawn_array_4km.shape[0])):
      if spawn_array_4km[y,x] != "--":
        if lon_all[x] < min_x or lat_all[y] > max_y:#subset by resolution extent
          for i in range(4*n_particles): #four sites for every 4km square
            coordx.append(x)
            coordy.append(y)
          
  n = len(coordx)   
  dx = (max(spawn_4km["lon"]) - min(spawn_4km["lon"])) / (spawn_array_4km.shape[1] * 2)
  dy = (max(spawn_4km["lat"]) - min(spawn_4km["lat"])) / (spawn_array_4km.shape[0] * 2)
  lon = spawn_4km["lon"][coordx].tolist()
  lat = spawn_4km["lat"][coordy].tolist()
  newlat_4km = (lat + np.random.uniform(low = -dy, high = dy, size = n)).tolist()
  newlon_4km = (lon + np.random.uniform(low = -dx, high = dx, size = n)).tolist()
        
  #randomly assign starting depths (4km)
  rand_depths_4km = []
  for i in range(n):
    try:
      rd = random.randrange(0, maxDepth+1, 1)
    except ValueError:
      rd = maxDepth #for opihi at surface
    depth = depth_4km["4kmdepth"][coordy[i], coordx[i]].tolist()
    if rd > depth:
      rd = depth
    rand_depths_4km.append(rd)
    
  newlat_1km.extend(newlat_4km)
  newlon_1km.extend(newlon_4km)
  rand_depths_1km.extend(rand_depths_4km)

  return(newlat_1km, newlon_1km, 
         rand_depths_1km)
  

def getSpawningEvents(monthSpawner, spawnMonths, timedSpawner, timeSpawnBegin, \
  timeSpawnEnd, peakSpawner, peakSpawnMonths, moonSpawner, spawnPhase):
    
  startTimestamp = getDate(0) #model starts
  endTimestamp = getDate(851) #model stops
  
  calendar = [startTimestamp + datetime.timedelta(days=x) for x in range((endTimestamp - startTimestamp).days + 1)]
  
  spawningdays = []
  peakspawningdays = []
  addSpawn = True
  addPeak = True
  
  for tmpDay in calendar:
    if monthSpawner and tmpDay.month not in spawnMonths: #if we donʻt spawn this month
      addSpawn = False
    if (not peakSpawner) or (peakSpawner and tmpDay.month not in peakSpawnMonths): #if we donʻt peak spawn this month
      addPeak = False
    if (not moonSpawner) or (addSpawn and moonSpawner and getMoon(tmpDay) not in spawnPhase): #if we donʻt spawn this moon phase
      addSpawn = False
      addPeak = False
      
    if addSpawn == True:
      spawningdays.append(tmpDay)
    if addPeak == True:
      peakspawningdays.append(tmpDay)
      
    addSpawn = True
    addPeak = True
    
  #month, no peak, no time of day
  #C. exarata
  if monthSpawner and not peakSpawner and not timedSpawner:
    #generate one random spawning event per site per spawning day
    spawn_starts = [x + datetime.timedelta(seconds = random.randrange(0, 86401, 1)) for x in spawningdays]
    spawn_stops = [x + datetime.timedelta(hours = 1) for x in spawn_starts]
    
  #month, no peak, time of day
  #P. meandrina
  if monthSpawner and not peakSpawner and timedSpawner:
    #add time to spawning days
    spawn_starts = [x + timeSpawnBegin for x in spawningdays]
    spawn_stops = [x + timeSpawnEnd for x in spawningdays]
      
  #build list of [spawning duration, no spawning duration, spawning duration, etc.]
  duration_spawning = []
  lastStop = ""
    
  for n in range(len(spawn_starts)):
    this_start = spawn_starts[n]
    this_stop = spawn_stops[n]
    
    if lastStop != "":
      #add no spawning duration
      duration_spawning.append(this_start - lastStop)
      
    #add spawning duration
    duration_spawning.append(this_stop - this_start)
    
    lastStop = this_stop
      
  #add extra chunk of time to the end
  starttime = (spawningdays[0] - startTimestamp).total_seconds()
  duration_spawning.append((endTimestamp - startTimestamp) - sum(duration_spawning, datetime.timedelta()))
        
  return(starttime, duration_spawning)
