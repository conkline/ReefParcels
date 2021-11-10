from parcels import ParcelsRandom as ParcelsRandom
import math

def HandleBoundaries(particle, fieldset, time):
  
  if particle.isDead == 0 and not particle.stage == 3:
  
    plat = particle.lat
    plon = particle.lon

    #check if we're out of bounds and remove if we are
    if plat <= fieldset.latmin or plat >= fieldset.latmax or plon <= fieldset.lonmin or plon >= fieldset.lonmax:
      particle.OOB = 1
    
    #if not OOB, check for beaching
    #maybe just if last coords = this coords; could also check UVW
    #run this as very last step
    #maybe a while loop?
    else:
      if ((particle.lon == particle.prev_lon) and (particle.lat == particle.prev_lat)) or particle.isBeached == 1:
      
        pdp = particle.depth
        plb = particle.lastBeached
  
        Nc = fieldset.F_N[time, pdp, plat, plon]
        Wc = fieldset.F_W[time, pdp, plat, plon]
  
        #if beached, give it a gentle "push" away from land and towards the sea
        #same scale as diffusion
        D1 = 0.2; R1 = 6378137
        push_lat = (((ParcelsRandom.uniform(0., 1.) * math.sqrt(D1 * particle.dt)) / R1) * (180 / math.pi))
        push_lon = (((ParcelsRandom.uniform(0., 1.) * math.sqrt(D1 * particle.dt)) / (R1 * math.cos(math.pi * particle.lon / 180))) * (180 / math.pi))
        
        if Nc != 0: #sea to the north
          particle.lat += push_lat
        else:
          Sc = fieldset.F_S[time, pdp, plat, plon] #test for sea to the south
          if Sc != 0: 
            particle.lat -= push_lat
        if Wc != 0: #sea to the west
          particle.lon += push_lon
        else:
          Ec = fieldset.F_E[time, pdp, plat, plon] #test for sea to the east
          if Ec != 0: 
            particle.lon += push_lon
  
        #update beaching status
        #keep track of when + where particle last beached
        particle.lastBeached = time
        particle.whereBeachedLat = plat
        particle.whereBeachedLon = plon
  
        #track time beached
        if (time - plb) == particle.dt:
          particle.timeBeached += particle.dt
        else:
          particle.timeBeached = particle.dt
  
        #test to see if we were successful
        (u, v, w) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        if u == 0 and v == 0:
          particle.isBeached = 1
        else:
          particle.isBeached = 0
  
      else:
        particle.isBeached = 0
    
def AdvectionRK4_3D_Lanaau(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration including vertical velocity.
    Function needs to be converted to Kernel object before execution"""
    
    #only run for living, not settled, non beached fish
    if particle.isDead == 0 and not particle.stage == 3 and particle.isBeached == 0 and particle.OOB == 0:
    #   
    #   advect = True
    #   
    #   #have to force interpolation on finer grid around coasts - hacky but works
    #   if particle.lon > 200 and particle.lon < 205.25 and particle.lat > 18.5 and particle.lat < 21.65:
    # 
    #     u1 = fieldset.U1[time, particle.depth, particle.lat, particle.lon]
    #     v1 = fieldset.V1[time, particle.depth, particle.lat, particle.lon]
    #     w1 = fieldset.W1[time, particle.depth, particle.lat, particle.lon]
    #     
    #     lon1 = particle.lon + u1*.5*particle.dt
    #     lat1 = particle.lat + v1*.5*particle.dt
    #     dep1 = particle.depth + w1*.5*particle.dt
    #     
    #     u2 = fieldset.U1[time + .5 * particle.dt, dep1, lat1, lon1]
    #     v2 = fieldset.V1[time + .5 * particle.dt, dep1, lat1, lon1]
    #     w2 = fieldset.W1[time + .5 * particle.dt, dep1, lat1, lon1]
    #     lon2 = particle.lon + u2*.5*particle.dt
    #     lat2 = particle.lat + v2*.5*particle.dt
    #     dep2 = particle.depth + w2*.5*particle.dt
    #     
    #     u3 = fieldset.U1[time + .5 * particle.dt, dep2, lat2, lon2]
    #     v3 = fieldset.V1[time + .5 * particle.dt, dep2, lat2, lon2]
    #     w3 = fieldset.W1[time + .5 * particle.dt, dep2, lat2, lon2]
    #     lon3 = particle.lon + u3*particle.dt
    #     lat3 = particle.lat + v3*particle.dt
    #     dep3 = particle.depth + w3*particle.dt
    #     
    #     u4 = fieldset.U1[time + particle.dt, dep3, lat3, lon3]
    #     v4 = fieldset.V1[time + particle.dt, dep3, lat3, lon3]
    #     w4 = fieldset.W1[time + particle.dt, dep3, lat3, lon3]
    #             
    #     #test for beaching
    #     if u4 == 0 and v4 == 0 and w4 == 0:
    #       advect = False
    #       particle.isBeached = 1
    #    
    #  #otherwise, use nested fields as usual
    #  else:
    
      advect = True
      n_valid_u = 0.
      n_valid_v = 0.
      n_valid_w = 0.

      (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
      lon1 = particle.lon + u1*.5*particle.dt
      lat1 = particle.lat + v1*.5*particle.dt
      dep1 = particle.depth + w1*.5*particle.dt
      if u1 != 0:
        n_valid_u += 1.
      if v1 != 0:
        n_valid_v += 1.
      if w1 != 0:
        n_valid_w += 1.
      
      (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1]
      lon2 = particle.lon + u2*.5*particle.dt
      lat2 = particle.lat + v2*.5*particle.dt
      dep2 = particle.depth + w2*.5*particle.dt
      if u2 != 0:
        n_valid_u += 1.
      if v2 != 0:
        n_valid_v += 1.
      if w2 != 0:
        n_valid_w += 1.
      
      (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2]
      lon3 = particle.lon + u3*particle.dt
      lat3 = particle.lat + v3*particle.dt
      dep3 = particle.depth + w3*particle.dt
      if u3 != 0:
        n_valid_u += 1.
      if v3 != 0:
        n_valid_v += 1.
      if w3 != 0:
        n_valid_w += 1.      
      
      (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3]
      if u4 != 0:
        n_valid_u += 1.
      if v4 != 0:
        n_valid_v += 1.
      if w4 != 0:
        n_valid_w += 1.
        
      #test for beaching
      if u4 == 0 and v4 == 0:
        advect = False
        particle.isBeached = 1
      
      #if we're not beached
      if advect:
        
        #exclude 0-velocity (land) from calculation
        particle.lon += ((u1 + 2*u2 + 2*u3 + u4) / n_valid_u * particle.dt)
        particle.lat += ((v1 + 2*v2 + 2*v3 + v4) / n_valid_v * particle.dt)
        particle.depth += ((w1 + 2*w2 + 2*w3 + w4) / n_valid_w * particle.dt)
        
        #if move would place particle through surface or seafloor, don't make that move
        #test_dep = particle.depth + ((w1 + 2*w2 + 2*w3 + w4) / n_valid_w * particle.dt)
        #move = True
      
        #if test_dep <= 0:
        #  move = False
        #  particle.depth = 0.1
        #else:
        #  ocean_dep = fieldset.seafloor[time, particle.depth, particle.lat, particle.lon]
        #  if test_dep > ocean_dep:
        #    move = False
        #    particle.depth = ocean_dep - 0.1
          
        ##set new depth
        #if move == True:
        #  particle.depth += ((w1 + 2*w2 + 2*w3 + w4) / n_valid_w * particle.dt)

def Behavior(particle, fieldset, time):

    #only run for living, not egg, not settled, non beached fish
    if particle.isDead == 0 and not particle.stage == 0 and not particle.stage == 3 and particle.isBeached == 0 and particle.OOB == 0:
      if fieldset.diel == 1:
        #diel vertical migration - either positive or negative phototaxis
        if particle.stage == 1:
          d1 = ParcelsRandom.uniform(fieldset.larvaeOneDiel, 0) * particle.dt
        elif particle.stage == 2:
          d1 = ParcelsRandom.uniform(fieldset.larvaeOneDiel, 0) * particle.dt
          
        #test for time of day
        secondsThisDay = time % 86400
        if secondsThisDay >= 21600 and secondsThisDay < 64800:
          particle.depth += d1
        else:
          particle.depth -= d1
      
        if fieldset.vertical == 1:
          #vertical swimming - if larvae is below max preferred depth, swim up
          if particle.depth >= fieldset.maxLarvalDepth:
            if particle.stage == 1:
              s1 = ParcelsRandom.uniform(fieldset.larvaeOneVert, 0) * particle.dt
            elif particle.stage == 2:
              s1 = ParcelsRandom.uniform(fieldset.larvaeOneVert, 0) * particle.dt  
            particle.depth += s1
      
      else:
        if fieldset.diel == 1:
          #vertical swimming - if larvae is below max preferred depth or above min preferred, swim up/down
          if particle.depth >= fieldset.maxLarvalDepth:
            if particle.stage == 1:
              s1 = ParcelsRandom.uniform(fieldset.larvaeOneVert, 0) * particle.dt
            elif particle.stage == 2:
              s1 = ParcelsRandom.uniform(fieldset.larvaeOneVert, 0) * particle.dt  
            particle.depth += s1
            
          elif particle.depth <= fieldset.minLarvalDepth:
            if particle.stage == 1:
              s1 = ParcelsRandom.uniform(fieldset.larvaeOneVert, 0) * particle.dt
            elif particle.stage == 2:
              s1 = ParcelsRandom.uniform(fieldset.larvaeOneVert, 0) * particle.dt  
            particle.depth -= s1
            
      #if fieldset.horizontal == 1:
      #  #horizontal swimming - random walk swimming
        
def Buoyancy(particle, fieldset, time):
    
    #only run for living, not settled, non beached fish
    if particle.isDead == 0 and not particle.stage == 3 and particle.isBeached == 0 and particle.OOB == 0:
   
      #add buoyancy
      if particle.stage == 0:
          b1 = ParcelsRandom.uniform(fieldset.eggBuoy, 0)
      elif particle.stage == 1:
          b1 = ParcelsRandom.uniform(fieldset.larvaeOneBuoy, 0)
      elif particle.stage == 2:
          b1 = ParcelsRandom.uniform(fieldset.larvaeTwoBuoy, 0)
          
      particle.depth += (b1 * particle.dt)
      
      #test for breaking the surface or ocean floor
      if particle.depth <=0:
        particle.depth = 0.1
      else:
        ocean_dep = fieldset.seafloor[time, particle.depth, particle.lat, particle.lon]
        if particle.depth >= ocean_dep:
          particle.depth = ocean_dep - 0.1
      
      #if move would place particle through surface or seafloor, don't make that move
      #test_dep = particle.depth + (b1 * particle.dt)
      #move = True
      #
      #if test_dep <= 0:
      #    move = False
      #else:
      #  ocean_dep = fieldset.seafloor[time, particle.depth, particle.lat, particle.lon]
      #  if test_dep > ocean_dep:
      #    move = False
          
      #set new depth
      #if move == True:
      #  particle.depth += (b1 * particle.dt)
         
def RWDiffusion_Lanaau(particle, fieldset, time):
  #only run for living, not settled, non beached fish
  if particle.isDead == 0 and not particle.stage == 3 and particle.isBeached == 0 and particle.OOB == 0:
  
    D = 0.2; R = 6378137 #D = diffusion constant, R = radius of earth in meters
    particle.lat += ((ParcelsRandom.uniform(-1., 1.) * math.sqrt(D * particle.dt)) / R) * (180 / math.pi)
    particle.lon += ((ParcelsRandom.uniform(-1., 1.) * math.sqrt(D * particle.dt)) / (R * math.cos(math.pi * particle.lon / 180))) * (180 / math.pi)

def deleteParticle(particle, fieldset, time):
    particle.OOB = 1
    particle.delete()
    
def submergeParticle(particle, fieldset, time):
    if particle.depth <= 0:
      particle.depth = 0.1
    particle.time = time + particle.dt #so we don't get stuck in infinite loop

def KillFish(particle, fieldset, time):
  
  #only run for living larval fish
  if not particle.stage == 3 and particle.isDead == 0:
    
    #check if weʻre OOB
    if particle.OOB == 1:
      particle.isDead = 1
      particle.howDead = 1
  
    #check if we're over PLD
    if particle.settleStage == 3:
      particle.isDead = 1
      particle.howDead = 4
  
    #check if we've been beached for a long time (24 hrs)
    if particle.timeBeached >= (24*60*60) and particle.isDead == 0:
      particle.isDead = 1
      particle.howDead = 2
  
    #kill by random instantaneous mortality
    if particle.isDead == 0:
        if particle.stage == 0: #egg
            mort = fieldset.eggMort
        elif particle.stage == 1: #larvae1
            mort = fieldset.larvaeOneMort
        elif particle.stage == 2: #larvae2
            mort = fieldset.larvaeTwoMort
          
    if ParcelsRandom.random() <= mort:
        particle.isDead = 1
        particle.howDead = 3
  
    #if particle is dead, remove particle
    if particle.isDead:
      particle.delete() 

def SettleFish(particle, fieldset, time):
  if particle.isDead == 0 and particle.settleStage == 1: #if particle is living and currently settling

    #see which grid we're interpolating on
    #TODO: just make this a nested field 
    #f = fieldset.F[time, particle.depth, particle.lat, particle.lon]
    #
    #if f == 1:
    #  settleTest = fieldset.km1settle[time, particle.depth, particle.lat, particle.lon]
    #elif f == 2:
    #  settleTest = fieldset.km4settle[time, particle.depth, particle.lat, particle.lon]
    
    settleTest = fieldset.settle[time, particle.depth, particle.lat, particle.lon]
  
    #test if we're somewhere we can settle
    #try out coin-flip settlement?
    if (settleTest) > 0 and (ParcelsRandom.uniform(0., 1.) >= 0.5):
        particle.settleStage = 2 #if we do, mark as settled and note location
        particle.final_lon = particle.lon
        particle.final_lat = particle.lat
        particle.final_depth = particle.depth
        particle.settleday = time
        
def GrowFish(particle, fieldset, time):

    #only run for living, not settled fish
    if particle.isDead == 0 and not particle.stage == 3:

        #advance age (in days)
        mydt = particle.dt
        particle.age += mydt / (60*60*24)
    
        #check if we develop
        if particle.stage == 0 and particle.age >= fieldset.hatchTime:
            #egg hatches
            particle.stage = 1
        if particle.stage == 1 and particle.age >= fieldset.flexionTime:
            #larvae advances to stage two
            #this will happen immediately for species without flexion
            particle.stage = 2
      
        #check if we've settled
        if particle.settleStage == 2:
            particle.stage = 3
      
        #check if we need to update settleStage (past PLD)
        if particle.age >= fieldset.maxPLD and particle.settleStage != 2:
            particle.settleStage = 3
        elif particle.settleStage == 0 and fieldset.minPLD <= particle.age and particle.age <= fieldset.maxPLD:
            particle.settleStage = 1
      
        #advance approx distance travelled (Euclidian)
        lat_dist = (particle.lat - particle.prev_lat) * 1.11e2
        lon_dist = (particle.lon - particle.prev_lon) * 1.11e2 * math.cos(particle.lat * math.pi / 180)
        particle.distance += math.sqrt(math.pow(lon_dist, 2) + math.pow(lat_dist, 2))
        
        if not particle.settleStage == 2:
            particle.prev_lon = particle.lon  # Set the stored values for next iteration
            particle.prev_lat = particle.lat
            particle.prev_depth = particle.depth

# def GetLandPolys():
#     myfile = open("/home/emily/Desktop/connectivity_5_13_19/land_simplified.csv", "r")
# 
#     groups = []
#     points = []
# 
#     for line in myfile.readlines()[1::]:
#         parts = line.split(",")
#         x = float(parts[1]); y = float(parts[2])
#         g = parts[7].strip("\n")
# 
#         groups.append(g)
#         points.append(Point(x,y))
# 
#     polygons = []
#     last = ""
#     thispoly = []
# 
#     for point in range(len(points)):
#         if last!="":
#             if last!=groups[point]:
#                 polygons.append(thispoly)
#                 thispoly = []
# 
#         thispoly.append(points[point])
#         last = groups[point]
# 
#     polygons.append(thispoly)
#     landpoly = [Polygon(sum(map(list, (p.coords for p in pts)), [])) for pts in polygons]
# 
#     return(landpoly)
# 
# def IsOnLand(particle):
#     test_point = Point(particle.lon, particle.lat)
#     onland = any([x.contains(test_point) for x in landpoly])
#     return(onland)
# 
# def MoveOffLand_old(particle):
#     if IsOnLand(particle):
#         if particle.prev_lon == 0:
#             dist_idx = findNearestPoint([particle.lon, particle.lat])
#                     new_point = gridPointsMasked[dist_idx[1]]
# 
#                     #set new lon and lat
#                     D = 0.2; R = 6378137
#             new_lon = new_point[0] + ((random.uniform(-1, 1)*math.sqrt(D*particle.dt)) / R) * (180/3.14)
#                     new_lat = new_point[1] + ((random.uniform(-1, 1)*math.sqrt(D*particle.dt)) / (R*math.cos(3.14*new_point[1]/180)) * (180/3.14))
#                     #particle.update()
#         else:
#             D = 0.2; R = 6378137
#             new_lon = particle.prev_lon + ((random.uniform(-1, 1)*math.sqrt(D*particle.dt)) / R) * (180/3.14)
#             new_lat = particle.prev_lat + ((random.uniform(-1, 1)*math.sqrt(D*particle.dt)) / (R*math.cos(3.14*particle.prev_lat/180)) * (180/3.14))
# 
#         particle.lon = new_lon
#         particle.lat = new_lat
# 
#         loopCount = 0
#         while IsOnLand(particle):
#             particle.lon += ((random.uniform(-1, 1)*math.sqrt(D*particle.dt)) / R) * (180/3.14)
#             particle.lat += ((random.uniform(-1, 1)*math.sqrt(D*particle.dt)) / (R*math.cos(3.14*particle.prev_lat/180)) * (180/3.14))
#             loopCount += 1
#             if loopCount > 100:
#                 print("Looks like particle " + str(particle.id) + "is stuck - deleting")
#                 particle.delete()
#                 break
# 
#     else:
#         print("Particle " + str(particle.id) + " stuck but not on land: deleting")
#         particle.delete()
# 
# 
# def MoveOffLand(particle, fieldset, time):
# 
#     D = 0.2; R = 6378137
# 
#     if particle.prev_lon == 0:
#         dist_idx = findNearestPoint([particle.lon, particle.lat])
#         new_point = gridPointsMasked[dist_idx[1]]
#         new_lon = new_point[0] + ((random.uniform(-1, 1)*math.sqrt(D*particle.dt)) / R) * (180/3.14)
#         new_lat = new_point[1] + ((random.uniform(-1, 1)*math.sqrt(D*particle.dt)) / (R*math.cos(3.14*new_point[1]/180)) * (180/3.14))
# 
#     else:
#         new_lon = particle.prev_lon + ((random.uniform(-1, 1)*math.sqrt(D*particle.dt)) / R) * (180/3.14)
#         new_lat = particle.prev_lat + ((random.uniform(-1, 1)*math.sqrt(D*particle.dt)) / (R*math.cos(3.14*particle.prev_lat/180)) * (180/3.14))
# 
#     #set new lon and lat
#     particle.lon = new_lon
#     particle.lat = new_lat
#     particle.depth = 0
