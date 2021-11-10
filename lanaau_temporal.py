import datetime
import time
from pyorbital_moonphase import moon_phase

#get model day from date
def getModelDay(modeldate):
  startDate = datetime.datetime(2011, 4, 1)
  dayOut = (modeldate - startDate).days
  return(dayOut)

#get month from model day
def getDate(modelday):
  #modelday is 1:852
  #starting on March 30th 2011 and ending on July 30th 2013
  
  startDate = datetime.datetime(2011, 4, 1)
  timeSince = datetime.timedelta(days = modelday)
  dateOut = startDate + timeSince
  
  #returns date object
  return dateOut

#get moon phase from model date
def getMoon(modeldate):
  #options: full, waning gibbous, waning half, waning crescent, new, waxing crescent, waxing half, waxing gibbous
  currentMoon = moon_phase(modeldate)
    
  #figure out phase
  if currentMoon >= 0.95:
    phase = "Full"
  elif currentMoon <= 0.02:
    phase = "New"
  else:
    yesterdayMoon = moon_phase(modeldate - datetime.timedelta(days = 1))
  
    #figure out if waxing or waning
    if yesterdayMoon > currentMoon:
      modif = "Waning"
    else:
      modif = "Waxing"
    
    #figure out phase
    if 0.02 < currentMoon <= 0.4:
      tmpphase = " crescent"
    elif 0.4 < currentMoon <= 0.65:
      tmpphase = " half"
    elif 0.65 < currentMoon < 0.95:
      tmpphase = " gibbous"
      
    phase = modif + tmpphase
    
  return(phase)
    
#def time of day from timestep
def getTimeOfDay(seconds):
  secondsThisDay = datetime.timedelta(seconds = seconds).seconds
  currTime = time.strftime('%H:%M:%S', time.gmtime(secondsThisDay))
  hm = [int(x) for x in currTime.split(":")[0:2]]
  
  #return [hours, minutes]
  return hm
