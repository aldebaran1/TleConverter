#!/usr/bin/env python
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
import ephem,  datetime
from pytz import timezone
import time
from math import sin, cos, asin, atan2, radians, pi, sqrt, degrees
import numpy as np
import yaml
import urllib 
import linecache

def getTLE():
	"""
	Function goes on CelesTrak web pages and greturns TLE data for wanted satellite.
	Function certainly works only for this satellite since web page is .txt file and not
	html. Txt file contains excess spaces (maybe for each column random), so file search 
	is maybe custom for each satellite!
	"""
	#Read .yaml configuration file and parse it into dictionary
	stream = open("inputs.yaml", 'r')
	inputs = yaml.load(stream)
	satellite = str(inputs.get('satellite'))
	#link = str(inputs.get('link'))
	
	#Get tle lines from celestrak page
	source = urllib.urlopen('http://www.celestrak.com/norad/elements/geo.txt').read()
	buff = source.split("\n")
	for num, line in enumerate(buff, 1):
		if satellite in line:
			ix = num
	
	line1 = buff[ix]
	line2 = buff[ix+1]
	
	return line1, line2

def readInputParameters():
	"""
	Function reads the input parametes from selected .yaml file and parses
	them into proper types and reutrns reciever location as list and 
	boolean value of presence WSG84 error in further calculations
	"""
	#Read .yaml configuration file and parse it into dictionary
	stream = open("inputs.yaml", 'r')
	inputs = yaml.load(stream)
	#Set types and store values from dictionary to individual variables
	receiver_lan = float(inputs.get('input lat'))
	receiver_lon = float(inputs.get('input lon'))
	receiver_alt = float(inputs.get('input alt'))
	WSG = bool(inputs.get('WSG84'))
	#Put receiver input parameters into list with proper units
	receiver_location = [receiver_lan, receiver_lon, receiver_alt * 0.001]
	return receiver_location, WSG

def llt2ecef(receiver_location, WSG):
	"""
	Function converts receiver location set in llt format(latitide [degrees], longitude [degrees], altitude [km])
	to ecef (x [km], y [km], z [km]) format and returns it.
	"""
	#Set global variables of receiver position
	global fi, lam, h
	fi = radians(receiver_location[0])
	lam = radians(receiver_location[1])
	h = receiver_location[2]
	#Convert to ecef with WSG84 error parameters of geoid calculations
	if WSG:
		#WSG parameters
		a = 6378.137
		b = 6356.7523142
		esq = 6.69437999014 * 0.001
		e1sq = 6.73949674228 * 0.001
		f = 1 / 298.257223563
		X1 = (a / (sqrt(1 - (2*f - f**2) * sin(fi) * sin(fi) ) ) + h) * cos(fi)
		#Coordinate system convertion
		x = X1 * cos(lam)
		y = X1 * sin(lam)
		z =  ((a * (1 - f)**2 ) / (sqrt(1 - (2*f - f**2) * sin(fi) * sin(fi) ) ) + h) * sin(fi)
		#Put new coordinats in list
		R = [x, y, z]
	#Convert to ecef with fixed earth radious
	else:
		#Earth radious
		a = 6378.137
		#Coordinate system convertion
		x = a * cos(lam) * cos(fi)
		y = a * sin(lam) * cos(fi)
		z = a* sin(fi)
		#Put new coordinats in list
		R = [x, y, z]
	return R

def tle2eci(l1, l2):
	"""
	Function takes first and second line of tle data and returns
	ECI position (x, y, z) of the satellite at current UTC time.
	For parsing TLE data and getting current location I am using
	spg4 library
	"""
	satellite = twoline2rv(l1, l2, wgs84)
	ljubljana = timezone('Europe/Ljubljana')
	lj_time = datetime.datetime.now(ljubljana)
	a_date = ephem.Date(lj_time.strftime('%Y/%m/%d %H:%M:%S'))
	#print "Time in Ljubljana: ", a_date
	
	(year, month, day, hour, minute, secunde) = a_date.tuple()
	position, velocity = satellite.propagate(year, month, day, hour-2, minute, secunde)
	#print "ECI position: ", position

	return position

def getGST(year, month, day, hour, minute, secunde):
	"""
	Function calculates the GST0 in degrees for a current year. This is a constans
	needed for eci2ecef function.
	"""
	
	stream = open("inputs.yaml", 'r')
	inputs = yaml.load(stream)
	longitude = float(inputs.get('input lon'))
	
	julian_date = ephem.julian_date(str(year) + '/' + str(month) + '/' + str(day))
	julian_centuries = (julian_date - 2451545) / 36525
	GST0 = 100.4606184 + (36000.77005361 * julian_centuries) + \
		(0.00038793 * julian_centuries**2) - (2.583*10**-8 * julian_centuries**3)
	GST0 = GST0 % 360
	#print "Julian date: ", julian_date
	#print "GMST0: ", GST0

	#Greenwich sidereal time
	GST = GST0 + 0.25068447733746215 * (hour*60 + minute +secunde/60)
	GST = GST % 360
	#print "GST: ", GST

	#Local sidereal time 
	LST = GST + longitude
	LST = LST % 360
	#print "LST: ", LST

	return GST

def eci2ecef(x, y, z):
	"""
	Function takes ECI location of the satellite (x, y, z) and returns
	its ECEF (x, y, z) location. For convertion between ECI and ECEF
	GST is needed. To get GST, I am using equation peveviously described
	in README file.
	"""

	ljubljana = timezone('Europe/Ljubljana')
	lj_time = datetime.datetime.now(ljubljana)
	a_date = ephem.Date(lj_time.strftime('%Y/%m/%d %H:%M:%S'))
	(year, month, day, hour, minute, secunde) = a_date.tuple()

	#Put input parameters into list
	position = [x, y, z]

	# #Convert GST in radians
	GST = radians(getGST(year, month, day, hour-2, minute, secunde))
	#Rotation matrix for ECI to ECEF twist
	rot = [[cos(GST), sin(GST), 0],
		[-sin(GST), cos(GST), 0],
		[0, 0, 1]]
	#Matrix product of rotation matrix and ECI coordinats to get ECEF coordinats
	ecef_sat = np.dot(rot, position)
	r = sqrt(position[0]**2 + position[1]**2 + position[2]**2)
	#print "ECEF position: ", ecef_sat

	return ecef_sat

def ecef2sez(receiver_position_ecef, satellite_position_ecef):
	"""
	Function takes receiver and satellite locations, both in ECEF format and returns
	satellite location relating to receiver position in ecef format. 
	Function first substract distances to get new satellite position but this is
	calculated to parallel coordinate system. We need to rotate moved coordinate
	system to be rectangluar to tangential surface to get right azimuth and
	elevation angles.
	"""
	position = list(satellite_position_ecef)
	R = list(receiver_position_ecef)
	#Substract position points to get satellite position relating
	#to receiver postion on earth 
	pos = [0, 0, 0]
	pos[0] = ((position[0]) - (R[0]))
	pos[1] = ((position[1]) - (R[1]))
	pos[2] = ((position[2]) - (R[2]))
	#Rotation matrix for coordinate axis twist
	rot = np.array([[sin(fi)*cos(lam), sin(fi)* sin(lam), -cos(fi)], 
		[-sin(lam), cos(lam), 0],
		[cos(fi)*cos(lam), cos(fi)*sin(lam), sin(fi)]])
	#Matrix product to get coordinates rectangular to tangetial surface on 
	#receiver location
	sez = np.dot(rot, pos)
	s = sez[0]
	e = sez[1] 
	z = sez[2]
	#print "SEZ position", s, e, z
	
	return s, e, z

def azimuthAndEvevation(s, e, z):
	"""
	Function takes twisted coordinates and returns azimuth and
	elevation angle of the satellite relating to receviver 
	"""
	distance = sqrt(s**2 + e**2 + z**2)
	El = degrees(asin(z / distance))
	Az = degrees(atan2(e, -s))
	#Azimuth angle is in range [0, 360] degrees
	Az = Az % 360
	print "AER position: ", Az, El, distance
	return Az, El

def getazel():
	"""
	Function starts the calculation cicle and returns final results as a list
	"""
	lines = getTLE()
	inputParameters = readInputParameters()	
	receiver_location_ecef = llt2ecef(inputParameters[0], inputParameters[1])
	sat_eci = tle2eci(lines[0], lines[1])
	satellite_location_ecef = eci2ecef(sat_eci[0], sat_eci[1], sat_eci[2])	
	local_coordinate_system_satellite_position = ecef2sez(receiver_location_ecef, satellite_location_ecef)	
	final_pos = azimuthAndEvevation(local_coordinate_system_satellite_position[0], local_coordinate_system_satellite_position[1],
		local_coordinate_system_satellite_position[2])

	return final_pos

if __name__ == '__main__':
	ae = getazel() 
	