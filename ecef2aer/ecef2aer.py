from math import  sqrt, asin, atan2, degrees, radians, atan
from scipy import sin, cos
import numpy as np
from numpy import linalg as la
import yaml

def getInputs():
	"""
	Function reads, parses, converts and returns all input parameters from .yaml configuration file
	"""
	stream = open("inputs.yaml", 'r')
	inputs = yaml.load(stream)
	input_list_file = inputs.get('input file name')
	output_list_file = inputs.get('output file name')
	receiver_lan = float(inputs.get('input lat'))
	receiver_lon = float(inputs.get('input lon'))
	receiver_alt = float(inputs.get('input alt'))
	WSG = bool(inputs.get('WSG84'))
	return input_list_file, output_list_file, receiver_lan, receiver_lon, receiver_alt, WSG

def ecef2aer(input_list_file, output_list_file, receiver_lan, receiver_lon, receiver_alt, WSG):
	"""
	Function takes all input parameters, does a conversion from ecef to aer and stores
	all converted data to file.
	"""
	#Put receiver parameters into list
	Rin = [receiver_lan, receiver_lon, receiver_alt * 0.001]
	#Read, edit and store oem data in list of lists as input parameter
	f = open(input_list_file, 'r')
	src_data = []
	src_data = [line.split() for line in f] 
	for i in src_data:
		#print i
		i.pop(6)
		i.pop(5)
		i.pop(4)
		i.pop(0)
	f.close()

	#base ecef vectors def. and rotation matrix definition
	fi = radians(Rin[0])
	lam = radians(Rin[1])
	h = Rin[2]
	ijk = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
	rot = np.array([[sin(fi)*cos(lam), sin(fi)* sin(lam), -cos(fi)], 
		[-sin(lam), cos(lam), 0],
		[cos(fi)*cos(lam), cos(fi)*sin(lam), sin(fi)]])

	#WSG84 parameters needed to calculate exact earth radious
	#depends to receiver position
	#a = 6378.137
	a = 6378.135
	b = 6356.7523142
	esq = 6.69437999014 * 0.001
	e1sq = 6.73949674228 * 0.001
	#f = 1 / 298.257223563
	f = 1 / 298.26

	#llt to xyz convetion of raceviver location R, including WSG84 parameters in 
	#geodetic to cartesian coo. system conversion.
	#X1 is WSG84 pre-calculated parameter further used in x and y equation
	if WSG:
		X1 = (a / (sqrt(1 - (2*f - f**2) * sin(fi) * sin(fi) ) ) + h) * cos(fi)

		x = X1 * cos(lam)
		y = X1 * sin(lam)
		z =  ((a * (1 - f)**2 ) / (sqrt(1 - (2*f - f**2) * sin(fi) * sin(fi) ) ) + h) * sin(fi)
		#Put new coordinats in list
		R = [x, y, z]
	else:
		x = a * cos(lam) * cos(fi)
		y = a * sin(lam) * cos(fi)
		z = a* sin(fi)
		#Put new coordinats in list
		R = [x, y, z]

	#Calculation of satelite vector between receiver and satelite position - difference of ecef basic sat. position and ecef basic receiver position
	ro = src_data
	roABS = []
	for i in ro:
		i[0] = float(i[0]) - (R[0])
		i[1] = float(i[1]) - (R[1])
		i[2] = float(i[2]) - (R[2])
		abs1 = sqrt(i[0]**2 + i[1]**2 + i[2]**2)
		roABS.append(abs1)

	#New basic vectors for receiver based coo. system
	sezABS = []
	s = []
	e = []
	z = []
	for i in range(len(ro)):
		sez = np.dot(rot, ro[i])
		s.append(sez[0])
		e.append(sez[1])
		z.append(sez[2])
		sezABS.append(sqrt(sez[0]**2 + sez[1]**2 + sez[2]**2))

	#Calculate azimuth and elevation angles and stor them to list
	El = []
	Az = []
	OutData = []
	for i in range(len(s)):
		El1 = degrees(asin(z[i]/roABS[i]))
		Az1 = degrees(atan2(e[i], -s[i]))
		El.append(El1)
		Az.append(360 + Az1)
		a = [Az1, El1]
		OutData.append(a)
	#print max(roABS), min(roABS)

	#print max(Az), min(Az)
	#print max(El), min(El)
	
	#Write data to file
	g = open(output_list_file, "w")
	for i  in OutData:
		for k in i:
			g.write(str(k))
			g.write(" ")
		g.write("\n")
	g.close()

if __name__ == '__main__':
	inputs = getInputs()
	ecef2aer(inputs[0], inputs[1], inputs[2], inputs[3], inputs[4], inputs[5])


